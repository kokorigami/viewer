var createContext = require('gl-context');
var createShader = require('gl-shader');
var createFBO = require('gl-fbo');
var createCamera = require('3d-view');

var createFaceBufferInfo = require('./createFaceBufferInfo.js');
var createFoldBufferInfo = require('./createFoldBufferInfo.js');
var createQuadBufferInfo = require('./createQuadBufferInfo.js');

var glslify = require('glslify');
var faceFs = glslify('./face-fs.glsl');
var faceVs = glslify('./face-vs.glsl');
var foldFs = glslify('./fold-fs.glsl');
var foldVs = glslify('./fold-vs.glsl');
var depthnormalFs = glslify('./depthnormal-fs.glsl');
var ssaoFs = glslify('./ssao-fs.glsl');
var ssaoVs = glslify('./ssao-vs.glsl');

var m4 = require('gl-mat4');

var Renderer = function (canvas, data) {
  this.render = this.render.bind(this);
  this.camera = createCamera({
    center: [0, 0.5, 0],
    eye: [0, 1, -3],
    distanceLimits: [1, 100],
    up: [0, 0, 1],
    mode: 'orbit'
  });

  this.initialize(canvas);

  if (data) this.data(data);
  return this;
};

module.exports = Renderer;

Renderer.prototype.initialize = function (canvas) {
  var gl = this.gl = createContext(canvas);
  gl.enable(gl.DEPTH_TEST);

  this.shaders = {};
  this.shaders.face = createShader(gl, faceVs, faceFs);
  this.shaders.depth = createShader(gl, faceVs, depthnormalFs);
  this.shaders.fold = createShader(gl, foldVs, foldFs);
  this.shaders.ssao = createShader(gl, ssaoVs, ssaoFs);

  this.fbo = createFBO(gl, [0, 0]);

  this.onMouseDown = getOnMouseDown(gl, this.camera);

  this.planeBufferInfo = createQuadBufferInfo(gl);
};

Renderer.prototype.teardown = function () {
  this.gl = null;
  if (this.fbo) {
    this.fbo.dispose();
  }
  var shaders = Object.keys(this.shaders);
  shaders.forEach(function (name) {
    this.shaders[name].dispose(); // TODO: check this line
  }.bind(this));
};

Renderer.prototype.data = function (data) {
  var gl = this.gl;
  var facesPerLayer = data.layers;
  var foldsPerLayer = data.folds;
  // var points = data.points;

  this.faceBufferInfo = createFaceBufferInfo(gl, facesPerLayer);
  this.foldBufferInfo = createFoldBufferInfo(gl, foldsPerLayer);
  return this;
};

Renderer.prototype.play = function () {
  this.stop();
  this.gl.canvas.addEventListener('mousedown', this.onMouseDown);
  this.frame = requestAnimationFrame(this.render);
  return this;
};

Renderer.prototype.render = function () {
  this.resize();
  var gl = this.gl;

  var uniforms = this.getUniforms();
  var fbo = this.getFramebuffer();

  this.bindFramebuffer(fbo);
  renderPass(gl, this.shaders.depth, this.faceBufferInfo, uniforms, 'TRIANGLES');

  this.bindFramebuffer(null);
  renderPass(gl, this.shaders.ssao, this.planeBufferInfo, uniforms, 'TRIANGLES');
  // renderPass(gl, this.shaders.face, this.faceBufferInfo, uniforms, 'TRIANGLES');
  // renderPass(gl, this.shaders.fold, this.foldBufferInfo, uniforms, 'LINES');

  //renderPoints(points, pointDivs);
  requestAnimationFrame(this.render);
};

Renderer.prototype.resize = function (resolution) {
  resolution = resolution || 1;
  resolution = Math.max(1, resolution);
  var canvas = this.gl.canvas;
  var width  = canvas.clientWidth  * resolution | 0;
  var height = canvas.clientHeight * resolution | 0;
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
    return true;
  }

  this.fbo.shape = [width, height];
  return false;
};

Renderer.prototype.getUniforms = function () {
  if (!this.uniforms) {
    this.uniforms = this.createUniforms();
  }
  var gl = this.gl;
  var t = Date.now();
  var camera = this.camera;

  camera.idle(t - 20);
  camera.flush(t - 100);
  camera.recalcMatrix(t - 25);

  // Update camera uniforms.
  var view = camera.computedMatrix;
  var projection = m4.perspective(
      [],
      Math.PI/4.0,
      gl.drawingBufferWidth/gl.drawingBufferHeight,
      1, 5
  );

  var width = gl.drawingBufferWidth;
  var height = gl.drawingBufferHeight;

  this.uniforms.u_screen = [width, height];
  this.uniforms.u_near = 1;
  this.uniforms.u_far = 5;
  this.uniforms.u_view = view;
  this.uniforms.u_camera = m4.invert([], view);
  this.uniforms.u_worldView = m4.multiply([], view, this.uniforms.u_world);
  this.uniforms.u_worldViewProjection = m4.multiply([], projection, this.uniforms.u_worldView);

  return this.uniforms;
};

Renderer.prototype.createUniforms = function () {
  var uniforms = {
    u_lightWorldPos: [-3, 3, -8],
    u_lightColor: [1, 0.9, 0.8, 1],
    u_ambient: [0, 0, 0, 1],
    u_shininess: 70,
    u_world: m4.identity([])
  };
  return uniforms;
};

Renderer.prototype.bindFramebuffer = function (fbo) {
  var gl = this.gl;

  if (fbo) fbo.bind();
  else gl.bindFramebuffer(gl.FRAMEBUFFER, null);

  gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
};

Renderer.prototype.stop = function () {
  cancelAnimationFrame(this.frame);
  this.gl.canvas.removeEventListener('mousedown', this.onMouseDown);
  return this;
};

function renderPass (gl, shader, geom, uniforms, drawType) {
  shader.bind();
  shader.uniforms = uniforms;
  geom.bind(shader);
  geom.draw(gl[drawType]);
  geom.unbind();
}

function getOnMouseDown (gl, camera) {
  var onMouseDown = function (ev) {
    var canvas = ev.currentTarget;
    canvas.addEventListener('mousemove', onMouseMove);
    canvas.addEventListener('mouseup', onMouseUp);

    var startX = ev.offsetX;
    var startY = ev.offsetY;
    var currentX = startX;
    var currentY = startY;

    function onMouseMove (mev) {
      var dx = (mev.offsetX - currentX) * 4/ gl.drawingBufferWidth;
      var dy = (mev.offsetY - currentY) * 4/ gl.drawingBufferHeight;

      camera.rotate(Date.now(), -dx, dy);
      currentX = mev.offsetX;
      currentY = mev.offsetY;
    }

    function onMouseUp () {
      canvas.removeEventListener('mousemove', onMouseMove);
      canvas.removeEventListener('mouseup', onMouseUp);
    }
  };
  return onMouseDown;
}

// var pointDivs = createPointDivs(points);

// function toPixelClipSpace (gl, point) {
//   var pixel = [];
//   pixel[0] = (point[0] *  0.5 + 0.5) * gl.canvas.width;
//   pixel[1] = (point[1] * -0.5 + 0.5) * gl.canvas.height;
//   return pixel;
// }

// function createPointDivs (points) {
//   return _.map(points, createPoint);
// }

// function renderPoints(points, divs) {
//   _.each(points, function (point, i) {
//     renderPoint(point, divs[i]);
//   });
// }

// function createPoint (point, id) {
//   var pointDiv = document.createElement('div');
//   pointDiv.classList.add('floating');
//   pointDiv.textContent = id;
//   document.body.appendChild(pointDiv);
//   return pointDiv;
// }

// function renderPoint (point, div) {
//   var adjustedPoint = twgl.v3.create();
//   adjustedPoint[0] = point[0];
//   adjustedPoint[1] = point[1];
//   adjustedPoint[2] = point[2];

//   twgl.m4.transformPoint(worldViewProjection, adjustedPoint, adjustedPoint);

//   var pixelPoint = toPixelClipSpace(gl, adjustedPoint);

//   div.style.left = Math.floor(pixelPoint[0]) + 'px';
//   div.style.top = Math.floor(pixelPoint[1]) + 'px';
// }
