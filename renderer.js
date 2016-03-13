var createCamera = require('3d-view');

var glslify = require('glslify');
var faceFs = glslify('./face-fs.glsl');
var faceVs = glslify('./face-vs.glsl');
var foldFs = glslify('./fold-fs.glsl');
var foldVs = glslify('./fold-vs.glsl');
var depthnormalFs = glslify('./depthnormal-fs.glsl');

var twgl = require('twgl.js');
var m4 = require('gl-mat4');

var createFaceBufferInfo = require('./createFaceBufferInfo.js');
var createFoldBufferInfo = require('./createFoldBufferInfo.js');

var Renderer = function (canvas, data) {
  var gl = this.gl = twgl.getWebGLContext(canvas);
  this.render = this.render.bind(this);

  var faceProgram = twgl.createProgramFromSources(gl, [faceVs, faceFs]);
  var depthProgram = twgl.createProgramFromSources(gl, [faceVs, depthnormalFs]);
  var foldProgram = twgl.createProgramFromSources(gl, [foldVs, foldFs]);

  this.faceProgramInfo = twgl.createProgramInfoFromProgram(gl, faceProgram);
  this.depthProgramInfo = twgl.createProgramInfoFromProgram(gl, depthProgram);
  this.foldProgramInfo = twgl.createProgramInfoFromProgram(gl, foldProgram);

  this.camera = createCamera({
    center: [0, 0.5, 0],
    eye: [0, 1, -3],
    distanceLimits: [1, 100],
    up: [0, 0, 1],
    mode: 'orbit'
  });

  //this.rotation = [0, 0];
  this.onMouseDown = getOnMouseDown(gl, this.camera);
  if (data) {
    this.data(data);
  }
  return this;
};

module.exports = Renderer;

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
  // Setup listeners
  this.stop();
  this.gl.canvas.addEventListener('mousedown', this.onMouseDown);
  this.frame = requestAnimationFrame(this.render);
  return this;
};

Renderer.prototype.render = function () {
  var gl = this.gl;

  twgl.resizeCanvasToDisplaySize(gl.canvas);
  gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);

  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  gl.enable(gl.DEPTH_TEST);

  var uniforms = this.getUniforms();

  //renderPass(gl, depthProgramInfo, faceBufferInfo, uniforms, 'TRIANGLES');
  renderPass(gl, this.faceProgramInfo, this.faceBufferInfo, uniforms, 'TRIANGLES');
  renderPass(gl, this.foldProgramInfo, this.foldBufferInfo, uniforms, 'LINES');

  //renderPoints(points, pointDivs);
  requestAnimationFrame(this.render);
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
      1,
      100
  );

  var viewProjection = m4.multiply([], projection, view);
  var worldViewProjection = m4.multiply([], viewProjection, this.uniforms.u_world);

  this.uniforms.u_near = 1;
  this.uniforms.u_far = 100;
  this.uniforms.u_view = view;
  this.uniforms.u_camera = m4.invert([], view);
  this.uniforms.u_worldViewProjection = worldViewProjection;

  return this.uniforms;
};

Renderer.prototype.createUniforms = function () {
  var world = m4.identity([]);

  var uniforms = {
    u_lightWorldPos: [-3, 3, -8],
    u_lightColor: [1, 0.9, 0.8, 1],
    u_ambient: [0, 0, 0, 1],
    u_specular: [1, 1, 0.8, 1],
    u_shininess: 70,
    u_specularFactor: 0.8,
    u_worldRotation: world,
    u_world: world
  };
  return uniforms;
};

Renderer.prototype.stop = function () {
  // TODO: remove listeners animation frames
  cancelAnimationFrame(this.frame);
  this.gl.canvas.removeEventListener('mousedown', this.onMouseDown);
  return this;
};

//var pointDivs = createPointDivs(points);

function renderPass (gl, programInfo, bufferInfo, uniforms, drawType) {
  gl.useProgram(programInfo.program);
  twgl.setBuffersAndAttributes(gl, programInfo, bufferInfo);
  twgl.setUniforms(programInfo, uniforms);
  twgl.drawBufferInfo(gl, gl[drawType], bufferInfo);
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

    function onMouseMove (ev) {
      var dx = (ev.offsetX - currentX) * 4/ gl.drawingBufferWidth;
      var dy = (ev.offsetY - currentY) * 4/ gl.drawingBufferHeight;

      camera.rotate(Date.now(), -dx, dy);
      currentX = ev.offsetX;
      currentY = ev.offsetY;
    }

    function onMouseUp () {
      canvas.removeEventListener('mousemove', onMouseMove);
      canvas.removeEventListener('mouseup', onMouseUp);
    }
  };
  return onMouseDown;
}

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
