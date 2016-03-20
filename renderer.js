var glslify = require('glslify');
var faceFs = glslify('./face-fs.glsl');
var faceVs = glslify('./face-vs.glsl');
var foldFs = glslify('./fold-fs.glsl');
var foldVs = glslify('./fold-vs.glsl');
var depthnormalFs = glslify('./depthnormal-fs.glsl');
var ssaoFs = glslify('./ssao-fs.glsl');
var ssaoVs = glslify('./ssao-vs.glsl');

var _ = require('underscore')._;
var twgl = require('twgl.js');
var v3 = twgl.v3;
var m4 = twgl.m4;

var createFaceBufferInfo = require('./createFaceBufferInfo.js');
var createFoldBufferInfo = require('./createFoldBufferInfo.js');

var Renderer = function (canvas, data) {
  var gl = this.gl = twgl.getWebGLContext(canvas);
  gl.enable(gl.DEPTH_TEST);
  this.render = this.render.bind(this);

  var faceProgram = twgl.createProgramFromSources(gl, [faceVs, faceFs]);
  var depthProgram = twgl.createProgramFromSources(gl, [faceVs, depthnormalFs]);
  var ssaoProgram = twgl.createProgramFromSources(gl, [ssaoVs, ssaoFs]);
  var foldProgram = twgl.createProgramFromSources(gl, [foldVs, foldFs]);

  this.faceProgramInfo = twgl.createProgramInfoFromProgram(gl, faceProgram);
  this.foldProgramInfo = twgl.createProgramInfoFromProgram(gl, foldProgram);
  this.depthProgramInfo = twgl.createProgramInfoFromProgram(gl, depthProgram);
  this.ssaoProgramInfo = twgl.createProgramInfoFromProgram(gl, ssaoProgram);

  this.planeBufferInfo = twgl.primitives.createXYQuadBufferInfo(gl);

  this.framebuffers = [];
  this.fbIndex = 0;

  this.rotation = [0, 0];
  this.onMouseDown = getOnMouseDown(this.rotation);
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

  var uniforms = this.getUniforms();

  //this.setFramebuffer(null);
  this.swapFramebuffer();
  renderPass(gl, this.depthProgramInfo, this.faceBufferInfo, uniforms, 'TRIANGLES');

  this.setFramebuffer(null);
  renderPass(gl, this.ssaoProgramInfo, this.planeBufferInfo, uniforms, 'TRIANGLES');
  //renderPass(gl, this.faceProgramInfo, this.faceBufferInfo, uniforms, 'TRIANGLES');
  // renderPass(gl, this.foldProgramInfo, this.foldBufferInfo, uniforms, 'LINES');


  //renderPoints(points, pointDivs);
  requestAnimationFrame(this.render);
};

Renderer.prototype.getUniforms = function () {
  var gl = this.gl;
  var eye = [0, 0.5, -4];
  var target = [0, 0.5, 0];
  var up = [0, 1, 0];
  var near = 1;
  var far = 6;

  var width = gl.canvas.clientWidth;
  var height = gl.canvas.clientHeight;
  var projection = m4.perspective(30 * Math.PI / 180, width / height, near, far);
  var camera = m4.lookAt(eye, target, up);
  var view = m4.inverse(camera);
  var viewProjection = m4.multiply(view, projection);

  var rotateX = m4.rotationX(this.rotation[0]);
  var rotateY = m4.rotationY(this.rotation[1]);
  var world = m4.multiply(rotateX, rotateY);

  var worldView = m4.multiply(world, view);
  var worldViewProjection = m4.multiply(world, viewProjection); // projection, view, world (of object), point (of object)

  var uniforms = {
    u_near: near,
    u_far: far,
    u_lightWorldPos: [-3, 3, -8],
    u_lightColor: [1, 0.9, 0.8, 1],
    u_ambient: [0, 0, 0, 1],
    u_specular: [1, 1, 0.8, 1],
    u_shininess: 70,
    u_specularFactor: 0.8,
    u_camera: camera,
    u_view: view,
    u_viewProjection: viewProjection,
    u_worldViewProjection: worldViewProjection,
    u_worldRotation: worldView,
    u_world: world,
    u_screen: [width, height]
  };
  return uniforms;
};

Renderer.prototype.swapFramebuffer = function () {
  var gl = this.gl;
  var attachments = [
    {format: gl.RGBA, type: gl.UNSIGNED_BYTE, min: gl.LINEAR, wrap: gl.CLAMP_TO_EDGE}
  ];

  if (!this.framebuffers[this.fbIndex]) {
    this.framebuffers[this.fbIndex] = twgl.createFramebufferInfo(gl, attachments);
  }

  var framebuffer = this.framebuffers[this.fbIndex];
  this.setFramebuffer(framebuffer);
  return this.fbIndex = 1 - this.fbIndex;
};

Renderer.prototype.setFramebuffer = function (fbo) {
  var gl = this.gl;
  var framebuffer = fbo && fbo.framebuffer || null;
  gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
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

function getOnMouseDown (rotation) {
  var onMouseDown = function (e) {
    var canvas = e.currentTarget;
    canvas.addEventListener('mousemove', onMouseMove);
    canvas.addEventListener('mouseup', onMouseUp);

    var startX = e.offsetX;
    var startY = e.offsetY;
    var currentX = startX;
    var currentY = startY;

    function onMouseMove (e) {
      var nextX = e.offsetX;
      var nextY = e.offsetY;
      pan(rotation, nextX - currentX, nextY - currentY, canvas);
      currentX = nextX;
      currentY = nextY;
    }

    function onMouseUp () {
      canvas.removeEventListener('mousemove', onMouseMove);
      canvas.removeEventListener('mouseup', onMouseUp);
    }
  };
  return onMouseDown;
}

function pan (rotation, dX, dY, canvas) {
  rotation[0] += 10 * dY / canvas.clientHeight;
  rotation[1] += 10 * dX / canvas.clientWidth;
  return rotation;
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
