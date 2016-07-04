var createContext = require('gl-context');
var createShader = require('gl-shader');
var createGeometry = require('gl-geometry');
var createCamera = require('3d-view');
var createTexture = require('gl-texture2d');

var glslify = require('glslify');
var faceFs = glslify('./face-fs.glsl');
var faceVs = glslify('./face-vs.glsl');
var foldFs = glslify('./fold-fs.glsl');
var foldVs = glslify('./fold-vs.glsl');
var depthnormalFs = glslify('./depthnormal-fs.glsl');

var getImage = require('./get-image.js');
var white = new Uint8ClampedArray([255, 255, 255, 255]);
var block = new ImageData(white, 1, 1);

var m4 = require('gl-mat4');

var Renderer = function (canvas) {
  this.gl = null;
  this.canvas = null;
  this.raf = null;
  this.texture = null;
  this.glTexture = null;
  this.render = this.render.bind(this);
  this.camera = createCamera({
    center: [0.5, 0.5, 0],
    eye: [0, 1, -3],
    distanceLimits: [1, 6],
    up: [0, 0, 1],
    mode: 'orbit'
  });
  this.uniforms = {
    u_lightWorldPos: [0, 0, 6],
    u_lightColor: [1, 0.9, 0.8, 1],
    u_ambient: [0, 0, 0, 1],
    u_shininess: 70,
    u_near: 0,
    u_far: 10,
    u_view: m4.identity([])
  };
  this.shaders = {};
  this.buffers = {};
  this.initialize(canvas);
  return this;
};

module.exports = Renderer;

Renderer.prototype.initialize = function (canvas) {
  var gl;
  try {
    gl = createContext(canvas);
    this.gl = gl;
    this.canvas = canvas;
    // gl.enable(gl.DEPTH_TEST);
    this.shaders.face = createShader(gl, faceVs, faceFs);
    this.shaders.depth = createShader(gl, faceVs, depthnormalFs);
    this.shaders.fold = createShader(gl, foldVs, foldFs);
    this.texturize(this.texture);
  } catch (e) {
    // No WebGL context
    console.log(e.message);
  }
};

Renderer.prototype.data = function (data) {
  if (!this.gl) return;

  this.buffers.face = createGeometry(this.gl)
    .attr('position', data.position)
    .attr('normal', data.normal)
    .attr('texcoord', data.texcoord);

  /*
    for folds, we'll need:
      foldType (1 float per vertex)
      lengthSoFar (1 float per vertex)
      position
  */
  this.updated = true;
  return this;
};

Renderer.prototype.play = function () {
  this.stop();

  if (!this.gl) return;
  this.render();
  return this;
};

Renderer.prototype.stop = function () {
  cancelAnimationFrame(this.raf);
  return this;
};

Renderer.prototype.render = function () {
  var gl = this.gl;
  this.update();

  if (this.updated) {
    gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    renderPass(gl, this.shaders.face, this.buffers.face, this.uniforms, 'TRIANGLES');
    this.updated = false;
  }

  this.raf = requestAnimationFrame(this.render);
};

Renderer.prototype.update = function (t) {
  t = t || Date.now();

  var gl = this.gl;
  var uniforms = this.uniforms;
  var camera = this.camera;

  camera.idle(t - 20);
  camera.flush(t - 100);
  camera.recalcMatrix(t - 25);

  var updated = this.resize() || !equivalent(camera.computedMatrix, uniforms.u_world);

  if (updated) {
    var projection = m4.perspective(
        [],
        Math.PI/4.0,
        gl.drawingBufferWidth/gl.drawingBufferHeight,
        uniforms.u_near,
        uniforms.u_far
    );

    uniforms.u_camera = m4.invert([], uniforms.u_view);
    uniforms.u_world = m4.copy([], camera.computedMatrix);
    uniforms.u_worldView = m4.multiply([], uniforms.u_view, uniforms.u_world);
    uniforms.u_worldViewProjection = m4.multiply([], projection, uniforms.u_worldView);

    this.updated = true;
  }
  return this;
};

Renderer.prototype.resize = function (resolution) {
  resolution = resolution || 1;
  resolution = Math.max(1, resolution);
  var canvas = this.canvas;
  var width  = canvas.clientWidth  * resolution | 0;
  var height = canvas.clientHeight * resolution | 0;
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
    return true;
  }
  return false;
};

Renderer.prototype.texturize = function (image) {
  var updateTexture = function (name, resource) {
    if (this.glTexture) this.glTexture.dispose();
    this.texture = name;
    this.glTexture = createTexture(this.gl, resource);
    this.uniforms.u_texture = this.glTexture.bind();
    this.updated = true;
  }.bind(this);

  // TODO: extend this block to handle resources covered by gl-texture2d without get-pixels
  if (!image) {
    updateTexture(null, block);
    return Promise.resolve();
  }

  return getImage(image)
    .then(function (pixels) { updateTexture(image, pixels); })
    .catch(function (err) { console.log(err); });
};

function renderPass(gl, shader, geom, uniforms, drawType) {
  shader.bind();
  shader.uniforms = uniforms;
  geom.bind(shader);
  geom.draw(gl[drawType]);
  geom.unbind();
}

var EPSILON = 0.000001;
function equivalent(a, b, delta) {
  delta = delta || EPSILON;
  if (a == null && b == null) return true;
  if ((a instanceof Float64Array || a instanceof Float32Array
    || a instanceof Array) && (b instanceof Float64Array
    || b instanceof Float32Array || b instanceof Array)) {

    if (a.length != b.length) return false;

    var check = true;
    for (var i = 0; i < a.length; i++) {
      check = check && equivalent(a[i], b[i], delta);
    }
    return check;
  }
  return Math.abs(a - b) <= delta;
}
