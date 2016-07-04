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
  this.uniforms = {};
  this.shaders = {};
  this.buffers = {};
  this.loading = Promise.reject('Renderer has not been initialized.');
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
    this.loading = this.texturize(this.texture);
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
  return this;
};

Renderer.prototype.play = function () {
  this.stop();

  if (!this.gl) return;
  this.loading
    .then(this.render)
    .catch(function (err) { console.log(err); });
  return this;
};

Renderer.prototype.stop = function () {
  cancelAnimationFrame(this.raf);
  return this;
};

Renderer.prototype.render = function () {
  this.resize();

  var gl = this.gl;
  gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  var uniforms = updateUniforms(gl, this.camera, this.glTexture, this.uniforms);
  renderPass(gl, this.shaders.face, this.buffers.face, uniforms, 'TRIANGLES');

  this.raf = requestAnimationFrame(this.render);
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
  return getImage(image)
    .then(function (pixels) {
      this.texture = image;
      this.glTexture = createTexture(this.gl, pixels);
    }.bind(this))
    .catch(function (err) { console.log(err); });
};

function updateUniforms(gl, camera, texture, uniforms) {
  // Updates uniforms in place.
  if (!Object.keys(uniforms).length) {
    uniforms.u_lightWorldPos = [0, 0, 6];
    uniforms.u_lightColor = [1, 0.9, 0.8, 1];
    uniforms.u_ambient = [0, 0, 0, 1];
    uniforms.u_shininess = 70;
  }

  var t = Date.now();
  camera.idle(t - 20);
  camera.flush(t - 100);
  camera.recalcMatrix(t - 25);

  // Update camera uniforms.
  var projection = m4.perspective(
      [],
      Math.PI/4.0,
      gl.drawingBufferWidth/gl.drawingBufferHeight,
      0,
      10
  );

  uniforms.u_near = 0;
  uniforms.u_far = 10;
  uniforms.u_view = m4.identity([]);
  uniforms.u_camera = m4.invert([], uniforms.u_view);
  uniforms.u_world = camera.computedMatrix;
  uniforms.u_worldView = m4.multiply([], uniforms.u_view, uniforms.u_world);
  uniforms.u_worldViewProjection = m4.multiply([], projection, uniforms.u_worldView);
  uniforms.u_texture = texture.bind();

  return uniforms;
}

function renderPass(gl, shader, geom, uniforms, drawType) {
  shader.bind();
  shader.uniforms = uniforms;
  geom.bind(shader);
  geom.draw(gl[drawType]);
  geom.unbind();
}
