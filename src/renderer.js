var createContext = require('gl-context');
var createShader = require('gl-shader');
var createGeometry = require('gl-geometry');
var createCamera = require('3d-view');
var createTexture = require('gl-texture2d');
var getPixels = require('get-pixels');

var glslify = require('glslify');
var faceFs = glslify('./face-fs.glsl');
var faceVs = glslify('./face-vs.glsl');
var foldFs = glslify('./fold-fs.glsl');
var foldVs = glslify('./fold-vs.glsl');
var depthnormalFs = glslify('./depthnormal-fs.glsl');

var m4 = require('gl-mat4');

var Renderer = function (canvas) {
  this.gl = null;
  this.canvas = null;
  this.raf = null;
  this.textures = [null, null];
  this.glTextures = [];
  this.render = this.render.bind(this);

  this.camera = createCamera({
    center: [0.5, 0.5, 0],
    eye: [0, 1, -3],
    distanceLimits: [1, 6],
    up: [0, 0, 1],
    mode: 'orbit'
  });

  this.uniforms = {
    u_color: [1, 1, 1, 1],
    u_lightWorldPos: [0, 0, 6],
    u_lightColor: [1, 0.9, 0.8, 1],
    u_ambient: [0.2, 0.2, 0.3, 1],
    u_shininess: 10,
    u_near: 0,
    u_far: 10,
    u_view: m4.create(),
    u_camera: m4.create(),
    u_world: m4.create(),
    u_worldView: m4.create(),
    u_worldViewProjection: m4.create(),
    u_texture: 0,
    u_textureBack: 1
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
    this.texturize(this.textures);
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
  this.update(true);
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

  if (this.update() || this.resize() || this.snap()) {
    gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    var projection = m4.create();
    m4.perspective(
      projection,
      Math.PI/4.0,
      gl.drawingBufferWidth/gl.drawingBufferHeight,
      this.uniforms.u_near,
      this.uniforms.u_far
    );

    m4.invert(this.uniforms.u_camera, this.uniforms.u_view);
    m4.copy(this.uniforms.u_world, this.camera.computedMatrix);
    m4.multiply(this.uniforms.u_worldView, this.uniforms.u_view, this.uniforms.u_world);
    m4.multiply(this.uniforms.u_worldViewProjection, projection, this.uniforms.u_worldView);

    this.uniforms.u_texture = this.glTextures[0].bind(0);
    this.uniforms.u_textureBack = this.glTextures[1].bind(1);

    renderPass(gl, this.shaders.face, this.buffers.face, this.uniforms, 'TRIANGLES');
    this.update(false);
  }

  this.raf = requestAnimationFrame(this.render);
};

Renderer.prototype.update = function (force) {
  if (force !== undefined) this._update = force;
  return this._update;
};

Renderer.prototype.snap = function () {
  var t = Date.now();
  this.camera.idle(t - 20);
  this.camera.flush(t - 100);
  this.camera.recalcMatrix(t - 25);

  return !equivalent(this.camera.computedMatrix, this.uniforms.u_world);
};

Renderer.prototype.resize = function () {
  var resolution = window.devicePixelRatio || 1;
  resolution = Math.max(1, resolution);
  var canvas = this.canvas;
  var width  = canvas.offsetWidth  * resolution || 0;
  var height = canvas.offsetHeight * resolution || 0;
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
    return true;
  }
  return false;
};

Renderer.prototype.texturize = function (sources) {
  var updateTexture = function (name, resource, i) {
    if (this.glTextures[i]) {
      this.glTextures[i].dispose();
    }
    this.textures[i] = name;
    this.glTextures[i] = createTexture(this.gl, resource);
    this.update(true);
  }.bind(this);

  var promises = sources.map(function (source, i) {
    if (!source) {
      updateTexture(null, createPixel([255, 255, 255, 255]), i);
      return Promise.resolve();
    }
    return getImage(source)
      .then(function (pixels) { updateTexture(source, pixels, i); })
      .catch(function (err) { console.log(err); });
  });

  return Promise.all(promises);
};

function createPixel(rgba) {
  var color = new Uint8ClampedArray(rgba);
  var pixel = new ImageData(color, 1, 1);
  return pixel;
}

function getImage(source) {
  return new Promise(function (res, rej) {
    getPixels(source, function (err, pixels) {
      if (err) rej(err);
      else res(pixels);
    });
  });
}

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
