{
  var createContext = require('gl-context');
  var createShader = require('gl-shader');
  var createGeometry = require('gl-geometry');
  var createCamera = require('3d-view');
  var createTexture = require('gl-texture2d');
  var createFBO = require('gl-fbo');
  var getPixels = require('get-pixels');
  var getQuad = require('gl-big-quad');

  var glslify = require('glslify');
  var faceFs = glslify('./face-fs.glsl');
  var faceVs = glslify('./face-vs.glsl');
  var foldFs = glslify('./fold-fs.glsl');
  var foldVs = glslify('./fold-vs.glsl');
  var drawFs = glslify('./draw-fs.glsl');
  var drawVs = glslify('./draw-vs.glsl');
  var ssaoFs = glslify('./ssao-fs.glsl');
  var depthnormalFs = glslify('./depthnormal-fs.glsl');

  var m4 = require('gl-mat4');
}

var Renderer = function (canvas) {
  this.gl = null;
  this.canvas = null;
  this.raf = null;
  this.textures = [null, null];
  this.glTextures = [null, null];
  this.render = this.render.bind(this);

  this.camera = createCamera({
    center: [0.5, 0.5, 0],
    eye: [0, 1, -3],
    distanceLimits: [2, 6],
    up: [0, 0, 1],
    mode: 'orbit'
  });

  this.uniforms = {
    u_lightWorldPos: [0, 0, 6],
    u_lightColor: [1, 0.9, 0.8, 1],
    u_ambient: [0.2, 0.2, 0.3, 1],
    u_shininess: 10,
    u_near: 1,
    u_far: 10,
    u_view: m4.create(),
    u_camera: m4.create(),
    u_world: m4.create(),
    u_worldView: m4.create(),
    u_worldViewProjection: m4.create(),
    u_origami: 0,
    u_sampler: 1,
    u_texture0: 2,
    u_texture1: 3
  };

  this.shaders = {};
  this.buffers = {};
  this.fboIndex = 0;
  this.framebuffers = [null, null];
  this.initialize(canvas);
  return this;
};

module.exports = Renderer;

Renderer.prototype.initialize = function (canvas) {
  var gl;
  try {
    gl = createContext(canvas);
    gl.enable(gl.DEPTH_TEST);
    this.gl = gl;
    this.canvas = canvas;
    this.shaders.face = createShader(gl, faceVs, faceFs);
    this.shaders.depth = createShader(gl, faceVs, depthnormalFs);
    this.shaders.fold = createShader(gl, foldVs, foldFs);
    this.shaders.draw = createShader(gl, drawVs, drawFs);
    this.shaders.ssao = createShader(gl, drawVs, ssaoFs);
    this.buffers.quad = getQuad(gl);
    this.framebuffers[0] = createFBO(gl, [canvas.width, canvas.height]);
    this.framebuffers[1] = createFBO(gl, [canvas.width, canvas.height]);
    this.texturize(this.textures);
    this.resize();
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
  var shaders = this.shaders;
  var buffers = this.buffers;
  var uniforms = this.uniforms;

  if (this.update() || this.resize() || this.snap()) {
    uniforms.u_sampler = this.renderSwap(function () {
      renderPass(gl, shaders.depth, buffers.face, uniforms, 'TRIANGLES');
    });

    uniforms.u_sampler = this.renderSwap(function () {
      renderPass(gl, shaders.ssao, buffers.quad, uniforms, 'TRIANGLES');
    });

    uniforms.u_origami = this.renderSwap(function () {
      renderPass(gl, shaders.face, buffers.face, uniforms, 'TRIANGLES');
    });

    this.renderBuffer(null, function () {
      renderPass(gl, shaders.draw, buffers.quad, uniforms, 'TRIANGLES');
    });

    this.update(false);
  }

  this.raf = requestAnimationFrame(this.render);
};

Renderer.prototype.renderSwap = function (render) {
  var index = this.fboIndex;
  this.fboIndex = (index + 1) % 2;
  return this.renderBuffer(index, render);
};

Renderer.prototype.renderBuffer = function (fboIndex, render) {
  var gl = this.gl;
  var isFBO = fboIndex != null;
  var fbo = null;

  if (!isFBO) gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  else {
    fbo = this.framebuffers[fboIndex];
    fbo.bind();
  }

  gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  render();

  if (isFBO) return fbo.color[0].bind(fboIndex);
};

Renderer.prototype.update = function (force) {
  if (force !== undefined) this._update = force;
  return this._update;
};

Renderer.prototype.snap = function () {
  var t = Date.now();
  var gl = this.gl;
  var uniforms = this.uniforms;
  var camera = this.camera;

  camera.idle(t - 20);
  camera.flush(t - 100);
  camera.recalcMatrix(t - 25);

  var willUpdate = !equivalent(camera.computedMatrix, uniforms.u_world);

  if (willUpdate) {
    var projection = m4.create();
    m4.perspective(
      projection,
      Math.PI/4.0,
      gl.drawingBufferWidth/gl.drawingBufferHeight,
      uniforms.u_near,
      uniforms.u_far
    );

    m4.invert(uniforms.u_camera, uniforms.u_view);
    m4.copy(uniforms.u_world, camera.computedMatrix);
    m4.multiply(uniforms.u_worldView, uniforms.u_view, uniforms.u_world);
    m4.multiply(uniforms.u_worldViewProjection, projection, uniforms.u_worldView);
  }

  return willUpdate;
};

Renderer.prototype.resize = function () {
  var resolution = window.devicePixelRatio || 1;
  resolution = Math.max(1, resolution);
  var canvas = this.canvas;
  var width  = canvas.offsetWidth  * resolution || 1;
  var height = canvas.offsetHeight * resolution || 1;
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
    this.framebuffers[0].shape = [width, height];
    this.framebuffers[1].shape = [width, height];
    return true;
  }
  return false;
};

Renderer.prototype.texturize = function (sources) {
  var updateTexture = function (name, resource, i) {
    var offset = this.framebuffers.length;
    var uniform = 'u_texture' + i;

    if (this.glTextures[i]) {
      this.glTextures[i].dispose();
    }
    this.textures[i] = name;
    this.glTextures[i] = createTexture(this.gl, resource);

    this.uniforms[uniform] = this.glTextures[i].bind(offset + i);
    this.update(true);
  }.bind(this);

  var promises = sources.map(function (source, i) {
    if (!source || source instanceof Array) {
      source = source || [255, 255, 255, 255];
      updateTexture(source, createPixel(source), i);
      return Promise.resolve();
    }
    return getImage(source)
      .then(function (pixels) { updateTexture(source, pixels, i); })
      .catch(function (err) { console.log(err); });
  });

  return Promise.all(promises);
};

Renderer.prototype.ambient = function (rgba) {
  if (rgba !== undefined) {
    var normalized = rgba.map(function (c) { return c / 255; });
    this.uniforms.u_ambient = normalized;
    return rgba;
  }

  var ambient = this.uniforms.u_ambient;
  var unnormalized = ambient.map(function (c) { return Math.round(c * 255); });
  return unnormalized;
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
