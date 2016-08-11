var Renderer = require('./renderer.js');
var Model = require('./model.js');
var controls = require('./controls.js');
var EventEmitter = require('events');

var Viewer = function (el) {
  this.el = el || document.createElement('canvas');
  this.renderer = new Renderer(this.el);
  this._emitter = new EventEmitter();
  this._model = new Model();
  this._frame = 0;
  this.raf = null;
  this.controls = null;

  this.el.classList.add('kokorigami-viewer');
  return this;
};

Viewer.style = require('./style.js');
Viewer.ActionBar = require('./actionbar.js');
Viewer.prototype = {};

Viewer.prototype.render = function () {
  this.controls = controls.attach(this.el, this.renderer.camera);
  this.renderer.play();
};

Viewer.prototype.play = function (frame, fps) {
  fps = fps || this.model.fps;

  this.stop();
  if (typeof frame === 'undefined') {
    frame = this.model.frames.length;
  }

  var spf = 1000 / fps;
  var startFrame = this.frame;
  var animationTime = Math.abs(frame - startFrame) * spf;
  var start;

  var render = function (timestamp) {
    start = start || timestamp;

    var delta = (timestamp - start) / animationTime;
    delta = Math.min(delta, 1);
    delta = Math.max(delta, 0);
    this.frame = (frame - startFrame) * delta + startFrame;

    if (this.frame !== frame) {
      this.raf = requestAnimationFrame(render);
    } else {
      this._emitter.emit('stop');
    }
  }.bind(this);
  this.raf = requestAnimationFrame(render);
  this._emitter.emit('play');
};

Viewer.prototype.playStep = function (step, fps) {
  var frames = this.model.stepFrames(step);
  frames[0] = Math.max(frames[0], 0);
  frames[1] = Math.min(frames[1], this.model.frames.length);
  if (this.frame < frames[0] || this.frame > frames[0]) {
    this.frame = frames[0];
  }
  this.play(frames[1] + 1, fps);
};

Viewer.prototype.stop = function () {
  if (!this.raf) return;
  cancelAnimationFrame(this.raf);
  this.raf = null;
  this._emitter.emit('stop');
};

Viewer.prototype.destroy = function () {
  this.stop();
  this.renderer.stop();
  this._emitter.removeAllListeners();
  controls.remove(this.el, this.controls);
};

Viewer.prototype.on = function (event, callback) {
  this._emitter.on(event, callback);
};

Viewer.prototype.once = function (event, callback) {
  this._emitter.once(event, callback);
};

Viewer.prototype.off = function (event, callback) {
  this._emitter.removeListener(event, callback);
};

Object.defineProperty(Viewer.prototype, 'model', {
  enumerable: true,
  set: function (data) {
    this._model.set(data);
    this.frame = 0; // reset frame
    this._emitter.emit('swap', data);
    return this.model;
  },

  get: function () {
    return this._model;
  }
});

Object.defineProperty(Viewer.prototype, 'frame', {
  enumerable: true,
  set: function (frame) {
    if (typeof frame !== 'number' || isNaN(frame)) return this.frame;
    frame = Math.max(frame, 0);
    frame = Math.min(frame, this.model.frames.length);
    var renderFrame = Math.min(frame, this.model.final);
    this.renderer.data(this.model.frameGeometry(renderFrame));
    this._frame = frame;
    this._emitter.emit('update', frame);
    return frame;
  },

  get: function () {
    return this._frame;
  }
});

Object.defineProperty(Viewer.prototype, 'step', {
  enumerable: true,
  get: function () {
    return this.model.stepFromFrame(this.frame);
  }
});

Object.defineProperty(Viewer.prototype, 'textures', {
  enumerable: true,
  set: function (images) {
    if (!images instanceof Array || images.length !== 2) {
      images = [images, images];
    }
    return this.renderer.texturize(images);
  },
  get: function () {
    return this.renderer.textures;
  }
});

Object.defineProperty(Viewer.prototype, 'background', {
  enumerable: true,
  set: function (rgba) {
    var css = [rgba[0], rgba[1], rgba[2], rgba[3] / 255];
    this.el.style.backgroundColor = 'rgba(' + css.join(',') + ')';
    this.renderer.ambient(rgba);
    return rgba;
  },
  get: function () {
    return this.renderer.ambient();
  }
});

module.exports = Viewer;
