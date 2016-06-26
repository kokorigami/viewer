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

Viewer.prototype.play = function (frame) {
  this.stop();
  if (typeof frame === 'undefined') {
    frame = this.model.frames.length;
  }

  var spf = 1000 / this.model.fps;
  var startFrame = this.frame;
  var animationTime = Math.abs(frame - startFrame) * spf;
  var start;

  var render = function (timestamp) {
    if (!start) start = timestamp;

    var delta = (timestamp - start) / animationTime;
    delta = Math.min(delta, 1);
    delta = Math.max(delta, 0);
    this.frame = (frame - startFrame) * delta + startFrame;

    if (this.frame !== frame) {
      this.raf = requestAnimationFrame(render);
    }
  }.bind(this);
  this.raf = requestAnimationFrame(render);
};

Viewer.prototype.stop = function () {
  cancelAnimationFrame(this.raf);
};

Viewer.prototype.next = function () {
  var nextFrames = this.model.stepFrames(this.step + 1);
  var toFrame = Math.min(nextFrames[0], this.model.frames.length);
  this.play(toFrame);
};

Viewer.prototype.prev = function () {
  var prevFrames = this.model.stepFrames(this.step - 1);
  var toFrame = Math.max(prevFrames[0], 0);
  this.play(toFrame);
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
    this._frame = frame;

    if (frame <= this.model.final) {
      this.renderer.data(this.model.frameGeometry(frame));
    }

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

module.exports = Viewer;
