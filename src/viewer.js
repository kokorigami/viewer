var Renderer = require('./renderer.js');
var Model = require('./model.js');
var controls = require('./controls.js');

var Viewer = function (el) {
  this.el = el || document.createElement('canvas');
  this.renderer = new Renderer(this.el);
  this._model = new Model();
  this._frame = 0;
  this.raf = null;

  this.el.classList.add('kokorigami-viewer');
  return this;
};

Viewer.prototype = {};

Viewer.prototype.render = function () {
  controls.set(this.el, this.renderer);
  this.renderer.play();
};

Viewer.prototype.play = function (frame) {
  this.stop();
  if (typeof frame === 'undefined') {
    frame = this.model.lastFrame;
  }

  var spf = 1000 / this.model.fps;
  var start;
  var delta = frame < this.frame ? -1 : 1;

  var render = function (timestamp) {
    var currentFrame = this.frame;
    if (!start) start = timestamp;

    if (currentFrame !== frame) {
      if (timestamp - start > spf) this.frame += delta;
      this.raf = requestAnimationFrame(render);
    }
  }.bind(this);
  this.raf = requestAnimationFrame(render);
};

Viewer.prototype.stop = function () {
  cancelAnimationFrame(this.raf);
};

Viewer.prototype.next = function () {
  var stepFrames = this.model.stepFrames(this.step);
  var nextFrames = this.model.stepFrames(this.step + 1);
  nextFrames = this.frame < stepFrames[1] ? stepFrames : nextFrames;
  this.play(nextFrames[1]);
};

Viewer.prototype.prev = function () {
  var stepFrames = this.model.stepFrames(this.step);
  var prevFrames = this.model.stepFrames(this.step - 1);
  prevFrames = this.frame > stepFrames[0] ? stepFrames : prevFrames;
  this.play(prevFrames[0]);
};

Viewer.prototype.destroy = function () {
  this.stop();
  this.renderer.stop();
  controls.remove(this.el);
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
    frame = Math.max(frame, 0);
    frame = Math.min(frame, this.model.lastFrame);
    this._frame = frame;
    this.renderer.data(this.model.frameGeometry(frame));
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
