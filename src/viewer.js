var Renderer = require('./renderer.js');
var Model = require('./model.js');

var Viewer = function (el) {
  this.el = el || document.createElement('canvas');
  this.el.classList.add('kokorigami-viewer');
  this._ = {
    model: new Model(),
    frame: 0,
    renderer: new Renderer(this.el)
  };
  return this;
};

Viewer.prototype = {};

Viewer.prototype.render = function () {
  this._.renderer.play();
};

Viewer.prototype.play = function (frame) {
  if (this.raf) cancelAnimationFrame(this.raf);
  if (typeof frame === 'undefined') {
    frame = this.model.lastFrame;
  }

  var ANIMATION_TIME = 800;
  var startFrame = this.frame;
  var start;

  var render = function (timestamp) {
    if (!start) start = timestamp;
    var currentFrame = this.frame;

    if (currentFrame !== frame) {
      var delta = (timestamp - start) / ANIMATION_TIME;
      delta = Math.min(delta, 1);
      delta = Math.max(delta, 0);
      this.frame = (frame - startFrame) * delta + startFrame;
      this.raf = requestAnimationFrame(render);
    }
  }.bind(this);
  this.raf = requestAnimationFrame(render);
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

Object.defineProperty(Viewer.prototype, 'model', {
  enumerable: true,
  set: function (data) {
    this._.model.set(data);
    this.frame = 0; // reset frame
    return this.model;
  },

  get: function () {
    return this._.model;
  }
});

Object.defineProperty(Viewer.prototype, 'frame', {
  enumerable: true,
  set: function (frame) {
    frame = Math.max(frame, 0);
    frame = Math.min(frame, this.model.lastFrame);
    this._.frame = frame;
    this._.renderer.data(this.model.frameGeometry(frame));
    return frame;
  },

  get: function () {
    return this._.frame;
  }
});

Object.defineProperty(Viewer.prototype, 'step', {
  enumerable: true,
  get: function () {
    return this.model.stepFromFrame(this.frame);
  }
});

module.exports = Viewer;
