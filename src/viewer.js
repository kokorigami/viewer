var Renderer = require('./renderer.js');
var Model = require('./model.js');

/*
TODO: think about forward/reverse interactions rather than jumping
to previous or next steps.
*/

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
  set: function (step) {
    step = Math.max(step, 0);
    step = Math.min(step, this.model.steps);
    if (step !== this.step) {
      this.frame = this.model.stepFrames(step)[0];
    }
    return step;
  },

  get: function () {
    return this.model.stepFromFrame(this.frame);
  }
});

Viewer.prototype.render = function () {
  this._.renderer.play();
};

Viewer.prototype.play = function () {
  this.playTo(this.model.lastFrame);
};

Viewer.prototype.playStep = function (step) {
  if (typeof step === 'number' && !isNaN(step)) {
    this.step = step;
  }
  var frameRange = this.model.stepFrames(this.step);
  this.playTo(frameRange[1]);
};

Viewer.prototype.playTo = function (frame) {
  var spf = 1000 / this.model.fps;
  var start, raf;

  var render = function (timestamp) {
    var currentFrame = this.frame;
    if (!start) start = timestamp;

    if (currentFrame < frame) {
      if (timestamp - start > spf) this.frame++;
      raf = requestAnimationFrame(render);
    } else {
      cancelAnimationFrame(raf);
    }
  }.bind(this);
  raf = requestAnimationFrame(render);
};

Viewer.prototype.next = function () {
  this.step++;
};

Viewer.prototype.prev = function () {
  this.step--;
};

module.exports = Viewer;
