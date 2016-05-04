var Renderer = require('./renderer.js');
var Model = require('./model.js');
var viewerHTML = require('./viewer.html');

var Viewer = function (el) {
  this.el = el || document.createElement('div');
  this.el.innerHTML = viewerHTML;

  var canvas = this.el.querySelector('canvas');

  this._ = {
    model: new Model(),
    progress: 0,
    renderer: new Renderer(canvas)
  };

  return this;
};

Viewer.prototype = {};

Object.defineProperty(Viewer.prototype, 'model', {
  enumerable: true,
  set: function (data) {
    this._.model.set(data);
    this.progress = 0; // reset progress
    return this.model;
  },

  get: function () {
    return this._.model;
  }
});

Object.defineProperty(Viewer.prototype, 'progress', {
  enumerable: true,
  set: function (progress) {
    progress = Math.max(progress, 0);
    progress = Math.min(progress, this.model.frames.length);
    this._.progress = progress;
    this._.renderer.data(this.model.frameGeometry(progress));
    return this.progress;
  },

  get: function () {
    return this._.progress;
  }
});

Viewer.prototype.play = function () {
  this._.renderer.play();
};

Viewer.prototype.next = function () {
  this.progress++;
};

Viewer.prototype.prev = function () {
  this.progress--;
};

module.exports = Viewer;
