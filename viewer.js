var Renderer = require('./renderer.js');
var Model = require('./model.js');
var viewerHTML = require('./viewer.html');

var Viewer = function (el) {
  this.el = el || document.createElement('div');
  this.el.innerHTML = viewerHTML;

  var canvas = this.el.querySelector('canvas');

  this.model = new Model();
  this.renderer = new Renderer(canvas);
  return this;
};

Viewer.prototype = {};

Object.defineProperty(Viewer.prototype, 'data', {
  enumerable: true,
  set: function (data) {
    this.model.set(data);
    this.renderer.data(this.model.frameGeometry(0));
    this.renderer.play();
    return this.model;
  },

  get: function () {
    return this.model;
  }
});

module.exports = Viewer;
