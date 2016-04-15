var Renderer = require('./renderer.js');
var viewerHTML = require('./viewer.html');

var Viewer = function (el) {
  this.el = el || document.createElement('div');
  this.el.innerHTML = viewerHTML;

  var canvas = this.el.querySelector('canvas');
  this.renderer = new Renderer(canvas);
  return this;
};

Viewer.prototype = {};

Object.defineProperty(Viewer.prototype, 'data', {
  enumerable: true,
  set: function (data) {
    this.renderer.data(data);
    this.renderer.play();
    return this._data = data;
  },

  get: function () {
    return this._data;
  }
});

module.exports = Viewer;
