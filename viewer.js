var Component = require('./component.js');
var Renderer = require('./renderer.js');

var viewerHtml = require('./viewer.html');

var Viewer = Component('kokorigami-viewer', HTMLElement, {
  _: {
    writable: false,
    enumerable: false,
    value: {
      shadow: null,
      renderer: null
    }
  },
  data: 'data',
  createdCallback: function () {
    this._.shadow = this.createShadowRoot();
    this._.shadow.innerHTML = viewerHtml;

    var canvas = this._.shadow.querySelector('canvas');
    this._.renderer = new Renderer(canvas, this.data);
    return;
  },
  attributeChangedCallback: function () {
    this._.renderer.data(this.data);
    this._.renderer.play();
    return;
  }
});

module.exports = Viewer;
