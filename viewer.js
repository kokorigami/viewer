var Component = require('./component.js');
var Renderer = require('./renderer.js');

var Viewer = Component('kokorigami-viewer', HTMLElement, {
  _: {
    writable: false,
    enumerable: false,
    value: {
      shadow: null,
      canvas: null,
      renderer: null
    }
  },
  data: 'data',
  createdCallback: function () {
    this._.shadow = this.createShadowRoot();
    this._.canvas = document.createElement('canvas');
    this._.shadow.appendChild(this._.canvas);
    this._.renderer = new Renderer(this._.canvas, this.data);
    return;
  },
  attributeChangedCallback: function () {
    this._.renderer.data(this.data);
    this._.renderer.play();
    return;
  }
});

module.exports = Viewer;
