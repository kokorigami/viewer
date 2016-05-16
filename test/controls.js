var chai = require('chai');
var spies = require('chai-spies');
var expect = chai.expect;
chai.use(spies);

var controls = require('../src/controls.js');
var createCamera = require('3d-view');

describe('controls', function () {
  var el, rotateSpy;
  var camera = createCamera({
    center: [0, 0.5, 0],
    eye: [0, 1, -3],
    distanceLimits: [1, 100],
    up: [0, 0, 1],
    mode: 'orbit'
  });

  beforeEach(function () {
    el = document.createElement('canvas');
    rotateSpy = chai.spy.on(camera, 'rotate');
    document.body.appendChild(el);
  });

  afterEach(function () {
    document.body.removeChild(el);
  });

  it('attaches the mousedown controller to an element given a camera', function () {
    var ctrl = controls.attach(el, camera);
    dispatchMouseEvents(el);
    expect(ctrl.mousedown).to.be.a('function');
    expect(rotateSpy).to.have.been.called.once;
  });

  it('removes the mousedown controller from an element given the controls', function () {
    var ctrl = controls.attach(el, camera);
    controls.remove(el, ctrl);
    dispatchMouseEvents(el);
    expect(rotateSpy).to.not.have.been.called();
  });
});

function dispatchMouseEvents(el) {
  var mousedown = new MouseEvent('mousedown', {view: window});
  var mousemove = new MouseEvent('mousemove', {view: window});
  var mouseup   = new MouseEvent('mouseup'  , {view: window});

  el.dispatchEvent(mousedown);
  el.dispatchEvent(mousemove);
  el.dispatchEvent(mouseup);
  el.dispatchEvent(mousemove);
}
