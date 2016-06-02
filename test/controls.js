var chai = require('chai');
var spies = require('chai-spies');
var expect = chai.expect;
chai.use(spies);

var controls = require('../src/controls.js');
var createCamera = require('3d-view');

describe('controls', function () {
  var el, rotateSpy, panSpy;
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
    panSpy = chai.spy.on(camera, 'pan');
    document.body.appendChild(el);
  });

  afterEach(function () {
    document.body.removeChild(el);
  });

  it('attaches the mousedown controller to an element given a camera', function () {
    var ctrl = controls.attach(el, camera);
    dispatchEvents(el, ['mousedown','mousemove','mouseup','mousemove']);
    expect(ctrl.mousedown).to.be.a('function');
    expect(rotateSpy).to.have.been.called.once;
  });

  it('attaches the wheel controller to an element given a camera', function () {
    var ctrl = controls.attach(el, camera);
    dispatchEvents(el, ['wheel']);
    expect(ctrl.wheel).to.be.a('function');
    expect(panSpy).to.have.been.called.once;
  });

  it('stops rotating the camera when the mouse leaves the element', function () {
    var ctrl = controls.attach(el, camera);
    dispatchEvents(el, ['mousedown','mousemove','mouseleave','mousemove']);
    expect(ctrl.mousedown).to.be.a('function');
    expect(rotateSpy).to.have.been.called.once;
  });

  it('removes the mousedown controller from an element given the controls', function () {
    var ctrl = controls.attach(el, camera);
    controls.remove(el, ctrl);
    dispatchEvents(el, ['mousedown','mousemove','mouseup','mousemove']);
    expect(rotateSpy).to.not.have.been.called();
  });

  it('removes the wheel controller from an element given the controls', function () {
    var ctrl = controls.attach(el, camera);
    controls.remove(el, ctrl);
    dispatchEvents(el, ['wheel']);
    expect(panSpy).to.not.have.been.called();
  });
});

function dispatchEvents(el, events) {
  var mouse = {
    mousedown:  new MouseEvent('mousedown' , {view: window}),
    mousemove:  new MouseEvent('mousemove' , {view: window}),
    mouseup:    new MouseEvent('mouseup'   , {view: window}),
    mouseleave: new MouseEvent('mouseleave', {view: window}),
    wheel:      new WheelEvent('wheel', {deltaY: 10})
  };

  events.forEach(function (event) {
    el.dispatchEvent(mouse[event]);
  });
}
