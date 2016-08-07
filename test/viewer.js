
var chai = require('chai');
var spies = require('chai-spies');
var expect = chai.expect;
chai.use(spies);
var Viewer = require('../src/viewer.js');
var Model = require('../src/model.js');
var data = require('./data.json');

describe('viewer', function () {
  var model = new Model(data);

  it('can be created', function () {
    var viewer = new Viewer();
    expect(viewer).to.exist;
  });

  it('can set the background color', function () {
    var viewer = new Viewer();
    viewer.background = [0, 255, 255, 102];
    expect(viewer.el.style.backgroundColor).to.equal('rgba(0, 255, 255, 0.4)');
  });

  it('can update frames', function () {
    var viewer = new Viewer();
    viewer.model = model;
    viewer.frame = 1;
    expect(viewer.frame).to.equal(1);
  });

  it('emits update when the frame is updated', function () {
    var callback = chai.spy();
    var viewer = new Viewer();
    viewer.model = model;
    viewer.on('update', callback);
    viewer.frame = 1;
    expect(callback).to.have.been.called();
  });

  it('can determine the current step', function () {
    var viewer = new Viewer();
    viewer.model = model;
    viewer.frame = data.fps + 1;
    expect(viewer.step).to.equal(1);
  });

  it('can play to a specific frame', function (done) {
    var viewer = new Viewer();
    viewer.model = model;
    viewer.play(2, 100);
    expect(viewer.frame).to.equal(0);
    setTimeout(function () {
      expect(viewer.frame).to.equal(2);
      done();
    }, 40);
  });

  it('can play a specific step', function () {
    var viewer = new Viewer();
    viewer.model = model;
    chai.spy.on(viewer, 'play');

    viewer.playStep(1);
    expect(viewer.step).to.equal(1);
    expect(viewer.frame).to.equal(data.fps);
    expect(viewer.play).to.have.been.called.with(data.fps * 2);
  });
});
