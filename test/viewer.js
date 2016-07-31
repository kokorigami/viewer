
var chai = require('chai');
var spies = require('chai-spies');
var expect = chai.expect;
chai.use(spies);
var Viewer = require('../src/viewer.js');

describe('viewer', function () {
  it('can be created', function () {
    var viewer = new Viewer();
    expect(viewer).to.exist;
  });

  it('emits update when the frame is updated', function () {
    var viewer = new Viewer();
    var callback = chai.spy();

    viewer.on('update', callback);
    viewer.frame = 1;
    expect(callback).to.have.been.called();
  });

  it('can set the background color', function () {
    var viewer = new Viewer();
    viewer.background = [0, 255, 255, 102];
    expect(viewer.el.style.backgroundColor).to.equal('rgba(0, 255, 255, 0.4)');
  });
});
