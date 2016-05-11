var expect = require('chai').expect;
var Viewer = require('../src/viewer.js');

describe('viewer', function () {
  it('can be created', function () {
    var viewer = new Viewer();
    expect(viewer).to.exist;
  });
});
