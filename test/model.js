
var expect = require('chai').expect;
var Model = require('../src/model.js');
var interpolate = require('../src/interpolate');
var data = require('./data.json');
var EPSILON = 0.000001;

describe('model', function () {
  var model;

  beforeEach(function () {
    model = new Model(data);
  });

  it('can be created from an exported model', function () {
    expect(model.triangles).to.equal(data.triangles);
    expect(model.naturals).to.equal(data.naturals);
    expect(model.frames).to.equal(data.frames);
    expect(model.fps).to.equal(data.fps);
  });

  it('counts the number of steps', function () {
    expect(model.steps).to.equal(2);
  });

  it('provides the index of the final frame', function () {
    expect(model.final).to.equal(11);
  });

  it('interpolates between frames', function () {
    var interpolated = model.frameInterpolate(1.55);
    var start = data.frames[1];
    var end = data.frames[2];
    var points = interpolate(start.points, end.points, 0.55);
    var normals = interpolate(start.normals, end.normals, 0.55);
    expectPointsCloseTo(interpolated.points, points);
    expectPointsCloseTo(interpolated.normals, normals);
  });

  it('generates geometry for frames', function () {
    var geometry = model.frameGeometry(0);
    expectPointsCloseTo(geometry.position, [[0,0,0],[1,0,0],[1,1,0],[0,0,0],[1,1,0],[0,1,0]]);
    expectPointsCloseTo(geometry.normal, [[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]]);
  });

  it('returns the range of frames for a step', function () {
    var frameRange = model.stepFrames(1);
    expect(frameRange).to.deep.equal([6, 11]);
  });

  it('returns the step given a frame index', function () {
    var step = model.stepFromFrame(9);
    expect(step).to.equal(1);
  });
});

function expectPointsCloseTo(points, expected) {
  points.forEach(function (point, i) {
    point.forEach(function (value, j) {
      expect(value).to.be.closeTo(expected[i][j], EPSILON);
    });
  });
}
