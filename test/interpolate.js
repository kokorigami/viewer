
var expect = require('chai').expect;
var interpolate = require('../src/interpolate.js');
var EPSILON = 0.000001;

describe('interpolate', function () {
  it('interpolates sets of points', function () {
    var setA = [[1, 1, 1], [0, 0, 0]];
    var setB = [[0, 0, 0], [1, 1, 1]];

    var interpolated = interpolate(setA, setB, 0.6);
    expectPointsCloseTo(interpolated, [[0.4, 0.4, 0.4], [0.6, 0.6, 0.6]]);
  });
});

function expectPointsCloseTo(points, expected) {
  points.forEach(function (point, i) {
    point.forEach(function (value, j) {
      expect(value).to.be.closeTo(expected[i][j], EPSILON);
    });
  });
}
