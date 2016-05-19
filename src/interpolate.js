var v3 = require('gl-vec3');

function interpolatePoints(start, end, amt) {
  var set = [];
  start.forEach(function (pointA, i) {
    var point = v3.create();
    v3.lerp(point, pointA, end[i], amt);
    set.push(point);
  });
  return set;
}

module.exports = interpolatePoints;
