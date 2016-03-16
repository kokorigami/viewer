var createGeometry = require('gl-geometry');

function createFoldBufferInfo (gl, foldsPerLayer) {
  var foldAttributes = getFoldBufferArrays(foldsPerLayer, 1); // mountain
  var geom = createGeometry(gl)
    .attr('foldType', foldAttributes.foldType, {size: 1})
    .attr('lengthSoFar', foldAttributes.lengthSoFar, {size: 1})
    .attr('position', foldAttributes.position);

  return geom;
}

module.exports = createFoldBufferInfo;

function getFoldBufferArrays (layers, type) {
  var thickness = 0.005;
  type = type || 0.5;

  var types = [];
  var lengths = [];
  var positions = [];

  layers.forEach(function (lines, layer) {
    var z = layer * thickness;

    lines.forEach(function (line) {
      var length = Math.sqrt(square(line[1][0] - line[0][0]) + square(line[1][1] - line[0][1]));

      line[0][2] = z;
      line[1][2] = z;

      positions.push.apply(positions, line[0]);
      positions.push.apply(positions, line[1]);
      types.push(type);
      types.push(type);
      lengths.push(0);
      lengths.push(length);
    });
  });

  var attributes = {
    foldType: types,
    lengthSoFar: lengths,
    position: positions
  };

  attributes.foldType.numComponents = 1;
  attributes.lengthSoFar.numComponents = 1;
  return attributes;
}

function square (x) {
  return x*x;
}
