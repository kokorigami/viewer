var createGeometry = require('gl-geometry');

function createQuadBufferInfo (gl) {
  var geom = createGeometry(gl)
    .attr('position', [
      -1, -1, 0,
      -1,  1, 0,
       1,  1, 0,
       1,  1, 0,
       1, -1, 0,
      -1, -1, 0
    ]);
  return geom;
}

module.exports = createQuadBufferInfo;
