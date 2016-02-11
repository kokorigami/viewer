var twgl = require('twgl.js');

function createXYPlaneBufferInfo (gl, w, h) {
  w = w || 1;
  h = h || 1;

  var planeAttributes = {
    position: [
      0, 0,
      0, h,
      w, h,
      0, 0,
      w, h,
      w, 0
    ]
  };
  planeAttributes.position.numComponents = 2;

  return twgl.createBufferInfoFromArrays(gl, planeAttributes);
}

module.exports = createXYPlaneBufferInfo;
