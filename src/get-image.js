var getPixels = require('get-pixels');
var white = new Uint8ClampedArray([255, 255, 255, 255]);
var block = new ImageData(white, 1, 1);

module.exports = function (texture) {
  if (!texture) return Promise.resolve(block);
  return new Promise(function (res, rej) {
    getPixels(texture, function (err, pixels) {
      if (err) rej(err);
      else res(pixels);
    });
  });
};
