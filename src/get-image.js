var getPixels = require('get-pixels');

module.exports = function (texture) {
  return new Promise(function (res, rej) {
    getPixels(texture, function (err, pixels) {
      if (err) rej(err);
      else res(pixels);
    });
  });
};
