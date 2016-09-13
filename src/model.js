
var interpolate = require('./interpolate.js');
var getPlaneNormal = require('get-plane-normal');

function computeFrameNormals(frame) {
  var points = frame.points;
  var normals = [];

  for (var i = 0; i < points.length; i += 3) {
    var normal = getPlaneNormal([], points[i], points[i+1], points[i+2]);
    normals[i] = normals[i+1] = normals[i+2] = normal;
  }

  frame.normals = normals;
}

var Model = function (data) {
  /*
    Expected format:
      naturals
      frames:
        points
  */

  return this.set(data);
};

Model.prototype = {};

Model.prototype.set = function (data) {
  data = data || {};
  this.naturals = data.naturals || [];
  this.frames = data.frames || [];
  this.fps = data.fps || 1;
  this.frames.forEach(computeFrameNormals);
  return this;
};

Model.prototype.frameInterpolate = function (frameIndex) {
  if (Number.isInteger(frameIndex)) return this.frames[frameIndex];

  var start = parseInt(frameIndex);
  var end = start + 1;
  var amt = frameIndex - start;
  var startFrames = this.frames[start];
  var endFrames = this.frames[end];
  var points = interpolate(startFrames.points, endFrames.points, amt);
  var normals = interpolate(startFrames.normals, endFrames.normals, amt);
  return {
    points: points,
    normals: normals
  };
};

Model.prototype.frameGeometry = function (frameIndex) {
  var frame = this.frameInterpolate(frameIndex);
  return {
    position: frame.points,
    normal: frame.normals,
    texcoord: this.naturals
  };
};

Model.prototype.stepFrames = function (step) {
  return [step * this.fps, step * this.fps + this.fps - 1];
};

Model.prototype.stepFromFrame = function (frame) {
  return Math.floor(frame / this.fps);
};

Object.defineProperty(Model.prototype, 'steps', {
  enumerable: true,
  get: function () {
    return this.frames.length / this.fps;
  }
});

Object.defineProperty(Model.prototype, 'final', {
  enumerable: true,
  get: function () {
    return this.frames.length - 1;
  }
});

module.exports = Model;
