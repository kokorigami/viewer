
var interpolate = require('./interpolate.js');

var Model = function (data) {
  /*
    Expected format:
      triangles
      naturals
      frames:
        points
        normals
  */

  return this.set(data);
};

Model.prototype = {};

Model.prototype.set = function (data) {
  data = data || {};
  this.triangles = data.triangles || [];
  this.naturals = data.naturals || [];
  this.frames = data.frames || [];
  this.fps = data.fps || 1;
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
  var naturals = this.naturals;
  var frame = this.frameInterpolate(frameIndex);

  var position = [];
  var normal = [];
  var texcoord = [];

  this.triangles.forEach(function (triangle, triangleIndex) {
    triangle.forEach(function (pointIndex) {
      position.push(frame.points[pointIndex]);
      normal.push(frame.normals[triangleIndex]);
      texcoord.push(naturals[pointIndex]);
    });
  });

  return {
    position: position,
    normal: normal,
    texcoord: texcoord
  };
};

Model.prototype.stepFrames = function (step) {
  step = Math.max(step, 0);
  step = Math.min(step, this.steps - 1);
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

Object.defineProperty(Model.prototype, 'lastFrame', {
  enumerable: true,
  get: function () {
    return this.frames.length - 1;
  }
});

module.exports = Model;
