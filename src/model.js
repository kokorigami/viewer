
var v3 = require('gl-vec3');

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
  var points, normals;
  if (Number.isInteger(frameIndex)) {
    points = this.frames[frameIndex].points;
    normals = this.frames[frameIndex].normals;
  } else {
    var start = parseInt(frameIndex);
    var end = start + 1;
    var amt = frameIndex - start;
    var startFrames = this.frames[start];
    var endFrames = this.frames[end];
    points = interpolateV3(startFrames.points, endFrames.points, amt);
    normals = interpolateV3(startFrames.normals, endFrames.normals, amt);
  }
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

function interpolateV3(start, end, amt) {
  var set = [];
  start.forEach(function (pointA, i) {
    var point = v3.create();
    v3.lerp(point, pointA, end[i], amt);
    set.push(point);
  });
  return set;
}

