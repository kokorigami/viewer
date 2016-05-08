
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

Model.prototype.frameGeometry = function (frameIndex) {
  var points = this.frames[frameIndex].points;
  var normals = this.frames[frameIndex].normals;
  var naturals = this.naturals;

  var position = [];
  var normal = [];
  var texcoord = [];

  this.triangles.forEach(function (triangle, triangleIndex) {
    triangle.forEach(function (pointIndex) {
      position.push(points[pointIndex]);
      normal.push(normals[triangleIndex]);
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
