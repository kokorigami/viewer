var _ = require('underscore')._;
var twgl = require('twgl.js');
var v3 = twgl.v3;

function createFaceBufferInfo (gl, facesPerLayer) {
  var faceAttributes = getFaceBufferArrays(facesPerLayer);
  var faceBufferInfo = twgl.createBufferInfoFromArrays(gl, faceAttributes);
  return faceBufferInfo;
}

module.exports = createFaceBufferInfo;

function getFaceBufferArrays (layers) {
  var thickness = 0.005;
  var positions = [];
  var normals = [];

  layers.forEach(function (faces, layer) {
    var z = layer * thickness;

    faces.forEach(function (face) {
      var clockwise = isClockwise(face);

      if (!clockwise) {
        face = face.reverse();
      }

      getBase(face, z, thickness, positions, normals);
      getSides(face, thickness, positions, normals);
    });
  });

  return {
    position: positions,
    normal: normals
  };
}

function getBase (face, z, thickness, positions, normals) {
  var vertices = triangulate(face);
  var mirrorVertices = mirror(vertices);

  var normal = getNormal(vertices);
  var mirrorNormal = getNormal(mirrorVertices);

  vertices.forEach(function (vertex) {
    vertex[2] = z;
    normals.push.apply(normals, normal);
    positions.push.apply(positions, vertex);
  });

  mirrorVertices.forEach(function (vertex) {
    vertex[2] = z - thickness;
    normals.push.apply(normals, mirrorNormal);
    positions.push.apply(positions, vertex);
  });
}

function getSides (face, thickness, positions, normals) {
  face.forEach(function (vertex, i) {
    var next = i + 1;
    if (next === face.length) {
      next = 0;
    }

    // Compute for sides
    var upperVertex = upper(vertex, thickness);
    var nextVertex = face[next];
    var nextUpperVertex = upper(nextVertex);

    var normalSideA = getNormal([upperVertex, vertex, nextUpperVertex]);
    var normalSideB = getNormal([vertex, nextVertex, nextUpperVertex]);

    positions.push.apply(positions, upperVertex);
    positions.push.apply(positions, vertex);
    positions.push.apply(positions, nextUpperVertex);
    positions.push.apply(positions, vertex);
    positions.push.apply(positions, nextVertex);
    positions.push.apply(positions, nextUpperVertex);

    normals.push.apply(normals, normalSideA);
    normals.push.apply(normals, normalSideA);
    normals.push.apply(normals, normalSideA);
    normals.push.apply(normals, normalSideB);
    normals.push.apply(normals, normalSideB);
    normals.push.apply(normals, normalSideB);
  });
}

function isClockwise (vertices) {
  var vA = v3.subtract(vertices[1], vertices[0]);
  var vB = v3.subtract(vertices[2], vertices[0]);
  var vC = v3.cross(vB, vA);
  return vC[2] > 0;
}

function getNormal (face) {
  var vA = v3.subtract(face[1], face[0]);
  var vB = v3.subtract(face[2], face[0]);
  var cross = v3.cross(vA, vB);
  v3.normalize(cross, cross);
  return cross;
}

function upper (vertex, thickness) {
  var upper = _.clone(vertex);
  upper[2] = vertex[2] + thickness;
  return upper;
}

function mirror (vertices) {
  var upperVertices = _.clone(vertices).reverse();
  return upperVertices;
}

function triangulate (face) {
  if (face.length < 3) return false;
  else if (face.length === 3) return face;
  else {
    var vertices = [];
    var first = face[0];
    var second = face[1];
    for (var i = 2; i < face.length; i++) {
      vertices.push(first);
      vertices.push(second);
      vertices.push(face[i]);
      second = face[i];
    }
    return vertices;
  }
}
