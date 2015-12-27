(function () {
  var canvas = document.createElement('canvas');
  document.body.appendChild(canvas);
  canvas.addEventListener('mousedown', onMouseDown);

  var gl = twgl.getWebGLContext(canvas);

  var layers = [ [ [ [ 0.5, 0.5, 0 ],
          [ 0, 0.2928932011127472, 0 ],
          [ 0, 1, 0 ],
          [ 1, 1, 0 ] ] ],
      [ [ [ 0.5, 0.5, 0 ],
          [ 1, 1, 0 ],
          [ 0, 1, 0 ],
          [ 0, 0.5, 0 ] ] ],
      [ [ [ 0, 0.7071067690849304, 0 ],
          [ 0.5, 0.5, 0 ],
          [ 0, 0.5, 0 ] ] ],
      [ [ [ -0.2071067988872528, 0.5, 0 ],
          [ 0, 0.2928932011127472, 0 ],
          [ 0.5, 0.5, 0 ] ],
        [ [ -0.2071067988872528, 0.5, 0 ],
          [ 0.5, 0.5, 0 ],
          [ 0, 0.7071067690849304, 0 ] ] ] ];

  var points = _.chain(layers).flatten().flatten().unique().value();

  var layerLines = [
    [],
    [],
    [],
    [
      [[ -0.2071067988872528, 0.5, 0 ], [ 0.5, 0.5, 0 ]]
    ]
  ];

  var pointDivs = createPointDivs(points);
  var faceArributes = getFaceBufferArrays(layers);
  var faceBufferInfo = twgl.createBufferInfoFromArrays(gl, faceArributes);
  var faceProgramInfo = twgl.createProgramInfo(gl, ['vs-face', 'fs-face']);
  var depthProgramInfo = twgl.createProgramInfo(gl, ['vs-face', 'fs-face-depth']);

  var lineAttributes = getLineBufferArrays(layerLines, 0.5); // mountain
  var lineBufferInfo = twgl.createBufferInfoFromArrays(gl, lineAttributes);
  var lineProgramInfo = twgl.createProgramInfo(gl, ['vs-fold', 'fs-fold']);

  var rotationY = 0;
  var rotationX = 0;

  var eye = [0, 1, -4];
  var target = [0, 0.5, 0];
  var up = [0, 1, 0];
  var near = 1;
  var far = 100;

  var projection = twgl.m4.perspective(30 * Math.PI / 180, gl.canvas.clientWidth / gl.canvas.clientHeight, near, far);
  var camera = twgl.m4.lookAt(eye, target, up);
  var view = twgl.m4.inverse(camera);
  var viewProjection = twgl.m4.multiply(view, projection);

  var rotateX = twgl.m4.rotationX(rotationX);
  var rotateY = twgl.m4.rotationY(rotationY);
  var world = twgl.m4.multiply(rotateX, rotateY);

  var worldView = twgl.m4.multiply(world, view);
  var worldViewProjection = twgl.m4.multiply(world, viewProjection);

  requestAnimationFrame(render);

  function getUniforms () {
    rotateX = twgl.m4.rotationX(rotationX);
    rotateY = twgl.m4.rotationY(rotationY);
    world = twgl.m4.multiply(rotateX, rotateY);

    worldView = twgl.m4.multiply(world, view);
    worldViewProjection = twgl.m4.multiply(world, viewProjection);

    return {
      u_near: near,
      u_far: far,
      u_lightWorldPos: [-3, 3, -8],
      u_lightColor: [1, 0.9, 0.8, 1],
      u_ambient: [0, 0, 0, 1],
      u_specular: [1, 1, 0.8, 1],
      u_shininess: 70,
      u_specularFactor: 0.8,
      u_camera: camera,
      u_worldViewProjection: worldViewProjection,
      u_worldRotation: world,
      u_world: world,
      // u_texture: textures.origami
    };
  }

  function isClockwise (vertices) {
    var v3 = twgl.v3;
    var vA = v3.subtract(vertices[1], vertices[0]);
    var vB = v3.subtract(vertices[2], vertices[0]);
    var vC = v3.cross(vB, vA);
    return vC[2] > 0;
  }

  function getNormal (face) {
    var v3 = twgl.v3;
    var vA = v3.subtract(face[1], face[0]);
    var vB = v3.subtract(face[1], face[2]);
    var cross = v3.cross(vA, vB);
    v3.normalize(cross, cross);
    return cross;
  }

  function getFaceBufferArrays (faces) {
    var thickness = 0.005;

    var positions = [];
    var normals = [];

    function upper (vertex) {
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
        for (i = 2; i < face.length; i++) {
          vertices.push(first);
          vertices.push(second);
          vertices.push(face[i]);
          second = face[i];
        }
        return vertices;
      }
    }

    function getBase (face, z) {
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

    function getSides (face) {
      face.forEach(function (vertex, i) {
        var next = i + 1;
        if (next === face.length) {
          next = 0;
        }

        // Compute for sides
        var upperVertex = upper(vertex);
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
    };

    layers.forEach(function (faces, layer) {
      var z = layer * thickness;

      faces.forEach(function (face) {
        var clockwise = isClockwise(face);

        if (!clockwise) {
          face = face.reverse();
        }

        getBase(face, z);
        getSides(face);
      });
    });

    return {
      position: positions,
      normal: normals
    };
  }

  function square (x) {
    return x*x;
  }

  function getLineBufferArrays (layers, type) {
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

  function render (time) {
    time *= 0.001;
    twgl.resizeCanvasToDisplaySize(gl.canvas);
    gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);

    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);

    var uniforms = getUniforms(time);

    //renderPass(gl, depthProgramInfo, faceBufferInfo, uniforms, 'TRIANGLES');
    renderPass(gl, faceProgramInfo, faceBufferInfo, uniforms, 'TRIANGLES');
    renderPass(gl, lineProgramInfo, lineBufferInfo, uniforms, 'LINES');

    renderPoints(points, pointDivs);
    requestAnimationFrame(render);
  }

  function renderPass (gl, programInfo, bufferInfo, uniforms, drawType) {
    gl.useProgram(programInfo.program);
    twgl.setBuffersAndAttributes(gl, programInfo, bufferInfo);
    twgl.setUniforms(programInfo, uniforms);
    twgl.drawBufferInfo(gl, gl[drawType], bufferInfo);
  }

  function onMouseDown (e) {
      var canvas = e.currentTarget;
      canvas.addEventListener('mousemove', onMouseMove);
      canvas.addEventListener('mouseup', onMouseUp);

      var startX = e.offsetX;
      var startY = e.offsetY;
      var currentX = startX;
      var currentY = startY;

      function onMouseMove (e) {
        var nextX = e.offsetX;
        var nextY = e.offsetY;
        pan(nextX - currentX, nextY - currentY);
        currentX = nextX;
        currentY = nextY;
      }

      function onMouseUp (e) {
        canvas.removeEventListener('mousemove', onMouseMove);
        canvas.removeEventListener('mouseup', onMouseUp);
      }
  }

  function pan (dX, dY) {
    rotationY = 10 * dX / canvas.clientWidth + rotationY;
    rotationX = 10 * dY / canvas.clientHeight + rotationX;
  }

  function toPixelClipSpace (gl, point) {
    var pixel = [];
    pixel[0] = (point[0] *  0.5 + 0.5) * gl.canvas.width;
    pixel[1] = (point[1] * -0.5 + 0.5) * gl.canvas.height;
    return pixel;
  }

  function createPointDivs (points) {
    return _.map(points, createPoint);
  }

  function renderPoints(points, divs) {
    _.each(points, function (point, i) {
      renderPoint(point, divs[i]);
    });
  }

  function createPoint (point, id) {
    var pointDiv = document.createElement('div');
    pointDiv.classList.add('floating');
    pointDiv.textContent = id;
    document.body.appendChild(pointDiv);
    return pointDiv;
  }

  function renderPoint (point, div) {
    var adjustedPoint = twgl.v3.create();
    adjustedPoint[0] = point[0];
    adjustedPoint[1] = point[1];
    adjustedPoint[2] = point[2];

    twgl.m4.transformPoint(worldViewProjection, adjustedPoint, adjustedPoint);

    var pixelPoint = toPixelClipSpace(gl, adjustedPoint);

    div.style.left = Math.floor(pixelPoint[0]) + 'px';
    div.style.top = Math.floor(pixelPoint[1]) + 'px';
  }
})();
