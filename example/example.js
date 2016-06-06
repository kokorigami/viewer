var KokorigamiViewer = require('../src/viewer.js');
var heartModel = require('./heart.json');

var viewer = new KokorigamiViewer();
document.body.appendChild(viewer.el);

var reset = document.getElementById('reset');
var prev = document.getElementById('prev');
var next = document.getElementById('next');
var frame = document.getElementById('frame');

prev.addEventListener('click', function () { viewer.prev(); });
next.addEventListener('click', function () { viewer.next(); });
reset.addEventListener('click', function () { viewer.frame = 0; });

viewer.emitter.on('update', function (newFrame) {
  frame.value = newFrame;
});

viewer.model = heartModel;
viewer.render();
