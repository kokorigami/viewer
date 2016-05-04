var KokorigamiViewer = require('../src/viewer.js');
var heartModel = require('./heart.json');

var viewer = new KokorigamiViewer();
document.body.appendChild(viewer.el);
viewer.model = heartModel;
viewer.play();

viewer.stream.onValue(function (val) {
  console.log(val);
});

var prev = document.getElementById('prev');
var next = document.getElementById('next');
prev.addEventListener('click', viewer.prev);
next.addEventListener('click', viewer.next);
