var KokorigamiViewer = require('../src/viewer.js');
var heartModel = require('./heart.json');

var viewer = new KokorigamiViewer();
document.body.appendChild(viewer.el);
viewer.model = heartModel;
viewer.play();

var prev = document.getElementById('prev');
var next = document.getElementById('next');
prev.addEventListener('click', function () { viewer.prev(); });
next.addEventListener('click', function () { viewer.next(); });
