var KokorigamiViewer = require('../src/viewer.js');
var heartModel = require('./heart.json');

var viewer = new KokorigamiViewer();
document.body.appendChild(viewer.el);
viewer.model = heartModel;
viewer.render();

var reset = document.getElementById('reset');
var play = document.getElementById('play');
var prev = document.getElementById('prev');
var next = document.getElementById('next');

prev.addEventListener('click', function () { viewer.prev(); viewer.playStep(viewer.step); });
next.addEventListener('click', function () { viewer.next(); viewer.playStep(viewer.step); });
play.addEventListener('click', function () { viewer.playStep(viewer.step); });
reset.addEventListener('click', function () { viewer.frame = 0; });
