var KokorigamiViewer = require('../viewer.js');
var heartModel = require('./heart.json');

console.log(heartModel);

var viewer = new KokorigamiViewer();
document.body.appendChild(viewer.el);
viewer.data = heartModel;