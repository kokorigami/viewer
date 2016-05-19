
function attachControls(el, camera) {
  var controls = {
    mousedown: controlCamera(camera)
  };
  el.addEventListener('mousedown', controls.mousedown);
  return controls;
}

function removeControls(el, controls) {
  el.removeEventListener('mousedown', controls.mousedown);
}

function controlCamera(camera) {
  var onMouseDown = function (ev) {
    var canvas = ev.currentTarget;
    canvas.addEventListener('mousemove', onMouseMove);
    canvas.addEventListener('mouseup', onMouseUp);

    var startX = ev.offsetX;
    var startY = ev.offsetY;
    var currentX = startX;
    var currentY = startY;

    function onMouseMove (mev) {
      var dx = (mev.offsetX - currentX) * 4/ canvas.width;
      var dy = (mev.offsetY - currentY) * 4/ canvas.height;

      camera.rotate(Date.now(), -dx, dy);
      currentX = mev.offsetX;
      currentY = mev.offsetY;
    }

    function onMouseUp () {
      canvas.removeEventListener('mousemove', onMouseMove);
      canvas.removeEventListener('mouseup', onMouseUp);
    }
  };
  return onMouseDown;
}

module.exports = {
  attach: attachControls,
  remove: removeControls
};
