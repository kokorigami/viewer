
function attachControls(el, camera) {
  var controls = {
    mousedown: controlCamera(camera),
    wheel: controlZoom(camera)
  };
  el.addEventListener('wheel', controls.wheel);
  el.addEventListener('mousedown', controls.mousedown);
  return controls;
}

function removeControls(el, controls) {
  el.removeEventListener('mousedown', controls.mousedown);
  el.removeEventListener('wheel', controls.wheel);
}

function controlCamera(camera) {
  var onMouseDown = function(ev) {
    var canvas = ev.currentTarget;
    canvas.addEventListener('mousemove', onMouseMove);
    canvas.addEventListener('mouseup', stopListening);
    canvas.addEventListener('mouseleave', stopListening);

    var startX = ev.offsetX;
    var startY = ev.offsetY;
    var currentX = startX;
    var currentY = startY;

    function onMouseMove(mev) {
      var dx = (mev.offsetX - currentX) * 4/ canvas.width;
      var dy = (mev.offsetY - currentY) * 4/ canvas.height;

      camera.rotate(Date.now(), -dx, dy);
      currentX = mev.offsetX;
      currentY = mev.offsetY;
    }

    function stopListening() {
      canvas.removeEventListener('mousemove', onMouseMove);
      canvas.removeEventListener('mouseup', stopListening);
      canvas.removeEventListener('mouseleave', stopListening);
    }
  };
  return onMouseDown;
}

function controlZoom(camera) {
  var onWheel = function (e) {
    camera.pan(Date.now(), 0, 0, e.deltaY / 50);
  };
  return onWheel;
}

module.exports = {
  attach: attachControls,
  remove: removeControls
};
