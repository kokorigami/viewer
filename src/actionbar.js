
var template = `
  <button class="reset" type="button">reset</button>
  <button class="prev" type="button">&lt;</button>
  <button class="next" type="button">&gt;</button>
  <input class="frame" disabled></input>
`;

function Actionbar(viewer, el) {
  var actionbar = el || document.createElement('div');
  actionbar.innerHTML = template;
  actionbar.classList.add('action-bar');

  var reset = actionbar.querySelector('.reset');
  var prev = actionbar.querySelector('.prev');
  var next = actionbar.querySelector('.next');
  var frame = actionbar.querySelector('.frame');

  prev.addEventListener('click', function () { viewer.prev(); });
  next.addEventListener('click', function () { viewer.next(); });
  reset.addEventListener('click', function () { viewer.frame = 0; });

  viewer.on('update', function (newFrame) { frame.value = newFrame; });

  return actionbar;
}

module.exports = Actionbar;
