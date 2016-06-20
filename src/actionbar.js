
var template = `
  <button class="reset" type="button">reset</button>
  <button class="prev" type="button">&lt;</button>
  <button class="next" type="button">&gt;</button>
  <input class="frame" disabled></input>
`;

function ActionBar(viewer, el) {
  this.el = el || document.createElement('div');
  this.el.innerHTML = template;
  this.el.classList.add('action-bar');

  this.viewer = viewer;

  this._onClickPrev = function () { viewer.prev(); };
  this._onClickNext = function () { viewer.next(); };
  this._onClickReset = function () { viewer.frame = 0; };
  this._onUpdate = function (newFrame) { frame.value = newFrame; };

  var reset = this.el.querySelector('.reset');
  var prev = this.el.querySelector('.prev');
  var next = this.el.querySelector('.next');
  var frame = this.el.querySelector('.frame');

  prev.addEventListener('click', this._onClickPrev);
  next.addEventListener('click', this._onClickNext);
  reset.addEventListener('click', this._onClickReset);
  viewer.on('update', this._onUpdate);

  return this;
}

ActionBar.prototype = {};

ActionBar.prototype.teardown = function () {
  var reset = this.el.querySelector('.reset');
  var prev = this.el.querySelector('.prev');
  var next = this.el.querySelector('.next');

  prev.removeEventListener('click', this._onClickPrev);
  next.removeEventListener('click', this._onClickNext);
  reset.removeEventListener('click', this._onClickReset);
  this.viewer.off('update', this._onUpdate);
};

module.exports = ActionBar;
