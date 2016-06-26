
var template = `
  <button class="reset" type="button">reset</button>
  <button class="prev" type="button">&lt;</button>
  <button class="next" type="button">&gt;</button>
  <input class="frame" disabled></input>
  <div class="progress-bar">
    <div class="progress"></div>
  </div>
`;

function ActionBar(viewer, el) {
  this.el = el || document.createElement('div');
  this.el.innerHTML = template;
  this.el.classList.add('kokorigami-action-bar');

  this.viewer = viewer;

  this.reset = this.el.querySelector('.reset');
  this.prev = this.el.querySelector('.prev');
  this.next = this.el.querySelector('.next');
  this.frame = this.el.querySelector('.frame');
  this.bar = this.el.querySelector('.progress');

  this._onClickPrev = function () { viewer.prev(); };
  this._onClickNext = function () { viewer.next(); };
  this._onClickReset = function () { viewer.frame = 0; };
  this.update = this.update.bind(this);

  this.prev.addEventListener('click', this._onClickPrev);
  this.next.addEventListener('click', this._onClickNext);
  this.reset.addEventListener('click', this._onClickReset);
  viewer.on('update', this.update);

  this.update();
  return this;
}

ActionBar.prototype = {};

ActionBar.prototype.update = function () {
  var frame = this.viewer.frame;
  var total = this.viewer.model.frames.length;
  var progress = 100 * frame / total;
  this.frame.value = frame;
  this.bar.style.width = progress + '%';
};

ActionBar.prototype.teardown = function () {
  var reset = this.el.querySelector('.reset');
  var prev = this.el.querySelector('.prev');
  var next = this.el.querySelector('.next');

  prev.removeEventListener('click', this._onClickPrev);
  next.removeEventListener('click', this._onClickNext);
  reset.removeEventListener('click', this._onClickReset);
  this.viewer.off('update', this.update);
};

module.exports = ActionBar;
