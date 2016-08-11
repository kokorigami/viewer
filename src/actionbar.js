
var template = `
  <button class="reset" type="button">reset</button>
  <button class="prev" type="button">&lt;</button>
  <button class="next" type="button">&gt;</button>
  <input class="frame" disabled></input>
  <div class="progress-bar">
    <div class="progress"></div>
    <div class="steps"></div>
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
  this.steps = this.el.querySelector('.steps');

  this._onClickPrev = function () { viewer.playStep(viewer.step - 1); };
  this._onClickNext = function () { viewer.playStep(viewer.step); };
  this._onClickReset = function () { viewer.frame = 0; };
  this.swap = this.swap.bind(this);
  this.update = this.update.bind(this);

  this.prev.addEventListener('click', this._onClickPrev);
  this.next.addEventListener('click', this._onClickNext);
  this.reset.addEventListener('click', this._onClickReset);
  viewer.on('swap', this.swap);
  viewer.on('update', this.update);

  this.swap();
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

ActionBar.prototype.swap = function () {
  var total = this.viewer.model.frames.length;
  var steps = this.viewer.model.steps;
  var stepEl;
  var stepFrames;
  var stepStart;
  var stepLength;

  this.steps.innerHTML = '';
  for (var i = 0; i < steps; i++) {
    stepEl = document.createElement('div');
    stepFrames = this.viewer.model.stepFrames(i);
    stepStart = 100 * stepFrames[0] / total;
    stepLength = 100 * (stepFrames[1] - stepFrames[0] + 1) / total;
    stepEl.style.left = stepStart + '%';
    stepEl.style.width = stepLength + '%';
    this.steps.appendChild(stepEl);
  }
};

ActionBar.prototype.teardown = function () {
  var reset = this.el.querySelector('.reset');
  var prev = this.el.querySelector('.prev');
  var next = this.el.querySelector('.next');

  prev.removeEventListener('click', this._onClickPrev);
  next.removeEventListener('click', this._onClickNext);
  reset.removeEventListener('click', this._onClickReset);
  this.viewer.off('update', this.update);
  this.viewer.off('swap', this.swap);
};

module.exports = ActionBar;
