
var template = `
  <button class="reset" type="button">reset</button>
  <button class="play" type="button">&gt;</button>
  <input class="frame" disabled></input>
  <div class="progress-bar">
    <div class="bar"></div>
    <div class="steps"></div>
  </div>
`;

function ActionBar(viewer, el) {
  this.el = el || document.createElement('div');
  this.el.innerHTML = template;
  this.el.classList.add('kokorigami-action-bar');

  this.viewer = viewer;

  this.reset = this.el.querySelector('.reset');
  this.play = this.el.querySelector('.play');
  this.frame = this.el.querySelector('.frame');
  this.bar = this.el.querySelector('.bar');
  this.steps = this.el.querySelector('.steps');

  this._onClickPlay = function () { viewer.playStep(viewer.step); };
  this._onClickReset = function () { viewer.frame = 0; };
  this._onClickSteps = function (e) {
    var stepEl = e.target;
    var step = stepEl.dataset.index;
    var startFrame = viewer.model.stepFrames(step)[0];
    viewer.frame = startFrame;
  };
  this.swap = this.swap.bind(this);
  this.update = this.update.bind(this);

  this.play.addEventListener('click', this._onClickPlay);
  this.reset.addEventListener('click', this._onClickReset);
  this.steps.addEventListener('click', this._onClickSteps);
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
    stepEl.dataset.index = i;
    this.steps.appendChild(stepEl);
  }
};

ActionBar.prototype.teardown = function () {
  this.play.removeEventListener('click', this._onClickNext);
  this.reset.removeEventListener('click', this._onClickReset);
  this.steps.removeEventListener('click', this._onClickSteps);
  this.viewer.off('update', this.update);
  this.viewer.off('swap', this.swap);
};

module.exports = ActionBar;
