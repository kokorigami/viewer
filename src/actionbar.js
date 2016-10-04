
var template = `
  <button class="play" type="button">
    <svg height="100" width="100" viewBox="0 0 100 100">
      <polygon points="20,20 80,50 20,80"/>
    </svg>
  </button>
  <button class="pause" type="button">
    <svg width="100" height="100" viewBox="0 0 100 100">
      <rect x="20" y="20" width="15" height="60" />
      <rect x="55" y="20" width="15" height="60" />
    </svg>
  </button>
  <div class="progress">
    <div class="progress-bar"><div class="bar"></div></div>
    <div class="steps"></div>
  </div>
`;

function ActionBar(viewer, el) {
  this.el = el || document.createElement('div');
  this.el.innerHTML = template;
  this.el.classList.add('kokorigami-action-bar');

  this.viewer = viewer;

  this.play = this.el.querySelector('.play');
  this.pause = this.el.querySelector('.pause');
  this.progress = this.el.querySelector('.progress');
  this.bar = this.progress.querySelector('.bar');
  this.steps = this.progress.querySelector('.steps');

  this._onClickPlay = function () { viewer.playStep(viewer.step); };
  this._onClickPause = function () { viewer.stop(); };
  this._onScrub = function (r) {
    if (r.button) return;
    var target = r.currentTarget;
    var width = target.clientWidth;

    var scrub = function (e) {
      var amt = e.offsetX / width;
      var total = viewer.model.frames.length;
      var frame = Math.round(amt * total);
      viewer.frame = frame;
    };

    var endScrub = function () {
      target.removeEventListener('mousemove', scrub);
      target.removeEventListener('mouseup', endScrub);
      target.removeEventListener('mouseout', endScrub);
    };

    target.addEventListener('mousemove', scrub);
    target.addEventListener('mouseup', endScrub);
    target.addEventListener('mouseout', endScrub);

    viewer.stop();
    scrub(r);
  };

  this.swap = this.swap.bind(this);
  this.update = this.update.bind(this);
  this.onStop = this.onStop.bind(this);
  this.onPlay = this.onPlay.bind(this);

  this.play.addEventListener('click', this._onClickPlay);
  this.pause.addEventListener('click', this._onClickPause);
  this.progress.addEventListener('mousedown', this._onScrub);
  viewer.on('swap', this.swap);
  viewer.on('update', this.update);
  viewer.on('stop', this.onStop);
  viewer.on('play', this.onPlay);

  this.swap();
  this.update();
  return this;
}

ActionBar.prototype = {};

ActionBar.prototype.update = function () {
  var frame = this.viewer.frame;
  var total = this.viewer.model.frames.length;
  var progress = 100 * frame / total;
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
  this.onStop();
};

ActionBar.prototype.teardown = function () {
  this.play.removeEventListener('click', this._onClickNext);
  this.pause.removeEventListener('click', this._onClickPause);
  this.progress.removeEventListener('mousedown', this._onScrub);
  this.viewer.off('update', this.update);
  this.viewer.off('swap', this.swap);
  this.viewer.off('stop', this.onStop);
  this.viewer.off('play', this.onPlay);
};

ActionBar.prototype.onStop = function () {
  this.pause.style.display = 'none';
  this.play.style.display = 'block';
};

ActionBar.prototype.onPlay = function () {
  this.pause.style.display = 'block';
  this.play.style.display = 'none';
};

module.exports = ActionBar;
