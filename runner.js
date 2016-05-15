
var flags = require('flags');
flags.defineInteger('timeout', 0, 'Connection timeout to start tests.');
flags.defineBoolean('keep', false, 'Keep the test server up after tests are run.');
flags.defineString('reporter', 'spec', 'Test reporter.');
flags.parse();

var mocha = require('mocha');
var express = require('express');
var bodyParser = require('body-parser');
var events = require('events');
var app = express();
var started = false;
var Reporter = mocha.reporters[flags.get('reporter')];
var reporter, runner;

if (flags.get('timeout') > 0) {
  setTimeout(function () {
    if (!started) finish('Connection timed out.', 1);
  }, flags.get('timeout'));
}

app.use(express.static(__dirname));
app.use(bodyParser.json());

app.get('/', function (req, res) {
  started = true;
  res.sendFile(__dirname + '/runner.html');
});

app.post('/start', function (req, res) {
  res.sendStatus(202);
  runner = new events.EventEmitter(); // mock runner
  reporter = new Reporter(runner);
  runner.emit('start');
});

app.post('/progress', function (req, res) {
  res.sendStatus(202);
  var event = req.body.event;
  var data = req.body.data;
  var err = req.body.err;

  // mock test functions
  if (data && data.slow) {
    var slow = data.slow;
    data.slow = function () { return slow; };
  }
  if (data && data.fullTitle) {
    var fullTitle = data.fullTitle;
    data.fullTitle = function () { return fullTitle; };
  }

  reporter.stats = req.body.stats;
  runner.emit(event, data, err);
});

app.post('/end', function (req, res) {
  res.sendStatus(202);
  reporter.stats = req.body;
  runner.emit('end');
  finish(null, reporter.stats.failures > 0 ? 1 : 0);
});

app.listen(3000, function () {
  print('Testing at http://localhost:3000');
});

function print(message, error) {
  var output = error ? process.stderr : process.stdout;
  output.write(message + '\n');
}

function finish(message, code) {
  if (message) print(message, code);
  if (!flags.get('keep')) process.exit(code);
}

process.on('exit', function (code) {
  print('Exited with code ' + code + '\n');
});
