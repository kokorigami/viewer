
var express = require('express');
var bodyParser = require('body-parser');
var app = express();

var flag = require('flags');
flag.defineInteger('timeout', 0, 'Timeout to start tests.');
flag.defineBoolean('keep', false, 'Keep the test server up after tests are run.');
// TODO: actually keep the server running
flag.parse();

app.use(express.static(__dirname));
app.use(bodyParser.json());

app.post('/results', function (req, res) {
  res.sendStatus(202);
  var result = req.body;
  if (result.failures) {
    print(result.failures + ' tests failed.\n', 1)
  }
  else {
    print('All tests passed!\n');
    process.exit(0);
  }
});

app.get('/', function (req, res) {
  started = true;
  res.sendFile(__dirname + '/runner.html');
});

app.listen(3000, function () {
  print('Test runner server at http://localhost:3000\n');
});

var started = false;
if (flag.get('timeout') > 0) {
  setTimeout(function () {
    if (!started) {
      print('Test runner did not start within timeout.', 1);
      process.exit(1);
    }
  }, flag.get('timeout'));
}

function print(message, error) {
  var output = error ? process.stderr : process.stdout;
  output.write(message);
}

process.on('exit', function (code) {
  print('Exited: ' + code + '\n');
});
