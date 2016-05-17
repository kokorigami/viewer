// Test bundle generator
var fs = require('fs');
var tests = fs.readdirSync('test');
var contents = '';
contents += '// Test bundle generated on ' + new Date() + '\n';

tests.forEach(function (file) {
  var segments = file.split('.');
  if (segments[segments.length - 1] === 'js') {
    contents += 'require(\'./test/' + file + '\');\n';
  }
});

fs.writeFile('test.bundle.js', contents,'utf8');
