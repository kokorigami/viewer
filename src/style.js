
var css = require('./style.css');

function style() {
  var el = document.querySelector('style');
  if (!el) {
    el = document.createElement('style');
    document.body.appendChild(el);
  }

  el.innerHTML += css;
  return el;
}

module.exports = style;
