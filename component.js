function Component(name, parent, _proto, extend) {
  var proto = {};
  var protoKeys = Object.keys(_proto || {});

  protoKeys.forEach(function (key) {
    var attr = _proto[key];
    if (typeof attr === 'string') {
      attr = createElementAttribute(attr);
    } else if (typeof attr === 'function') {
      attr = createFunctionAttribute(attr);
    }

    proto[key] = attr;
  });

  var prototype = Object.create(parent.prototype, proto);
  return document.registerElement(name, {prototype: prototype, extends: extend});
}

module.exports = Component;

function createElementAttribute(name) {
  return {
    get: function() { return JSON.parse(this.getAttribute(name)); },
    set: function(val) { return this.setAttribute(name, JSON.stringify(val)); }
  };
}

function createFunctionAttribute(fn) {
  return {
    value: fn
  };
}
