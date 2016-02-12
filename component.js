function Component(name, _proto, parent, extend) {
  parent = parent || HTMLElement;

  var proto = {};
  var protoKeys = Object.keys(_proto || {});

  protoKeys.forEach(function (key) {
    var attr = _proto[key];
    if (typeof attr === 'function') {
      attr = createFunctionAttribute(attr);
    } else if (typeof attr === 'string') {
      attr = createStringElementAttribute(attr);
    } else if (!isDefinedProperty(attr)) {
      attr = createElementAttribute(attr);
    }

    proto[key] = attr;
  });

  var prototype = Object.create(parent.prototype, proto);
  return document.registerElement(name, {prototype: prototype, extends: extend});
}

module.exports = Component;

function isDefinedProperty(prop) {
  var isDataProp = 'value' in prop || 'writable' in prop;
  var isAccessProp = 'get' in prop || 'set' in prop;
  return isDataProp || isAccessProp;
}

function createElementAttribute(name) {
  return {
    get: function() { return JSON.parse(this.getAttribute(name)); },
    set: function(val) { return this.setAttribute(name, JSON.stringify(val)); }
  };
}

function createStringElementAttribute(name) {
  return {
    get: function() { return this.getAttribute(name); },
    set: function(val) { return this.setAttribute(name, val); }
  };
}

function createFunctionAttribute(fn) {
  return {
    value: fn
  };
}
