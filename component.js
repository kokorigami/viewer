function Component(name, _proto, parent, extend) {
  var proto = {};
  var protoKeys = Object.keys(_proto || {});

  protoKeys.forEach(function (key) {
    var attr = _proto[key];
    var attachment = attr;
    if (typeof attr === 'function') {
      attachment = { value: attr };
    } else if (attr.get == null && attr.set == null) {
      var isstring = typeof attr == 'string';
      var def = isstring ? attr : JSON.stringify(attr);
      attachment = {
        get: function() {
          var result = this.hasAttribute(key) ? this.getAttribute(key) : def;
          if (isstring) return result;
          return JSON.parse(result);
        }, set: function(val) {
          return this.setAttribute(key, isstring ? val : JSON.stringify(val));
        }
      };
    }

    proto[key] = attachment;
  });

  parent = parent || HTMLElement;
  var prototype = Object.create(parent.prototype, proto);
  return document.registerElement(name, {prototype: prototype, extends: extend});
}

module.exports = Component;

