var expect = require('chai').expect;
var Model = require('../src/model.js');

describe('model', function () {
  it('can be created', function () {
    var model = new Model();
    expect(model).to.exist;
  });
});
