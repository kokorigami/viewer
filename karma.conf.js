
module.exports = function(config) {
  config.set({
    basePath: '',
    frameworks: ['browserify', 'mocha'],
    files: ['test/*.js'],
    exclude: [],
    preprocessors: {
      'test/*.js': ['browserify']
    },
    reporters: ['mocha'],
    port: 9876,
    colors: true,
    logLevel: config.LOG_INFO,
    autoWatch: false,
    browsers: ['Firefox'],
    singleRun: true,
    concurrency: Infinity,
    client: {
      mocha: {
        reporter: 'html',
        ui: 'bdd'
      }
    }
  });
};
