
module.exports = function(config) {
  config.set({
    basePath: '',
    frameworks: ['mocha'],
    files: ['test.build.js'],
    exclude: [],
    preprocessors: {},
    reporters: ['progress'],
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
