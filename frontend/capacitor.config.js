const { CapacitorConfig } = require('@capacitor/cli');

const config = {
  appId: 'com.neurobotanica.dispensary',
  appName: 'NeuroBotanica Budtender',
  webDir: 'out',
  bundledWebRuntime: false,
  plugins: {
    SplashScreen: {
      launchShowDuration: 3000,
      launchAutoHide: true,
      backgroundColor: '#0f172a',
      androidSplashResourceName: 'splash',
      androidScaleType: 'CENTER_CROP',
      showSpinner: true,
      androidSpinnerStyle: 'large',
      iosSpinnerStyle: 'small',
      spinnerColor: '#1e293b',
    },
  },
  android: {
    buildOptions: {
      keystorePath: null,
      keystoreAlias: null,
      keystorePassword: null,
      keystoreAliasPassword: null,
    },
  },
  ios: {
    scheme: 'NeuroBotanica',
  },
};

module.exports = config;