# NeuroBotanica Budtender Application - Tablet Installation Guide

## Overview
This guide covers the installation and deployment of the NeuroBotanica Budtender Assistant on tablet devices for Nevada dispensary operations. The application supports multiple deployment methods: Progressive Web App (PWA), native Android/iOS apps via Capacitor, and Electron desktop applications.

## Prerequisites

### System Requirements
- **Tablets**: iPad (iOS 14+), Android tablets (Android 8.0+)
- **Storage**: Minimum 2GB free space
- **Network**: WiFi connectivity for initial setup and updates
- **Backend**: NeuroBotanica API server running (localhost:8000 or cloud deployment)

### Development Environment
- **Node.js**: 18.0+ with npm
- **Python**: 3.9+ with FastAPI backend
- **Capacitor CLI**: `npm install -g @capacitor/cli`
- **Android Studio**: For Android builds (optional)
- **Xcode**: For iOS builds (macOS only, optional)

## Installation Methods

### Method 1: Progressive Web App (PWA) - Recommended for Quick Deployment

#### Browser-Based Installation
1. **Start the Application**
   ```bash
   # Backend server
   cd neurobotanica_project
   python -m uvicorn backend.main:app --host 127.0.0.1 --port 8000 --reload

   # Frontend (separate terminal)
   cd frontend
   npm run build
   npm start
   ```

2. **Access on Tablet**
   - Open Safari (iOS) or Chrome (Android) on the tablet
   - Navigate to: `http://[SERVER_IP]:3000`
   - The app will detect tablet mode automatically

3. **Install as PWA**
   - **iOS Safari**: Tap share button → "Add to Home Screen"
   - **Android Chrome**: Tap menu (⋮) → "Add to Home screen" → "Add"
   - The app icon will appear on the home screen

4. **Offline Configuration**
   - The PWA includes service worker caching
   - Core functionality works offline
   - Data syncs when connection is restored

#### PWA Features
- ✅ Native app-like experience
- ✅ Offline functionality
- ✅ Automatic updates
- ✅ No app store approval required
- ✅ Cross-platform compatibility

### Method 2: Native Android App via Capacitor

#### Development Build
1. **Install Dependencies**
   ```bash
   cd frontend
   npm install
   npm install -g @capacitor/cli @capacitor/core @capacitor/android
   ```

2. **Build Web Assets**
   ```bash
   npm run build
   ```

3. **Initialize Capacitor**
   ```bash
   npx cap init "NeuroBotanica Budtender" "com.neurobotanica.dispensary"
   npx cap add android
   ```

4. **Configure Capacitor**
   - The `capacitor.config.js` is already configured
   - Web directory points to `out/` (Next.js export)

5. **Sync and Build**
   ```bash
   npx cap sync android
   npx cap open android
   ```

6. **Build APK**
   ```bash
   cd android
   ./gradlew assembleDebug
   ```

#### Production Build
1. **Generate Signed APK**
   ```bash
   # Configure keystore in capacitor.config.js
   ./gradlew assembleRelease
   ```

2. **Deploy to Tablets**
   - Transfer APK file to tablets
   - Enable "Install from Unknown Sources" in Android settings
   - Install the APK

#### Android-Specific Features
- Native Android notifications
- Camera integration for product photos
- Enhanced offline storage
- Google Play Store distribution option

### Method 3: Native iOS App via Capacitor

#### Prerequisites
- macOS computer with Xcode 14+
- Apple Developer Account (for distribution)
- iOS device for testing

#### Development Build
1. **Install Dependencies**
   ```bash
   cd frontend
   npm install
   npm install -g @capacitor/cli @capacitor/core @capacitor/ios
   ```

2. **Build Web Assets**
   ```bash
   npm run build
   ```

3. **Initialize Capacitor**
   ```bash
   npx cap add ios
   npx cap sync ios
   ```

4. **Open in Xcode**
   ```bash
   npx cap open ios
   ```

5. **Configure Signing**
   - Open the project in Xcode
   - Select team in Signing & Capabilities
   - Configure bundle identifier

6. **Build and Run**
   - Select target device/tablet
   - Build (⌘B) and run (⌘R)

#### Production Build
1. **Archive for App Store**
   ```bash
   # In Xcode: Product → Archive
   ```

2. **Distribute**
   - Upload to App Store Connect
   - Or create Ad Hoc distribution for direct installation

#### iOS-Specific Features
- Native iOS design integration
- Face ID/Touch ID authentication
- iCloud sync capabilities
- App Store distribution

### Method 4: Electron Desktop Application

#### Build Process
1. **Install Electron Builder**
   ```bash
   cd frontend
   npm install
   ```

2. **Build Application**
   ```bash
   npm run build:electron
   ```

3. **Package for Distribution**
   ```bash
   npx electron-builder --win --mac --linux
   ```

4. **Install on Tablets**
   - For Windows tablets: Install the .exe or .msi
   - For macOS: Install the .dmg
   - For Linux: Install the .AppImage

## Configuration

### Backend Connection
Update the API endpoint in the application:

```javascript
// src/utils/api.ts
const API_BASE_URL = process.env.NODE_ENV === 'production'
  ? 'https://your-api-domain.com'
  : 'http://localhost:8000';
```

### Environment Variables
Create `.env.local` for local development:
```
NEXT_PUBLIC_API_URL=http://localhost:8000
NEXT_PUBLIC_ENVIRONMENT=development
```

### Database Configuration
For production deployments, configure PostgreSQL:
```python
# backend/config/database.py
DATABASE_URL = "postgresql://user:password@host:port/database"
```

## Deployment Checklist

### Pre-Installation
- [ ] Backend server running and accessible
- [ ] Database initialized with patient tables
- [ ] API endpoints tested (4/4 passing)
- [ ] Frontend build successful
- [ ] Network connectivity verified

### Tablet Setup
- [ ] Device storage: 2GB+ free space
- [ ] Screen resolution: 1280x800 minimum
- [ ] Browser: Chrome 90+ or Safari 14+
- [ ] Permissions: Camera, storage access granted

### Post-Installation
- [ ] App launches successfully
- [ ] Customer profile creation works
- [ ] Recommendations generate properly
- [ ] Offline mode functions
- [ ] Data syncs with backend

## Troubleshooting

### Common Issues

#### PWA Won't Install
- Clear browser cache and cookies
- Ensure HTTPS in production (required for PWA)
- Check manifest.json is accessible

#### Capacitor Build Fails
- Verify Node.js and npm versions
- Clear node_modules: `rm -rf node_modules && npm install`
- Update Capacitor: `npm install @capacitor/cli@latest`

#### Backend Connection Issues
- Verify API server is running on correct port
- Check CORS configuration in FastAPI
- Test API endpoints with curl/Postman

#### Performance Issues
- Enable production build optimizations
- Implement lazy loading for components
- Optimize images and assets

## Security Considerations

### Data Protection
- All customer data is HIPAA-compliant
- End-to-end encryption for sensitive information
- Secure API authentication with tokens

### Network Security
- HTTPS required for production deployments
- API rate limiting implemented
- Input validation on all endpoints

### Device Security
- Regular security updates
- Remote wipe capability for lost devices
- Encrypted local storage

## Support and Maintenance

### Updates
- PWA: Automatic updates via service worker
- Native Apps: App store updates or manual APK/IPA deployment
- Backend: Rolling updates with zero downtime

### Monitoring
- API health checks: `/health` endpoint
- Error logging: Integrated with backend logging
- Performance metrics: Response time monitoring

### Backup and Recovery
- Database backups: Daily automated backups
- Configuration backups: Version controlled
- Disaster recovery: Multi-region deployment option

## Contact Information

For technical support:
- **Development Team**: neurobotanica@cloakandquill.org
- **Documentation**: [GitHub Repository](https://github.com/cloakandquill/neurobotanica)
- **Issues**: Create GitHub issue with tablet model and error details

---

*Installation Guide v1.0 | January 6, 2026 | NeuroBotanica Budtender Assistant*