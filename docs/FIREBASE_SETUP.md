# Firebase Setup for NeuroBotanica

This document describes the minimal steps to enable Firebase authentication for the NeuroBotanica MVP and how to add the required secrets to GitHub Actions.

Prerequisites
- You must be an owner or have permission to manage the Google Cloud project `neurobotanica-9ddf3`.
- If your organization enforces the `iam.disableServiceAccountKeyCreation` constraint, either temporarily disable it for this project or use JWKS-only verification (project ID only).

Steps
1. Enable Firebase for an existing GCP project
   - Open https://console.firebase.google.com
   - Click **Add project** → choose **Import Google Cloud project** or select the existing `neurobotanica-9ddf3` project.

2. Create (or use existing) service account and generate a key
   - Open Google Cloud Console → **IAM & Admin** → **Service Accounts**.
   - Locate the Firebase service account: `firebase-adminsdk-<id>@neurobotanica-9ddf3.iam.gserviceaccount.com` or create a new account (name: `neurobotanica-ci-sa`).
   - Click the service account → **Keys** → **Add Key** → **Create new key** → choose **JSON** → **Create**.
   - The browser will download a private JSON key. Save it securely and do NOT commit it to git.

3. If key creation is blocked by organization policy
   - Open Google Cloud Console → **IAM & Admin** → **Organization Policies**.
   - Search for `Disable service account key creation` (`constraints/iam.disableServiceAccountKeyCreation`).
   - Edit the policy to allow key creation for your project, or coordinate with your Org Admin to authorize keys for this project.

4. Add secrets to GitHub (recommended)
   - Repository → Settings → Secrets and variables → Actions → New repository secret.
     - `FIREBASE_PROJECT_ID` = `neurobotanica-9ddf3`
     - `FIREBASE_SERVICE_ACCOUNT_JSON` = *paste the full JSON file contents*
   - Alternatively, use GitHub CLI locally:
     ```bash
     gh secret set FIREBASE_PROJECT_ID --body "neurobotanica-9ddf3"
     gh secret set FIREBASE_SERVICE_ACCOUNT_JSON --body-file ./path/to/service-account.json
     ```

5. Use test mode for local development and CI until production keys are available
   - Set `FIREBASE_TEST_MODE=1` in your environment or repo secrets to enable the codebase test-token bypass.

Security notes
- Never commit the JSON key to source control.
- After verification, restrict service account roles to the least privilege required.
- Rotate keys immediately if exposed.

If you want, I can add the `FIREBASE_SERVICE_ACCOUNT_JSON` secret for you once you confirm the file path to the downloaded key.
  
See also: docs/FIREBASE_JWKS_SNIPPET.md - instructions to enable the optional JWKS refresher for staging. 
