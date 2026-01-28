# GitHub Push Instructions for GECKO

## Your GitHub Account
- **Username:** RishabhGautam22

---

## Step 1: Create Repository on GitHub

1. Go to: https://github.com/new
2. Fill in:
   - **Repository name:** `gecko-web` (or your preferred name)
   - **Description:** `GECKO-A Atmospheric Chemistry Web Interface`
   - **Visibility:** Public (or Private)
   - **DO NOT** initialize with README (you already have one)
3. Click **"Create repository"**

---

## Step 2: Push Code to GitHub

Open Terminal and run these commands:

```bash
# Navigate to project
cd /Users/rgpro/Desktop/GECKO

# Initialize git (if not already done)
git init

# Add all files
git add .

# Create commit
git commit -m "Initial commit: GECKO-A Web Interface v3.0.10

- FastAPI backend with atmospheric chemistry simulation
- GECKO-A and BOXMODEL4GECKO Fortran integration
- 150+ VOC compound database
- RDKit molecular visualization
- Fly.io deployment configuration
- Cross-platform support (macOS, Windows, Linux)

Co-Authored-By: Claude <noreply@anthropic.com>"

# Add GitHub remote
git remote add origin https://github.com/RishabhGautam22/gecko-web.git

# Push to GitHub
git push -u origin main
```

**When prompted for password:** Use your GitHub password or Personal Access Token (PAT)

---

## Step 3: If Authentication Fails

GitHub no longer accepts passwords. Create a Personal Access Token:

1. Go to: https://github.com/settings/tokens
2. Click **"Generate new token (classic)"**
3. Select scopes: `repo` (full control)
4. Copy the token
5. Use the token as your password when pushing

Or use GitHub CLI:
```bash
# Install GitHub CLI
brew install gh

# Login
gh auth login

# Push
git push -u origin main
```

---

## Step 4: Deploy to Fly.io

After pushing to GitHub:

```bash
# Install Fly CLI (if not installed)
brew install flyctl

# Login to Fly.io
fly auth login

# Launch app (uses existing fly.toml)
fly launch --no-deploy

# When prompted:
# - App name: gecko-web (or choose your own)
# - Region: iad (US East) or your preferred region
# - PostgreSQL: No
# - Redis: No
# - Deploy now: No

# Create persistent volume for job data
fly volumes create gecko_data --region iad --size 10

# Deploy (first time takes 45-90 minutes)
fly deploy

# Open in browser
fly open
```

---

## Step 5: Connect GitHub to Fly.io (Auto-Deploy)

### Get Fly.io API Token:
```bash
fly tokens create deploy -x 999999h
```
Copy the token.

### Add to GitHub Secrets:
1. Go to: `https://github.com/RishabhGautam22/gecko-web/settings/secrets/actions`
2. Click **"New repository secret"**
3. Name: `FLY_API_TOKEN`
4. Value: (paste your Fly.io token)
5. Click **"Add secret"**

Now every `git push` to `main` will auto-deploy to Fly.io!

---

## Expected URLs After Deployment

| Resource | URL |
|----------|-----|
| **GitHub Repo** | `https://github.com/RishabhGautam22/gecko-web` |
| **Live App** | `https://gecko-web.fly.dev` |
| **Fly Dashboard** | `https://fly.io/apps/gecko-web` |

---

## Quick Reference Commands

```bash
# Check git status
git status

# Add changes
git add .

# Commit
git commit -m "Your message"

# Push to GitHub (triggers auto-deploy)
git push origin main

# Manual deploy to Fly.io
fly deploy

# View Fly.io logs
fly logs

# Open app in browser
fly open

# Check app status
fly status
```

---

## Troubleshooting

### "remote origin already exists"
```bash
git remote remove origin
git remote add origin https://github.com/RishabhGautam22/gecko-web.git
```

### "failed to push some refs"
```bash
git pull origin main --rebase
git push origin main
```

### "authentication failed"
Use Personal Access Token instead of password, or:
```bash
gh auth login
git push origin main
```

---

## Files Being Uploaded

✅ **Included (~50 MB):**
- All Python code (`gecko_web/`)
- All Fortran source (492 `.f90` files)
- Reaction databases (`docker/gecko_source/DATA/`)
- Docker configuration
- Fly.io configuration
- Documentation

❌ **Excluded (saves ~1.1 GB):**
- Virtual environments (`.venv/`, `venv/`)
- Windows installer (`Gecko Windows Installation/`)
- Compiled Fortran objects (`*.o`, `*.mod`)
- Test data directories

---

*Last updated: 2026-01-28*
