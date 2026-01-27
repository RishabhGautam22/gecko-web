# GECKO-A Web Interface - Fly.io Deployment Guide

## Overview

This guide walks you through deploying GECKO to Fly.io with GitHub integration for automatic deployments.

**Estimated Cost:** $10-15/month
**Deployment Time:** ~60-90 minutes (first build), ~5 minutes (subsequent)

---

## Prerequisites

- GitHub account with your GECKO repository
- Credit card for Fly.io (free tier available, but GECKO needs paid resources)
- Terminal/Command Prompt access

---

## Step 1: Create GitHub Repository

If you haven't already pushed to GitHub:

```bash
cd /Users/rgpro/Desktop/GECKO

# Initialize git (if not done)
git init

# Add all files
git add .

# Commit
git commit -m "Initial commit: GECKO-A Web Interface v3.0.6"

# Create repo on GitHub (using GitHub CLI)
gh repo create gecko-web --public --source=. --push

# Or manually:
# 1. Go to https://github.com/new
# 2. Create repository named "gecko-web"
# 3. Push:
git remote add origin https://github.com/YOUR_USERNAME/gecko-web.git
git branch -M main
git push -u origin main
```

---

## Step 2: Install Fly.io CLI

### macOS
```bash
brew install flyctl
```

### Windows (PowerShell)
```powershell
powershell -Command "iwr https://fly.io/install.ps1 -useb | iex"
```

### Linux
```bash
curl -L https://fly.io/install.sh | sh
```

---

## Step 3: Create Fly.io Account & Login

```bash
# Sign up / Login (opens browser)
fly auth signup
# or
fly auth login
```

---

## Step 4: Launch Application

```bash
cd /Users/rgpro/Desktop/GECKO

# Launch (uses existing fly.toml)
fly launch --no-deploy

# When prompted:
# - App name: gecko-web (or choose your own)
# - Region: iad (US East) or choose closer to you
# - Would you like to set up a Postgresql database? No
# - Would you like to set up an Upstash Redis database? No
# - Would you like to deploy now? No (we need to create volume first)
```

---

## Step 5: Create Persistent Volume

**Important:** Create the volume BEFORE first deploy, otherwise job data won't persist.

```bash
# Create 10GB volume for job data
fly volumes create gecko_data --region iad --size 10

# Verify
fly volumes list
```

---

## Step 6: First Deployment

```bash
# Deploy (this will take 45-90 minutes for Fortran compilation)
fly deploy

# Watch the build logs
fly logs
```

**Note:** The first build is slow because it compiles:
- GECKO-A Fortran source (79 object files)
- BOXMODEL4GECKO (kinetic solver)

Subsequent deploys use cached layers and take ~5 minutes.

---

## Step 7: Verify Deployment

```bash
# Check status
fly status

# Open in browser
fly open

# View logs
fly logs --no-tail
```

Your app will be live at: `https://gecko-web.fly.dev`

---

## Step 8: Connect GitHub for Auto-Deploy

### Option A: GitHub Actions (Recommended)

1. **Get Fly.io API Token:**
   ```bash
   fly tokens create deploy -x 999999h
   ```
   Copy the token.

2. **Add Token to GitHub Secrets:**
   - Go to: `https://github.com/YOUR_USERNAME/gecko-web/settings/secrets/actions`
   - Click "New repository secret"
   - Name: `FLY_API_TOKEN`
   - Value: (paste your token)

3. **Push to Deploy:**
   Now every push to `main` branch will auto-deploy:
   ```bash
   git add .
   git commit -m "Update feature"
   git push origin main
   # GitHub Actions will deploy to Fly.io automatically
   ```

### Option B: Fly.io Dashboard Link

1. Go to: https://fly.io/dashboard
2. Select your app (gecko-web)
3. Go to Settings → Source
4. Connect GitHub repository
5. Select branch: `main`
6. Enable auto-deploy

---

## Configuration Details

### Current Settings (fly.toml)

| Setting | Value | Purpose |
|---------|-------|---------|
| Region | `iad` (US East) | Change based on your users' location |
| Memory | 4096 MB | Required for complex VOC simulations |
| CPUs | 2 shared | Adequate for typical usage |
| Min machines | 1 | Avoids cold starts |
| Volume | 10 GB | Stores simulation results |

### Change Region

Edit `fly.toml` and change `primary_region`:
- `iad` - US East (Virginia)
- `lax` - US West (Los Angeles)
- `lhr` - Europe (London)
- `fra` - Europe (Frankfurt)
- `sin` - Asia (Singapore)
- `syd` - Australia (Sydney)

```bash
# After changing fly.toml
fly deploy
```

### Scale Resources

```bash
# Increase memory (if complex VOCs fail)
fly scale memory 8192

# Add more CPU
fly scale vm shared-cpu-4x

# Check current scale
fly scale show
```

---

## Costs Breakdown

| Resource | Amount | Monthly Cost |
|----------|--------|--------------|
| VM (shared-cpu-2x, 4GB) | 1 | ~$12 |
| Persistent Volume (10GB) | 1 | ~$1.50 |
| Bandwidth | ~10 GB | Free (100GB included) |
| **Total** | | **~$13.50/month** |

### Reduce Costs

To reduce costs when not actively using:

```bash
# Stop the app (no charges while stopped)
fly scale count 0

# Restart when needed
fly scale count 1
```

Or enable auto-stop (saves ~50% but has cold starts):
```toml
# In fly.toml
[http_service]
  auto_stop_machines = true
  min_machines_running = 0
```

---

## Troubleshooting

### Build Fails / Timeout

```bash
# Check build logs
fly logs

# If Fortran compilation fails, ensure Dockerfile has:
# - gfortran installed
# - libnetcdf-dev, libnetcdff-dev installed
```

### Out of Memory

```bash
# Increase memory
fly scale memory 8192

# Check memory usage
fly status
```

### Volume Not Mounted

```bash
# Verify volume exists
fly volumes list

# Should show gecko_data attached to your app
```

### App Won't Start

```bash
# Check logs
fly logs

# SSH into container for debugging
fly ssh console

# Inside container:
ls -la /app
ls -la /data
python -c "import gecko_web"
```

### DNS Not Working

It can take 5-10 minutes for DNS to propagate after first deploy:
```bash
# Check if app is running
fly status

# Try direct IP access from status output
```

---

## Useful Commands

```bash
# Deploy
fly deploy

# View logs (live)
fly logs

# View logs (recent)
fly logs --no-tail

# Check status
fly status

# Open in browser
fly open

# SSH into container
fly ssh console

# Scale up/down
fly scale count 1
fly scale count 0

# Restart
fly apps restart gecko-web

# Delete app (careful!)
fly apps destroy gecko-web

# List all apps
fly apps list

# Check billing
fly billing
```

---

## Automatic Deployment Flow

```
GitHub Push (main branch)
    ↓
GitHub Actions triggered (.github/workflows/fly-deploy.yml)
    ↓
Fly.io CLI runs: flyctl deploy --remote-only
    ↓
Fly.io builds Docker image (uses cache if available)
    ↓
New version deployed with zero-downtime
    ↓
Health check passes (/api/environment)
    ↓
Live at: https://gecko-web.fly.dev
```

---

## Links

- **Your App Dashboard:** https://fly.io/apps/gecko-web
- **Fly.io Documentation:** https://fly.io/docs/
- **Fly.io Status:** https://status.fly.io/
- **Billing:** https://fly.io/dashboard/billing

---

## Support

- **Fly.io Community:** https://community.fly.io/
- **GECKO Issues:** https://github.com/YOUR_USERNAME/gecko-web/issues

---

*Last updated: 2026-01-28*
