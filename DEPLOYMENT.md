# Deployment Guide - Publishing Your Dashboard Online
This is the prototype version. See Bangladeshi version: https://github.com/mdabrarfaiyaj/bangladesh-dengue-variant-tracker

## Deploying to shinyapps.io (Free Tier)

shinyapps.io offers free hosting for Shiny applications with these limits:
- 5 applications
- 25 active hours/month
- Perfect for portfolio projects!

---

## Step-by-Step Deployment

### 1. Create shinyapps.io Account

1. Visit https://www.shinyapps.io/
2. Click "Sign Up"
3. Use GitHub, Google, or email to sign up
4. Choose the **Free** plan

### 2. Get Your Deployment Token

1. Log into shinyapps.io
2. Click your name (top-right) → **Tokens**
3. Click **Show** next to your token
4. Copy the displayed code (looks like):
   ```r
   rsconnect::setAccountInfo(name='yourname',
                             token='ABC123...',
                             secret='xyz789...')
   ```

### 3. Install rsconnect Package

In R console:
```r
install.packages('rsconnect')
```

### 4. Configure Deployment

In R console, paste your token code:
```r
rsconnect::setAccountInfo(
  name='yourname',
  token='YOUR_TOKEN_HERE',
  secret='YOUR_SECRET_HERE'
)
```

This only needs to be done once per computer.

### 5. Prepare for Deployment

Before deploying, ensure:
- All data files are in `data/processed/` folder
- Analysis has been run (`qc_analysis.R`)
- Dashboard works locally

**Important:** The free tier has a 1GB size limit. To reduce size:

```r
# Check app size
tools::dirSize(".")  

# If too large, remove unnecessary files:
unlink("plots/*.png")  # Plots regenerate from data
```

### 6. Deploy Your App!

```r
# From your project directory
setwd("path/to/viral_tracker_dashboard")

# Deploy
rsconnect::deployApp(
  appName = "dengue-variant-tracker",  # Custom name (no spaces)
  appTitle = "Dengue Variant Tracker Dashboard",
  account = "yourname"
)
```

**First deployment takes 5-10 minutes** as it installs packages on the server.

### 7. Success! 🎉

Once complete, you'll see:
```
Application successfully deployed to https://yourname.shinyapps.io/dengue-variant-tracker/
```

Your dashboard is now **live and shareable**!

---

## Managing Your Deployed App

### View App Settings
```r
rsconnect::showLogs(appName = "dengue-variant-tracker")
```

### Update Your App
After making changes locally:
```r
rsconnect::deployApp(appName = "dengue-variant-tracker")
```

### View Usage Stats
1. Log into shinyapps.io
2. Click your app name
3. View **Metrics** tab for:
   - Active hours used
   - Number of visits
   - Performance metrics

### Terminate App (to save hours)
```r
rsconnect::terminateApp(appName = "dengue-variant-tracker")
```

Or on website: Apps → (your app) → Settings → **Archive**

---

## Troubleshooting Deployment

### Error: "Deployment failed"
**Check logs:**
```r
rsconnect::showLogs()
```

**Common issues:**
1. **Missing packages:** Ensure all packages are in CRAN/Bioconductor
2. **Large file size:** Remove unnecessary files
3. **File paths:** Use relative paths only (e.g., `data/processed/`)

### Error: "Exceeded size limit"
```r
# Check what's large
file.info(list.files(recursive = TRUE))$size %>% sort(decreasing = TRUE)

# Remove large files not needed for deployment
unlink("data/raw/*.fasta")  # Keep only processed data
```

### Error: "Package installation failed"
Some Bioconductor packages may fail on shinyapps.io. 

**Solution:** Pre-process data locally, deploy only dashboard with results:
1. Run analysis completely locally
2. Keep only processed CSV files
3. Modify `app.R` to load CSVs only (skip Biostrings)

### Dashboard is slow
Free tier has limited resources. Optimize:
```r
# In app.R, add caching
data <- reactiveFileReader(10000, session, "data/processed/qc_summary.csv", read.csv)
```

---

## Alternative Deployment Options

### Option 1: GitHub Pages (Static Version)

For a non-interactive version:

1. **Convert to static HTML:**
   ```r
   # Create R Markdown report
   rmarkdown::render("analysis_report.Rmd", 
                     output_file = "index.html")
   ```

2. **Push to GitHub:**
   ```bash
   git add index.html
   git commit -m "Add static report"
   git push
   ```

3. **Enable GitHub Pages:**
   - Repository → Settings → Pages
   - Source: main branch, / (root)
   - Your site: https://yourusername.github.io/repo-name/

**Pros:** Unlimited usage, free forever  
**Cons:** No interactivity

### Option 2: Render.com (Free Tier)

Similar to shinyapps.io with 750 hours/month:

1. Create `Dockerfile`:
   ```dockerfile
   FROM rocker/shiny:latest
   RUN R -e "install.packages(c('shiny', 'shinydashboard', 'ggplot2', 'dplyr', 'plotly', 'DT'))"
   RUN R -e "BiocManager::install(c('Biostrings', 'ShortRead'))"
   COPY . /srv/shiny-server/
   EXPOSE 3838
   CMD ["/usr/bin/shiny-server"]
   ```

2. Deploy via Render.com interface

### Option 3: Hugging Face Spaces

Good for ML/AI focused projects:
- Visit https://huggingface.co/spaces
- Create new Space with R Shiny template
- Upload your code

---

## Optimization Tips

### Reduce Package Dependencies
```r
# Instead of loading entire library, import specific functions
plotly::ggplotly()  # Instead of library(plotly)
```

### Cache Data Loading
```r
# In app.R
data <- reactive({
  # Only reload if file changes
  invalidateLater(60000, session)  # Check every minute
  read.csv("data/processed/qc_summary.csv")
})
```

### Minimize File Size
```r
# Compress processed data
saveRDS(data, "data/processed/data.rds", compress = TRUE)

# Load compressed
data <- readRDS("data/processed/data.rds")
```

---

## Monitoring & Maintenance

### Check Active Hours
Free tier = 25 hours/month

**Tips to conserve:**
- Archive app when not showcasing
- Set app to sleep after 15 min idle
- Deploy only during active job search

### Update Regularly
```r
# When you improve the project
rsconnect::deployApp(
  appName = "dengue-variant-tracker",
  forceUpdate = TRUE
)
```

### Backup Deployment Config
```bash
# Save your deployment settings
cp -r rsconnect/ rsconnect_backup/
```

---

## Sharing Your Dashboard

### Professional Links

**LinkedIn:**
```
🔗 Live Dashboard: https://yourname.shinyapps.io/dengue-variant-tracker/
Built with R, Shiny, and Bioconductor. Analyzes dengue virus genomic variants.
```

**GitHub README:**
```markdown
## 🌐 Live Demo
👉 [View Interactive Dashboard](https://yourname.shinyapps.io/dengue-variant-tracker/)
```

**Resume/CV:**
```
Portfolio Projects:
• Dengue Variant Tracker: Interactive R Shiny dashboard for viral genomics
  Live: yourname.shinyapps.io/dengue-variant-tracker
  Code: github.com/yourname/dengue-tracker
```

### Generate QR Code

Create QR code for your dashboard URL:
- Visit https://www.qr-code-generator.com/
- Enter your shinyapps.io URL
- Download PNG
- Add to presentation slides!

---

## Analytics (Optional)

Track usage with Google Analytics:

1. Get GA tracking ID
2. Add to `app.R` in UI section:
   ```r
   tags$head(
     tags$script(async = NA, 
                src = "https://www.googletagmanager.com/gtag/js?id=GA_ID"),
     tags$script(HTML("
       window.dataLayer = window.dataLayer || [];
       function gtag(){dataLayer.push(arguments);}
       gtag('js', new Date());
       gtag('config', 'GA_ID');
     "))
   )
   ```

---

## Checklist for Deployment

- [ ] shinyapps.io account created
- [ ] Token configured locally
- [ ] rsconnect package installed
- [ ] App tested locally
- [ ] Unnecessary files removed
- [ ] App size < 1GB
- [ ] Deployment successful
- [ ] App accessible via URL
- [ ] No errors in logs
- [ ] Shared on LinkedIn/GitHub
- [ ] QR code created (optional)

---

## Support Resources

- **shinyapps.io Docs:** https://docs.posit.co/shinyapps.io/
- **rsconnect Package:** https://github.com/rstudio/rsconnect
- **Community Forum:** https://community.rstudio.com/

---

**Your dashboard is now live!** Share it with the world! 🚀
