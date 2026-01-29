# Quick Start Guide - Dengue Variant Tracker Dashboard

## 🚀 Complete Setup in 15 Minutes

This guide will help you get the Dengue Variant Tracker Dashboard running on your system, even with limited RAM (4GB).

---

## Prerequisites

### Required Software
1. **R** (version 4.0 or higher)
   - Download: https://cran.r-project.org/
   - Verify installation: Open terminal/command prompt, type `R --version`

2. **RStudio** (optional but recommended)
   - Download: https://posit.co/download/rstudio-desktop/

3. **Git** (for version control)
   - Download: https://git-scm.com/downloads
   - Verify: `git --version`

4. **Internet connection** (for downloading sequences and packages)

---

## Installation Steps

### Step 1: Clone or Download Repository
```bash
# If using Git:
git clone <your-repo-url>
cd viral_tracker_dashboard

# Or download ZIP and extract, then navigate to folder
```

### Step 2: Install R Packages (5-10 minutes)
```bash
# Make installation script executable (Mac/Linux)
chmod +x install_packages.R

# Run installation
Rscript install_packages.R
```

**OR in R/RStudio console:**
```r
source("install_packages.R")
```

**What this installs:**
- CRAN packages: shiny, shinydashboard, ggplot2, dplyr, tidyr, plotly, DT
- Bioconductor: Biostrings, ShortRead, BSgenome

**Troubleshooting:**
- If package installation fails, try installing individually:
  ```r
  install.packages("shiny")
  # Repeat for each package
  ```
- For Bioconductor packages:
  ```r
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("Biostrings")
  ```

### Step 3: Download Dengue Data (2-3 minutes)
```bash
# Make script executable (Mac/Linux)
chmod +x download_data.sh

# Run download
./download_data.sh
```

**Windows users:** Right-click `download_data.sh` and open with Git Bash, or:
```bash
bash download_data.sh
```

**What this does:**
- Creates data folders
- Downloads dengue virus sequences from NCBI (or creates test dataset)
- Limits to 50 sequences for low-RAM systems

**If download fails:**
The script includes sample test data. Alternatively:
1. Visit: https://www.ncbi.nlm.nih.gov/labs/virus/
2. Search: "Dengue virus"
3. Download 10-50 sequences as FASTA
4. Save to: `data/raw/dengue_sequences.fasta`

### Step 4: Run Analysis (2-3 minutes)
```bash
Rscript qc_analysis.R
```

**OR in R/RStudio:**
```r
source("qc_analysis.R")
```

**What this does:**
- Quality control filtering
- Motif pattern detection
- Generate visualizations
- Create processed datasets

**Output files created:**
- `data/processed/qc_summary.csv` - Quality metrics
- `data/processed/motif_matches.csv` - Detected motifs
- `data/processed/filtered_sequences.fasta` - Clean sequences
- `plots/*.png` - 4 visualization files

### Step 5: Launch Dashboard! 🎉
```bash
Rscript -e "shiny::runApp('app.R')"
```

**OR in R/RStudio:**
```r
shiny::runApp('app.R')
```

**What to expect:**
- Browser window opens automatically
- Dashboard runs on http://127.0.0.1:XXXX
- Explore 6 interactive tabs
- Stop server: Press `Ctrl+C` in terminal or click stop in RStudio

---

## Using the Dashboard

### Overview Tab
- See dataset statistics at a glance
- View sequence length and GC content distributions
- Check top motifs found

### Sequence Quality Tab
- Detailed quality control metrics
- Base composition analysis
- Correlation plots

### Motif Explorer Tab
- Filter by specific motifs
- View match frequencies
- Explore position details

### Custom Analysis Tab
- **Try this:** Search for custom DNA patterns
- Example motifs:
  - `ATG` - Start codon
  - `CACAG` - 5' UTR region
  - `AATAAA` - Poly-A signal
- Adjust mismatch tolerance for fuzzy matching

### Data Table Tab
- Browse complete dataset
- Filter and sort
- Export to CSV

### About Tab
- Project information
- Methodology
- Citations

---

## Low-RAM Tips (4GB Systems)

If you experience lag or crashes:

1. **Use test dataset only:**
   - In `qc_analysis.R`, ensure `USE_TEST_DATA <- TRUE` (line 13)

2. **Close other applications:**
   - Close browsers, editors before running analysis

3. **Process in batches:**
   - Comment out heavy visualizations temporarily

4. **Use Google Colab (free cloud alternative):**
   ```r
   # Upload files to Colab
   # Install packages in Colab environment
   # Run analysis there
   ```

5. **Reduce sequence count:**
   - Edit `download_data.sh`, change `seqnum=50` to `seqnum=10`

---

## Publishing Your Work

### GitHub
```bash
# Initialize repository (if not cloned)
git init
git add .
git commit -m "Initial commit - Dengue Variant Tracker"

# Create repository on GitHub, then:
git remote add origin https://github.com/yourusername/dengue-tracker.git
git push -u origin main
```

### Deploy Dashboard Online (Free!)

**Option 1: shinyapps.io** (Recommended)
1. Sign up: https://www.shinyapps.io/
2. Get token from account settings
3. In R console:
   ```r
   library(rsconnect)
   rsconnect::setAccountInfo(name='yourname', 
                             token='YOUR_TOKEN',
                             secret='YOUR_SECRET')
   rsconnect::deployApp()
   ```
4. Share your live dashboard URL!

**Option 2: Render as static HTML**
```r
# For non-interactive version
rmarkdown::render("analysis_report.Rmd")  # Create Rmd first
```

### LinkedIn Post Template
```
🦟 Excited to share my latest bioinformatics project!

I built an interactive R Shiny dashboard for tracking dengue virus variants 
using Bioconductor tools. The dashboard analyzes public NCBI sequences to 
identify mutation patterns relevant to public health in tropical regions.

🔬 Key features:
- Automated quality control pipeline
- Motif detection for known mutation sites  
- Interactive visualizations
- Custom pattern search

💻 Technologies: R, Bioconductor, Shiny, ggplot2
🌍 Real-world impact: Supports local disease surveillance

🔗 Live Dashboard: [your-shinyapps-url]
🔗 GitHub Code: [your-github-url]

#Bioinformatics #DataScience #PublicHealth #RStats #Genomics #DengueVirus
```

---

## Troubleshooting

### "Package not found" error
```r
# Reinstall missing package
install.packages("package_name")
# Or for Bioconductor:
BiocManager::install("package_name")
```

### "File not found" error
- Check you're in correct directory: `getwd()` in R
- Ensure scripts were run in order (download → analysis → dashboard)

### Dashboard won't load
- Check R console for error messages
- Verify all data files exist in `data/processed/`
- Try clearing browser cache

### Out of memory
- Reduce dataset size (use test_dengue.fasta only)
- Close other programs
- Restart R session: `.rs.restartR()` in RStudio

### Download script fails
- Use manual download instructions in script output
- Check internet connection
- Try sample test data (automatically created)

---

## Next Steps - Enhancing Your Project

1. **Add more analyses:**
   - Phylogenetic tree visualization
   - Mutation rate calculations
   - Serotype classification

2. **Improve visualizations:**
   - Add more plot types
   - Create animated transitions
   - Export high-res figures

3. **Expand dataset:**
   - Include more sequences (when RAM permits)
   - Add other arboviruses (Zika, Chikungunya)
   - Time series analysis

4. **Write documentation:**
   - Create detailed methodology
   - Add code comments
   - Write blog post about findings

5. **Present findings:**
   - Create presentation slides
   - Submit to conferences
   - Share on ResearchGate

---

## Getting Help

- **R Help:** `?function_name` in R console
- **Bioconductor:** https://support.bioconductor.org/
- **Shiny:** https://shiny.posit.co/r/getstarted/
- **This project:** Open issue on GitHub

---

## Success Checklist

- [ ] R and packages installed
- [ ] Data downloaded
- [ ] Analysis completed successfully
- [ ] Dashboard launches in browser
- [ ] All tabs working
- [ ] Custom motif search functional
- [ ] GitHub repository created
- [ ] Dashboard deployed online
- [ ] LinkedIn post published

---

**Congratulations!** You now have a professional bioinformatics portfolio project! 🎉

For questions or improvements, feel free to contribute via GitHub pull requests.
