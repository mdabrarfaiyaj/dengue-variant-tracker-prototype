# Dengue Variant Tracker Dashboard

## 🦟 Project Overview
An interactive R Shiny dashboard for tracking dengue virus variants using Bioconductor tools. This project analyzes public dengue virus sequences from NCBI to identify mutation patterns and motifs relevant to public health in tropical regions like Bangladesh.

## 🎯 Objectives
- Fetch and process public dengue virus sequences
- Perform quality control on genomic data
- Identify mutation patterns and motifs in viral genomes
- Visualize variant patterns through an interactive dashboard
- Demonstrate end-to-end bioinformatics workflow

## 🔬 Real-World Impact
This dashboard addresses practical public health needs:
- Monitoring vaccine escape mutations
- Tracking dengue serotype evolution
- Supporting local health surveillance efforts in Dhaka and Bangladesh

## 📊 Data Source
- **Source**: NCBI Virus Database (public, open-access data)
- **Virus**: Dengue virus (DENV) - all serotypes
- **Data Type**: Nucleotide sequences in FASTA format
- **Dataset Size**: ~10-50 sequences (optimized for low RAM environments)
- **Citation**: NCBI Virus Database. https://www.ncbi.nlm.nih.gov/labs/virus/

## 🔧 Technical Stack
- **Language**: R (4.0+)
- **Core Packages**: 
  - Bioconductor: `Biostrings`, `ShortRead`, `BSgenome`
  - Data Processing: `dplyr`, `tidyr`
  - Visualization: `ggplot2`, `plotly`
  - Dashboard: `shiny`, `shinydashboard`
- **Automation**: Bash shell scripting
- **Version Control**: Git/GitHub
- **Deployment**: shinyapps.io (free tier)

## 🚀 Quick Start

### Prerequisites
```r
# Install R packages
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "tidyr", "plotly", "DT"))

# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "ShortRead", "BSgenome"))
```

### Setup & Run
```bash
# 1. Clone repository
git clone <your-repo-url>
cd viral_tracker_dashboard

# 2. Download data (requires internet)
chmod +x download_data.sh
./download_data.sh

# 3. Run quality control and analysis
Rscript qc_analysis.R

# 4. Launch dashboard
Rscript -e "shiny::runApp('app.R')"
```

## 📁 Project Structure
```
viral_tracker_dashboard/
├── README.md                 # Project documentation
├── download_data.sh          # Data acquisition script
├── qc_analysis.R            # Quality control and motif analysis
├── app.R                    # Shiny dashboard application
├── data/                    # Data directory (gitignored)
│   ├── raw/                 # Raw downloaded sequences
│   └── processed/           # Processed analysis results
├── plots/                   # Generated visualizations
├── utils/                   # Helper functions
│   └── analysis_helpers.R   # Reusable analysis functions
└── docs/                    # Additional documentation
    └── methodology.md       # Detailed methods
```

## 🧬 Key Features

### 1. Data Processing Pipeline
- Automated sequence downloading from NCBI
- Quality filtering (sequence length, ambiguous bases)
- Sequence alignment and motif detection

### 2. Analysis Capabilities
- **Motif Scanning**: Identifies known dengue mutation hotspots
- **Variant Statistics**: Calculates mutation frequencies
- **Sequence Composition**: GC content, codon usage
- **Phylogenetic Markers**: Serotype-specific signatures

### 3. Interactive Dashboard
- **Overview Tab**: Dataset summary statistics
- **Motif Explorer**: Search custom or predefined motifs
- **Variant Visualization**: Interactive plots of mutation patterns
- **Data Table**: Browse and filter sequence information

## 🧪 Analysis Workflow

### Step 1: Data Acquisition
```bash
./download_data.sh
# Downloads dengue sequences from NCBI (limited to 50 for low RAM)
```

### Step 2: Quality Control
```r
source("qc_analysis.R")
# - Loads sequences
# - Filters by length (>500 bp)
# - Removes sequences with >5% ambiguous bases
# - Generates QC report
```

### Step 3: Motif Analysis
```r
# Searches for known dengue motifs:
# - ATG: Start codons
# - GAC: Common mutation site in E protein
# - AATAAA: Poly-A signal regions
# - Custom user-defined patterns
```

### Step 4: Dashboard Launch
```r
shiny::runApp("app.R")
# Opens interactive dashboard in browser
```

## 🔒 Ethical Considerations
- **Public Data Only**: Uses exclusively open-access, anonymized viral sequences
- **No Personal Health Information**: Complies with data privacy regulations
- **Proper Attribution**: All data sources are cited
- **Reproducibility**: Complete code and methodology shared openly

## 💾 Low-RAM Optimization
This project is optimized for systems with 4GB RAM:
- Processes sequences in small batches
- Limits dataset to 10-50 sequences
- Uses efficient Bioconductor data structures
- Includes memory cleanup (`gc()`) after intensive operations
- Alternative: Use Google Colab for R sessions (free cloud computing)

## 📈 Results & Outputs
- **motif_matches.csv**: Detected motif positions and frequencies
- **qc_summary.csv**: Quality control metrics
- **variant_plots.png**: Static visualizations for reports
- **Interactive Dashboard**: Real-time exploration via Shiny

## 🌐 Deployment
```r
# Deploy to shinyapps.io (requires free account)
library(rsconnect)
rsconnect::setAccountInfo(name="<ACCOUNT>", token="<TOKEN>", secret="<SECRET>")
rsconnect::deployApp()
```

## 🤝 Contributing
Suggestions and improvements welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with clear description

## 📚 References & Resources
- NCBI Virus Database: https://www.ncbi.nlm.nih.gov/labs/virus/
- Bioconductor: https://bioconductor.org/
- Dengue WHO Fact Sheet: https://www.who.int/news-room/fact-sheets/detail/dengue-and-severe-dengue
- R Shiny: https://shiny.posit.co/

## 📝 License
MIT License - Free to use with attribution

## 👤 Author
**Your Name**
- LinkedIn: [https://www.linkedin.com/in/md-abrar-faiyaj-559246381/]
- GitHub: [https://github.com/mdabrarfaiyaj]
- Email: faiyaj.mdabrar@gmail.com

## 🏆 Skills Demonstrated
- Bioinformatics data analysis (Bioconductor/R)
- Shell scripting automation
- Interactive data visualization (Shiny)
- Version control (Git)
- Public health data interpretation
- Low-resource computing optimization

---
*This project was developed as part of a bioinformatics portfolio showcasing real-world genomic data analysis skills applicable to tropical disease surveillance.*
# Dengue Variant Tracker
# Dengue Variant Tracker
