#!/bin/bash
# Dengue Variant Tracker - Data Download Script
# Downloads dengue virus sequences from NCBI Virus Database
# Optimized for low RAM environments (limited to 50 sequences)

set -e  # Exit on error

echo "=========================================="
echo "Dengue Virus Data Download Script"
echo "=========================================="
echo ""

# Create data directories
echo "[1/5] Creating data directories..."
mkdir -p data/raw
mkdir -p data/processed
mkdir -p plots

# Check for internet connectivity
echo "[2/5] Checking internet connection..."
if ! ping -c 1 google.com &> /dev/null; then
    echo "ERROR: No internet connection detected."
    echo "Please check your network and try again."
    exit 1
fi
echo "✓ Internet connection confirmed"

# Download dengue sequences from NCBI
echo "[3/5] Downloading dengue virus sequences from NCBI..."
echo "Note: This downloads a small subset (50 sequences) optimized for 4GB RAM systems"

# Method 1: Using NCBI Datasets (recommended if available)
# Uncomment if you have ncbi-datasets-cli installed
# datasets download virus genome taxon 12637 --include genome --limit 50 --filename data/raw/dengue_ncbi.zip
# unzip -o data/raw/dengue_ncbi.zip -d data/raw/

# Method 2: Direct download via wget (backup method)
# Note: NCBI's direct download URLs may change. This is a fallback approach.
# For production use, consider using NCBI E-utilities or Datasets API

# Alternative: Using a pre-selected small dataset URL
# This is a sample URL - you may need to update based on NCBI's current API
NCBI_URL="https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/Database/nph-select.cgi"

echo "Attempting to download dengue sequences..."
wget -O data/raw/dengue_sequences.fasta \
    "${NCBI_URL}?cmd=download&viruses=12637&flt=1&fmt=1&origin=vvr" \
    --post-data="seqnum=50" \
    --timeout=60 \
    --tries=3 \
    2>&1 || {
    echo ""
    echo "WARNING: Automated download failed. This is common with NCBI's dynamic URLs."
    echo ""
    echo "MANUAL DOWNLOAD INSTRUCTIONS:"
    echo "1. Visit: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Dengue%20virus,%20taxid:12637"
    echo "2. Click 'Download' button"
    echo "3. Select 'Nucleotide' and 'FASTA' format"
    echo "4. Limit to 50 sequences"
    echo "5. Save file as: data/raw/dengue_sequences.fasta"
    echo ""
    echo "OR use the sample test data provided in the next step..."
    echo ""
}

# Check if download was successful
if [ -f "data/raw/dengue_sequences.fasta" ] && [ -s "data/raw/dengue_sequences.fasta" ]; then
    echo "✓ Download successful"
    FILE_SIZE=$(wc -c < "data/raw/dengue_sequences.fasta")
    echo "  File size: $(numfmt --to=iec-i --suffix=B $FILE_SIZE)"
else
    echo "[4/5] Creating sample test dataset..."
    # Create a small test dataset with sample dengue sequences
    # This allows the project to run even without successful NCBI download
    cat > data/raw/dengue_sequences.fasta << 'EOF'
>NC_001477.1 Dengue virus 1, complete genome
AGTTGTTAGTCTACGTGGACCGACAAGAACAGTTTCGAATCGGAAGCTTGCTTAACGTAGTTCTAACAGT
TTTTATTTAGAGAGCAGATCTCTGATGAATAACCAACGGAAAAAGGCGAGAAATACACGCTTTTCAATATA
TGCGAAAAAAGCAGGAAACACATGGATAAAACAGGCTGAGAAGACCATCAGCCAGAAACCAATTCTCACTG
AAAGGGCTGAAACACGCGGCCAATAATCTAAAGCAATACATCTCAGTTGTGCTAGGCATTATTGGGAAGGA
GTTGGAAACACATCTTGGAGAAATGGAAAGGTTAAGGAAACACAAATTTTCAATTCTGCAGACATGGACCG
ATGAAACACATGGAGAGCTTCACATTCTCACCAAGAGCCGGAGGTGGGTGGACGCCAGAGTCATACCCTCC
ACCACAGAGAGAAGAGAGCCCGAGCAGCAAGACAAACCTGGAGGGTGAGGCTCAACTGGTATCTTGGGACA
CGAGGAGCAATCTGACCAATGGGTGGGTGATGGAGATCTATTTCAATGGGGGGCCACCAACATTGCTGGTT
GTAGGAGTTGACTCTGAAGGAACAATAGGGGAAGGAGTCGGAATCATTTTACAACAGGAGCCGGATTCACC
CATCTATGGTGGTTGGGAGGAGAACAGGGAGATTCTTCTATCCAAGGAAAAAGGAAGTGAAATAACAAGAG
GAGAGGACTGGGCAGGATGGCTCACAGCCATTTAAGCCAGCATTTGTGGGGCTGGGCAGGGAAACATGTCG
TTGGGTGGACAGGCAGGCAAGCCAGCTTCAGGGAATTTGATCAATAATAAGGATGACATTGAAAGTCAAGA
AGCACGGGAACTTGGTGGTTGTCATAACTGCTGCCTCTTTTCTGGCTCTCATGGTAGGTCAGGATGGATGT
>NC_001474.2 Dengue virus 2, complete genome
AGTTGTTAGTCTACGTGGACCGACAAGATACAGTTCGATTCCGGAGGCTGCGTCACGTAGTTCTAACAGTT
TTTATTGAAGAGCAGACATCACAATGGATAAAACAGAGTGATACGACCACACCACTGAAAGAGAAAAGAAA
ACTGAAGATGGACCCACTGGAAACAGGCTGAATGGCCAGGCAACAGCACCATTCATGGATATTCCCGAGGA
TCTGGCCAACATGGCTCAGCGGATGGTTCTCATACAGACCTCGGGCCATCCGCAGAGGCCGCACGCACAAG
GCAAGATGAACATCGACGTGAAACCAACGTCAACAATCGTTGAAGCGCCAAAGGGACAAAGCTTCACGTGG
AAGTCCATCACAACTCCATGGGACTTCACGATACCGGTCACTGGAAAGTAGGACTAGTGTTCGCAGGATCC
CTGCTAATACAAAGCAAGACATGGAAAGTGGTTCCTGGAAAGATGGAGACATCGGTCACCAGCCAGGATCA
GGGAAGGAGATTGATGGGGAAGCTGTTGGAGGAATGATTTGGAAAGGGAGGAGTTGTGACCTGTGCATATA
>NC_001475.2 Dengue virus 3, complete genome
AGTTGTTACTCTACGTGGACCGACAAGAACAGTTTCGAATCGGAAGCTTGCTTAACGTAGTTCTAACAGTT
TTTATTAGAGAGCAGGTCTCTGGAATTAACAGGAAAGAAGGCTGAAACGTAGTTTTCCATTCACCATCTGC
GGGGCTTGCAATAAGCATCTTCTAGGGATTTCCATGGATGAAACACATGGAGAAATTCACATTCTCACCAA
GAGCTGGGGCTGGGTGGACGCCAGAGTTATACCAATCACCACAGAGAGAAGAGAACCTGAGAAGCATGACA
AAAATGGAAGGGTGAGGCTTCATTGGTATCTTGGGGCCAGAAGAACAATCAGTCAATGGGTGGGTGATGGA
GATCTACTTCAATGGGGGACCACCAACATTGCTGGTTCTAGGAGTTGACTCTGAAGGAACAATAGGAGAAG
GAGTCGGAGTCATTTTACAGCAGGACCCAGATTCACCCATCTATGGTGGTTGGGAGGAGAACAGGGAGATT
EOF
    
    echo "✓ Sample test dataset created (3 dengue serotypes)"
    echo "  Note: For full analysis, please download real data from NCBI manually"
fi

# Create a smaller test subset for initial testing (saves RAM)
echo "[5/5] Creating test subset for low-RAM testing..."
head -n 40 data/raw/dengue_sequences.fasta > data/raw/test_dengue.fasta
echo "✓ Test subset created: data/raw/test_dengue.fasta"

# Summary
echo ""
echo "=========================================="
echo "Download Complete!"
echo "=========================================="
echo "Files created:"
echo "  - data/raw/dengue_sequences.fasta (full dataset)"
echo "  - data/raw/test_dengue.fasta (test subset)"
echo ""
echo "Next steps:"
echo "  1. Run quality control: Rscript qc_analysis.R"
echo "  2. Launch dashboard: Rscript -e \"shiny::runApp('app.R')\""
echo ""
echo "For RAM-constrained systems, start with test_dengue.fasta"
echo "=========================================="
