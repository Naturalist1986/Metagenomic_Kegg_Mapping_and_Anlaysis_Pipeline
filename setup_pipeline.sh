#!/bin/bash

################################################################################
# Pipeline Setup Script
#
# Run this once before submitting the pipeline to create necessary directories
#
# Usage: ./setup_pipeline.sh
################################################################################

echo "Setting up Metagenomics Pipeline..."

# Create logs directory
mkdir -p logs
echo "✓ Created logs/ directory"

# Create scripts directory if it doesn't exist
mkdir -p scripts
echo "✓ Verified scripts/ directory"

# Make all scripts executable
chmod +x run_metagenomics_pipeline.sh
chmod +x scripts/*.sh
chmod +x *.py
echo "✓ Made scripts executable"

echo ""
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Copy and edit config file:"
echo "     cp config.example.sh config.sh"
echo "     nano config.sh"
echo ""
echo "  2. Submit pipeline:"
echo "     sbatch run_metagenomics_pipeline.sh config.sh"
echo ""
echo "=========================================="
