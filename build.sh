#!/bin/bash
set -e  # exit on first error

# Copy required files into .package
cp LICENSE .package/LICENSE
cp README.md .package/README.md

# Run build
python -m build .package

# Clean up temporary files
rm .package/LICENSE
rm .package/README.md

# Move output to .dist (remove if exists)
rm -rf .dist
mv .package/dist .dist
rm -rf .package/dist
