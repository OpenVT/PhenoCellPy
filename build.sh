set -e
cp LICENSE .package/LICENSE
cp README.md .package/README.md
cp pyproject.toml .package/pyproject.toml

# Clean and symlink source dirs
rm -rf .package/phenocellpy .package/CC3D_examples .package/DeveloperTestZone .package/TF_examples .package/paper
cp -r phenocellpy .package/
cp -r CC3D_examples .package/phenocellpy
cp -r DeveloperTestZone .package/phenocellpy
cp -r TF_examples .package/phenocellpy
# Do not copy over the paper sub dir

python -m build .package

rm -rf .package/LICENSE .package/README.md .package/pyproject.toml \
       .package/phenocellpy

rm -rf .dist
mv .package/dist .dist
