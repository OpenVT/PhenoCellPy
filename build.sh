set -e
cp LICENSE .package/LICENSE
cp README.md .package/README.md

# Clean and symlink source dirs
rm -rf .package/phenocellpy .package/CC3D_examples .package/DeveloperTestZone .package/TF_examples
cp -r phenocellpy .package/
cp -r CC3D_examples .package/
cp -r DeveloperTestZone .package/
cp -r TF_examples .package/

python -m build .package

rm -rf .package/LICENSE .package/README.md \
       .package/phenocellpy .package/CC3D_examples .package/DeveloperTestZone .package/TF_examples

rm -rf .dist
mv .package/dist .dist
