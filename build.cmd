@echo off
copy LICENSE .package\LICENSE >nul
copy README.md .package\README.md >nul

rmdir /s /q .package\phenocellpy .package\CC3D_examples .package\DeveloperTestZone .package\TF_examples 2>nul
xcopy /e /i phenocellpy .package\phenocellpy
xcopy /e /i CC3D_examples .package\CC3D_examples
xcopy /e /i DeveloperTestZone .package\DeveloperTestZone
xcopy /e /i TF_examples .package\TF_examples

python -m build .package
if errorlevel 1 exit /b %errorlevel%

del .package\LICENSE
del .package\README.md
rmdir /s /q .package\phenocellpy .package\CC3D_examples .package\DeveloperTestZone .package\TF_examples

rmdir /s /q .dist 2>nul
move .package\dist .dist
