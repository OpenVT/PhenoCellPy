@echo off

copy LICENSE .package\LICENSE >nul
copy README.md .package\README.md >nul

python -m build .package
if errorlevel 1 exit /b %errorlevel%

del .package\LICENSE
del .package\README.md

:: Move output to .dist
rmdir /s /q .dist 2>nul
move .package\dist .dist
rmdir /s /q .package\dist 2>nul
