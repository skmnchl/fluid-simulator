#!/bin/bash

# Build project
emcc main.cpp FlipFluid.cpp FluidRenderer.cpp -sINITIAL_MEMORY=64MB -sUSE_SDL=2 -sFULL_ES2=1 -o dist/index.html --shell-file shell.html

# Commit and push to gh-pages
git checkout gh-pages
rm -rf *
cp -r dist/* .
git add .
git commit -m "Deploy update"
git push origin gh-pages
git checkout main