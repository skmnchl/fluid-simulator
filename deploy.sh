#!/bin/bash

# Build project
emcc main.cpp FlipFluid.cpp FluidRenderer.cpp MarchingCubes.cpp -sINITIAL_MEMORY=64MB -sUSE_SDL=2 -sFULL_ES2=1 -o index.html --shell-file shell.html

# Stash changes
git stash push --include-untracked

# Checkout to deploy branch
git checkout gh-pages
rm -rf *

# Pop stashed files
git stash pop

# Commit changes
git add .
git commit -m "Deploy update"

# Deploy
# git push origin gh-pages

# Checkout back to main
git checkout main