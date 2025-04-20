# Fluid Simulator
FLIP fluid simulator built with pure C++, featuring real-time particle-based fluid dynamics and rendering in the browser using WebGL-based WASM.

# Demo
You can see live demo at https://skmnchl.github.io/fluid-simulator/

# Setup
Install emscripten: https://emscripten.org/docs/getting_started/downloads.html

# Build
You can build this project either using wasm or just by console.
This project support multi threading, but to deploy as static page, using only single thread is default when conpiling with wasm.

WASM(multi threading)
```bash
emcc -pthread main.cpp FlipFluid.cpp FluidRenderer.cpp -sINITIAL_MEMORY=64MB -sPTHREAD_POOL_SIZE=5 -sUSE_SDL=2 -sFULL_ES2=1 -o index.html --shell-file shell.html
```

WASM(single threading) : 
```bash
emcc main.cpp FlipFluid.cpp FluidRenderer.cpp -sINITIAL_MEMORY=64MB -sUSE_SDL=2 -sFULL_ES2=1 -o index.html --shell-file shell.html
```

Console : 
```bash
g++ main.cpp FlipFluid.cpp -o main
```

# Run
Run a http server serving index.html with proper COOP and COEP headers (Required for multi threading) or just run the command below and check localhost:8888

WASM: 
```bash
emrun --no_browser --port 8888 --serve_after_close index.html
```

Console:
```bash
./main
```