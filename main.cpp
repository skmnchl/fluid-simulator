//
// Setup:
//     Install emscripten: https://emscripten.org/docs/getting_started/downloads.html
//
// Build:
//     WASM : emcc -pthread main.cpp FlipFluid.cpp FluidRenderer.cpp -sINITIAL_MEMORY=64MB -sPTHREAD_POOL_SIZE=5 -sUSE_SDL=2 -sFULL_ES2=1 -o index.html --shell-file shell.html
//     Console : g++ main.cpp FlipFluid.cpp -o main
//
// Run:
//     WASM: Run a http server serving index.html with proper COOP and COEP headers (Required for multi threading)
//     or just run the command below and check localhost:8888
//         emrun --no_browser --port 8888 --serve_after_close index.html
//     Console:
//         ./main
//

#include <chrono>

#include "FlipFluid.h"

#ifdef __EMSCRIPTEN__
#include "FluidRenderer.h"
#endif

int main() {
    // faster I/O
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    // use small grid but more particles instead of large grid for better performance
    int xSize = 25, ySize = 25, zSize = 40;

    // make basic box type grid
    std::vector<std::vector<std::vector<Grid::CellType>>> initCell;
    initCell.resize(
        xSize, std::vector<std::vector<Grid::CellType>>(
        ySize, std::vector<Grid::CellType>(
        zSize, Grid::CellType::AIR
    )));

    for(int i=0;i<xSize;i++) {
        for(int j=0;j<ySize;j++) {
            for(int k=0;k<zSize;k++) {
                if (i==0 || i==xSize-1
                 || j==0 || j==ySize-1
                 || k==0 || k==zSize-1) {
                    initCell[i][j][k] = Grid::CellType::SOLID;
                }
            }
        }
    }
    Grid grid = Grid(initCell);

    // make half filled atoms
    std::vector<float> pos;
    std::vector<float> vel;
    float spacing = 0.6;
    for(float i=2;i<xSize-2;i+=spacing) {
        for(float j=ySize/3;j<ySize-2;j+=spacing) {
            for(float k=zSize*2/3;k<zSize-2;k+=spacing) {
                float x = i;
                float y = j;
                float z = k;
                float u = 0;
                float v = 0;
                float w = 0;
                pos.push_back(x);
                pos.push_back(y);
                pos.push_back(z);
                vel.push_back(u);
                vel.push_back(v);
                vel.push_back(w);
            }
        }
    }
    Atom atom = Atom(pos, vel);

    // other configs
    float gravity = 50;
    int simulateStep = 5;
    float dt = 0.02;
    float minDist = 0.5;
    int pushApartIteration = 0;
    int incompressIteration = 5;
    float overRelaxation = 1.9;
    float stiff = 2.0;
    float restDensity = 3.0;
    int numThreads = 0; // single threading for github pages deploy

    // make fluid
    FlipFluid fluid = FlipFluid(
        grid,
        atom,
        gravity,
        simulateStep,
        dt,
        minDist,
        pushApartIteration,
        incompressIteration,
        overRelaxation,
        stiff,
        restDensity,
        numThreads
    );

#ifdef __EMSCRIPTEN__
    FluidRenderer renderer = FluidRenderer(&fluid, 512, 512);

    emscripten_set_main_loop_arg(FluidRendererMainLoop, &renderer, 0, true);
#else
    while(true) {
        auto start = std::chrono::high_resolution_clock::now(); // start time
    
        fluid.update();
        fluid.printFluid(100, 60, 50, 40, 0.4);

        auto end = std::chrono::high_resolution_clock::now(); // end time
        std::chrono::duration<double> duration = end - start; // compute duration
        std::cout << "Time elapsed: " << duration.count() << '\n';
    }
#endif
}