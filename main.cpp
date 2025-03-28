//
// Setup:
//     Install emscripten: https://emscripten.org/docs/getting_started/downloads.html
//
// Build:
//     emcc -pthread --bind main.cpp flip_fluid.cpp -sINITIAL_MEMORY=64MB -sPTHREAD_POOL_SIZE=5 -o main.js
//
// Run:
// Run a http server serving index.html with proper COOP and COEP headers (Required for multi threading)
// or just run the command below and check localhost:8888
//     emrun --no_browser --port 8888 --serve_after_close index.html
//

#include <chrono>

#ifdef __EMSCRIPTEN__
#include <cstdio> // usde for hex color manipulation
#include <emscripten.h>
#include <emscripten/val.h>
#endif

#include "flip_fluid.h"

// using global variable for simple implementation
FlipFluid *fluid = nullptr;

#ifdef __EMSCRIPTEN__

const float scale = 12.0;
const float offsetX = 600;
const float offsetY = 300;

// update fluid and draw it to "canvas" every frame
void animateOnHtmlCanvas() {
    const auto updateStart = std::chrono::high_resolution_clock::now(); // start time
    // update fluid status
    fluid->update();
    const auto updateEnd = std::chrono::high_resolution_clock::now(); // end time

    const auto drawStart = std::chrono::high_resolution_clock::now(); // start time
    emscripten::val document = emscripten::val::global("document");
    emscripten::val canvas = document.call<emscripten::val>("getElementById", emscripten::val("canvas"));
    emscripten::val ctx = canvas.call<emscripten::val>("getContext", emscripten::val("2d"));

    // clear canvas
    ctx.call<void>("clearRect", 0, 0, canvas["width"].as<int>(), canvas["height"].as<int>());

    for(int a=0;a<fluid->atom.cnt;a++) {
        // set color or atom
        int r = 128 + fluid->atom.vel[a*3 + 0]*10;
        int g = 128 + fluid->atom.vel[a*3 + 1]*10;
        int b = 128 + fluid->atom.vel[a*3 + 2]*10;
        r = (r>255) ? 255 : r; r = (r<0) ? 0 : r;
        g = (g>255) ? 255 : g; g = (g<0) ? 0 : g;
        b = (b>255) ? 255 : b; b = (b<0) ? 0 : b;

        char hexColor[10];
        std::sprintf(hexColor, "#%02X%02X%02X80", r, g, b);
        ctx.set("fillStyle", emscripten::val(std::string(hexColor)));

        // draw atom
        float i = fluid->atom.pos[a*3 + 0];
        float j = fluid->atom.pos[a*3 + 1];
        float k = fluid->atom.pos[a*3 + 2];                

        float y = -0.5*(i+k) + j;
        float x = -0.86602540378*(k-i); // sqrt(3)/2 * (k-i)
        y = offsetY + y*scale;
        x = offsetX + x*scale;
        ctx.call<void>("fillRect", x, y, 5, 5);
    }

    auto drawEnd = std::chrono::high_resolution_clock::now(); // end time

    // compute duration
    std::chrono::duration<double> updateDuration = updateEnd - updateStart; 
    std::chrono::duration<double> drawDuration = drawEnd - drawStart;
    
    // write text
    ctx.call<void>("fillText", emscripten::val("Flip Fluid simulator: "), 30, 30);
    ctx.call<void>("fillText", emscripten::val("Update duration: "), 30, 50);
    ctx.call<void>("fillText", emscripten::val(std::to_string(updateDuration.count())), 110, 50);
    ctx.call<void>("fillText", emscripten::val("Draw duration: "), 30, 70);
    ctx.call<void>("fillText", emscripten::val(std::to_string(drawDuration.count())), 110, 70);
}
#endif

int main() {
    // faster I/O
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    int xSize = 60, ySize = 40, zSize = 100;
    // int xSize = 30, ySize = 20, zSize = 50;

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
    for(int i=2;i<xSize*3/4;i++) {
        for(int j=2;j<ySize*3/4;j++) {
            for(int k=2;k<zSize/2;k++) {
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
    float dt = 0.03;
    float minDist = 0.5;
    int pushApartIteration = 0;
    int incompressIteration = 5;
    float overRelaxation = 1.9;
    float stiff = 2.0;
    float restDensity = 3.0;
    int numThreads = 5;

    // make fluid
    fluid = new FlipFluid(
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
    emscripten_set_main_loop(animateOnHtmlCanvas, 0, 0);
#else
    while(true) {
        auto start = std::chrono::high_resolution_clock::now(); // start time
    
        fluid->update();
        fluid->printFluid(100, 60, 50, 40, 0.4);

        auto end = std::chrono::high_resolution_clock::now(); // end time
        std::chrono::duration<double> duration = end - start; // compute duration
        std::cout << "Time elapsed: " << duration.count() << '\n';
    }
#endif
}