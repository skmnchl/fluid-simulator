#include <chrono>

#include "flip_fluid.h"

int main() {
    // faster I/O
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    int xSize = 60, ySize = 40, zSize = 100;

    // make basic box type grid
    std::vector<std::vector<std::vector<Grid::CellType>>> init_cell;
    init_cell.resize(
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
                    init_cell[i][j][k] = Grid::CellType::SOLID;
                }
            }
        }
    }
    Grid grid = Grid(init_cell);

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

    while(true) {
        auto start = std::chrono::high_resolution_clock::now(); // start time
    
        fluid.update();
        fluid.printFluid(100, 60, 50, 40, 0.4);

        auto end = std::chrono::high_resolution_clock::now(); // end time
        std::chrono::duration<double> duration = end - start; // compute duration
        std::cout << "Time elapsed: " << duration.count() << '\n';
    }
}