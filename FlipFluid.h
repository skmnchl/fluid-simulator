#pragma once

#include <iostream>
#include <vector>
#include <cassert>

#include <thread>

struct Atom {
    std::vector<float> pos; // {x1, y1, z1, x2, y2, z2, ...}
    std::vector<float> vel; // {u1, v1, w1, u2, v2, w2, ...}
    int cnt; // number of atoms

    Atom(std::vector<float> pos, std::vector<float> vel);
};

class Grid {
    public:
    enum CellType {
        SOLID,
        LIQUID,
        AIR,
    };

    int xSize;
    int ySize;
    int zSize;
    std::vector<CellType> cell;
    std::vector<float> vel[3]; // velocity value for each grid cell
    std::vector<float> prevVel[3]; // previous velocity value for each grid cell
    std::vector<bool> s; // true if fluid can be filled (=false if cell type is solid)
    std::vector<float> r[3]; // total sum of weights
    std::vector<CellType> screen;

    Grid(std::vector<std::vector<std::vector<Grid::CellType>>> cell);

    // print grid cross-section in the cut direction according to the given direction
    // direction: 0 -> x, 1 -> y, 2 -> z, 3 -> all
    void printGridSlice(int direction);

    // print 3d atoms projected to 2d screen with given screen size, offset and scale
    void printGridView(int screenWidth, int screenHeight, int offsetX, int offsetY, float scale);
};

class FlipFluid {
    public:
    int xSize;
    int ySize;
    int zSize;
    Atom atom = Atom({},{}); // particles in fluid
    Grid grid = Grid({{{Grid::CellType::AIR}}}); // grid for flip method
    std::vector<float> density; // grid density
    std::vector<int> incompressCounter; // used for multi threading while solving incompressibility

    float gravity;
    int simulateStep;
    float dt;
    float minDist;
    int pushApartIteration;
    int incompressIteration;
    float overRelaxation;
    float stiff;
    float restDensity;
    int numThreads;

    FlipFluid(
        Grid grid, 
        Atom atom,
        float gravity,
        int simulateStep,
        float dt,
        float minDist,
        int pushApartIteration,
        int incompressIteration,
        float overRelaxation,
        float stiff,
        float restDensity,
        int numThreads);

    private:
    // simulate atoms in given range
    void simulatePartialAtoms(float gravity, int step, float dt, int firstAtomIndex, int lastAtomIndex);

    // simulate entire atoms using multi-threading    
    void simulateAtoms(float gravity, int step, float dt);

    // efficient collision detection using spacial hash
    void pushAtomsApart(float minDist, int numIter);

    // update each cell's type in grid
    void updateCellType();

    // transfer velocity of given range of atoms to grid
    void transferPartialAtomVelocityToGrid(int dim, int firstAtomIndex, int lastAtomIndex);

    // transfer entire velocity of atoms to grid using multi-threading
    void transferAtomVelocityToGrid();

    // update density of cells in grid
    void updateDensity();

    // make given section of fluid incompressible using Gauss-Seidal method
    void solvePartialGridIncompressibility(int threadIndex, float overRelaxation, float stiff, float restDensity);

    // make fluid incompressible using multi-threading 
    // thread pipelined per layer
    // TODO: handle numIter properly using thread pooling
    void solveGridIncompressibility(int numIter, float overRelaxation, float stiff, float restDensity);

    // single thread version of solveGridIncompressibility
    void solveGridIncompressibilitySingleThreaded(int numIter, float overRelaxation, float stiff, float restDensity);

    // transfer velocity of grid to given range of atoms 
    void transferGridVelocityToPartialAtoms(int dim, int firstAtomIndex, int lastAtomIndex);

    // transfer velocity of entire grid to atom
    void transferGridVelocityToAtoms();

    // update fluid grid and atom state using multiple threads
    void updateMultiThreaded();

    // update fluid grid and atom state using single thread
    void updateSingleThreaded();

    public:
    // update fluid grid and atom state
    // use updateMultiThreaded or updateSingleThreaded
    void update();

    // print fluid to screen using printGridView
    void printFluid(int screenWidth, int screenHeight, int offsetX, int offsetY, float scale);
};