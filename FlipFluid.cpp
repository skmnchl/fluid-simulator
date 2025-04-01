#include "FlipFluid.h"

Atom::Atom(std::vector<float> pos, std::vector<float> vel) {
    // each atom should have three position and velocity
    assert(pos.size() == vel.size());
    assert(pos.size() % 3 == 0);

    this->pos = pos;
    this->vel = vel;
    this->cnt = pos.size() / 3;
}

Grid::Grid(std::vector<std::vector<std::vector<Grid::CellType>>> cell) {
    // grid should be 3D nonempty vector
    assert(cell.size()>0 && cell[0].size()>0 && cell[0][0].size()>0);

    // init grid size
    this->xSize = cell.size();
    this->ySize = cell[0].size();
    this->zSize = cell[0][0].size();
    this->cell.resize(xSize * ySize * zSize);

    // init each cell in grid
    // flatten input cell data to 1d vector for cache locality
    for(int i=0;i<xSize;i++) {
        for(int j=0;j<ySize;j++) {
            for(int k=0;k<zSize;k++) {
                this->cell[(i*ySize + j)*zSize + k] = cell[i][j][k];
            }
        }
    }
    
    // each vel represent following direction
    // prevVel points same direction as vel
    // (u, v, w => vel[0], vel[1], vel[2])
    //
    //       y    v
    //       |    |
    //       +----|------+
    //      /|    |     /|
    //     / |         / |
    //    +-----------+  |
    //    |  |        | ------ w
    //    |  +--------|--+-- z
    //    | /   /     | /
    //    |/   /      |/
    //    +---/-------+
    //   /   u
    //  x
    //
    this->vel[0].resize((xSize+1) * ySize * zSize);
    this->vel[1].resize(xSize * (ySize+1) * zSize);
    this->vel[2].resize(xSize * ySize * (zSize+1));
    this->prevVel[0].resize((xSize+1) * ySize * zSize);
    this->prevVel[1].resize(xSize * (ySize+1) * zSize);
    this->prevVel[2].resize(xSize * ySize * (zSize+1));

    // init others
    this->r[0].resize((xSize+1) * (ySize+1) * (zSize+1));
    this->r[1].resize((xSize+1) * (ySize+1) * (zSize+1));
    this->r[2].resize((xSize+1) * (ySize+1) * (zSize+1));
    this->s.resize(xSize * ySize * zSize);
    for(int i=0;i<xSize;i++) {
        for(int j=0;j<ySize;j++) {
            for(int k=0;k<zSize;k++) {
                if (cell[i][j][k] == Grid::CellType::SOLID)
                    s[(i*ySize + j)*zSize + k] = false;
                else
                    s[(i*ySize + j)*zSize + k] = true;
            }
        }
    }
}

// print grid cross-section in the cut direction according to the given direction
// direction: 0 -> x, 1 -> y, 2 -> z, 3 -> all
void Grid::printGridSlice(int direction) {
    if (direction == 0 || direction == 4) {
        int i=xSize/2;
        for(int j=0;j<ySize;j++) {
            for(int k=0;k<zSize;k++) {
                switch (cell[(i*ySize + j)*zSize + k]) {
                    case Grid::CellType::SOLID:
                        std::cout << "■ ";
                        break;
                    case Grid::CellType::LIQUID:
                        std::cout << "□ ";
                        break;
                    case Grid::CellType::AIR:
                        std::cout << "  ";
                        break;
                }
            }
            std::cout << '\n';
        }
        std::cout << "\n\n";
    }
    
    if (direction == 1 || direction == 4) {
        // int j=ySize/2;
        int j=ySize-4;
        for(int i=0;i<xSize;i++) {
            for(int k=0;k<zSize;k++) {
                switch (cell[(i*ySize + j)*zSize + k]) {
                    case Grid::CellType::SOLID:
                        std::cout << "■ ";
                        break;
                    case Grid::CellType::LIQUID:
                        std::cout << "□ ";
                        break;
                    case Grid::CellType::AIR:
                        std::cout << "  ";
                        break;
                }
            }
            std::cout << '\n';
        }
        std::cout << "\n\n";
    }
    
    if (direction == 2 || direction == 4) {
        int k=zSize/2;
        for(int i=0;i<xSize;i++) {
            for(int j=0;j<ySize;j++) {
                switch (cell[(i*ySize + j)*zSize + k]) {
                    case Grid::CellType::SOLID:
                        std::cout << "■ ";
                        break;
                    case Grid::CellType::LIQUID:
                        std::cout << "□ ";
                        break;
                    case Grid::CellType::AIR:
                        std::cout << "  ";
                        break;
                }
            }
            std::cout << '\n';
        }
        std::cout << "\n\n";
    }
}

// print 3d atoms projected to 2d screen with given screen size, offset and scale
void Grid::printGridView(int screenWidth, int screenHeight, int offsetX, int offsetY, float scale) {
    // resize screen
    if (screenWidth * screenHeight != screen.size()) {
        screen.resize(screenWidth * screenHeight);
    }

    // init screen
    for(int i=0;i<screen.size();i++) {
        screen[i] = Grid::CellType::AIR;
    }

    for(int i=0;i<xSize;i++) {
        for(int j=0;j<ySize;j++) {
            for(int k=0;k<zSize;k++) {
                float y = -0.5*(i+k) + j;
                float x = -0.86602540378*(k-i); // sqrt(3)/2 * (k-i)
                int yi = offsetY + floor(y*scale);
                int xi = offsetX + floor(x*scale);
                if (((i==0 || i==xSize-1) && (j==0 || j==ySize-1))
                    || ((j==0 || j==ySize-1) && (k==0 || k==zSize-1))
                    || ((k==0 || k==zSize-1) && (i==0 || i==xSize-1))) {
                    screen[yi*screenWidth + xi] = Grid::CellType::SOLID;
                }

                if (cell[(i*ySize + j)*zSize + k] == Grid::CellType::LIQUID) {
                    if (screen[yi*screenWidth + xi] == Grid::CellType::AIR)
                        screen[yi*screenWidth + xi] = Grid::CellType::LIQUID;
                }
            }
        }
    }

    for(int i=0;i<screenHeight;i++) {
        for(int j=0;j<screenWidth;j++) {
            switch (screen[i*screenWidth+j]) {
                case Grid::CellType::SOLID:
                    std::cout << "■ ";
                    break;
                case Grid::CellType::LIQUID:
                    std::cout << "□ ";
                    break;
                case Grid::CellType::AIR:
                    std::cout << "  ";
                    break;
            }
        }
        std::cout << '\n';
    }
}



FlipFluid::FlipFluid(
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
    int numThreads) {

    // init grid
    this->xSize = grid.xSize;
    this->ySize = grid.ySize;
    this->zSize = grid.zSize;
    this->grid = grid;
    
    // init atoms
    this->atom = atom;

    // init others
    this->density.resize(xSize * ySize * zSize);
    this->incompressCounter.resize(xSize);

    this->gravity = gravity;
    this->simulateStep = simulateStep;
    this->dt = dt;
    this->minDist = minDist;
    this->pushApartIteration = pushApartIteration;
    this->incompressIteration = incompressIteration;
    this->overRelaxation = overRelaxation;
    this->stiff = stiff;
    this->restDensity = restDensity;
    this->numThreads = numThreads;
}

// simulate atoms in given range
void FlipFluid::simulatePartialAtoms(float gravity, int step, float dt, int firstAtomIndex, int lastAtomIndex) {
    for(int a=firstAtomIndex;a<lastAtomIndex;a++) {
        atom.vel[a*3 + 1] += gravity * dt; // add gravity to y direction

        float dx = atom.vel[a*3 + 0] * dt / step;
        float dy = atom.vel[a*3 + 1] * dt / step;
        float dz = atom.vel[a*3 + 2] * dt / step;
        float x = atom.pos[a*3 + 0];
        float y = atom.pos[a*3 + 1];
        float z = atom.pos[a*3 + 2];
        // avoid atoms to be too close to solid type cells        
        float rx = (dx >= 0) ? 0.5 : -0.5;
        float ry = (dy >= 0) ? 0.5 : -0.5; 
        float rz = (dz >= 0) ? 0.5 : -0.5;
        bool xHit = false;
        bool yHit = false;
        bool zHit = false;
        for(int i=0;i<step;i++) {
            // atoms that go outside of grid are not simulated
            if (x<0 || x>=xSize 
                || y<0 || y>=ySize
                || z<0 || z>=zSize ) {
                break;
            }

            // if atom is too close to solid type cells, stop moving
            if (!xHit && grid.cell[(floor(x+dx+rx)*ySize + floor(y))*zSize + floor(z)] != Grid::CellType::SOLID) {
                x += dx;
            } else {
                xHit = true;
            }
            if (!yHit && grid.cell[(floor(x)*ySize + floor(y+dy+ry))*zSize + floor(z)] != Grid::CellType::SOLID) {
                y += dy;
            } else {
                yHit = true;
            }
            if (!zHit && grid.cell[(floor(x)*ySize + floor(y))*zSize + floor(z+dz+rz)] != Grid::CellType::SOLID) {
                z += dz;
            } else {
                zHit = true;
            }
        }

        // update atom position and velocity
        atom.pos[a*3 + 0] = x;
        atom.pos[a*3 + 1] = y;
        atom.pos[a*3 + 2] = z;
        if (xHit) { atom.vel[a*3 + 0] = 0; }
        if (yHit) { atom.vel[a*3 + 1] = 0; }
        if (zHit) { atom.vel[a*3 + 2] = 0; }
    }
}

// simulate entire atoms using multi-threading    
void FlipFluid::simulateAtoms(float gravity, int step, float dt) {
    std::vector<std::thread> threads;
    for(int t=0;t<numThreads;t++) {
        int firstAtomIndex = atom.cnt * t / numThreads;
        int lastAtomIndex = atom.cnt * (t+1) / numThreads;
        threads.emplace_back(&FlipFluid::simulatePartialAtoms, this, gravity, step, dt, firstAtomIndex, lastAtomIndex);
    }

    for (auto& t : threads) {
        t.join();
    }
}

// efficient collision detection using spacial hash
void FlipFluid::pushAtomsApart(float minDist, int numIter) {
    return; // TODO: add proper seperation function for 3d grid
}

// update each cell's type in grid
void FlipFluid::updateCellType() {
    // reset cell types
    for(int i=0;i<grid.cell.size();i++) {
        if (grid.cell[i] != Grid::CellType::SOLID)
            grid.cell[i] = Grid::CellType::AIR;
    }

    // set cell that has atom as liquid
    for(int a=0;a<atom.cnt;a++) {
        int x_cell = floor(atom.pos[a*3 + 0]); // x coordinate of current atom cell
        int y_cell = floor(atom.pos[a*3 + 1]); // y coordinate of current atom cell
        int z_cell = floor(atom.pos[a*3 + 2]); // z coordinate of current atom cell
        if (x_cell < 0 || x_cell >= xSize 
            || y_cell < 0 || y_cell >= ySize
            || z_cell < 0 || z_cell >= zSize) // out of grid boundary
            continue;
        if (grid.cell[(x_cell*ySize + y_cell)*zSize + z_cell] == Grid::CellType::SOLID) // do not modify solid type cell
            continue;

        grid.cell[(x_cell*ySize + y_cell)*zSize + z_cell] = Grid::CellType::LIQUID;
    }
}

// transfer velocity of given range of atoms to grid
void FlipFluid::transferPartialAtomVelocityToGrid(int dim, int firstAtomIndex, int lastAtomIndex) {
    // assume cell height is 1, offset is 0.5 for each cells
    float offset[3][3] = {
        {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5},
        {0.5, 0.5, 0.0}
    };

    int dimSize[3][3] = {
        {xSize+1, ySize, zSize},
        {xSize, ySize+1, zSize},
        {xSize, ySize, zSize+1}
    };

    for(int a=firstAtomIndex;a<lastAtomIndex;a++) {
        // assume cell height is 1
        float x_p = atom.pos[a*3 + 0] - offset[dim][0];
        float y_p = atom.pos[a*3 + 1] - offset[dim][1];
        float z_p = atom.pos[a*3 + 2] - offset[dim][2];
        int x_cell = floor(x_p); // x direction cell number of atom
        int y_cell = floor(y_p); // y direction cell number of atom
        int z_cell = floor(z_p); // z direction cell number of atom
        float dx = x_p - x_cell; // atom x in particular cell
        float dy = y_p - y_cell; // atom y in particular cell
        float dz = z_p - z_cell; // atom z in particular cell
        
        // calculate weight for cell corner weight
        // each weight represent following corner
        //
        //       y
        //       | w[0][1][0]  w[0][1][1]
        //       +-----------+
        //      /|          /| 
        //  w[1][1][0]     w[1][1][1]
        //    +-----------+  |
        //    | w[0][0][0]|  | w[0][0][1]
        //    |  +--------|--+-- z
        //    | /         | /
        //    |/          |/
        //    +-----------+
        //   / w[1][0][0]  w[1][0][1]
        //  x
        float weight[2][2][2] = { 
            {
                { (1-dx)*(1-dy)*(1-dz), (1-dx)*(1-dy)*dz },
                { (1-dx)*dy*(1-dz), (1-dx)*dy*dz }
            },
            {
                { dx*(1-dy)*(1-dz), dx*(1-dy)*dz },
                { dx*dy*(1-dz), dx*dy*dz }
            }
        };

        for(int i=0;i<2;i++) {
            for(int j=0;j<2;j++) {
                for(int k=0;k<2;k++) {
                    if (x_cell+i<0 || x_cell+i>=xSize 
                        || y_cell+j<0 || y_cell+j>=ySize
                        || z_cell+k<0 || z_cell+k>=zSize)
                        continue;
                    
                    grid.vel[dim][((x_cell+i)*dimSize[dim][1] + y_cell+j)*dimSize[dim][2] + z_cell+k] += weight[i][j][k] * atom.vel[a*3 + dim];
                    grid.r[dim][((x_cell+i)*(ySize+1) + y_cell+j)*(zSize+1) + z_cell+k] += weight[i][j][k];
                }
            }
        }
    }
}

// transfer entire velocity of atoms to grid using multi-threading
void FlipFluid::transferAtomVelocityToGrid() {
    std::vector<std::thread> threads;
    for(int d=0;d<3;d++) {
        // reset grid velocity components
        for(int i=0;i<grid.vel[d].size();i++)
            grid.vel[d][i] = 0;

        // reset r
        for(int i=0;i<grid.r[d].size();i++)
            grid.r[d][i] = 0;

        threads.emplace_back(&FlipFluid::transferPartialAtomVelocityToGrid, this, d, 0, atom.cnt);

        // !race condition!
        // for(int t=0;t<3;t++) {
        //     int firstAtomIndex = atom.cnt * t / numThreads;
        //     int lastAtomIndex = atom.cnt * (t+1) / numThreads;
        //     threads.emplace_back(&FlipFluid::transferPartialAtomVelocityToGrid, this, d, firstAtomIndex, lastAtomIndex);
        // }
    }

    for (auto& t : threads) {
        t.join();
    }

    int dimSize[3][3] = {
        {xSize+1, ySize, zSize},
        {xSize, ySize+1, zSize},
        {xSize, ySize, zSize+1}
    };

    // normalize velocity to max 1
    for(int d=0;d<3;d++) {
        for(int i=0;i<dimSize[d][0];i++) {
            for(int j=0;j<dimSize[d][1];j++) {
                for(int k=0;k<dimSize[d][2];k++) {
                    float totalWeight = grid.r[d][(i*(ySize+1) + j)*(zSize+1) + k];
                    if (totalWeight == 0.0)
                        continue;
                    grid.vel[d][(i*dimSize[d][1] + j)*dimSize[d][2] + k] /= totalWeight;
                }
            }
        }
    }
}

// update density of cells in grid
void FlipFluid::updateDensity() {
    // reset density
    for(int i=0;i<density.size();i++)
        density[i] = 0;

    // assume cell height is 1, offset is 0.5 for each cells
    float offset = 0.5;
    for(int a=0;a<atom.cnt;a++) {
        float x_p = atom.pos[a*3 + 0] - offset;
        float y_p = atom.pos[a*3 + 1] - offset;
        float z_p = atom.pos[a*3 + 2] - offset;
        int x_cell = floor(x_p); // x direction cell number of atom
        int y_cell = floor(y_p); // y direction cell number of atom
        int z_cell = floor(z_p); // z direction cell number of atom
        float dx = x_p - x_cell; // atom x in particular cell
        float dy = y_p - y_cell; // atom y in particular cell
        float dz = z_p - z_cell; // atom z in particular cell
                
        // calculate weight for cell corner weight
        float weight[2][2][2] = { 
            {
                { (1-dx)*(1-dy)*(1-dz), (1-dx)*(1-dy)*dz },
                { (1-dx)*dy*(1-dz), (1-dx)*dy*dz }
            },
            {
                { dx*(1-dy)*(1-dz), dx*(1-dy)*dz },
                { dx*dy*(1-dz), dx*dy*dz }
            }
        };

        for(int i=0;i<2;i++) {
            for(int j=0;j<2;j++) {
                for(int k=0;k<2;k++) {
                    if (x_cell+i<0 || x_cell+i>=xSize 
                        || y_cell+j<0 || y_cell+j>=ySize
                        || z_cell+k<0 || z_cell+k>=zSize)
                        continue;

                    density[((x_cell+i)*ySize + y_cell+j)*zSize + z_cell+k] += weight[i][j][k];
                }
            }
        }
    }
}

// make given section of fluid incompressible using Gauss-Seidal method
void FlipFluid::solvePartialGridIncompressibility(int threadIndex, float overRelaxation, float stiff, float restDensity) {
    for(int i=0;i<=xSize;i++) {
        // count last layer
        if (i == xSize) {
            incompressCounter[i-1] = threadIndex + 1;
            break;
        }

        // check if the calculation was completed by the previous thread and wait
        while (incompressCounter[i] != threadIndex) {
            std::this_thread::yield(); // give cpu time to other threads
        }

        for(int j=0;j<ySize;j++) {
            for(int k=0;k<zSize;k++) {
                if (grid.cell[(i*ySize + j)*zSize + k] != Grid::CellType::LIQUID) // only consider liquid type cell
                    continue;

                float divergence = grid.vel[0][((i+1)*ySize + j)*zSize + k] - grid.vel[0][(i*ySize + j)*zSize + k]
                                    + grid.vel[1][(i*(ySize+1) + j+1)*zSize + k] - grid.vel[1][(i*(ySize+1) + j)*zSize + k]
                                    + grid.vel[2][(i*ySize + j)*(zSize+1) + k+1] - grid.vel[2][(i*ySize + j)*(zSize+1) + k];
                divergence *= overRelaxation;
                if (density[(i*ySize + j)*zSize + k] > 0)
                    divergence -= stiff*(density[(i*ySize + j)*zSize + k] - restDensity);

                bool s1 = (i-1>=0)     ? grid.s[((i-1)*ySize + j)*zSize + k] : false;
                bool s2 = (i+1<xSize) ? grid.s[((i+1)*ySize + j)*zSize + k] : false;
                bool s3 = (j-1>=0)     ? grid.s[(i*ySize + (j-1))*zSize + k] : false;
                bool s4 = (j+1<ySize) ? grid.s[(i*ySize + (j+1))*zSize + k] : false;
                bool s5 = (k-1>=0)     ? grid.s[(i*ySize + j)*zSize + (k-1)] : false;
                bool s6 = (k+1<zSize) ? grid.s[(i*ySize + j)*zSize + (k+1)] : false;
                int s = s1+s2+s3+s4+s5+s6;

                if (s == 0)
                    continue;
                grid.vel[0][((i+0)*ySize + j)*zSize + k] += divergence * s1 / s;
                grid.vel[0][((i+1)*ySize + j)*zSize + k] -= divergence * s2 / s;
                grid.vel[1][(i*(ySize+1) + (j+0))*zSize + k] += divergence * s3 / s;
                grid.vel[1][(i*(ySize+1) + (j+1))*zSize + k] -= divergence * s4 / s;
                grid.vel[2][(i*ySize + j)*(zSize+1) + k+0] += divergence * s5 / s;
                grid.vel[2][(i*ySize + j)*(zSize+1) + k+1] -= divergence * s6 / s;
            }
        }

        // check done layer so other thread can proceed
        if (i-1>=0) {
            incompressCounter[i-1] = threadIndex + 1;
        }
    }
}

// make fluid incompressible using multi-threading 
// thread pipelined per layer
// TODO: handle numIter properly using thread pooling
void FlipFluid::solveGridIncompressibility(int numIter, float overRelaxation, float stiff, float restDensity) {
    // reset incompressCounter
    for(int i=0;i<incompressCounter.size();i++)
        incompressCounter[i] = 0;

    std::vector<std::thread> threads;
    for(int t=0;t<numThreads;t++) {
        threads.emplace_back(&FlipFluid::solvePartialGridIncompressibility, this, t, overRelaxation, stiff, restDensity);
    }

    for (auto& t : threads) {
        t.join();
    }
}

// single thread version of solveGridIncompressibility
void FlipFluid::solveGridIncompressibilitySingleThreaded(int numIter, float overRelaxation, float stiff, float restDensity) {
    for(int iter=0;iter<numIter;iter++) {
        for(int i=0;i<xSize;i++) {
            for(int j=0;j<ySize;j++) {
                for(int k=0;k<zSize;k++) {
                    if (grid.cell[(i*ySize + j)*zSize + k] != Grid::CellType::LIQUID) // only consider liquid type cell
                        continue;

                    float divergence = grid.vel[0][((i+1)*ySize + j)*zSize + k] - grid.vel[0][(i*ySize + j)*zSize + k]
                                        + grid.vel[1][(i*(ySize+1) + j+1)*zSize + k] - grid.vel[1][(i*(ySize+1) + j)*zSize + k]
                                        + grid.vel[2][(i*ySize + j)*(zSize+1) + k+1] - grid.vel[2][(i*ySize + j)*(zSize+1) + k];
                    divergence *= overRelaxation;
                    if (density[(i*ySize + j)*zSize + k] > 0)
                        divergence -= stiff*(density[(i*ySize + j)*zSize + k] - restDensity);

                    bool s1 = (i-1>=0)    ? grid.s[((i-1)*ySize + j)*zSize + k] : false;
                    bool s2 = (i+1<xSize) ? grid.s[((i+1)*ySize + j)*zSize + k] : false;
                    bool s3 = (j-1>=0)    ? grid.s[(i*ySize + (j-1))*zSize + k] : false;
                    bool s4 = (j+1<ySize) ? grid.s[(i*ySize + (j+1))*zSize + k] : false;
                    bool s5 = (k-1>=0)    ? grid.s[(i*ySize + j)*zSize + (k-1)] : false;
                    bool s6 = (k+1<zSize) ? grid.s[(i*ySize + j)*zSize + (k+1)] : false;
                    int s = s1+s2+s3+s4+s5+s6;

                    if (s == 0)
                        continue;
                    grid.vel[0][((i+0)*ySize + j)*zSize + k] += divergence * s1 / s;
                    grid.vel[0][((i+1)*ySize + j)*zSize + k] -= divergence * s2 / s;
                    grid.vel[1][(i*(ySize+1) + (j+0))*zSize + k] += divergence * s3 / s;
                    grid.vel[1][(i*(ySize+1) + (j+1))*zSize + k] -= divergence * s4 / s;
                    grid.vel[2][(i*ySize + j)*(zSize+1) + k+0] += divergence * s5 / s;
                    grid.vel[2][(i*ySize + j)*(zSize+1) + k+1] -= divergence * s6 / s;
                }
            }
        }
    }
}

// transfer velocity of grid to given range of atoms 
void FlipFluid::transferGridVelocityToPartialAtoms(int dim, int firstAtomIndex, int lastAtomIndex) {
    // assume cell height is 1, offset is 0.5 for each cells
    float offset[3][3] = {
        {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5},
        {0.5, 0.5, 0.0}
    };

    int dimSize[3][3] = {
        {xSize+1, ySize, zSize},
        {xSize, ySize+1, zSize},
        {xSize, ySize, zSize+1}
    };

    // d points vel
    for(int a=firstAtomIndex;a<lastAtomIndex;a++) {
        // assume cell height is 1
        float x_p = atom.pos[a*3 + 0] - offset[dim][0];
        float y_p = atom.pos[a*3 + 1] - offset[dim][1];
        float z_p = atom.pos[a*3 + 2] - offset[dim][2];
        int x_cell = floor(x_p); // x direction cell number of atom
        int y_cell = floor(y_p); // y direction cell number of atom
        int z_cell = floor(z_p); // z direction cell number of atom
        float dx = x_p - x_cell; // atom x in particular cell
        float dy = y_p - y_cell; // atom y in particular cell
        float dz = z_p - z_cell; // atom z in particular cell
        
        // calculate weight for cell corner weight
        float weight[2][2][2] = { 
            {
                { (1-dx)*(1-dy)*(1-dz), (1-dx)*(1-dy)*dz },
                { (1-dx)*dy*(1-dz), (1-dx)*dy*dz }
            },
            {
                { dx*(1-dy)*(1-dz), dx*(1-dy)*dz },
                { dx*dy*(1-dz), dx*dy*dz }
            }
        };

        float totalVelocity=0, totalWeight=0;
        for(int i=0;i<2;i++) {
            for(int j=0;j<2;j++) {
                for(int k=0;k<2;k++) {
                    if (x_cell+i<0 || x_cell+i>=xSize 
                        || y_cell+j<0 || y_cell+j>=ySize
                        || z_cell+k<0 || z_cell+k>=zSize)
                        continue;
                    int idx = ((x_cell+i)*dimSize[dim][1] + y_cell+j)*dimSize[dim][2] + z_cell+k;
                    // add difference of component to atom
                    totalVelocity += (grid.vel[dim][idx]-grid.prevVel[dim][idx])*weight[i][j][k];
                    totalWeight += weight[i][j][k];
                }
            }
        }
        atom.vel[a*3 + dim] += totalVelocity / totalWeight;
    }
}

// transfer velocity of entire grid to atom
void FlipFluid::transferGridVelocityToAtoms() {
    std::vector<std::thread> threads;
    for(int d=0;d<3;d++) {
        threads.emplace_back(&FlipFluid::transferGridVelocityToPartialAtoms, this, d, 0, atom.cnt);
    }

    for (auto& t : threads) {
        t.join();
    }
}

// update fluid grid and atom state using multiple threads
void FlipFluid::updateMultiThreaded() {
    // update atom velocity by gravity
    // update atom position by it's velocity
    // simulate particles
    simulateAtoms(gravity, simulateStep, dt);

    pushAtomsApart(minDist, pushApartIteration);

    // update cell state
    updateCellType();

    // transfer atom velocity to grid
    transferAtomVelocityToGrid();

    updateDensity();

    // save current grid velocity
    for(int d=0;d<3;d++)
        grid.vel[d].swap(grid.prevVel[d]);
    
    // make fluid incompressible
    // using Gauss-Seidal method
    solveGridIncompressibility(incompressIteration, overRelaxation, stiff, restDensity);

    // add change to atoms
    transferGridVelocityToAtoms();
}

// update fluid grid and atom state using single thread
void FlipFluid::updateSingleThreaded() {
    // update atom velocity by gravity
    // update atom position by it's velocity
    // simulate particles
    simulatePartialAtoms(gravity, simulateStep, dt, 0, atom.cnt);

    pushAtomsApart(minDist, pushApartIteration);

    // update cell state
    updateCellType();

    int dimSize[3][3] = {
        {xSize+1, ySize, zSize},
        {xSize, ySize+1, zSize},
        {xSize, ySize, zSize+1}
    };

    // transfer atom velocity to grid
    for(int d=0;d<3;d++) {
        // reset grid velocity components
        for(int i=0;i<grid.vel[d].size();i++)
            grid.vel[d][i] = 0;

        // reset r
        for(int i=0;i<grid.r[d].size();i++)
            grid.r[d][i] = 0;

        transferPartialAtomVelocityToGrid(d, 0, atom.cnt);

        // normalize velocity to max 1
        for(int i=0;i<dimSize[d][0];i++) {
            for(int j=0;j<dimSize[d][1];j++) {
                for(int k=0;k<dimSize[d][2];k++) {
                    float totalWeight = grid.r[d][(i*(ySize+1) + j)*(zSize+1) + k];
                    if (totalWeight == 0.0)
                        continue;
                    grid.vel[d][(i*dimSize[d][1] + j)*dimSize[d][2] + k] /= totalWeight;
                }
            }
        }
    }

    updateDensity();

    // save current grid velocity
    for(int d=0;d<3;d++)
        grid.vel[d].swap(grid.prevVel[d]);
    
    // make fluid incompressible
    // using Gauss-Seidal method
    solveGridIncompressibilitySingleThreaded(incompressIteration, overRelaxation, stiff, restDensity);

    // add change to atoms
    for(int d=0;d<3;d++)
        transferGridVelocityToPartialAtoms(d, 0, atom.cnt);
}

void FlipFluid::update() {
    if (numThreads > 0) {
        updateMultiThreaded();
    } else {
        updateSingleThreaded();
    }
}

void FlipFluid::printFluid(int screenWidth, int screenHeight, int offsetX, int offsetY, float scale) {
    grid.printGridView(screenWidth, screenHeight, offsetX, offsetY, scale);
}