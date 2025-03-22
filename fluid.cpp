#include <iostream>
#include <vector>
#include <cassert>

#include <chrono>
#include <thread>

struct Atom {
    std::vector<float> pos; // {x1, y1, z1, x2, y2, z2, ...}
    std::vector<float> vel; // {u1, v1, w1, u2, v2, w2, ...}
    int cnt;

    Atom(std::vector<float> pos, std::vector<float> vel) {
        // each atom should have three position and velocity
        assert(pos.size() == vel.size());
        assert(pos.size() % 3 == 0);

        this->pos = pos;
        this->vel = vel;
        this->cnt = pos.size() / 3;
    }
};

class Grid {
    public:
    enum CellType {
        SOLID,
        LIQUID,
        AIR,
    };

    int x_size;
    int y_size;
    int z_size;
    std::vector<CellType> cell;
    std::vector<float> vel[3]; // velocity value for each grid cell
    std::vector<float> prev_vel[3]; // previous velocity value for each grid cell
    std::vector<bool> s; // true if fluid can be filled (=false if cell type is solid)
    std::vector<float> r[3]; // total sum of weights
    std::vector<CellType> screen;

    Grid(std::vector<std::vector<std::vector<Grid::CellType>>> cell) {
        // grid should be 3D nonempty vector
        assert(cell.size()>0 && cell[0].size()>0 && cell[0][0].size()>0);

        // init grid size
        this->x_size = cell.size();
        this->y_size = cell[0].size();
        this->z_size = cell[0][0].size();
        this->cell.resize(x_size * y_size * z_size);

        // init each cell in grid
        // flatten input cell data to 1d vector for cache locality
        for(int i=0;i<x_size;i++) {
            for(int j=0;j<y_size;j++) {
                for(int k=0;k<z_size;k++) {
                    this->cell[(i*y_size + j)*z_size + k] = cell[i][j][k];
                }
            }
        }
        
        // each vel represent following direction
        // prev_vel points same direction as vel
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
        this->vel[0].resize((x_size+1) * y_size * z_size);
        this->vel[1].resize(x_size * (y_size+1) * z_size);
        this->vel[2].resize(x_size * y_size * (z_size+1));
        this->prev_vel[0].resize((x_size+1) * y_size * z_size);
        this->prev_vel[1].resize(x_size * (y_size+1) * z_size);
        this->prev_vel[2].resize(x_size * y_size * (z_size+1));

        // init others
        this->r[0].resize((x_size+1) * (y_size+1) * (z_size+1));
        this->r[1].resize((x_size+1) * (y_size+1) * (z_size+1));
        this->r[2].resize((x_size+1) * (y_size+1) * (z_size+1));
        this->s.resize(x_size * y_size * z_size);
        for(int i=0;i<x_size;i++) {
            for(int j=0;j<y_size;j++) {
                for(int k=0;k<z_size;k++) {
                    if (cell[i][j][k] == Grid::CellType::SOLID)
                        s[(i*y_size + j)*z_size + k] = false;
                    else
                        s[(i*y_size + j)*z_size + k] = true;
                }
            }
        }
    }

    // print grid cross-section in the cut direction according to the given direction
    // direction: 0 -> x, 1 -> y, 2 -> z, 3 -> all
    void print_grid(int direction) {
        if (direction == 0 || direction == 4) {
            int i=x_size/2;
            for(int j=0;j<y_size;j++) {
                for(int k=0;k<z_size;k++) {
                    switch (cell[(i*y_size + j)*z_size + k]) {
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
            // int j=y_size/2;
            int j=y_size-4;
            for(int i=0;i<x_size;i++) {
                for(int k=0;k<z_size;k++) {
                    switch (cell[(i*y_size + j)*z_size + k]) {
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
            int k=z_size/2;
            for(int i=0;i<x_size;i++) {
                for(int j=0;j<y_size;j++) {
                    switch (cell[(i*y_size + j)*z_size + k]) {
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
    void print_screen(int screen_width, int screen_height, int x_offset, int y_offset, float scale) {
        // resize screen
        if (screen_width * screen_height != screen.size()) {
            screen.resize(screen_width * screen_height);
        }

        // init screen
        for(int i=0;i<screen.size();i++) {
            screen[i] = Grid::CellType::AIR;
        }

        for(int i=0;i<x_size;i++) {
            for(int j=0;j<y_size;j++) {
                for(int k=0;k<z_size;k++) {
                    float y = -0.5*(i+k) + j;
                    float x = -0.86602540378*(k-i); // sqrt(3)/2 * (k-i)
                    int yi = y_offset + floor(y*scale);
                    int xi = x_offset + floor(x*scale);
                    if (((i==0 || i==x_size-1) && (j==0 || j==y_size-1))
                     || ((j==0 || j==y_size-1) && (k==0 || k==z_size-1))
                     || ((k==0 || k==z_size-1) && (i==0 || i==x_size-1))) {
                        screen[yi*screen_width + xi] = Grid::CellType::SOLID;
                    }

                    if (cell[(i*y_size + j)*z_size + k] == Grid::CellType::LIQUID) {
                        if (screen[yi*screen_width + xi] == Grid::CellType::AIR)
                            screen[yi*screen_width + xi] = Grid::CellType::LIQUID;
                    }
                }
            }
        }

        for(int i=0;i<screen_height;i++) {
            for(int j=0;j<screen_width;j++) {
                switch (screen[i*screen_width+j]) {
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
};

class FlipFluid {
    private:
    int x_size;
    int y_size;
    int z_size;
    Atom atom = Atom({},{}); // particles in fluid
    Grid grid = Grid({{{Grid::CellType::AIR}}}); // grid for flip method
    std::vector<float> density; // grid density
    std::vector<int> incompress_counter; // used for multi threading while solving incompressibility

    float gravity;
    int simulate_step;
    float dt;
    float min_dist;
    int pushapart_iteration;
    int incompress_iteration;
    float over_relaxation;
    float stiff;
    float rest_density;
    int num_threads;

    public:
    FlipFluid(
        Grid grid, 
        Atom atom,
        float gravity,
        int simulate_step,
        float dt,
        float min_dist,
        int pushapart_iteration,
        int incompress_iteration,
        float over_relaxation,
        float stiff,
        float rest_density,
        int num_threads) {

        // init grid
        this->x_size = grid.x_size;
        this->y_size = grid.y_size;
        this->z_size = grid.z_size;
        this->grid = grid;
        
        // init atoms
        this->atom = atom;

        // init others
        this->density.resize(x_size * y_size * z_size);
        this->incompress_counter.resize(x_size);

        this->gravity = gravity;
        this->simulate_step = simulate_step;
        this->dt = dt;
        this->min_dist = min_dist;
        this->pushapart_iteration = pushapart_iteration;
        this->incompress_iteration = incompress_iteration;
        this->over_relaxation = over_relaxation;
        this->stiff = stiff;
        this->rest_density = rest_density;
        this->num_threads = num_threads;
    }

    // simulate atoms in given range
    void simulate_atom_partial(float gravity, int step, float dt, int atom_start, int atom_end) {
        for(int a=atom_start;a<atom_end;a++) {
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
            bool x_hit = false;
            bool y_hit = false;
            bool z_hit = false;
            for(int i=0;i<step;i++) {
                // atoms that go outside of grid are not simulated
                if (x<0 || x>=x_size 
                 || y<0 || y>=y_size
                 || z<0 || z>=z_size ) {
                    break;
                }

                // if atom is too close to solid type cells, stop moving
                if (!x_hit && grid.cell[(floor(x+dx+rx)*y_size + floor(y))*z_size + floor(z)] != Grid::CellType::SOLID) {
                    x += dx;
                } else {
                    x_hit = true;
                }
                if (!y_hit && grid.cell[(floor(x)*y_size + floor(y+dy+ry))*z_size + floor(z)] != Grid::CellType::SOLID) {
                    y += dy;
                } else {
                    y_hit = true;
                }
                if (!z_hit && grid.cell[(floor(x)*y_size + floor(y))*z_size + floor(z+dz+rz)] != Grid::CellType::SOLID) {
                    z += dz;
                } else {
                    z_hit = true;
                }
            }

            // update atom position and velocity
            atom.pos[a*3 + 0] = x;
            atom.pos[a*3 + 1] = y;
            atom.pos[a*3 + 2] = z;
            if (x_hit) { atom.vel[a*3 + 0] = 0; }
            if (y_hit) { atom.vel[a*3 + 1] = 0; }
            if (z_hit) { atom.vel[a*3 + 2] = 0; }
        }
    }

    // simulate entire atoms using multi-threading    
    void simulate_atom(float gravity, int step, float dt) {
        std::vector<std::thread> threads;
        for(int t=0;t<num_threads;t++) {
            int a_start = atom.cnt * t / num_threads;
            int a_end = atom.cnt * (t+1) / num_threads;
            threads.emplace_back(&FlipFluid::simulate_atom_partial, this, gravity, step, dt, a_start, a_end);
        }

        for (auto& t : threads) {
            t.join();
        }
    }

    // efficient collision detection using spacial hash
    void push_atoms_apart(float min_dist, int num_iter) {
        return; // TODO: add proper seperation function for 3d grid
    }

    // update each cell's type in grid
    void update_cell() {
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
            if (x_cell < 0 || x_cell >= x_size 
             || y_cell < 0 || y_cell >= y_size
             || z_cell < 0 || z_cell >= z_size) // out of grid boundary
                continue;
            if (grid.cell[(x_cell*y_size + y_cell)*z_size + z_cell] == Grid::CellType::SOLID) // do not modify solid type cell
                continue;

            grid.cell[(x_cell*y_size + y_cell)*z_size + z_cell] = Grid::CellType::LIQUID;
        }
    }

    // transfer velocity of given range of atoms to grid
    void transfer_velocity_to_grid_partial(int dim, int atom_start, int atom_end) {
        // assume cell height is 1, offset is 0.5 for each cells
        float offset[3][3] = {
            {0.0, 0.5, 0.5},
            {0.5, 0.0, 0.5},
            {0.5, 0.5, 0.0}
        };

        int dim_size[3][3] = {
            {x_size+1, y_size, z_size},
            {x_size, y_size+1, z_size},
            {x_size, y_size, z_size+1}
        };

        for(int a=atom_start;a<atom_end;a++) {
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
                        if (x_cell+i<0 || x_cell+i>=x_size 
                            || y_cell+j<0 || y_cell+j>=y_size
                            || z_cell+k<0 || z_cell+k>=z_size)
                            continue;
                        
                        grid.vel[dim][((x_cell+i)*dim_size[dim][1] + y_cell+j)*dim_size[dim][2] + z_cell+k] += weight[i][j][k] * atom.vel[a*3 + dim];
                        grid.r[dim][((x_cell+i)*(y_size+1) + y_cell+j)*(z_size+1) + z_cell+k] += weight[i][j][k];
                    }
                }
            }
        }
    }

    // transfer entire velocity of atoms to grid using multi-threading
    void transfer_velocity_to_grid() {
        std::vector<std::thread> threads;
        for(int d=0;d<3;d++) {
            // reset grid velocity components
            for(int i=0;i<grid.vel[d].size();i++)
                grid.vel[d][i] = 0;

            // reset r
            for(int i=0;i<grid.r[d].size();i++)
                grid.r[d][i] = 0;

            threads.emplace_back(&FlipFluid::transfer_velocity_to_grid_partial, this, d, 0, atom.cnt);

            // !race condition!
            // for(int t=0;t<3;t++) {
            //     int atom_start = atom.cnt * t / num_threads;
            //     int atom_end = atom.cnt * (t+1) / num_threads;
            //     threads.emplace_back(&FlipFluid::transfer_velocity_to_grid_partial, this, d, atom_start, atom_end);
            // }
        }

        for (auto& t : threads) {
            t.join();
        }

        int dim_size[3][3] = {
            {x_size+1, y_size, z_size},
            {x_size, y_size+1, z_size},
            {x_size, y_size, z_size+1}
        };

        // normalize velocity to max 1
        for(int d=0;d<3;d++) {
            for(int i=0;i<dim_size[d][0];i++) {
                for(int j=0;j<dim_size[d][1];j++) {
                    for(int k=0;k<dim_size[d][2];k++) {
                        float weight_sum = grid.r[d][(i*(y_size+1) + j)*(z_size+1) + k];
                        if (weight_sum == 0.0)
                            continue;
                        grid.vel[d][(i*dim_size[d][1] + j)*dim_size[d][2] + k] /= weight_sum;
                    }
                }
            }
        }
    }

    // update density of cells in grid
    void update_density() {
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
                        if (x_cell+i<0 || x_cell+i>=x_size 
                         || y_cell+j<0 || y_cell+j>=y_size
                         || z_cell+k<0 || z_cell+k>=z_size)
                            continue;

                        density[((x_cell+i)*y_size + y_cell+j)*z_size + z_cell+k] += weight[i][j][k];
                    }
                }
            }
        }
    }

    // make given section of fluid incompressible using Gauss-Seidal method
    void solve_incompressibility_partial(int thread_idx, float over_relaxation, float stiff, float rest_density) {
        for(int i=0;i<=x_size;i++) {
            // count last layer
            if (i == x_size) {
                incompress_counter[i-1] = thread_idx + 1;
                break;
            }

            // check if the calculation was completed by the previous thread and wait
            while (incompress_counter[i] != thread_idx) {
                std::this_thread::yield(); // give cpu time to other threads
            }
    
            for(int j=0;j<y_size;j++) {
                for(int k=0;k<z_size;k++) {
                    if (grid.cell[(i*y_size + j)*z_size + k] != Grid::CellType::LIQUID) // only consider liquid type cell
                        continue;

                    float divergence = grid.vel[0][((i+1)*y_size + j)*z_size + k] - grid.vel[0][(i*y_size + j)*z_size + k]
                                        + grid.vel[1][(i*(y_size+1) + j+1)*z_size + k] - grid.vel[1][(i*(y_size+1) + j)*z_size + k]
                                        + grid.vel[2][(i*y_size + j)*(z_size+1) + k+1] - grid.vel[2][(i*y_size + j)*(z_size+1) + k];
                    divergence *= over_relaxation;
                    if (density[(i*y_size + j)*z_size + k] > 0)
                        divergence -= stiff*(density[(i*y_size + j)*z_size + k] - rest_density);

                    bool s1 = (i-1>=0)     ? grid.s[((i-1)*y_size + j)*z_size + k] : false;
                    bool s2 = (i+1<x_size) ? grid.s[((i+1)*y_size + j)*z_size + k] : false;
                    bool s3 = (j-1>=0)     ? grid.s[(i*y_size + (j-1))*z_size + k] : false;
                    bool s4 = (j+1<y_size) ? grid.s[(i*y_size + (j+1))*z_size + k] : false;
                    bool s5 = (k-1>=0)     ? grid.s[(i*y_size + j)*z_size + (k-1)] : false;
                    bool s6 = (k+1<z_size) ? grid.s[(i*y_size + j)*z_size + (k+1)] : false;
                    int s = s1+s2+s3+s4+s5+s6;

                    if (s == 0)
                        continue;
                    grid.vel[0][((i+0)*y_size + j)*z_size + k] += divergence * s1 / s;
                    grid.vel[0][((i+1)*y_size + j)*z_size + k] -= divergence * s2 / s;
                    grid.vel[1][(i*(y_size+1) + (j+0))*z_size + k] += divergence * s3 / s;
                    grid.vel[1][(i*(y_size+1) + (j+1))*z_size + k] -= divergence * s4 / s;
                    grid.vel[2][(i*y_size + j)*(z_size+1) + k+0] += divergence * s5 / s;
                    grid.vel[2][(i*y_size + j)*(z_size+1) + k+1] -= divergence * s6 / s;
                }
            }

            // check done layer so other thread can proceed
            if (i-1>=0) {
                incompress_counter[i-1] = thread_idx + 1;
            }
        }

    }

    // make fluid incompressible using multi-threading 
    // thread pipelined per layer
    void solve_incompressibility(int num_iter, float over_relaxation, float stiff, float rest_density) {
        // reset incompress_counter
        for(int i=0;i<incompress_counter.size();i++)
            incompress_counter[i] = 0;

        std::vector<std::thread> threads;
        for(int t=0;t<num_threads;t++) {
            threads.emplace_back(&FlipFluid::solve_incompressibility_partial, this, t, over_relaxation, stiff, rest_density);
        }

        for (auto& t : threads) {
            t.join();
        }
    }

    // transfer velocity of grid to given range of atoms 
    void transfer_velocity_to_atom_partial(int dim, int atom_start, int atom_end) {
        // assume cell height is 1, offset is 0.5 for each cells
        float offset[3][3] = {
            {0.0, 0.5, 0.5},
            {0.5, 0.0, 0.5},
            {0.5, 0.5, 0.0}
        };

        int dim_size[3][3] = {
            {x_size+1, y_size, z_size},
            {x_size, y_size+1, z_size},
            {x_size, y_size, z_size+1}
        };

        // d points vel
        for(int a=0;a<atom.cnt;a++) {
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

            float numerator=0, denominator=0;
            for(int i=0;i<2;i++) {
                for(int j=0;j<2;j++) {
                    for(int k=0;k<2;k++) {
                        if (x_cell+i<0 || x_cell+i>=x_size 
                            || y_cell+j<0 || y_cell+j>=y_size
                            || z_cell+k<0 || z_cell+k>=z_size)
                            continue;
                        int idx = ((x_cell+i)*dim_size[dim][1] + y_cell+j)*dim_size[dim][2] + z_cell+k;
                        // add difference of component to atom
                        numerator += (grid.vel[dim][idx]-grid.prev_vel[dim][idx])*weight[i][j][k];
                        denominator += weight[i][j][k];
                    }
                }
            }
            float quantity = numerator / denominator;
            atom.vel[a*3 + dim] += quantity;
        }
    }

    // transfer velocity of entire grid to atom
    void transfer_velocity_to_atom() {
        std::vector<std::thread> threads;
        for(int d=0;d<3;d++) {
            threads.emplace_back(&FlipFluid::transfer_velocity_to_atom_partial, this, d, 0, atom.cnt);
        }

        for (auto& t : threads) {
            t.join();
        }
    }

    void update() {
        // update atom velocity by gravity
        // update atom position by it's velocity
        // simulate particles
        simulate_atom(gravity, simulate_step, dt);

        push_atoms_apart(min_dist, pushapart_iteration);

        // update cell state
        update_cell();

        // transfer atom velocity to grid
        transfer_velocity_to_grid();

        update_density();

        // save current grid velocity
        for(int d=0;d<3;d++)
            grid.vel[d].swap(grid.prev_vel[d]);
        
        // make fluid incompressible
        // using Gauss-Seidal method
        solve_incompressibility(incompress_iteration, over_relaxation, stiff, rest_density);

        // add change to atoms
        transfer_velocity_to_atom();
    }

    void print_fluid(int screen_width, int screen_height, int offset_x, int offset_y, float scale) {
        grid.print_screen(screen_width, screen_height, offset_x, offset_y, scale);
    }
};

int main() {
    // faster I/O
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    // int x_size = 61, y_size = 41, z_size = 101; // odd dim for easy solve_incompressibility
    int x_size = 60, y_size = 40, z_size = 100;

    // make basic box type grid
    std::vector<std::vector<std::vector<Grid::CellType>>> init_cell;
    init_cell.resize(
        x_size, std::vector<std::vector<Grid::CellType>>(
        y_size, std::vector<Grid::CellType>(
        z_size, Grid::CellType::AIR
    )));

    for(int i=0;i<x_size;i++) {
        for(int j=0;j<y_size;j++) {
            for(int k=0;k<z_size;k++) {
                if (i==0 || i==x_size-1
                 || j==0 || j==y_size-1
                 || k==0 || k==z_size-1) {
                    init_cell[i][j][k] = Grid::CellType::SOLID;
                }
            }
        }
    }
    Grid grid = Grid(init_cell);

    // make half filled atoms
    std::vector<float> pos;
    std::vector<float> vel;
    for(int i=2;i<x_size*3/4;i++) {
        for(int j=2;j<y_size*3/4;j++) {
            for(int k=2;k<z_size/2;k++) {
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
    int simulate_step = 5;
    float dt = 0.03;
    float min_dist = 0.5;
    int pushapart_iteration = 0;
    int incompress_iteration = 5;
    float over_relaxation = 1.9;
    float stiff = 2.0;
    float rest_density = 3.0;
    int num_threads = 5;

    // make fluid
    FlipFluid fluid = FlipFluid(
        grid,
        atom,
        gravity,
        simulate_step,
        dt,
        min_dist,
        pushapart_iteration,
        incompress_iteration,
        over_relaxation,
        stiff,
        rest_density,
        num_threads
    );

    while(true) {
        auto start = std::chrono::high_resolution_clock::now(); // start time
    
        fluid.update();
        fluid.print_fluid(100, 60, 50, 40, 0.4);

        auto end = std::chrono::high_resolution_clock::now(); // end time
        std::chrono::duration<double> duration = end - start; // compute duration
        std::cout << "Time elapsed: " << duration.count() << '\n';
    }
}