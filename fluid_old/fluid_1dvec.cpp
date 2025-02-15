#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <chrono>
#include <thread>

// return random float between 0 and 1
class RNG {
    public:
    RNG() {
        srand (static_cast <unsigned> (time(0)));
    }

    float get_rand_float() {
        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        return r;
    }

};

class Matrix {
    public:
    std::vector<float> data;  // 1D storage
    size_t rows, cols;

    public:
    Matrix() 
        : rows(0), cols(0), data(0, 0) {}

    Matrix(size_t rows, size_t cols) 
        : rows(rows), cols(cols), data(rows * cols, 0) {}

    // Accessor for (row, col) indexing
    float& operator()(size_t i, size_t j) {
        return data[i * cols + j];  // Convert 2D index to 1D
    }

    const float& operator()(size_t i, size_t j) const {
        return data[i * cols + j];
    }

    void resize(size_t rows, size_t cols) {
        this->rows = rows;
        this->cols = cols;
        this->data.resize(rows * cols);
    }

    void reset() {
        for(float &d : data)
            d = 0;
    }

    // Get raw 1D data pointer (useful for optimizations)
    float* raw_data() { return data.data(); }

    // Get size information
    size_t get_rows() const { return rows; }
    size_t get_cols() const { return cols; }
};

struct Atom {
    float pos[2];
    float vel[2];

    Atom(float x, float y, float u, float v) {
        this->pos[0] = x;
        this->pos[1] = y;
        this->vel[0] = u;
        this->vel[1] = v;
    }
};

class Grid {
    public:
    enum CellType {
        SOLID,
        LIQUID,
        AIR,
    };

    int height;
    int width;
    std::vector<std::vector<CellType>> cell;
    Matrix vel[2]; // grid velocity
    Matrix prev_vel[2]; // grid previous velocity
    Matrix density; // grid density
    // std::vector<std::vector<bool>> s;
    Matrix r;
    float rest_density = 0;

    Grid(std::vector<std::vector<Grid::CellType>> cell) {
        assert(cell.size() > 0 && cell[0].size() > 0);
        
        // init grid size
        this->height = cell.size();
        this->width = cell[0].size();
        this->cell = cell;
        
        // init grid velocity
        this->vel[0].resize(height, width+1);
        this->vel[1].resize(height+1, width);
        this->prev_vel[0].resize(height, width+1);
        this->prev_vel[1].resize(height+1, width);

        // init others
        this->density.resize(height, width);
        this->r.resize(height+1, width+1);

        // this->s.resize(height, std::vector<float>(width));
    }

    void reset_cell() {
        for(int i=0;i<height;i++) {
            for(int j=0;j<width;j++) {
                if (cell[i][j] != Grid::CellType::SOLID)
                    cell[i][j] = Grid::CellType::AIR;
            }
        }
    }

    void reset_velocity() {
        for(int d=0;d<2;d++)
            vel[d].reset();
    }

    void reset_density() {
        density.reset();
    }

    void reset_r() {
        r.reset();
    }

    void print_grid() {
        for(int i=0;i<height;i++) {
            for(int j=0;j<width;j++) {
                switch (cell[i][j]) {
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

    void print_grid_status(std::string msg) {
        return;

        std::cout << msg;
        for(int i=0;i<height*2+1;i++) {
            for(int j=0;j<width+1;j++) {
                if (i%2==0 && j<width) {
                    if (vel[1](i/2, j) >= 0) {
                        std::cout.precision(2);
                        std::cout << std::fixed;
                    } else {
                        std::cout.precision(1);
                        std::cout << std::fixed;
                    }
                    std::cout << "    " << vel[1](i/2, j);
                }
                else if (i%2==1 && i/2<height) {
                    if (vel[0](i/2, j) >= 0) {
                        std::cout.precision(2);
                        std::cout << std::fixed;
                    } else {
                        std::cout.precision(1);
                        std::cout << std::fixed;
                    }
                    std::cout << vel[0](i/2, j);
                    if (i/2>=height || j>=width)
                        continue;
                    if (cell[i/2][j] == Grid::CellType::LIQUID) std::cout << "□ □ ";
                    else if (cell[i/2][j] == Grid::CellType::SOLID) std::cout << "■ ■ ";
                    else std::cout << "    ";
                }
            }
            std::cout << '\n';
        }
    }
};

class FlipFluid {
    private:
    int height;
    int width;
    std::vector<Atom> atom;
    Grid grid = Grid({{Grid::CellType::AIR}});
    std::vector<std::vector<std::vector<int>>> cell_atom;

    public:
    FlipFluid(Grid grid, std::vector<Atom> atom) {
        // init grid
        this->height = grid.height;
        this->width = grid.width;
        this->grid = grid;
        
        // init atoms
        this->atom = atom;

        // init others
        this->cell_atom.resize(height, std::vector<std::vector<int>>(width, std::vector<int>()));
    }

    void reset_cell_atom() {
        for(int i=0;i<cell_atom.size();i++)
            for(int j=0;j<cell_atom[0].size();j++)
                cell_atom[i][j].clear();
    }

    void transfer_velocity_to_atom() {
        // assume cell height is 1
        float offset[2][2] = {
            {0, 0.5},
            {0.5, 0}
        };

        // d points vel
        for(int d=0;d<2;d++) {    
            for(Atom &a : atom) {
                // assume cell height is 1
                float x_p = a.pos[0] - offset[d][0];
                float y_p = a.pos[1] - offset[d][1];
                int x_cell = floor(x_p); // x direction cell number of atom
                int y_cell = floor(y_p); // y direction cell number of atom
                float dx = x_p - x_cell; // atom x in particular cell
                float dy = y_p - y_cell; // atom y in particular cell
                
                // calculate weight for cell corner weight
                float weight[2][2] = {
                    { (1-dx)*(1-dy), dx*(1-dy) },
                    { (1-dx)*dy, dx*dy },
                };

                float numerator=0, denominator=0;
                for(int i=0;i<2;i++) {
                    for(int j=0;j<2;j++) {
                        if (y_cell+i<0 || y_cell+i>=height || x_cell+j<0 || x_cell+j>=width)
                            continue;
                        float q_curr = grid.vel[d](y_cell+i, x_cell+j);
                        float q_prev = grid.prev_vel[d](y_cell+i, x_cell+j);
                        // add difference of component to atom
                        numerator += (q_curr-q_prev)*weight[i][j];
                        denominator += weight[i][j];
                    }
                }
                float quantity = numerator / denominator;
                a.vel[d] += quantity;
            }
        }
    }

    void transfer_velocity_to_grid() {
        // reset grid velocity components
        grid.reset_velocity();
        grid.print_grid_status("\n\nreset grid component\n");
        
        // assume cell height is 1
        float offset[2][2] = {
            {0, 0.5},
            {0.5, 0}
        };

        // d points vel
        for(int d=0;d<2;d++) {
            grid.reset_r();

            for(Atom &a : atom) {
                // assume cell height is 1
                float x_p = a.pos[0] - offset[d][0];
                float y_p = a.pos[1] - offset[d][1];
                int x_cell = floor(x_p); // x direction cell number of atom
                int y_cell = floor(y_p); // y direction cell number of atom
                float dx = x_p - x_cell; // atom x in particular cell
                float dy = y_p - y_cell; // atom y in particular cell
                
                // calculate weight for cell corner weight
                float weight[2][2] = {
                    { (1-dx)*(1-dy), dx*(1-dy) },
                    { (1-dx)*dy, dx*dy },
                };

                for(int i=0;i<2;i++) {
                    for(int j=0;j<2;j++) {
                        if (y_cell+i<0 || y_cell+i>=height || x_cell+j<0 || x_cell+j>=width)
                            continue;
                        
                        grid.vel[d](y_cell+i, x_cell+j) += weight[i][j] * a.vel[d];
                        grid.r(y_cell+i, x_cell+j) += weight[i][j];
                    }
                }
            }

            grid.print_grid_status("\n\nbefore dividing with r\n");

            for(int i=0;i<grid.vel[d].get_rows();i++) {
                for(int j=0;j<grid.vel[d].get_cols();j++) {
                    if (grid.r(i, j) == 0)
                        continue;
                    grid.vel[d](i, j) /= grid.r(i, j);
                }
            }

            grid.print_grid_status("\n\nafter dividing with r\n");
        }
    }

    void force_incompressible(int num_iter, float over_relaxation, float stiff) {
        for(int iter=0;iter<num_iter;iter++) {
            for(int i=0;i<height;i++) {
                for(int j=0;j<width;j++) {
                    if (grid.cell[i][j] != Grid::CellType::LIQUID)
                        continue;
                    float divergence = over_relaxation*(grid.vel[0](i, j+1) - grid.vel[0](i, j) + grid.vel[1](i+1, j) - grid.vel[1](i, j));
                    if (stiff*(grid.density(i, j) - grid.rest_density) > 0)
                        divergence -= stiff*(grid.density(i, j) - grid.rest_density);
                    int s1 = (j-1<1 || grid.cell[i][j-1] == Grid::CellType::SOLID) ? 0 : 1;
                    int s2 = (j+1>=width-1 || grid.cell[i][j+1] == Grid::CellType::SOLID) ? 0 : 1;
                    int s3 = (i-1<1 || grid.cell[i-1][j] == Grid::CellType::SOLID) ? 0 : 1;
                    int s4 = (i+1>=height-1 || grid.cell[i+1][j] == Grid::CellType::SOLID) ? 0 : 1;
                    int s = s1+s2+s3+s4;

                    if (s == 0)
                        continue;
                    grid.vel[0](i, j) += divergence * s1 / s;
                    grid.vel[0](i, j+1) -= divergence * s2 / s;
                    grid.vel[1](i, j) += divergence * s3 / s;
                    grid.vel[1](i+1, j) -= divergence * s4 / s;
                }
            }
        }
    }

    void push_atoms_apart(float min_dist, int num_iter) {
        RNG rng;
        // use grid method
        // divide each atom into grid
        for(int iter=0;iter<num_iter;iter++) {
            reset_cell_atom();
            // consider near 9 cells
            for(int i=0;i<atom.size();i++) {
                Atom a = atom[i];
                int a_x = floor(a.pos[0]);
                int a_y = floor(a.pos[1]);
                int cell_x[3] = {a_x-1, a_x, a_x+1};
                int cell_y[3] = {a_y-1, a_y, a_y+1};
                for(int i=0;i<3;i++) {
                    for(int j=0;j<3;j++) {
                        if (cell_x[j]<0 || cell_x[j]>=width || cell_y[i]<0 || cell_y[i]>=height)
                            continue;
                        cell_atom[cell_y[i]][cell_x[j]].push_back(i);
                    }
                }
            }

            // consider near 9 cells
            for(int i=0;i<height;i++) {
                for(int j=0;j<width;j++) {
                    // see if atoms collide in each cell
                    for(int p=0;p<cell_atom[i][j].size();p++) {
                        for(int q=p+1;q<cell_atom[i][j].size();q++) {
                            Atom &a1 = atom[p];
                            Atom &a2 = atom[q];
                            float dist = sqrt(pow(a1.pos[0]-a2.pos[0],2)+pow(a1.pos[1]-a2.pos[1],2));
                            if (dist < min_dist) {
                                float mid_x = (a1.pos[0]+a2.pos[0]) / 2;
                                float mid_y = (a1.pos[1]+a2.pos[1]) / 2;

                                if (dist == 0) {
                                    a1.pos[0] += min_dist / 2;
                                    a2.pos[0] -= min_dist / 2;
                                    // a1.pos[1] += min_dist / 2;
                                    // a2.pos[1] -= min_dist / 2;
                                    continue;
                                }

                                a1.pos[0] = (a1.pos[0]-mid_x) * min_dist / dist;
                                a2.pos[0] = (a2.pos[0]-mid_x) * min_dist / dist;
                                a1.pos[1] = (a1.pos[1]-mid_y) * min_dist / dist;
                                a2.pos[1] = (a2.pos[1]-mid_y) * min_dist / dist;
                            }
                        }
                    }

                }
            }
        }
    }

    void handle_atom_collisions() {
        for(Atom &a : atom) {
            if (a.pos[0] <= 1 || a.pos[0] >= width-2) {
                a.pos[0] = a.pos[0]>width/2 ? width-2 : 1;
                a.vel[0] = 0;
            }
            if (a.pos[1] <= 1 || a.pos[1] >= height-2) {
                a.pos[1] = a.pos[1]>height/2 ? height-2 : 1;
                a.vel[1] = 0;
            }
        }
    }

    void update_atom_density() {
        grid.reset_density();
        // assume cell height is 1
        float offset = 0.5;
        for(Atom a : atom) {
            float x_p = a.pos[0]-offset;
            float y_p = a.pos[1]-offset;
            int x_cell = floor(x_p); // x direction cell number of atom
            int y_cell = floor(y_p); // y direction cell number of atom
            float dx = x_p - x_cell; // atom x in particular cell
            float dy = y_p - y_cell; // atom y in particular cell
                    
            // calculate weight for cell corner weight
            float weight[2][2] = {
                { (1-dx)*(1-dy), dx*(1-dy) },
                { (1-dx)*dy, dx*dy },
            };

            for(int i=0;i<2;i++) {
                for(int j=0;j<2;j++) {
                    if (y_cell+i<0 || y_cell+i>=height || x_cell+j<0 || x_cell+j>=width)
                        continue;
                    grid.density(y_cell+i, x_cell+j) += weight[i][j];
                }
            }
        }

        if (grid.rest_density == 0) {
            float sum = 0;
            int num_fluid_cells = 0;

            for (int i=0;i<height;i++) {
                for(int j=0;j<width;j++) {
                    if (grid.cell[i][j] == Grid::CellType::LIQUID) {
                        sum += grid.density(i, j);
                        num_fluid_cells+=1;
                    }
                }
            }

            if (num_fluid_cells > 0)
                grid.rest_density = sum / num_fluid_cells;
        }
    }

    void update(float gravity, float dt, int num_iter) {
        // update atom velocity by gravity
        // update atom position by it's velocity
        // simulate particles
        for(Atom &a : atom) {
            a.vel[1] += gravity * dt;
            a.pos[0] += a.vel[0] * dt;
            a.pos[1] += a.vel[1] * dt;
        }

        // push_atoms_apart(1, 2);
        handle_atom_collisions();

        grid.print_grid_status("\n\ninit grid:\n");

        // reset cell types
        grid.reset_cell();

        // set cell that has atom as liquid
        for(Atom a : atom) {
            int x_cell = floor(a.pos[0]);
            int y_cell = floor(a.pos[1]);
            if (x_cell < 0 || x_cell >= width || y_cell < 0 || y_cell >= height)
                continue;
            if (grid.cell[y_cell][x_cell] == Grid::CellType::SOLID)
                continue;

            grid.cell[y_cell][x_cell] = Grid::CellType::LIQUID;
        }

        // transfer atom velocity to grid
        transfer_velocity_to_grid();
        grid.print_grid_status("\n\natom->grid:\n");

        update_atom_density();

        // save current grid velocity
        for(int d=0;d<2;d++)
            grid.prev_vel[d].data.assign(grid.vel[d].data.begin(), grid.vel[d].data.end());
        
        // make fluid incompressible
        // using Gauss-Seidal method
        force_incompressible(num_iter, 1.9, 1);
        grid.print_grid_status("\n\nincompressible:\n");

        // add change to atoms
        transfer_velocity_to_atom();
    }

    void print_fluid() {
        grid.print_grid();
    }
};

int main() {
    // faster I/O
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    // int height=40, width=50, num_atoms=100000;
    // int height = 12, width = 12, num_atoms = 5;
    int height = 80, width = 120, num_atoms = 5;

    // make basic box type grid
    std::vector<std::vector<Grid::CellType>> init_cell;
    init_cell.resize(height, std::vector<Grid::CellType>(width, Grid::CellType::AIR));
    for(int i=0;i<height;i++) {
        init_cell[i][0] = Grid::CellType::SOLID;
        init_cell[i][width-1] = Grid::CellType::SOLID;
    }
    for(int i=0;i<width;i++) {
        init_cell[0][i] = Grid::CellType::SOLID;
        init_cell[height-1][i] = Grid::CellType::SOLID;
    }

    Grid grid = Grid(init_cell);

    // make half filled atoms
    RNG rng;
    std::vector<Atom> atom;
    for(int i=2;i<(height-2)*3/4;i++) {
        for(int j=2;j<width*3/4;j++) {
            float x = j;//+width/4;
            float y = height/4+i-2;
            float u=0, v=0;
            Atom a = Atom(x, y, u, v);
            atom.push_back(a);
        }
    }


    // make fluid
    FlipFluid fluid = FlipFluid(grid, atom);

    while(true) {
        auto start = std::chrono::high_resolution_clock::now(); // Start time
    
        fluid.update(9.81*2, 0.04, 50);
        fluid.print_fluid();
        // int temp;
        // std::cin >> temp;
        // std::this_thread::sleep_for(std::chrono::milliseconds(20));

        auto end = std::chrono::high_resolution_clock::now(); // End time
        std::chrono::duration<double> duration = end - start; // Compute duration
        std::cout << "Time elapsed: " << duration.count() << '\n';
    }
}