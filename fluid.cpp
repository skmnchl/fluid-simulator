#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <chrono>
#include <thread>

struct Atom {
    std::vector<float> pos; // {x1, y1, x2, y2, ...}
    std::vector<float> vel; // {u1, v1, u2, v2, ...}
    int cnt;

    Atom(std::vector<float> pos, std::vector<float> vel) {
        assert(pos.size() == vel.size());
        assert(pos.size() % 2 == 0);

        this->pos = pos;
        this->vel = vel;
        this->cnt = pos.size()/2;
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
    std::vector<CellType> cell;
    std::vector<float> vel[2]; // grid velocity
    std::vector<float> prev_vel[2]; // grid previous velocity
    std::vector<bool> s; // true if fluid can be filled
    std::vector<float> r; // total sum of weights

    Grid(std::vector<std::vector<Grid::CellType>> cell) {
        assert(cell.size() > 0 && cell[0].size() > 0);
        
        // init grid size
        this->height = cell.size();
        this->width = cell[0].size();
        this->cell.resize(width * height);
        for(int i=0;i<height;i++) {
            for(int j=0;j<width;j++) {
                this->cell[i * width + j] = cell[i][j];
            }
        }
        
        // init grid velocity
        this->vel[0].resize(height * (width+1));
        this->vel[1].resize((height+1) * width);
        this->prev_vel[0].resize(height * (width+1));
        this->prev_vel[1].resize((height+1) * width);

        // init others
        this->r.resize((height+1) * (width+1));
        this->s.resize(height * width);
        for(int i=0;i<height;i++) {
            for(int j=0;j<width;j++) {
                if (cell[i][j] == Grid::CellType::SOLID)
                    s[i*width+j] = false;
                else
                    s[i*width+j] = true;
            }
        }

    }

    void reset_cell() {
        for(int i=0;i<height;i++) {
            for(int j=0;j<width;j++) {
                if (cell[i*width+j] != Grid::CellType::SOLID)
                    cell[i*width+j] = Grid::CellType::AIR;
            }
        }
    }

    void reset_velocity() {
        for(int d=0;d<2;d++)
            for(int i=0;i<vel[d].size();i++)
                vel[d][i] = 0;
    }

    void reset_r() {
        for(int i=0;i<r.size();i++)
            r[i] = 0;
    }

    void print_grid() {
        for(int i=0;i<height;i++) {
            for(int j=0;j<width;j++) {
                switch (cell[i*width+j]) {
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
        std::cout << msg;
        for(int i=0;i<height*2+1;i++) {
            for(int j=0;j<width+1;j++) {
                if (i%2==0 && j<width) {
                    if (vel[1][i/2*width+j] >= 0) {
                        std::cout.precision(2);
                        std::cout << std::fixed;
                    } else {
                        std::cout.precision(1);
                        std::cout << std::fixed;
                    }
                    std::cout << "    " << vel[1][i/2*width+j];
                }
                else if (i%2==1 && i/2<height) {
                    if (vel[0][i/2*(width+1)+j] >= 0) {
                        std::cout.precision(2);
                        std::cout << std::fixed;
                    } else {
                        std::cout.precision(1);
                        std::cout << std::fixed;
                    }
                    std::cout << vel[0][i/2*(width+1)+j];
                    if (i/2>=height || j>=width)
                        continue;
                    if (cell[i/2*width+j] == Grid::CellType::LIQUID) std::cout << "□ □ ";
                    else if (cell[i/2*width+j] == Grid::CellType::SOLID) std::cout << "■ ■ ";
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
    Atom atom = Atom({},{});
    Grid grid = Grid({{Grid::CellType::AIR}});
    std::vector<int> num_cell_atom;
    std::vector<int> cell_atom_idx;
    std::vector<float> density; // grid density

    float gravity;
    float dt;
    float min_dist;
    int pushapart_iteration;
    int incompress_iteration;
    float over_relaxation;
    float stiff;
    float rest_density;

    public:
    FlipFluid(
        Grid grid, 
        Atom atom,
        float gravity,
        float dt,
        float min_dist,
        int pushapart_iteration,
        int incompress_iteration,
        float over_relaxation,
        float stiff,
        float rest_density) {

        // init grid
        this->height = grid.height;
        this->width = grid.width;
        this->grid = grid;
        
        // init atoms
        this->atom = atom;

        // init others
        this->num_cell_atom.resize(height * width + 1);
        this->cell_atom_idx.resize(atom.cnt);
        this->density.resize(height * width);
        this->gravity = gravity;
        this->dt = dt;
        this->min_dist = min_dist;
        this->pushapart_iteration = pushapart_iteration;
        this->incompress_iteration = incompress_iteration;
        this->over_relaxation = over_relaxation;
        this->stiff = stiff;
        this->rest_density = rest_density;
    }

    void reset_num_cell_atom() {
        for(int i=0;i<num_cell_atom.size();i++)
            num_cell_atom[i] = 0;
    }
    
    void reset_cell_atom_idx() {
        for(int i=0;i<cell_atom_idx.size();i++)
            cell_atom_idx[i] = 0;
    }

    void reset_density() {
        for(int i=0;i<density.size();i++)
            density[i] = 0;
    }

    void simulate_atom(float gravity, float dt) {
        int step = 10;
        for(int a=0;a<atom.cnt;a++) {
            atom.vel[a*2 + 1] += gravity * dt; // add gravity

            float dx = atom.vel[a*2] * dt / step;
            float dy = atom.vel[a*2 + 1] * dt / step;
            float x0 = atom.pos[a*2];
            float y0 = atom.pos[a*2 + 1];
            float x = x0;
            float y = y0;
            float rx = dx>0 ? min_dist : -min_dist;
            float ry = dy>0 ? min_dist : -min_dist;
            for(int i=0;i<step;i++) {
                if (x>=width || x<=0 || y>=height || y<=0) {
                    break;
                }

                if (grid.cell[floor(y+ry) * width + floor(x+rx)] == Grid::CellType::SOLID) {
                    if (dx == 0) {
                        atom.vel[a*2 + 1] = 0;
                        break;
                    } else if (dy == 0) {
                        atom.vel[a*2] = 0;
                        break;
                    }

                    float corner_x = floor(x) + (dx>0 ? -min_dist : 1+min_dist);
                    float corner_y = floor(y) + (dy>0 ? -min_dist : 1+min_dist);
                    float cross_product = (corner_x - x0) * (y - y0) - (corner_y - y0) * (x - x0);
                    float direction = -cross_product * dx * dy;
                    if (direction > 0) {
                        atom.vel[a*2] = 0;
                    } else if (direction < 0) {
                        atom.vel[a*2 + 1] = 0;
                    } else {
                        atom.vel[a*2] = 0;
                        atom.vel[a*2 + 1] = 0;
                    }

                    break;
                }

                x += dx;
                y += dy;
            }

            atom.pos[a*2] = x;
            atom.pos[a*2 + 1] = y;
        }
    }

    void push_atoms_apart(float min_dist, int num_iter) {
        // use grid method
        // divide each atom into grid
        reset_num_cell_atom();
        reset_cell_atom_idx();
        
        // count number of atoms in each cell 
        for(int a=0;a<atom.cnt;a++) {
            float xi = floor(atom.pos[a*2]);
            float yi = floor(atom.pos[a*2 + 1]);
            if (yi>=height || yi<0 || xi>=height || xi<0)
                continue;
            num_cell_atom[yi * width + xi] += 1;
        }

        // make partial sum of number of atoms in each cell
        int atom_cnt = 0;
        for(int i=0;i<num_cell_atom.size();i++) {
            atom_cnt += num_cell_atom[i];
            num_cell_atom[i] = atom_cnt;
        }

        // fill atom idx in seperate vector
        for(int a=0;a<atom.cnt;a++) {
            float xi = floor(atom.pos[a*2]);
            float yi = floor(atom.pos[a*2 + 1]);
            if (yi>=height || yi<0 || xi>=height || xi<0)
                continue;
            num_cell_atom[yi * width + xi] -= 1;
            cell_atom_idx[num_cell_atom[yi * width + xi]] = a;
        }

        // push atoms apart
        for(int iter=0;iter<num_iter;iter++) {
            for(int a1=0;a1<atom.cnt;a1++) {
                // consider near 9 cells
                float px = atom.pos[2*a1];
                float py = atom.pos[2*a1 + 1];
                int xi = floor(px);
                int yi = floor(py);
                int x0 = (xi-1>=0) ? xi-1 : 0;
                int x1 = (xi+1<width) ? xi+1 : width-1;
                int y0 = (yi-1>=0) ? yi-1 : 0;
                int y1 = (yi+1<height) ? yi+1 : height-1;

                for(int y=y0;y<=y1;y++) {
                    for(int x=x0;x<=x1;x++) {
                        int cell_num = y * width + x;
                        int first = num_cell_atom[cell_num];
                        int last = num_cell_atom[cell_num + 1];
                        for (int atom_id = first; atom_id < last; atom_id++) {
                            int a2 = cell_atom_idx[atom_id];
                            if (a2 == a1)
                                continue;
                            float qx = atom.pos[2*a2];
                            float qy = atom.pos[2*a2 + 1];

                            float dx = qx - px;
                            float dy = qy - py;
                            // float dist = sqrt(pow(dx, 2) + pow(dy, 2)); -> this is double slower
                            float dist = std::sqrt(dx*dx + dy*dy);
                            if (dist > min_dist || dist == 0.0)
                                continue;
                            float separate = 0.5 * (min_dist - dist) / dist;
                            dx *= separate;
                            dy *= separate;
                            
                            atom.pos[2*a1] -= dx;
                            atom.pos[2*a1 + 1] -= dy;
                            atom.pos[2*a2] += dx;
                            atom.pos[2*a2 + 1] += dy;
                        }
                    }
                }
            }
        }
    }

    void handle_atom_collisions(float min_dist) {
        for(int a=0;a<atom.cnt;a++) {
            if (atom.pos[a*2] <= min_dist+1 || atom.pos[a*2] >= width-min_dist-1) {
                atom.pos[a*2] = atom.pos[a*2]>width/2 ? width-min_dist-1 : min_dist+1;
                atom.vel[a*2] = 0;
            }
            if (atom.pos[a*2 + 1] <= min_dist+1 || atom.pos[a*2 + 1] >= height-min_dist-1) {
                atom.pos[a*2 + 1] = atom.pos[a*2 + 1]>height/2 ? height-min_dist-1 : min_dist+1;
                atom.vel[a*2 + 1] = 0;
            }
        }
    }

    void update_cell() {
        // reset cell types
        grid.reset_cell();

        // set cell that has atom as liquid
        for(int a=0;a<atom.cnt;a++) {
            int x_cell = floor(atom.pos[a*2]);
            int y_cell = floor(atom.pos[a*2 + 1]);
            if (x_cell < 0 || x_cell >= width || y_cell < 0 || y_cell >= height)
                continue;
            if (grid.cell[y_cell*width+x_cell] == Grid::CellType::SOLID)
                continue;

            grid.cell[y_cell*width+x_cell] = Grid::CellType::LIQUID;
        }
    }

    void transfer_velocity_to_grid() {
        // reset grid velocity components
        grid.reset_velocity();
        // grid.print_grid_status("\n\nreset grid component\n");
        
        // assume cell height is 1
        float offset[2][2] = {
            {0, 0.5},
            {0.5, 0}
        };

        // d points vel
        for(int d=0;d<2;d++) {
            grid.reset_r();

            for(int a=0;a<atom.cnt;a++) {
                // assume cell height is 1
                float x_p = atom.pos[a*2] - offset[d][0];
                float y_p = atom.pos[a*2 + 1] - offset[d][1];
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
                        
                        grid.vel[d][(y_cell+i)*(width+1-d)+x_cell+j] += weight[i][j] * atom.vel[a*2 + d];
                        grid.r[(y_cell+i)*(width+1)+x_cell+j] += weight[i][j];
                    }
                }
            }

            // grid.print_grid_status("\n\nbefore dividing with r\n");

            for(int i=0;i<height+d;i++) {
                for(int j=0;j<width+1-d;j++) {
                    if (grid.r[i*(width+1)+j] == 0)
                        continue;
                    grid.vel[d][i*(width+1-d)+j] /= grid.r[i*(width+1)+j];
                }
            }

            // grid.print_grid_status("\n\nafter dividing with r\n");
        }
    }

    void update_density() {
        reset_density();
        // assume cell height is 1
        float offset = 0.5;
        for(int a=0;a<atom.cnt;a++) {
            float x_p = atom.pos[a*2]-offset;
            float y_p = atom.pos[a*2 + 1]-offset;
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
                    density[(y_cell+i)*width+x_cell+j] += weight[i][j];
                }
            }
        }
    }

    void solve_incompressibility(int num_iter, float over_relaxation, float stiff, float rest_density) {
        for(int iter=0;iter<num_iter;iter++) {
            for(int i=0;i<height;i++) {
                for(int j=0;j<width;j++) {
                    if (grid.cell[i*width+j] != Grid::CellType::LIQUID)
                        continue;
                    float divergence = grid.vel[0][i*(width+1)+j+1] - grid.vel[0][i*(width+1)+j] + grid.vel[1][(i+1)*width+j] - grid.vel[1][i*width+j];
                    divergence *= over_relaxation;
                    if (density[i*width+j] > 0)
                        divergence -= stiff*(density[i*width+j] - rest_density);
                    int s1 = grid.s[i*width+j-1];
                    int s2 = grid.s[i*width+j+1];
                    int s3 = grid.s[(i-1)*width+j];
                    int s4 = grid.s[(i+1)*width+j];
                    int s = s1+s2+s3+s4;

                    if (s == 0)
                        continue;
                    grid.vel[0][i*(width+1)+j] += divergence * s1 / s;
                    grid.vel[0][i*(width+1)+j+1] -= divergence * s2 / s;
                    grid.vel[1][i*width+j] += divergence * s3 / s;
                    grid.vel[1][(i+1)*width+j] -= divergence * s4 / s;
                }
            }
        }
    }

    void transfer_velocity_to_atom() {
        // assume cell height is 1
        float offset[2][2] = {
            {0, 0.5},
            {0.5, 0}
        };

        // d points vel
        for(int d=0;d<2;d++) {    
            for(int a=0;a<atom.cnt;a++) {
                // assume cell height is 1
                float x_p = atom.pos[a*2] - offset[d][0];
                float y_p = atom.pos[a*2 + 1] - offset[d][1];
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
                        float q_curr = grid.vel[d][(y_cell+i)*(width+1-d)+x_cell+j];
                        float q_prev = grid.prev_vel[d][(y_cell+i)*(width+1-d)+x_cell+j];
                        // add difference of component to atom
                        numerator += (q_curr-q_prev)*weight[i][j];
                        denominator += weight[i][j];
                    }
                }
                float quantity = numerator / denominator;
                atom.vel[a*2 + d] += quantity;
            }
        }
    }

    void update() {
        // update atom velocity by gravity
        // update atom position by it's velocity
        // simulate particles
        simulate_atom(gravity, dt);

        handle_atom_collisions(min_dist);
        push_atoms_apart(min_dist, pushapart_iteration);

        // grid.print_grid_status("\n\ninit grid:\n");

        // update cell state
        update_cell();

        // transfer atom velocity to grid
        transfer_velocity_to_grid();
        // grid.print_grid_status("\n\natom->grid:\n");

        update_density();

        // save current grid velocity
        for(int d=0;d<2;d++)
            grid.vel[d].swap(grid.prev_vel[d]);
        
        // make fluid incompressible
        // using Gauss-Seidal method
        solve_incompressibility(incompress_iteration, over_relaxation, stiff, rest_density);
        // grid.print_grid_status("\n\nincompressible:\n");

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

    // int height=40, width=50;
    // int height = 12, width = 12;
    int height = 80, width = 120;

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
    for(int i=height*3/4;i<height*15/16;i++) {
        init_cell[i][width*3/4] = Grid::CellType::SOLID;
        init_cell[i][width/4] = Grid::CellType::SOLID;
    }
    for(int i=width/4;i<width*3/4;i++) {
        init_cell[height*15/16][i] = Grid::CellType::SOLID;
    }
    Grid grid = Grid(init_cell);

    // make half filled atoms
    std::vector<float> pos;
    std::vector<float> vel;
    for(int i=2;i<(height-2)*3/4;i++) {
        for(int j=2;j<width*3/4;j++) {
            float x = j;//+width/4;
            float y = height/4+i-2;
            float u = 0;
            float v = 0;
            pos.push_back(x);
            pos.push_back(y);
            vel.push_back(u);
            vel.push_back(v);
        }
    }
    Atom atom = Atom(pos, vel);

    // other configs
    float gravity = 50;
    float dt = 0.03;
    float min_dist = 0.5;
    int pushapart_iteration = 0;
    int incompress_iteration = 50;
    float over_relaxation = 1.9;
    float stiff = 2.0;
    float rest_density = 3.0;

    // make fluid
    FlipFluid fluid = FlipFluid(
        grid,
        atom,
        gravity,
        dt,
        min_dist,
        pushapart_iteration,
        incompress_iteration,
        over_relaxation,
        stiff,
        rest_density
    );

    while(true) {
        auto start = std::chrono::high_resolution_clock::now(); // start time
    
        fluid.update();
        fluid.print_fluid();
        // std::this_thread::sleep_for(std::chrono::milliseconds(2000));

        auto end = std::chrono::high_resolution_clock::now(); // end time
        std::chrono::duration<double> duration = end - start; // compute duration
        std::cout << "Time elapsed: " << duration.count() << '\n';
    }
}