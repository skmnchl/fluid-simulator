#include <iostream>
#include <vector>
#include <cassert>

#include <chrono>

struct Atom {
    std::vector<float> pos; // {x1, y1, z1, x2, y2, z2, ...}
    std::vector<float> vel; // {u1, v1, w1, u2, v2, w2, ...}
    int cnt;

    Atom(std::vector<float> pos, std::vector<float> vel) {
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
    std::vector<float> vel[3]; // grid velocity
    std::vector<float> prev_vel[3]; // grid previous velocity
    std::vector<bool> s; // true if fluid can be filled (=false if cell type is solid)
    std::vector<float> r; // total sum of weights
    std::vector<bool> screen;

    Grid(std::vector<std::vector<std::vector<Grid::CellType>>> cell) {
        assert(cell.size()>0 && cell[0].size()>0 && cell[0][0].size()>0);

        // init grid size
        this->x_size = cell.size();
        this->y_size = cell[0].size();
        this->z_size = cell[0][0].size();
        this->cell.resize(x_size * y_size * z_size);
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
        //       z    w
        //       |    |
        //       +----|------+
        //      /|    |     /|
        //     / |         / |
        //    +-----------+  |
        //    |  |        | ------ v
        //    |  +--------|--+-- y
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
        this->r.resize((x_size+1) * (y_size+1) * (z_size+1));
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
        this->screen.resize(x_size * y_size * 2);
    }

    void print_grid() {
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

    void print_screen() {
        // init screen
        for(int i=0;i<screen.size();i++) {
            screen[i] = false;
        }

        // offset at screen
        int x_offset = x_size;
        int y_offset = y_size;

        for(int i=0;i<x_size;i++) {
            for(int j=0;j<y_size;j++) {
                for(int k=0;k<z_size;k++) {
                    if (cell[(i*y_size + j)*z_size + k] == Grid::CellType::LIQUID) {
                        float y = -0.5*(i+j) + k;
                        float x = -0.86602540378*(j-i);
                        int yi = y_offset + floor(y);
                        int xi = x_offset + floor(x);
                        screen[yi*2*z_size + xi] = true;
                    }
                }
            }
        }

        for(int i=0;i<x_size*2;i++) {
            for(int j=0;j<y_size*2;j++) {
                if (screen[i*2*x_size + j])
                    std::cout << "■ ";
                else
                    std::cout << "  ";
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
    Atom atom = Atom({},{});
    Grid grid = Grid({{{Grid::CellType::AIR}}});
    std::vector<int> num_cell_atom;
    std::vector<int> cell_atom_idx;
    std::vector<float> density; // grid density

    float gravity;
    int simulate_step;
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
        int simulate_step,
        float dt,
        float min_dist,
        int pushapart_iteration,
        int incompress_iteration,
        float over_relaxation,
        float stiff,
        float rest_density) {

        // init grid
        this->x_size = grid.x_size;
        this->y_size = grid.y_size;
        this->z_size = grid.z_size;
        this->grid = grid;
        
        // init atoms
        this->atom = atom;

        // init others
        this->num_cell_atom.resize(x_size * y_size * z_size + 1);
        this->cell_atom_idx.resize(atom.cnt);
        this->density.resize(x_size * y_size * z_size);
        this->gravity = gravity;
        this->simulate_step = simulate_step;
        this->dt = dt;
        this->min_dist = min_dist;
        this->pushapart_iteration = pushapart_iteration;
        this->incompress_iteration = incompress_iteration;
        this->over_relaxation = over_relaxation;
        this->stiff = stiff;
        this->rest_density = rest_density;
    }
    
    // void reset_cell_atom_idx() {
    //     for(int i=0;i<cell_atom_idx.size();i++)
    //         cell_atom_idx[i] = 0;
    // }

    void simulate_atom(float gravity, int step, float dt) {
        for(int a=0;a<atom.cnt;a++) {
            atom.vel[a*3 + 1] += gravity * dt; // add gravity

            float dx = atom.vel[a*3 + 0] * dt / step;
            float dy = atom.vel[a*3 + 1] * dt / step;
            float dz = atom.vel[a*3 + 2] * dt / step;
            float x = atom.pos[a*3 + 0];
            float y = atom.pos[a*3 + 1];
            float z = atom.pos[a*3 + 2];
            bool x_hit = false;
            bool y_hit = false;
            bool z_hit = false;
            for(int i=0;i<step;i++) {
                if (x<0 || x>=x_size 
                 || y<0 || y>=y_size
                 || z<0 || z>=z_size ) {
                    break;
                }

                if (!x_hit && grid.cell[(floor(x+dx)*y_size + floor(y))*z_size + floor(z)] != Grid::CellType::SOLID) {
                    x += dx;
                } else {
                    x_hit = true;
                }
                if (!y_hit && grid.cell[(floor(x)*y_size + floor(y+dy))*z_size + floor(z)] != Grid::CellType::SOLID) {
                    y += dy;
                } else {
                    y_hit = true;
                }
                if (!z_hit && grid.cell[(floor(x)*y_size + floor(y))*z_size + floor(z+dz)] != Grid::CellType::SOLID) {
                    z += dz;
                } else {
                    z_hit = true;
                }
            }

            atom.pos[a*3 + 0] = x;
            atom.pos[a*3 + 1] = y;
            atom.pos[a*3 + 2] = z;
            if (x_hit) { atom.vel[a*3 + 0] = 0; }
            if (y_hit) { atom.vel[a*3 + 1] = 0; }
            if (z_hit) { atom.vel[a*3 + 2] = 0; }
        }
    }

    // efficient collision detection using spacial hash
    void push_atoms_apart(float min_dist, int num_iter) {
        return; // TODO: add proper seperation function for 3d grid

        // // use grid method
        // // divide each atom into grid
        // for(int i=0;i<num_cell_atom.size();i++)
        //     num_cell_atom[i] = 0;
        // reset_cell_atom_idx();
        
        // // count number of atoms in each cell 
        // for(int a=0;a<atom.cnt;a++) {
        //     float xi = floor(atom.pos[a*2]);
        //     float yi = floor(atom.pos[a*2 + 1]);
        //     if (yi>=height || yi<0 || xi>=height || xi<0)
        //         continue;
        //     num_cell_atom[yi * width + xi] += 1;
        // }

        // // make partial sum of number of atoms in each cell
        // int atom_cnt = 0;
        // for(int i=0;i<num_cell_atom.size();i++) {
        //     atom_cnt += num_cell_atom[i];
        //     num_cell_atom[i] = atom_cnt;
        // }

        // // fill atom idx in seperate vector
        // for(int a=0;a<atom.cnt;a++) {
        //     float xi = floor(atom.pos[a*2]);
        //     float yi = floor(atom.pos[a*2 + 1]);
        //     if (yi>=height || yi<0 || xi>=height || xi<0)
        //         continue;
        //     num_cell_atom[yi * width + xi] -= 1;
        //     cell_atom_idx[num_cell_atom[yi * width + xi]] = a;
        // }

        // // push atoms apart
        // for(int iter=0;iter<num_iter;iter++) {
        //     for(int a1=0;a1<atom.cnt;a1++) {
        //         // consider near 9 cells
        //         float px = atom.pos[2*a1];
        //         float py = atom.pos[2*a1 + 1];
        //         int xi = floor(px);
        //         int yi = floor(py);
        //         int x0 = (xi-1>=0) ? xi-1 : 0;
        //         int x1 = (xi+1<width) ? xi+1 : width-1;
        //         int y0 = (yi-1>=0) ? yi-1 : 0;
        //         int y1 = (yi+1<height) ? yi+1 : height-1;

        //         for(int y=y0;y<=y1;y++) {
        //             for(int x=x0;x<=x1;x++) {
        //                 int cell_num = y * width + x;
        //                 int first = num_cell_atom[cell_num];
        //                 int last = num_cell_atom[cell_num + 1];
        //                 for (int atom_id = first; atom_id < last; atom_id++) {
        //                     int a2 = cell_atom_idx[atom_id];
        //                     if (a2 == a1)
        //                         continue;
        //                     float qx = atom.pos[2*a2];
        //                     float qy = atom.pos[2*a2 + 1];

        //                     float dx = qx - px;
        //                     float dy = qy - py;
        //                     // float dist = sqrt(pow(dx, 2) + pow(dy, 2)); -> this is double slower
        //                     float dist = std::sqrt(dx*dx + dy*dy);
        //                     if (dist > min_dist || dist == 0.0)
        //                         continue;
        //                     float separate = 0.5 * (min_dist - dist) / dist;
        //                     dx *= separate;
        //                     dy *= separate;
                            
        //                     atom.pos[2*a1] -= dx;
        //                     atom.pos[2*a1 + 1] -= dy;
        //                     atom.pos[2*a2] += dx;
        //                     atom.pos[2*a2 + 1] += dy;
        //                 }
        //             }
        //         }
        //     }
        // }
    }

    void handle_atom_collisions(float min_dist) {
        for(int a=0;a<atom.cnt;a++) {
            int dim_size[3] = {x_size, y_size, z_size};
            for(int d=0;d<3;d++) {
                if (atom.pos[a*3 + d] <= min_dist+1) {
                    atom.pos[a*3 + d] = min_dist+1;
                    atom.vel[a*3 + d] = 0;
                } else if (atom.pos[a*3 + d] >= dim_size[d]-min_dist-1) {
                    atom.pos[a*3 + d] = dim_size[d]-min_dist-1;
                    atom.vel[a*3 + d] = 0;
                }
            }
        }
    }

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

    void transfer_velocity_to_grid() {
        // reset grid velocity components
        for(int d=0;d<3;d++)
            for(int i=0;i<grid.vel[d].size();i++)
                grid.vel[d][i] = 0;
        // grid.print_grid_status("\n\nreset grid component\n");
        
        // assume cell height is 1
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
        for(int d=0;d<3;d++) {
            // reset r
            for(int i=0;i<grid.r.size();i++)
                grid.r[i] = 0;

            for(int a=0;a<atom.cnt;a++) {
                // assume cell height is 1
                float x_p = atom.pos[a*3 + 0] - offset[d][0];
                float y_p = atom.pos[a*3 + 1] - offset[d][1];
                float z_p = atom.pos[a*3 + 2] - offset[d][2];
                int x_cell = floor(x_p); // x direction cell number of atom
                int y_cell = floor(y_p); // y direction cell number of atom
                int z_cell = floor(z_p); // z direction cell number of atom
                float dx = x_p - x_cell; // atom x in particular cell
                float dy = y_p - y_cell; // atom y in particular cell
                float dz = z_p - z_cell; // atom z in particular cell
                
                // calculate weight for cell corner weight
                // each weight represent following corner
                //
                //       z
                //       | w[0][0][1]  w[0][1][1]
                //       +-----------+
                //      /|          /|
                //  w[1][0][1]     w[1][1][1]
                //    +-----------+  |
                //    | w[0][0][0]|  | w[0][1][0]
                //    |  +--------|--+-- y
                //    | /         | /
                //    |/          |/
                //    +-----------+
                //   / w[1][0][0]  w[1][1][0]
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
                            
                            grid.vel[d][((x_cell+i)*dim_size[d][1] + y_cell+j)*dim_size[d][2] + z_cell+k] += weight[i][j][k] * atom.vel[a*3 + d];
                            grid.r[((x_cell+i)*(y_size+1) + y_cell+j)*(z_size+1) + z_cell+k] += weight[i][j][k];
                        }
                    }
                }
            }

            // grid.print_grid_status("\n\nbefore dividing with r\n");

            for(int i=0;i<dim_size[d][0];i++) {
                for(int j=0;j<dim_size[d][1];j++) {
                    for(int k=0;k<dim_size[d][2];k++) {
                        float weight_sum = grid.r[(i*(y_size+1) + j)*(z_size+1) + k];
                        if (weight_sum == 0.0)
                            continue;
                        grid.vel[d][(i*dim_size[d][1] + j)*dim_size[d][2] + k] /= weight_sum;
                    }
                }
            }

            // grid.print_grid_status("\n\nafter dividing with r\n");
        }
    }

    void update_density() {
        // reset density
        for(int i=0;i<density.size();i++)
            density[i] = 0;

        // assume cell height is 1
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

    void solve_incompressibility(int num_iter, float over_relaxation, float stiff, float rest_density) {
        for(int iter=0;iter<num_iter;iter++) {
            for(int i=0;i<x_size;i++) {
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
            }
        }
    }

    void transfer_velocity_to_atom() {
        // assume cell height is 1
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
        for(int d=0;d<3;d++) {    
            for(int a=0;a<atom.cnt;a++) {
                // assume cell height is 1
                float x_p = atom.pos[a*3 + 0] - offset[d][0];
                float y_p = atom.pos[a*3 + 1] - offset[d][1];
                float z_p = atom.pos[a*3 + 2] - offset[d][2];
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
                            int idx = ((x_cell+i)*dim_size[d][1] + y_cell+j)*dim_size[d][2] + z_cell+k;
                            // add difference of component to atom
                            numerator += (grid.vel[d][idx]-grid.prev_vel[d][idx])*weight[i][j][k];
                            denominator += weight[i][j][k];
                        }
                    }
                }
                float quantity = numerator / denominator;
                atom.vel[a*3 + d] += quantity;
            }
        }
    }

    void update() {
        // update atom velocity by gravity
        // update atom position by it's velocity
        // simulate particles
        simulate_atom(gravity, simulate_step, dt);

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
        for(int d=0;d<3;d++)
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
        // grid.print_screen();
    }
};

int main() {
    // faster I/O
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    int x_size = 30, y_size = 40, z_size = 50;

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
    int simulate_step = 10;
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
        simulate_step,
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
        // std::this_thread::sleep_for(std::chrono::milliseconds(20));

        auto end = std::chrono::high_resolution_clock::now(); // end time
        std::chrono::duration<double> duration = end - start; // compute duration
        std::cout << "Time elapsed: " << duration.count() << '\n';
    }
}