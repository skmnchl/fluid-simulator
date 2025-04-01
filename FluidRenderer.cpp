//
// Setup:
//     Install emscripten: https://emscripten.org/docs/getting_started/downloads.html
//
// Build:
//     WASM : emcc -pthread --bind main.cpp flip_fluid.cpp -sINITIAL_MEMORY=64MB -sPTHREAD_POOL_SIZE=5 -sUSE_SDL=2 -sFULL_ES2=1 -o index.html --shell-file shell.html
//     Console : g++ main.cpp flip_fluid.cpp -o main
//
// Run:
//     WASM: Run a http server serving index.html with proper COOP and COEP headers (Required for multi threading)
//     or just run the command below and check localhost:8888
//         emrun --no_browser --port 8888 --serve_after_close index.html
//     Console:
//         ./main
//

#include <chrono>

#ifdef __EMSCRIPTEN__
#include <cstdio> // usde for hex color manipulation
#include <emscripten.h>
#include <emscripten/val.h>
#include <SDL.h>
#include <SDL_opengles2.h>
#endif

#include "flip_fluid.h"

// using global variables for brevity
FlipFluid *fluid = nullptr;

#ifdef __EMSCRIPTEN__
std::vector<GLfloat> vertexData;

// Vertex shader
const GLchar* vertexSource =
    "attribute vec4 position;                      \n"
    "attribute vec3 vertexColor;                   \n"
    "varying vec3 color;                           \n"
    "void main()                                   \n"
    "{                                             \n"
    "    gl_PointSize = 3.0;                       \n"
    "    gl_Position = vec4(position.xyz, 1.0);    \n"
    "    color = vertexColor;                      \n"
    "}                                             \n";

// Fragment/pixel shader
const GLchar* fragmentSource =
    "precision mediump float;                      \n"
    "varying vec3 color;                           \n"
    "void main()                                   \n"
    "{                                             \n"
    "    vec2 coord = gl_PointCoord - vec2(0.5);   \n"
    "    float dist = length(coord) * 2.0;         \n"
    "                                              \n"
    "    float alpha = smoothstep(1.0, 0.5, dist); \n"
    "    gl_FragColor = vec4 ( color, alpha );     \n"
    "}                                             \n";

GLuint initShader()
{
    // Enable blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Default blending function

    // Create and compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);

    // Create and compile fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    glCompileShader(fragmentShader);

    // Link vertex and fragment shader into shader program and use it
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    return shaderProgram;
}

void initGeometry(GLuint shaderProgram)
{
    GLuint vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Interleaved vertex position and color data
    for(int a=0;a<fluid->atom.cnt * 6;a++) {
        vertexData[a] = 0.0f;
    }

    glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(GLfloat), vertexData.data(), GL_STREAM_DRAW);

    // Link position attribute
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);

    // Link color attribute
    GLint colorAttrib = glGetAttribLocation(shaderProgram, "vertexColor");
    glEnableVertexAttribArray(colorAttrib);
    glVertexAttribPointer(colorAttrib, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
}

// udpate fluid and draw to "canvas" using WebGL
void mainLoop(void* mainLoopArg) 
{
    // update fluid status
    fluid->update();

    // Resize on every frame for brevity, normally resizing is done on resize event
    SDL_Window* pWindow = (SDL_Window*)mainLoopArg;
    int winWidth, winHeight;
    SDL_GL_GetDrawableSize(pWindow, &winWidth, &winHeight);
    glViewport(0, 0, winWidth, winHeight);   

    // Clear screen
    glClear(GL_COLOR_BUFFER_BIT);

    // Update vertex buffer
    // Interleaved vertex position and color data
    float scale = 1.0f;
    for(int a=0;a<fluid->atom.cnt;a++) {
        // calculate view projection
        float i = fluid->atom.pos[a*3 + 0];
        float j = fluid->atom.pos[a*3 + 1];
        float k = fluid->atom.pos[a*3 + 2];                
        float x = -0.86602540378*(k-i); // sqrt(3)/2 * (k-i)
        float y = 0.5*(i+k) - j;
        x = x * scale;
        y = y * scale;

        // set color for each atom
        float cameraX = -fluid->xSize * 0.1;
        float cameraY = -fluid->ySize * 0.1;
        float cameraZ = -fluid->zSize * 0.1;
        float depth = sqrt((cameraX - fluid->atom.pos[a*3 + 0])*(cameraX - fluid->atom.pos[a*3 + 0])
                         + (cameraY - fluid->atom.pos[a*3 + 1])*(cameraY - fluid->atom.pos[a*3 + 1])
                         + (cameraZ - fluid->atom.pos[a*3 + 2])*(cameraZ - fluid->atom.pos[a*3 + 2]));
        float vSquare = sqrt((fluid->atom.vel[a*3 + 0])*(fluid->atom.vel[a*3 + 0])
                           + (fluid->atom.vel[a*3 + 1])*(fluid->atom.vel[a*3 + 1])
                           + (fluid->atom.vel[a*3 + 2])*(fluid->atom.vel[a*3 + 2]));
        float r = 1.5f - depth/50 + vSquare/50;
        float g = 1.7f - depth/50;
        float b = 2.0f - depth/50;
        r = (r>1.0f) ? 1.0f : r; r = (r<0.0f) ? 0.0f : r;
        g = (g>1.0f) ? 1.0f : g; g = (g<0.0f) ? 0.0f : g;
        b = (b>1.0f) ? 1.0f : b; b = (b<0.0f) ? 0.0f : b;

        int maxSize = fluid->xSize;
        maxSize = maxSize > fluid->ySize ? maxSize : fluid->ySize;
        maxSize = maxSize > fluid->zSize ? maxSize : fluid->zSize;
        vertexData[a*6 + 0] = x / maxSize;
        vertexData[a*6 + 1] = y / maxSize;
        vertexData[a*6 + 2] = 0.0f;
        vertexData[a*6 + 3] = r;
        vertexData[a*6 + 4] = g;
        vertexData[a*6 + 5] = b;
    }
    glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(GLfloat), vertexData.data(), GL_STREAM_DRAW);

    // Draw the vertex buffer
    glDrawArrays(GL_POINTS, 0, fluid->atom.cnt);

    // Swap front/back framebuffers
    SDL_GL_SwapWindow(pWindow);
}
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
    for(float i=2;i<xSize-2;i+=0.5) {
        for(float j=ySize/2;j<ySize-2;j+=0.5) {
            for(float k=zSize*2/3;k<zSize-2;k+=0.5) {
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
    vertexData.resize(fluid->atom.cnt * 6);
    int winWidth = 768, winHeight = 768;

    // Create SDL window
    SDL_Window *pWindow = 
        SDL_CreateWindow("Flip Fluid", 
                         SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                         winWidth, winHeight, 
                         SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_SHOWN);

    // Create OpenGLES 2 context on SDL window
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
    SDL_GL_SetSwapInterval(1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GLContext glc = SDL_GL_CreateContext(pWindow);

    // Set clear color to black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // Get actual GL window size in pixels, in case of high dpi scaling
    SDL_GL_GetDrawableSize(pWindow, &winWidth, &winHeight);
    glViewport(0, 0, winWidth, winHeight);   

    // Initialize shader and geometry
    GLuint shaderProgram = initShader();
    initGeometry(shaderProgram);

    // Start the main loop
    void* mainLoopArg = pWindow;

    emscripten_set_main_loop_arg(mainLoop, mainLoopArg, 0, true);
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