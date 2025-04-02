#include "FluidRenderer.h"

FluidRenderer::FluidRenderer(FlipFluid *fluid, int winWidth, int winHeight) {
    this->fluid = fluid;
    this->vertexData.resize(fluid->atom.cnt * 6);
    this->winWidth = winWidth;
    this->winHeight = winHeight;

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
    this->pWindow = pWindow;
}

GLuint FluidRenderer::initShader()
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

void FluidRenderer::initGeometry(GLuint shaderProgram)
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
void FluidRendererMainLoop(void* RendererInstance)
{
    // get arguments from renderer instance
    FluidRenderer *renderer = static_cast<FluidRenderer*>(RendererInstance);
    FlipFluid *fluid = renderer->fluid;
    int winWidth = renderer->winWidth;
    int winHeight = renderer->winHeight;
    float &azimuthalAngle = renderer->azimuthalAngle;
    float &polarAngle = renderer->polarAngle;
    float &scale = renderer->scale;
    SDL_Window *pWindow = renderer->pWindow;

    auto updateStart = std::chrono::high_resolution_clock::now(); // start time
        
    // update fluid status
    fluid->update();

    // rotate the fluid
    azimuthalAngle += 0.008f;

    auto updateEnd = std::chrono::high_resolution_clock::now(); // end time
    std::chrono::duration<double> updateDuration = updateEnd - updateStart; // compute duration


    auto drawStart = std::chrono::high_resolution_clock::now(); // start time

    // Resize on every frame for brevity, normally resizing is done on resize event
    SDL_GL_GetDrawableSize(pWindow, &winWidth, &winHeight);
    glViewport(0, 0, winWidth, winHeight);   

    // Clear screen
    glClear(GL_COLOR_BUFFER_BIT);

    // Update vertex buffer
    // Interleaved vertex position and color data
    const float PI = 3.1416; 
    // polarAngle = cos(azimuthalAngle) / 3;
    for(int a=0;a<fluid->atom.cnt;a++) {
        // calculate view projection
        float i = fluid->atom.pos[a*3 + 0] - fluid->xSize/2;
        float j = fluid->atom.pos[a*3 + 1] - fluid->ySize/2;
        float k = fluid->atom.pos[a*3 + 2] - fluid->zSize/2;
        
        float i_y = cos(azimuthalAngle)*i - sin(azimuthalAngle)*k;
        float j_y = -j;
        float k_y = sin(azimuthalAngle)*i + cos(azimuthalAngle)*k;

        float i_yz = i_y;
        float j_yz = cos(polarAngle)*j_y - sin(polarAngle)*k_y;
        float k_yz = sin(polarAngle)*j_y + cos(polarAngle)*k_y;
        k_yz += fluid->zSize;

        if (k == 0)
            continue;
        float x = i_yz * scale / k_yz;
        float y = j_yz * scale / k_yz;

        // set color for each atom
        float cameraX = -fluid->xSize * 0.1;
        float cameraY = -fluid->ySize * 0.1;
        float cameraZ = -fluid->zSize * 0.1;
        float depth = sqrt((cameraX - i)*(cameraX - i)
                         + (cameraY - j)*(cameraY - j)
                         + (cameraZ - k)*(cameraZ - k));
        float vSquare = sqrt((fluid->atom.vel[a*3 + 0])*(fluid->atom.vel[a*3 + 0])
                           + (fluid->atom.vel[a*3 + 1])*(fluid->atom.vel[a*3 + 1])
                           + (fluid->atom.vel[a*3 + 2])*(fluid->atom.vel[a*3 + 2]));
        float r = 1.0f - depth/50 + vSquare/50;
        float g = 1.2f - depth/50;
        float b = 1.5f - depth/50;
        r = (r>1.0f) ? 1.0f : r; r = (r<0.0f) ? 0.0f : r;
        g = (g>1.0f) ? 1.0f : g; g = (g<0.0f) ? 0.0f : g;
        b = (b>1.0f) ? 1.0f : b; b = (b<0.0f) ? 0.0f : b;

        int maxSize = fluid->xSize;
        maxSize = maxSize > fluid->ySize ? maxSize : fluid->ySize;
        maxSize = maxSize > fluid->zSize ? maxSize : fluid->zSize;
        renderer->vertexData[a*6 + 0] = x / maxSize;
        renderer->vertexData[a*6 + 1] = y / maxSize;
        renderer->vertexData[a*6 + 2] = 0.0f;
        renderer->vertexData[a*6 + 3] = r;
        renderer->vertexData[a*6 + 4] = g;
        renderer->vertexData[a*6 + 5] = b;
    }
    glBufferData(GL_ARRAY_BUFFER, renderer->vertexData.size() * sizeof(GLfloat), renderer->vertexData.data(), GL_STREAM_DRAW);

    // Draw the vertex buffer
    glDrawArrays(GL_POINTS, 0, fluid->atom.cnt);

    // Swap front/back framebuffers
    SDL_GL_SwapWindow(pWindow);

    auto drawEnd = std::chrono::high_resolution_clock::now(); // end time
    std::chrono::duration<double> drawDuration = drawEnd - drawStart; // compute duration

    // output renderer status
    EM_ASM({
        document.getElementById("update-duration").innerText = "Update duration: "+$0.toFixed(5)+"s";
        document.getElementById("draw-duration").innerText = "Draw duration: "+$1.toFixed(5)+"s";
        document.getElementById("azimuthal-angle").innerText = "Azimuthal Angle: "+$2.toFixed(2)+"°";
        document.getElementById("polar-angle").innerText = "Polar Angle: "+$3.toFixed(2)+"°";
        document.getElementById("scale").innerText = "Scale: "+$4.toFixed(2);
    }, 
    updateDuration.count(),
    drawDuration.count(),
    azimuthalAngle/PI*180,
    polarAngle/PI*180,
    scale
    );

    // output fluid status
    EM_ASM({
        document.getElementById("grid-size").innerText = "Grid Size: "+$0+" * "+$1+" * "+$2;
        document.getElementById("num-particles").innerText = "Particle Count: "+$3;
        document.getElementById("gravity").innerText = "Gravity: "+$4;
        document.getElementById("dt").innerText = "dt: "+$5.toFixed(2);
        document.getElementById("simulation-step").innerText = "Particle Simulation Step: "+$6;
        document.getElementById("min-dist").innerText = "Particle Minum Distance: "+$7;
        document.getElementById("push-apart-iteration").innerText = "Push Apart Iteration: "+$8;
        document.getElementById("incompress-iteration").innerText = "Solve Incompressibility Iteration: "+$9;
        document.getElementById("over-relaxation").innerText = "Over Relaxation: "+$10.toFixed(2);
        document.getElementById("stiff").innerText = "Stiffness: "+$11;
        document.getElementById("rest-density").innerText = "Rest Density: "+$12;
        document.getElementById("num-threads").innerText = "Threads Count: "+$13;
    },
    fluid->xSize,
    fluid->ySize,
    fluid->zSize,
    fluid->atom.cnt,
    fluid->gravity,
    fluid->dt,
    fluid->simulateStep,
    fluid->minDist,
    fluid->pushApartIteration,
    fluid->incompressIteration,
    fluid->overRelaxation,
    fluid->stiff,
    fluid->restDensity,
    fluid->numThreads
    );
}