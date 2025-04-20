#pragma once 

#include <vector>
#include <chrono>
#include <math.h>

#include <emscripten.h>
#include <emscripten/html5.h>
#include <SDL.h>
#include <SDL_opengles2.h>

#include "FlipFluid.h"

class FluidRenderer {
    public:
    // fluid to render
    FlipFluid *fluid;

    // vertexData that contains position and rgb values of each atom
    std::vector<GLfloat> vertexData;

    // window
    SDL_Window* pWindow;
    int winWidth;
    int winHeight;

    // mouse event
    bool isMouseDown = false;
    float mouseX = 0.0f;
    float mouseY = 0.0f;
    float moveX = 0.0f;
    float moveY = 0.0f;

    // user options
    bool autoRotate = true;
    bool playSimulation = true;
    int gravityMode = 0;

    // view
    float azimuthalAngle = 0.0f;
    float polarAngle = 0.1f;
    float scale = 50.0f;
    
    FluidRenderer(FlipFluid *fluid, int winWidth, int winHeight);

    // Vertex shader
    const GLchar* vertexSource =
        "attribute vec4 position;                      \n"
        "attribute vec3 vertexColor;                   \n"
        "varying vec3 color;                           \n"
        "void main()                                   \n"
        "{                                             \n"
        "    gl_PointSize = 2.0;                       \n"
        "    gl_Position = vec4(position.xyz, 1.0);    \n"
        "    color = vertexColor;                      \n"
        "}                                             \n";

    // Fragment/pixel shader
    const GLchar* fragmentSource =
        "precision mediump float;                      \n"
        "varying vec3 color;                           \n"
        "void main()                                   \n"
        "{                                             \n"
        "    gl_FragColor = vec4 ( color, 0.8 );       \n"
        "}                                             \n";

    // inits shader using vertexSource and fragmentSource
    GLuint initShader();

    // inits fluid geometry
    void initGeometry(GLuint shaderProgram);

};

// main animation loop for emscripten_main_loop_arg
void FluidRendererMainLoop(void* RendererInstance);

// mouse event callback function
bool FluidRendererMouseCallback(int eventType, const EmscriptenMouseEvent *e, void* RendererInstance);

// touch event callback function
bool FluidRendererTouchCallback(int eventType, const EmscriptenTouchEvent *e, void* RendererInstance);