#pragma once 

#include <vector>
#include <chrono>
#include <math.h>

#include <emscripten.h>
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

    // others
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