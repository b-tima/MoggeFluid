/*********************
 *      INCLUDES
 *********************/

#include "drawing.h"
#include <stdio.h>

/*********************
 *      DEFINES
 *********************/

/**********************
 *      TYPEDEFS
 **********************/

/**********************
 *  STATIC PROTOTYPES
 **********************/

static void set_color(tColor color);

static void set_color_a(tColor_a color);

static void draw_circle(int32_t centreX, int32_t centreY, int32_t radius);

/**********************
 *  STATIC VARIABLES
 **********************/

static SDL_Renderer* drawing_Renderer;

/**********************
 *      MACROS
 **********************/

/**********************
 *   GLOBAL FUNCTIONS
 **********************/

void drawing_init(){
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        printf("SDL_Init Error: %s\n", SDL_GetError());
        exit(1);
    }

    SDL_Window* window = SDL_CreateWindow("Hello SDL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, 0);
    if (window == NULL) {
        printf("SDL_CreateWindow Error: %s\n", SDL_GetError());
        exit(1);
    }

    drawing_Renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

	SDL_SetRenderDrawBlendMode(drawing_Renderer, SDL_BLENDMODE_BLEND);
}

void drawing_clearScreen(){
    SDL_SetRenderDrawColor(drawing_Renderer, BACKGROUND_COLOR.r, BACKGROUND_COLOR.g, BACKGROUND_COLOR.b, 255);
    SDL_RenderClear(drawing_Renderer);
}

void drawing_drawScreen(){
	SDL_RenderPresent(drawing_Renderer);
}

void drawing_delay(uint32_t delay){
	SDL_Delay(delay);
}

void drawing_drawCircle(tVector2_int pos, int32_t radius, tColor color){
	set_color(color);
	draw_circle(pos.x, pos.y, radius);
}

void drawing_drawPixel(tVector2_int pos, tColor_a color){
	set_color_a(color);
	SDL_RenderDrawPoint(drawing_Renderer, pos.x, pos.y);
}

void drawing_drawRect(tVector2_int start, tVector2_int end, tColor_a color){
	set_color_a(color);
	SDL_Rect rect;
	rect.h = end.y - start.y;
	rect.w = end.x - start.x;
	rect.x = start.x;
	rect.y = start.y;
	SDL_RenderFillRect(drawing_Renderer, &rect);
}

/**********************
 *   STATIC FUNCTIONS
 **********************/

static void set_color(tColor color){
	SDL_SetRenderDrawColor(drawing_Renderer, color.r, color.g, color.b, 255);
}

static void set_color_a(tColor_a color){
	SDL_SetRenderDrawColor(drawing_Renderer, color.r, color.g, color.b, color.a);
}

static void draw_circle(int32_t centreX, int32_t centreY, int32_t radius){
	for (int w = 0; w < radius * 2; w++)
    {
        for (int h = 0; h < radius * 2; h++)
        {
            int dx = radius - w; // horizontal offset
            int dy = radius - h; // vertical offset
            if ((dx*dx + dy*dy) <= (radius * radius))
            {
                SDL_RenderDrawPoint(drawing_Renderer, centreX + dx, centreY + dy);
            }
        }
    }
}