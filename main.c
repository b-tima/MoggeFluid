
#include <SDL.h>
#include "drawing.h"
#include "fluidSim.h"
#include "common.h"
#include <stdbool.h>

static bool running = true;

int main() {
	drawing_init();
	fluidSim_init();

	while(running) {
		SDL_Event event;
		bool mouse_hyst = false;
		while(SDL_PollEvent(&event)){
			switch(event.type){
				case SDL_QUIT:
					running = false;
					break;
				case SDL_MOUSEBUTTONDOWN:
					if(mouse_hyst){
						break;
					}
					mouse_hyst = true;
					tVector2_int mouse_pos;
					mouse_pos.x = event.button.x;
					mouse_pos.y = event.button.y;
					fluidSim_onClick(mouse_pos);
					break;
				case SDL_MOUSEBUTTONUP:
					mouse_hyst = false;
					break;
				default:
					break;
			}
		}

		drawing_clearScreen();

		fluidSim_update();

		drawing_drawScreen();

		drawing_delay(TIME_DELTA);
	}
	
    return 0;
}