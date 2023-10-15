
#include <SDL.h>
#include "drawing.h"
#include "fluidSim.h"
#include "common.h"
#include <stdbool.h>

static bool running = true;

int main() {
	drawing_init();
	fluidSim_init();

	bool mouse_hyst = false;
	while(running) {
		tVector2_int mouse_pos;

		SDL_Event event;
		while(SDL_PollEvent(&event)){
			switch(event.type){
				case SDL_QUIT:
					running = false;
					break;
				case SDL_MOUSEBUTTONDOWN:
					mouse_hyst = true;
					mouse_pos.x = event.button.x;
					mouse_pos.y = event.button.y;
					break;
				case SDL_MOUSEMOTION:
					if(!mouse_hyst){
						break;
					}
					mouse_pos.x = event.button.x;
					mouse_pos.y = event.button.y;
					break;
				case SDL_MOUSEBUTTONUP:
					mouse_hyst = false;
					break;
				default:
					break;
			}
		}

		if(mouse_hyst){
			fluidSim_onClick(mouse_pos);
		}

		drawing_clearScreen();

		fluidSim_update();

		drawing_drawScreen();

		drawing_delay(TIME_DELTA);
	}
	
	fluidSim_destroy();

    return 0;
}