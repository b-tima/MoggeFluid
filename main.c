
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
		while(SDL_PollEvent(&event)){
			switch(event.type){
				case SDL_QUIT:
					running = false;
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