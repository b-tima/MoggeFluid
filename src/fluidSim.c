/*********************
 *      INCLUDES
 *********************/

#include <math.h>

#include "fluidSim.h"
#include "drawing.h"

/*********************
 *      DEFINES
 *********************/

#define PARTICLE_RADIUS (10)

#define NUM_PARTICLES (100)

#define PARTICLE_SPACING (20)

/**
 * 	If Set to 1, the starting positions will be randomized in the bounding box.
 * 	If set to 0, the starting positions will be in a grid in the center.
 */
#define RANDOMIZE_STARTING_POSITIONS (1)

/**
 * 	If set to 1, the simulation will not be started after initialiation.
 * 	If set to 0, the simulation will start
 */
#define DONT_START_SIMULARITON (0)

/**********************
 *      TYPEDEFS
 **********************/

/**********************
 *  STATIC PROTOTYPES
 **********************/

static tVector2 translate_position(tVector2 vector);
static void draw_particle(tVector2 vector);

static void resolve_edge_collision(tVector2* local_position, tVector2* local_velocity);

/**********************
 *  STATIC VARIABLES
 **********************/

static float gravity = 1000;
static float collision_damping = 0.6;

static tVector2 position[NUM_PARTICLES];
static tVector2 velocity[NUM_PARTICLES];

const tVector2 bounding_box = {
	.x = SCREEN_WIDTH,
	.y = SCREEN_HEIGHT
};

/**********************
 *      MACROS
 **********************/

/**********************
 *   GLOBAL FUNCTIONS
 **********************/

void fluidSim_init() {
#if !RANDOMIZE_STARTING_POSITIONS
	int particles_per_row = sqrt(NUM_PARTICLES);
	int particles_per_column = ((NUM_PARTICLES - 1) / particles_per_row) + 1;
	double spacing = PARTICLE_RADIUS * 2 + PARTICLE_SPACING;

	for(int i = 0; i < NUM_PARTICLES; i++){
		double x = (i % particles_per_row - particles_per_row / 2 + 0.5) * spacing;
		double y = (i / particles_per_row - particles_per_column / 2 + 0.5) * spacing;
		position[i].x = x;
		position[i].y = y;
	}
#else
	for(int i = 0; i < NUM_PARTICLES; i++){
		double x = bounding_box.x * (double)(((double)rand()) / RAND_MAX) - (bounding_box.x / 2);
		double y = bounding_box.y * (double)(((double)rand()) / RAND_MAX) - (bounding_box.y / 2);
		position[i].x = x;
		position[i].y = y;
	}
#endif
}

void fluidSim_update() {
#if !DONT_START_SIMULARITON
	for(int i = 0; i < NUM_PARTICLES; i++){
		velocity[i].y -= gravity * TIME_DELTA_S;
		position[i].y += velocity[i].y * TIME_DELTA_S;

		resolve_edge_collision(&position[i], &velocity[i]);

		draw_particle(position[i]);
	}
#else
	for(int i = 0; i < NUM_PARTICLES; i++){
		draw_particle(position[i]);
	}
#endif
}

/**********************
 *   STATIC FUNCTIONS
 **********************/

static tVector2 translate_position(tVector2 vector) {
    tVector2 translated_vector = {.x = vector.x + (SCREEN_WIDTH / 2),
                                  .y = (SCREEN_HEIGHT / 2) - vector.y};
    return translated_vector;
}

static void draw_particle(tVector2 vector) {
    tVector2 tranpos = translate_position(vector);
    drawing_drawCircle(VECTOR2_TO_INT(tranpos), PARTICLE_RADIUS,
                       COLOR_LIGHTBLUE);
}

static void resolve_edge_collision(tVector2* local_position, tVector2* local_velocity){
	tVector2 half_bounds = VECTOR2_SCALED(bounding_box, 0.5);

	if(abs(local_position->x) + PARTICLE_RADIUS > half_bounds.x){
		local_position->x = (half_bounds.x - PARTICLE_RADIUS) * NUM_SIGN(local_position->x);
		local_velocity->x *= -1 * collision_damping;
	}

	if(abs(local_position->y) + PARTICLE_RADIUS > half_bounds.y) {
		local_position->y = (half_bounds.y - PARTICLE_RADIUS) * NUM_SIGN(local_position->y);
		local_velocity->y *= -1 * collision_damping;
	}
}
