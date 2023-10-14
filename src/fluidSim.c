/*********************
 *      INCLUDES
 *********************/

#include <math.h>

#include "fluidSim.h"
#include "drawing.h"

/*********************
 *      DEFINES
 *********************/

#define PARTICLE_RADIUS (0.1)

#define SCREEN_SCALING (100)

#define NUM_PARTICLES (100)

#define PARTICLE_SPACING (1)

#define SMOOTHING_RADIUS (1)

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

/**
 *	If set to 1, the simulation will draw the density
 * 	If set to 0, the simulation will not draw the density
 */
#define DRAW_DENSITY (0)

/**********************
 *      TYPEDEFS
 **********************/

/**********************
 *  STATIC PROTOTYPES
 **********************/

static tVector2 translate_global_position(tVector2_int vector);
static tVector2 translate_local_position(tVector2 vector);
static void draw_particle(tVector2 vector);

static void resolve_edge_collision(tVector2* local_position, tVector2* local_velocity);

static double calculate_density(tVector2 sample_point);
static double smoothing_kernel(double radius, float distance);

/**********************
 *  STATIC VARIABLES
 **********************/

static float gravity = 10;
static float collision_damping = 0.6;

static tVector2 position[NUM_PARTICLES];
static tVector2 velocity[NUM_PARTICLES];

const tVector2 bounding_box = {
	.x = SCREEN_WIDTH / SCREEN_SCALING,
	.y = SCREEN_HEIGHT / SCREEN_SCALING
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
#if DRAW_DENSITY
	// Calculate density field...
	for(int x = 0; x < SCREEN_WIDTH; x++){
		for(int y = 0; y < SCREEN_HEIGHT; y++){
			tVector2_int global = {
				.x = x,
				.y = y
			};

			tVector2 pos = translate_global_position(global);

			double density = calculate_density(pos);
			
			double a = 255 * NUM_MIN(1, density / 3);
			tColor_a color = COLOR_WITH_A(COLOR_LIGHTBLUE, (uint8_t)a);
			if (a < 20){
				continue;
			}
			drawing_drawPixel(global, color);
		}
	}
#endif

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

		tVector2 sample_point = {
			.x = 0,
			.y = 0
		};
	}
#endif
}

void fluidSim_onClick(tVector2_int pos){
	tVector2 local_pos = translate_global_position(pos);
	double density = calculate_density(local_pos);
	printf("density = %f\n", density);
}

/**********************
 *   STATIC FUNCTIONS
 **********************/

static tVector2 translate_global_position(tVector2_int vector){
	tVector2 translated_vector = {
		.x = ((double)vector.x) / SCREEN_SCALING - bounding_box.x / 2,
		.y = -((double)vector.y) / SCREEN_SCALING + bounding_box.y / 2
	};

	return translated_vector;
}

static tVector2 translate_local_position(tVector2 vector) {
    tVector2 translated_vector = {.x = vector.x * SCREEN_SCALING + (SCREEN_WIDTH / 2),
                                  .y = (SCREEN_HEIGHT / 2) - vector.y * SCREEN_SCALING};
    return translated_vector;
}

static void draw_particle(tVector2 vector) {
    tVector2 tranpos = translate_local_position(vector);
    drawing_drawCircle(VECTOR2_TO_INT(tranpos), PARTICLE_RADIUS * SCREEN_SCALING,
                       COLOR_LIGHTBLUE);
}

static void resolve_edge_collision(tVector2* local_position, tVector2* local_velocity){
	tVector2 half_bounds = VECTOR2_SCALED(bounding_box, 0.5);

	if(NUM_ABS(local_position->x) + PARTICLE_RADIUS > half_bounds.x){
		local_position->x = (half_bounds.x + PARTICLE_RADIUS) * NUM_SIGN(local_position->x);
		local_velocity->x *= -1 * collision_damping;
	}

	if(NUM_ABS(local_position->y) + PARTICLE_RADIUS > half_bounds.y) {
		local_position->y = (half_bounds.y - PARTICLE_RADIUS) * NUM_SIGN(local_position->y);
		local_velocity->y *= -1 * collision_damping;
	}
}

static double calculate_density(tVector2 sample_point){
	double density = 0;
	const double mass = 1;

	for(int i = 0; i < NUM_PARTICLES; i++){
		tVector2 diff = VECTOR2_SUB(position[i], sample_point);
		double distance = VECTOR2_MAG(diff);

		double influence = smoothing_kernel(SMOOTHING_RADIUS, distance);
		density += mass * influence;
	}

	return density;
}

static double smoothing_kernel(double radius, float distance){
	double value = NUM_MAX(0, radius * radius - distance * distance);
	const double volume = NUM_PI * pow(radius, 8) / 4;
	return value * value * value / volume;
}
