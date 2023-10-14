/*********************
 *      INCLUDES
 *********************/

#include <math.h>

#include "fluidSim.h"
#include "drawing.h"

/*********************
 *      DEFINES
 *********************/

#define PARTICLE_RADIUS (0.02)

#define SCREEN_SCALING (100)

#define NUM_PARTICLES (200)

#define PARTICLE_SPACING (0.1)

#define SMOOTHING_RADIUS (1)

/**
 * 	If Set to 1, the starting positions will be randomized in the bounding box.
 * 	If set to 0, the starting positions will be in a grid in the center.
 */
#define RANDOMIZE_STARTING_POSITIONS (0)

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

static tVector2 get_random_direction();

static void update_densites();
static double calculate_density(tVector2 sample_point);
static double convert_density_to_pressure(double density);
static tVector2 calculate_pressure_force(int particle_index);
static double smoothing_kernel(double radius, float distance);
static double smoothing_kernel_derivative(double radius, double distance);

/**********************
 *  STATIC VARIABLES
 **********************/

const double gravity = 10;
const double collision_damping = 0.6;
static double target_density = 10;
static double pressure_multiplier = 1;

static tVector2 position[NUM_PARTICLES];
static tVector2 velocity[NUM_PARTICLES];
static double particle_density[NUM_PARTICLES];

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
	for(int i = 0; i < NUM_PARTICLES; i++){
		particle_density[i] = 0;
	}

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
	update_densites();

	for(int i = 0; i < NUM_PARTICLES; i++){
		tVector2 pressure_force = calculate_pressure_force(i);
		tVector2 pressure_acceleration = VECTOR2_SCALED(pressure_force, 1/particle_density[i]);
		velocity[i] = VECTOR2_SCALED(pressure_acceleration, TIME_DELTA_S);
	}

	for(int i = 0; i < NUM_PARTICLES; i++){
		position[i] = VECTOR2_ADD(position[i], velocity[i]);

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

static tVector2 get_random_direction(){
	double x = (double)((double)rand()) / RAND_MAX;
	double y = (double)((double)rand()) / RAND_MAX;
	tVector2 direction = VECTOR2_CREATE(x, y);
	return VECTOR2_SCALED(direction, VECTOR2_MAG(direction));
}

static void update_densites(){
	for(int i = 0; i < NUM_PARTICLES; i++){
		particle_density[i] = calculate_density(position[i]);
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

static double convert_density_to_pressure(double density){
	double density_error = density - target_density;
	double pressure = density_error * pressure_multiplier;
	return pressure;
}

static tVector2 calculate_pressure_force(int particle_index){
	tVector2 pressure_gradient = VECTOR2_ZERO;
	const double mass = 1;

	for(int other_particle_index = 0; other_particle_index < NUM_PARTICLES; other_particle_index++){
		if(other_particle_index == particle_index){
			continue;
		}
		tVector2 diff = VECTOR2_SUB(position[other_particle_index], position[particle_index]);
		double distance = VECTOR2_MAG(diff);
		tVector2 direction;
		if(distance == 0){
			direction = get_random_direction();
		}
		else{
			direction = VECTOR2_SCALED(VECTOR2_SUB(position[other_particle_index], position[particle_index]), 1/distance);
		}
		double slope = smoothing_kernel_derivative(distance, SMOOTHING_RADIUS);
		double density = calculate_density(position[other_particle_index]);
		double pressure_gradient_scaler = -convert_density_to_pressure(density) * slope * mass / density;

		pressure_gradient = VECTOR2_ADD(pressure_gradient, VECTOR2_SCALED(direction, -pressure_gradient_scaler));
	}

	return pressure_gradient;
}

static double smoothing_kernel(double radius, float distance){
	double value = NUM_MAX(0, radius * radius - distance * distance);
	const double volume = NUM_PI * pow(radius, 8) / 4;
	return value * value * value / volume;
}

static double smoothing_kernel_derivative(double radius, double distance){
	if(distance >= radius){
		return 0;
	}	

	float f = radius * radius - distance * distance;
	float scale = -24 / (NUM_PI * pow(radius, 8));
	return scale * distance * f * f;
}
