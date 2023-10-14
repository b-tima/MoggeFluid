/*********************
 *      INCLUDES
 *********************/

#include <math.h>
#include <stdlib.h>

#include "fluidSim.h"
#include "drawing.h"

/*********************
 *      DEFINES
 *********************/

#define PARTICLE_RADIUS (0.04)

#define NUM_PARTICLES (1000)

#define PARTICLE_SPACING (0.1)

#define SMOOTHING_RADIUS (1)

#define SCREEN_SCALING (100)

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

typedef struct{
	uint32_t key;
	int index;
} cell_key_pair;

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
static double calculate_shared_pressure(double d1, double d2);
static tVector2 calculate_pressure_force(int particle_index);
static double smoothing_kernel(double radius, float distance);
static double smoothing_kernel_derivative(double radius, double distance);

static int get_cell_key_from_particle(int particle_index);
static int get_cell_key_from_grid_index(tVector2_int grid);
static tVector2_int position_to_cell_coordinate(int particle_index);
static int cell_sorting_compare(const void* c1, const void* c2);
static void update_cell_keys();

/**********************
 *  STATIC VARIABLES
 **********************/

const double gravity = 10;
const double collision_damping = 0.6;
const double target_density = 1;
const double pressure_multiplier = 1;
const double mass = 1;

static tVector2 position[NUM_PARTICLES];
static tVector2 velocity[NUM_PARTICLES];
static double particle_density[NUM_PARTICLES];

const tVector2 bounding_box = {
	.x = SCREEN_WIDTH / SCREEN_SCALING,
	.y = SCREEN_HEIGHT / SCREEN_SCALING
};

static cell_key_pair* spacial_LUT;
static int spacial_LUT_rows;
static int spacial_LUT_cols;
static int num_spacial_LUT;

static int* spacial_LUT_start_indeces;

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

	spacial_LUT_rows = ceil(bounding_box.x / SMOOTHING_RADIUS);
	spacial_LUT_cols = ceil(bounding_box.y / SMOOTHING_RADIUS);
	num_spacial_LUT = spacial_LUT_rows * spacial_LUT_cols;

	spacial_LUT = (cell_key_pair*) malloc(NUM_PARTICLES * sizeof(cell_key_pair));
	spacial_LUT_start_indeces = (int*)malloc(num_spacial_LUT * sizeof(int));
 
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
	for(int x = 0; x < SCREEN_WIDTH / 10; x++){
		for(int y = 0; y < SCREEN_HEIGHT / 10; y++){
			tVector2_int global = {
				.x = 10*x + 5,
				.y = 10*y + 5
			};

			tVector2 pos = translate_global_position(global);

			double density = calculate_density(pos)/10000;

			if(density < 0){
				density = 0;
			}

			double a = 255 * NUM_MIN(1, density);
			tColor_a color = COLOR_WITH_A(COLOR_LIGHTBLUE, (uint8_t)a);

			tVector2_int start = {
				.x = global.x - 5,
				.y = global.y - 5
			};
			tVector2_int end = {
				.x = global.x + 5,
				.y = global.y + 5
			};
			drawing_drawRect(start, end, color);
		}
	}
#endif

#if !DONT_START_SIMULARITON
	update_cell_keys();
	update_densites();

	for(int i = 0; i < NUM_PARTICLES; i++){
		tVector2 pressure_force = calculate_pressure_force(i);
		tVector2 pressure_acceleration = VECTOR2_SCALED(pressure_force, 1/particle_density[i]);
		velocity[i] = VECTOR2_ADD(velocity[i], VECTOR2_SCALED(pressure_acceleration, TIME_DELTA_S));
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
	double pressure = convert_density_to_pressure(calculate_density(local_pos));
	printf("pressure = %f\n", pressure);
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
		local_position->x = (half_bounds.x - PARTICLE_RADIUS) * NUM_SIGN(local_position->x);
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

static double calculate_shared_pressure(double d1, double d2){
	double pressure1 = convert_density_to_pressure(d1);
	double pressure2 = convert_density_to_pressure(d2);
	return (pressure1 + pressure2) / 2;
}

static tVector2 calculate_pressure_force(int particle_index){
	tVector2_int particle_cell_coord = position_to_cell_coordinate(particle_index);
	double sqr_radius = SMOOTHING_RADIUS * SMOOTHING_RADIUS;

	tVector2 pressure_gradient = VECTOR2_ZERO;

	for(int x = -1; x < 2; x++){
		for(int y = -1; y < 2; y++){
			tVector2_int grid = (tVector2_int){
				.x = particle_cell_coord.x + x,
				.y = particle_cell_coord.y + y
			};
			int key = get_cell_key_from_grid_index(grid);
			int cell_start_index = spacial_LUT_start_indeces[key];

			if(cell_start_index < 0){
				continue;
			}

			for(int i = cell_start_index; i < NUM_PARTICLES; i++){
				if(spacial_LUT[i].key != key){
					break;
				}

				if(spacial_LUT[i].index == particle_index){
					continue;
				}

				int other_particle_index = spacial_LUT[i].index;
				tVector2 diff = VECTOR2_SUB(position[other_particle_index], position[particle_index]);
				double sqr_distance = VECTOR2_SQR_MAG(VECTOR2_SUB(position[other_particle_index], position[particle_index]));

				if(sqr_distance < sqr_radius){
					double distance = sqrt(sqr_distance);
					tVector2 direction;
					if(__glibc_unlikely(distance == 0)){
						direction = get_random_direction();
					}
					else{
						direction = VECTOR2_SCALED(diff, 1/distance);
					}
					double slope = smoothing_kernel_derivative(distance, SMOOTHING_RADIUS);
					double density = particle_density[other_particle_index];
					double shared_pressure = calculate_shared_pressure(density, particle_density[particle_index]);
					double pressure_gradient_scaler = shared_pressure * slope * mass / density;
					pressure_gradient = VECTOR2_ADD(pressure_gradient, VECTOR2_SCALED(direction, pressure_gradient_scaler));
				}
			}
		}
	}

	return pressure_gradient;
}

static double smoothing_kernel(double distance, float radius){
	if(distance >= radius){
		return 0;
	}

	float volume = (NUM_PI * pow(radius, 4)) / 6;
	return (radius - distance) * (radius - distance) / volume;
}

static double smoothing_kernel_derivative(double distance, double radius){
	if(distance >= radius){
		return 0;
	}

	float scale = 12 / (pow(radius, 4) * NUM_PI);
	return (distance - radius) * scale;
}

static int get_cell_key_from_particle(int particle_index){
	// Map particle position to a cell
	tVector2_int coord = position_to_cell_coordinate(particle_index);
	uint32_t key = get_cell_key_from_grid_index(coord);
	return key;
}

static int get_cell_key_from_grid_index(tVector2_int grid){
	return combinging_hash(grid.x, grid.y) % NUM_PARTICLES;
}

static tVector2_int position_to_cell_coordinate(int particle_index){
	uint32_t x = (uint32_t)floor(position[particle_index].x / SMOOTHING_RADIUS);
	uint32_t y = (uint32_t)floor(position[particle_index].y / SMOOTHING_RADIUS);

	return (tVector2_int){
		.x = x,
		.y = y
	};
}

static int cell_sorting_compare(const void* c1, const void* c2){
	const cell_key_pair* p1 = (const cell_key_pair*)c1;
	const cell_key_pair* p2 = (const cell_key_pair*)c2;

	if(p1->key < p2->key){
		return -1;
	}
	else if(p1->key > p2->key){
		return 1;
	}
	return 0;
}

static void update_cell_keys(){
	for(int i = 0; i < NUM_PARTICLES; i++){
		int cell_key = get_cell_key_from_particle(i);
		spacial_LUT[i] = (cell_key_pair){
			.key = cell_key,
			.index = i
		};
	}

	qsort(spacial_LUT, NUM_PARTICLES, sizeof(cell_key_pair), cell_sorting_compare);

	for(int i = 0; i < NUM_PARTICLES; i++){
		uint32_t prev_key = i == 0 ? -1 : spacial_LUT[i - 1].key;
		if(spacial_LUT[i].key != prev_key){
			prev_key = spacial_LUT[i].key;
			spacial_LUT_start_indeces[spacial_LUT[i].key] = i;
		}
		else{
			spacial_LUT_start_indeces[i] = -1;
		}
	}
}
