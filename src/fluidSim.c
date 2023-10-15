/*********************
 *      INCLUDES
 *********************/

#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include <omp.h>

#include "fluidSim.h"
#include "drawing.h"

/*********************
 *      DEFINES
 *********************/

#define PARTICLE_RADIUS (0.04)

#define NUM_PARTICLES (2000)

#define PARTICLE_SPACING (0.04)

#define SMOOTHING_RADIUS (0.3)
#define VISCOSITY_RADIUS (1)

#define MAX_VISCOSITY (10)

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

/**
 * 	If set to 1, the cells around the mouse will be drawn on mouse click
 * 	If set to 0, the cells will not be drawn
 */
#define DRAW_CELLS_ON_MOUSE_CLICK (0)

/**********************
 *      TYPEDEFS
 **********************/

typedef struct{
	int key;
	int index;
} cell_key_pair;

/**********************
 *  STATIC PROTOTYPES
 **********************/

static tVector2 translate_global_position(tVector2_int vector);
static tVector2_int translate_local_position(tVector2 vector);
static void draw_particle(tVector2 vector);

static void resolve_edge_collision(tVector2* local_position, tVector2* local_velocity);

static tVector2 get_random_direction();

static void update_densites();
static tVector2 interaction_force(tVector2 input_pos, double radius, double strength, int particle_index);
static double calculate_density(tVector2 sample_point);
static double convert_density_to_pressure(double density);
static double calculate_shared_pressure(double d1, double d2);
static void get_close_particles(int particle_index, int *particles, int* num_particles);
static tVector2 calculate_pressure_force(int particle_index);
static tVector2 calculate_viscosity_force(int particle_index);
static double smoothing_kernel(double radius, double distance);
static double smoothing_kernel_derivative(double radius, double distance);
static double viscosity_kernel(double radius, double distance);

static int get_cell_key_from_particle(int particle_index);
static int get_cell_key_from_grid_index(tVector2_int grid);
static tVector2_int particle_to_cell_coordinate(int particle_index);
static tVector2_int position_to_cell_coordinate(tVector2 sample_point);
static int cell_sorting_compare(const void* c1, const void* c2);
static void update_cell_keys();

/**********************
 *  STATIC VARIABLES
 **********************/

const double gravity = 0;
const double collision_damping = 0.1;
const double target_density = 10;
const double pressure_multiplier = 2.3;
const double mass = 0.1;

const double interaction_strenght = 0.08;
const double interaction_radius = 2;

const double viscosity_strenght = 1;

static tVector2 position[NUM_PARTICLES];
static tVector2 predicted_position[NUM_PARTICLES];
static tVector2 velocity[NUM_PARTICLES];
static tVector2 additional_force[NUM_PARTICLES];
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


static int rect_timeout = 0;
static tVector2_int rect_start[1000];
static tVector2_int rect_end[1000];
static int num_rects = 0;

/**********************
 *      MACROS
 **********************/

/**********************
 *   GLOBAL FUNCTIONS
 **********************/

void fluidSim_init() {
	for(int i = 0; i < NUM_PARTICLES; i++){
		particle_density[i] = 0;
		velocity[i] = VECTOR2_ZERO;
		additional_force[i] = VECTOR2_ZERO;
	}

	spacial_LUT_rows = ceil(bounding_box.y / SMOOTHING_RADIUS);
	spacial_LUT_cols = ceil(bounding_box.x / SMOOTHING_RADIUS);
	num_spacial_LUT = NUM_PARTICLES * 10;

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

void fluidSim_destroy(){
	free(spacial_LUT);
	free(spacial_LUT_start_indeces);
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

			double a = 255 * NUM_MIN(1, density/100);
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

	// Predict future position
	for(int i = 0; i < NUM_PARTICLES; i++){
		velocity[i] = VECTOR2_ADD(velocity[i], VECTOR2_SCALED(VECTOR2_DOWN, gravity * TIME_DELTA_S));
		predicted_position[i] = VECTOR2_ADD(position[i], VECTOR2_SCALED(velocity[i], 120*TIME_DELTA_S));
	}

	// Update densities and cell keys
	update_cell_keys();
	update_densites();

	double largest_pressure = 0;
	double largest_vis = 0;
	for(int i = 0; i < NUM_PARTICLES; i++){
		tVector2 pressure_force = calculate_pressure_force(i);
		tVector2 viscosity_force = calculate_viscosity_force(i);
		tVector2 particle_acceleration = VECTOR2_SCALED(VECTOR2_ADD(pressure_force, viscosity_force), 1/particle_density[i]);
		velocity[i] = VECTOR2_ADD(velocity[i], VECTOR2_SCALED(particle_acceleration, TIME_DELTA_S));
		velocity[i] = VECTOR2_ADD(velocity[i], additional_force[i]);

		double mag_p = VECTOR2_MAG(pressure_force);
		double mag_v = VECTOR2_MAG(viscosity_force);
		if(mag_p > largest_pressure){
			largest_pressure = mag_p;
		}

		if(mag_v > largest_vis){
			largest_vis = mag_v;
		}

		additional_force[i] = VECTOR2_ZERO;
	}

	for(int i = 0; i < NUM_PARTICLES; i++){
		position[i] = VECTOR2_ADD(position[i], velocity[i]);

		resolve_edge_collision(&position[i], &velocity[i]);

		draw_particle(position[i]);
	}

	if(rect_timeout > 0){
		for(int i = 0; i < num_rects; i++){
			drawing_drawHollowRect(rect_start[i], rect_end[i], COLOR_RED);
		}
		rect_timeout -= TIME_DELTA;
	}
	else{
		num_rects = 0;
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

	for(int i = 0; i < NUM_PARTICLES; i++){
		additional_force[i] = interaction_force(local_pos, interaction_radius, interaction_strenght, i);
	}

#if DRAW_CELLS_ON_MOUSE_CLICK
	// Get all of the cells in a 3x3 grid around the local position
	num_rects = 0;
	tVector2_int cell = position_to_cell_coordinate(local_pos);
	for(int x = -1; x < 2; x++){
		for(int y = -1; y < 2; y++){
			tVector2_int new_grid = {
				.x = cell.x + x,
				.y = cell.y + y
			};

			printf("(%d, %d)\n", new_grid.x, new_grid.y);
			
			int cell_key = get_cell_key_from_grid_index(new_grid);

			int num_particles = 0;
			for(int i = spacial_LUT_start_indeces[cell_key]; i < NUM_PARTICLES; i++){
				if(spacial_LUT[i].key != cell_key){
					break;
				}
				num_particles += 1;
			}

			// Find all cells with this cell key
			for(int i = 0; i < spacial_LUT_cols; i++){
				for(int j = 0; j < spacial_LUT_rows; j++){
					tVector2 cell_coord = {
						.x = i,
						.y = j
					};

					if(cell_key == get_cell_key_from_grid_index(VECTOR2_TO_INT(cell_coord))){
						// Draw a rectangle!
						cell_coord = VECTOR2_ADD(cell_coord, VECTOR2_SCALED(bounding_box, -0.5));
						//cell_coord = VECTOR2_ADD(cell_coord, VECTOR2_CREATE(SMOOTHING_RADIUS / 2, SMOOTHING_RADIUS / 2));

						tVector2 cell_end = VECTOR2_ADD(cell_coord, VECTOR2_CREATE(SMOOTHING_RADIUS, SMOOTHING_RADIUS));

						tVector2_int global_start = translate_local_position(cell_coord);
						tVector2_int global_end = translate_local_position(cell_end);
						memcpy(&rect_start[num_rects], &global_start, sizeof(tVector2_int));
						memcpy(&rect_end[num_rects], &global_end, sizeof(tVector2_int));
						num_rects += 1;
					}
				}
			}
		}
	}

	rect_timeout = 2000;
#endif
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

static tVector2_int translate_local_position(tVector2 vector) {
    tVector2_int translated_vector = {.x = vector.x * SCREEN_SCALING + (SCREEN_WIDTH / 2),
                                  .y = (SCREEN_HEIGHT / 2) - vector.y * SCREEN_SCALING};
    return translated_vector;
}

static void draw_particle(tVector2 vector) {
    tVector2_int tranpos = translate_local_position(vector);
    drawing_drawCircle(tranpos, PARTICLE_RADIUS * SCREEN_SCALING,
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

static tVector2 interaction_force(tVector2 input_pos, double radius, double strength, int particle_index){
	tVector2 force = VECTOR2_ZERO;
	tVector2 offset = VECTOR2_SUB(input_pos, position[particle_index]);
	double sqr_distance = VECTOR2_SQR_MAG(offset);

	if(sqr_distance < radius * radius){
		double distance = sqrt(sqr_distance);

		tVector2 dir_to_input_point = VECTOR2_ZERO;
		if(distance > 0.00001){
			dir_to_input_point= VECTOR2_SCALED(offset, 1/distance);
		}

		double centre_t = 1 - distance / radius;

		force = VECTOR2_ADD(force, 
			VECTOR2_SUB(
				VECTOR2_SCALED(dir_to_input_point, strength),
				velocity[particle_index]
			)
		);

		force = VECTOR2_SCALED(force, centre_t);
	}

	return force;
}

static tVector2 get_random_direction(){
	double x = (double)((double)rand()) / RAND_MAX;
	double y = (double)((double)rand()) / RAND_MAX;
	tVector2 direction = VECTOR2_CREATE(x, y);
	return VECTOR2_SCALED(direction, VECTOR2_MAG(direction));
}

static void update_densites(){
	#pragma omp parallel for num_threads(20)
	for(int i = 0; i < NUM_PARTICLES; i++){
		particle_density[i] = calculate_density(predicted_position[i]);
	}
}

static double calculate_density(tVector2 sample_point){
	double density = 0;

	#pragma omp parallel for reduction(+:density) num_threads(20)
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

static void get_close_particles(int particle_index, int *particles, int* num_particles){
	tVector2_int particle_cell_coord = particle_to_cell_coordinate(particle_index);
	double sqr_radius = SMOOTHING_RADIUS * SMOOTHING_RADIUS;

	*num_particles = 0;

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
				tVector2 diff = VECTOR2_SUB(predicted_position[other_particle_index], predicted_position[particle_index]);
				double sqr_distance = VECTOR2_SQR_MAG(diff);

				if(sqr_distance < sqr_radius){
					particles[*num_particles] = other_particle_index;
					*num_particles += 1;
				}
			}
		}
	}
}

static tVector2 calculate_pressure_force(int particle_index){
	int close_particles[NUM_PARTICLES];
	int num_particles;
	get_close_particles(particle_index, close_particles, &num_particles);

	tVector2 pressure_gradient = VECTOR2_ZERO;

	for(int i = 0; i < num_particles; i++){
		int other_particle_index = close_particles[i];

		tVector2 diff = VECTOR2_SUB(predicted_position[other_particle_index], predicted_position[particle_index]);
		double distance = VECTOR2_MAG(diff);
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

	return pressure_gradient;
}

static double smoothing_kernel(double distance, double radius){
	if(distance >= radius){
		return 0;
	}

	double radius_4 = radius * radius * radius * radius;

	double volume = (NUM_PI * radius_4) / 6;
	return NUM_MAX(0, (radius - distance) * (radius - distance) / volume);
}

static double smoothing_kernel_derivative(double distance, double radius){
	if(distance >= radius){
		return 0;
	}

	double radius_4 = radius * radius * radius * radius;

	double scale = 12 / (radius_4 * NUM_PI);
	return (distance - radius) * scale;
}

static tVector2 calculate_viscosity_force(int particle_index){
	int close_particles[NUM_PARTICLES];
	int num_particles;
	get_close_particles(particle_index, close_particles, &num_particles);

	tVector2 force = VECTOR2_ZERO;

	for(int i = 0; i < num_particles; i++){
		int other_particle_index = close_particles[i];

		tVector2 diff = VECTOR2_SUB(predicted_position[other_particle_index], predicted_position[particle_index]);
		double distance = VECTOR2_MAG(diff);

		double influence = viscosity_kernel(distance, VISCOSITY_RADIUS);

		tVector2 f = VECTOR2_ADD(force, 
			VECTOR2_SUB(
				velocity[other_particle_index],
				velocity[particle_index]
			)
		);
		force = VECTOR2_SCALED(f, influence);
	}

	if(VECTOR2_MAG(force) > MAX_VISCOSITY){
		force = VECTOR2_SCALED(force, MAX_VISCOSITY/VECTOR2_MAG(force));
	}

	return VECTOR2_SCALED(force, viscosity_strenght);
}

static double viscosity_kernel(double distance, double radius){
	if(distance > radius){
		return 0;
	}

	double radius_8 = radius * radius * radius * radius * radius * radius * radius * radius;

	double value = radius * radius - distance * distance;
	double volume = (NUM_PI * radius_8) / 4;
	return value * value * value / volume;
}

static int get_cell_key_from_particle(int particle_index){
	// Map particle position to a cell
	tVector2_int coord = particle_to_cell_coordinate(particle_index);
	int key = get_cell_key_from_grid_index(coord);
	return key;
}

static int get_cell_key_from_grid_index(tVector2_int grid){
	return combinging_hash(grid.x, grid.y) % num_spacial_LUT;
}

static tVector2_int particle_to_cell_coordinate(int particle_index){
	return position_to_cell_coordinate(predicted_position[particle_index]);
}

static tVector2_int position_to_cell_coordinate(tVector2 sample_point){
	uint32_t x = (uint32_t)floor(sample_point.x / SMOOTHING_RADIUS) + bounding_box.x / 2;
	uint32_t y = (uint32_t)floor(sample_point.y / SMOOTHING_RADIUS) + bounding_box.y / 2;

	return (tVector2_int){
		.x = x,
		.y = y
	};
}

static int cell_sorting_compare(const void* c1, const void* c2){
	const cell_key_pair* p1 = (const cell_key_pair*)c1;
	const cell_key_pair* p2 = (const cell_key_pair*)c2;

	return p1->key - p2->key;
}

static void update_cell_keys(){
	for(int i = 0; i < NUM_PARTICLES; i++){
		int cell_key = get_cell_key_from_particle(i);
		spacial_LUT[i] = (cell_key_pair){
			.key = cell_key,
			.index = i
		};
		spacial_LUT_start_indeces[i] = -1;
	}

	qsort((void*)spacial_LUT, NUM_PARTICLES, sizeof(cell_key_pair), cell_sorting_compare);

	for(int i = 0; i < NUM_PARTICLES; i++){
		int key = spacial_LUT[i].key;
		int prev_key = i == 0 ? -1 : spacial_LUT[i - 1].key;
		if(key != prev_key){
			spacial_LUT_start_indeces[spacial_LUT[i].key] = i;
		}
	}
}
