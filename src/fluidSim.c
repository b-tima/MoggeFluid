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

/**********************
 *      TYPEDEFS
 **********************/

/**********************
 *  STATIC PROTOTYPES
 **********************/

static tVector2 translate_position(tVector2 vector);
static void draw_particle(tVector2 vector);

static void resolve_collision();

/**********************
 *  STATIC VARIABLES
 **********************/

static float gravity = 100;
static float collision_damping = 0.6;

static tVector2 position;
static tVector2 velocity;

const tVector2 bounding_box = {
	.x = SCREEN_HEIGHT,
	.y = SCREEN_HEIGHT
};

/**********************
 *      MACROS
 **********************/

/**********************
 *   GLOBAL FUNCTIONS
 **********************/

void fluidSim_init() {
    position = VECTOR2_ZERO;
    velocity = VECTOR2_ZERO;
}

void fluidSim_update() {
    velocity.y -= gravity * TIME_DELTA_S;
    position.y += velocity.y * TIME_DELTA_S;

	printf("y = %f\n", velocity.y);

	resolve_collision();

    draw_particle(position);
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
    tVector2 tranpos = translate_position(position);
    drawing_drawCircle(VECTOR2_TO_INT(tranpos), PARTICLE_RADIUS,
                       COLOR_LIGHTBLUE);
}

static void resolve_collision(){
	tVector2 half_bounds = VECTOR2_SCALED(bounding_box, 0.5);

	if(abs(position.x) + PARTICLE_RADIUS > half_bounds.x){
		position.x = (half_bounds.x - PARTICLE_RADIUS) * NUM_SIGN(position.x);
		velocity.x *= -1 * collision_damping;
	}

	if(abs(position.y) + PARTICLE_RADIUS > half_bounds.y) {
		position.y = (half_bounds.y - PARTICLE_RADIUS) * NUM_SIGN(position.y);
		velocity.y *= -1 * collision_damping;
	}
}
