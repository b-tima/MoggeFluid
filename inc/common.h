#ifndef COMMON_H_
#define COMMON_H_

/*********************
 *      INCLUDES
 *********************/

#include <stdint.h>

/*********************
 *      DEFINES
 *********************/

#define TIME_DELTA (20)
#define TIME_DELTA_S (((double)TIME_DELTA) / 1000)

#define SCREEN_WIDTH ((double)1270)
#define SCREEN_HEIGHT ((double)720)

/**********************
 *      TYPEDEFS
 **********************/

typedef struct {
    double x;
    double y;
} tVector2;

typedef struct {
    int32_t x;
    int32_t y;
} tVector2_int;

typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
} tColor;

typedef struct{
	uint8_t r;
	uint8_t g;
	uint8_t b;
	uint8_t a;
} tColor_a;

/**********************
 *      MACROS
 **********************/

#ifndef VECTOR2_ZERO
#define VECTOR2_ZERO ((tVector2){0, 0})
#endif

#ifndef VECTOR2_CREATE
#define VECTOR2_CREATE(x, y) ((tVector2){x,y})
#endif

#ifndef VECTOR2_ADD
#define VECTOR2_ADD(v1, v2) ((tVector2){(v1.x + v2.x), (v1.y + v2.y)})
#endif

#ifndef VECTOR2_SUB
#define VECTOR2_SUB(v1, v2) ((tVector2){((v1.x) - (v2.x)), ((v1.y) - (v2.y))})
#endif

#ifndef VECTOR2_SCALED
#define VECTOR2_SCALED(v, s) ((tVector2){(v.x * s), (v.y * s)})
#endif

#ifndef VECTOR2_MAG
#define VECTOR2_MAG(v) (sqrt(((v.x) * (v.x) + (v.y) * (v.y))))
#endif

#ifndef VECTOR2_TO_INT
#define VECTOR2_TO_INT(v) ((tVector2_int){(int32_t)v.x, (int32_t)v.y})
#endif

#ifndef BACKGROUND_COLOR
#define BACKGROUND_COLOR ((tColor){10, 10, 10})
#endif

#ifndef NUM_SIGN
#define NUM_SIGN(n) (n > 0 ? 1 : -1)
#endif

#ifndef NUM_ABS
#define NUM_ABS(n) (n > 0 ? n : -n)
#endif

#ifndef NUM_MAX
#define NUM_MAX(n1, n2) (n1 > n2 ? n1 : n2)
#endif

#ifndef NUM_MIN
#define NUM_MIN(n1, n2) (n1 < n2 ? n1 : n2)
#endif

#ifndef NUM_PI
#define NUM_PI (3.141592653589793)
#endif

#ifndef COLOR_WITH_A
#define COLOR_WITH_A(color, a) ((tColor_a){(color.r), (color.g), (color.b), (a)})
#endif

/**********************
 *      CONSTANTS
 **********************/

#ifndef VECTOR2_UP
#define VECTOR2_UP VECTOR2_CREATE(0, 1)
#endif

#ifndef VECTOR2_DOWN
#define VECTOR2_DOWN VECTOR2_CREATE(0, -1)
#endif

#ifndef VECTOR2_RIGHT
#define VECTOR2_RIGHT VECTOR2_CREATE(1, 0)
#endif

#ifndef VECTOR2_LEFT
#define VECTOR2_LEFT VECTOR2_CREATE(-1, 0)
#endif

#ifndef COLOR_RED
#define COLOR_RED ((tColor){255, 0, 0})
#endif

#ifndef COLOR_GREEN
#define COLOR_GREEN ((tColor){0, 255, 0})
#endif

#ifndef COLOR_BLUE
#define COLOR_BLUE ((tColor){0, 0, 255})
#endif

#ifndef COLOR_LIGHTBLUE
#define COLOR_LIGHTBLUE ((tColor){30,144,255})
#endif

#endif /* COMMON_H_ */