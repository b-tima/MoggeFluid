#ifndef COMMON_H_
#define COMMON_H_

/*********************
 *      INCLUDES
 *********************/

#include <stdint.h>
#include <math.h>

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

inline tVector2 VECTOR2_CREATE(double x, double y){
	return (tVector2){
		.x = x,
		.y = y
	};
}

inline tVector2 VECTOR2_ADD(tVector2 v1, tVector2 v2){
	return (tVector2){
		.x = v1.x + v2.x,
		.y = v1.y + v2.y
	};
}

inline tVector2 VECTOR2_SUB(tVector2 v1, tVector2 v2){
	return (tVector2){
		.x = v1.x - v2.x,
		.y = v1.y - v2.y
	};
}

inline tVector2 VECTOR2_SCALED(tVector2 v, double scaler){
	return (tVector2){
		.x = v.x * scaler,
		.y = v.y * scaler
	};
}

inline double VECTOR2_SQR_MAG(tVector2 v){
	return (v.x*v.x + v.y*v.y);
}

inline double VECTOR2_MAG(tVector2 v){
	return sqrt(v.x*v.x + v.y*v.y);
}

inline tVector2_int VECTOR2_TO_INT(tVector2 v){
	return (tVector2_int){
		.x = (int32_t)v.x,
		.y = (int32_t)v.y
	};
}

inline double NUM_SIGN(double n){
	return n > 0 ? 1 : -1;
}

inline double NUM_ABS(double n){
	return n > 0 ? n : -n;
}

inline double NUM_MAX(double n1, double n2){
	return n1 > n2 ? n1 : n2;
}

inline double NUM_MIN(double n1, double n2){
	return n1 < n2 ? n1 : n2;
}

inline uint32_t combinging_hash(uint32_t n1, uint32_t n2){
	uint32_t a = n1 * 15823;
	uint32_t b = n2 * 9737333;
	uint32_t x = (a + b) * (a + b + 1) / 2 + a;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
	x = ((x >> 16) ^ x) * 0x45d9f3b;
	x = (x >> 16) ^ x;
	return x;
}

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

#ifndef VECTOR2_ZERO
#define VECTOR2_ZERO VECTOR2_CREATE(0, 0)
#endif

#ifndef BACKGROUND_COLOR
#define BACKGROUND_COLOR ((tColor){10, 10, 10})
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

#ifndef NUM_PI
#define NUM_PI (3.141592653589793)
#endif

#endif /* COMMON_H_ */