#ifndef DRAWING_H_
#define DRAWING_H_

#ifdef __cplusplus
extern "C" {
#endif


/*********************
 *      INCLUDES
 *********************/

#include <SDL.h>
#include <stdint.h>

#include "common.h"

/*********************
 *      DEFINES
 *********************/

/**********************
 *      TYPEDEFS
 **********************/

/**********************
 * GLOBAL PROTOTYPES
 **********************/

void drawing_init();

void drawing_clearScreen();

void drawing_drawScreen();

void drawing_delay(uint32_t delay);

void drawing_drawCircle(tVector2_int pos, int32_t radius, tColor color);

void drawing_drawPixel(tVector2_int pos, tColor_a color);

/**********************
 *  GLOBAL VARIABLES
 **********************/

/**********************
 *      MACROS
 **********************/

/**********************
 *   POST INCLUDES
 *********************/

#ifdef __cplusplus
} /*extern "C"*/
#endif

#endif /* DRAWING_H_ */