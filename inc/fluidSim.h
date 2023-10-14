#ifndef FLUIDSIM_H_
#define FLUIDSIM_H_

#ifdef __cplusplus
extern "C" {
#endif


/*********************
 *      INCLUDES
 *********************/

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

void fluidSim_init();

void fluidSim_update();

void fluidSim_onClick(tVector2_int pos);

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

#endif /* FLUIDSIM_H_ */