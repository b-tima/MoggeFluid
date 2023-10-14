/*********************
 *      INCLUDES
 *********************/

#include "example_lib.h"
#include <stdio.h>

/*********************
 *      DEFINES
 *********************/

/**********************
 *      TYPEDEFS
 **********************/

/**********************
 *  STATIC PROTOTYPES
 **********************/

/**********************
 *  STATIC VARIABLES
 **********************/

/**********************
 *      MACROS
 **********************/

/**********************
 *   GLOBAL FUNCTIONS
 **********************/

void example_function(char* hello_world){
	sprintf(hello_world, "Hello, %s\n", "World!");
}

int example_testable_function(int a){
	return a + 1;
}

/**********************
 *   STATIC FUNCTIONS
 **********************/
