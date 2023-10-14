#include <stdio.h>
#include <stdlib.h>

#include "example_lib.h"

int main(){
	char* buffer = malloc(sizeof(char)*100);
	example_function(buffer);
	printf("%s", buffer);
	free(buffer);
	return 0;
}
