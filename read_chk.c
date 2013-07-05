#include <stdio.h>
#include <stdlib.h>
#include "fdl.h"
#include "fdl.c"
#include "htable.h"
#include "htable.c"

int main(int argc, char **argv){
	FDL_CTX *fdl;
	int an_int, another_int;
	double a_double, another_double;

#ifdef FDL_HINCLUDED
	printf("FDL included\n");
#endif	
#ifdef HTABLE_HINCLUDED
	printf("HTABLE included\n");
#endif

	/*This function modifies checkpoint files
	format:
	read_chk [FILENAME] [VARIABLE_NAME] [VALUE] [IDENT]
	
	filename = name of chk. file
	variable_name = name of variable to write 
	e.g. 'version', or 'dExtraStore' ->see master.c for full list
	ident = type of variable: n, i = integer, d = double*/
	fdl = FDL_modify(argv[1]);
	if(argv[4][0] == 'n' || argv[4][0] == 'i'){
		FDL_read(fdl,argv[2], &an_int);
		printf("Old value of %s is: %i\n", argv[2] ,an_int);
		another_int = atoi(argv[3]);
		printf("Writing value %i to checkpoint\n", another_int);
		FDL_write(fdl, argv[2], &another_int);
		FDL_read(fdl,argv[2], &an_int);
		printf("New value of %s is: %i\n", argv[2] ,an_int);
	}
	else if(argv[4][0] == 'd'){
		FDL_read(fdl, argv[2], &a_double);
		printf("Old value of %s is: %E\n", argv[2] ,a_double);
		another_double = atof(argv[3]);
		printf("Writing value %E to checkpoint\n", another_double);
		FDL_write(fdl, argv[2], &another_double);
		FDL_read(fdl, argv[2], &a_double);
		printf("New value of %s is: %E\n", argv[2] ,a_double);
	}
	
	return(0);
}