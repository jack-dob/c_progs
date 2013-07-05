#include <stdio.h>
#include "sshead.h"
#include "orb_params.h"
#include <math.h>


void bin_ecc(SSHEAD sshead, SSDATA *ssdata, int *ecc_bins, int nbins);
void print_bins_i(int nbins, double start, double stop, int *bin_array);

int main(int argc, char **argv){
	int i, j, n_eject;
	SSHEAD sshead;
	SSDATA *ssdata;
	char *infile;
	
	//assume that each argument is a file to read header of
	//printf("arg count: %d\n", argc);
	//exit(0);
	if(argc < 1){
		puts("Error: No input files passed, exiting...");
		exit(0);
	}

	
	printf("Years\tp_num\tejected\torb_num\n");
	for(i=1;i<argc;i++){
		infile = argv[i]; //1st argument is the input file
		ssdata = readin_ss(infile, &sshead); //read in header
		n_eject=0;
		for(j=0;j<sshead.n_data;j++){
			if(sma(ssdata[j]) < 0)
				n_eject++;
		}
		printf("%G\t%d\t%d\t%d\n", sshead.time/(2.0*M_PI), sshead.n_data, n_eject, sshead.n_data-n_eject);
		free(ssdata); //have to free to avoid memory leaks, maybe use realloc in readin_ss?
	}
	


}








