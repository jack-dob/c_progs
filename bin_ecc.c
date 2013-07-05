#include <stdio.h>
#include "sshead.h"
#include "orb_params.h"
#include <math.h>

void bin_ecc(SSHEAD sshead, SSDATA *ssdata, int *ecc_bins, int nbins);
void print_bins_i(int nbins, double start, double stop, int *bin_array);

int main(int argc, char **argv){
	int i, nbins;
	int *ecc_bins;
	SSHEAD sshead;
	SSDATA *ssdata;
	char *infile;
	
	if(argc < 1){
		puts("Error: No input file passed, exiting...");
		exit(0);
	}
	infile = argv[1]; //1st argument is the input file
	
	ssdata = readin_ss(infile, &sshead); //read in data
	//puts("A");
	nbins = 10000; //number of bins
	ecc_bins = (int*) malloc( nbins*sizeof(int));
	//puts("B");
	bin_ecc(sshead, ssdata, ecc_bins, nbins);
	//puts("C");
	print_bins_i(nbins, 0, 1, ecc_bins);


	free(ecc_bins);
	free(ssdata);
}

void bin_ecc(SSHEAD sshead, SSDATA *ssdata, int *ecc_bins, int nbins){
	//bins are always from 0 to 1 (how eccentricity works)
	int i, bin_num;
	double ecc1, sma1;
	double sma_cutoff = 1.5;
	
	for(i=0;i<nbins;i++){ //make sure this is zero, don't want surprises
		ecc_bins[i] = 0;
	}
	
	for(i=0;i<sshead.n_data;i++){
		sma1 = sma(ssdata[i]);
		//printf("i: %d, sma : %G, N: %d\n", i, sma1, sshead.n_data);
		if((sma1 > sma_cutoff) || (sma1 < 0)) //if outside our range or ejected
			continue; //next object
		ecc1 = ecc(ssdata[i]); //find eccentricity
		//puts("D");
		bin_num = (int)floor(ecc1*(double)nbins); //eccentricity is always a fraction from 0 to 1
		//puts("E");
		ecc_bins[bin_num] += 1; //add one to bin counter
	}
	return;
}


void print_bins_i(int nbins, double start, double stop, int *bin_array){
	int i;
	double x, dx;
	
	dx = (stop - start)/(double)nbins;
	printf("x\tN(x)\n");
	for(i=0;i<nbins;i++){
		x = start + dx*(double)i;
		printf("%G\t%d\n", x, bin_array[i]);
	}
	return;
}








