#include "sshead.h"
#include "vector.h"
#include "orb_params.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct rad_bin{
	double ecc;
	double h_over_r;
	int count;
} RAD_BIN;


#define N_BINS 1000.0
#define MAX_RAD 4.0
#define VERBOSE 0

void time_fix(SSHEAD *head, SSDATA *data);
void mass_fix(SSHEAD *head, SSDATA *data);
void delete_rand(SSHEAD *head, SSDATA *data, int nout);

int main(int argc, char *argv[]){
	int i;
	SSHEAD head;
	SSDATA *data;
	PARTICLE *particle;
	char *infile=NULL, *outfile=NULL;
	FILE *f;
	double *params;
	int n_bins;
	RAD_BIN *radbins;
	
	
	if(argc >= 2){
		infile = argv[1];
		outfile = argv[2];
	}
	if(infile==NULL) infile="stdin";
	if(outfile==NULL) {
		outfile="stdout";
		f = stdout;
	}
	else f = fopen(outfile, "w");
	
	printf("Fetching from file: %s\n", infile);
	printf("Sending to file: %s\n", outfile);
	
	
	if(infile=="stdin"){
		return(0);
	}
	
	data = readin_ss(infile, &head);
	//have ssdata read in by this point
	
	
	//MASS FIX
	//mass_fix(&head, data);
	
	//TIME_FIX
	//time_fix(&head, data);
	
	// NUMBER REDUCE
	//delete_rand(&head, data, 100000);
	
	//writeout_ss(infile, &head, data);
	
	
	
	
	particle = (PARTICLE*) malloc( (int)head.n_data*sizeof(PARTICLE) );
	if(particle==NULL){
		fprintf(stderr, "Error: Could not allocate memory in ftn 'main'\n");
		exit(0);
	}
	
	orbital_params(&head, data, particle);
	
	if(!VERBOSE){
		int i = 0;
		printf("Properties of particle %i:\n", i);
		printf("sma: %lf, ecc: %lf, inc: %lf\n",particle[i].sma, particle[i].ecc, particle[i].inc);
		printf("rx: %lf, ry: %lf, rz: %lf\n",particle[i].pos[0], particle[i].pos[1], particle[i].pos[2]);
		printf("vx: %lf, vy: %lf, vz: %lf\n", particle[i].vel[0], particle[i].vel[1], particle[i].vel[2]);
		printf("Mass: %E, Radius: %E\n", particle[i].mass, particle[i].radius);

		return(0);
	}

	printf("Binning data...\n", outfile);
	
	/*want to (sma) bin ecc and inc*/
	radbins = (RAD_BIN*) malloc( N_BINS*sizeof(RAD_BIN) );
	
	/*want to output in radial bins with <i> and <e>*/
	for(i=0; i<head.n_data; i++){
		int whichbin;
		
		whichbin = (int)floor(particle[i].sma*N_BINS/MAX_RAD);
		if( whichbin< N_BINS){
			/*if within correct sma range, bin the data*/
			/*do <e^2>^0.5, <i^2>^0.5 etc*/
			radbins[whichbin].ecc+=(particle[i].ecc*particle[i].ecc);
			radbins[whichbin].h_over_r+=(tan(particle[i].inc)*tan(particle[i].inc));
			radbins[whichbin].count++;
		}
	}
	
	
	for(i=0;i<N_BINS;i++){
		double whatrad = (double)i*MAX_RAD/N_BINS;
		double av_h_over_r = sqrt(radbins[i].h_over_r/radbins[i].count);
		double av_ecc = sqrt(radbins[i].ecc/radbins[i].count);
		
		if(radbins[i].count==0){
			av_h_over_r = 0.0;
			av_ecc = 0.0;
		}
		fprintf(f, "%E,%E,%E\n", whatrad, av_h_over_r, av_ecc);
	
	}

	close(f);
	free(radbins);
	free(particle);
	free(data);
	return(0);
}

void time_fix(SSHEAD *head, SSDATA *data){
	//fixes the time to line up with slower growth curve
	double gamma = 1.0/(1000*2*M_PI)*log(0.001/4.5E-5);
	head->time = log(data[0].mass/4.5E-5)*(1.0/gamma);
	return;
}

void mass_fix(SSHEAD *head, SSDATA *data){
	double new_mass, alpha, gamma, A, B, ta, tb;
	double TIME;
	
	//TIME = 8.0/10.0*6283.185307;
	TIME = head->time;
	A = 1.556937E-04;
	ta = 2513.274123;
	B = 1.001206E-03;
	tb = 6283.185307;
	gamma = (ta - tb)/log(A/B);
	alpha = A/exp(ta/gamma);
	new_mass = alpha*exp(TIME/gamma);
	//head->time = 0.0;
	data[0].mass = new_mass;
	return;
}

void swap_ss(SSDATA *data1, SSDATA *data2){
	SSDATA temp;
	temp = *data1; 
	*data1 = *data2;
	*data2 = temp;
	return;
}

void delete_rand(SSHEAD *head, SSDATA *data, int nout){
	int i, r;
	
	
	srand(time(NULL)); //seed random number generator
	for(i=1;i<nout;i++){
		r = rand()%(head->n_data-nout) + nout; //find random number
		swap_ss(&(data[i]),&(data[r])); //swap the two entries
	}
	head->n_data = nout; //will only write nout entries
	return;
}












