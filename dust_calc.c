#ifndef DEBUG
	#define DEBUG 0
#endif
#ifndef STDIO_INCLUDED
	#include <stdio.h>
	#define STDIO_INCLUDED
#endif
#ifndef SSHEAD_INCLUDED
	#include "sshead.h"
	#define SSHEAD_INCLUDED
#endif
#ifndef VECTOR_INCLUDED
	#include "vector.h"
	#define VECTOR_INCLUDED
#endif
#ifndef MATH_INCLUDED
	#include <math.h>
	#define MATH_INCLUDED
#endif
#ifndef FITSIO_INCLUDED
	#include "fitsio.h" /*Included in 'cfitsio' package*/
#endif


#ifndef MALLOC_ERROR
	#define MALLOC_ERROR "Error: %s could not allocate memory to array %s\n"
#endif

#define PI 3.14159265359 /* approximate value of pi */
#define SQRT_PI 1.77245385091 /* approximate value of sqrt(pi) */
#define SOLAR_MASS 1.99E30 /*kg per solar mass*/
#define AU 1.5E11 /*meters per astronomical unit*/

/*Yay, pre-processor*/
typedef struct dust_blob{
	double *pos;
	double mass;
	int n_particles;
	int org_idx;
}DUST_BLOB;

struct parameters { /*Holds all the parameters for the routines, can pass around as needed*/
	double image_scale; /*defines the scale of the image*/
	double cube_res; /*resolution in x, y direction*/
	double cube_scale; /*the size of the image with respect to the data*/
	double cube_side; /*the side length of a cube-cell*/
	double cube_z_scale; /*the size of the image in z-direction*/
	double cube_z_res; /*the resolution of the z-direction, fixed so that grid is square*/
	double cube_min; /*minimum value in x-y direction*/
	double cube_max; /*maximum value in x-y direction*/
	double cube_z_min; /* minimum value in z-direction*/
	double cube_z_max; /*maximum value in z-direction*/
	double orb_res;/*number of 'blobs' to split planetesimals into*/
	double smallest_mass; /*the smallest mass a planetesimal can be*/
	double smallest_dust; /*want the mass of a 10nm dust grain in kg: 4/3*pi*r^3*density = mass*/
	double destruction_radius; /*dust destruction radius in AU*/
	char *fitsname; /*default name for fits file*/
	char *gridname; /*default name for grid output*/
	char *dataname; /*default name for data output*/
	double dust_rad_min; //minimum radius of dust
	double dust_rad_max; //maximum radius of dust
	double debris_res; //maximum size of stuff that isn't a particle
	double debris_blowout; //blowout size due to radiation pressure
	double dust_factor; /*factor either side of centre*/
	double dust_min; /*minimum mass in dust bin*/
	double dust_max; /*maximum mass in dust bin*/
	double dust_norm; //mass to normalise the amount of dust to, all dust will add up to this mass
	//double ****scalars; //array of cubes to be point-wise multiplied with the dust cubes
	//int n_scalars; //how many of those scaling cubes are there?
	int n_rad_bins; //number of radial bins to use when azimuthally smoothing
	double rad_bins_inner; //where to start the bins (in AU, will be continuous)
	double rad_bins_outer; //where to end the bins (int AU)
	char *rad_bins_file; //name default name of resulting fits file
	char *dustopac; //path of dustopac.inp file
};

#ifndef DUST_PARAMS_INCLUDED
	#include "dust_params.h"
	#define DUST_PARAMS_INCLUDED
#endif
#include "radmc3d_conv.h" //routines to convert stuff to radmc3d formats

/*function declarations*/
double blob_mass(struct parameters *PAR, DUST_ELLIPSES dust_ellipse, double constant);
double surf_area(struct parameters *PAR, double mass, double constant);
double number_of_particles(double mass, double constant);
double sum_of_particles(struct parameters *PAR, double mass, double constant);
void create_dust_cube(struct parameters *PAR, DUST_ELLIPSES *dust_ellipses, SSHEAD head, double ****dust_cube);
void add_to_cube(struct parameters *PAR, DUST_BLOB *dust_blobs, double ****dust_cube);
void filter_dust_cube(struct parameters *PAR, double ****dust_cube, SSDATA *data, SSHEAD *head);
void flatten_cube(struct parameters *PAR, int nx, int ny, int nz, double ****dust_cube);
void az_smooth_square(struct parameters *PAR, int nx, int ny, double ****dust_cube);
void write_dust_cube(struct parameters *PAR, char* outfile, int nx, int ny, int nz, double ***dust_cube, char* infile);
void txt_dust_cube(struct parameters *PAR, double **dust_cube);

void init_par(struct parameters *PAR); /*Holds initialisation for PAR structure, see line 815*/

/*Start Code*/
int main(int argc, char *argv[]){
	int i, j, k;
	SSHEAD head;
	SSDATA *data;
	DUST_ELLIPSES *dust_ellipses;
	char *infile=NULL;
	FILE *f, *dustopac;
	double ***dust_cube;
	struct parameters *PAR;
	
	PAR = (struct parameters*) malloc( (int)sizeof(struct parameters) );
	
	init_par(PAR);
	
	printf("PAR->cube_z_res: %G, make sure this is an integer\n", PAR->cube_z_res);
	if (PAR->cube_z_res < 1.0){
		printf("ERROR: Need to increase resolution, or do in 2D\n");
		exit(1);
	}
	
	/*READ ARGUMENTS*/
	if(argc >= 3){
		infile = argv[1];
		PAR->fitsname = argv[2];
	}
	else{
		printf("Please pass an input and output file...\n");
		exit(1);
	}
	if(argc >= 4){
		PAR->rad_bins_file = argv[3];
	}
	
	//if(infile==NULL) {
	//	printf("Please pass an input file as 1st argument\n");
	//	exit(1);
	//}
	//if(PAR->fitsname==NULL) {
	//	printf("Please pass an output file as 2nd argument\n");
	//	exit(1);
	//}
	printf("Fetching from file: %s\n", infile);
	printf("Sending to file: %s\n", PAR->fitsname);
	printf("Smoothing to file: %s\n", PAR->rad_bins_file);
	/*make sure we're using decent files*/
	if(infile=="stdin"){
		printf("Error: Please enter valid input file\n");
		exit(1);
	}
	if(PAR->fitsname=="stdout"){
		printf("Error: Please enter a valid output file as second argument\n");
		exit(1);
	}
	/* END READ FILES */
	
	/*read in data*/
	data = readin_ss(infile, &head);
	
	/*create container for data*/
	dust_ellipses = (DUST_ELLIPSES*) malloc( (int)head.n_data*sizeof(DUST_ELLIPSES) );
	if(dust_ellipses==NULL){
		fprintf(stderr, "Error: Could not allocate memory in to array 'dust_ellipses' in ftn 'main'\n");
		exit(1);
	}
	if(DEBUG) printf("DEBUG: created dust ellipses array...\n");
	/*get orbital parameters of all particles in array 'particle'*/
	dust_params(&head, data, dust_ellipses);
	

	if(0) {
		i = 535;
		printf("DEBUG: outputting data for %ith particle, eject_flag: %i...\n", i, dust_ellipses[i].ejected_flag);
		printf("org_idx: %i, max_mass: %E, ecc: %E, inc: %E\n", dust_ellipses[i].org_idx, dust_ellipses[i].max_mass, dust_ellipses[i].ecc, dust_ellipses[i].inc);
		printf("sma: %E, psi: %E, dr: %E\n", dust_ellipses[i].sma, dust_ellipses[i].psi, dust_ellipses[i].dr);
		printf("pos: %G, %G, %G\n", data[i].pos[0], data[i].pos[1], data[i].pos[2]);
		printf("vel: %G, %G, %G\n", data[i].vel[0], data[i].vel[1], data[i].vel[2]);
		printf("Inverse Transform:\n");
		for(j=0;j<3;j++){
			printf("%E\t%E\t%E\n", dust_ellipses[i].inv_transform[j][0],dust_ellipses[i].inv_transform[j][1],dust_ellipses[i].inv_transform[j][2]);
		}
	}		
	/*close files and free arrays*/
	close(f);
	//free(data);
	//data = NULL;
	/*START OUTPUT ROUTINES*/
	/* Put stuff here that will create pictures from data*/
	/*Best idea is to split up x-y plane orbit into blobs
	and inverse transform them one at a time, onto a 3d grid and replace with appropriate amount of dust
	May need to sample orbit at high resolution to avoid gaps after transform*/
	
	/*Allocate memory for dust_cube*/
	printf("Allocating dust_cube array...\n");

	dust_cube = (double***) malloc( (int)(PAR->cube_res*sizeof(double**)) );
	if(dust_cube==NULL){
		printf(MALLOC_ERROR, __FUNCTION__, "dust_cube");
		exit(1);
	}
	for(i=0;i<PAR->cube_res;i++){
		dust_cube[i] = (double**) malloc((int)PAR->cube_res*sizeof(double*));
		if(dust_cube[i]==NULL){
			char *str;
			sprintf(str, "dust_cube[%i]", i);
			printf(MALLOC_ERROR, __FUNCTION__, str);
			exit(1);
		}
		for(j=0;j<PAR->cube_res;j++){
			dust_cube[i][j] = (double*) calloc(PAR->cube_z_res, sizeof(double));
			if(dust_cube[i][j] == NULL){
				char *str;
				sprintf(str, "dust_cube[%i][%i]", i, j);
				printf(MALLOC_ERROR, __FUNCTION__, str);
				exit(1);
			}
		}
	}

	
	printf("dust_cube array allocated.\n");
	/*finished memory allocation*/
	

	/*printf("Inverse Transform again: %i\n", dust_ellipses[0].inv_transform);
	for(j=0;j<3;j++){
		printf("%E\t%E\t%E\n", dust_ellipses[0].inv_transform[j][0],dust_ellipses[0].inv_transform[j][1],dust_ellipses[0].inv_transform[j][2]);
	}*/
	printf("create the dust cube..\n");
	create_dust_cube(PAR, dust_ellipses, head, &dust_cube); /*corrected for 3d case*/
	
	printf("Applying filters to 'dust_cube'...\n");
	filter_dust_cube(PAR, &dust_cube, data, &head);
	
	printf("write the dust cube...\n");
	/*dust_cube[500][500] = 10000.0;*/
	write_dust_cube(PAR, PAR->fitsname, PAR->cube_res, PAR->cube_res, PAR->cube_z_res, dust_cube, infile);
	
	puts("Checking for dust opacity input file...");
	dustopac = fopen(PAR->dustopac, "r");
	if(dustopac!=NULL){
		puts("File found, using size and proportion data...");
		radmc3d_grid(PAR); //write the dust grid
		radmc3d_dust(PAR, dust_cube, dustopac); //write the dust file
	}
	else{
		puts("Could not find file, using default parameters");
		printf("writing radmc3d input files...\n");
		if(PAR->gridname == "amr_grid.binp"){
			radmc3d_grid_b(PAR);
		}
		else{
			radmc3d_grid(PAR);
		}
		if(PAR->dataname == "dust_density.binp"){
			radmc3d_data_b(PAR, dust_cube);
		}
		else{
			radmc3d_data(PAR, dust_cube);
		}
	}
	//txt_dust_cube(dust_cube, PAR->fitsname);
	puts("flattening cube...");
	flatten_cube(PAR, PAR->cube_res, PAR->cube_res, PAR->cube_z_res, &dust_cube);
	puts("smoothing square...");
	az_smooth_square(PAR, PAR->cube_res, PAR->cube_res, &dust_cube);
	puts("writing dust square...");
	write_dust_cube(PAR, PAR->rad_bins_file, PAR->cube_res, PAR->cube_res, 1, dust_cube, infile);	
	
	
	/*free arrays*/
	
	printf("freeing arrays...\n");
	
	for(i=0;i<PAR->cube_res;i++){
		for(j=0;j<PAR->cube_res;j++){
			if(dust_cube[i][j]==NULL) continue;
			free(dust_cube[i][j]);
		}
	}
	for(i=0;i<PAR->cube_res;i++){
		if(dust_cube[i] == NULL) continue;
		free(dust_cube[i]);
	}
	free(data);
	data = NULL;
	free(dust_cube);
	free(dust_ellipses);
	printf("arrays freed, ending program...\n");
	return(0);
}


void create_dust_cube(struct parameters *PAR, DUST_ELLIPSES *dust_ellipses, SSHEAD head, double ****dust_cube){
	int i, j, k;
	DUST_BLOB *dust_blobs;
	double *temp_pos, *temp_pos2;
	
	/*Allocate Memory*/
	printf("Allocating memory to temporary storage array\n");
	temp_pos = (double*) calloc(3, sizeof(double));
	if(temp_pos==NULL){
		printf(MALLOC_ERROR, __FUNCTION__, "temp_pos");
	}
	printf("temp_pos allocated.\nAllocating dust_blobs array...\n");
	dust_blobs = (DUST_BLOB*) malloc( (int)(PAR->orb_res*sizeof(DUST_BLOB)) );
	if(dust_blobs == NULL){
		printf(MALLOC_ERROR, "dust_blobs");
		exit(1);
	}
	for(j=0;j<PAR->orb_res;j++){
		dust_blobs[j].pos = (double*) malloc((int)(3*sizeof(double*)));
		if(dust_blobs[j].pos==NULL){
			char *str;
			sprintf(str, "dust_blobs[%i].pos", j);
			printf(MALLOC_ERROR, __FUNCTION__, str);
			exit(1);
		}
	}
	
	/*Main processing loop*/
	printf("dust_blob array allocated.\nStarting loop...\n");
	for(i=1;i<head.n_data;i++){ //start at 1, avoid Jupiter mass planet
		/*find semi-latus rectum*/
		double P = dust_ellipses[i].r_a*dust_ellipses[i].r_p/dust_ellipses[i].sma;
		double smia;
		double phi;
		double constant;
		/*START PARTICLE FILTERING
		ignore particle if it has been ejected
		Put any other particle filtering stuff here also*/
		if(dust_ellipses[i].ejected_flag){
			continue;
		}
		/*if(dust_ellipses[i].max_mass > (100.0*PAR->smallest_mass)){
			continue;
		}/**/
		/*if(dust_ellipses[i].ecc < 0.03){
			continue;
		}*/
		
		/*if(dust_ellipses[i].sma < 0.5){
			continue;
		}*/
		/*END PARTICLE FILTERING*/
		
		printf("Calculating position of dust from particle %i\r", i);
		/*find semi-minor axis*/
		smia = sqrt(dust_ellipses[i].sma*P);
		/*	1 find position of mini-particles in xy plane using step in theta
			2 store locations in dust_blobs array
			3 transform positions and place in dust_cube array
			--See 'Making images' 11/06/2012 in lab-book--
		*/
		
		phi=0;
		/*calculate constant that goes into 'blob_mass' funtion*/
		constant = pow((dust_ellipses[i].max_mass*SOLAR_MASS)/PAR->orb_res, 11.0/6.0); /*this is for dust of ALL SIZES*/
		//printf("in 'created_dust_cube': constant: %G\n", constant);
		/*constant*=pow(4.0*PI, 1.0/3.0);
		constant*=pow(3.0, 2.0/3.0);
		constant/=4.0;/*at the moment this is for the surface area function*/
		/*start creating mini-particles from the apogee*/
		for(j=0;j<PAR->orb_res;j++){
			double r = P/(1.0 - dust_ellipses[i].ecc*cos(phi));
			double theta = phi + dust_ellipses[i].psi;
			double d_phi = 2*PI*dust_ellipses[i].sma*smia/(r*r*PAR->orb_res);
			double x = r*sin(theta);
			double y = r*cos(theta);
			double z = 0.0;

			temp_pos[0] = x;
			temp_pos[1] = y;
			temp_pos[2] = z;

			col_vec_multi(dust_ellipses[i].inv_transform, temp_pos, &dust_blobs[j].pos, 3, 3);

			dust_blobs[j].mass = blob_mass(PAR, dust_ellipses[i], constant);
			
			phi+=d_phi;
		}
		if(0){
			/*FOR DEBUGGING*/
			printf("\nBlobs at positions:\n");
			for(j=0;j<PAR->orb_res;j++){
				/*
				printf("Blob %i: x: %E, y: %E, z: %E\n", j, 
					dust_blobs[j].pos[0], dust_blobs[j].pos[1], dust_blobs[j].pos[2]);
				*/
				printf("%E,%E,%E\n", dust_blobs[j].pos[0], dust_blobs[j].pos[1], dust_blobs[j].pos[2]);
			}
			/*END FOR DEBUGGING*/
		}
		//puts("nowhere");
		add_to_cube(PAR, dust_blobs, dust_cube);
		
		//puts("bloop");
	}
	printf("\n");
	return;
}


void add_to_cube(struct parameters *PAR, DUST_BLOB *dust_blobs, double ****dust_cube){
	int i;
	
	for(i=0;i<PAR->orb_res;i++){
		int xbin, ybin, zbin;
		xbin = (int)floor((dust_blobs[i].pos[0] - PAR->cube_min)*PAR->cube_res/PAR->cube_scale);
		ybin = (int)floor((dust_blobs[i].pos[1] - PAR->cube_min)*PAR->cube_res/PAR->cube_scale);
		zbin = (int)floor((dust_blobs[i].pos[2] - PAR->cube_z_min)*PAR->cube_z_res/PAR->cube_z_scale);
		
		if( (xbin >= PAR->cube_res) || (ybin >= PAR->cube_res) || (zbin >= PAR->cube_z_res) || (xbin<0) || (ybin<0) || (zbin<0) ){
			/*printf("Warning: could not write dust blob, %i %i %i is outside range %i %i...\r",
				xbin, ybin, zbin, 0, PAR->cube_res);
			*/
			continue;		
		}
		//printf("\ndust_blob[%i] bins: %i, %i, %i, cube_res: %g, cube_z_res: %g \n",i, xbin, ybin, zbin, PAR->cube_res, PAR->cube_z_res);
		/*cfits uses fortrans style ordering, will have to re-format this data later*/

		(*dust_cube)[xbin][ybin][zbin] += dust_blobs[i].mass;
	}
	
	return;
}


double blob_mass(struct parameters *PAR, DUST_ELLIPSES dust_ellipse, double constant){
	/*
	put algorithm here!
	This finds the mass of the blobs of dust that are scattered around the particle's 
	initial orbit
	*/
	/*this is the maximum mass of a power law distribution*/
	//double max_mass = (dust_ellipse.max_mass*SOLAR_MASS)/PAR->orb_res;
	/*max_mass is now in kilograms*/
	//double rad = 0.00001; /*radius of particles we are interested in (in meters)*/
	//double density = 2000.0; /* assumed density of rock (in kg m^-3)*/
	//double mass_d = 4.0/3.0*PI*rad*rad*rad*density;
	//double dust_mass, m_big = 1E-2, m_small = 1E-4;
	//dust_mass = constant*(pow(PAR->dust_max, 1.0/6.0)-pow(PAR->dust_min,1.0/6.0))*dust_ellipse.t_coll;
	//printf("in blob_mass: constant: %G, dust_mass: %G, density: %G\n", constant, dust_mass, density);
	//return(dust_mass);
	return(dust_ellipse.max_mass*dust_ellipse.t_coll/PAR->orb_res);
	//return(1.0);
	/*return(mass_d*number_of_particles(mass_d, constant));*/
	/*return(mass_d*(sum_of_particles(0.99*mass_d, constant) - sum_of_particles(1.01*mass_d, constant)));/**/
	/*gives the mass of particles with a mass 1% around mass_d - the mass of a dust particle with radius = rad*/
	
	/*return( pow(mass_d, (1.0/6.0)) );*/ 
	/*return(max_mass);/*debug*/
	/*return(surf_area(max_mass, constant)); /*leave this for now, fix later*/
}

double surf_area(struct parameters *PAR, double mass, double constant){
	/*this should give surface area  of particles with mass
	greater than 'mass'*/
	double small_SA = (1.0/pow((PAR->smallest_dust), 1.0/6.0))*constant;
	double large_SA = (1.0/pow(mass, 1.0/6.0))*constant;
	double diff_SA = small_SA - large_SA;

	return(diff_SA);	/*maybe use surf area of circle as can only see 1/2 of sphere
						 at any one time*/

}

double number_of_particles(double mass, double constant){
	/*this gives the number of particles with a given mass
	  according to the Donyhani power law distribution*/
	double num = pow(mass,(-11.0/6.0))*constant;
	return(num);
}

double sum_of_particles(struct parameters *PAR, double mass, double constant){
	/*this gives the number of particles with mass greater than the 
	  given mass*/
	double sum_small = pow((PAR->smallest_dust), -5.0/6.0)*constant;
	double sum = pow(mass, -5.0/6.0)*constant;
	double sum_diff = sum_small - sum;
	return(sum_diff);
}


void filter_dust_cube(struct parameters *PAR, double ****dust_cube, SSDATA *data, SSHEAD *head){
	/*Use this function to apply various filters that are necessary*/
	int i, j, k, l;
	double max_SA = (PAR->cube_scale*AU/PAR->cube_res)*(PAR->cube_scale*AU/PAR->cube_res);
	double sum=0.0, sum2=0.0;
	double constant;
	/*more filtering*/
	
	//vel_disp_map(PAR, head, data, dust_cube);
	
	for(i=0;i<PAR->cube_res;i++){
		double x2 = (PAR->cube_scale/PAR->cube_res)*(PAR->cube_res/2.0 - (double)i)*(PAR->cube_scale/PAR->cube_res)*(PAR->cube_res/2.0 - (double)i);
		for(j=0;j<PAR->cube_res;j++){
			double y2 = (PAR->cube_scale/PAR->cube_res)*(PAR->cube_res/2.0 - (double)j)*(PAR->cube_scale/PAR->cube_res)*(PAR->cube_res/2.0 - (double)j);
			for(k=0;k<PAR->cube_z_res;k++){
				double z2 = (PAR->cube_z_scale/PAR->cube_z_res)*(PAR->cube_z_res/2.0 - (double)k)*(PAR->cube_z_scale/PAR->cube_z_res)*(PAR->cube_z_res/2.0 - (double)k);
				double r2 = x2+y2+z2;
				if(r2 < (PAR->destruction_radius*PAR->destruction_radius)){ /* no dust inside dust destruction radius*/
					(*dust_cube)[i][j][k]=0.0;
				}
				else{
					(*dust_cube)[i][j][k] /= sqrt(sqrt(r2)); //scale by light spreadyoutyness (for testing)
				}
				sum+=(*dust_cube)[i][j][k];
			}
		}
	}

		
	if(sum != PAR->dust_norm){
		constant = PAR->dust_norm/sum;
		for(i=0;i<PAR->cube_res;i++){
			for(j=0;j<PAR->cube_res;j++){
				for(k=0;k<PAR->cube_z_res;k++){
					(*dust_cube)[i][j][k]*=constant;
					sum2+=(*dust_cube)[i][j][k];
				}
			}
		}
	}
	printf("dust mass before: %G\nnomalised dust mass: %G\n", sum, sum2);
	return;
}

//dust cube operations

void flatten_cube(struct parameters *PAR, int nx, int ny, int nz, double ****dust_cube){
	//always flattens to the [0] layer in the z direction
	int i, j, k;
	double sum;
	
	for(i=0;i<nx;i++){
		for(j=0;i<ny;i++){
			sum = 0;
			for(k=0;k<nz;k++){
				sum += (*dust_cube)[i][j][k];
			}
			(*dust_cube)[i][j][0] = sum;
		}
	}
	return;
}

void az_smooth_square(struct parameters *PAR, int nx, int ny, double ****dust_cube){
	//will assume that stuff is on z[0] dimension of cube
	int i, j, r;
	double r2;
	double *rad_bins;
	double dr = (PAR->rad_bins_outer - PAR->rad_bins_inner)/((double)PAR->n_rad_bins);
	rad_bins = (double*) calloc( PAR->n_rad_bins, sizeof(double) );
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			r2 = sqrt(((double)i-(double)nx/2.0)*((double)i-(double)nx/2.0) + ((double)j-(double)ny/2.0)*((double)j-(double)ny/2.0)); //finds radius
			r2 *= PAR->cube_scale/nx;//now in AU
			r = (int)floor((r2 - PAR->rad_bins_inner)/dr); //which radial bin is it in?
			if( (r >= 0) && (r < PAR->n_rad_bins) ){ //make sure it's within bin ranges
				rad_bins[r] += (*dust_cube)[i][j][0]; //add to radial bins
			}
		}
	}
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			r2 = sqrt(((double)i-(double)nx/2.0)*((double)i-(double)nx/2.0) + ((double)j-(double)ny/2.0)*((double)j-(double)ny/2.0)); //finds radius
			//have number of pixels away from centre
			r2 *= PAR->cube_scale/nx; //now in AU
			r = (int)floor((r2 - PAR->rad_bins_inner)/dr); //which radial bin is it in?
			if( (r >= 0) && (r < PAR->n_rad_bins) ){ //make sure it's within bin ranges
				(*dust_cube)[i][j][0] = rad_bins[r]/(2.0*M_PI*r*dr); //smooth with radial bins (take account of averaging?)
				//(*dust_cube)[i][j][0] = r;
			}
			else{
				(*dust_cube)[i][j][0] = 0.0; //nothing outside
			}
		}
	}
	free(rad_bins);
	return;
}

/* write_dust_cube */
void write_dust_cube(struct parameters *PAR, char* outfile, int nx, int ny, int nz, double ***dust_cube, char* infile){
	long naxis = 3; /*dimension of image*/
	int i, j, k, l;
	long first_pixel;
	long nelements=(nx*ny*nz);
	fitsfile *fitsptr;
	int status=0;
	int bitpix = DOUBLE_IMG; /*image is in double format*/

	long fpixel[3] = {1, 1, 1};
	long naxes[3] = {nx, ny, nz};

	double *out_array; /*need to put into 1d format as cfits likes to be annoying*/
	char *bunit="mJy", *ctype="offset (arcseconds)";
	float crpix1=nx/2.0, crval=0.0, cdelt1=(PAR->cube_scale*PAR->image_scale)/nx;
	float crpix2=ny/2.0, cdelt2=(PAR->cube_scale*PAR->image_scale)/ny;

	out_array = (double*) malloc( (int)(nx*ny*nz*sizeof(double)) );

	if(out_array==NULL){
		printf(MALLOC_ERROR, "out_arrray");
		exit(1);
	}
	
	k=0;
	for(l=0;l<nz;l++){
		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
					out_array[k] = dust_cube[i][j][l];
					k++;			
			}
		}
	}
	printf("Creating FITS file...\n");
	if( !fits_create_file(&fitsptr, outfile, &status) ){
		/*Write the fits file here,*/
		
		first_pixel = 1; 							/*first pixel to write*/
		
		printf("Constants found, creating FITS image...\n");
		fits_create_img(fitsptr,  bitpix, naxis, naxes, &status);
		if(status){
			printf("Error: There was a problem creating the image to the fits file.\n");
			return;
		}
		printf("Writing FITS image...\n");
		/*fits_write_img(fitsptr, TDOUBLE, first_pixel, nelements, dust_cube[0], &status);
		*/
		fits_write_pix(fitsptr, TDOUBLE, fpixel, nelements, out_array, &status);
		if(status){
			printf("Error: There was a problem writing the image to the fits file.\n");
			return;
		}
		printf("Writing header data...\n");
		fits_update_key(fitsptr, TSTRING, "INPUT_FILE", infile,
			"File this cube was created from", &status);
		fits_update_key(fitsptr, TSTRING, "BUNIT", bunit,
			"", &status);
		fits_update_key(fitsptr, TSTRING, "CTYPE1", &ctype,
			"", &status);
		fits_update_key(fitsptr, TSTRING, "CTYPE2", &ctype,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CRPIX1", &crpix1,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CRPIX2", &crpix2,
			"", &status);
			
		fits_update_key(fitsptr, TFLOAT, "CRVAL1", &crval,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CRVAL2", &crval,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CDELT1", &cdelt1,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CDELT2", &cdelt2,
			"", &status);
         if(status){
			printf("Error: There was a problem writing the origin file to the header.\n");
			return;
		}
		printf("Closing FITS file...\n");
		if(fits_close_file(fitsptr, &status)){
			printf("Warning: Could not close file successfully\n");
		}
	}
	if(status){
		printf("Error: There was a problem writing the fits file, it may be corrupted.\n");
	}
	return;
}

void txt_dust_cube(struct parameters *PAR, double **dust_cube){
	/*Only really use this for debugging*/
	FILE *f;
	int i, j, k;
	
	f = fopen(PAR->fitsname, "w");
	
	for(i=0;i<PAR->cube_res;i++){
		for(j=0;j<PAR->cube_res;j++){
			fprintf(f, "%G ", dust_cube[i][j]);
		}
		fprintf(f, "\n");
	}

	close(f);

}



/*initialise the values of the parameter array, do any error checking here if possible*/
void init_par(struct parameters *PAR){
	PAR->image_scale = 170.0/3.6; /*defines the scale of the image*/
	PAR->cube_res = 100.0; /*resolution in x, y direction*/
	PAR->cube_scale = 10.0; /*the size of the image in the x-y direction (AU)*/
	PAR->cube_side = PAR->cube_scale/PAR->cube_res; /*the side length of a cube-cell*/
	PAR->cube_z_scale = 0.1; /*the size of the image in z-direction (AU)*/
	PAR->cube_z_res = PAR->cube_z_scale/PAR->cube_side; /*the resolution of the z-direction, fixed so that grid is square*/
	/*here I have to make sure that cube_z_res is an integer, and at least 1*/
	if(PAR->cube_z_res < PAR->cube_z_scale){
		PAR->cube_z_res = PAR->cube_z_scale;
		PAR->cube_z_scale = PAR->cube_side;
	}
	else if(PAR->cube_z_res != floor(PAR->cube_z_res)){
		PAR->cube_z_res = floor(PAR->cube_z_res);
		PAR->cube_z_scale = PAR->cube_z_res*PAR->cube_side;
	}
	PAR->cube_min = PAR->cube_scale/2.0 - PAR->cube_scale; /*minimum value in x-y direction*/
	PAR->cube_max = PAR->cube_scale/2.0; /*maximum value in x-y direction*/
	PAR->cube_z_min = PAR->cube_z_scale/2.0 - PAR->cube_z_scale; /* minimum value in z-direction*/
	PAR->cube_z_max = PAR->cube_z_scale/2.0; /*maximum value in z-direction*/
	PAR->orb_res = 10000.0; /*number of 'blobs' to split planetesimals into*/
	PAR->smallest_mass = (4.7*5.97E24/2.0E30)/1E6; /*the smallest mass a planetesimal can be*/
	/*want the mass of a 10nm dust grain in kg: 4/3*pi*r^3*density = mass*/
	PAR->smallest_dust = ((4.0/3.0)*3.142*(1E-8*1E-8*1E-8)*2000.0);
	PAR->destruction_radius = 0.5; /*dust destruction radius in AU*/
	PAR->fitsname = "outtest.fits"; /*default name for fits file*/
	PAR->gridname = "amr_grid.inp"; /*default name for grid output*/
	PAR->dataname = "dust_density.inp"; /*default name for data output*/
	PAR->dust_rad_min = 1E-7; /*Minimum radius of dust to use in meters*/
	PAR->dust_rad_max = 1E-3; //Maximum radius of dust to use in meters
	PAR->debris_res = 1E3; //Maximum size of stuff that isn't a particle
	PAR->debris_blowout = 1E-7; //radiation pressure blow out limit (smallest size)
	PAR->dust_factor = 10.0; /*factor either side that makes dust bin*/
	//PAR->dust_min = ((4.0/3.0)*3.142*((PAR->dust_size/PAR->dust_factor)*(PAR->dust_size/PAR->dust_factor)*(PAR->dust_size/PAR->dust_factor))*2000.0); /*smallest dust in bin*/
	//if(PAR->dust_min < PAR->smallest_dust) 
	//	PAR->dust_min = PAR->smallest_dust; /*cannot be smaller than grain blow-out size*/
	//PAR->dust_max = ((4.0/3.0)*3.142*(PAR->dust_size*PAR->dust_factor*PAR->dust_size*PAR->dust_factor*PAR->dust_size*PAR->dust_factor)*2000.0); /*largest dust in bin*/
	PAR->dust_norm = SOLAR_MASS*2E-7; //mass to normalise dust output to, use 0.0000002 solar masses as default
	
	//PAR->scalars = NULL; //this is a array of cubes with same dimensions as dust-cube, to be point-wise multiplied i.e. scales the output
	//PAR->n_scalars = 0; //this is the number of different cubes 
	PAR->n_rad_bins = 15; //number of radial bins to use when azimuthally smoothing
	PAR->rad_bins_inner = 0.7; //where to start the bins (in AU, will be continuous)
	PAR->rad_bins_outer = 3.9; //where to end the bins (int AU)
	PAR->rad_bins_file = "!radtest.fits"; //name default name of resulting fits file
	PAR->dustopac = "/Users/glzjd/Documents/data/d_opac_store/dustopac.inp"; //path to dustopac.inp file (default is same folder)
	//PAR->dustopac = ""; //path to dustopac.inp file, set to nothing to turn off
	puts("Parameter structure initialised...");
	
	printf("\nimage_scale: %G\ncube_res : %G\ncube_scale: %G\ncube_side: %G\ncube_z_scale: %G\n\
	\rcube_z_res: %G\ncube_min: %G\ncube_max: %G\ncube_z_min: %G\ncube_z_max: %G\norb_res: %G\n\
	\rsmallest_mass: %G\ndestruction_radius: %G\nfitsname: %s\ngridname: %s\ndataname: %s\n\
	\rdust_rad_min: %G\ndust_rad_max: %G\ndust_factor: %G\ndust_min: %G\ndust_max: %G\ndust_norm: %G\n\n",\
	PAR->image_scale, PAR->cube_res, PAR->cube_scale, PAR->cube_side, PAR->cube_z_scale, \
	PAR->cube_z_res, PAR->cube_min, PAR->cube_max, PAR->cube_z_min, PAR->cube_z_max,\
	PAR->orb_res, PAR->smallest_mass, PAR->destruction_radius, PAR->fitsname, PAR->gridname,\
	PAR->dataname, PAR->dust_rad_min, PAR->dust_rad_max, PAR->dust_factor, PAR->dust_min, PAR->dust_max, PAR->dust_norm);
	
	return;
}
