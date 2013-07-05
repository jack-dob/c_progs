#ifndef RADMC3D_CONV_INCLUDED
#define RADMC3D_CONV_INCLUDED
#include "mass_bin_dust.h" //make sure we can read in dustopac.inp files

void radmc3d_grid(struct parameters *PAR);
void radmc3d_data(struct parameters *PAR, double ***dust_cube);
void radmc3d_grid_b(struct parameters *PAR);
void radmc3d_data_b(struct parameters *PAR, double ***dust_cube);
void radmc3d_dust(struct parameters *PAR, double ***dust_cube, FILE *dustopac);
void unique_sizes(double *dust_rad, double *dust_factor, int n_dust, double **udust_rad, double **udust_reps, int *n_urad);
void find_size_bins(double **dust_limits, double* udust_rad, double* udust_reps, int n_urad, int n_dust);

void radmc3d_grid(struct parameters *PAR){
	int i;
	int iformat, grid_style, co_ord, grid_info, inclx, incly, inclz, nx, ny, nz;
	double *wall_x, *wall_y, *wall_z;
	FILE *f;
	/*"amr_grid.inp"; /*default name*/

	if (PAR->dataname[9]=='b'){
		radmc3d_grid_b(PAR);
		return;
	}

	iformat = 1; /*the first format they've come up with*/
	grid_style = 0; /*0-> rectangular, 1-> oct tree, 2-> layer*/
	co_ord = 0; /*000 ->cartesian, 100->sperical polar, 200->cylindircal polars*/
	grid_info = 0; /*is extra information included? almost always =0*/
	inclx = 1; /*include x dimension?*/
	incly = 1; /*include y dimension?*/
	inclz = 1; /*include z dimension?*/
	nx = PAR->cube_res; /*number of grid points in x*/
	ny = PAR->cube_res; /*number of grid points in y*/
	nz = PAR->cube_z_res; /*number of grid points in z*/
	
	
	wall_x = (double*) malloc( (int)(PAR->cube_res+1)*sizeof(double) ); /*position of cell walls in x*/
	wall_y = (double*) malloc( (int)(PAR->cube_res+1)*sizeof(double) ); /*position of cell walls in y*/
	wall_z = (double*) malloc( (int)(PAR->cube_z_res+1)*sizeof(double) ); /*position of cell walls in z*/
	
	/*Have to change from AU to cm as radmc uses CGS*/
	for(i=0;i<(PAR->cube_res+1);i++){
		wall_x[i] = PAR->cube_min*100.0*AU + (double)i*100.0*AU*PAR->cube_scale/PAR->cube_res;
		wall_y[i] = PAR->cube_min*100.0*AU + (double)i*100.0*AU*PAR->cube_scale/PAR->cube_res;
	}
	for(i=0;i<(PAR->cube_z_res+1);i++){
		wall_z[i] = PAR->cube_z_min*100.0*AU + (double)i*100.0*AU*PAR->cube_z_scale/PAR->cube_z_res;
	}
	
	/***	Print Header Data	***/
	f = fopen(PAR->gridname, "w");
	
	fprintf(f, "%i\n", iformat);
	fprintf(f, "%i\n", grid_style);
	fprintf(f, "%i\n", co_ord);
	fprintf(f, "%i\n", grid_info);
	fprintf(f, "%i %i %i\n", inclx, incly, inclz);
	fprintf(f, "%i %i %i\n", nx, ny, nz);
	
	/***	Print Data	***/
	for(i=0;i<(PAR->cube_res+1);i++){
		fprintf(f, "%E ", wall_x[i]);
	}
	fprintf(f, "\n");
	for(i=0;i<(PAR->cube_res+1);i++){
		fprintf(f, "%E ", wall_y[i]);
	}
	fprintf(f, "\n");
	for(i=0;i<(PAR->cube_z_res+1);i++){
		fprintf(f, "%E ", wall_z[i]);
	}
	fprintf(f, "\n");
	
	close(f);
	free(wall_x);
	free(wall_y);
	free(wall_z);
	
	return;
}

void radmc3d_data(struct parameters *PAR, double ***dust_cube){
	//Things get tricky here. This file, the non-binary file needs actual densities
	//so we have to divide the mass of dust in dust_cube by the volume of each cube
	//and make sure we are in CGS!!
	int i,j,k;
	int iformat, nrcells, nrspec;
	double vol=(100.0*AU*PAR->cube_side)*(100.0*AU*PAR->cube_side)*(100.0*AU*PAR->cube_side);
	FILE *f;
	/*"dust_density.inp"; /*default name*/
	
	if (PAR->dataname[13]=='b'){
		radmc3d_data_b(PAR, dust_cube);
		return;
	}
	
	
	iformat = 1; /*Format of file*/
	nrcells = (int)(PAR->cube_res*PAR->cube_res*PAR->cube_z_res); /*number of cells in total*/
	nrspec = 1; /*number of dust species (i.e. size, composition)*/
	
	/***	Print Header Data	***/
	f = fopen(PAR->dataname, "w");
	
	fprintf(f, "%i\n", iformat);
	fprintf(f, "%i\n", nrcells);
	fprintf(f, "%i\n", nrspec);
	
	/***	Print Data	***/
	for(k=0;k<PAR->cube_z_res;k++){
		for(j=0;j<PAR->cube_res;j++){
			for(i=0;i<PAR->cube_res;i++){
				dust_cube[i][j][k] *= 1000.0/vol;
				fprintf(f, "%E\n", dust_cube[i][j][k]);
			}
		}
	}
	close(f);
	return;
}

void radmc3d_grid_b(struct parameters *PAR){
	/*Binary version of above function*/
	/*radmc uses CGS!! i've been using SI, must convert*/
	int i;
	long int iformat, grid_style, co_ord, grid_info, inclx, incly, inclz, nx, ny, nz;
	double *wall_x, *wall_y, *wall_z;
	FILE *f;
	/*"amr_grid.binp"; /*default name*/
	
	if (PAR->gridname[9]=='i'){
		radmc3d_grid(PAR);
		return;
	}
	
	printf("sizeof long int: %i, size of double: %i\n", sizeof(long int), sizeof(double));
	if((int)sizeof(long int)!=8){
		printf("ERROR: Have to have 8 byte long integers, will need to re-write this bit in ftn %s\n", __FUNCTION__);
		exit(1);
	}
	if((int)sizeof(double)!=8){
		printf("ERROR: Have to have 8 byte doubles, will need to re-write this bit in ftn %s\n", __FUNCTION__);
		exit(1);
	}
	
	puts("In 'radmc3d_grid_b', assigning constants");
	
	iformat = 1; /*the first format they've come up with*/
	grid_style = 0; /*0-> rectangular, 1-> oct tree, 2-> layer*/
	co_ord = 0; /*000 ->cartesian, 100->sperical polar, 200->cylindircal polars*/
	grid_info = 0; /*is extra information included? almost always =0*/
	inclx = 1; /*include x dimension?*/
	incly = 1; /*include y dimension?*/
	inclz = 1; /*include z dimension?*/
	nx = PAR->cube_res; /*number of grid points in x*/
	ny = PAR->cube_res; /*number of grid points in y*/
	nz = PAR->cube_z_res; /*number of grid points in z*/
	
	puts(" Constants assigned, allocating memory... ");
	
	wall_x = (double*) malloc( (int)(PAR->cube_res+1)*sizeof(double) ); /*position of cell walls in x*/
	wall_y = (double*) malloc( (int)(PAR->cube_res+1)*sizeof(double) ); /*position of cell walls in y*/
	wall_z = (double*) malloc( (int)(PAR->cube_z_res+1)*sizeof(double) ); /*position of cell walls in z*/
	
	for(i=0;i<(PAR->cube_res+1);i++){
		/*factor of 100 is due to conversion from meters to centimeters*/
		wall_x[i] = PAR->cube_min*100.0 + (double)i*100.0*PAR->cube_scale/PAR->cube_res;
		wall_y[i] = PAR->cube_min*100.0 + (double)i*100.0*PAR->cube_scale/PAR->cube_res;
	}
	for(i=0;i<(PAR->cube_z_res+1);i++){
		wall_z[i] = PAR->cube_z_min*100.0*AU + (double)i*100.0*PAR->cube_z_scale/PAR->cube_z_res;
	}
	
	puts("Memory allocated and assigned, opening file and writing output");
	
	/***	Output Header Data	***/
	f = fopen(PAR->gridname, "w");
	
	fwrite(&iformat, sizeof(long int), 1, f);
	fwrite(&grid_style, sizeof(long int), 1, f);
	fwrite(&co_ord, sizeof(long int), 1, f);
	fwrite(&grid_info, sizeof(long int), 1, f);
	fwrite(&inclx, sizeof(long int), 1, f);
	fwrite(&incly, sizeof(long int), 1, f);
	fwrite(&inclz, sizeof(long int), 1, f);
	fwrite(&nx, sizeof(long int), 1, f);
	fwrite(&ny, sizeof(long int), 1, f);
	fwrite(&nz, sizeof(long int), 1, f);
	
	/***	Output Data	***/
	fwrite(wall_x, sizeof(double), (PAR->cube_res+1), f);
	fwrite(wall_y, sizeof(double), (PAR->cube_res+1), f);
	fwrite(wall_x, sizeof(double), (PAR->cube_z_res+1), f);
	
	puts("Data output, freeing arrays...");

	close(f);
	free(wall_x);
	free(wall_y);
	free(wall_z);

	puts("Arrays freed, returning.");

	return;
}

void radmc3d_data_b(struct parameters *PAR, double ***dust_cube){
	/*Binary version of above function*/
	/*radmc uses CGS!! i've been using SI, must convert*/
	/*this time the data needs to be in MASS format, NOT DENSITY!!*/
	int i,j,k;
	long int iformat, precis, nrcells, nrspec;
	double dist;
	FILE *f;
	/*"dust_density.binp"; /*default name*/
	
	if (PAR->dataname[13]=='i'){
		radmc3d_data(PAR, dust_cube);
		return;
	}
	
	
	if(sizeof(long int)!=8){
		printf("ERROR: Have to have 8 byte long integers, will need to re-write this bit in ftn %s\n", __FUNCTION__);
		exit(1);
	}
	if(sizeof(double)!=8){
		printf("ERROR: Have to have 8 byte doubles, will need to re-write this bit in ftn %s\n", __FUNCTION__);
		exit(1);
	}
	
	puts("In 'radmc3d_data_b' assigning constants");
	
	iformat = 1; /*Format of file*/
	precis = 8; /*precision of numbers*/
	nrcells = PAR->cube_res*PAR->cube_res*PAR->cube_z_res; /*number of cells in total*/
	nrspec = 1; /*number of dust species (i.e. size, composition)*/
	
	puts("Constants assigned, writing constants");
	
	/***	Print Header Data	***/
	f = fopen(PAR->dataname, "w");
	
	fwrite(&iformat, sizeof(long int), 1, f);
	fwrite(&precis, sizeof(long int), 1, f);
	fwrite(&nrcells, sizeof(long int), 1, f);
	fwrite(&nrspec, sizeof(long int), 1, f);

	puts("Constants written, writing data...");
	
	/***	Print Data	***/
	/*fwrite(dust_cube, sizeof(double), nrcells, f); /*Does this work?*/
	/* position of star is CUBE_REZ/2.0, PAR->cube_res/2.0, PAR->cube_z_res/2.0*/
	/* want to make sure there's not loads of dust right next to star*/
	/* set exclusion zone of 0.3AU (inner orbit of mercury)*/
	/* PAR->cube_scale/PAR->cube_res = size of cube */
	/* scale/res(res/2 - co-ord) = dist in a direction */
	for(k=0;k<PAR->cube_z_res;k++){ /* z-coord */
		for(j=0;j<PAR->cube_res;j++){ /* y-coord */
			for(i=0;i<PAR->cube_res;i++){ /* x-coord */
				/*printf("address of data [%i][%i][%i]: %i\n", i, j, k, &(dust_cube[i][j][k]));*/
				if( isnan(dust_cube[i][j][k]) || isinf(dust_cube[i][j][k]) ){
					printf("ERROR: Data is corrupted, dust_cube[%i][%i][%i]: %E\n", i, j, k, dust_cube[i][j][k]);
					exit(1);
				}
				/*we have mass in SI, want density in CGS*/
				dust_cube[i][j][k] *= 1000.0;//dispite what file is called, is actually just mass
				fwrite(&(dust_cube[i][j][k]), sizeof(double), 1, f);
			}
		}
	}
	close(f);
	puts("Data written, returning.");
	
	return;
}

void radmc3d_dust(struct parameters *PAR, double ***dust_cube, FILE *dustopac){
	//use mass_bin_dust.h fuctions to get the data needed to find the correct mass 
	//to put into each dust species, then write each species to file
	int n_dust, i, j, k, n;
	double *dust_rad, *dust_factor, **dust_limits, *udust_rad, *udust_reps;
	double species_factor; //what I need to multiply the normalised dust by to get correct mass for each species
	double constant; //needed for species_factor
	int iformat, nrcells, nrspec;
	int n_urad;
	double fac;
	double vol=(100.0*AU*PAR->cube_side)*(100.0*AU*PAR->cube_side)*(100.0*AU*PAR->cube_side);
	FILE *f;
	//"dust_density.inp"; //default name
	
	read_dustopac(dustopac, &n_dust, &dust_rad, &dust_factor);
	close(dustopac);
	puts("read dustopac.inp file...");
	if(n_dust==1){
		puts("Only one dust species, do it manually...");
		exit(0);
	}
	// problem, what if I have two different species at the same size e.g. carbon and sillicon?
	//I will assume that any entry with a dust_factor below 1 will have things the same size as
	//it following it
	unique_sizes(dust_rad, dust_factor, n_dust, &udust_rad, &udust_reps, &n_urad);
	
	puts("Finding dust size bin limits...");
	dust_limits = (double**) malloc( n_dust*sizeof(double*) );
	for(i=0;i<n_dust;i++){
		dust_limits[i] = (double*) malloc( 2*sizeof(double) ); //this will hold upper and lower limits of dust size bins
	}
	
	find_size_bins(dust_limits, udust_rad, udust_reps, n_urad, n_dust);
	
	puts("Writing radmc data file...");
	//Things get tricky here. This file, the non-binary file needs actual densities
	//so we have to divide the mass of dust in dust_cube by the volume of each cube
	//and make sure we are in CGS!!
	
	
	iformat = 1; //Format of file
	nrcells = (int)(PAR->cube_res*PAR->cube_res*PAR->cube_z_res); //number of cells in total
	nrspec = n_dust; //number of dust species (i.e. size, composition)
	
	//***	Print Header Data	***
	f = fopen(PAR->dataname, "w");
	
	fprintf(f, "%i\n", iformat);
	fprintf(f, "%i\n", nrcells);
	fprintf(f, "%i\n", nrspec);
	
	//constant should be defined in CGS!
	constant = sqrt(PAR->debris_res*1E2) - sqrt(PAR->debris_blowout*1E2);
	puts("Header written, dumping data");
	//***	Print Data	***
	//printf("Writing dust species");
	for(n=0;n<n_dust;n++){
		//make sure all is in CGS!!
		//printf(", %d", n);
		// (1000.0/vol) factor makes sure we're in CGS
		//this is where the dust mass is worked out from the size distribution
		//dust limits is in MICROMETERS!! need to change to METERS!! then to CGS!!
		species_factor = (1000.0/vol)*dust_factor[n]*(sqrt(1E-6*dust_limits[n][1]*1E2) - sqrt(1E-6*dust_limits[n][0]*1E2))/constant;
		printf("dust_limits[%d]: %G to %G\n", n, dust_limits[n][0], dust_limits[n][1]);
		printf("constant: %G\n", constant);
		printf("species_factor[%d]: %G\n", n, species_factor);
		for(k=0;k<PAR->cube_z_res;k++){
			for(j=0;j<PAR->cube_res;j++){
				for(i=0;i<PAR->cube_res;i++){
					//put dust_cube into GCS mass density!
					fprintf(f, "%E\n", species_factor*dust_cube[i][j][k]);
				}
			}
		}	
	}
	printf("\n");
	close(f);
	//free malloc'd stuff
	free(dust_rad);
	free(dust_factor);
	for(i=0;i<n_dust;i++){
		free(dust_limits[i]);
	}
	free(dust_limits);
	free(udust_rad);
	free(udust_reps);
	return;
}

void unique_sizes(double *dust_rad, double *dust_factor, int n_dust, double **udust_rad, double **udust_reps, int *n_urad){
	int i, j, k;
	double fac_sum=0.0;
	*n_urad = 0; //assume no unique sizes to start
	for(i=0;i<n_dust;i++){
		fac_sum += dust_factor[i];
		if(fac_sum >1.01){
			printf("Error: Dust factors greater than one within one percent tolerance\n");
			for(k=0;k<n_dust;k++){
				printf("i: %d, dust_factor[%d]: %G\n", i, k, dust_factor[k]);
			}
			exit(0);
		}
		if(fac_sum > 0.99){
			(*n_urad)++; //increase number of unique sizes found
			fac_sum = 0.0; //reset sum
		}	
	}
	(*udust_rad) = (double*) malloc( *n_urad*sizeof(double)); //allocate storage for unique sizes
	(*udust_reps) = (double*) calloc( *n_urad, sizeof(double) );//storage for number of repeats
	j=0;
	for(i=0;i<n_dust;i++){
		fac_sum += dust_factor[i];
		(*udust_reps)[j]++; //have an extra one of these
		printf("fac_sum: %G, udust_reps[%d]: %d\n", fac_sum, i, (*udust_reps)[i]);
		//if(fac_sum >1.01){
		//	printf("Error: Dust factors greater than one within one percent tolerance\n");
		//	for(k=0;k<n_dust;k++){
		//		printf("i: %d, dust_factor[%d]: %G\n", i, k, dust_factor[k]);
		//	}
		//	exit(0);
		//}
		if(fac_sum > 0.99){
			(*udust_rad)[j] = dust_rad[i]; //record unique size
			j++; //move to space for new record
			fac_sum = 0.0; //reset sum
		}
	}
	return;
}

void find_size_bins(double **dust_limits, double* udust_rad, double* udust_reps, int n_urad, int n_dust){
	int i=0, j=0, k=0;
	double low_lim=0.0, up_lim=0.0;
	k=0;
	for(i=0;i<n_urad;i++){//loop over all the unique dust sizes
		printf("udust_rad[%d]: %G, udust_reps[%d]: %G\n", i, udust_rad[i], i, udust_reps[i]);
		if(i==0){ //set the first one's lower limit
			low_lim = udust_rad[i]/2.0; //half the size of the smallest one
		}
		else{
			low_lim = (udust_rad[i-1]+udust_rad[i])/2.0;
		}
		if(i==(n_urad-1)){//set the last one's upper limit
			up_lim = udust_rad[i]*2.0; //twice the size of the biggest one
		}
		else{//set the other one's limits
			up_lim = (udust_rad[i]+udust_rad[i+1])/2.0;
		}
		for(j=0;j<udust_reps[i];j++){
			dust_limits[k][0] = low_lim;
			dust_limits[k][1] = up_lim;
			k++;
		}
	}
	return;
}

#endif //RADMC3D_CONV_INCLUDED