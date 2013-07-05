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
#ifndef FREE
	#define FREE(pointer) free(pointer); pointer = NULL;
#endif
#ifndef MALLOC_CHECK
	#define MALLOC_CHECK(pointer)\
		if( (pointer)==NULL ){\
			printf("Malloc Error: Line %d, Function %s, File %s\n", __LINE__, __func__, __FILE__);\
			exit(1);\
		}
#endif	
//#define MALLOC_ERROR(string) fprintf(stderr, "Error: Could not allocate '%s' in file: %s, function: %s, line: %s\n", string, __FILE__, __func__, __LINE__)

typedef struct particle{
	double	mass;
	double	radius;
	double	pos[N_DIM];
	double	vel[N_DIM];
	double	spin[N_DIM];
	int		ident;
	int		colour;
} PARTICLE;
#define TREE_TYPE PARTICLE //controls which type of object mem_tree.h operates on
#define TREE_SORT pos
#define N_DIM 3
#ifndef MEM_TREE_INCLUDED
	#include "mem_tree.h"
	#define MEM_TREE_INCLUDED
#endif

// function prototypes //
void step_2_au(double* point, int i, int j, int size, double scale);
void nn_density(NODE **n_near, double *n_dist2, double scale_len, double* location, int n_dim, int n);
void nn_num_density(NODE **n_near, double *n_dist2, double sigmas, double *location, int n_dim, int n);
double max(double *array, int n);

int main(int argc, char **argv){
	NODE *root = NULL, **n_near = NULL;
	SSDATA *ss_data = NULL;
	SSHEAD head;
	PARTICLE **data = NULL, *point = NULL;
	double *n_dist2 = NULL;
	double **density = NULL; //maybe make this 3d?
	double scale = 10.0; //size of each side in AU
	double sigmas = 5.0; //affects the gaussian dropoff of the density calculation, is used as the number of standard divs the farthest neighbour is away
	int size = 255, n = 32, N;//size = resolution, n=number of neighbours, N=number of particles
	int i=0, j=0, k=0;
	FILE* f;
	char* out_f_name = "data.txt";
	
	ss_data = readin_ss(argv[1], &head);
	N = head.n_data-1;
	data = (PARTICLE**) malloc( (int)sizeof(PARTICLE*)*N );
	
	//copy all the information to 'data' array
	MALLOC_CHECK(data);
	
	for(i=1;i<head.n_data;i++){//start at 1 to cut out jupiter mass planet
		data[i-1] = (PARTICLE*) malloc( (int)sizeof(PARTICLE) );
		MALLOC_CHECK(data[i-1]);
		for(j=0;j<N_DIM;j++){
			data[i-1]->pos[j] = ss_data[i].pos[j];
			data[i-1]->vel[j] = ss_data[i].vel[j];
		}
		data[i-1]->mass = ss_data[i].mass;
		data[i-1]->ident = i;
		data[i-1]->colour = ss_data[i].color;
	}
	//free unused storage
	FREE(ss_data);
	puts("ss_data freed");
	
	root = (NODE*) malloc((int)sizeof(NODE));
	MALLOC_CHECK(root);
	
	density = (double**) malloc( (int)sizeof(double*)*size );
	MALLOC_CHECK(density);
	
	for(i=0;i<size;i++){
		density[i] = (double*) malloc( (int)sizeof(double)*size );
		MALLOC_CHECK(density[i]);
	}
	
	//for(i=0;i<head.n_data;i++){ //gets data into useful format
	//	data[i] = ss_data[i].pos;
	//}
	
	create_tree(&root, data, N, N_DIM); //create tree of the data
	
	puts("Allocating memory for query point.");
	point = (PARTICLE*) malloc( (int)sizeof(PARTICLE) );
	MALLOC_CHECK(point);

	n_dist2 = (double*) malloc( (int)sizeof(double)*n );
	MALLOC_CHECK(n_dist2);
	n_near = (NODE**) malloc( (int)sizeof(NODE*)*n );
	MALLOC_CHECK(n_near);
	
	f = fopen(out_f_name, "w");
	
	init_nn_list(n_near, n_dist2, n);
	puts("finding density on grid.");
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			step_2_au(point->pos, i, j, size, scale);
			//init_nn_list(n_near, n_dist2, n);
			//recompute_nn_dist(n_near, n_dist2, point, n, N_DIM);
			nearest_nn(root, n_near, n_dist2, point, n, 0, N_DIM);
			//if((i==0 && j==0) || (i==0 && j==50)){ //for debugging
			//	printf("i: %d, j: %d\n", i, j);
			//	for(k=0;k<n;k++){
			//		printf("%G\t%G\t%G\n", n_near[k]->val[0], n_near[k]->val[1], n_dist2[k]);
			//	}
			//}
			nn_num_density(n_near, n_dist2, sigmas, &(density[i][j]), N_DIM, n);
			fprintf(f, "%E\t", density[i][j]);
		}
		fprintf(f, "\n");
	}
	printf("finished grid density loop, output file: %s\n", out_f_name);
	//print_tree(root, 0, 0);
	//free all the memory
	puts("burning tree...");
	destroy_tree(root);
	puts("..tree burnt.");
	puts("freeing density.");
	for(i=0;i<size;i++){
		FREE(density[i]);
	}
	FREE(density);
	puts("density freed.");
	FREE(ss_data);
	puts("ss_data freed");
	for(i=0;i<N;i++){
		FREE(data[i]);
	}
	FREE(data);
	puts("data freed.");
	FREE(point);
	puts("point freed.");
}

void step_2_au(double* point, int i, int j, int size, double scale){
	(point)[0] = ((double)i/(double)size)*scale - scale/2.0;
	(point)[1] = ((double)j/(double)size)*scale - scale/2.0;
	(point)[2] = 0.0;//((double)j/(double)size)*scale - scale/2.0;
	return;
}

void nn_density(NODE **n_near, double *n_dist2, double sigmas, double* location, int n_dim, int n){
	//n_near: an array of nearest neighbours
	//n_dist2: an array of the squared distances of those nearest neighbours, in order
	//scale_len: the length scale of the smoothing, approximately the size of a cell
	//location: where the final answer will be stored (usually a cell in a 2d array)
	//n_dim: how many dimensions we have
	//n: how many nearest neighbours we have
	int i;
	double sum = 0.0, var2, factor;
	var2 = max(n_dist2, n)/(sigmas*sigmas); //assume 'max(n_dist2)' is at 'sigmas' sigma level
	factor = 1.0/(sqrt(var2)*sqrt(2.0*M_PI));
	for(i=0;i<n;i++){
		double gauss = factor*exp(-n_dist2[i]/(2.0*var2));
		sum += n_near[i]->val->mass*gauss; //assumes mass is in 'n_dim' slot of 'val' attribute
		//sum += particle[(n_near[i]->val[n_dim])]->mass*gauss //this could be a possible re-write
	}
	*location = sum;
}

void nn_num_density(NODE **n_near, double *n_dist2, double sigmas, double *location, int n_dim, int n){
	//This should find the number density associated with the n_near list
	int i;
	double sum = 0.0, var2, factor;
	var2 = max(n_dist2, n)/(sigmas*sigmas); //assume 'max(n_dist2) is at a 'sigmas' sigma level
	factor = 1.0/(sqrt(var2)*sqrt(2.0*M_PI));
	for(i=0;i<n;i++){
		sum += factor*exp(-n_dist2[i]/(2.0*var2)); //want number density, no scaling
	}
	*location = sum;
}

void nn_vel_disp(NODE **n_near, double n_dist2, double *ref_vel, double *location, int n_dim, int n){
	//finds the velocity dispersion of the nearest neighbours
	int i;
	double sum=0.0;
	for(i=0;i<n;i++){
		sum += sqrt(dot_prod(ref_vel,ref_vel) - dot_prod(n_near[i]->val->vel, n_near[i]->val->vel)); //assumes where the velocity starts
		//sum += sqrt(dot_prod(ref_vel, ref_vel) - dot_prod(particle[(n_near->val[n_dim])]->vel, particle[(n_near->val[n_dim])]->vel)); //possible re-write
	}
	*location = sum/(double)n;
}

double max(double *array, int n){
	int i;
	double max = 0.0;
	for(i=0;i<n;i++){
		if(array[i] > max){
			max = array[i];
		}
	}
	return(max);
}


