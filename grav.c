#ifndef STDIO_INCLUDED
	#include <stdio.h>
	#define STDIO_INCLUDED
#endif
#ifndef P_TREE_INCLUDED
	#include "p_tree.h"
	#define P_TREE_INCLUDED
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

typedef struct g_cell{
	double pos_bbl[N_DIM]; //back bottom left position
	double pos_tfr[N_DIM]; //top front right position
	double pos_com[N_DIM]; //centre of mass
	double eff_m; //effective mass (usually the sum)
	PARTICLE **obj; //list of the objects inside it
	int n_obj; //number of objects inside it
} G_CELL;

create_gcells(NODE* root, G_CELL **g_cells, int N);

int main(int argc, char **argv){
	SSDATA *ssdata;
	SSHEAD *sshead;
	PARTICLE **data;
	NODE *root;
	G_CELL *g_cells;
	int N;
	
	sshead = (SSHEAD*) malloc( (int)sizeof(SSHEAD) );
	ssdata = readin_ss(argv[1], sshead);
	N = sshead->n_data
	data = (PARTICLE**) malloc( (int)sizeof(PARTICLE*)*sshead->n_data );
	
	//copy all the information to 'data' array
	MALLOC_CHECK(data);
	
	for(i=0;i<sshead->n_data;i++){//copy data over
		data[i] = (PARTICLE*) malloc( (int)sizeof(PARTICLE) );
		MALLOC_CHECK(data[i]);
		for(j=0;j<N_DIM;j++){
			data[i]->pos[j] = ssdata[i].pos[j];
			data[i]->vel[j] = ssdata[i].vel[j];
			data[i]->spin[j] = ssdata[i].spin[j];
		}
		data[i]->mass = ss_data[i].mass;
		data[i]->radius = ss_data[i].radius;
		data[i]->ident = i;
	}
	//free unused storage
	FREE(ssdata);
	FREE(sshead);
	puts("ss_data freed");
	
	//allocate the root node, the rest will be allocated later
	root = (NODE*) malloc( sizeof(NODE) );
	create_tree(&root, data, N, N_DIM);
	create_gcells(root, &g_cells, N);
	
	return;	
}

create_gcells(NODE* node, G_CELL **g_cells, int N){
	int max_n = 8;
	int n_cells = N/max_n +1; // 10/8=1 in integer division but will need 2 cells...
	
	
	*g_cells = (G_CELL*) malloc( n_cells*sizeof(G_CELL) );
	//move through the cells from the left. 
	while(node->left!=NULL){
		node = node->left;
	}
	

}






