#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#ifndef TREE_TYPE
	#error pre-processor macro "TREE_TYPE" not defined, see line __LINE__ in __FILE__.
	//TREE_TYPE tells p_tree.h what type of data structure will be in the tree
	//e.g. #define TREE_TYPE double, or #define TREE_TYPE struct TREE_TYPE etc...
#endif

#ifndef TREE_SORT
	#error pre-processos macro "TREE_SORT" not defined, see line __LINE__ in __FILE__.
	//TREE_SORT tells p_tree.h which element of TREE_TYPE to sort over
	//e.g. if TREE_TYPE is 'TREE_TYPE', and TREE_SORT is 'TREE_SORT' then will sort over the 'TREE_SORT' attribute
#endif
//#ifndef N_DIM //have put this into the functions
//	#error pre-processos macro "NDIM" not defined, see line __LINE__ in __FILE__.
	//N_DIM tells p_tree.h how many members the sorting variable has
	//e.g. #define NDIM 3 means we have a 3 element array
//#endif

#define P_TREE_INCLUDED
#ifndef NULL_CHECK
	#define NULL_CHECK(array, str, name, i, N) \
		for(i=0;i<N;i++){\
			if(array[i]==NULL){\
			fprintf(stderr, "%sWarning: %s[%d] is NULL\n", str, name, i);\
			exit(0);\
			}\
		}
#endif
#ifndef MALLOC_CHECK
	#define MALLOC_CHECK(pointer)\
		if( (pointer)==NULL ){\
			printf("Malloc Error: Line %d, Function %s, File %s\n", __LINE__, __func__, __FILE__);\
			exit(1);\
		}
#endif
#ifndef FREE
	#define FREE(pointer)\
		free(pointer); pointer=NULL;
#endif


/*typedef struct TREE_TYPE{
	double	mass;
	double	radius;
	double	TREE_SORT[N_DIM];
	double	vel[N_DIM];
	double	spin[N_DIM];
	int		ident;
} TREE_TYPE;*/

//creates the structure 'node', represents a node in a tree
struct node{
	TREE_TYPE		*val; //value of the node
	struct node	*left, *right; //pointer to the children
	struct node *parent; //pointer to the parent node
};

typedef struct node NODE; //typedef the struct, easier to type

//*****function declarations*****//
void init_node(NODE* node);
void init_nn_list(NODE **nn_list, double *dist2_list, int n);
void recompute_nn_dist(NODE **nn_list, double *nn_dist2, double *point, int n, int n_dim);
void qsort_nd2(TREE_TYPE** data, int N, int dim);
int median_array(TREE_TYPE **array, int N, int dim);
int order_lt_gt(TREE_TYPE** array, int m, int N, int dim);
NODE* build_tree_array(TREE_TYPE **point_array, int N, int K, int depth);
void printbits(unsigned int v, unsigned int depth);
void print_tree(NODE *node, unsigned int TREE_SORT, unsigned int depth);
void link_parents(NODE *node);
double dist2(TREE_TYPE *a, double *b, int n_dim);
double dist2_d(double *a, double *b, int n_dim);
int any_gt(double* list, double test, int n);
int is_in(NODE** list, NODE* test, int n);
void nearest_nn(NODE* node, NODE** guess, double *g_dist2, TREE_TYPE* point, int n, int depth, int n_dim);
void nearest_n(NODE* node, NODE** guess, TREE_TYPE* point, double *g_dist2, int depth, int n_dim);
void create_tree(NODE** root, TREE_TYPE** data, int N, int K);
void destroy_tree(NODE* node);
void tree_insert(NODE* node, TREE_TYPE* point, int K, int depth);
void tree_insert_array(NODE* root, TREE_TYPE** data, int N, int K);
int find_in_array(TREE_TYPE **point_array, TREE_TYPE *point, int N);

//*****BEGIN CODE*****//

void init_node(NODE* node){
	//printf("\tin 'init_node' val %u, left %u, right %u, parent %u\n", node->val, node->left, node->right, node->parent);
	node->val = NULL;
	node->left = NULL;
	node->right = NULL;
	node->parent = NULL;
	//printf("\tin 'init_node' val %u, left %u, right %u, parent %u\n", node->val, node->left, node->right, node->parent);
	return;
} 

void init_nn_list(NODE **nn_list, double *nn_dist2, int n){
	int i;
	for(i=0;i<n;i++){ //set the nearest neighbours initially
		nn_dist2[i] = DBL_MAX; //maximum value a double can hold, should use infinity, but can't
		nn_list[i] = NULL; //nothing
	}
}

void recompute_nn_dist(NODE **nn_list, double *nn_dist2, double *point, int n, int n_dim){
	//re-computes the distances in an nearest-neighbour list for a new point
	int i;
	double temp;
	//puts("in 'recompute_nn_dist'");
	//printf("point: {%G, %G}\n", point[0], point[1]);
	for(i=0;i<n;i++){
		if(nn_list[i]==NULL){ //make sure null values are treated properly
			//printf("nn_list[%d] is NULL\n", i);
			nn_dist2[i]=DBL_MAX;
			continue;
		}
		temp = nn_dist2[i];
		//printf("nn_list[%d]: {%G, %G}, old_dist2: %G ", i, nn_list[i]->val[0], nn_list[i]->val[1], nn_dist2[i]);
		//printf("nn_list[%d]: %x\told_dist2: %G\t", i, nn_list[i], nn_dist2[i]);
		nn_dist2[i] = dist2_d(nn_list[i]->val->TREE_SORT, point, n_dim);
		//printf("new_dist: %G\tdiff: %G\n", nn_dist2[i], temp - nn_dist2[i]);
	}
}

void qsort_nd2(TREE_TYPE** data, int N, int dim){
	//need to pass a data array and a storage array. will junk the contents of both arrays
	//so must ensure that data array is sorted correctly
	int i, j, pivot, piv_id, temp_i;
	TREE_TYPE *temp;
	if(N<2) return; //if we have 1 or less TREE_TYPEs then the list is already sorted
	//printf("In qsort_nd2\narguments: N=%d, dim=%d\n", N, dim);
	//puts("data array:");
	//for(i=0;i<N;i++){
	//	printf("data[%d]->TREE_SORT[%d] = %G\t%u\n", i, dim, data[i]->TREE_SORT[dim], data[i]);
	//}
	
	pivot = N/2; //choose middle integer number as the pivot
	piv_id = data[pivot]->ident; //remember the pivot id
	//printf("pivot: data[%d]->TREE_SORT[%d] = %G \t %u\n", pivot, dim, data[pivot]->TREE_SORT[dim], data[pivot]);
	
	//j=0;
	for(i=0;i<N;i++){
		//printf("%07d, median %d: %G, test: %G\r", i, pivot, data[pivot]->TREE_SORT[dim], data[i]->TREE_SORT[dim]);
		if(data[i]->ident == piv_id)//don't sort the pivot against its self
			continue;
		if((data[i]->TREE_SORT[dim] < data[pivot]->TREE_SORT[dim]) && (i > pivot))//if data is smaller than pivot and we are after the pivot
		{
			temp = data[i]; //swap the data and pivot
			data[i] = data[pivot];
			data[pivot] = temp;
			temp_i = pivot; //remember old pivot TREE_SORTition
			pivot = i; //our pivot is now in a new TREE_SORTition
			i = temp_i; //must check the stuff after the old pivot TREE_SORTition again
		}
		else if((data[i]->TREE_SORT[dim] > data[pivot]->TREE_SORT[dim]) && (i < pivot))//if data is greater than or equal to pivot and we are before the pivot
		{
			temp = data[i]; //swap the data and pivot
			data[i] = data[pivot];
			data[pivot] = temp;
			pivot = i; //pivot is now at new TREE_SORTition, can set off from here as we moved the pivot to earlier in the list
		}
		else{ //if not last two then are equal so sort by a secondary metric
			continue; //secondary sorting goes here
		}
	}
	
	j = pivot+1; //j should be the element after the pivot
	
	//puts("storage array:");
	//for(i=0;i<N;i++){
	//	printf("storage[%d][%d] = %G\t%u\n", i, dim, storage[i][dim], storage[i]);
	//}
	
	//don't include the pivot in the new arrays, it's already in the correct place
	//puts("sort 'left' array.");
	if(pivot>1)//if more than one value in 'less' array
		qsort_nd2(data, pivot, dim); //sort the values and put in 'data' array
	
	//puts("sort 'right' array");
	//printf("sorted arrays, %d, %d\n", j, N);
	if(j<(N-1))// if more than one value is in 'greater' array
		qsort_nd2( (data+j), N-j, dim);//sort the values and put into 'data' array
		
	//now the arrays should be sorted
	//puts("SORTED array:");
	//for(i=0;i<N;i++){
	//	printf("data[%d][%d] = %G\t%u\n", i, dim, data[i][dim], data[i]);
	//}
	//puts("return from 'qsort_nd2'");
	
	return;
}

int median_array(TREE_TYPE **array, int N, int dim){
	int i;
	//returns the TREE_SORTition of the median value
	//puts("In 'median_array':");
	for(i=N/2; (array[i]->TREE_SORT[dim]-array[i-1]->TREE_SORT[dim])==0.0; i--){
		if(i==1){
			i=0;
			break;
		}
	}
	//printf("Median is array[%d]: {%G %G}\t N is %d\n", i, array[i][dim], array[i][dim], N);
	return(i);

}

int order_lt_gt(TREE_TYPE **array, int m, int N, int dim){
	//put 'array' into order of: stuff less than 'array[m]', 'array[m]', stuf greater than 'array[m]'
	//return TREE_SORTition of 'val[dim]'
	int i=0, j=N-1, k=0;
	int temp_i;
	TREE_TYPE **temp_p;
	temp_p = (TREE_TYPE**) malloc (N*sizeof(TREE_TYPE*));
	//puts("starting ordering loop\n");
	//printf("k:%d, j:%d, N:%d\n", k, j, N);
	for(i=0;i<N;i++){
		//if(array[i]==NULL) printf("\n\tNULL\n");
		//printf("%07d %% complete, median at %d: %G, test: %G\n", i, m, array[m]->TREE_SORT[dim], array[i]->TREE_SORT[dim]);
		if(array[i] == array[m]){ //if these point to the same place
			//printf("\tignoring value ident[%d] = %d, ident[%d] = %d\n", i, array[i]->ident, m, array[m]->ident);
			continue; //ignore it
		}
		//puts("\t1");
		else if((array[i]->TREE_SORT[dim] < array[m]->TREE_SORT[dim]) ){ //if data is less than value and we are after the median
			//puts("\tmoving earlier");
			temp_p[k] = array[i];//put at start of array
			k++; //next 'lower than' space
		}
		//puts("\t2");
		else if((array[i]->TREE_SORT[dim] >= array[m]->TREE_SORT[dim])){ //if data is greater than or equal to value and we are before the median
			//puts("\tmoving later");
			temp_p[j] = array[i];//put at end of array
			j--; //next 'higher than' space
		}
		else exit(0);
		//if(i==m+100) exit(0); //for debugging
	}
	//printf("k:%d, j:%d, N:%d\n", k, j, N);
	temp_p[k] = array[m]; //put median where it belongs
	//NULL_CHECK(temp_p, "", "temp_p", i, N);
	for(i=0;i<N;i++){
		array[i] = temp_p[i]; //copy to correct array
	}
	FREE(temp_p); //free array we created
	//puts("finished ordering loop");
	//puts("returning");
	return(k);
}

NODE* build_tree_array(TREE_TYPE **point_array, int N, int K, int depth){
	//This builds a k-dimensional tree from the point_array 2d array.
	
	int dim, m=0, i; //counters
	int approx_limit = 10000; //if more points than this, find approximate median
	int approx_N = 1000; //use this many points to find approximate median
	NODE *node = NULL; //node that we're finding
	TREE_TYPE **larray, **rarray;// array to send to left and right child
	
	//printf("In 'build_tree_array', depth = %d, N = %d, K = %d\n", depth, N, K);

	//NULL_CHECK(point_array, "", "point_array", i, N);

	if(N==0){ // if no data, return
		return(NULL);
	}
	
	//NULL_CHECK(point_array, "", "point_array", i, N);
	
	//puts("Allocating memory to node");
	node = (NODE*) malloc( (int)sizeof(NODE)); // allocate the memory for the node
	//printf("node: %u, node->left: %u\n", node, node->left);
	MALLOC_CHECK(node);
	//puts("Allocation successful");
	init_node(node); //initialise the node
	//NULL_CHECK(point_array, "", "point_array", i, N);

	dim = depth % K; // find which dimension we're looking at

	//puts("Finding median of 'point_array'.");
	if( (N > 1) && (N < approx_limit) ){ // if 1<N<approx_limit, sort them and find median, this gives an exactly balanced tree
		qsort_nd2(point_array, N, dim);
		m = median_array(point_array, N, dim);
	}
	else if(N >= approx_limit){ //if N<approx_limit, find an approximate value of the median, this gives pretty well balanced trees
		TREE_TYPE **approx_array; //pointer
		//puts("N is larger than approx_limit, allocating space for approximation arrays");
		approx_array = (TREE_TYPE**) malloc( (int)sizeof(TREE_TYPE*)*approx_N ); //allocate stuff
		srand(time(NULL)); //seed the random number generator
		//puts("seeding approximate array");
		for(i=0;i<approx_N;i++){
			int j = rand() % N; //make sure that j is a random number between 0 and N-1
			approx_array[i] = point_array[j]; //populate the approximation array
		}
		//puts("sort reduced array");
		qsort_nd2(approx_array, approx_N, dim); //sort our reduced array
		//puts("find median of reduced array");
		m = median_array(approx_array, approx_N, dim); //find it's median
		//printf("m: %d\n", m);
		//puts("order the real array to be split up: lt, median, gt.");
		m = find_in_array(point_array, approx_array[m], N); //where is it in the big array?
		//printf("median in big array, m: %d\n", m);
		m = order_lt_gt(point_array, m, N, dim);	//change order of 'point_array' to {[<m], m, [>m]}
		//printf("median now at, m: %d\n", m);		//and tell us 'm'
		FREE(approx_array); //free stuff and set to NULL
		//continue as normal
	}
	//puts("Found median");
	//printf("m=%d\n", m);
	//NULL_CHECK(point_array, "", "point_array", i, N);
	node->val = point_array[m]; //store the value of the median data point in this node
	//puts("Finding left array");
	larray = point_array; //left array is from the beginning, to the value before the median
	node->left = build_tree_array(larray, m, K, depth+1);
	//puts("Finding right array");
	rarray = &(point_array[m+1]); // right array is from the value after the median to the end
	node->right = build_tree_array(rarray, (N-(m+1)), K, depth+1);
	return(node);
}

void printbits(unsigned int v, unsigned int depth){
	int i, j=10;
  	for(i = depth; i >= 0; i--) putchar('0' + ((v >> i) & 1));
  	for(i=depth-j; i<0;i++) putchar(' '); //line up everything
}

void print_tree(NODE *node, unsigned int TREE_SORT, unsigned int depth){
	//puts("printing tree structure:");
	//puts("1 is a left turn, 0 is a right turn, root starts at 0");
	printf("TREE_SORT, depth(%d): ", depth);
	printbits(TREE_SORT, depth-1);
	printf("\t{%G, %G, %G}\t%u\n", node->val[0], node->val[1], node->val[2], node);

	if(node->left != NULL){
		TREE_SORT = TREE_SORT << 1;
		TREE_SORT++;
		print_tree(node->left, TREE_SORT, depth+1);
		TREE_SORT--;
		TREE_SORT = TREE_SORT >> 1;
	}
	if(node->right != NULL){
		TREE_SORT = TREE_SORT << 1;
		print_tree(node->right, TREE_SORT, depth+1);
		TREE_SORT = TREE_SORT >> 1;
	}
	return;
}

void link_parents(NODE *node){
	//puts("Linking parent nodes");
	
	if(node->left!=NULL){
		(node->left)->parent = node;
		link_parents(node->left);
	}
	if(node->right!=NULL){
		(node->right)->parent = node;
		link_parents(node->right);
	}
	return;
}

double dist2(TREE_TYPE *a, double *b, int n_dim){
	//this decides what is called close, change this if need to do strange things
	double dist2 = 0.0; //make sure this is zero
	int i;
	for(i=0;i<n_dim;i++){
		dist2 += (a->TREE_SORT[i]-b[i])*(a->TREE_SORT[i]-b[i]); // dx^2+dy^2+dz^2=dist2
	}
	return(dist2);
}

double dist2_d(double *a, double *b, int n_dim){
	//this decides what is called close, change this if need to do strange things
	double dist2 = 0.0; //make sure this is zero
	int i;
	for(i=0;i<n_dim;i++){
		dist2 += (a[i]-b[i])*(a[i]-b[i]); // dx^2+dy^2+dz^2=dist2
	}
	return(dist2);
}

int any_gt(double* list, double test, int n){
	//check if any of 'list' is greater than 'test'
	//returns the index of the value which is most greater than 'test'
	double trial=0.0, best=0.0;
	int largest = -1;
	int i;
	for(i=0;i<n;i++){
		//printf("\t\tlargest %d, test %G, best %G, trial %G, list[i] %G\n", largest, test, best, trial,i, list[i]);
		if(list[i] > test){
			trial = list[i] - test;
		}
		if (trial > best){
			best = trial;
			largest = i;
		}
	}
	return(largest);
}

int is_in(NODE** list, NODE* test, int n){
	int i;
	for(i=0;i<n;i++){
		if(list[i] == test){
			return(i); //'test' is a member of 'list'
		}
	}
	return(-1);//'test' is not a member of 'list'
}

void nearest_nn(NODE* node ,NODE** guess, double *g_dist2, TREE_TYPE* point, int n, int depth, int n_dim){
	//'node' = root of the tree to search
	//'guess' = an array that can hold all 'n' nearest neighbours
	//'point' = the point that we want to find the neighbours of
	//'g_dist2' = an array that can hold all 'n' nearest neighbour distances
	//'n' = the number of nearest neighbours to look for
	//'depth' = how far down the tree we are (root is at depth 0)
	//'n_dim' = now many dimensions the data in the tree has
	double t_dist2, diff;
	double hyp_dist2;
	NODE *first, *second;
	int dim = depth % n_dim;
	int change, present;

	if(depth==0){ //this will only run right at the first call, for inits and checks
		while(node->parent!=NULL){
			node = node->parent; //always start from the root node
		}
		recompute_nn_dist(guess, g_dist2, point->TREE_SORT, n, n_dim); //make sure any previous guesses are correct
	}
	if(node->val->ident != point->ident){ //if the node we're testing isn't the point we're looking around
		present = is_in(guess, node, n);
		if((present == -1)){ //if the node we're testing is not in the guess list already
			t_dist2 = dist2(node->val, point->TREE_SORT, n_dim); //find it's distance
			//puts("\tsending to 'any_gt'");
			change = any_gt(g_dist2, t_dist2, n); // find if any current guesses are farther away than 'node'
			//printf("\tchange = %d\n", change);
			if(change != -1){ //if we found a closer one
				g_dist2[change] = t_dist2; //replace it
				guess[change] = node;
			}
		}
	}
	hyp_dist2 = dist2_d(&(node->val->TREE_SORT[dim]), &(point->TREE_SORT[dim]), 1); //find distance to hyper-plane
	if(point->TREE_SORT[dim] < node->val->TREE_SORT[dim]){
		first = node->left;
		second = node->right;
	}
	else{
		first = node->right;
		second = node->left;
	}
	
	if(first!=NULL){ //if can go left, do so
		nearest_nn(first, guess, g_dist2, point, n, depth+1, n_dim);
	}
	change = any_gt(g_dist2, hyp_dist2, n); // is the hyper-plane closer than any of the guesses?
	if((second != NULL) && (change!=-1)){ //if can go right, and there could be a closer node there, do so
		nearest_nn(second, guess, g_dist2, point, n, depth+1, n_dim);
	}
	//puts("\tReturning.");
	return;
}

void nearest_n(NODE* node, NODE** guess, TREE_TYPE* point, double *g_dist2, int depth, int n_dim){
	double t_dist2;
	double hyp_dist2;
	int dim = depth % n_dim;
	NODE *first, *second;
	
	//printf("\t'nearest_n', node: %u, guess %u, g_dist2: %G, depth: %d\n", node, guess, *g_dist2, depth);
	if (point->ident != node->val->ident){ //if the node we're looking at isn't the point we're looking around
		t_dist2 = dist2(node->val, point->TREE_SORT, n_dim);
		if(t_dist2 < *g_dist2){//if this point is closer
			*g_dist2 = t_dist2;//it becomes new guess
			*guess = node;
		}
	}
	hyp_dist2 = dist2_d(&(node->val->TREE_SORT[dim]), &(point->TREE_SORT[dim]), 1); //what is distance to hyper plane?
	if(point->TREE_SORT[dim] < node->val->TREE_SORT[dim]){
		first = node->left;
		second = node->right;
	}
	else{
		first = node->right;
		second = node->left;
	}
	if(first != NULL){// if there is a left node go there
		nearest_n(first, guess, point, g_dist2, depth+1, n_dim);
	}
	if( (second != NULL) && (hyp_dist2 < *g_dist2) ){ //if there is a right node, only go there if there might be closer points
		nearest_n(second, guess, point, g_dist2, depth+1, n_dim);
	}
	
	return;
}

void create_tree(NODE** root, TREE_TYPE** data, int N, int K){
	int i;
	//TREE_TYPE** storage=NULL; //storage array for 'build_tree_array' function
	//root = pointer to a pointer to a NODE structure
	//data = array in format { pointer_to_TREE_TYPE0, ptr2TREE_TYPE1, ptr2TREE_TYPE2,... ptr2TREE_TYPEN }
	// N   = Number of data points
	// K   = Number of dimensions (in above example have three, x, y, and z)
	// THIS MUST BE RUN BEFORE EVERYTHING ELSE!
	printf("In 'create_tree': root = %u, data = %u, N = %d, K = %d\n", root, data, N, K);
	puts("passing to 'build_tree_array'.");
	*root = build_tree_array(data, N, K, 0);
	
	//link up parents and children in tree
	puts("Linking up parents.");
	link_parents(*root);
	
	//exit
	puts("Returning from 'create_tree'.");
	return;
}

void destroy_tree(NODE* node){
	//have to assume that each 'val' pointer in a tree is unique
	//this WILL NOT free the 'data' array associated with the tree,

	if(node->left != NULL){
		//printf("at node: {%G, %G}\n", node->val[0], node->val[1]);
		//puts("this node has a left node, kill that first, then un-link");
		destroy_tree(node->left);
		//free(node->left);
		node->left=NULL;
	}
	if(node->right != NULL){
		//printf("at node: {%G, %G}\n", node->val[0], node->val[1]);
		//puts("this node has a right node, kill that first, then un-link");
		destroy_tree(node->right);
		//free(node->right);
		node->right=NULL;
	}
	if(node->left==NULL && node->right==NULL){ //if we are a child node
		//printf("at node: {%G, %G}, location: %X\n", node->val[0], node->val[1], node->val);
		//puts("I am a leaf node, kill myself\nun-link my value...");
		//free(node->val); //free my value
		node->val = NULL; //set the pointer to NULL
		//puts("un-link from parent...");
		node->parent=NULL;
		//puts("free myself...");
		free(node);
		node = NULL;
		//puts("all done, now returning.\n");
		
		return;
	}
}

void tree_insert(NODE* node, TREE_TYPE* point, int K, int depth){
	//inserts a node into the tree that 'node' is a part of, best to start with the root node
	int dim = depth % K;
	
	if(point->TREE_SORT[dim] < node->val->TREE_SORT[dim]){ //if should go left
		if(node->left == NULL){ //and there isn't a left node
			node->left = (NODE*) malloc( sizeof(NODE) ); //create one
			init_node(node->left); //initialise it
			node->left->val = point; //give it a value
			node->left->parent = node; //link to parent
		}
		else{ //otherwise go to that node
			tree_insert(node->left, point, K, depth+1);
		}
	}
	if(point->TREE_SORT[dim] >= node->val->TREE_SORT[dim]){ //if should go right
		if(node->right == NULL){ //and there isn't a right node
			node->right = (NODE*) malloc( sizeof(NODE) ); //create one
			init_node(node->right); //initialise it
			node->right->val = point; //give it a value
			node->right->parent = node; //link to parent
		}
		else{ //otherwise go to that node
			tree_insert(node->right, point, K, depth+1);
		}
	}
	return;
}

void tree_insert_array(NODE* root, TREE_TYPE** data, int N, int K){
	//routine for inserting an array of data, will always go to the root node first
	int i;
	
	while(root->parent!=NULL){
		root->parent = root;
	}
	for(i=0;i<N;i++){
		tree_insert(root, data[i], K, 0);
	}
	return;
}

int find_in_array(TREE_TYPE **point_array, TREE_TYPE *point, int N){
	int i;
	for(i=0;i<N;i++){
		//printf("check id: %d against id: %d\n", point->ident, point_array[i]->ident);
		if(point == point_array[i]){ //if these are the same is in array (works on pointers)
			return(i);
		}
	}
	return(-1); //if the object is not in the array
}