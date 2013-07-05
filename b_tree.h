#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#ifndef NULL_CHECK
	#define NULL_CHECK(array, str, name, i, N) \
		for(i=0;i<N;i++){\
			if(array[i]==NULL)\
			fprintf(stderr, "%sWarning: %s[%d] is NULL\n", str, name, i);\
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

//creates the structure 'node', represents a node in a tree
struct node{
	double		*val; //value of the node
	struct node	*left, *right; //pointer to the children
	struct node *parent; //pointer to the parent node
};

typedef struct node NODE;

//*****function declarations*****//
void init_node(NODE* node);
void init_nn_list(NODE **nn_list, double *dist2_list, int n);
void recompute_nn_dist(NODE **nn_list, double *nn_dist2, double *point, int n, int n_dim);
void qsort_nd2(double** data, double** storage, int N, int n_dim, int dim);
int median_array(double **array, int N, int dim);
int order_lt_gt(double** array, double** store, double* val, int N, int K, int dim);
NODE* build_tree_array(double **point_array, double **storage, int N, int K, int depth);
void printbits(unsigned int v, unsigned int depth);
void print_tree(NODE *node, unsigned int pos, unsigned int depth);
void link_parents(NODE *node);
double dist2(double *a, double *b, int n_dim);
int any_gt(double* list, double test, int n);
int is_in(NODE** list, NODE* test, int n);
void nearest_nn(NODE* node, NODE** guess, double *g_dist2, double* point, int n, int depth, int n_dim);
void nearest_n(NODE* node, NODE** guess, double* point, double *g_dist2, int depth, int n_dim);
void create_tree(NODE** root, double** data, int N, int K);
void destroy_tree(NODE* node);
void tree_insert(NODE* node, double* point, int K, int depth);
void tree_insert(NODE* node, double* point, int K, int depth);
void tree_insert_array(NODE* root, double** data, int N, int K);

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
		nn_dist2[i] = dist2(nn_list[i]->val, point, n_dim);
		//printf("new_dist: %G\tdiff: %G\n", nn_dist2[i], temp - nn_dist2[i]);
	}
}

void qsort_nd2(double** data, double** storage, int N, int n_dim, int dim){
	//need to pass a data array and a storage array. will junk the contents of both arrays
	//so must ensure that data array is sorted correctly
	int i, j, k, pivot;
	
	if(N<=1) return;
	//printf("In qsort_nd2\narguments: N=%d, n_dim=%d, dim=%d\n", N, n_dim, dim);
	//puts("data array:");
	for(i=0;i<N;i++){
	//	printf("data[%d][%d] = %G\t%u\n", i, dim, data[i][dim], data[i]);
		storage[i] = data[i]; //store the data
	}
	
	pivot = N/2; //choose middle number as the pivot
	
	//printf("pivot: data[%d][%d] = %G \t %u\n", pivot, dim, data[pivot][dim], data[pivot]);
	
	j=0;
	k=(N-1) ;
	for(i=0;i<N;i++){
		//printf("loop count %d, looking at %G\n", i, storage[i][dim]);
		if(i==pivot)//don't sort the pivot yet
			continue;
		if(storage[i][dim]<storage[pivot][dim])//if smaller than pivot
		{
			//printf("put %G left\n", storage[i][dim]);
			data[j] = storage[i]; //put value at beginning of data array
			j++;
		}
		if(storage[i][dim]>=storage[pivot][dim])//if greater than or equal to pivot
		{	
			//printf("put %G right\n", storage[i][dim]);
			data[k] = storage[i]; //put value at end of data array
			k--;
		}
	}
	
	data[j] = storage[pivot]; //pivot is now in correct position in data array
	k++; //add one to k, want to start at beginning of 'greater' array
	//printf("j:%d, k:%d\n", j ,k);
	
	//puts("storage array:");
	//for(i=0;i<N;i++){
	//	printf("storage[%d][%d] = %G\t%u\n", i, dim, storage[i][dim], storage[i]);
	//}
	
	//puts("sort 'left' array.");
	if(j>1)//if more than one value in 'less' array
		qsort_nd2(data, storage, j, n_dim, dim); //sort the values and put in 'data' array
	
	//puts("sort 'right' array");
	//printf("%d, %d\n", k, N-1);
	if(k<(N-1))// if more than one value is in 'greater' array
		qsort_nd2( (data+k), (storage+k), N-k, n_dim, dim);//sort the values and put into 'data' array
		
	//now the arrays should be sorted
	//puts("SORTED array:");
	//for(i=0;i<N;i++){
	//	printf("data[%d][%d] = %G\t%u\n", i, dim, data[i][dim], data[i]);
	//}
	//puts("return from 'qsort_nd2'");
	
	return;
}

int median_array(double **array, int N, int dim){
	int i;
	//returns the position of the median value
	//puts("In 'median_array':");
	for(i=N/2; (array[i][dim]-array[i-1][dim])==0.0; i--){
		if(i==1){
			i=0;
			break;
		}
	}
	//printf("Median is array[%d]: {%G %G}\t N is %d\n", i, array[i][dim], array[i][dim], N);
	return(i);

}

int order_lt_gt(double** array, double** store, double* val, int N, int K, int dim){
	//put 'array[dim]' into order of: stuff less than 'val[dim]', 'val[dim]', stuf greater than 'val[dim]'
	//return position of 'val[dim]'
	int i, j, k;
	//puts("copying array into storage");
	
	//NULL_CHECK(array, "\t", "array", i, N);
	
	for(i=0;i<N;i++){
		store[i] = array[i]; //copy array into storage
	}
	//NULL_CHECK(store, "\t", "store", i, N);
	
	j=0; //start at beginning of array
	k=N-1; //start at end of array
	//puts("starting ordering loop");
	for(i=0;i<N;i++){
		//if(array[i]==NULL) printf("\t NULL");
		//printf("\tin loop, store[%d] = %u, val = %u\n", i, store[i], val);
		if(store[i] == val){ //if these point to the same place
			//puts("\tignoring value");
			continue; //ignore it
		}
		//puts("\t1");
		if(store[i][dim] < val[dim]){ //if data is less than value
			//puts("\tmoving to start");
			array[j] = store[i]; //put in left of array
			j++;
		}
		//puts("\t2");
		if(store[i][dim] >= val[dim]){ //if data is greater than or equal to value
			//puts("\tmoving to end");
			array[k] = store[i]; //put in right of array
			k--;
		}
	}
	//puts("finished ordering loop");
	array[j] = val; // put 'val[dim]' between the 'less than' and 'greater than' arrays
	//puts("returning");
	return(j);
}

NODE* build_tree_array(double **point_array, double **storage, int N, int K, int depth){
	//This builds a k-dimensional tree from the point_array 2d array.
	
	int dim, m=0, i; //counters
	int approx_limit = 10000; //if more points than this, find approximate median
	int approx_N = 1000; //use this many points to find approximate median
	NODE *node = NULL; //node that we're finding
	double **larray, **lstore; // array to send to left child
	double **rarray, **rstore;// array to send to right child

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
		qsort_nd2(point_array, storage, N, K, dim);
		m = median_array(point_array, N, dim);
	}
	else if(N >= approx_limit){ //if N<approx_limit, find an approximate value of the median, this gives pretty well balanced trees
		double **approx_array, **approx_store; //pointers
		//puts("N is larger than approx_limit, allocating space for approximation arrays");
		approx_array = (double**) malloc( (int)sizeof(double*)*approx_N ); //allocate stuff
		approx_store = (double**) malloc( (int)sizeof(double*)*approx_N );
		srand(time(NULL)); //seed the random number generator
		//puts("seeding approximate array");
		for(i=0;i<approx_N;i++){
			int j = rand() % N; //make sure that j is a random number between 0 and N-1
			approx_array[i] = point_array[j]; //populate the approximation array
		}
		//puts("sort reduced array");
		qsort_nd2(approx_array, approx_store, approx_N, K, dim); //sort our reduced array
		//puts("find median of reduced array");
		m = median_array(approx_array, approx_N, dim); //find it's median
		//puts("order the real array to be split up: lt, median, gt.");
		m = order_lt_gt(point_array, storage, approx_array[m], N, K, dim);	//change order of 'point_array' to {[<m], m, [>m]}
																			//and tell us 'm'
		free(approx_array); //free stuff and set to NULL
		approx_array = NULL;
		free(approx_store);
		approx_store = NULL;
		//continue as normal
	}
	//puts("Found median");
	node->val = point_array[m]; //store the value of the median data point in this node
	//printf("m=%d\n", m);
	//puts("Finding left array");
	larray = point_array; //left array is from the beginning, to the value before the median
	lstore = storage; //similarly with storage array
	node->left = build_tree_array(larray, lstore, m, K, depth+1); // recurse

	//puts("Finding right array");
	rarray = &(point_array[m+1]); // right array is from the value after the median to the end
	rstore = &(storage[m+1]); //similarly with storage
	node->right = build_tree_array(rarray, rstore, (N-(m+1)), K, depth+1); // recurse

	return(node);
}

void printbits(unsigned int v, unsigned int depth){
	int i, j=10;
  	for(i = depth; i >= 0; i--) putchar('0' + ((v >> i) & 1));
  	for(i=depth-j; i<0;i++) putchar(' '); //line up everything
}

void print_tree(NODE *node, unsigned int pos, unsigned int depth){
	//puts("printing tree structure:");
	//puts("1 is a left turn, 0 is a right turn, root starts at 0");
	printf("pos, depth(%d): ", depth);
	printbits(pos, depth-1);
	printf("\t{%G, %G, %G}\t%u\n", node->val[0], node->val[1], node->val[2], node);

	if(node->left != NULL){
		pos = pos << 1;
		pos++;
		print_tree(node->left, pos, depth+1);
		pos--;
		pos = pos >> 1;
	}
	if(node->right != NULL){
		pos = pos << 1;
		print_tree(node->right, pos, depth+1);
		pos = pos >> 1;
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

double dist2(double *a, double *b, int n_dim){
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

void nearest_nn(NODE* node ,NODE** guess, double *g_dist2, double* point, int n, int depth, int n_dim){
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
	//after this = debugging
	//int i;
	//printf("\tIn 'nearest_nn': node: {%G, %G, %G}, location %x, depth: %d\n", node->val[0], node->val[1], node->val[2], node, depth);
	//puts("\tguesses:");
	//for(i=0;i<n;i++){
	//	if(guess[i]!=NULL){
	//		printf("\t%x: {%G, %G, %G}\n", guess[i], guess[i]->val[0], guess[i]->val[1], guess[i]->val[2]);
	//	}
	//	else printf("\tNULL\n");
	//}
	//puts("\tdistances:");
	//for(i=0;i<n;i++){
	//	printf("\t%G\n", g_dist2[i]);
	//}
	if(depth==0){ //this will only run right at the first call, for inits and checks
		while(node->parent!=NULL){
			node = node->parent; //always start from the root node
		}
		recompute_nn_dist(guess, g_dist2, point, n, n_dim); //make sure any previous guesses are correct
	}
	present = is_in(guess, node, n);
	if(present == -1){ //if the node we're testing is not in the guess list already
		t_dist2 = dist2(node->val, point, n_dim);
		//puts("\tsending to 'any_gt'");
		change = any_gt(g_dist2, t_dist2, n); // find if any current guesses are farther away than 'node'
		//printf("\tchange = %d\n", change);
		if(change != -1){ //if we found a closer one
			g_dist2[change] = t_dist2; //replace it
			guess[change] = node;
		}
	}
	hyp_dist2 = dist2(&(node->val[dim]), &(point[dim]), 1); //find distance to hyper-plane
	if(point[dim] < node->val[dim]){
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

void nearest_n(NODE* node, NODE** guess, double* point, double *g_dist2, int depth, int n_dim){
	double t_dist2;
	double hyp_dist2;
	int dim = depth % n_dim;
	NODE *first, *second;
	
	//printf("\t'nearest_n', node: %u, guess %u, g_dist2: %G, depth: %d\n", node, guess, *g_dist2, depth);
	t_dist2 = dist2(node->val, point, n_dim);
	if(t_dist2 < *g_dist2){//if this point is closer
		*g_dist2 = t_dist2;//it becomes new guess
		*guess = node;
	}
	hyp_dist2 = dist2(&(node->val[dim]), &(point[dim]), 1); //what is distance to hyper plane?
	if(point[dim] < node->val[dim]){
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

void create_tree(NODE** root, double** data, int N, int K){
	int i;
	double** storage=NULL; //storage array for 'build_tree_array' function
	//root = pointer to a pointer to a NODE structure
	//data = array in format { {x0,y0,z0}, {n1,y1,z1}, ... {xn, yn, zn} }
	// N   = Number of data points
	// K   = Number of dimensions (in above example have three, x, y, and z)
	printf("In 'create_tree': root = %u, data = %u, N = %d, K = %d\n", root, data, N, K);
	//NULL_CHECK(data, "\t", "data", i, N);
	//allocate memory, only going to store pointers to doubles in it
	puts("Allocating 'storage' array.");
	storage = (double**) malloc( (int)sizeof(double*)*N );
	MALLOC_CHECK(storage);
	//make the tree
	//'root' is points to where we can store the pointer to the root of the array
	//'build_tree_array()' returns the location of the root node of the tree it just built.
	puts("passing to 'build_tree_array'.");
	*root = build_tree_array(data, storage, N, K, 0);
	
	//free the memory
	puts("Freeing 'storage' array.");
	FREE(storage);
	
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

void tree_insert(NODE* node, double* point, int K, int depth){
	//inserts a node into the tree that 'node' is a part of, best to start with the root node
	int dim = depth % K;
	
	if(point[dim] < node->val[dim]){ //if should go left
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
	if(point[dim] >= node->val[dim]){ //if should go right
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

void tree_insert_array(NODE* root, double** data, int N, int K){
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