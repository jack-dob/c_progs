#ifndef SSHEAD_INCLUDED
	#include "sshead.h"
	#define SSHEAD_INCLUDED
#endif
#ifndef VECTOR_INCLUDED
	#include "vector.h"
	#define VECTOR_INCLUDED
#endif
#ifndef B_TREE_INCLUDED
	#include "b_tree.h"
	#define B_TREE_INCLUDED
#endif

typedef struct particle{
	double 	time;
	double	mass;
	double	radius;
	double	hill_rad;
	double	pos[N_DIM];
	double	vel[N_DIM];
	double	spin[N_DIM];
	double	ecc;
	double	inc;
	double	sma;
	int		color;
	int		org_idx;
} PARTICLE;

typedef struct dust_ellipses{
	/*	this is separate from the particle structure
		contains all necessary data to plot orbit of 
		a particle and it's associated dust distribution
		ra, rp give size of orbit
		psi and angle_vec give orientation of orbit
	*/
	int org_idx; /*so we can cross-compare the dust and particles if needed, although they should be in same order*/
	int ejected_flag; /*flags whether the particle has been ejected, also used for corruption or if could not allocate memory*/
	double dr; /*dr at apohelion, this should be the biggest*/
	double max_mass; /*the mass of the largest object in the cascade, is mass of parent particle*/
	double ecc;/*orbit of dust, is instantaneous eccentricity of particle*/
	double sma;/*orbit of dust, is instantaneous semi-major axis of particle*/
	double inc;/*orbit of dust, is instantaneous inclination of particle, not sure if will use this*/
	double r_a; /*radial distance of the particle's apogee, will be in +ve x-direction*/
	double r_p; /*radial distance of the particle's perigee, will be in -ve x-direction*/
	double psi; /*value of the angular mis-match between x-axis and line connecting apogee and perigee of orbit*/
	double angle_vec[3]; /*values of the alpha, beta, gamma angles needed to rotate data to ecliptic*/
	double **inv_transform; /*Stores the transform from plane of ecliptic to original orbit*/
	double t_coll; /*collision timescale, for use in dust mass approximation later*/
}DUST_ELLIPSES;

void vel_disp_map(struct parameters *PAR, SSHEAD *head, SSDATA *data, double ****dust_cube);
double t_coll(DUST_ELLIPSES *dust_ellipses, SSHEAD *head, SSDATA *data, int i, NODE *root, int N);
void rot_to_ecliptic(SSDATA *data, double **angle_vec, double ***inv_transform);
void dust_params(SSHEAD *head, SSDATA *data, DUST_ELLIPSES *dust_ellipses);
double sma(SSDATA particle);
double hill_rad(SSDATA particle);
double ecc(SSDATA particle);
double r_a(SSDATA particle);
double r_p(SSDATA particle);
int grab_planets(SSHEAD *head, SSDATA *data, SSDATA **planet_list, int *N_planets);
int delete_flagged(SSHEAD *head, SSDATA **data);
int orbital_params(SSHEAD *head, SSDATA* data, PARTICLE *particle);
void nn_density_SS(NODE **n_near, double *n_dist2, SSDATA *data, double sigmas, double* location, int n_dim, int n);
void nn_num_density_SS(NODE **n_near, double *n_dist2, SSDATA *data, double sigmas, double *location, int n_dim, int n);
void nn_vel_disp_SS(NODE **n_near, double *n_dist2, SSDATA *data, double *ref_vel, double *location, int n_dim, int n);
void nn_vel_disp_point(NODE **n_near, double *n_dist2, SSDATA *data, double *location, int n_dim, int n);
double max(double *array, int n);
int min(double *array, int n);
void find_vel(double r2, double *ref_point, double *re_vel);

void vel_disp_map(struct parameters *PAR, SSHEAD *head, SSDATA *data, double ****dust_cube){
	int e, f, g, i, j;
	int n=128, N_dim=3, N;
	NODE *root;
	NODE **n_near;
	double **tree_data;
	double *n_dist2;
	double r2;
	double vel_disp;
	double ref_vel[3]; //need to make this a vector
	double ref_point[3];
	double z, y, x;
	
	root = (NODE*) malloc( (int)sizeof(NODE) ); //have our root node
	//now create data to go into the tree
	//only include non-ejected particles
	printf("number of particles to analyse = %i\n", head->n_data);
	N=0;
	for(i=0;i<head->n_data;i++){
		if(sma(data[i]) >= 0.0){
			N++;
		}
	}
	/*Create tree outside loop*/
	printf("Number of particles that are not ejected: %d\n", N);
	tree_data = (double**) malloc( N*sizeof(double*) );
	j=0;
	for(i=0;i<head->n_data;i++){
		if(sma(data[i]) >= 0.0){
			tree_data[j] = (double*) malloc( (N_dim+1)*sizeof(double) ); //don't need to do this, space already there
			tree_data[j][0] = data[i].pos[0];
			tree_data[j][1] = data[i].pos[1];
			tree_data[j][2] = data[i].pos[2];
			tree_data[j][3] = (double)i;
			j++;
		}
	}
	create_tree(&root, tree_data, N, N_dim);
	//want to store the velocity dispersion of each grid point, will have to be relative to kepler velocity (circular orbit)
	//put all this stuff outside this function!! that way it is only done once like it's supposed to.
	printf("starting dust_cube scaling...\n");
		
	n_near = (NODE**) malloc( n*sizeof(NODE*) );
	n_dist2 = (double*) malloc( n*sizeof(double) );
	//puts("1");
	init_nn_list(n_near, n_dist2, n);
	
	for(e=0;e<PAR->cube_res;e++){
		x = (PAR->cube_scale/PAR->cube_res)*(PAR->cube_res/2.0 - (double)e);
		//printf("making scalar: x %g            \r",x);
		for(f=0;f<PAR->cube_res;f++){
			y = (PAR->cube_scale/PAR->cube_res)*(PAR->cube_res/2.0 - (double)f);
			for(g=0;g<PAR->cube_z_res;g++){
				z = (PAR->cube_z_scale/PAR->cube_z_res)*(PAR->cube_z_res/2.0 - (double)g);
				r2 = x*x+y*y+z*z;
				ref_point[0] = x;
				ref_point[1] = y;
				ref_point[2] = z;
				find_vel(sqrt(r2), ref_point, ref_vel);
				nearest_nn(root, n_near, n_dist2, ref_point, n, 0, N_dim); //last number is number of dimensions to find nearest neighbour in
				nn_vel_disp_SS(n_near, n_dist2, data, ref_vel, &vel_disp, N_dim, n);
				//nn_vel_disp_point(n_near, n_dist2, data, &vel_disp, N_dim, n);
				(*dust_cube)[e][f][g] = vel_disp;
				printf("making scalar: x %g, y %g, z %g            \r",x, y, z);
			}
		}
	}
	
	//free everything
	destroy_tree(root);
	free(tree_data);

}

double t_coll(DUST_ELLIPSES *dust_ellipses, SSHEAD *head, SSDATA *data, int i, NODE *root, int N){
	//annoyingly this doesn't do what I'd like it too. I think I have to include an orbit crossing
	//term in my nearest neighbour algorithm
	NODE **n_near;
	double *n_dist2;
	double density, num_density, vel_disp;
	double av_m_coll, R1, R2, p_0, v2_esc, x_section, t_coll, N_coll, t_orb, r_orb;
	double l_coll, h_rad;
	double *search_data;
	int n = 128; //more neighbours gives better results but can be very slow
	int j=0, N_dim = 1, k, l, condition=1;
	int e=0, f=0, g=0;
	
	search_data = (double*) malloc(2*sizeof(double) );
	search_data[0] = sma(data[i]);
	if(search_data[0] <= 0.0){ // if particle is ejected, it does not count towards dust.
		free(search_data);
		return(0.0);
	}
	search_data[1] = ecc(data[i]);
	//puts("in 't_coll'");
	
	n_near = (NODE**) malloc( n*sizeof(NODE*) );
	n_dist2 = (double*) malloc( n*sizeof(double) );
	//puts("1");
	init_nn_list(n_near, n_dist2, n);
	//puts("2");
	
	nearest_nn(root, n_near, n_dist2, search_data, n, 0, N_dim); //last number is number of dimensions to find nearest neighbour in
	//puts("3");
	k=0;
	l = min(n_dist2, n); //find the closest 'neighbour' to our particle. It is the actual particle not a neighbour.
	for(j=0;j<n;j++){ //check that each neighbour crosses the orbit of our particle
		int p_num = (int)n_near[j]->val[N_dim];
		if(j==l){ // don't count the actual particle.
			continue;
		}
		if( sma(data[p_num]) < 0.0 ){
			n_dist2[j] = DBL_MAX; //make lose all weighting if neighbour is ejected
		}
		//condition = condition && ( (r_a(data[i]) >= r_a(data[p_num])) && (r_p(data[i]) <= r_p(data[p_num])) );
		//condition = condition && ( (r_p(data[i]) <= r_a(data[p_num])) && (r_p(data[i]) >= r_a(data[p_num])) );
		//condition = condition && ( (r_p(data[i]) >= r_p(data[p_num])) && (r_a(data[i]) <= r_a(data[p_num])) );
		//condition = condition && ( (r_a(data[i]) >= r_p(data[p_num])) && (r_p(data[i]) <= r_a(data[p_num])) );
		if(cross_orbit(dust_ellipses[i], dust_ellipses[p_num])){
		//if(condition){
			k++;
			continue;
		}
		else{
			n_dist2[j] = DBL_MAX; //cannot put at infinity, this will have to do. Will lose all weighting
		}
	}
	

	
	//puts("4");
	nn_density_SS(n_near, n_dist2, data, 5.0, &density, N_dim, n);
	nn_num_density_SS(n_near, n_dist2, data, 5.0, &num_density, N_dim, n);
	nn_vel_disp_SS(n_near, n_dist2, data, data[i].vel, &vel_disp, N_dim, n);
	
	//puts("5");
	
	if(num_density==0.0){
		av_m_coll = DBL_MAX;
	}
	else if(density==0.0 && num_density==0.0){
		av_m_coll=1.0;
	}
	else{
		av_m_coll = density/num_density;
	}
	p_0 = 3.36627518E6; //2000.0 kg m^-3 in M_sol AU^-3
	R1 = cbrt(3.0*data[i].mass/(4.0*p_0*M_PI)); //in AU
	R2 = cbrt(3.0*av_m_coll/(4.0*p_0*M_PI)); //in AU
	v2_esc = 2.0*1.0*(data[i].mass + av_m_coll)/(R1+R2); //G=1 in system units, escape velocity for target + collider
	x_section = M_PI*(R1+R2)*(R1+R2)*(1.0+v2_esc/(vel_disp*vel_disp));
	
	if(num_density==0){
		t_coll = DBL_MAX;
		l_coll = DBL_MAX;
	}
	else{
		t_coll = 1.0/(sqrt(dot_prod(data[i].vel, data[i].vel))*num_density*x_section); //gives time for 1 collision
		l_coll = 1.0/(x_section*num_density); //gives average distance between collisions
	}
	h_rad = hill_rad(data[i]); //l_coll/h_rad is how likely something is to collide with an object in it's hill radius
	N_coll = data[i].mass/av_m_coll; //number of collisions = target mass / average collision mass
	r_orb = sqrt(dot_prod(data[i].pos, data[i].pos));
	t_orb = cbrt(4.0*M_PI*(r_orb*r_orb*r_orb));
	
	
	
	//printf("\n\tdensity: %G, num_density: %G, vel_disp: %G, av_m_coll: %G, R1: %G, R2: %G, p_0: %G, v2_esc: %G, x_section: %G, t_coll: %G, N_coll: %G",
	//	density, num_density, vel_disp, av_m_coll, R1, R2, p_0, v2_esc, x_section, t_coll, N_coll);
	free(n_near);
	free(n_dist2);
	free(search_data);
	return(t_orb/(N_coll*t_coll)); //works also
	//return(t_orb/N_coll) //works!, bright ring co-orbital with planet. Dark patches either side.
	// t_orb/t_coll is number of collisions per orbit, c_orb
	// N_coll/c_orb is number of orbits needed to blow up our mass, n_orb
	// n_coll*t_orb/t_coll = n_orb
	// dust_mass = max_mass/n_orb = max_mass*t_coll/(N_coll*t_orb)
	// finding number density and average collision mass is the thing I need to get right.
	//return(N_coll*t_coll/t_orb);
	//return(N_coll*t_orb/t_coll);
	//return(h_rad/l_coll);
	//return(t_coll*h_rad/(l_coll*t_orb*N_coll)); //kinda works
	//return(1.0/t_coll);
	//return(r_a(data[i]) - r_p(data[i]));
	//return(k);
}

void rot_to_ecliptic(SSDATA *data, double **angle_vec, double ***inv_transform){
	double alpha, beta, gamma;
	double x2y2 = sqrt(data->pos[0]*data->pos[0] + data->pos[1]*data->pos[1]);
	double mag_pos = sqrt(dot_prod(data->pos, data->pos));
	double pos3[3] = {data->pos[0]*mag_pos/x2y2, data->pos[1]*mag_pos/x2y2, 0.0};
	double *vel2;
	double **transform;
	double angles[2];
	int axis[2];
	double *ptr;
	int i, j;

	/*
	printf("INFO: data.pos: [0] = %G, [1] = %G, [2] = %G\n", data->pos[0], data->pos[1], data->pos[2]);
	printf("INFO: data.vel: [0] = %G, [1] = %G, [2] = %G\n", data->vel[0], data->vel[1], data->vel[2]);
	*/
	alpha = atan(data->pos[0]/data->pos[1]);
	beta = atan(-data->pos[2]/x2y2);
	if(DEBUG) printf("DEBUG: alpha and beta found...\n");
	
	transform = (double**) malloc( (int)3*sizeof(double*) );
	if(transform==NULL){
		printf("Error: rot_to_ecliptic could not allocate memory for transform\n");
		exit(1);
	}
	for(i=0;i<3;i++){
		transform[i] = (double*) malloc( (int)3*sizeof(double) );
		if(transform[i]==NULL){
			printf("Error: rot_to_ecliptic could not allocate memory for transform[%i]\n", i);
			exit(1);
		}
		for(j=0;j<3;j++){
			if(i!=j) transform[i][j] = 0.0;
			else transform[i][j] = 1.0;
		}
	}
	if(DEBUG) printf("DEBUG: transform array initialised...\n");
	if(DEBUG){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("DEBUG: transform[%i][%i] = %G\n", i, j, transform[i][j]);
			}
		}
	}
	
	vel2 = (double*) malloc( (int)3*sizeof(double) );
	if(vel2==NULL){
		printf("Error: rot_to_ecliptic could not allocate memory for vel2\n");
		exit(1);
	}
	for(i=0;i<3;i++){
		vel2[i] = 0.0;
	}
	if(DEBUG) printf("DEBUG: vel2 array initialised\n");
	/*printf("INFO: vel2: [0] = %G, [1] = %G, [2] = %G\n", vel2[0], vel2[1], vel2[2]);*/
	angles[0] = alpha;
	angles[1] = beta;
	axis[0] = 2;
	axis[1] = 0;
	if(DEBUG) printf("DEBUG: finding half-way rotation transformation, addr: %i...\n", &transform);
	rotate(&transform, angles, axis, 2);
	
	if(DEBUG){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("DEBUG: transform[%i][%i] = %G\n", i, j, transform[i][j]);
			}
		}
	}
	
	if(DEBUG) printf("DEBUG: applying half-way transformation to velocity, storing in vel2 %i...\n", &vel2);
	if(DEBUG) printf("DEBUG: vel2[0] = %lf, vel2[1] = %lf, vel2[2] = %lf\n", vel2[0], vel2[1], vel2[2]);
	col_vec_multi(transform, data->vel, &(vel2), 3, 3);
	/*printf("INFO: vel2: [0] = %G, [1] = %G, [2] = %G\n", vel2[0], vel2[1], vel2[2]);*/
	gamma = atan(vel2[2]/vel2[0]);
	
	angles[0] = gamma;
	angles[1] = -alpha;
	axis[0] = 1;
	axis[1] = 2;
	if(DEBUG) printf("DEBUG: angle gamma calculated, finding complete transfrom...\n");
	rotate(&transform, angles, axis, 2);
	
	for(i=0;i<N_DIM;i++){
		data->pos[i] = pos3[i];
	}
	
	if(DEBUG) printf("DEBUG: transformed position stored, applying complete transform to velocity, storing in vel2..\n");
	col_vec_multi(transform, data->vel, &(vel2), 3, 3);
	
	for(i=0;i<N_DIM;i++){
		data->vel[i] = vel2[i];
	}
		
	(*angle_vec)[0] = alpha;
	(*angle_vec)[1] = beta;
	(*angle_vec)[2] = gamma;
	if(DEBUG) printf("DEBUG: transformed velocity stored, transform angles stored (retain orientation of original orbit)\n");
	
	/*printf("BLAH: %G\n", transform[1][2]);*/
	
	invert(transform, inv_transform, 3);
	
	if(DEBUG) printf("DEBUG: Inverted transformation...\n");
	if(0){
		/*FOR DEBUGGING*/
		double **check;
		check = (double**) malloc( (int)(3*sizeof(double*)) );
		for(i=0;i<3;i++){
			check[i] = (double*) malloc( (int)(3*sizeof(double)) );
		}
		printf("\nChecking Transforms...\n");
		for(i=0;i<3;i++){
			printf("|%012G  %012G  %012G|\t\t|%012G  %012G  %012G|\n",transform[i][0],transform[i][1],transform[i][2],
			(*inv_transform)[i][0],(*inv_transform)[i][1],(*inv_transform)[i][2]);
		}
		mat_multi(transform, (*inv_transform), &check, 3,3,3,3,2);
		for(i=0;i<3;i++){
			printf("|%012G  %012G  %012G|\t|%012G  %012G  %012G|\t\t|%012G  %012G  %012G|\n",transform[i][0],transform[i][1],transform[i][2],
			(*inv_transform)[i][0],(*inv_transform)[i][1],(*inv_transform)[i][2],
			check[i][0],check[i][1],check[i][2]);
		}
		for(i=0;i<3;i++){
			free(check[i]);
		}
		free(check);
		/*exit(0);
		/*END FOR DEBUGGING*/
	}
	for(i=0;i<3;i++){
		free(transform[i]);
	}
	free(transform);
	if(DEBUG) printf("DEBUG: transform array freed\n");
	/*printf("INFO: vel2: [0] = %G, [1] = %G, [2] = %G\n", vel2[0], vel2[1], vel2[2]);*/
	free(vel2);
	/*
	printf("INFO: data.pos: [0] = %G, [1] = %G, [2] = %G\n", data->pos[0], data->pos[1], data->pos[2]);
	printf("INFO: data.vel: [0] = %G, [1] = %G, [2] = %G\n", data->vel[0], data->vel[1], data->vel[2]);
	*/
}

void dust_params(SSHEAD *head, SSDATA *data, DUST_ELLIPSES *dust_ellipses){
	int i, j, k, order, N;
	double *angle_vec;
	double **inv_transform;
	NODE* root;
	//NODE* root_cart;//for cartesian tree
	double **tree_data;
	//double **tree_data_cart; //cartesian tree
	double testing[5000];
	//DEBUGGING
	FILE* dat = fopen("dat.txt", "w");
	for(i=0;i<5000;i++){
		testing[i]=0;
	}
	order = 3;/*set dimensions of system, here 3d*/
	
	inv_transform = (double**) malloc( (int)order*sizeof(double*) );
	if(inv_transform==NULL){
		printf(MALLOC_ERROR, __FUNCTION__, "inv_transform");
		exit(1);
	}
	for(i=0;i<3/*order*/;i++){
		inv_transform[i] = (double*) malloc( (int)order*sizeof(double) );
		if(inv_transform[i]==NULL){
			char *str;
			sprintf(str, "inv_transform[%i]", i);
			printf( MALLOC_ERROR, __FUNCTION__, str );
			exit(1);
		}
	}
	angle_vec = (double*) malloc( (int)3*sizeof(double) );
	if(angle_vec==NULL){
		printf("Error: dust_params could not allocate memory to angle_vec\n");
		exit(1);
	}
	printf("number of particles to analyse = %i\n", head->n_data);
	N=0;
	for(i=0;i<head->n_data;i++){
		if(sma(data[i]) >= 0.0){
			N++;
		}
	}
	/*Create tree outside loop*/
	printf("Number of particles that are not ejected: %d\n", N);
	tree_data = (double**) malloc( N*sizeof(double*) );
	//tree_data_cart = (double**) malloc( head->n_data*sizeof(double*) );
	j=0;
	for(i=0;i<head->n_data;i++){ //loop over everything
		double semi_ma;
		//tree_data_cart[i] = (double*) malloc( 3*sizeof(double) );
		//tree_data_cart[i][0] = data[i].pos[0];
		//tree_data_cart[i][1] = data[i].pos[1];
		//tree_data_cart[i][2] = data[i].pos[2];
		//tree_data[i][3] = (double)i; //hopefully this won't be a problem
		//re-do data
		//printf("i: %d, j: %d\n", i, j);
		
		semi_ma = sma(data[i]); // to avoid computing this twice
		if(semi_ma >= 0.0){
			tree_data[j] = (double*) malloc( 2*sizeof(double) );
			
			tree_data[j][0] = semi_ma;
			// need not include ejected particles as it buggers up nearest-neighbour search
			//tree_data[j][1] = ecc(data[i]);
			tree_data[j][1] = (double)i;
			//fprintf(dat, "%d\t%G\t%G\n", j, tree_data[j][0], tree_data[j][1]);
			j++;
		}
	}
	puts("Assigned data, creating tree...");
	root = (NODE*) malloc( sizeof(NODE) );
	//root_cart = (NODE*) malloc( sizeof(NODE) );
	//create_tree(&root_cart, tree_data_cart, N, N_DIM); //sort by x, y, z
	//vel_disp_map(PAR, head, data, root_cart);
	create_tree(&root, tree_data, N, 1); //sort by sma
	for(i=0;i<(head->n_data);i++){
		double h[3], unit_h[3];
		double unit_z[3] = {(double)0.0, (double)0.0, (double)1.0};
		double h2, v2, r2;
		double theta, p, r, r_h;

		
		/*need to transform particle to ecliptic plane to be sure the psi calculation is correct*/
		/*might as well do it first, see lab-book for explanation of maths, will also comment ftn*/
		printf("analysing particle %i\r", i);
		dust_ellipses[i].ejected_flag = 0;/*set to false initially*/
		
		
		dust_ellipses[i].inv_transform = (double**) malloc((int)(order*sizeof(double*)) );
		if(dust_ellipses[i].inv_transform==NULL){
			char *str;
			sprintf(str, "dust_ellipses[%i].inv_transform", i);
			printf(MALLOC_ERROR, __FUNCTION__, str);
			/*ignore bad allocations, but continue with rest of program*/
			dust_ellipses[i].ejected_flag = 1;
		}
		for(j=0;j<order;j++){
			dust_ellipses[i].inv_transform[j] = (double*) malloc( (int)(order*sizeof(double)) );
			if(dust_ellipses[i].inv_transform[j]==NULL){
				char *str;
				sprintf(str, "dust_ellipses[%i].inv_transform[%i]", i, j);
				printf(MALLOC_ERROR, __FUNCTION__, str);
				/*ignore bad allocations, but continue with rest of program*/
				dust_ellipses[i].ejected_flag = 1;
			}
		}
		/*catch bad allocation of memory, continue to next particle*/
		if(dust_ellipses[i].ejected_flag){
			continue;
		}
		
		if(DEBUG) printf("DEBUG: Rotating to ecliptic plane...\n");
		rot_to_ecliptic(&data[i], &angle_vec, &inv_transform);
		if(DEBUG) printf("DEBUG: Rotation finished...\n");
		
		if(DEBUG) printf("DEBUG: Storing inverse transformation...\n");
		for(j=0;j<order;j++){
			for(k=0;k<order;k++){
				dust_ellipses[i].inv_transform[j][k] = inv_transform[j][k];
			}
		}
		r2 = dot_prod(data[i].pos, data[i].pos);
		v2 = dot_prod(data[i].vel, data[i].vel);
		r = sqrt(r2);
		cross_prod(data[i].pos, data[i].vel, h);
		norm_vec(h, unit_h);
		r_h = hill_rad(data[i]);
		h2 = dot_prod(h, h);
		if(DEBUG) printf("DEBUG: Constants computed\n");
		/*
		printf("INFO: (in dust_params) data.pos: [0] = %G, [1] = %G, [2] = %G\n", data[i].pos[0], data[i].pos[1], data[i].pos[2]);
		printf("INFO: (in dust_params) data.vel: [0] = %G, [1] = %G, [2] = %G\n", data[i].vel[0], data[i].vel[1], data[i].vel[2]);
		printf("INFO: (in dust_params) r = %G, v2 = %G, mass = %G\n", r, v2, data[i].mass);
		printf("INFO: (in dust_params) sma[%i] = %G\n", i ,1.0/(2.0/r - v2/(1.0+data[i].mass)));
		printf("INFO: (in dust_params) alpha = %G, beta = %g, gamma = %G\n", angle_vec[0], angle_vec[1], angle_vec[2]);
		*/
		/*store rotation angles*/
		for(j=0;j<3;j++){
			dust_ellipses[i].angle_vec[j] = angle_vec[j];
		}
				
		/*Compute semi-major axis*/
		dust_ellipses[i].sma = 1.0/(2.0/r - v2/(1.0+data[i].mass));
		if(dust_ellipses[i].sma <= 0.0){
			dust_ellipses[i].ejected_flag = 1;
			continue;
		}
		else{
			dust_ellipses[i].ejected_flag = 0;
		}
		
		/*compute eccentricity*/
		dust_ellipses[i].ecc = sqrt(1.0 - h2/( (1.0+data[i].mass)*dust_ellipses[i].sma ));
		
		/*compute inclination,
		Assume that ecliptic is z=0 plane
		acos() outputs in RADIANS, the trig functions use RADIANS!!
		if(i==0) printf("h = %E %E %E, \nunit_h = %E %E %E\n", h[0],h[1],h[2], unit_h[0],unit_h[1],unit_h[2]);*/
		dust_ellipses[i].inc = acos( dot_prod(unit_h, unit_z) );

		/*store the maximum mass*/
		dust_ellipses[i].max_mass = data[i].mass;

		/*remember particle index*/
		dust_ellipses[i].org_idx = data[i].org_idx;
		
		/*compute apogee*/
		dust_ellipses[i].r_a = dust_ellipses[i].sma*(1.0 + dust_ellipses[i].ecc);
		
		/*compute perigee*/
		dust_ellipses[i].r_p = dust_ellipses[i].sma*(1.0 - dust_ellipses[i].ecc);
		
		
		/*compute psi*/
		/*assume is in ecliptic plane as co-ord transform was performed earlier*/
		theta = atan(data[i].pos[0]/data[i].pos[2]);
		p = dust_ellipses[i].sma*(1 - dust_ellipses[i].ecc*dust_ellipses[i].ecc);
		dust_ellipses[i].psi = theta - acos( (p - r)/(r*dust_ellipses[i].ecc) );
		
		/*compute dr*/
		dust_ellipses[i].dr = sqrt((24.0*r_h*r_h*r_h)/r);
		/*printf("i at end: %i\n", i);*/
		
		/*Use tree inside t_coll*/
		dust_ellipses[i].t_coll = t_coll(dust_ellipses, head, data, i, root, N);
		
		//DEBUGGING:
		if(1){//(dot_prod(data[i].pos, data[i].pos) > 4.0) &&(dot_prod(data[i].pos, data[i].pos) < 25.0)){
			//for(j=0;j<5000;j+=2){
			//	if(sma(data[i]) < 0.0) break;
			//	if (((double)j/(double)1000 <= r_a(data[i])) && ((double)j/(double)1000 >= r_p(data[i])) ){
			//		testing[j]+=data[i].mass;
			//		testing[j+1]+=1;
			//	}
			//}
			fprintf(dat, "%d\t%G\t%G\n", i, sma(data[i]), dust_ellipses[i].t_coll );
			//exit(0);
		}
	}
	//for(i=0;i<5000;i+=2){
	//	fprintf(dat, "%G\t%G\n", (double)i/(double)500, testing[i+1]/testing[i] );
	//}
	//DEBUGGING
	fclose(dat);
	//free stuff
	destroy_tree(root);
	//destroy_tree(root_cart);
	for(i=0;i<N;i++){
		free(tree_data[i]);
		//free(tree_data_cart[i]);
	}
	free(tree_data);
	//free(tree_data_cart);
	for(i=0;i<order;i++){
		free(inv_transform[i]);
	}
	free(inv_transform);
}

double sma(SSDATA particle){
	double r2, v2, sma;
	
	r2 = dot_prod(particle.pos, particle.pos);
	v2 = dot_prod(particle.vel, particle.vel);
	
	sma = 1.0/(2.0/sqrt(r2) - v2/(1.0+particle.mass));
	return(sma);
}

double hill_rad(SSDATA particle){
	/*units in solar masses, so M_sun = 1
	*/
	double hill_rad;
	
	hill_rad = cbrt(particle.mass/3.0)*sma(particle);

	return(hill_rad);
}

double ecc(SSDATA particle){
	double h[3], h2, ecc;
	cross_prod(particle.pos, particle.vel, h);
	h2 = dot_prod(h, h);
	ecc = sqrt(1.0 - h2/( (1.0+particle.mass)*sma(particle) ));
	return(ecc);
}

double r_a(SSDATA particle){
	return((1.0+ecc(particle))*sma(particle));
}

double r_p(SSDATA particle){
	return((1.0-ecc(particle))*sma(particle));
}

int grab_planets(SSHEAD *head, SSDATA *data, SSDATA **planet_list, int *N_planets){
	/*atm detect colours != 3, later add specific planet support
	*/
	int i, j;
	
	*N_planets=0;
	for(i=0;i<head->n_data;i++){ /*loop through and count planets*/
		if(data[i].color!=3) (*N_planets)++;
	}
	
	*planet_list = (SSDATA*) malloc( (*N_planets)*sizeof(SSDATA) );
	if((*planet_list) == NULL){
		fprintf(stderr, "Error: Could not allocate memory, aborting...\n");
		abort();
	}
	j=0;
	for(i=0;i<head->n_data;i++){ /*add planets to array*/
		if(data[i].color!=3) (*planet_list)[j] = data[i];
	}
	return(0);
}

int delete_flagged(SSHEAD *head, SSDATA **data){
	/*	Flag for deletion by setting SSDATA.org_idx = -1
		This modifies the data array, it needs a pointer to an array
		of data. e.g. if "array" is a data array, you would pass
		"&array", the address of "array".
	*/
	int i=0, j=head->n_data;
	SSDATA *temp=NULL;
	
	for(i=0;i<head->n_data;i++){ /*count how many will be left*/
		if((*data)[i].org_idx ==-1) j--;
	}
	
	temp = (SSDATA*) malloc( j*sizeof(SSDATA) );
	if(temp==NULL){
		fprintf(stderr, "Error: Could not allocate memory, aborting...\n");
		abort();
	}
	
	j=0;
	for(i=0;i<head->n_data;i++){
		/*printf("*data[i].org_idx: %i, j: %i\n", (*data)[i].org_idx, j);*/
		if((*data)[i].org_idx == -1) continue;
		temp[j] = (*data)[i];
		j++;
	}

	head->n_data = j;
	free(*data);
	*data=temp;
	
	return(0);
}

int orbital_params(SSHEAD *head, SSDATA* data, PARTICLE *particle){
	/*Computes the three main orbital prameters for a particle
	in simulation units G=1, M_sun=1*/
	int i, j;

	for(i=0;i<head->n_data;i++){
		double h[3], unit_h[3];
		double unit_z[3] = {(double)0.0, (double)0.0, (double)1.0};
		double h2, v2, r2;
		
		r2 = dot_prod(data[i].pos, data[i].pos);
		v2 = dot_prod(data[i].vel, data[i].vel);
		
		cross_prod(data[i].pos, data[i].vel, h);
		norm_vec(h, unit_h);
		
		h2 = dot_prod(h, h);
		
		/*Compute semi-major axis*/
		particle[i].sma = 1.0/(2.0/sqrt(r2) - v2/(1.0+data[i].mass));
		
		/*compute eccentricity*/
		particle[i].ecc = sqrt(1.0 - h2/( (1.0+data[i].mass)*particle[i].sma ));
		
		/*compute inclination,
		Assume that ecliptic is z=0 plane
		acos() outputs in RADIANS, the trig functions use RADIANS!!*/
		if(i==0) printf("h = %E %E %E, \nunit_h = %E %E %E\n", h[0],h[1],h[2], unit_h[0],unit_h[1],unit_h[2]);
		particle[i].inc = acos( dot_prod(unit_h, unit_z) );
		
		particle[i].mass = data[i].mass;
		particle[i].radius = data[i].radius;
		particle[i].hill_rad = hill_rad(data[i]);
		for(j=0;j<N_DIM;j++){
			particle[i].pos[j] = data[i].pos[j];
			particle[i].vel[j] = data[i].vel[j];
			particle[i].spin[j] = data[i].spin[j];
		}
		particle[i].mass = data[i].mass;
		particle[i].color = data[i].color;
		particle[i].org_idx = data[i].org_idx;
	}


	return(0);
}


void nn_density_SS(NODE **n_near, double *n_dist2, SSDATA *data, double sigmas, double* location, int n_dim, int n){
	//n_near: an array of nearest neighbours
	//n_dist2: an array of the squared distances of those nearest neighbours, in order
	//scale_len: the length scale of the smoothing, approximately the size of a cell
	//location: where the final answer will be stored (usually a cell in a 2d array)
	//n_dim: how many dimensions we have
	//n: how many nearest neighbours we have
	int i, j = min(n_dist2, n);
	double sum = 0.0, var2, factor;
	var2 = max(n_dist2, n)/(sigmas*sigmas); //assume 'max(n_dist2)' is at 'sigmas' sigma level
	if(var2==0.0){ //then everything is at DBL_MAX, i.e. not included
		*location = 0.0; //prob. can't set to zero
		return;
	}
	factor = 1.0/(sqrt(var2)*sqrt(2.0*M_PI));
	for(i=0;i<n;i++){
		double gauss = factor*exp(-n_dist2[i]/(2.0*var2));
		int p_num = (int)(n_near[i]->val[n_dim]);
		if(i == j){ //very likely to not be a neighbour but the particle
			continue;
		}
		sum += data[p_num].mass*gauss; //this could be a possible re-write
	}
	*location = sum;
}

void nn_num_density_SS(NODE **n_near, double *n_dist2, SSDATA *data, double sigmas, double *location, int n_dim, int n){
	//This should find the number density associated with the n_near list
	int i, j = min(n_dist2, n);
	double sum = 0.0, var2, factor;
	var2 = max(n_dist2, n)/(sigmas*sigmas); //assume 'max(n_dist2) is at a 'sigmas' sigma level
	if(var2==0.0){ //then everything is at DBL_MAX, i.e. not included
		*location = 0.0; //prob. can't set to zero
		return;
	}
	factor = 1.0/(sqrt(var2)*sqrt(2.0*M_PI));
	for(i=0;i<n;i++){
		if(i == j){ //very likely to not be a neighbour but the particle
			continue;
		}
		sum += factor*exp(-n_dist2[i]/(2.0*var2)); //want number density, no scaling
	}
	*location = sum;
}

void nn_vel_disp_SS(NODE **n_near, double *n_dist2, SSDATA *data, double *ref_vel, double *location, int n_dim, int n){
	//finds the velocity dispersion of the nearest neighbours
	int i, j = min(n_dist2, n);
	double sum=0.0;
	for(i=0;i<n;i++){
		double dv[3];
		int p_num = (int)(n_near[i]->val[n_dim]);
		if(i == j){ //very likely to not be a neighbour but the particle
			continue;
		}
		v_minus(ref_vel, data[p_num].vel, dv);
		sum += sqrt(dot_prod(dv, dv));
	}
	*location = sum/(double)n;
}

void nn_vel_disp_point(NODE **n_near, double *n_dist2, SSDATA *data, double *location, int n_dim, int n){
	int i, pnum;
	double sum[3]={0.0, 0.0, 0.0}, sum2=0.0, mean2;
	for(i=0;i<n;i++){
		pnum = (int)(n_near[i]->val[n_dim]);
		v_plus(sum, data[i].vel, sum);
		sum2 += dot_prod(data[i].vel, data[i].vel);
	}
	scalar_multi(sum, 1.0/(double)n, sum);
	mean2 = dot_prod(sum, sum); //square of the mean
	sum2 /= (double)n; //mean of the square
	*location = sum2 - mean2; //mean of the squares - square of the mean
}

double max(double *array, int n){
	int i, j=0;
	double max = 0.0;
	for(i=0;i<n;i++){ //DBL_MAX is what I set things to if I want them ignored
		if((array[i] > max) && (array[i]!=DBL_MAX)){
			max = array[i];
		}
	}
	return(max);
}

int min(double *array, int n){
	int i, j;
	double min = DBL_MAX;
	for(i=0;i<n;i++){
		if(array[i] < min){
			min = array[i];
			j = i;
		}
	}
	return(j);
}

int cross_orbit(DUST_ELLIPSES obj_1, DUST_ELLIPSES obj_2){
	double P1, P2, a, b, c;
	P1 = obj_1.sma*(1.0-obj_1.ecc*obj_1.sma);
	P2 = obj_2.sma*(1.0-obj_2.ecc*obj_2.sma);
	a = P1 - P2;
	b = obj_2.ecc*P1*cos(obj_2.psi) - obj_1.ecc*P2*cos(obj_1.psi);
	c = obj_2.ecc*P1*sin(obj_2.psi) - obj_1.ecc*P2*sin(obj_1.psi);
	
	if(a*a > (b*b + c*c)){
		return(0); //orbits do not cross
	}
	else{
		return(1); //orbits do cross
	}
}

void find_vel(double r, double *ref_point, double *ref_vel){
	//assume circular orbit
	//G=1, solar mass = 1
	double speed = sqrt(1.0/r);
	double cphi, ctheta, stheta, rho;
	cphi = (ref_point[2]/r);
	rho = r*sin(acos(cphi));
	ctheta = (ref_point[0]/rho);
	stheta = ref_point[1]/rho;
	//may have to play with the signs to get the right co-ord system
	ref_vel[2] = speed*cphi;
	ref_vel[0] = speed*rho*ctheta;
	ref_vel[1] = speed*rho*stheta;
	return;
}