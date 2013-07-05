#include <stdio.h>
#include <math.h>
#include "p_tree.h"
#include "sshead.h"
#include "dust_params_p.h"

#define XPIX 1000.0
#define YPIX 1000.0
#define XSCALE 10.0
#define YSCALE 10.0
#define XMIN -5.0
#define XMAX 5.0
#define YMIN -5.0
#define YMAX 5.0
#define ORB_REZ 1000

int main(int *argc, char **argv){
	char *infile, *outfile;
	SSDATA *ssdata;
	SSHEAD sshead;
	PARTICLE **p_data;
	NODE *root;
	double **t_coll_map;
	double **N_part_map;

	if(argc>1)
		infile =  argv[1];
	if(argc>2)
		outfile = argv[2];
	
	//want to write out a map of t_coll for a given ss-file
	//read in data
	readin_sshead(infile, &sshead);
	ssdata = readin_ssdata(infile, &head);
	
	//change to correct format for p_tree
	ssdata2ptree(&ssdata, sshead, &p_data); // translates and frees previous version
	
	//allocate the root node
	root = (NODE*) malloc((int)sizeof(NODE));
	create_tree(&root, p_data, sshead.n_data, NDIM);
	
	//get t_coll map
	get_t_coll(p_data, &t_coll_map, &N_part_map, sshead, root);
	
	//write out t_coll map
	return(0);
}

void get_t_coll(PARTICLE **p_data, double ***t_coll_map, double ***N_part_map, SSHEAD sshead, NODE *root){
	int i, j;
	double inv_transform[3][3];
	double angle_vec[3];
	double **t_coll_m;
	double **N_part_m;
	double t_coll;
	
	t_coll_m = (double**) malloc(YPIX*sizeof(double*));
	N_part_m = (double**) malloc(YPIX*sizeof(double*));
	for(i=0;i<XPIX;i++){ //allocate maps
		t_coll_m[i] = (double*) malloc(YPIX*sizeof(double));
		N_part_m[i] = (double*) malloc(YPIX*sizeof(double));
		for(j=0;j<YPIX;j++){	
			t_coll_m[i][j] = 0.0;
			N_part_m[i][j] = 0;
		}
	}
	
	for(i=0;i<sshead.n_data;i++){
		rot_to_ecliptic_p(p_data[i], &angle_vec, &inv_transform); //rotate each particle to the ecliptic
		t_coll = t_coll_p(&sshead, p_data, root, head.n_data, i);
		//write t_coll and number of objects per pixel to maps
		write_maps(&t_coll_m, N_part_m, p_data[i], t_coll);
		
	}
	//divide t_coll_map by N-part_map to get average t_coll per pixel
	div_grid(&t_coll_m, &N_part_m);
	t_coll_map = t_coll_m;
	N_part_map = N_part_m;
	return;
}

void div_grid(double ***t_coll_m, double ***N_part_m){
	//div arg1 by arg2
	int i, j;
	for(i=0;i<XPIX;i++){
		for(j=0;j<YPIX;j++){
			(*t_coll_m)[i][j] /= (*N_part_m)[i][j];
		}
	}
	return;
}

void write_maps(double ***t_coll_m, double ***N_part_m, PARTICLE *p_dat, double t_coll){
	int i;
	double phi=0.0;
	double P = r_a_p(p_dat)*r_p_p(p_dat)/sma_p(p_dat);
	double smia = sqrt(sma_p(p_dat)*P);
	double temp_pos[3];
	for(i=0;i<ORB_RES;i++){
		double r = P/(1.0 - ecc_p(p_dat)*cos(phi));
		double theta = phi + psi_p(p_dat);
		double d_phi = 2.0*M_PI*sma_p(p_dat)*smia/(r*r*(double)ORB_RES);
		double x = r*sin(theta);
		double y = r*cos(theta);
		double z = 0.0;

		temp_pos[0] = x;
		temp_pos[1] = y;
		temp_pos[2] = z;
		add2grid(t_coll_m, temp_pos, t_coll);
		add2grid(N_part_m, temp_pos, 1.0);
		phi+=d_phi;
	}
	return;
}

void add2grid(double ***grid, double *pos, double value){
	int i;
	
	for(i=0;i<ORB_RES;i++){
		int xbin, ybin, zbin;
		xbin = (int)floor((pos[0] - XMIN)*XPIX/(XMAX-XMIN));
		ybin = (int)floor((pos[1] - YMIN)*YPIX/(YMAX-YMIN));

		if( (xbin >= XPIX) || (ybin >= YPIX) || (xbin<0) || (ybin<0) ){
			/*printf("Warning: could not write dust blob, %i %i %i is outside range %i %i...\r",
				xbin, ybin, zbin, 0, PAR->cube_res);
			*/
			continue;		
		}
		//printf("\ndust_blob[%i] bins: %i, %i, %i, cube_res: %g, cube_z_res: %g \n",i, xbin, ybin, zbin, PAR->cube_res, PAR->cube_z_res);
		/*cfits uses fortrans style ordering, will have to re-format this data later*/

		(*grid)[xbin][ybin] += value;
	}
	
	return;
}

void ssdata2ptree(SSDATA **ssdata_out, SSHEAD sshead, PARTICLE ***p_data_out){
	int i, j;
	SSDATA *ssdata;
	PARTICLE **p_data;
	
	ssdata = *ssdata_out; //for ease of use
	p_data = (PARTICLE**) malloc( sshead.n_data*sizeof(PARTICLE*));
	for(i=0;i<sshead.n_data;i++){
		p_data[i] = (PARTICLE*) malloc( sizeof(PARTICLE) );
		p_datum[i]->mass = ssdata[i].mass;
		for(j=0;j<NDIM;j++)
			p_data[i]->pos[j] = ssdata[i].pos[j];
		for(j=0;j<NDIM;j++)
			p_data[i]->vel[j] = ssdata[i].vel[j];
		for(j=0;j<NDIM;j++)
			p_data[i]->spin[j] = ssdata[i].spin[j];
		p_data[i]->ident = i;
	}
	FREE(ssdata);
	*ssdata_out = ssdata;
	*p_data_out = p_data;
	return;
}


double t_coll_p(SSHEAD *head, PARTICLE **data, NODE *root, int N, int i){
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
	search_data[0] = sma_p(data[i]);
	if(search_data[0] <= 0.0){ // if particle is ejected, it does not count towards dust.
		free(search_data);
		return(0.0);
	}
	search_data[1] = ecc_p(data[i]);
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
		int p_num = (int)n_near[j]->val->ident;
		if(j==l){ // don't count the actual particle.
			continue;
		}
		if( sma_p(data[p_num]) < 0.0 ){
			n_dist2[j] = DBL_MAX; //make lose all weighting if neighbour is ejected
		}
		if(cross_orbit_p(data[p_num], data[p_num])){
		//if(condition){
			k++;
			continue;
		}
		else{
			n_dist2[j] = DBL_MAX; //cannot put at infinity, this will have to do. Will lose all weighting
		}
	}
	

	
	//puts("4");
	nn_density_p(n_near, n_dist2, data, 5.0, &density, N_dim, n);
	nn_num_density_p(n_near, n_dist2, data, 5.0, &num_density, N_dim, n);
	nn_vel_disp_p(n_near, n_dist2, data, data[i]->vel, &vel_disp, N_dim, n);
	
	//puts("5");
	
	if(num_density==0.0){ //fixed
		av_m_coll = DBL_MAX;
	}
	else if(density==0.0 && num_density==0.0){ //fixed
		av_m_coll=1.0;
	}
	else{
		av_m_coll = density/num_density;
	}
	p_0 = 3.36627518E6; //2000.0 kg m^-3 in M_sol AU^-3
	R1 = cbrt(3.0*data[i]->mass/(4.0*p_0*M_PI)); //in AU
	R2 = cbrt(3.0*av_m_coll/(4.0*p_0*M_PI)); //in AU
	v2_esc = 2.0*1.0*(data[i]->mass + av_m_coll)/(R1+R2); //G=1 in system units, escape velocity for target + collider
	x_section = M_PI*(R1+R2)*(R1+R2)*(1.0+v2_esc/(vel_disp*vel_disp));
	
	if(num_density==0){ //fixed
		t_coll = DBL_MAX;
		l_coll = DBL_MAX;
	}
	else{
		t_coll = 1.0/(sqrt(dot_prod(data[i]->vel, data[i]->vel))*num_density*x_section); //gives time for 1 collision
		l_coll = 1.0/(x_section*num_density); //gives average distance between collisions
	}
	h_rad = hill_rad_p(data[i]); //l_coll/h_rad is how likely something is to collide with an object in it's hill radius
	N_coll = data[i]->mass/av_m_coll; //number of collisions = target mass / average collision mass
	r_orb = sqrt(dot_prod(data[i]->pos, data[i]->pos));
	t_orb = cbrt(4.0*M_PI*(r_orb*r_orb*r_orb));

	free(n_near);
	free(n_dist2);
	free(search_data);
	return(t_orb/(N_coll*t_coll)); //works also

}



void rot_to_ecliptic_p(PARTICLE *data, double **angle_vec, double ***inv_transform){
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

void nn_density_p(NODE **n_near, double *n_dist2, PARTICLE **data, double sigmas, double* location, int n_dim, int n){
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
		int p_num = (int)(n_near[i]->val->ident);
		if(i == j){ //very likely to not be a neighbour but the particle
			continue;
		}
		sum += data[p_num]->mass*gauss; //this could be a possible re-write
	}
	*location = sum;
}

void nn_num_density_p(NODE **n_near, double *n_dist2, PARTICLE **data, double sigmas, double *location, int n_dim, int n){
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

void nn_vel_disp_p(NODE **n_near, double *n_dist2, PARTICLE **data, double *ref_vel, double *location, int n_dim, int n){
	//finds the velocity dispersion of the nearest neighbours
	int i, j = min(n_dist2, n);
	double sum=0.0;
	for(i=0;i<n;i++){
		double dv[3];
		int p_num = (int)(n_near[i]->val->ident);
		if(i == j){ //very likely to not be a neighbour but the particle
			continue;
		}
		v_minus(ref_vel, data[p_num]->vel, dv);
		sum += sqrt(dot_prod(dv, dv));
	}
	*location = sum/(double)n;
}

int cross_orbit_p(PARTICLE obj_1, PARTICLE obj_2){
	double P1, P2, a, b, c;
	P1 = sma_p(obj_1)*(1.0-ecc_p(obj_1)*sma_p(obj_1));
	P2 = sma_p(obj_2)*(1.0-ecc_p(obj_2)*sma_p(obj_2));
	a = P1 - P2;
	b = ecc_p(obj_2)*P1*cos(psi_p(obj_2)) - ecc_p(obj_1)*P2*cos(psi-p(obj_1);
	c = ecc_p(obj_2)*P1*sin(psi_p(obj_2)) - ecc_p(obj_1)*P2*sin(psi_p(obj_1));
	
	if(a*a > (b*b + c*c)){
		return(0); //orbits do not cross
	}
	else{
		return(1); //orbits do cross
	}
}

double sma_p(PARTICLE *particle){
	double r2, v2, sma;
	
	r2 = dot_prod(particle->pos, particle->pos);
	v2 = dot_prod(particle->vel, particle->vel);
	
	sma = 1.0/(2.0/sqrt(r2) - v2/(1.0+particle->mass));
	return(sma);
}

double hill_rad_p(PARTICLE *particle){
	/*units in solar masses, so M_sun = 1
	*/
	double hill_rad;
	
	hill_rad = cbrt(particle->mass/3.0)*sma(particle);

	return(hill_rad);
}

double ecc_p(PARTICLE *particle){
	double h[3], h2, ecc;
	cross_prod(particle->pos, particle->vel, h);
	h2 = dot_prod(h, h);
	ecc = sqrt(1.0 - h2/( (1.0+particle->mass)*sma(particle) ));
	return(ecc);
}

double r_a_P(PARTICLE *particle){
	return((1.0+ecc(particle))*sma(particle));
}

double r_p_p(PARTICLE *particle){
	return((1.0-ecc(particle))*sma(particle));
}

double psi_p(PARTICLE *particle){
	double theta, p, psi, i, r=0.0;
	/*assume is in ecliptic plane as co-ord transform was performed earlier*/
	for (i=0;i<NDIM;i++){
		r += particle->pos[i]*particle->pos[i];
	}
	theta = atan(particle->pos[0]/particle->pos[1]);
	p = sma_p(particle)*(1.0 - ecc_p(particle)*ecc_p(particle));
	psi = theta - acos( (p - r)/(r*ecc_p(particle)) );
	return(psi);
}


