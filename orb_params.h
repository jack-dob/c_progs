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