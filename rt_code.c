#include <stdio.h>

#define SB_CONST 6E-8
#define PI 3.142
#define AU 1.5E11
#define N_BLOBS 10
#define C_REF 0.05
#define C_TRANS 0.9

typedef struct star{
	/*assume that star is at the centre*/
	double mass;
	double radius;
	double temp;
	double surf_flux;
} STAR;

typedef struct blob{
	double inner_radius;
	double outer_radius;
	double height;
	double c_abs; /*absorption coefficient*/
	double c_ref; /*reflection coefficient*/
	double c_trans; /*transmission coefficient*/
	double inner_flux;
	double outer_flux;
	double temp;
} BLOB;

double stefan_boltzmann(double temp, double area);
void blob_temp(STAR *star, BLOB blob[10]);

int main(int argc, char **argv){
	
	STAR sun;
	BLOB blob[N_BLOBS];
	int i;
	
	/*all units in SI*/
	sun.mass = 2E30;
	sun.radius = 7E8;
	sun.temp = 5778;
	sun.surf_flux = SB_CONST*sun.temp*sun.temp*sun.temp*sun.temp;
	
	printf("Sun:\nMass: %G, Radius: %G, Temp: %G, surf_flux: %G\n",
			sun.mass, sun.radius, sun.temp, sun.surf_flux);
	
	/*setup of blobs*/
	for(i=0;i<(N_BLOBS);i++){
		blob[i].inner_radius = 0.1*((double)i)*AU+0.5*AU;
		blob[i].outer_radius = 0.1*((double)i)*AU+0.6*AU;
		blob[i].height = 0.1*AU;
		blob[i].c_ref = C_REF;
		blob[i].c_trans = C_TRANS;
		blob[i].c_abs = 1.0 - (C_REF+C_TRANS);
		printf("Blob %i:\nr_i: %G, r_o: %G, h: %G, c_abs: %G, c_ref: %G\n",
				i, blob[i].inner_radius/AU, blob[i].outer_radius/AU, blob[i].height/AU, blob[i].c_abs, blob[i].c_ref);
	}
	
	blob_temp(&sun, blob);
	
	for(i=0;i<(N_BLOBS);i++){
		printf("Blob %i:\nr_i: %G, r_o: %G, h: %G, c_abs: %G, c_ref: %G, temp: %G\n",
				i, blob[i].inner_radius/AU, blob[i].outer_radius/AU, blob[i].height/AU,
				blob[i].c_abs, blob[i].c_ref, blob[i].temp);
	}
	
	return(0);
}

void blob_temp(STAR *sun, BLOB blob[N_BLOBS]){
	int i;
	double flux_at_start;
	
	flux_at_start = sun->surf_flux/(4.0*PI*blob[0].inner_radius*blob[0].inner_radius);
	printf("%G\n", flux_at_start);
	for(i=0;i<(N_BLOBS);i++){
		printf("Blob %i:\nr_i: %G, r_o: %G, h: %G, c_abs: %G, c_trans: %G, c_ref: %G\n",
				i, blob[i].inner_radius/AU, blob[i].outer_radius/AU, blob[i].height/AU, blob[i].c_abs,
				blob[i].c_trans, blob[i].c_ref);
		
	}
}


double stefan_boltzmann(double temp, double area){
	double power;
	power = SB_CONST*temp*temp*temp*temp*area;
	return(power);
}