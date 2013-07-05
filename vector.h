#ifndef VECTOR_INCLUDED
#define VECTOR_INCLUDED

#ifndef DEBUG
	#define DEBUG 0
#endif /*DEBUG*/

#define MALLOC_ERROR "Error: %s could not initialise array %s\n"

#include <math.h>

double mag(double a[3]);
void scalar_multi(double a[3], double scalar, double out[3]);
int norm_vec(double a[3], double out[3]);
int cross_prod(double a[3], double b[3], double out[3]);
double dot_prod(double a[3], double b[3]);
double displacement(double a[3], double b[3]);
int v_minus(double a[3], double b[3], double out[3]);
int v_plus(double a[3], double b[3], double out[3]);
void rotate_x(double ***rot_mat, double angle);
void rotate_y(double ***rot_mat, double angle);
void rotate_z(double ***rot_mat, double angle);
void mat_multi(double **a,double **b, double ***output, int a_n, int a_m, int b_n, int b_m, int malloc_flag);
void rotate(double ***a, double *angles, int *axis, int n_ops);
void col_vec_multi(double **a, double *b, double **output,  int a_n, int a_m);
void invert(double **a, double ***out, int order);
double cofactor(double **a, int i, int j, int order);
void minor_matrix(double **a, double ***out, int row, int col, int order);
double det(double **a, int order);
void transpose(double **a, double ***out, int order);

double mag(double a[3]){
	double dp;
	dp = dot_prod(a, a);
	return(sqrt(dp));
}

void scalar_multi(double a[3], double scalar, double out[3]){
	int i;
	for(i=0; i<3; i++){
		out[i] = scalar*a[i];
	}
	return;
}

void transpose(double **a, double ***out, int order){
	int i, j;
	double **temp;	
	
	for(i=0;i<order;i++){
		for(j=0;j<order;j++){
			(*out)[i][j] = a[j][i];
		}
	}
}

double det(double **a, int order){
	int i, j;
	double **minor;
	double sum=0;
	
	if(order==2){
		return(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	}
	else if(order!=3){
		printf("Error: %s only works for order 2 or 3 matricies at the momen...\n");
		exit(0);
	}
	
	minor = (double**) malloc((int)(order*sizeof(double*) ) );
	if(minor==NULL){
		printf(MALLOC_ERROR, __FUNCTION__, "minor");
	}
	for( i=0;i<order;i++){
		minor[i] = (double*) calloc( order, sizeof(double) );
		if(minor[i]==NULL){
			char *str;
			sprintf(str, "minor[%i]", i);
			printf(MALLOC_ERROR, __FUNCTION__, str);
		}
	}
	
	for(i=0;i<order;i++){
		for(j=0;j<order;j++){
			minor_matrix(a, &minor, i, j, order-1);
			sum+=(a[i][j]*det(minor, order-1));
		}
	}
	
	for(i=0;i<order;i++){
		free(minor[i]);
	}
	free(minor);
	return(sum);
}

void minor_matrix(double **a, double ***out, int row, int col, int order){
	/*row, col are the row and column respectively of the element who's minor matrix is
	being computed*/
	int i,j,k,l;
	k=0;
	l=0;
	
	if(DEBUG) printf("a:\n");
	for(i=0;i<order;i++){
			if(DEBUG) printf("%G  %G  %G\n", a[i][0], a[i][1], a[i][2]);
	}
	k=0;
	for(i=0;i<order;i++){
		if(i==row){
			continue;
		}
		l=0;
		for(j=0;j<order;j++){
			if(j!=col){
				(*out)[k][l] = a[i][j];
				l++;
			}
		}
		k++;
	}
}

double cofactor(double **a, int row, int col, int order){
	/*This computes the co-factor of the row-col element of matrix a and returns it as a double*/
	/*this code is for 3X3 matrices only, can expand if necessary*/
	double **minor;
	double determinant=0.0, cofactor=1.0, answer=0.0, det_a = 0.0;
	int k;
	int i, j;
	det_a = det(a, order);
	minor = (double**) malloc( (int)(order-1)*sizeof(double*) );
	if(minor==NULL){
		/*printf(MALLOC_ERROR("minor"));*/
		exit(1);
	}
	for(k=0;k<(order-1);k++){
		minor[k] = (double*) malloc( (int)(order-1)*sizeof(double) );
		if(minor[k]==NULL){
			/*printf(MALLOC_ERROR( sprintf("minor[%i]", k) ));*/
			exit(1);
		}
	}
	
	minor_matrix(a, &minor, row, col, order);
	determinant = det(minor, (order-1) );

	for(k=0;k<(row+col);k++){
		cofactor*=-1.0;
	}
	
	if(DEBUG) printf("Minor matrix:\n");
	for(k=0;k<(order-1);k++){
		if(DEBUG) printf("%G  %G\n", minor[k][0],minor[k][1]);
		free(minor[k]);
	}
	free(minor);
	answer = cofactor*determinant/det_a;
	if(DEBUG) printf("cofactor: %G, determinant: %G, det_a: %G, answer: %G\n", cofactor, determinant, det_a, answer);
	return(answer);
}

void invert(double **a, double ***out, int order){
	int k, i, j;
	double **temp;

	if(DEBUG) printf("Allocating temporary storage...\n");
	temp = (double**) malloc( (int)(order)*sizeof(double*) );
	if(temp==NULL){
		printf(MALLOC_ERROR, __FUNCTION__, "temp");
		exit(1);
	}
	for(k=0;k<(order);k++){
		temp[k] = (double*) calloc( (int)(order), sizeof(double) );
		if(temp[k]==NULL){
			printf(MALLOC_ERROR, __FUNCTION__, "temp");
			exit(1);
		}
	}
	if(DEBUG) printf("\n");
	if(DEBUG) printf("computing cofactors...\n");
	for(i=0;i<order;i++){
		for(j=0;j<order;j++){
			if(DEBUG) printf("\t\t\ta[%i][%i] = %G\n", i, j, a[i][j]);
			temp[i][j] = cofactor(a, i, j, order);
			if(DEBUG) printf("temp[%i][%i] = %G\n", i, j, temp[i][j]);
		}
	}	
	if(DEBUG) printf("transposing temporary array...\n");
	transpose(temp, out, order);
	
	if(DEBUG) printf("Freeing temporary storage...\n");
	for(k=0;k<order;k++){
		free(temp[k]);
	}
	free(temp);
	if(DEBUG) printf("Temporary storage freed...\n");
}

void col_vec_multi(double **a, double *b, double **output, int a_n, int a_m){
	/* 	a = transform, an a_n by a_m matrix
		b = original column vector, an a_n by 1 matrix
		output = transformed vector, a 1 by a_n matrix
		a_n, a_m = dimensions of a*/
	int i, j;
	/*printf("addr of output = %i\n", output);
	printf("*output[0] = %lf, *output[1] = %lf, *output[2] = %lf\n", (*output)[0], (*output)[1], (*output)[2]);
	*/
	/*printf("Before loop\n");*/
	/*coming from 'create_dust_cube' the a array is screwed for some reason,
	try copying it into a different array and see if that helps*/
	/*printf("%E\n", a[0][0]);*/
	for(i=0;i<a_n;i++){
		double sum=0;
		/*printf("In 1st loop\n");*/
		for(j=0;j<a_m;j++){
			/*printf("in 2nd loop\n");*/
			if(DEBUG) printf("DEBUG: a[%i][%i]  = %G ,b[%i] = %G\n", i, j ,a[i][j], j, b[j]);
			sum += a[i][j]*b[j];
			
		}
		/*printf("sum  = %G\n", sum);*/
		(*output)[i] = sum;
	}
}

void mat_multi(double **a,double **b, double ***output, int a_n, int a_m, int b_n, int b_m, int malloc_flag){
	int i, j, k, duplicate_flag;
	double sum;
	double **temp;

	
	/*printf("IN MAT_MULTI\n");*/
	if(malloc_flag==2) printf("\n input matricies a, b:\n");
	if(malloc_flag==2) {
		for(i=0;i<3;i++){
			printf("|%012G  %012G  %012G|\t\t|%012G  %012G  %012G|\n",a[i][0],a[i][1],a[i][2],
			b[i][0],b[i][1],b[i][2]);
		}
	}	
	if(!malloc_flag){/*allocate memory to output array*/
		*output = (double**) malloc( a_n*sizeof(double*) );
		if (*output==NULL){
			printf("Error: mat_multi could not allocate memory to output\n");
			exit(1);
		}
		for(i=0;i<a_n;i++){
			(*output)[a_n] = (double*) malloc( b_m*sizeof(double) );
			if((*output)[a_n]==NULL){
				printf("Error: mat_multi could not allocate memory to output[%i]\n", i);
				exit(1);
			}
		}
	}
	/*	a = left hand matrix,
		b = right hand matrix,
		n = number of rows (of a or b)
		m = number of columns (of a or b)
		to make a row matrix to a column matrix, pass pointer to it
	*/
	if(a_m != b_n){
		printf("Error: your matrix sizes are not consistent, a is %i X %i, b is %i X %i.\n", a_n, a_m, b_n, b_m);
		exit(1);
	}
	if(malloc_flag==2) printf("HERE\n");
	/*printf("START MULTIPLICATION\n");*/
	for(k=0;k<a_n;k++){/*rows of a*/
		sum=0.0;
		for(j=0;j<b_n;j++){/*columns of b*/
			sum=0.0;
			for(i=0;i<b_m;i++){/*sum over rows of b*/
				if(malloc_flag==2) printf("%E * %E + ",a[k][i], b[i][j]);
				sum+=a[k][i]*b[i][j];
			}
			if(malloc_flag==2) printf("\nout[%i][%i] = %E\n", k, j, sum);
			(*output)[k][j] = sum;
		}
	}
	/*printf("END MULTIPLICATION, RETURN\n");*/
}

/*this creates the x-axis rotation matrix for a given angle*/
void rotate_x(double ***rot_mat, double angle){
	int i,j;
	double s_angle = sin(angle);
	double c_angle = cos(angle);
	double temp[3][3]= {	{1.0, 0.0, 0.0},
							{0.0, c_angle, -s_angle},
							{0.0, s_angle, c_angle}
					    };
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			(*rot_mat)[i][j] = temp[i][j];
		}
	}
}

/*this creates the y-axis rotation matrix for a given angle*/
void rotate_y(double ***rot_mat, double angle){
	int i,j;
	double s_angle = sin(angle);
	double c_angle = cos(angle);
	double temp[3][3] = {	{c_angle, 0.0, s_angle},
							{0.0, 1.0, 0.0},
							{-s_angle, 0.0, c_angle}
						};
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			(*rot_mat)[i][j] = temp[i][j];
		}
	}
}

/*this creates the z-axis rotation matrix for a given angle*/
void rotate_z(double ***rot_mat, double angle){
	int i,j;
	double s_angle = sin(angle);
	double c_angle = cos(angle);
	/*printf("INFO: (in rotate_z) cos(%G) = %G, sin(%G) = %G\n", angle, c_angle, angle, s_angle);*/
	double temp[3][3] =	{	{c_angle, -s_angle, 0.0},
							{s_angle, c_angle, 0.0},
							{0.0, 0.0, 1.0}
					};

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			/*printf("temp[%i][%i] = %G\n", i, j, temp[i][j]);*/
			(*rot_mat)[i][j] = temp[i][j];
			/*printf("(*rot_mat)[%i][%i] = %G\n", i, j, (*rot_mat)[i][j]);*/
		}
	}
}

void rotate(double ***a, double *angles, int *axis, int n_ops){
	/*  'a' is a pointer to a 2d matrix that will contain final transformation
		'angles' are the angles you want to rotate by in order of operation
		'axis' are the axis you want to rotate around in order of operation, 0=x, 1=y, 2=z
		'n-ops' is the number of operations
		NOTE: this does not check matrix dimension, make sure they are correct.
		angles must be in RADIANS*/
	int i, j, k;
	double **temp;
	double **rot_mat;
	/****ALLOCATE MEMORY****/
	rot_mat = (double**) malloc( (int)3*sizeof(double*) );
	if(rot_mat==NULL){
		printf("Error: rotate could not allocate rot_mat\n");
		exit(1);
	}
	for(i=0;i<3;i++){
		rot_mat[i] = (double*) malloc( (int)3*sizeof(double) );
		if(rot_mat[i]==NULL){
			printf("Error: rotate could not allocate rot_mat[%i]\n", i);
			exit(0);
		}
	}
	
	temp = (double**) malloc( (int)3*sizeof(double*) );
	if(temp==NULL){
		printf("Error: rotate could not allocate temp\n");
		exit(1);
	}
	for(k=0;k<3;k++){
		temp[k] = (double*) malloc( (int)3*sizeof(double) );
		if(temp[k]==NULL){
			printf("Error: rotate could not allocate temp[%i]\n", k);
			exit(0);
		}
		for(j=0;j<3;j++){
			temp[k][j] = (*a)[k][j];
		}
	}
	/****END MEMORY ALLOCATION****/
	if(DEBUG) printf("DEBUG: (in rotate) address of transform %i\n", a);
	
	if(DEBUG){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("DEBUG: (in rotate, at start) transform[%i][%i] = %G\n", i, j, (*a)[i][j]);
			}
		}
	}
	/*printf("INFO: n_ops = %i\n", n_ops);*/
	for(i=0;i<n_ops;i++){
		/*printf("INFO: i = %i\n", i);*/
		switch (axis[i]){
			case 0:
				/*printf("axis %i, angle = %G\n", axis[i], angles[i]);*/
				/*find rotation matrix*/
				rotate_x(&rot_mat, angles[i]);
				/*output debug info if needed*/
				if(DEBUG){
					for(k=0;k<3;k++){
						for(j=0;j<3;j++){
							printf("DEBUG: (in rotate) rot_mat[%i][%i] = %G\n", k, j, rot_mat[k][j]);
						}	
					}
				}
				/*perform rotation*/
				mat_multi(rot_mat, temp, a, 3, 3, 3, 3, 1);
				/*set temporary array to reflect changes*/
				for(k=0;k<3;k++){
					for(j=0;j<3;j++){
						temp[k][j] = (*a)[k][j];
					}
				}
				break;
			case 1:
				/*printf("axis %i, angle = %G\n", axis[i], angles[i]);*/
				/*find rotation matrix*/
				rotate_y(&rot_mat, angles[i]);
				/*output debug info if needed*/
				if(DEBUG){
					for(k=0;k<3;k++){
						for(j=0;j<3;j++){
							printf("DEBUG: (in rotate) rot_mat[%i][%i] = %G\n", k, j, rot_mat[k][j]);
						}	
					}
				}
				/*perform rotation*/
				mat_multi(rot_mat, temp, a, 3, 3, 3, 3, 1);
				/*set temporary array to reflect changes*/
				for(k=0;k<3;k++){
					for(j=0;j<3;j++){
						temp[k][j] = (*a)[k][j];
					}
				}
				break;
			case 2:
				/*printf("axis %i, angle = %G\n", axis[i], angles[i]);*/
				/*find rotation matrix*/
				rotate_z(&rot_mat, angles[i]);
				/*output debug info if needed*/
				if(DEBUG){
					for(k=0;k<3;k++){
						for(j=0;j<3;j++){
							printf("DEBUG: (in rotate) rot_mat[%i][%i] = %G\n", k, j, rot_mat[k][j]);
						}	
					}
				}
				/*perform rotation*/
				mat_multi(rot_mat, temp, a, 3, 3, 3, 3, 1);
				/*set temporary array to reflect changes*/
				for(k=0;k<3;k++){
					for(j=0;j<3;j++){
						temp[k][j] = (*a)[k][j];
					}
				}
				break;
		}
		if(DEBUG){
			for(k=0;k<3;k++){
				for(j=0;j<3;j++){
					printf("DEBUG: (in rotate, loop %i) transform[%i][%i] = %G\n", i, k, j, (*a)[k][j]);
				}
			}
		}
		/*printf("restart loop...\n");*/
	}
	/*printf("THERE\n");*/
	for(i=0;i<3;i++){
		free(rot_mat[i]);
	}
	free(rot_mat);
	if(DEBUG){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				printf("DEBUG: (in rotate, at end) transform[%i][%i] = %G\n", i, j, (*a)[i][j]);
			}
		}
	}
}

int v_plus(double a[3], double b[3], double out[3]){
	out[0] = b[0] + a[0];
	out[1] = b[1] + a[1];
	out[2] = b[2] + a[2];
	return(0);
}

int v_minus(double a[3], double b[3], double out[3]){
	
	out[0] = b[0] - a[0];
	out[1] = b[1] - a[1];
	out[2] = b[2] - a[2];
	return(0);
}

double displacement(double a[3], double b[3]){
	double out;
	double diff[3];
	
	v_minus(a, b, diff);
	
	out = sqrt(dot_prod(diff, diff));

	return(out);
}

int cross_prod(double a[3], double b[3], double out[3]){
	/* |out| = |a||b|sin(angle)*/
	
	out[0] = a[1]*b[2]-b[1]*a[2];
	out[1] = b[0]*a[2]-a[0]*b[2];
	out[2] = a[0]*b[1]-b[0]*a[1];
	
	return(0);
}

double dot_prod(double a[3], double b[3]){
	/* c = |a||b|cos(angle) */
	double c=0.0;
	int i;
	
	for(i=0;i<3;i++){
		c+=a[i]*b[i];
	}
	
	return(c);
}

int norm_vec(double a[3], double out[3]){
	int i;
	double mag_a = sqrt(dot_prod(a,a));
	
	for(i=0;i<3;i++){
		out[i] = a[i]/mag_a;
	}

	return(0);
}

#endif /*VECTOR_INCLUDED*/