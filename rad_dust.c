#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "fitsio.h" //must use compiler flag '-l cfitsio'
//	Make a load of functions to read in .fits files

#define BUFF_SIZE 256

void bin2pic(double **data, long *naxes, double *iscale, int nbins, double *d_bins);
int str_end(char *str, char *end);
int str_len(char *str);
void split_str(char *str, char delim, int **pos);
void read_dust(char* fname, double **d_bins, int *nbins);
void rad_smooth(double *data, long *naxes, double *iscale);
void min_max(double *data, long *naxes);
void get_out_file(char out_file[128], char* filename, char* prefix);
void flatten_cube(double *data, long int *naxes);
int sizeof_bitpix(int bitpix);
void printerror( int status);
void vec_wrap(long int *fpixel, long int *naxes, int naxis, int add_num);
void write_ds9_cube(char* outfile, long *naxes, double *iscale, double *ds9_cube);

int main(int argc, char **argv){
	//use main to test the reading of a 3d data cube
	//int i, j, k, pnum;
	int i;
	double *d_bins=NULL;
	int n_bins = 15;
	fitsfile *fitsptr;
	char *filename=NULL;
	char *out_file=NULL;
	int status=0;
	//int bitpix, bpix_size; //the type of data stored in the image
	//int max_elem;
	//int naxis=8;
	long naxes[3];
	//long fpixel[3]={1, 1, 1};
	//long fpixel = 1;
	//long nelements;
	double *data=NULL;
	double *iscale=NULL;
	//int anynul;
	char buf[100];
	
	filename = argv[1]; //first argument is the file name
	out_file = argv[2]; //second argument is output file
	puts("A");
	read_dust(filename, &d_bins, &n_bins);
	
	
	//do stuff here with the whole array
	
	// RADIAL SMOOTH
	//get_out_file(out_file, filename, "rad_");
	iscale = (double*) malloc(3*(int)sizeof(double));
	iscale[0] = 10;
	iscale[1] = 10;
	iscale[2] = 0.01; //placeholder
	naxes[0]=1000;
	naxes[1]=1000;
	naxes[2]=1;
	puts("A");
	for(i=0;i<n_bins;i++){
		printf("d_bins[%d]: %G\n", i, d_bins[i]);
	}
	bin2pic(&data, naxes, iscale, n_bins, d_bins);
	puts("A");
	//get_out_file(out_file, filename, "radbins_");
	puts("A");
	write_ds9_cube(out_file, naxes, iscale, data);
	puts("A");
	free(iscale);
	
	//close files
	fits_close_file(fitsptr, &status);
	//free arrays
	free(data);
	free(d_bins);
	return(0);
}

void read_dust(char* fname, double **d_bins, int *nbins){
	int i;
	int *pos=NULL;
	FILE *f=NULL;
	char buffer[BUFF_SIZE]; //hopefully lines won't be over 128 characters long
	
	if (!str_end(fname, ".dust")){
		puts("Error: Not correct file extension '.dust'");
		return;
	}
	f = fopen(fname, "r"); //open file for reading
	//puts("R");
	*nbins=0;
	while(fgets(buffer, BUFF_SIZE, f)!=NULL){ //NULL indicates problem or end of file
		//assume this is one line
		(*nbins)++;
	}
	rewind(f); //go back to start of file
	//puts("R");
	(*d_bins) = (double*) realloc((*d_bins), (int)(*nbins)*sizeof(double) );
	//printf("nbins: %d\n",(*nbins));
	for(i=0;i<(*nbins);i++){
		//printf("i:%d\n", i);
		fgets(buffer, BUFF_SIZE, f);
		//puts("L");
		split_str(buffer, ' ', &pos); //split the buffer by spaces and store split positions in pos		
		//printf("pos[1]: %d, &buffer[pos[1]]: %s\n", pos[1], &(buffer[pos[1]]));
		//printf("pos[1]: %d\n", pos[1]);
		//exit(0);
		(*d_bins)[i] = strtod(&(buffer[pos[0]]), NULL); //get double from string
		printf("i: %d, d_bins[%d]: %G\n", i, i, (*d_bins)[i]);
	}
	//puts("R");
	return;
}

void split_str(char *str, char delim, int **pos){
	//returns an integer array of starting points for each string
	//replaces the delimiter with '\0'
	
	int i, j, n;
	
	for(n=0;str[n]!='\0';n++); //finds length of string
	for(i=0;i<n;i++){
		if(str[i]==delim){
			j++;
		}
	}
	(*pos) = (int*) realloc((*pos), j*sizeof(int));
	j=0;
	for(i=0;i<n;i++){
		if(str[i]==delim){
			str[i]='\0';
			(*pos)[j] = i+1;
			j++;
		}
	}
	return;
}

int str_len(char *str){ //find the length of the string
	int i=0;
	for(i=0;str[i]!='\0';i++);
	return(i);
}

int str_end(char *str, char *end){
	int i, j, k;
	j = str_len(end); //length of suffix string
	i = str_len(str) - j; //where the suffix should start in test string

	for(k=0;k<j;k++){
		if(str[i] != end[k]){
			break;
		}
		i++;
	}
	if(k<j)
		return(0);
	else if(k==j)
		return(1);
	else{
		fprintf(stderr, "Error: %s did not perform correctly\n", __FUNCTION__);
		exit(1);
	}
}

void bin2pic(double **data, long *naxes, double *iscale, int nbins, double *d_bins){
	int i, j, k;
	int nelements=1, nbin=(nbins-1), bin; //nbins-1 as the last is the trash bin
	double x, y, z, r;
	double binstart = 0.7, binend = 3.9;
	double binstep, area, test1;
	double xoff = -5, yoff = -5; //bottom left corner
	
	binstep = (binend-binstart)/(double)nbin;
	
	printf("binstep: %g, binstart: %g, binend: %g\n", binstep, binstart, binend);
	
	for(i=0;i<3;i++){
		nelements*=(int)naxes[i];
	}
	
	(*data) = (double*) realloc((*data), nelements*sizeof(double)); //allocate data array
	//for(i=0;i<nbin;i++){
	//	printf("bin_array[%d]: %G\n", i, bin_array[i]);
	//}
	puts("creating picture");
	for(i=0;i<nelements;i++){
		printf("progrsss: %4.2lf%%\r", 100.0*(double)i/(double)nelements);
		x = (double)(i%naxes[0]);
		//printf("xi: %g", x);
		x = (iscale[0]*x/(double)naxes[0]) - iscale[0]/2.0; //should be in AU from centre now
		y = (double)floor(((double)i)/(double)naxes[0]);
		//printf("\tyi: %g\n", y);
		y = (iscale[1]*y/(double)naxes[1]) - iscale[1]/2.0; //should be in AU from centre now
		//z = floor((float)i/(naxes[0]*naxes[1]));
		//z = (z/iscale[2]) - iscale[2]/2.0;
		z = 0.0; //can do this for thin disk
		r = sqrt(x*x + y*y + z*z);
		bin = (int)floor((r - binstart)/binstep);
		area = M_PI*((double)bin*binstep+binstart)*binstep;
		//printf("i: %d, x: %g, y: %g, z: %g, r: %g, bin: %d, binstep: %g, nbins: %d, nbin: %d\n",\
		//i, x, y, z, r, bin, binstep, nbins, nbin);
		if((bin >= 0) && (bin < nbin)){
			(*data)[i] = d_bins[bin];///area;
		}
		else
			(*data)[i] =0.0;
	}
	return;
}

void min_max(double *data, long *naxes){
	int i, nelements=1;
	double min = DBL_MAX;
	double max = DBL_MIN;
	//double temp;
	
	for(i=0;i<3;i++){
		nelements*=naxes[i];
	}
	
	for(i=0;i<nelements;i++){
		//temp = data[i];
		if (data[i] > max){
			max = data[i];
		}
		if(data[i] < min){
			min = data[i];
		}
	}
	//printf("min: %E\t\tmax: %E\n", min, max);
	printf("%E\n",max);
	return;
}

void get_out_file(char out_file[128], char* filename, char* prefix){
	int i, len;
	char temp[128]; 
	
	strcpy(temp, prefix);//decide on prefix for file
	
	for(len=0;filename[len]!='\0';len++); //find the lenght of the filename
	//printf("filename: %s\n", filename);
	//printf("out_file: %s\n", out_file);
	//printf("temp: %s\n", temp);
	//printf("len: %d\n", len);
	
	//find where the last '/' character is, signals start of file
	//if no '/' character use the start of the string
	for(i=len;filename[i]!='/';i--){
		if(i==0)
			break;
	} 
	if(i!=0) //if we're not at the start add one on (only want the bit after the slash)
		i++;
	
	//printf("filename[%d]> : %s\n", i, &(filename[i]));
	strcat(temp, &(filename[i])); //concatenate prefix and filename
	//printf("temp: %s\n", temp);
	
	strcpy(out_file, filename); //copy the filename to the output file
	out_file[i] = '\0'; //make the string end at the slash (now is only the path)
	strcat(out_file, temp); //copy out prefixed filename after the slash
	
	//printf("out_file: %s\n", out_file);
	return;
}

void flatten_cube(double *data, long *naxes){
	//this will flatten a cube onto one plane, will destroy some data!!
	int i, j, square, nelements=1;
	
	puts("Flattening cube...");
	
	for(i=0;i<3;i++){ //add the 
		nelements*=naxes[i];
		printf("naxes[%d]: %d, nelements: %d\n", i, naxes[i], (int)nelements);
	}
	square = naxes[0]*naxes[1];
	//return;
	for(i=square;i<nelements;i++){
		j = i%square;
		//printf("i: %d, j: %d\n", i, j);
		if(j==0){
			printf("progress: %0.2lf\r", (double)i/(double)nelements);
		}
		data[j] += data[i]; //add the other squares to the first one
	}
	
	naxes[2] = 1; //have flattened image
	return;
}

int sizeof_bitpix(int bitpix){
	//returns the machine dependent size of the type of data that bitpix tells us to use
	switch(bitpix){
		case BYTE_IMG:
			return(sizeof(char));
		case SHORT_IMG:
			return(sizeof(short int));
		case LONG_IMG:
			return(sizeof(int));
		case LONGLONG_IMG:
			return(sizeof(long int));
		case FLOAT_IMG:
			return(sizeof(float));
		case DOUBLE_IMG:
			return(sizeof(double));
		default:
			fprintf(stderr, "Error: Could not determine data type from bitpix=%d\n", bitpix);
			return(0);
	
	}
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

void vec_wrap(long int *fpixel, long int *naxes, int naxis, int add_num){
	int i, nboxes=1;
	int *temp;
	int current_pix=(fpixel[0]-1);
	int max_pix=naxes[0];
	
	for(i=1;i<naxis;i++){
		nboxes *= naxes[i-1];
		current_pix += (fpixel[i]-1)*nboxes;
		max_pix *= naxes[i];
	}
	//fprintf(stderr, "add_num: %d current_pix: %d, max_pix: %d\n", add_num, current_pix, max_pix);
	if(current_pix+add_num > max_pix){
		fprintf(stderr, "ERROR:\tadding %d will take the array out of bounds\n", add_num);
		fprintf(stderr, "\tcurrent_pix: %d, max_pix: %d\n", current_pix, max_pix);
		exit(1);
	}
	
	temp = (int*) calloc( naxis, sizeof(int) );
	//FIX THIS!!!!
	add_num += current_pix; //start from zero, is the best way
	temp[0] = (add_num % (int)naxes[0]); //number to add to the 1st axis
	//printf("add_num %d, naxes[0] %d, temp[0] %d\n", add_num, naxes[0], temp[0]);
	nboxes=1;
	for(i=1;i<(naxis-1);i++){ 
		nboxes *= naxes[i-1]; //C indexed from 0
		temp[i] = (add_num/nboxes) % (int)naxes[i]; //number to add to the 2nd, 3rd, ... n-1th axis
	}
	
	nboxes *= naxes[naxis-2];//naxis-2 because C is indexed from 0
	temp[naxis-1] = add_num/nboxes; //hopefully integer division will work nice here
	
	for(i=0;i<naxis;i++){ //add the 
		fpixel[i] = (long int)(temp[i]+1);
		//printf("fpixel[%d]: %d\t", i, fpixel[i]);
		//printf("naxes[%d]: %d\n", i, naxes[i]);
	}
	//puts("");
	return;
}

void write_ds9_cube(char* outfile, long *naxes, double *iscale, double *ds9_cube){
	//outfile: string containing output file name
	//naxes: array containing the number of elements in each dimension
	//iscale: array containing the expanse of each dimension
	//ds9_cube: double array containing data to go into fits image
	long naxis = 3; /*dimension of image*/
	int i, j, k, l;
	long first_pixel;
	long nelements=1;
	fitsfile *fitsptr;
	int status=0;
	int bitpix = DOUBLE_IMG; /*image is in double format*/

	long fpixel[3] = {1, 1, 1};
	//long naxes[3] = {nx, ny, nz};

	double *out_array; /*need to put into 1d format as cfits likes to be annoying*/
	char *bunit="mJy", *ctype="offset (arcseconds)";
	float crpix1=((float)naxes[0])/2.0, crval=0.0, cdelt1=((float)iscale[0])/((float)naxes[0]);
	float crpix2=((float)naxes[1])/2.0, cdelt2=((float)iscale[1])/((float)naxes[1]);
	float crpix3=((float)naxes[2])/2.0, cdelt3=((float)iscale[2])/((float)naxes[2]);
	
	for(i=0;i<naxis;i++){
		nelements*=naxes[i];
	}
	
	
	printf("Creating FITS file...\n");
	if( !fits_create_file(&fitsptr, outfile, &status) ){
		//Write the fits file here
		
		first_pixel = 1; //first pixel to write
		
		printf("Constants found, creating FITS image...\n");
		fits_create_img(fitsptr,  bitpix, naxis, naxes, &status);
		if(status){
			printf("Error: There was a problem creating the image to the fits file.\n");
			return;
		}
		printf("Writing FITS image...\n");
		/*fits_write_img(fitsptr, TDOUBLE, first_pixel, nelements, dust_cube[0], &status);
		*/
		fits_write_pix(fitsptr, TDOUBLE, fpixel, nelements, ds9_cube, &status);
		if(status){
			printf("Error: There was a problem writing the image to the fits file.\n");
			return;
		}
		printf("Writing header data...\n");
		//fits_update_key(fitsptr, TSTRING, "INPUT_FILE", infile,
		//	"File this cube was created from", &status);
		fits_update_key(fitsptr, TSTRING, "BUNIT", bunit,
			"", &status);
		fits_update_key(fitsptr, TSTRING, "CTYPE1", &ctype,
			"", &status);
		fits_update_key(fitsptr, TSTRING, "CTYPE2", &ctype,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CRPIX1", &crpix1,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CRPIX2", &crpix2,
			"", &status);
			
		fits_update_key(fitsptr, TFLOAT, "CRVAL1", &crval,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CRVAL2", &crval,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CDELT1", &cdelt1,
			"", &status);
		fits_update_key(fitsptr, TFLOAT, "CDELT2", &cdelt2,
			"", &status);
         if(status){
			printf("Error: There was a problem writing the origin file to the header.\n");
			return;
		}
		printf("Closing FITS file...\n");
		if(fits_close_file(fitsptr, &status)){
			printf("Warning: Could not close file successfully\n");
		}
	}
	if(status){
		printf("Error: There was a problem writing the fits file, it may be corrupted.\n");
	}
	return;
}
