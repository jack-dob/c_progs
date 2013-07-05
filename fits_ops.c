#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "fitsio.h" //must use compiler flag '-l cfitsio'
//	Make a load of functions to read in .fits files


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
	int i, j, k, pnum;
	fitsfile *fitsptr;
	char *filename;
	char out_file[128];
	int status=0;
	int bitpix, bpix_size; //the type of data stored in the image
	int max_elem;
	int naxis=8;
	long naxes[3];
	long fpixel[3]={1, 1, 1};
	//long fpixel = 1;
	long nelements;
	double *data;
	double *iscale;
	int anynul;
	char buf[100];
	
	filename = argv[1]; //first argument is the file name
	
	//get_out_file(out_file, filename, "flat_");
	//return; //for testing
	if (fits_open_file(&fitsptr, filename, READONLY, &status) )
		printerror(status);
	if(fits_read_key(fitsptr, TINT, "NAXIS", &naxis, NULL, &status))
		printerror(status);
	//printf("naxis: %d, %s\n", naxis, buf);
	for(i=1;i<=naxis;i++){
		sprintf(buf, "NAXIS%d", i);
		if(fits_read_key(fitsptr, TINT, buf, &(naxes[i-1]), NULL, &status))
			printerror(status);
		//printf("naxes[%d]: %d\n", i-1, naxes[i-1]);
	}
	
	//if(naxis != maxdim){
	//	printf("Error: reading a maximum of %d axes, found %d axes. Too many, exiting...\n", maxdim, naxis);
	//	exit(0);
	//}
	//bpix_size = sizeof_bitpix(bitpix);
	//allocate container array
	nelements=1;
	for(i=0;i<naxis;i++){
		nelements *= naxes[i];
	}
	//printf("nelements: %d\n", nelements);
	if(argc>=3){
		max_elem=atoi(argv[2]);
	}
	else{
		max_elem = nelements;
	}
	//printf("max_elem: %d\n", max_elem);
	//nelements=100;
	data = (double*) malloc(max_elem*sizeof(double));
	pnum=0;
	//puts("Starting loop");
	k=0;
	//still breaks if do not read the whole file at once!!
	while((int)nelements > 0){
		//puts("Blaap");
		if(fits_read_pix(fitsptr, TDOUBLE, fpixel, max_elem, NULL, data, &anynul, &status))
			printerror(status);
		// do stuff
		nelements-=max_elem;
		//puts("Bleep");
		vec_wrap(fpixel, naxes, naxis, max_elem);
		pnum+=max_elem;
		//printf("pnum: %d, nelements: %d\n", pnum, nelements);
		//puts("Bloop");
		if(k==2)
			exit(0);
		k++;		
	}
	
	/*
	// FIND MINIMUM AND MAXIMUM
	//printf("%s:\t", filename);
	min_max(data, naxes);
	*/
	
	
	// FLATTEN CUBE
	if(naxis==3){
		flatten_cube(data, naxes);
	}
	iscale = (double*) malloc(3*(int)sizeof(double));
	iscale[0] = 10;
	iscale[1] = 10;
	iscale[2] = 0.1; //placeholder
	get_out_file(out_file, filename, "flat_");
	write_ds9_cube(out_file, naxes, iscale, data);
	free(iscale);
	
	
	//do stuff here with the whole array
	
	// RADIAL SMOOTH
	//get_out_file(out_file, filename, "rad_");
	iscale = (double*) malloc(3*(int)sizeof(double));
	iscale[0] = 10;
	iscale[1] = 10;
	iscale[2] = 0.1; //placeholder
	rad_smooth(data, naxes, iscale);
	get_out_file(out_file, filename, "rad_");
	write_ds9_cube(out_file, naxes, iscale, data);
	free(iscale);
	
	//close files
	fits_close_file(fitsptr, &status);
	//free arrays
	free(data);
	
	return(0);
}

void rad_smooth(double *data, long *naxes, double *iscale){
	int i, j, k;
	int nelements=1, nbin=15, bin;
	double x, y, z, r;
	double binstart = 0.7, binend = 3.9;
	double binstep, area, test1;
	double xoff = -5, yoff = -5;
	double *bin_array;
	
	bin_array = (double*) malloc( nbin*(int)sizeof(double) );
	
	binstep = (binend-binstart)/(double)nbin;
	
	printf("binstep: %g, binstart: %g, binend: %g\n", binstep, binstart, binend);
	
	for(i=0;i<3;i++){
		nelements*=(int)naxes[i];
		printf("naxis[%d]: %d, iscale[%d]: %g, nelements: %d\n", i, naxes[i], i, iscale[i], nelements);
	}
	puts("creating bin array");
	for(i=0;i<nelements;i++){
		printf("progrsss: %0.2lf\r", (double)i/(double)nelements);
		x = (double)(i%naxes[0]);
		//printf("xi: %g", x);
		x = (iscale[0]*x/(double)naxes[0]) - iscale[0]/2.0; //should be in AU from centre now
		y = (double)floor(((double)i)/(double)naxes[0]);
		//printf("\tyi: %g\n", y);
		y = (iscale[1]*y/(double)naxes[1]) - iscale[1]/2.0; //should be in AU from centre now
		//z = floor((float)i/(naxes[0]*naxes[1]));
		//z = (z/iscale[2]) - iscale[2]/2.0;
		z = 0.0; //can assume this for thin disk
		r = sqrt(x*x + y*y + z*z);
		//test1 = floor((r-binstart)/binstep);
		//printf("test1: %g\n", test1);
		bin = (int)floor((r - binstart)/binstep);
		//printf("i: %d, x: %g, y: %g, z: %g, r: %g, bin: %d, test1: %g, binstep: %g\n",\
		//i, x, y, z, r, bin, test1, binstep);
		//may want to check all this
		//if(i==2)
		//	exit(0); //for testing
		if((bin >= 0) && (bin < nbin)){
			bin_array[bin] += data[i]; //add mass to the correct bin
		}
	}
	//for(i=0;i<nbin;i++){
	//	printf("bin_array[%d]: %G\n", i, bin_array[i]);
	//}
	puts("creating picture");
	for(i=0;i<nelements;i++){
		printf("progrsss: %0.2lf\r", (double)i/(double)nelements);
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
		//printf("i: %d, x: %g, y: %g, z: %g, r: %g, bin: %d, binstep: %g\n",\
		//i, x, y, z, r, bin, binstep);
		if((bin >= 0) && (bin < nbin)){
			data[i] = bin_array[bin]/area;
		}
		else
			data[i] = 0.0;
	}
	free(bin_array);
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
