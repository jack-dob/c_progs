#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#ifndef MASS_BIN_DUST_INCLUDED

typedef struct opacity_dat{
	char *tag;
	double rad;
	double factor;
}OPAC_DAT;

void read_dustopac(FILE *dp, int *n_dust, double **dust_rad, double **dust_factor);
int str_split(char *bigstr, char ***splitstr, char delim);
void parse_tags(char **lines, OPAC_DAT *opac_dat, int i, int j);
void sort_opac_dat(OPAC_DAT *opac_dat, int n_dust);

void read_dustopac(FILE *dp, int *n_dust, double **dust_rad, double **dust_factor){
	//dp = dustopac.inp file handle
	//n_dust = number of dust species
	//dust_rad = array of radii of dust species
	//dust_factor = array of mass proportion in a dust bin of a given dust species
	int max_char = 80, i, j; //maximum length of a line
	char *f_data, **lines;
	char tags[2];
	char ch;
	int f_len=0, n_lines; //length of file
	OPAC_DAT *opac_dat;
	
	//dp = fopen("dustopac.inp", "r");
	for(ch=fgetc(dp); ch != EOF; ch=fgetc(dp)){
		f_len++;
	} //read to the end of the file
	
	f_data = (char*) malloc( f_len*sizeof(char)); //allocate space for the file

	if(f_data ==NULL){
		puts("Could not allocate memory to read in dustopac.inp file, exiting...");
		exit(1);
	}

	rewind(dp); //go back to beginning of file
	fread(f_data, sizeof(char), f_len, dp); // read in the file to f_data
	close(dp); //don't need file any more

	n_lines = str_split(f_data, &lines, '\n');
	*n_dust = atoi(lines[1]); //convert string to integer
	(*dust_rad) = (double*) malloc(*n_dust*sizeof(double));
	(*dust_factor) = (double*) malloc(*n_dust*sizeof(double));
	opac_dat = (OPAC_DAT*) malloc(*n_dust*sizeof(OPAC_DAT));
	
	j=0;
	for(i=5;i<n_lines;i++){
		if( !((i-5)%4) ){ //chooses only the species names: <tag>_<size>
			//printf("i: %d\n", i);
			//printf("%s\n", lines[i]);
			parse_tags(lines, opac_dat, i, j);
			j++;
		}
	}
	//sort opac_dat so that the different sizes are next to eachother, irrespective of tag
	sort_opac_dat(opac_dat, *n_dust);
	//assign correct stuff to arrays
	for(i=0;i<*n_dust;i++){
		(*dust_rad)[i] = opac_dat[i].rad;
		(*dust_factor)[i] = opac_dat[i].factor;
	}
	puts("Done.");
	
	//free malloc'd stuff
	free(f_data);
	free(lines);
	free(opac_dat);
	return;
}

int str_split(char *bigstr, char ***splitstr, char delim){
	//splits bigstr into chunks separated by delim and stores i splitstr array
	int n_chunks=0, max_len=0, test=0, i, j, len=0;
	char ch;
	
	for(i=0;bigstr[i]!='\0';i++){
		if( bigstr[i] == delim){
			n_chunks++;
		}
		len++;
	}
	
	(*splitstr) = (char**) malloc( n_chunks*sizeof(char*)); //make space for pointers
	(*splitstr)[0] = bigstr; //first string starts at beginning
	j=1;
	for(i=0; i<(len-1); i++){
		if( bigstr[i]==delim){ //if we are at a delimiter
			bigstr[i] = '\0'; //split the string
			(*splitstr)[j] = &(bigstr[i+1]); //start new string at the next line
			j++;
		}
	}
	return(n_chunks);
}

void parse_tags(char **lines, OPAC_DAT *opac_dat, int i, int j){
	char **parts;
	int n_parts;
	//puts("HERE 1");
	n_parts = str_split(lines[i], &parts, '_'); //n_parts should be 2
	//puts("HERE 2");
	opac_dat[j].rad = strtod(parts[1], NULL);
	opac_dat[j].tag = parts[0];
	//puts("HERE 3");
	//printf("%s\t%s\n", parts[0], parts[1]);
	if(!strcmp(opac_dat[j].tag,"AC1")){ //amorphous carbon makes up ~20% of dust
		//puts("BLOB");
		opac_dat[j].factor = 0.2;
	}
	else if(!strcmp(opac_dat[j].tag,"sil")){ //silicates make up ~80% of dust
		//puts("EUGH");
		opac_dat[j].factor = 0.8;
	}
	else{ //if none of the above, assume there is only one type of dust
		//puts("HURG");
		opac_dat[j].factor = 1.0;
	}
	return;
}

void sort_opac_dat(OPAC_DAT *opac_dat, int n_dust){
	//uses quicksort algorithm
	int i, j, pivot, temp_i;
	OPAC_DAT temp;
	
	if(n_dust<2){
		return;
	}

	pivot = n_dust/2; //start with middle number
	for(i=0;i<n_dust;i++){
		if( (opac_dat[i].rad == opac_dat[pivot].rad) && \
			(opac_dat[i].factor == opac_dat[pivot].factor) && \
			(opac_dat[i].tag == opac_dat[pivot].tag) ){
			continue; //don't sort pivot against it's self	
		}
		if( (opac_dat[i].rad < opac_dat[pivot].rad) && (i > pivot) ){ //if less than pivot
			//swap data and pivot
			temp = opac_dat[i];
			opac_dat[i] = opac_dat[pivot];
			opac_dat[pivot] = temp;
			//remember old pivot position
			temp_i = pivot;
			// pivot is now in at current position
			pivot = i;
			// must check stuff after old pivot position again
			i = temp_i;
		}
		if( (opac_dat[i].rad > opac_dat[pivot].rad) && (i < pivot) ){
			//swap data and pivot
			temp = opac_dat[i];
			opac_dat[i] = opac_dat[pivot];
			opac_dat[pivot] = temp;
			//pivot is at current position
			pivot = i;
		}
	}
	j = pivot+1;
	if(pivot>1){
		sort_opac_dat(opac_dat, pivot);
	}
	if(j<(n_dust-1)){
		sort_opac_dat(opac_dat+j, n_dust-j);
	}
	return;
}

#endif //MASS_BIN_DUST_INCLUDED









