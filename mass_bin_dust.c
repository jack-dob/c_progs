#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void read_dustopac(int *n_dust, double **dust_rad, double **dust_factor);
int str_split(char *bigstr, char ***splitstr, char delim);
void parse_tags(char **lines, double **dust_rad, double **dust_factor, int i, int j); 

int main(int argc, char **argv){
	int n_dust, i;
	double *dust_rad, *dust_factor;
	//n_dust = number of dust species
	//dust_rad = list of radii for each dust species (used to work out mass)
	//dust_factor = list of the proportion of the total mass (in a given mass bin) that this species makes up
	read_dustopac(&n_dust, &dust_rad, &dust_factor);
	
	printf("Number of dust species: %d\n", n_dust);
	for(i=0;i<n_dust;i++){
		printf("Species radius: %G, Species factor: %G\n", dust_rad[i], dust_factor[i]);
	}
	return(0);
}

void read_dustopac(int *n_dust, double **dust_rad, double **dust_factor){
	//n_dust = number of dust species
	//dust_rad = array of radii of dust species
	//dust_factor = array of mass proportion in a dust bin of a given dust species
	FILE* dp;
	int max_char = 80, i, j; //maximum length of a line
	char *f_data, **lines;
	char tags[2];
	char ch;
	int f_len=0, n_lines; //length of file

	dp = fopen("dustopac.inp", "r");
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
	j=0;
	for(i=5;i<n_lines;i++){
		if( !((i-5)%4) ){ //chooses only the species names: <tag>_<size>
			//printf("i: %d\n", i);
			//printf("%s\n", lines[i]);
			parse_tags(lines, dust_rad, dust_factor, i, j);
			j++;
		}
	}
	
	puts("Done.");
	
	//free malloc'd stuff
	free(f_data);
	free(lines);
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

void parse_tags(char **lines, double **dust_rad, double **dust_factor, int i, int j){
	char **parts;
	int n_parts;
	//puts("HERE 1");
	n_parts = str_split(lines[i], &parts, '_'); //n_parts should be 2
	//puts("HERE 2");
	(*dust_rad)[j] = strtod(parts[1], NULL);
	//puts("HERE 3");
	//printf("%s\t%s\n", parts[0], parts[1]);
	if(!strcmp(parts[0],"AC1")){ //amorphous carbon makes up ~20% of dust
		//puts("BLOB");
		(*dust_factor)[j] = 0.2;
	}
	else if(!strcmp(parts[0],"sil")){ //silicates make up ~80% of dust
		//puts("EUGH");
		(*dust_factor)[j] = 0.8;
	}
	else{ //if none of the above, assume there is only one type of dust
		//puts("HURG");
		(*dust_factor)[j] = 1.0;
	}
	return;
}











