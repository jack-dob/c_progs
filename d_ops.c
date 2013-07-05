#include <stdio.h>
#include <stdlib.h>

#define BUFF_SIZE 128

typedef struct grid_info{
	char *filename;
	int iformat;
	int blank;
	int coord_sys;
	int extra_info;
	int incl[3];
	int num[3];
	double ***walls;
	char *remainder;
}GRID_INFO;

typedef struct dust_info{
	char *filename
	int iformat, nrcells, nrspec;
	double *density;
} DUST_INFO;

int* split_str(char *str, char delim);
void read_amr_grid(GRID_INFO *grid_info);

int main(int argc, char **argv){
	GRID_INFO grid_info;
	DUST_INFO dust_info;
	
	if(argc>1){
		dust_info.filename = argv[1];
	}
	else{
		puts("No dust_density.inp file given");
		exit(1);
	}
	if(argc>2){
		grid_info.filename = argv[2];
	}
	else{
		puts("No amr_grid.inp file given");
		exit(1);
	}
	
	read_amr_grid(&grid_info);
	
	return(0);
}

void read_amr_grid(GRID_INFO *grid_info){
	FILE *f;
	char buff[BUFF_SIZE];
	int *split_pos;
	
	f = fopen(grid_info->filename, "r");
	
	fgets(buff, BUFF_SIZE, f); //read line
	grid_info->iformat = atoi(buff);
	fgets(buff, BUFF_SIZE, f); //read line
	grid_info->blank = atoi(buff);
	fgets(buff, BUFF_SIZE, f); //read line
	grid_info->coord_sys = atoi(buff);
	fgets(buff, BUFF_SIZE, f); //read line
	grid_info->extra_info = atoi(buff);
	fgets(buff, BUFF_SIZE, f); //read line, this consists of three numbers

	

}


int* split_str(char *str, char delim){
	//returns an integer array of starting points for each string
	//replaces the delimiter with '\0'
	
	int i, j, n;
	int *pos;
	
	for(n=0;str[n]!='\0';n++); //finds length of string
	for(i=0;i<n;i++){
		if(str[i]==delim){
			j++;
		}
	}
	pos = (int*) malloc( j*sizeof(int));
	j=0;
	for(i=0;i<n;i++){
		if(str[i]==delim){
			str[i]='\0';
			pos[j] = i+1;
			j++;
		}
	}
	return(pos);
}






