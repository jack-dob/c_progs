#ifndef SSIO_HINCLUDED
#define SSIO_HINCLUDED

/*
 ** ssio.h -- DCR 98-09-16
 ** ======
 ** Header file specific to Solar System data I/O routines.
 */
#ifndef STDIO_INCLUDED
	#include <stdio.h>
	#define STDIO_INCLUCED
#endif
#ifndef STRING_INCLUDED
	#include <string.h>
	#define STRING_INCLUDED
#endif
#ifndef STDLIB_INCLUDED
	#include <stdlib.h>
	#define STDLIB_INCLUDED
#endif
#ifndef PARAM_INCLUDED
	#include <sys/param.h>	/*for MAXPATHLEN */
	#define PARAM_INCLUDED
#endif
#ifndef RPC_INCLUDED
	#include <rpc/rpc.h>	/* for XDR routines */
	#define RPC_INCLUDED
#endif
#ifndef ASSERT_INCLUDED
	#include <assert.h>
	#define ASSERT_INCLUDED
#endif
#ifndef MATH_INCLUDED
	#include <math.h>
	#define MATH_INCLUDED
#endif

#define FLAG_TO_DELETE(a) (a).org_idx = -1


#ifndef VERBOSE
#define VERBOSE 0
#endif /*VERBOSE*/

#ifndef DEBUG
	#define DEBUG 0
#endif

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

#ifndef N_DIM
#define N_DIM 3
#endif

#define MALLOC_ERROR "Error: %s could not initialise array %s\n"

/*
 ** Following structures intended for use with xdr I/O routines. Note that
 ** internal representations may differ from format written to disk. As a
 ** result, for now use hardwired sizes for seeking. Is there an xdr routine
 ** to figure this out automatically?...
 */

typedef struct ss_head {
	double	time;
	int		n_data;
	int		pad;	/* unused: pad to 8-byte boundary */
	} SSHEAD;

#define SSHEAD_SIZE 16	/* XDR assumes 8-byte doubles and 4-byte ints */

typedef struct ss_data {
	double	mass;
	double	radius;
	double	pos[N_DIM];
	double	vel[N_DIM];
	double	spin[N_DIM];
	int		color;
	int		org_idx;
	} SSDATA;

#define SSDATA_SIZE 96	/* XDR assumes 8-byte doubles and 4-byte ints */

typedef struct ssio {
	FILE *fp;
	XDR xdrs;
	} SSIO;

#define SSIO_READ	0
#define SSIO_WRITE	1
#define SSIO_UPDATE	2

#define SS_EXT ".ss"



int writeout_ss(char* filename, SSHEAD *head, SSDATA *data);
SSDATA* readin_ss(char* filename, SSHEAD* head);
int ssioNewExt(const char *infile, const char *inext,
			   char *outfile, const char *outext);
int ssioOpen(const char *filename, SSIO *ssio, u_int mode);
int ssioHead(SSIO *ssio, SSHEAD *head);
int ssioData(SSIO *ssio, SSDATA *data);
int ssioClose(SSIO *ssio);
int ssioSetPos(SSIO *ssio, u_int pos);
void ssioRewind(SSIO *ssio);


int writeout_ss(char* filename, SSHEAD *head, SSDATA *data){
	SSIO ssio;
	int i;
	if(ssioOpen(filename, &ssio, SSIO_WRITE)){
		fprintf(stderr, "Error: Could not open file for writing: %s", filename);
		return(1);
	}
	
	if(ssioHead(&ssio, head)){
		fprintf(stderr, "Error: Could not write header data for some reason...");
		return(1);
	}
	
	for(i=0;i<head->n_data;i++){
		if(ssioData(&ssio, &data[i])){
			fprintf(stderr, "Error: Could not write particle data for some reason...");
			return(1);
		}
	}
	
	if(ssioClose(&ssio)){
		fprintf(stderr, "Error: Could not close the file for some reason...");
		return(1);
	}
	return (0);
}

void readin_sshead(char* filename, SSHEAD *head){
	SSIO ssio;
	
	if(VERBOSE)	printf("Opening file: %s\n", filename);
	
	if(ssioOpen(filename, &ssio, SSIO_READ)){
		fprintf(stderr, "Error: readin_ss() could not open file %s...\n", filename);
		return;
	}
	
	if(ssioHead(&ssio, head)){
		fprintf(stderr, "Error: readin_ss() could read header data for some reason...");
		return;
	}
	
	if(VERBOSE) printf("Grabbed header data: Time = %lf, N = %i, pad = %i\n",
			head->time, head->n_data, head->pad);
	
	if(VERBOSE) printf("Header read, closing file and returning...\n");
	if(ssioClose(&ssio)){
		fprintf(stderr, "Error: readin_ss() could not close file %s for some reason...", filename);
		return;
	}
	return;
}

SSDATA* readin_ss(char* filename, SSHEAD *head){
	/*
	WARNING!!!
	This contains a malloc statement, you will need to free
	your array after using this!!!!
	*/
	SSIO ssio;
	SSDATA *temp;
	int i;
	
	if(VERBOSE)	printf("Opening file: %s\n", filename);
	
	if(ssioOpen(filename, &ssio, SSIO_READ)){
		fprintf(stderr, "Error: readin_ss() could not open file %s...\n", filename);
		return(NULL);
	}
	
	if(ssioHead(&ssio, head)){
		fprintf(stderr, "Error: readin_ss() could read header data for some reason...");
		return(NULL);
	}
	
	if(VERBOSE) printf("Grabbed header data: Time = %lf, N = %i, pad = %i\n",
			head->time, head->n_data, head->pad);
	
	if(VERBOSE) printf("Allocating memory...\n");
	temp = (SSDATA*) malloc( head->n_data*sizeof(SSDATA) );

	if(temp==NULL){
		fprintf(stderr, "Error: readin_ss() could not allocate memory...\n");
		return(NULL);
	}
	if(VERBOSE) printf("Memory allocated...\n");
	
	if(VERBOSE) printf("Reading in data...\n");
	for(i=0; i<head->n_data ; i++){
		if(ssioData(&ssio, &temp[i])){
			fprintf(stderr, "Error: readin_ss() could not read particle data for some reason...");
			return(NULL);
		}
		if(VERBOSE) printf("%i\r", temp[i].org_idx);
	}
	if(VERBOSE) printf("Data read, closing file and returning...\n");
	if(ssioClose(&ssio)){
		fprintf(stderr, "Error: readin_ss() could not close file %s for some reason...", filename);
		return(NULL);
	}
	return(temp);
}


static const char *
ssioBasename(const char *path)
{
	char *p;

	assert(path != NULL);
	p = strrchr(path,'/');
	if (p) return p + 1;
	else return path;
	}

int
ssioNewExt(const char *infile,const char *inext,
		   char *outfile,const char *outext)
{
	const char *basename;
	char *c;
	size_t n;

	assert(infile != NULL && inext != NULL && outfile != NULL && outext != NULL);
	basename = ssioBasename(infile);
	if ((c = strrchr(basename,'.')) && strstr(c,inext))
		n = c - basename;
	else
		n = strlen(basename);
	if (n + strlen(outext) >= (size_t) MAXPATHLEN)
		return 1;
	(void) strncpy(outfile,basename,n); /* not null terminated */
	(void) strcpy(outfile + n,outext);
	return 0;
	}

int
ssioOpen(const char *filename,SSIO *ssio,const u_int mode)
{
	const char type[][3] = {"r","w","r+"};
	const enum xdr_op op[] = {XDR_DECODE,XDR_ENCODE,XDR_ENCODE};

	assert(filename != NULL && ssio != NULL);
	assert(mode == SSIO_READ || mode == SSIO_WRITE || mode == SSIO_UPDATE);
	if (!(ssio->fp = fopen(filename,type[mode]))){
		return 1;
	}
	
	xdrstdio_create(&ssio->xdrs,ssio->fp,op[mode]);
	return 0;
	}

int
ssioHead(SSIO *ssio,SSHEAD *head)
{
	assert(ssio != NULL && head != NULL);
	if (!xdr_double(&ssio->xdrs,&head->time)) return 1;
	if (!xdr_int(&ssio->xdrs,&head->n_data)) return 1;
	if (!xdr_int(&ssio->xdrs,&head->pad)) return 1;
	return 0;
	}

int
ssioData(SSIO *ssio,SSDATA *data)
{
	int i;

	assert(ssio != NULL && data != NULL);
	if (!xdr_double(&ssio->xdrs,&data->mass)) return 1;
	if (!xdr_double(&ssio->xdrs,&data->radius)) return 1;
	for (i=0;i<N_DIM;i++)
		if (!xdr_double(&ssio->xdrs,&data->pos[i])) return 1;
	for (i=0;i<N_DIM;i++)
		if (!xdr_double(&ssio->xdrs,&data->vel[i])) return 1;
	for (i=0;i<N_DIM;i++)
		if (!xdr_double(&ssio->xdrs,&data->spin[i])) return 1;
	if (!xdr_int(&ssio->xdrs,&data->color)) return 1;
	if (!xdr_int(&ssio->xdrs,&data->org_idx)) return 1;
	return 0;
	}

int
ssioClose(SSIO *ssio)
{
	assert(ssio != NULL);
	xdr_destroy(&ssio->xdrs);
	if (fclose(ssio->fp)) return 1;
	return 0;
	}

int
ssioSetPos(SSIO *ssio,const u_int pos)
{
	assert(ssio != NULL);
	if (!xdr_setpos(&ssio->xdrs,pos)) return 1;
	return 0;
	}

void
ssioRewind(SSIO *ssio)
{
	assert(ssio != NULL);
	rewind(ssio->fp);
	}
	
#endif /*SSHEAD_INCLUDED*/
