#include <stdio.h>

#define TYPE THAT
#define SORT a

#define STR(x) #x

#define A_MEMBER(kind, thing) (((kind*)thing)->a)

#define a_MEMBER(thing) (((TYPE*)thing)->a)

#define MEMBER(thing) (((TYPE*)thing)->SORT)

typedef struct this{
	int a;
	int b;
}THIS;

typedef struct that{
	double fill;
	int a;
	double b;
}THAT;


void printstruct_a(void *blob, char* type);

int main(int *argc, char **argv){
	THIS this1 = {1, 2};
	THAT that1 = {100, 3, 4};
	printf("%s", STR(this1));
	printstruct_a(&this1, "THIS");
	//#define TYPE THAT
	printstruct_a(&that1, "THAT");
	return(0);
}

void printstruct_a(void *blob, char* kind){
	printf("%d\t", MEMBER(blob) );
	printf("%lf\n", MEMBER(blob) );
	return;
}