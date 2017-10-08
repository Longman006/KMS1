
#include<stdio.h>
#include<stdlib.h>
#define NN 0
#define MM 1
#define EE 2
#define RR 3
#define FF 4
#define LL 5
#define AA 6
#define T_0 7
#define TAU 8
#define S_O 9
#define S_D 10
#define S_OUT 11
#define S_XYZ 12

typedef struct  {
	double* paramsTab;
	int size;
}Params;
typedef struct {
	double* xyz;
	int size;
}Vector;



Vector* newVector(double* tab){
	int i =0;
	Vector* retValue = malloc(sizeof( Vector));

	if(retValue == NULL)
	return NULL;

	retValue->size = 3;
	retValue->xyz = malloc(sizeof(double)*retValue->size);
	if (retValue->xyz == NULL) {
        free (retValue);
        return NULL;
    }
    for( i = 0 ; i< retValue->size ; i++){
		retValue->xyz[i] = tab[i];
	}

	return retValue;
}
	Vector* vectorMultiply(Vector* v1,double a){
	double tab[3];
	int i=0;
	for(i = 0 ; i < v1->size ; i++){
		tab[i] = v1->xyz[i]*a;
	}
	return newVector(tab);
}
void printVector(Vector* v1){
	int i =0;
	printf("Vector \n");
	for(i = 0 ; i< v1->size ;i++){
		printf("%lf",v1->xyz[i]);
	}
}
 Params* newParameters(char* plikParams){

	Params* retValue = malloc(sizeof( Params));

	if(retValue == NULL)
		return NULL;

	retValue->paramsTab = malloc(sizeof(double)*13);
	if (retValue->paramsTab == NULL) {
        free (retValue);
        return NULL;
    }

	FILE *fp;
	fp=fopen(plikParams, "r");


	if (fp == NULL) {
		fprintf(stderr, "Can't open input file !\n");
		exit(1);
	}

	retValue->size = 13;
	int i=0;
	while (fscanf(fp, "%lf", &retValue->paramsTab[i]) == 1){
		printf("%d %lf \n",i,&retValue->paramsTab[i]);
		i++;
	}

	fclose(fp);
	return retValue;
}


int main(void){

	double tab[] = {1,2,3};
	Params* params = newParameters("PlikWejsciowy.txt");
	Vector* vector = newVector(tab);
	Vector* v2 = newVector(tab);
	printVector(v2);
	printVector(vectorMultiply(v2,2));

	printf("Lab1 \n");

	return 1;
}
