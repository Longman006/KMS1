
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>

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

double randomDouble(double max){
	static int check = 0;
	double retValue = 0;
	if(!check){
		srand(time(NULL));
		check = 1;
	}
	return (double)rand()/(double)RAND_MAX*max;
}
double randomSign(void){
	if(randomDouble(2)>1)
		return 1;
	return -1;
}

/**
 * nie ingeruje w tab, alokuje pamięć i przyrównuje
 */
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

/**
 * Mnozy v1 przez a i zwraca nowy wektor (alokuje)
 */
Vector* vectorMultiply(Vector* v1,double a){
	double tab[3];
	int i=0;
	for(i = 0 ; i < v1->size ; i++){
		tab[i] = v1->xyz[i]*a;
	}
	return newVector(tab);
}
/**
 * wypisuje współrzędne
 */
void printVector(Vector* v1){
	int i =0;
	printf("Vector \n");
	for(i = 0 ; i< v1->size ;i++){
		printf("%lf",v1->xyz[i]);
	}

}


/**
 * Zwraca nowy vector ktory jest sumą
 */
Vector* vectorSum(Vector** vvv,int size){
	double tab[] = {0,0,0};
	int iii,xyz = 0;

	//po Vektorach
	for(iii = 0 ; iii<size ; iii++ ){
		//Po współrzędnych
		for(xyz = 0 ; xyz<3 ; xyz++){
			tab[xyz]+= vvv[iii]->xyz[xyz];
		}
	}

	return newVector(tab);
}


/**
 * Wczytuje parametry alokuje, nazwa pliku w argumencie
 */
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
		printf("%d %lf \n",i,retValue->paramsTab[i]);
		i++;
	}

	fclose(fp);
	return retValue;
}

/**
 * Ustawiamy polozenia poczatkowe
 */
Vector** initPosition(Params* params){

	double tab[] = {0,0,0};
	Vector* temprr;
	double n = params->paramsTab[NN];
	double N = pow(n,3);
	//Vector* bb0,bb1,bb2;
	Vector* bb[3];
	Vector** rr;
	int wsp;
	int iii[3] = {0};
	int index=0;

	rr = malloc(sizeof(Vector*)*N);
	if(rr == NULL)
		return NULL;
	/**
	 * b_0
	 */
	tab[0] = params->paramsTab[AA];
	bb[0] = newVector(tab);

	/**
	 * b_1
	 */
	tab[0] = tab[0]/2;
	tab[1] = sqrt(3)*tab[0];
	bb[1] = newVector(tab);

	/**
	 * b_2
	 */
	tab[1] = tab[1]/2;
	tab[2] = params->paramsTab[AA]*sqrt(2/3);
	bb[2] = newVector(tab);

	for(iii[0]=0;iii[0]<params->paramsTab[NN];iii[0]++){
		for(iii[1] = 0 ; iii[1]<params->paramsTab[NN]; iii[1]++){
			for(iii[2] = 0 ; iii[2]<params->paramsTab[NN] ; iii[2]++){

				/**
				 * Numer petli
				 */
				for(wsp=0,index = 0 ; wsp<3 ; wsp++){

					vectorMultiply(bb[wsp],iii[wsp]-((n-1)/2) );
					index +=iii[wsp]*pow(n,wsp);
				}
				rr[index] = vectorSum(bb, 3);
			}
		}
	}
	return rr;

}

Vector** initEnergy(Params* params){
	double N = pow(params->paramsTab[NN],3);
	double k = 1.3806*pow(10,-2);
	double T0 = params->paramsTab[T_0];
	double lambda[3];
	double stala = -1.0/2.0*k*T0;
	int iii,wsp=0;
	Vector** energie;

	printf("stala %lf \n",stala);

	energie = malloc(sizeof(Vector*)*N);
	if(energie == NULL)
		return NULL;

	for(iii = 0 ; iii < N ; iii++){
		for(wsp = 0 ; wsp<3 ; wsp++){
			lambda[wsp] =log(randomDouble(1));
			printf("lambda[%d],%lf",wsp,lambda[wsp]);
		}
		energie[iii] = vectorMultiply(newVector(lambda),stala);
	}
	return energie;
}

Vector** initMomentum(Params* params,Vector** energies){
	double mmm = params->paramsTab[MM];
	double nnn = pow(params->paramsTab[NN],3);
	double tab[3] = {0};
	int wsp = 0;
	int iii = 0 ;
	Vector** retValue;

	retValue = malloc(sizeof(Vector*)*nnn);
	if(retValue== NULL)
		return NULL;

	for(iii = 0 ; iii < nnn ; iii++ ){
		for(wsp = 0 ; wsp <3;wsp++){
			tab[wsp] = randomSign()* sqrt(2*mmm*energies[iii]->xyz[wsp]);
		}
		retValue[iii] = newVector(tab);
	}
	return retValue;
}

int main(void){

	int iii=0;
	Params* params = newParameters("PlikWejsciowy.txt");
	Vector** vvv;
	Vector** eeenergie;
	Vector** pedy;

	vvv = initPosition(params);
	eeenergie = initEnergy(params);
	pedy = initMomentum(params, eeenergie);

	for(iii=0; iii<params->paramsTab[NN] ; iii++ ){
		printVector(vvv[iii]);
		printVector(eeenergie[iii]);
		printVector(pedy[iii]);
	}

	return 1;
}
