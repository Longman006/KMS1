
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

typedef struct {
	double t;
	double V;
	double P;
	double H;
	double T;
}State;

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

	for(i = 0 ; i< v1->size ;i++){
		printf("%lf ",v1->xyz[i]);
	}
	printf("Vector \n");

}


/**
 * Zwraca nowy vector ktory jest sumą
 */
Vector* vectorSum(Vector** vvv,int size){
	double tab[] = {0,0,0};
	int iii,xyz = 0;
	tab[0] = 0;
	//po Vektorach
	for(iii = 0 ; iii<size ; iii++ ){
		//printf(" vector sum iii = %d\n",iii);
		//printVector(vvv[iii]);
		//Po współrzędnych
		for(xyz = 0 ; xyz<3 ; xyz++){
			//printf(" vector sum xyz = %d\n",xyz);
			//printf("tab[xyz] = %lf\n",tab[xyz]);
			//printf("vvv(xyz) = %lf\n",vvv[iii]->xyz[xyz]);
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
	Vector* temprr[3];
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
	tab[2] = params->paramsTab[AA]*sqrt(2.0/3.0);
	bb[2] = newVector(tab);

	printf("Vectory bb : \n");
	for(iii[0] = 0 ; iii[0]<3 ; iii[0]++){
		printVector(bb[iii[0]]);
	}

	for(iii[0]=0;iii[0]<n;iii[0]++){
		for(iii[1] = 0 ; iii[1]<n; iii[1]++){
			for(iii[2] = 0 ; iii[2]<n ; iii[2]++){

				/**
				 * Numer petli
				 */
				for(wsp=0,index = 0 ; wsp<3 ; wsp++){

					temprr[wsp] = vectorMultiply(bb[wsp],iii[wsp]-((n-1)/2) );
					index +=iii[wsp]*pow(n,wsp);
				}
				rr[index] = vectorSum(temprr, 3);
			}
		}
	}
	return rr;

}

Vector** initEnergy(Params* params){
	double N = pow(params->paramsTab[NN],3);
	double k = 8.31*pow(10,-3);
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
void writeToFileRawData(char* nazwaPliku,Vector** tablicaVect,int size){
	FILE *f = fopen(nazwaPliku, "w");
	int iii,jjj= 0 ;
	if (f == NULL){
	    printf("Error opening file!\n");
	    exit(1);
	}

	for(iii =  0 ; iii<size ; iii++){

		for(jjj = 0 ; jjj<3 ; jjj++){
			fprintf(f,"%lf ",tablicaVect[iii]->xyz[jjj]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
void writeToFile(char* nazwaPliku,Vector ** tablicaVect,int size){
	FILE *f = fopen(nazwaPliku, "w");
	int iii,jjj= 0 ;
	if (f == NULL){
	    printf("Error opening file!\n");
	    exit(1);
	}
	fprintf(f,"%d\n\n",size);

	for(iii =  0 ; iii<size ; iii++){
		fprintf(f,"Ar ");
		for(jjj = 0 ; jjj<3 ; jjj++){
			fprintf(f,"%lf ",tablicaVect[iii]->xyz[jjj]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
Vector* addVector(Vector* v1, Vector* v2){
	double tab[3];
	int iii;
	for(iii = 0 ; iii< 3 ; iii++){
		tab[iii] = v1->xyz[iii]+v2->xyz[iii];
	}
	return newVector(tab);
}
Vector* subtractVector(Vector* v1,Vector* v2){
	double tab[3];
	int iii;
	for(iii = 0 ; iii< 3 ; iii++){
		tab[iii] = v1->xyz[iii]-v2->xyz[iii];
	}
	return newVector(tab);
}
double vectorLengthSquare(Vector* vvv){
	int iii = 0 ;
	double retValue = 0 ;
	for(iii = 0 ; iii<3 ; iii++){
		retValue += vvv->xyz[iii]*vvv->xyz[iii];
	}
	return retValue;

}
double vectorLength(Vector* vvv){
	return sqrt(vectorLengthSquare(vvv));
}
void momentumCorrection(Vector** pedy,int size){
	int iii = 0 ;
	Vector* poprawka = vectorSum(pedy, size);
	poprawka = vectorMultiply(poprawka, -1.0/(double)(size));
	for(iii = 0 ; iii<size ; iii++){
		pedy[iii] = addVector(pedy[iii], poprawka);
	}
}
Vector** initForce(Vector** polozenia,Params* params){
	int N = pow(params->paramsTab[NN],3);
	double rij=0;
	double ri = 0 ;
	double epsilon = params->paramsTab[EE];
	double R = params->paramsTab[RR];
	double L = params->paramsTab[LL];
	double f = params->paramsTab[FF];
	double stalaP,stalaS = 0 ;
	int iii,jjj,wsp=0;
	double zeroVect[] = {0,0,0};

	Vector** forceP = malloc(sizeof(Vector*)*N);
	Vector** forcePjte = malloc(sizeof(Vector*)*N);
	Vector** forceS = malloc(sizeof(Vector*)*N);
	Vector** retValue = malloc(sizeof(Vector*)*N);

	if(retValue == NULL || forceP==NULL || forcePjte == NULL || forceS == NULL ){
		return NULL;
	}


	for(iii = 0 ; iii<N ; iii++){
		//printf("\n");
		for(jjj = 0 ;jjj <N ; jjj++){
			if(jjj==iii){
				forcePjte[jjj] = newVector(zeroVect);
				continue;
			}

			rij = vectorLength(subtractVector(polozenia[iii], polozenia[jjj]));
			//printf("rij = %lf\n",rij);
			stalaP = 12*epsilon/rij/rij*(pow(R/rij,12)-pow(R/rij,6));
			forcePjte[jjj] = vectorMultiply(subtractVector(polozenia[iii], polozenia[jjj]), stalaP);
		}

		ri = vectorLength(polozenia[iii]);
		forceP[iii]  = vectorSum(forcePjte, N);
		stalaS = f*(L-ri)/ri;
		if(ri<L)
			forceS[iii] =newVector(zeroVect);
		else
			forceS[iii] = vectorMultiply(polozenia[iii], stalaS);

		retValue[iii] = addVector(forceP[iii], forceS[iii]);

	}
	return retValue;

}
double initPot(Vector** polozenia,Params* params){
	int N = pow(params->paramsTab[NN],3);
	double rij=0;
	double ri = 0 ;
	double epsilon = params->paramsTab[EE];
	double R = params->paramsTab[RR];
	double L = params->paramsTab[LL];
	double f = params->paramsTab[FF];
	double stalaP,stalaS = 0 ;
	int iii,jjj;

	double potP = 0;
	double potPjte = 0;
	double potS = 0;
	double retValue =0;

	for(iii = 0 ; iii<N ; iii++){
		for(jjj = iii ;jjj <N ; jjj++){
			if(jjj==iii){
				continue;
			}
			rij = vectorLength(subtractVector(polozenia[iii], polozenia[jjj]));
			stalaP = epsilon*((pow(R/rij,12)-2*pow(R/rij,6)));
			potP+=stalaP;
		}

		ri = vectorLength(polozenia[iii]);
		stalaS = f*(ri-L)*(ri-L)/2;

		if(ri<L)
			potS+=0;
		else
			potS +=stalaS;

	}
	retValue+= potS;
	retValue+=potP;
	return retValue;

}

double calculateHamiltonian(Vector** pedy,double pot, Params* params){
	double Hamiltonian = 0 ;
	double m = params->paramsTab[MM];
	int N = pow(params->paramsTab[NN],3);
	int iii = 0 ;
	for(iii = 0 ; iii<N ; iii++){
		Hamiltonian+=vectorLengthSquare(pedy[iii])/2/m;
	}
	Hamiltonian+=pot;
	return Hamiltonian;
}
double calculatePressure(Vector** force,Params* param){
	double L = param->paramsTab[LL];
	double size = pow(param->paramsTab[NN],3);
	double retValue =  0;
	int iii = 0 ;

	for(iii = 0 ; iii<size ; iii++){
		retValue+=vectorLengthSquare(force[iii]);
	}
	return retValue;
}
double* calculateEnergy(Vector** pedy,Params* params){
	int iii = 0 ;
	int size = pow(params->paramsTab[NN],3);
	double m = params->paramsTab[MM];
	double* tabEn = malloc(sizeof(double)*size);

	for(iii = 0 ; iii < size ; iii++){
		tabEn[iii] = vectorLength(pedy[iii])/2/m;
	}
	return tabEn;
}
double sumTable(double* tab,int size){
	int iii = 0 ;
	double retValue = 0 ;
	for(iii = 0 ; iii<size ; iii++){
		retValue+=tab[iii];
	}
	return retValue;
}
double calculateTemperature(double* tabEn,Params* param){
	double k = 8.31*pow(10,-3);
	int size = pow(param->paramsTab[NN],3);
	int iii = 0 ;
	double retValue = 0 ;

	for(iii = 0 ; iii<size ; iii++){
		retValue+=tabEn[iii];
	}
	retValue= retValue*2/3/size/k;
	return retValue;
}
void newMomentumHalf(Vector** momentum,Vector** force,Params* params){
	int iii = 0 ;
	double tau = params->paramsTab[TAU];
	double m = params->paramsTab[MM];
	double N = pow(params->paramsTab[NN]);

	for(iii = 0 ; iii<N ; iii++){
		momentum[iii] = addVector(momentum[iii],vectorMultiply(force[iii], tau/2.0));
	}
}
void newPosition(Vector** polozenia,Vector** momentumHalf,Params* params){
	int iii = 0 ;
	double tau = params->paramsTab[TAU];
	double m = params->paramsTab[MM];
	double N = pow(params->paramsTab[NN]);

	for(iii =0 ; iii<N ;iii++){
		polozenia[iii] = vectorMultiply(momentumHalf[iii], tau/m);
	}
}
void newForce(Vector** polozenia,Vector** force,Params* params){
	//free(force);
	force = initForce(polozenia, params);
}
void newMomentum(Vector** momentum,Vector** force,Vector** momentumHalf,Params* params){
	int iii = 0 ;
	double tau = params->paramsTab[TAU];
	double m = params->paramsTab[MM];
	double N = pow(params->paramsTab[NN]);

	for(iii = 0 ; iii<N ; iii++){
		momentum[iii] = addVector(momentumHalf[iii], vectorMultiply(force[iii], tau/2.0));
	}
}
void nextStep(Vector** polozenia,Vector** pedy,Vector** sily,Vector** pedyHalf,Params* params){
	newMomentumHalf(pedy, sily, params);
	newPosition(polozenia, pedyHalf, params);
	newForce(polozenia, sily, params);
	newMomentum(pedy, sily, pedyHalf, params);
}
void calculateState(double t ,State* state,Vector** polozenia,Vector** pedy,Vector** sily,Vector** pedyHalf,Params* params ){
	state->t= t;
	state->V = initPot(polozenia, params);
	state->H = calculateHamiltonian(pedy, state->V, params);
	state->T = calculateTemperature(calculateEnergy(pedy, params), params);
	state->P = calculatePressure(sily, params);
}
void writeToFileState(State* state,char* nazwaPliku){
	FILE *f = fopen(nazwaPliku, "a");
		int iii,jjj= 0 ;
		if (f == NULL){
		    printf("Error opening file!\n");
		    exit(1);
		}
		fprintf(f,"%lf,%lf,%lf,%lf,%lf\n",state->t,state->H,state->V,state->T,state->P);
		fclose(f);
}
void simulate(State* state,Vector** polozenia,Vector** pedy,Vector** sily,Vector** pedyHalf,Params* params){
	int Sxyz = params->paramsTab[S_XYZ];
	int Sout = params->paramsTab[S_OUT];
	int S0 = params->paramsTab[S_O];
	int koniec = 100;
	int iii,jjj;
	for(iii = 0 ; iii<S0 ; iii++){
		nextStep(polozenia, pedy, sily, pedyHalf, params);
	}
	for(jjj = 0 ; jjj<koniec ; jjj++){
		if(jjj%Sout){
			writeToFileState(state, "plikStanowy.txt");
		}

	}
}

int main(void){

	int iii=0;
	int N =0;
	double ham = 0 ;
	Params* params = newParameters("PlikWejsciowy.txt");
	Vector** vvv;
	Vector** eeenergie;
	Vector** pedy;
	Vector** sily;
	double pot,pressure,temp;
	double* tabEn;

	N =  pow(params->paramsTab[NN] ,3);
	vvv = initPosition(params);
	pot = initPot(vvv,params);
	eeenergie = initEnergy(params);
	pedy = initMomentum(params, eeenergie);
	sily = initForce(vvv, params);
	printf("en\n");
	tabEn = calculateEnergy(pedy, params);
	printf("ham\n");
	ham = calculateHamiltonian(pedy, pot, params);
	printf("press\n");
	pressure = calculatePressure(sily, params);
	printf("temp\n");
	temp = calculateTemperature(tabEn, params);

	printf("Ham = %lf\nPress = %lf\n temp = %lf\n",ham,pressure,temp);
	writeToFile("plikPol.txt", vvv, pow(params->paramsTab[NN],3));
	writeToFileRawData("pedy.txt",pedy,N );
	return 1;
}
