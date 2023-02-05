#include<math.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>


#define D 2 
#define L 100
#define Q 2
#define extf 0
#define PI 3.14159265

const int stat=10000;

struct sito {
	int s;
	int f;
	int flag;
};
struct sito lattice[L][L];

void initialize(){
	int i, j;
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
			lattice[i][j].s=0;
			lattice[i][j].flag=0;
		}
	}
	return;
}

int spinew (int i, int j){
	int z;
	int sold;
	int snew;
	sold=lattice[i][j].s;
	z=(int) rand()%Q;
	while (z==sold){
		z= (int) rand()%Q;
	}
	snew=z;
	return snew;
}
void geometry (int *n, int *m){
 int i;
	for (i=0;i<L;i++){
		n[i]=i+1;
		m[i]=i-1;
		
	}
	n[L-1]=0;
	m[0]=L-1;
	return;
}
double magnetizzazione(){
	double Re,Im;
	double m;
	int i,j;
	int F[Q];
	Re=0;
	Im=0;
	for(i=0;i<Q;i++)	F[i]=0;
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
			F[lattice[i][j].s]++;
			}
		}
	for (i=0;i<Q;i++){
		Re+=F[i]*cos(2*PI*i/Q);
		Im+=F[i]*sin(2*PI*i/Q);
	}
	m=sqrt((pow(Re,2)+pow(Im,2)));
	m=m/(L*L);
	return m;
	}

double energia(int *n, int *m){
	double e;
	int i, j;
	int ip,im,jp,jm;
	e=0;
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
			ip=n[i];
			im=m[i];
			jp=n[j];
			jm=m[j];
			if (lattice[i][j].s==lattice[i][jp].s)	e=e+0.5;
			if (lattice[i][j].s==lattice[i][jm].s)	e=e+0.5;
			if (lattice[i][j].s==lattice[ip][j].s)	e=e+0.5;
			if (lattice[i][j].s==lattice[im][j].s)	e=e+0.5;
		}
	}
	e=e/(L*L);
	return e;
}
void esploro (int *n, int *m,int i,int j, int snew, double h){
	double y;
	lattice[i][j].flag=1;
	
	
	if(lattice[i][n[j]].flag==0 && lattice[i][j].s==lattice[i][n[j]].s){
		y=(double) rand()/RAND_MAX;
		if(y<1-h)	esploro (n,m,i,n[j],snew,h);
	}
	else	lattice[i][n[j]].flag=1;

	if(lattice[i][m[j]].flag==0 && lattice[i][j].s==lattice[i][m[j]].s){
		y=(double) rand()/RAND_MAX;
		if(y<1-h)	esploro (n,m,i,m[j],snew,h);	
	}
	else lattice[i][m[j]].flag=1;

	if(lattice[n[i]][j].flag==0 && lattice[i][j].s==lattice[n[i]][j].s){
		y=(double) rand()/RAND_MAX;
		if(y<1-h)	esploro (n,m,n[i],j,snew,h);
	}
	else	lattice[n[i]][j].flag=1;

	if(lattice[m[i]][j].flag==0 && lattice[i][j].s==lattice[m[i]][j].s){
		y=(double) rand()/RAND_MAX;
		if(y<1-h)	esploro (n,m,m[i],j,snew,h);
	}
	else lattice[m[i]][j].flag=1;

	lattice[i][j].s=snew;
	return;
}
void wolff (int *n, int *m, double h){
	int i,j;
	int snew,sold;
	i=rand()%L;
	j=rand()%L;
	sold=lattice[i][j].s;
	snew=spinew(i,j);
	esploro (n,m,i,j, snew,h);
	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
			lattice[i][j].flag=0;
		}
	}
	return;
}


int main (void){
	double beta;
	double h;
	FILE *fp;
	int i,j,b;
	int npp[L],nmm[L];
	char datafile[100];
	double M,E,E2;
	beta=0.8;
	h=exp(-beta);
	srand(time(NULL));
 	sprintf(datafile, "2wolfprova_%d.txt", Q);
 	fp=fopen(datafile, "w");
		if(fp==NULL) {
		printf("Errore nell'apertura del data'");
		exit(1);
	}
		h=exp(-beta);
		initialize();
 		geometry(npp,nmm);
	 	E=0;
	 	M=0;
	 	E2=0;
	 	
	 	for(b=0;b<stat;b++){
		 	wolff(npp,nmm,h);
			M=magnetizzazione();
		 	E=energia(npp,nmm);
		fprintf(fp,"%.4lf	",M);
		fprintf(fp,"%.4lf\n",E);
	}		
	fclose (fp);
	return 0;
}
