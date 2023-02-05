#include<math.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>


#define D 2 
#define L 50
#define deco 30
#define extf 0
#define PI 3.14159265

const int stat=50000;
const int vicini[5]={4,8,12,16,20};
int N[2];

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
	N[0]=L*L;
	N[1]=0;
	return;
}

int spinew (int i, int j){
	int z;
	int sold;
	int snew;
	sold=lattice[i][j].s;
	z=(int) rand()%2;
	while (z==sold){
		z= (int) rand()%2;
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


double magnetizzazione(double m){
	double n;
	int i,j;
	n=abs(N[0]-N[1]);
	m=(double) n/(L*L);
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
			if (lattice[i][j].s==lattice[i][jp].s)	e=e-0.5;
			else e=e+0.5;
			if (lattice[i][j].s==lattice[i][jm].s)	e=e-0.5;
			else e=e+0.5;
			if (lattice[i][j].s==lattice[ip][j].s)	e=e-0.5;
			else e=e+0.5;
			if (lattice[i][j].s==lattice[im][j].s)	e=e-0.5;
			else e=e+0.5;
		}
	}
	e=e/(L*L);
	return e;
}
void esploro (int *n, int *m,int i,int j, int snew, double h){
	double y;
	lattice[i][j].flag=1;
	N[lattice[i][j].s]-=1;
	
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
	N[lattice[i][j].s]+=1;
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



void mediadev(double *vet, double *d){
	int i;
	double m,m2;
	m=0;
	m2=0;
	
	for (i=0;i<stat;i++){
		m+=vet[i];
		m2+=vet[i]*vet[i];
	}
	m=m/stat;
	m2=m2/stat;
	d[0]=m;
	d[1]=sqrt(m2-(m*m)/stat);
	return;
}


void bootstrap (double *vet, double *m, int k) {
	int i,j,z;
	double medie[k],medie2[k];
	m[0]=0;
	m[1]=0;
		
	for (i=0;i<k;i++){
		medie[k]=0;
		medie2[k]=0;
		for(j=0;j<stat;j++){
			z=rand()%stat;
			medie[k]+=vet[z];
			medie2[k]+=vet[z]*vet[z];
		}
		medie[k]=medie[k]/stat;
		medie2[k]=medie2[k]/stat;
		m[0]+=medie[k];
		m[1]+=medie2[k];
	}
	m[0]=(double) m[0]/k;
	m[1]=(double) sqrt(m[1]-m[0]*m[0]);
	return;
}

void bootstrap2 (double *vet1, double *vet2, double *x, int k, double t){
	int i,j,z;
	double m1[k],m2[k],d1[k],d2[k];
	double y,y2;
	y=0;
	y2=0;
	
	for(i=0;i<k;i++){
		m1[i]=0;
		m2[i]=0;
		for(j=0;j<stat;j++){
			z=rand()%stat;
			m1[i]+=vet1[z];
			m2[i]+=vet2[z];
		}
		m1[i]=m1[i]/stat;
		m2[i]=m2[i]/stat;
		y+=(m2[i]-m1[i]*m1[i])/t;
		y2+=pow(((m2[i]-m1[i]*m1[i])/t),2);
		
	}
	x[0]=(y/k)*L*L;
	x[1]=L*L*sqrt(y2/k-pow((y/k),2));
	

	return;
}

void bootstrappo (double *vet, double *m, int k) {
	int i,j,z;
	double medie[k],medie2[k];
	m[0]=0;
	m[1]=0;
		
	for (i=0;i<k;i++){
		medie[k]=0;
		medie2[k]=0;
		for(j=0;j<stat;j++){
			z=rand()%stat;
			medie[k]+=vet[z];
			medie2[k]+=vet[z]*vet[z];
		}
		medie[k]=medie[k]/stat;
		medie2[k]=medie2[k]/stat;
		m[0]+=medie[k];
		m[1]+=medie2[k];
	}
	m[0]=(double) m[0]/k;
	m[1]=(double) m[1]/k;
	m[1]=(double) m[1]-m[0]*m[0];
	return;
}
double mean(double* array, int length) {
    double sum = 0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
    return sum / length;
}


int main (void){
	double t,beta;
	double h;
	FILE *fp;
	int i,j,a,b,g;
	int npp[L],nmm[L];
	char datafile[100];
	double mana[stat], ene[stat], ene2[stat],mana2[stat];
	double dm[2],de[2],dene[2],dx[2],dh[2],dmana[2];
	double M;
	int z,tran;
	beta=0.35;
	t=1/beta;
	tran=40;
	h=exp(-2*beta);
	srand(time(NULL));
 	sprintf(datafile, "wolf_%d.txt", L);
 	fp=fopen(datafile, "w");
		if(fp==NULL) {
		printf("Errore nell'apertura del data'");
		exit(1);
	}
	fprintf(fp, "#t	m	dm	e	de	x	dx	h	dh\n");
	for (g=0;g<50;g++){
		printf("%d\n", g);
		beta=beta+0.003;
		t=1/beta;
		h=exp(-2*beta);
		initialize();
 		geometry(npp,nmm);
 		for (z=0;z<tran;z++){
 			wolff(npp,nmm,h);
		 }
	 	for(a=0;a<stat;a++){
	 		for(b=0;b<deco;b++){
	 			wolff(npp,nmm,h);
			}
			
			M=magnetizzazione(M);
	 		ene[a]=energia(npp,nmm);
	 		ene2[a]=ene[a]*ene[a];
			mana[a]=M;
			mana2[a]=M*M;
			
		}
		fprintf(fp,"%.4lf	",beta);
		bootstrappo(mana,dm,40);
		bootstrappo(ene,de,40);
		bootstrappo(ene2,dh,40);
		bootstrappo(mana2,dx,40);
		fprintf(fp,"%.4lf	",dm[0]);
		fprintf(fp,"%.4lf	",dm[1]);
		fprintf(fp,"%.4lf	",de[0]);
		fprintf(fp,"%.4lf	",de[1]);
		fprintf(fp,"%f	",dx[0]);
		fprintf(fp,"%f	",dx[1]);
		fprintf(fp,"%f	",dh[0]);
		fprintf(fp,"%f\n",dh[1]);
		
	
}
	fclose (fp);
	return 0;
}

