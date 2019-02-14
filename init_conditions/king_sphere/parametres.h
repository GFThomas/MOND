#define _FILE_OFFSET_BITS 64 
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

/**************************
 ** Constantes physiques **
 **************************/
#define pi 3.141592653589793
#define G 4.302113489*pow(10.0,-6.0)		// kpc km2 s-2 M⊙-1
#define a0 3456.790123		// km2 s-2 kpc-1

/*****************************
 ** Paramrtres du programme **
 *****************************/
#define N_max 1000000		// Nombre d'intervalle maximal
#define delta_r	0.0001		// intervalle de rayon (kpc)
#define JMAX 40

/****************
 ** Structures **
 ****************/
struct king_model{
	double psi,rho,mass;
};

struct particules{
	double x,y,z;
	double vx,vy,vz;
	double mass;
};

/****************
 ** Prototypes **
 ****************/
struct king_model * king(double sigma, double r_king, double W0, char name_nu[], double fact_G, double *rho1, long *N_lim, double *Mtot);
double potentiel(double R2, double R1, double psi_NM2, double psi_NM1, double psi1, double sigma , double rho1, char name_nu[], double fact_G, double *psi_NM);
double densite(double sigma, double rho1, double psi);
void particules(double sigma, double rho1, long N_lim, long Nb_part, double Mtot,struct king_model *para_king);
double * quaterion(double w1, double w2 ,double w3, double x, double y, double z);
double * multiplication(double *u, double d, double a ,double b, double c);
double * velocity(double psi, double rho1, double sigma);
double frand(void);
double f_energie(double v, double sigma, double psi, double a);
double rtbis(double x1, double x2, double xacc, double sigma, double psi, double a);
double integrale(double v_max, double sigma, double phi);










