#include "parametres.h"

struct king_model * king(double sigma, double r_king, double W0, char name_nu[], double fact_G, double *rho1, long *N_lim, double *Mtot){


	long i,tmp;
	double psi_NM[N_max],*pot;				// Potentiel calcule avec poisson en Newton
	double R,R1,R2;
	double rh,mh;
   double C,rt;
	struct king_model *para;
   	*Mtot=0.0;
	para=malloc(N_max*sizeof(struct king_model));		// Alloue la mémoire

	// Calcul du potentiel et densite centraux
	 para[1].psi=W0*pow(sigma,2.0);								// Calcul du potentiel relatif central à partir du ratio W0 et de sigma
	 para[1].rho=(9.0*sigma*sigma)/(4.0*pi*G*fact_G*pow(r_king,2.0));	// Calcul de la densite centrale à partir de r_king et sigma ( def du rayon de King)
	// Grad psi=0 au centre
	 para[2].psi=para[1].psi;
	 para[2].rho=para[1].rho;
	// Potentiel Mondien = potentiel Newtonien au centre
	 psi_NM[1]=para[1].psi;
	 psi_NM[2]=para[2].psi;

	// Calcul de la normalisation rho1
		*rho1=para[1].rho/(exp(para[1].psi/pow(sigma,2.0))*erf(sqrt(para[1].psi)/sigma)-sqrt((4.0*para[1].psi)/(pi*pow(sigma,2.0)))*(1.0+(2.0*para[1].psi)/(3.0*pow(sigma,2.0))));

	// Calcul du potentiel de la densite et de la masse pour tout les intervalle jusque rt
	for(i=1;i<N_max;i++){
		// Definition des rayons des intervalles :
		 R=i*delta_r;
		 R1=(i-1.0)*delta_r;
		 R2=(i-2.0)*delta_r;
		if(i>2){ // Pour les intervalles autres que les deux premiers centraux
			
			para[i].psi=potentiel(R2, R1, psi_NM[i-2], psi_NM[i-1], para[i-1].psi, sigma , para[i-2].rho, name_nu,fact_G, &psi_NM[i]); // Calcul du potentiel

			para[i].rho=densite(sigma, *rho1, para[i].psi);
			
		}
        if(para[i].rho <=0.0 || para[i].psi <= 0.0){
            *N_lim=i-1;
	    tmp=i-1;
	
            rt=R1;
            C=log10(rt/r_king);
            break;
        }

        // Calcul de la masse inscrite et de la masse totale
         	para[i].mass=*Mtot+para[i].rho*(4.0/3.0*pi*(pow(R,3.0)-pow(R1,3.0)));
		 *Mtot=*Mtot+para[i].rho*(4.0/3.0*pi*(pow(R,3.0)-pow(R1,3.0)));
	}
/*
	// Calcul le vrai potentiel
	pot=malloc(tmp*sizeof(double));
	pot[tmp]=-G*(*Mtot)/rt;
	for(i=1;i<tmp;i++){
		pot[i]=pot[tmp]-para[i].psi;
	}
*/
for(i=1;i<=(*N_lim);i++){
	if(para[i].mass>=(*Mtot/2.0)) break;
	rh=i*delta_r;
	mh=para[i].mass;
	
}

	// Creer le fichier king.out contenant les paramètres de King en fonction du rayon
	FILE * king_file=NULL;
	king_file=fopen("king.out","w");
	if(king_file!=NULL){
		fprintf(king_file,"#R\trho\tpsi\tacceleration\tmass\n");
		for(i=1;i<=(*N_lim);i++){
			R=i*delta_r;
			fprintf(king_file,"%lf\t%g\t%lf\t%lf\t%g\n",R,para[i].rho,para[i].psi,-(para[i].psi-para[i-1].psi)/delta_r,para[i].mass);
		}
	}
	fclose(king_file);

    //Entre les données dans table_of_parameters.out
    FILE * table_file=NULL;
    table_file=fopen("table_of_parameters.out","a+");
    if(table_file!=NULL){
        fprintf(table_file,"sigma= %lf ro = %lf W0= %lf --> rt= %lf C= %lf mass= %.4g rh= %lf Mh= %.4g\n Gnorm=G*%lf",sigma,r_king,W0,rt,C,*Mtot,rh,mh,fact_G);		
	}
	fclose(table_file);
  


	return para;
}


















/*************************************
 ** Fonction de calcul du potentiel **
 *************************************/
double potentiel(double R2, double R1, double psi_NM2, double psi_NM1, double psi1, double sigma , double rho1, char name_nu[], double fact_G, double *psi_NM){
	double psi,nu,y;
	double gamma=4.0;
	
	// Calcul du potentiel par Poisson classique
	 *psi_NM=delta_r*pow(R2/R1,2.0)*((psi_NM1-psi_NM2)/delta_r-4.0*pi*G*fact_G*delta_r*rho1)+psi_NM1;

	// Choix de la fonction nu
	 y=fabs((*psi_NM-psi_NM1)/delta_r)/a0;
	 if(strcmp(name_nu,"simple")==0){
		nu=0.5*sqrt(1.0+4.0/y)+0.5;
	}else if (strcmp(name_nu,"std")==0 || strcmp(name_nu,"standard")==0){
		nu=pow((1.0+pow(1.0+4.0*pow(y,-2.0),0.5))/2.0,0.5);
	}else if (strcmp(name_nu,"exp")==0 || strcmp(name_nu,"exponentiel")==0){
		nu=pow(1.0-exp(-pow(y,(gamma/2.0))),-1.0/gamma)+(1.0-1.0/gamma)*exp(-pow(y,(gamma/2.0)));

	}else if (strcmp(name_nu,"newton")==0 || strcmp(name_nu,"pot_ext")==0){
		nu=1.0;
	}

	// Calcul du potentiel relatif MONDien
	 psi=nu*(*psi_NM-psi_NM1)+psi1;
	return psi;
}


/**************************************
 ** Fonction de calcul de la densite **
 **************************************/
double densite(double sigma, double rho1, double psi){
	double rho;

        rho=rho1*(exp(psi/pow(sigma,2.0))*erf(sqrt(psi)/sigma)-sqrt((4.0*psi)/(pi*pow(sigma,2.0)))*(1.0+(2.0*psi)/(3.0*pow(sigma,2.0))));
	return rho;
}
