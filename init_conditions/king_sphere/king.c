#include "parametres.h"

int main(int argc, char *argv[]){
	long Nb_part;
	long N_lim;				// N° de l'intervalle final du modele de King (défini rt)
	double sigma,r_king,W0;
	double rho1,Mtot,fact_G,g_ext,mue,Le,x; 			// Facteur de normalisation et masse totale et Facteur de renormalisation de G si pot_ext majoritaire
	char name_nu[256];
	struct king_model *para_king;

	// ENTREE
	scanf("%lf",&sigma);	// Dispersion de vitesse
	scanf("%lf",&r_king);	// Rayon de king
	scanf("%lf",&W0);		// Rapport psi0/sigma2
	scanf("%ld",&Nb_part);	// Nombre de particules souhaitees
	scanf("%s",name_nu);	// Choix de quel nu function on souhaite
	if (strcmp(name_nu,"pot_ext")==0){
		scanf("%lf",&g_ext); // g ext in km2.s-2.kpc-1
		x=g_ext/a0;
		mue=x/(1.0+x);	
		Le=(1.0-x/(1.0+x));
		fact_G= 27.5;//1.0/(mue*(1.0+Le));
		printf("Fact G = %lf\n",fact_G);
	}else	fact_G=1.0;

	para_king=king(sigma,r_king,W0,name_nu,fact_G,&rho1,&N_lim,&Mtot);	// Calcul la densite, le potentiel relatif et la mass suivant un profil de King (1966)
	particules(sigma,rho1,N_lim,Nb_part,Mtot,para_king);

	return 0;
}




