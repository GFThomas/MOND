#include "parametres.h"


void particules(double sigma, double rho1, long N_lim, long Nb_part, double Mtot,struct king_model *para_king){


	long i,j=1,part_ID=0,nb_part_tot=0,k,tmp;
	double R,R1,M_shell;
	double *res;
	double w1,w2,w3;
	double part_mass=Mtot/Nb_part; // Masse d'une particules

	double *nb_part_shell;	
	nb_part_shell=malloc((N_lim+1)*sizeof(double));
	double *rest;
	rest=malloc((N_lim+1)*sizeof(double));
	double M_rest=0.0,pos=0.0; // Nombre de particules dans un intervalle et masse restante dnas une coquille
	double dist,omega,theta;
	struct particules *part;
	part=malloc(Nb_part*sizeof(struct particules));



	// Création de nombre aléatoires
	 long idum=-(long)time(NULL)*100^getpid();
	 srand(time(NULL)*1000^getpid());	

	// Calcul du nombre de particules par cellule
	for(i=1;i<=N_lim;i++){
		R=i*delta_r;
		R1=(i-1.0)*delta_r;

		// Calcul du nombre de particules dans un intervalle
		 M_shell=para_king[i].rho*(4.0/3.0*pi*(pow(R,3.0)-pow(R1,3.0)));
		 nb_part_shell[i]=(long)(M_shell/part_mass);
		 nb_part_tot+=nb_part_shell[i];


		//Calcul de la masse restante pour rajouter une part au particule du reste.
		rest[i]=M_shell-(nb_part_shell[i]*part_mass);
		M_rest= M_rest+rest[i];
		if(i==N_lim && nb_part_tot<Nb_part) M_rest=part_mass;
 
		if(M_rest>=part_mass){ // if rest of mass is greater than mass of a particule, we determine it position
			M_rest = M_rest-part_mass;
			for(k=j;k<=i;k++){
				if(k==i) rest[i]=(rest[i]-M_rest);
				 pos=pos+rest[k]*k*delta_r;
			}
				pos=pos/part_mass;
				k=(long)(pos/delta_r);
				j = i;
				rest[i]=M_rest;
				pos = 0.0;
				nb_part_shell[k] = nb_part_shell[k] + 1;
				nb_part_tot=nb_part_tot + 1;
		}
	}

	// Attribution de la position et de la vitesses a chauque particule
	for(i=1;i<=N_lim;i++){
		R=i*delta_r;
		R1=(i-1.0)*delta_r;
		j=1;
		// Give position and velocity for each particules of the shell 
		while(j<=nb_part_shell[i]){

			// Random position of particules
			w1=frand();
			w2=frand();
			w3=frand();

			//spherical coordinates of the particule		
			dist=R1+(R-R1)*w1;
			theta=acos(2.0*w2-1.0);
			omega=2.0*pi*w3;
			// Spherical coordinates to carthesian coordinates
			part[part_ID].x=dist*sin(theta)*cos(omega);
			part[part_ID].y=dist*sin(theta)*sin(omega);
			part[part_ID].z=dist*cos(theta);

			// Attribution des vitesses
			 res=velocity(para_king[i].psi, rho1, sigma);
			 part[part_ID].vx=res[1];
			 part[part_ID].vy=res[2];				
			 part[part_ID].vz=res[3];
			 part[part_ID].mass=part_mass;			
		
			part_ID=part_ID+1; // ID of a particules in the table "part"
			j=j+1;
		}		
	}

printf("Nombre de particules : %ld\n",part_ID);
	// Cree Le fichier ic_part
	FILE * part_file=NULL;
   	part_file=fopen("ic_part","w");
    	if(part_file!=NULL){
        	for(i=0;i<Nb_part;i++){
    			fprintf(part_file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",part[i].x,part[i].y,part[i].z,part[i].vx,part[i].vy,part[i].vz,part[i].mass);
		}
	}
	fclose(part_file);

	return;
}


/**********************************************
 ** Attribution de la vitesse des particules **
 **********************************************/
double * velocity(double psi, double rho1, double sigma){
	double *vitesse;
	double w1,w2,w3,w4;
	double v,fv,energie;
	
	w1=frand();
	w2=frand();
	w3=frand();
	w4=frand();


	v=rtbis(0.0000001, sqrt(2.0*psi), 0.00001, sigma, psi, w4);
	


	vitesse=quaterion(w1, w2 ,w3, 1.0, 0.0, 0.0);

	vitesse[1]=vitesse[1]*v;
	vitesse[2]=vitesse[2]*v;
	vitesse[3]=vitesse[3]*v;

	return vitesse;
}






















/******************************
***    Partie Quaterion    ****
******************************/
double * multiplication(double *u, double d, double a ,double b, double c){
	double *mul;
	mul=malloc(4*sizeof(double));	

	mul[0]=u[0]*d-(u[1]*a+u[2]*b+u[3]*c);
	mul[1]=(u[0]*a+u[2]*c-u[3]*b+u[1]*d);
	mul[2]=(u[0]*b-u[1]*c+u[3]*a+u[2]*d);
	mul[3]=(u[0]*c+u[1]*b-u[2]*a+u[3]*d);
return mul;
}



double * quaterion(double w1, double w2 ,double w3, double x, double y, double z){
	//x,y,z coordonées x,y,z du point que l'on fait tourner.
	double u[4];
	double * fin;
	fin=malloc(4*sizeof(double));
	

// Calcul le quaternion unitaire u
// http://planning.cs.uiuc.edu/node198.html#eqn:shoemake
// http://mathworld.wolfram.com/Quaternion.html
	u[0]=sqrt(w1)*cos(2*pi*w3);
	u[1]=sqrt(1.0-w1)*sin(2*pi*w2);
	u[2]=sqrt(1.0-w1)*cos(2*pi*w2);
	u[3]=sqrt(w1)*sin(2*pi*w3);
	

	fin=multiplication(u,0.0,x,y,z);
	fin=multiplication(fin,u[0],-u[1],-u[2],-u[3]);
return fin;
}









/*******************************
 ** rtbis et fonction energie **
 *******************************/
double f_energie(double v, double sigma, double psi, double a){
	return (integrale(v,sigma,psi)/integrale(sqrt(2.0*psi),sigma,psi) - a);

	
}

double rtbis(double x1, double x2, double xacc, double sigma, double psi, double a)
{
	int j;
	double dx,f,fmid,xmid,rtb;

	f=f_energie(x1, sigma, psi,a);
	fmid=f_energie(x2, sigma, psi,a);
	
	if (f*fmid >= 0.0){
		printf("Root must be bracketed for bisection in rtbis\n");
		
	}
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=f_energie((xmid=rtb+(dx *= 0.5)), sigma, psi,a);
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0){ return rtb;}
		
	}
	printf("Too many bisections in rtbis\n");
	return 0.0;
}








double integrale(double v_max, double sigma, double psi){
	return( -pow(v_max,3.0)/3.0 + 0.5*exp(psi/pow(sigma,2.0)) *pow(sigma,2.0)* ( -2.0*exp(-pow(v_max,2.0)/(2.0*pow(sigma,2.0)))*v_max + sqrt(2.0*pi ) * sigma *erf(v_max/(sqrt(2.0)*sigma))));
}












/****************************
*** Flottants aléatoires  ***
****************************/
double frand(void){
    return ( rand()/(double)RAND_MAX );
}
