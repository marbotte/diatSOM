#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//***********FONCTION PRINCIPALE***********
void som(
  int *roughLength,
  int *finetuneLength,
  int *nb_lignes_donnees,
  int *nb_colonnes_donnees,
  int *nb_hexagones,
  int *dist,
  int *cellsite,
  int *relat,
  double *sigmaDebRough,
  double *sigmaFinRough,
  double *sigmaFinFinetune,
  double *alphadebrough,
  double *alphadebfinetune,
  double *sigmaVec,
  double *alphaVec,
  int *givenInd,
  double *MatDistVect,
  double *VirtualUnitsVect,
  double *DataMatrixVect,
  double *ab)

{
  //************prototypes******************
  int bmu(double *donnees_site,int *nb_hexagones,int *dist,int *nb_colonnes_donnees,double *ab);
  void updateCodebook(double *donnees_site, double *ab, int meilleur_noeud, double sigma, double alpha, double *MatDistVect,  int *nb_colonnes_donnees, int *nb_hexagones);
  void relatCodebook(double *ab, int  *nb_colonnes_donnees, int*nb_hexagones);
  //****************************************
  double alpha,
  sigma,
  facSigma,
  *donnees_site,
  b,
  a;
  
  int t,
  tt,
  i,
  meilleur_noeud,
  compteurRough=0;
  
  srand(time(NULL));
  
  //mode automatique ou manuel pour les indices sigma et alpha
  if(*givenInd==0)
  {
    Rprintf("SOM indices alpha and sigma are calculated during SOM processing\n");
  }
  else
  {
    Rprintf("SOM indices alpha and sigma were given\n");
  }
  //initialisation des valeurs du codebook ab avec les valeurs de VirtualUnitsVect
  for (i=0;i<*nb_hexagones*(*nb_colonnes_donnees);i++){
    ab[i]=VirtualUnitsVect[i];
  }
  
  //******************************************************
  //*****************ROUGH PHASE**************************
  //******************************************************
  
  //pour calculer alpha rough
  b=(double)((double)(*roughLength)-1.00)/25.000;
  a=b*(*alphadebrough);
  //pour calculer sigma rough
  facSigma=(*sigmaFinRough - *sigmaDebRough)/((double)*roughLength - 1.00);
  for (t=1;t<=*roughLength;t++)
  {
    alpha=a/(b+t-1);
    sigma=(*sigmaDebRough)+(t-1)*facSigma;
    //calcul de l'indice aleatoire du site de debut de boucle (cette etape permet d'eviter que ce soit toujour le meme site qui soit en dernier dans la rough phase et donc ait une importance trop forte)
    if(t%5==0){
      Rprintf("Rough phase ... epoch %d \n",t);
    }
    int indice=(int)(rand() / (double)RAND_MAX * (*nb_lignes_donnees ));//renvoie un chiffre aleatoire de 0 à nb_lignes_donnees -1
    //boucle 1 sur les sites (récursivement de l'indice à la fin du tableau)
    for(tt=indice;tt<(*nb_lignes_donnees);tt++)
    {
      compteurRough++;
      if(*givenInd==0)
      {
        //calcul alpha
        alphaVec[compteurRough-1]=alpha;
        //calcul sigma
        sigmaVec[compteurRough-1]=sigma;
      }
      else
      {
        sigma=sigmaVec[compteurRough-1];
        alpha=alphaVec[compteurRough-1];
      }
      //initialisation du pointeur sur la première donnée du site
      donnees_site=(DataMatrixVect+(tt*(*nb_colonnes_donnees)));
      //recherche du BMU
      meilleur_noeud=bmu(donnees_site,nb_hexagones,dist,nb_colonnes_donnees,ab);
      cellsite[tt]=meilleur_noeud+1;
      //mise a jour du codebook
      updateCodebook(donnees_site, ab, meilleur_noeud, sigma, alpha, MatDistVect, nb_colonnes_donnees, nb_hexagones);
      //remise en abondances relatives
      if(*relat==1){
        relatCodebook(ab,nb_colonnes_donnees,nb_hexagones);
      }
    }
    
    //boucle 2 sur les sites (recursivement du début à l'indice)
    for(tt=0;tt<indice;tt++){
      compteurRough++;
      if(*givenInd==0)
      {
        //calcul alpha
        alphaVec[compteurRough-1]=alpha;
        //calcul sigma
        sigmaVec[compteurRough-1]=sigma;
      }
      else
      {
        sigma=sigmaVec[compteurRough-1];
        alpha=alphaVec[compteurRough-1];
      }
      //initialisation du pointeur sur la première donnée du site
      donnees_site=(DataMatrixVect+(tt*(*nb_colonnes_donnees)));
      //recherche du BMU
      meilleur_noeud=bmu(donnees_site,nb_hexagones,dist,nb_colonnes_donnees,ab);
      cellsite[tt]=meilleur_noeud+1;
      //mise a jour du codebook
      updateCodebook(donnees_site, ab, meilleur_noeud, sigma, alpha, MatDistVect, nb_colonnes_donnees, nb_hexagones);
      //remise en abondances relatives
      if(*relat==1){
        relatCodebook(ab,nb_colonnes_donnees,nb_hexagones);
      }
    }
  }
  //******************************************
  //***********FINETUNING PHASE***************
  //******************************************
  //pour calculer alpha finetune
  b=(*finetuneLength-1)/20.00;
  a=b*(*alphadebfinetune);
  //pour calculer sigma finetune
  facSigma=(*sigmaFinFinetune - *sigmaFinRough)/((double)*finetuneLength-1.00);
  for(t=1;t<=*finetuneLength;t++)
    {
      if(t%(10*(*nb_lignes_donnees))==0)
	{ 
	  Rprintf("finetune phase %d / %d \n",t,*finetuneLength);
	}
      tt=(int)(rand() / (double)RAND_MAX * (*nb_lignes_donnees));
      if(*givenInd==0)
      {
	//calcul alpha
	alpha=a/(b+t-1);
	alphaVec[compteurRough+t-1]=alpha;
	//calcul sigma
	sigma=(*sigmaFinRough)+(t-1)*facSigma;
	sigmaVec[compteurRough+t-1]=sigma;
      }
      else
      {
	sigma=sigmaVec[compteurRough+t-1];
	alpha=alphaVec[compteurRough+t-1];
      }
      //initialisation du pointeur sur la première donnée du site
      donnees_site=(DataMatrixVect+tt*(*nb_colonnes_donnees));
      //recherche du BMU
      meilleur_noeud=bmu(donnees_site,nb_hexagones,dist,nb_colonnes_donnees,ab);
      cellsite[tt]=meilleur_noeud+1;
      //mise a jour du codebook
      updateCodebook(donnees_site, ab, meilleur_noeud, sigma, alpha, MatDistVect, nb_colonnes_donnees, nb_hexagones);
      //remise en abondances relatives
      if(*relat==1)
      {
	relatCodebook(ab,nb_colonnes_donnees,nb_hexagones);
      }
      
    }
  
}
int bmu(double *donnees_site, int *nb_hexagones,int *dist,int *nb_colonnes_donnees,double *ab)
{
  double distance=0.0,//distance avec le site calculée sur chaque hexagone
  Numerateur=0.0,//utilisé dans le calcul de la distance de BC
  mindist=10000000.00,//on range la distance minimale dans mindist pour la comparer avec distance
  SomCell=0.0,//somme totale de la cellule pour le calcul de BC
  SomSite=0.0;//somme totale du site pour le calcul de BC
  int jj,//boucle sur colonnes
  ii,//boucle sur hexagones
  meilleur_noeud;//resultat de la fonction à retourner
  
  
  for (ii=0;ii<*nb_hexagones;ii++)
  {
    distance=0.0;
    Numerateur=0.0;
    SomCell=0.0;
    SomSite=0.0;
    
    switch (*dist)
    {
      //euclidienne
      case 1:
	for (jj=0;jj<*nb_colonnes_donnees;jj++)
	{
	  distance += pow((ab[ii*(*nb_colonnes_donnees)+jj] - donnees_site[jj]),2);
	}
	distance =sqrt(distance);
	
	break;
	//Bray Curtis
      case 2:
	for (jj=0;jj<*nb_colonnes_donnees;jj++)
	{
	  SomCell+=ab[ii*(*nb_colonnes_donnees)+jj];
	  SomSite+=donnees_site[jj];
	  Numerateur+=fabs(ab[ii*(*nb_colonnes_donnees)+jj]-donnees_site[jj]);
	}
	distance=Numerateur/(SomCell+SomSite);
	break;
	//Manhattan
      case 3:
	for (jj=0;jj<*nb_colonnes_donnees;jj++)
	{
	  distance+=fabs(ab[ii*(*nb_colonnes_donnees)+jj]-donnees_site[jj]);
	}
	break;
      default:
	Rprintf("Distance non prise en compte\n");
    }
    //test de la distance et affectation du BMU
    if (distance<mindist){
      mindist=distance;
      meilleur_noeud=ii;
    }
  }
  
  return(meilleur_noeud);
}
void updateCodebook(double *donnees_site, double *ab, int meilleur_noeud, double sigma, double alpha, double *MatDistVect,  int *nb_colonnes_donnees, int *nb_hexagones)
{
  double coef;//partie de la mise a jour qui ne changera pas au sein d'un hexagone
  int ii,//boucle sur les cellules
  jj;//boucle sur les colonnes
  
  for (ii=0; ii<*nb_hexagones;ii++){
    coef=exp(-pow((MatDistVect[ii*(*nb_hexagones) + meilleur_noeud]),2)/(0.5*pow(sigma,2)));//a verifier: que l'on prenne la bonne distance
    //Rprintf(" %f %d ",coef, sigma);
    for (jj=0;jj<*nb_colonnes_donnees;jj++){
      ab[ii*(*nb_colonnes_donnees)+jj] += alpha*coef*(donnees_site[jj]-(ab[ii*(*nb_colonnes_donnees)+jj]));
    }
  }
}
void relatCodebook(double *ab, int  *nb_colonnes_donnees, int *nb_hexagones)
{
  int ii,//boucle
  jj;
  double SomCell;
  for(ii=0;ii<*nb_hexagones;ii++)
  {
    SomCell=0.0;
    for(jj=0;jj<*nb_colonnes_donnees;jj++){
      SomCell+=ab[ii*(*nb_colonnes_donnees)+jj];
    }
    for(jj=0;jj<*nb_colonnes_donnees;jj++){
      ab[ii*(*nb_colonnes_donnees)+jj]/=SomCell;
    }
  }
}

