####Attention le nombre minimale d'étape a été modifié en 10 rough 20 finetune dans ce script'


################################################################
################################################################
#############Fonctions accessoires##############################
################################################################
################################################################
################################################################

################################################################
#####Pour l'initialisation linéaire du codebook#################
################################################################

normVec=function(vec)
##permet de calculer la norme d'un vecteur
{
  return(sqrt(vec%*%vec))
}

normValMainEig=function(tab)
##Fonction qui sert dans l'initialisation linéaire du codebook: traduction exacte de la somtoolbox##PAS SUR: sert à calculer la norme des vecteurs propres
{
  AutoCov=cov(tab)
  A=eigen(AutoCov)
  res=apply(A$vectors[,1:2],2,function(vec)vec/normVec(vec))
  for(i in 1:2){res[,i]=res[,i]*sqrt(A$values[i])}
  return(res)
}

linInit=function(matrice,nbHex,coord)
##initiation linéaire des valeurs du codebook : traduction exacte de la som toolbox
{
  me=colMeans(matrice)
  D=t(t(matrice)-me)
  vecMain=normValMainEig(D)
  codebook= matrix(rep(me,each=nbHex),ncol=ncol(matrice))
  CoordNorm=(apply(coord,2,function(vec)(vec-min(vec))/(max(vec)-min(vec)))-.5)*2
  for(i in 1:nbHex)
  {
    vecMainFac=t(CoordNorm[i,1:2]*t(vecMain))
    codebook[i,]=codebook[i,]+vecMain[,1]-vecMainFac[,2]
  }
  return(codebook)
}


###############################################################
#########Calcul des erreurs et BMU#############################
###############################################################
searchBMU_1_2=function(data,codebook,distUsed)
{

  BMU1_2=t(switch(distUsed,
                  
                 "euclidean"=apply(data,1,function(Vec,Tab) 
                 {
                   sort(apply(Tab,1,function(Vec2,Vec1)
                             { sqrt(sum((Vec2-Vec1)^2)) }
                             ,Vec1=Vec),index.return=T)$ix[1:2] 
                 },Tab=codebook),
                 
                 "BrayCurtis"=apply(data,1,function(Vec,Tab)
                 {
                   sort(apply(Tab,1,function(Vec2,Vec1)
                             { sum(abs(Vec1-Vec2))/(sum(Vec1)+sum(Vec2)) }
                             ,Vec1=Vec),index.return=T)$ix[1:2]
                 },Tab=codebook),
                 
                 "manhattan"=apply(data,1,function(Vec,Tab)
                 {
                   sort(apply(Tab,1,function(Vec2,Vec1)
                               {sum(abs(Vec1-Vec2)) } 
                               ,Vec1=Vec),index.return=T)$ix[1:2]
                 },Tab=codebook)
         ))
  return(BMU1_2)
}

################################################################
################################################################
###Chargement des librairies dynamiques#########################
#####et gestion des OS##########################################
################################################################
################################################################
################################################################
if(Sys.info()['sysname']=='Windows')
{
sample.int=function(n,m,...){return(sample(1:n,m,...))}
}

################################################################
##Fonction proprement dite
###############################################################
SOM=function(tableau,#data
        nbRowSom,#output's number of rows
        nbColSom,#output's number of cols
        dista,#distance measure in SOM processing
        roughLength=NA,#number of steps of the rough phase
        finetuneLength=NA,#number of steps of the finetuning phase
        randomize=T,# Should data be randomized before SOM processing
        codebookInit=NULL,#manually initialised codebook, if null lininit is assumed
        lininit=is.null(codebookInit),#linear intialisation of codebook
        randinit=!lininit&is.null(codebookInit),#random initialisation of codebook (random sampling into data)
        relat = F,#should data be transformed into relative data
        keepData = F ,#should all be saved during SOMprocess?
        keepCdbIni=F,#should initialized codebook be returned
        alphadebrough=0.8,#initial value of alpha in rough phase
        alphadebfinetune=0.2,#initial value of alpha in finetune phase (and final value of alpha in rough phase)
        sigmaDebRough=max(3,sqrt((0.5*sqrt(3)*nbRowSom)^2+(nbColSom)^2)),#initial value of sigma in rough phase
											       # (default to maximal distance on the map)
        sigmaFinRough=max(1.5,sigmaDebRough/2),#final value of sigma in rough phase
        sigmaFinFinetune=1.1,#final value of sigma in finetuning phase
        sigmaVec=NULL,
        alphaVec=NULL,
        keepIndices=F,#should all alpha and sigma indices be returned
        onlyPos=F)
{
########INITIALISATION DES DIFFERENTS PARAMETRES ET VALEURS DE LA SOM#############
nb_lignes_donnees=nrow(tableau)
nb_colonnes_donnees=ncol(tableau)
nb_hexagones=nbRowSom*nbColSom

if(length(rownames(tableau))==0)
{
rownames(tableau)=as.character(1:nrow(tableau))
warnings("data matrix without rownames, using 1:nrow(tableau) instead")
}
rowOrd=rownames(tableau)
if(randomize)
  {
    tableau=tableau[sample.int(nrow(tableau),nrow(tableau)),] 
  }
  distancesCodees=c('euclidean','BrayCurtis','manhattan')

  cellsite=rep(0,nb_lignes_donnees)
  
  if(relat)
  {
    tableau=t(apply(tableau,1,function(x)x/sum(x)))
  }
  
  tableau=as.matrix(tableau)
  DataMatrixVect=as.vector(t(tableau))
  distUse=match.arg(dista,choices=distancesCodees)
  
  ######calcul des coordonnées des hexagones sur la carte et des distances qui les séparent
  xmat=matrix(1:nbColSom,nrow=nbRowSom,ncol=nbColSom,byrow=T);
  xmat[which(1:nbRowSom%%2==0),]=xmat[which(1:nbRowSom%%2==0),]+.5#xmat et ymat calculés ainsi permettent d'avoir
  #une distance de 1 entre les hexagones qui se touchent
  ymat=matrix(seq(from=0,to=(0.5*sqrt(3))*(nbRowSom-1),by=.5*sqrt(3)),nrow=nbRowSom,ncol=nbColSom)
  CentrHex=cbind(as.double(t(xmat)),as.vector(t(ymat)))#coordonnées complètes sous forme de tableau
  MatDistVect=as.vector(t(as.matrix(dist(CentrHex))))#calcul des distances entre hexagones et mise en forme pour C
  
  ######initialisation linéaire du codebook grâce aux vecteurs propres normalisés

  if(lininit&is.null(codebookInit))
  {
    cat('initalisation linéaire du codebook\t...\t')
    codebook=linInit(tableau,nb_hexagones,CentrHex)
    cat('done\n')
  }
  ######initialisation aléatoire du codebook####
  if(randinit&!lininit&is.null(codebookInit))
  {
    codebook=as.matrix(tableau[sample(1:nb_lignes_donnees,nb_hexagones),])
  }

  ######CodebookInit#############################
  if(!is.null(codebookInit))
    {
      codebook=codebookInit
    }

  ####codebook
  if(relat)
  {
    codebook=t(apply(codebook,1,function(x)x/sum(x)))
  }
  
  if(onlyPos)
  {
    codebook[codebook<0]=0
  }
  
  rownames(codebook)=1:nrow(codebook)
  colnames(codebook)=colnames(tableau)
  VirtualUnitsVect=as.vector(t(codebook))#codebook initial (ne sera pas changé dans la fonction en C
  ab=double(length(VirtualUnitsVect))#codebook final (modifié en C)
  
  ##calcul du nombre d'étape dans les différentes phases##
  if(is.na(roughLength)|is.na(finetuneLength))
  {
    mpd=(nb_hexagones*log(nb_colonnes_donnees))/nb_lignes_donnees
    if(is.na(roughLength)){roughLength=max(round(10*mpd),10)}
    if(is.na(finetuneLength)){finetuneLength=max(round(40*mpd)*nb_lignes_donnees,20*nb_lignes_donnees) }
  }
  cat('rough phase : ',roughLength,' epoch \n')
  cat('finetune phase : ',finetuneLength,' real steps\n')
  
  ##Cas ou l'on règle "à la main" les indices sigma et alpha
  givenIndices=(!is.null(alphaVec)|!is.null(sigmaVec))
  
  if(!givenIndices){sigmaVec<-alphaVec<-double(roughLength*nb_lignes_donnees+finetuneLength)}
  
  if(givenIndices&&(is.null(alphaVec)|is.null(sigmaVec)))
  {
    stop('You must fill in ',ifelse(is.null(alphavec),'alpha','sigma'),
         ' if you wish to run SOM algorithm in manual mode')
  }
  #############APPEL DE LA FONCTION C###################

  result=.C("som",
            as.integer(roughLength),
            as.integer(finetuneLength),
            as.integer(nb_lignes_donnees),
            as.integer(nb_colonnes_donnees),
            as.integer(nb_hexagones),
            as.integer(factor(distUse,levels=distancesCodees)),
            as.integer(cellsite),
            as.integer(relat),
            as.double(sigmaDebRough),
            as.double(sigmaFinRough),
            as.double(sigmaFinFinetune),
            as.double(alphadebrough),
            as.double(alphadebfinetune),
            as.double(sigmaVec),
            #vecteur qui contiendra les facteurs sigma à la fin, ou qui les contient dès le début 
            as.double(alphaVec),
            #vecteur qui contiendra les facteurs alpha,ou qui les règle pour l'algo
            as.integer(givenIndices),
            as.double(MatDistVect),
            as.double(VirtualUnitsVect),
            as.double(DataMatrixVect),
            as.double(ab),
					 PACKAGE="diatSOM")
  
  ######RECUPERATION DES RESULTATS DE LA FONCTION C
  names(result)=c('roughLength',
                  'finetuneLength',
                  'nbRowData',
                  'nbColData',
                  'nbCell',
                  'distUse',
                  'BMUproc',
                  'relat',
                  'sigmaDebRough',
                  'sigmaFinRough',
                  'sigmaFinFinetune',
                  'alphaDebRough',
                  'alphaDebFinetune',
                  'sigmaVec',
                  'alphaVec',
                  'givenIndices',
                  'DistCell',
                  'codebookIni',
                  'Data',
                  'codebook')
  result$sigma=c(result$sigmaDebRough,result$sigmaFinRough,result$sigmaFinRough,
                 result$sigmaFinFinetune)
  result$distUse=distancesCodees[result$distUse]
  result$Data=matrix(result$Data,nrow=nrow(tableau),ncol=ncol(tableau),byrow=T)
  attributes(result$Data)=attributes(tableau)
	result$Data=result$Data[rowOrd,]
  result$alpha=c(result$alphaDebRough,result$alphaDebFinetune)
  result$codebookIni=matrix(result$codebookIni,nrow=nrow(codebook),ncol=ncol(codebook),byrow=T)
  attributes(result$codebookIni)=attributes(codebook)
  result$codebook=matrix(result$codebook,nrow=nrow(codebook),ncol=ncol(codebook),byrow=T)
  attributes(result$codebook)=attributes(codebook)
  #####RESULTATS ADDITIONNELS
  result$CoordHex=CentrHex
  colnames(result$CoordHex)=c('x','y')
  names(result$BMUproc)=rownames(tableau)
result$BMUproc=result$BMUproc[rowOrd]
  result$DistCell=matrix(result$DistCell,nrow=result$nbCell,ncol=result$nbCell,byrow=T)
  ##recherche des BMU et calcul des erreurs
  cat("Calcul des erreurs et recherche des BMUdef\t...")
    #BMU 1 & 2
  result$BMU1_2=searchBMU_1_2(data=result$Data,codebook=result$codebook,distUsed=result$distUse)
    #erreur quant
  result$error=list()
  result$error$quant=sum(sapply(rownames(result$BMU1_2),function(Site,SOM)
    {switch(SOM$distUse,
      "euclidean"=sqrt(sum((SOM$Data[Site,]-SOM$codebook[SOM$BMU1_2[Site,1],])^2)),
      "BrayCurtis"=sum(abs(SOM$Data[Site,]-SOM$codebook[SOM$BMU1_2[Site,1],]))/(sum(SOM$Data[Site,])+
  sum(SOM$codebook[SOM$BMU1_2[Site,1],])),
      "manhattan"=sum(abs(SOM$Data[Site,]-SOM$codebook[SOM$BMU1_2[Site,1],]))
      )
    },SOM=result)
  )              
    #erreur topo
  result$error$topo=length(which(apply(result$BMU1_2,1,function(bmus,distCell)distCell[bmus[1],bmus[2]],
  distCell=result$DistCell)>1.001))/result$nbRowData
  cat("done\n") 
  result=result[c('roughLength',"finetuneLength",'distUse','sigma','alpha',if(keepIndices)'sigmaVec',if(keepIndices)'alphaVec','CoordHex',if(keepData)'Data',
  if(keepCdbIni)'codebookIni','BMUproc' , 'BMU1_2' ,'codebook','error')]

  result$taille=c(nbRowSom,nbColSom)
  result$call=match.call()
  ##########LIBERATION DE LA MEMOIRE ET RESULTATS FINAUX
  rm(list=c('ab','cellsite','CentrHex','codebook','DataMatrixVect','MatDistVect',
  'tableau','VirtualUnitsVect','xmat','ymat'));gc()
  return(result)
}
