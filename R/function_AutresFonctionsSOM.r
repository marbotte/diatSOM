hexagone=function(xcenter,ycenter,d,col=NULL,...)
  {xpoints=xcenter+c(0,d,d,0,-d,-d)
  ypoints=ycenter+c(d*2/sqrt(3),.5*d,-.5*d,-d*2/sqrt(3),-.5*d,.5*d)
  polygon(xpoints,ypoints,col=col,...)
  }

SOMshowColscale=function(col,coordHex,tit=NA,d=F,placeLegendRight=0)#reste à coder le cas où les hexagones font des
	#tailles différentes
  {
  if(ncol(coordHex)!=2)
    {stop('Les coordonnées des cellules doivent être données dans un tableau à deux colonnes\n("x","y")\n')
    }
  plot(0,xlim=range(coordHex[,'x'])+c(-1,1)+c(0,placeLegendRight),ylim=range(coordHex[,'y'])+c(-1,1),ylab="",xlab="",
	type="n",axes=F,main=tit)
  if(all(!d))
    {apply(cbind(coordHex,col),1,function(vecCoordCol)
      {hexagone(xcenter=as.numeric(vecCoordCol[1]),ycenter=as.numeric(vecCoordCol[2]),d=.5,col=vecCoordCol[3])
      })
    }else{
    tabCoord=cbind(coordHex,col,d)
    apply(tabCoord,1,function(vecCoordCol)
    {hexagone(xcenter=as.numeric(vecCoordCol[1]),ycenter=as.numeric(vecCoordCol[2]),d=as.numeric(vecCoordCol[4]),col=
    vecCoordCol[3])
    })
  }
  }
minMaxScaleCut=function(scaleCut,index)#fonction qui permet de donner la première ou la deuxième valeur d'un interval
	# retourné par cut
  {if(!any(c(1,2)==index))
    {stop('l\'index doit être 1 (min) ou 2 (max)')}
  scaleCut=as.character(scaleCut)
  val=sapply(strsplit(scaleCut,','),function(x)x[index])
  res=switch(index,{substr(val,2,nchar(val))},{substr(val,1,nchar(val)-1)})
  return(res)
  }

show.specie=function(specie,SOMresult,nbscale=10,colscale=NA,cexlegend=1)
  {
  if(is.na(colscale[1]))
    {colscale=rev(gray.colors(nbscale))
    }
  facCodebook=cut(SOMresult$codebook[,specie],nbscale,ordered_result=T)
  colCodebook=colscale[as.numeric(facCodebook)]
  SOMshowColscale(colCodebook,coordHex=SOMresult$CoordHex,tit=specie,d=F,placeLegendRight=1)
  legend('right',fill=colscale,legend=minMaxScaleCut(as.character(levels(facCodebook)),2),cex=cexlegend)
  }

show.species=function(SOMresult,nbscale=10,colscale=NA,rowPerPage=3,colPerPage=3,PDF=F,nomPDF=ifelse(PDF,paste('SOM',Sys.Date(),'_ShowSp.pdf'),NA),sep='_')
  {
  if(!PDF)
    {x11(w=9,h=9)
    }else
    {pdf(nomPDF,9,9)
    }
  par(mfrow=c(rowPerPage,colPerPage))
  for(i in 1:ncol(SOMresult$codebook))
    {
    show.specie(colnames(SOMresult$codebook)[i],SOMresult,nbscale,colscale,cexlegend=2/colPerPage)
    if(i%%(rowPerPage*colPerPage)==0&&!PDF)
      {
      x11(w=9,h=9);par(mfrow=c(rowPerPage,colPerPage))
      }
    }
  if(PDF)dev.off()
  }

nbSUCell=function(SOMresult,plot=T,colEmpty='red')
  {justFilled=tapply(SOMresult$BMU1_2[,1],SOMresult$BMU1_2[,1],length)
  res=numeric(nrow(SOMresult$codebook))
  res[as.numeric(names(justFilled))]=justFilled
  if(plot)
    {facRes=cut(res,5,ordered_result=T)
    Col=rev(gray.colors(5))[facRes]
    if(!is.na(colEmpty))
      Col[res==0]=colEmpty
    SOMshowColscale(Col,SOMresult$CoordHex,tit='Sampling Units by Cell',placeLegendRight=1)
    if(is.na(colEmpty)){fillLegend=rev(gray.colors(5))}else{fillLegend=c(colEmpty,rev(gray.colors(5)))}
    if(is.na(colEmpty)){legendLegend=minMaxScaleCut(as.character(levels(facRes)),2)}else{legendLegend=c(
	0,minMaxScaleCut(as.character(levels(facRes)),2))}
    legend('right',fill=fillLegend,legend=legendLegend)
    }
  return(res)
  }

show.numCell=function(SOMresult)
  {SOMshowColscale('white',SOMresult$CoordHex)
  for(i in 1:nrow(SOMresult$CoordHex))
    {text(SOMresult$CoordHex[i,1],SOMresult$CoordHex[i,2],i)
    }
  }

SOMclassify=function(SOMresult,Dist="none",algohclust="none",plot=T,test=F,nomtest=NA)
  {
  #appel=match.call()
  #browser()
  if(Dist=="none")
    {Dist=SOMresult$distUse
    }else
    {Dist=match.arg(Dist,choices=c('euclidean','BrayCurtis','manhattan'))
    }
  dd=switch(Dist,
    "euclidean"=dist(SOMresult$codebook),
    'BrayCurtis'=vegdist(SOMresult$codebook,method='bray'),
    'manhattan'=vegdist(SOMresult$codebook,method='manhattan')
    )
  if(algohclust=="none")
    {algohclust=ifelse(Dist=="euclidean",'ward','average')
    }
  dendro=hclust(dd,method=algohclust)
  if(plot)
    {plot(dendro)
    }
	return(dendro)
  }


SOMclusterize=function(SOMresult,nbGroup,Dist="none",algohclust="none",plot=T,plotTOT=T,cexlegend=1)
  {if(Dist=="none")
    {Dist=SOMresult$distUse
    }else
    {Dist=match.arg(Dist,choices=c('euclidean','BrayCurtis','manhattan'))
    }
  dd=switch(Dist,
    "euclidean"=dist(SOMresult$codebook),
    'BrayCurtis'=vegdist(SOMresult$codebook,method='bray'),
    'manhattan'=vegdist(SOMresult$codebook,method='manhattan')
    )
  if(algohclust=="none")
    {algohclust=ifelse(Dist=="euclidean",'ward','average')
    }
  colscale=rainbow(nbGroup)
  dendro=hclust(dd,method=algohclust)
  if(plot&&plotTOT)
    {par(mfrow=c(1,2))
    plot(dendro)
    colClustRect(dendro,k=nbGroup)
    }
  gp=cutree(dendro,k=nbGroup)
  if(plot)
    {
    SOMshowColscale(col=colscale[gp],coordHex=SOMresult$CoordHex,tit="Clusters",placeLegendRight=ifelse(plotTOT,2,1))
    legend('right',fill=colscale,legend=levels(as.factor(gp)),cex=cexlegend)
    }
  return(gp)
  }

plot.SOM_gp=function(SOMresult,vecnbGroup=c(3,5),tit=paste(SOMresult$distUse,'[',paste(SOMresult$taille,collapse=' ')
,']'))
    {
	par(mfrow=c(length(vecnbGroup)+1,2),oma=c(1,2,1,1))
	show.numCell(SOMresult)
	nbSUCell(SOMresult)
	for(i in vecnbGroup)
	    {
		Dist=SOMresult$distUse
	      	dd=switch(Dist,
		    "euclidean"=dist(SOMresult$codebook),
		    'BrayCurtis'=vegdist(SOMresult$codebook,method='bray'),
		    'manhattan'=vegdist(SOMresult$codebook,method='manhattan')
		)
		algohclust=ifelse(Dist=="euclidean",'ward','average')
		dendro=hclust(dd,method=algohclust)
		plot(dendro)
		colClustRect(dendro,k=i,legend=F)
		SOMclusterize(SOMresult,nbGroup=i,plotTOT=F)
	    }
	title(main=tit,outer=T)
    }

clustDisjoint=function(SOMresult,nbGpes=2:(SOMresult$taille[1]*SOMresult$taille[2]-1),Dist="none",algohclust="none")
    {
  if(Dist=="none")
    {Dist=SOMresult$distUse
    }else
    {Dist=match.arg(Dist,choices=c('euclidean','BrayCurtis','manhattan'))
    }
  dd=switch(Dist,
    "euclidean"=dist(SOMresult$codebook),
    'BrayCurtis'=vegdist(SOMresult$codebook,method='bray'),
    'manhattan'=vegdist(SOMresult$codebook,method='manhattan')
    )
  if(algohclust=="none")
    {algohclust=ifelse(Dist=="euclidean",'ward','average')
    }
	#distGeo=round(as.matrix(dist(SOMresult$CoordHex)),3)
	dendro=hclust(dd,method=algohclust)
	tabclust=cutree(dendro,k=nbGpes)
	resComp=apply(tabclust,2,function(gps,coord)
	    {sapply(
		lapply(
		    by(coord,gps,FUN=dist)
		,as.matrix)
	    ,clustDisjoint1gp)
		},coord=as.data.frame(SOMresult$CoordHex))
	#browser()
	AnyPb_gp=(!sapply(resComp,all))
	AnyPb_tot=(any(AnyPb_gp))
	nbGpesPb=names(which(AnyPb_gp))
	WhichPb=lapply(resComp[nbGpesPb],function(logic)which(!logic))
	if(AnyPb_tot){res=list(AnyPb_tot=AnyPb_tot,AnyPb_gp=AnyPb_gp,nbGpesPb=nbGpesPb,WhichPb=WhichPb)
		}else{res=list(AnyPb_tot=AnyPb_tot)}
	if(AnyPb_tot){cat('Attention : cette SOM comprends des clusters disjoints \n')}
	return(res)
	
    }

clustDisjoint1gp=function(distHexGp)
    {
	#1 if ici nous permet de régler le problème des groupes avec 1 seul 
	cell1=1
	cellOK=c(cell1,which(round(distHexGp[,cell1],3)==1.000))
	nb1=0
	while(nb1!=length(cellOK)&&length(cellOK)!=nrow(distHexGp)&&length(cellOK)!=1)
	    {
		nb1=length(cellOK)
		search=apply(as.matrix(distHexGp[,cellOK]),2,function(cells){
		  which(round(cells,3)==1.000)})
		for(i in 1:length(search))
		    {
			cellOK=unique(c(cellOK,search[[i]]))
		    }

	    }
	res=length(cellOK)==nrow(distHexGp)
	return(res)

}
colClustRect=function(dendro,k,col=NULL,legend=T)
    {
	if(is.null(col)){col=rainbow(k)}
	clusts=cutree(dendro,k)[dendro$order]
	height=ifelse(k>2,mean(dendro$height[length(dendro$height)-(k-c(1,2))]),1.02*max(dendro$height))
	indMax=which(clusts[2:length(clusts)]!=clusts[2:length(clusts)-1])
	indMin=c(0.6,(indMax+1)-0.4)
	indMax=c(indMax,length(clusts))+0.4
	rect(indMin,rep(0,k),indMax,height,border=col[unique(clusts)])
	if(legend){legend('topright',legend=1:k,fill=col)}
    }

SOM_listeVR=function(SOMresult)#renvoie une liste ayant pour nom les noms de ligne du codebook et à l'interieur des
#éléments les noms des lignes du tableau
    {return(
	lapply(
	    by(rownames(SOMresult$BMU1_2),factor(SOMresult$BMU1_2[,1],levels=as.character(1:nrow(SOMresult$codebook))),
		function(vec)return(vec))
	    ,as.character)
	)
    }
showVRnames=function(SOMresult,nb=5,cexVR=0.8)
    {
	SOMshowColscale('white',SOMresult$CoordHex)
	namesVR=lapply(SOM_listeVR(SOMresult),function(vecChar){
		tailleSupNb=(length(vecChar)>nb)
		vecChar=na.omit(vecChar[1:nb])
		vecChar=paste(vecChar,collapse='\n')
		if(tailleSupNb){vecChar=paste(vecChar,'ETC.',sep='\n')}
		return(vecChar)
		})
	 for(i in 1:nrow(SOMresult$CoordHex))
	    {text(SOMresult$CoordHex[i,1],SOMresult$CoordHex[i,2],namesVR[[i]],cex=cexVR,col='red')
	    }
	invisible(namesVR)
    }

show.extVar=function(SOMresult,varExt,addClust=T,nbGp=5,ncolPage=3,nrowPage=3,PDF=F,
	nomPDF=paste('SOM',Sys.Date(),'_varExt.pdf'))
    {
	varExtScale=scale(as.matrix(varExt))
	VRnames=SOM_listeVR(SOMresult)
	vecApp=character(nrow(varExt))
	names(vecApp)=rownames(varExt)
	for(i in 1:length(VRnames)){vecApp[VRnames[[i]]]=names(VRnames)[i]}
	vecApp=factor(vecApp,levels=as.character(1:prod(SOMresult$taille)))
	meanVV=aggregate(varExt,list(vecApp),FUN=mean)
	rownames(meanVV)=meanVV[,1];meanVV=meanVV[,2:ncol(meanVV)]
	meanVV=meanVV[rownames(SOMresult[['codebook']]),];rownames(meanVV)=rownames(SOMresult[['codebook']])
	sdVV=aggregate(varExtScale,list(vecApp),FUN=sd)
	rownames(sdVV)=sdVV[,1];sdVV=sdVV[,2:ncol(sdVV)]
	sdVV=sdVV[rownames(SOMresult[['codebook']]),];rownames(sdVV)=rownames(SOMresult[['codebook']])
	tailHex=.5-sdVV/4
	tailHex[tailHex<0]=0;tailHex[is.na(tailHex)]=.5
	colHexGP=lapply(meanVV,cut,10,ordered_result=T)
	colHex=lapply(colHexGP,as.numeric)
	colHex=lapply(colHex,function(vec)rev(gray.colors(10))[vec])
	nbSU=nbSUCell(SOMresult,plot=F)
	colHex=lapply(colHex,function(vec,nb){vec[nb==0]='red';
		#vec[nb==1]='orange';
		return(vec)}
		,nb=nbSU)
	if(!PDF)
	    {x11(w=9,h=9)
	    }else
	    {pdf(nomPDF,9,9)
	    }
	par(mfrow=c(nrowPage,ncolPage))
	for(i in 1:length(colHex)){
		SOMshowColscale(col=colHex[[i]],coordHex=SOMresult$CoordHex,d=tailHex[,i],tit=names(colHex)[i],
			placeLegendRight=3)
		legend('right',fill=c('red',rev(gray.colors(10))),legend=c('empty',
			minMaxScaleCut(as.character(levels(colHexGP[[i]])),2)))
		if(any(nbSU==1)){for(j in which(nbSU==1))
			{hexagone(SOMresult$CoordHex[j,1],SOMresult$CoordHex[j,2],col=NA,d=.5,density=20)}}
		if(addClust){###ADDCLUST
			addclust(SOMresult,nbGp=nbGp)
		    }
		 if(i%%(nrowPage*ncolPage)==0&&!PDF)
		    {
			cat(i,'\n')
			x11(w=9,h=9);par(mfrow=c(nrowPage,ncolPage))
		    }
	    }
    }
addclust=function(SOMresult,nbGp)
{
clust=SOMclusterize(SOMresult,nbGp,plot=F)
colclust=rainbow(nbGp)
for(i in 1:nrow(SOMresult$CoordHex))
{
hexagone(SOMresult$CoordHex[i,1],SOMresult$CoordHex[i,2],d=.47,col=NA,border=colclust[clust[i]],lwd=1.5)
}
}
siteClust=function(SOMresult,gp1=F){
	if(gp1){
		tabClustCell=SOMclusterize(SOMresult,nbGroup=1:prod(SOMresult$taille),plot=F)
	}else{
		tabClustCell=SOMclusterize(SOMresult,nbGroup=2:prod(SOMresult$taille),plot=F)
	}
  tabClustSite=tabClustCell[as.character(SOMresult$BMU1_2[,1]),]
  rownames(tabClustSite)=rownames(SOMresult$BMU1_2)
  return(tabClustSite)
}

siteInClust=function(SOMresult,nbGroup){
	clustCell=SOMclusterize(SOMresult,nbGroup=nbGroup,plot=F)
	cellSite=as.character(SOMresult$BMU1_2[,1]);names(cellSite)=rownames(SOMresult$BMU1_2)
	clustSite=clustCell[cellSite];names(clustSite)=names(cellSite)
	return(clustSite)
}

tailleSOM=function(tab)
{
	tailleOpt=5*nrow(tab)^0.54321#calcul de la taille optimale
	rowcolpos=combn(1:(2*sqrt(tailleOpt)),2)#
	taillePossibles=apply(rowcolpos,2,prod)
	rowcolpos=apply(rowcolpos,2,range)
	rowcolpos=rowcolpos[,which(taillePossibles>(.7*tailleOpt)&taillePossibles<(1.2*tailleOpt))]
	taillePossibles=apply(rowcolpos,2,prod)
	rowcolpos=rowcolpos[,apply(rowcolpos,2,function(x){x[2]<=1.8*x[1]})]
	rowcolpos=rowcolpos[,apply(rowcolpos,2,function(x){x[2]>=round(1.3*x[1])})]
	return(rowcolpos)
}
testTailleSOM=function(tab,nb=5,...)
{
	listTaille=apply(tailleSOM(tab),2,as.list)
		res=lapply(listTaille,function(taille,Nb,data,...){
			listSom=list()
			for(i in 1:Nb)
			{
				listSom[[i]]=SOM(data,taille[[2]],taille[[1]],...)
			}
			return(listSom)
		},Nb=nb,data=tab,...)
	return(res)
}
