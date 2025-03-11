######## FUNCTION mBiplext ####################
# Compositional biplot, colored by population,
# with additional real variables incorporated.
# Programmed by J.J. Egozcue (2020) based on
# previous function mBPPOP (Dec. 2014)
##### draws a CoDa-biplot with data coming from
# populations coded by a number xpop indicating
# color from a color sequence in colist.
# carries out clr of data set.
# centres the clr and the additional real variables
# in extr. Then they are added to the centered clr.
# carries out svd of centred clr and added real
# extra variables in extr
# plots biplot with data colored by xpop (number of color)
##### input:
# x compositional data by rows (matrix)
# xpop factor indicating the population of each row in x
# extr a vector or a matrix of real variables to be added
#    to the biplot.
# namextr name of a variable (only for single variable added)
# biscale    = 1 covariance biplot
#            = 0 form biplot
# circ       = FALSE (default) ; TRUE plots a unit circle
#              in form biplots
# punch      = 0 do not plot data points
#            = 1 plot symbols for data points
#              following pchpoint
#            = 2 plot numbers for data points following
#              pchpoints
#            = 3 lines between consecutive data points
#              are drawn
# choice[1:2] = PC's to be plotted, eg c(1,2), c(1,3)...
# colist a color sequence for populations.
# pchpoint integer sequence determining the plot pch symbol
#     of each point; when punch=2 is the number to be plotted.
# optpdf = 1  prints on pdf
# filename (optional) defines the name of the output pdf file.
#         by default formBPPOP.pdf or covBPPOP.pdf are used
#### output: a list containing
# the svd matrices: U, V, and singular values in D
# explained variance in the biplot explvar
# additionally
# a pdf file containing the biplot is printed in working dir.
###################################################

mBiplext <- function(x,xpop=NULL,extr=NULL,choice=c(1,2),
                     biscale=1,punch=1,colist=1:10,
                     circ=FALSE,colcirc="grey70",
                     optpdf=0,filename=NULL,
                     namextr=c("tot"),
                     colray="red",colextr="darkgreen",
                     cextext=1,lwdray=1,pchpoint=1){
  # point colors
  colpoint = rep(1,length(colist))
  if(!is.null(xpop)){
    colpoint = colist
  }
  # clr of x
  logx=log(x)
  xclr = logx - outer(rowMeans(logx),rep(1,length=ncol(logx)))
  # centring xclr
  cxclr = xclr - outer(rep(1,nrow(xclr)), colMeans(xclr))  
  # centering real variables extr, if any
  extrue=FALSE
  nextr=NULL
  if(!is.null(extr)){
    if(is.vector(extr)){
      cextr = extr - mean(extr)
      mextr = matrix(cextr,nrow=length(extr),ncol=1)
      colnames(mextr)=namextr
      nextr=1
      extrue = TRUE
    }
    if(is.matrix(extr)){
      namextr = colnames(extr)
      mextr = extr-outer(rep(1,nrow(extr)),colMeans(extr))
      nextr=ncol(mextr) }
    # append real variables in extr  
    cxclr1= cbind(cxclr,mextr)
    colnames(cxclr1)=c(colnames(xclr),namextr)
    cxclr=cxclr1
    extrue = TRUE
  }
  # svd (cxclr)
  SVDxclr = svd(cxclr)
  U = SVDxclr$u
  V = SVDxclr$v
  rownames(V)=colnames(cxclr)
  D = SVDxclr$d
  # scores and loadings
  ## covariance biplot
  if(biscale==1){
    ld=t(diag(D)%*%t(V))/sqrt(nrow(cxclr))
    mainT="covariance biplot"
    fileT="covBiplext"
    mld=max(abs(ld[,choice]))
    msc=max(abs(U[,choice]))
    sc=U *(mld/msc)  # scaling scores
    xylimit=c(-mld,mld)
  }
  
  ## form biplot
  ## scaling: unit norm of V-vectors
  if(biscale==0){
    sc=U%*%diag(D)
    ld=V
    mainT="form biplot"
    fileT="formBiplext"
    mld = max(abs(ld[,choice]))     # scaling basis vectors
    msc = max(abs(sc[,choice]))
    sc = sc*(mld/msc)
    xylimit = c(-mld,mld)
  }
  
  # numeric output
  variances = D^2
  totvar=sum(variances)/(nrow(x)-1)
  extrvar=0
  if(extrue==TRUE){
    extrvar = var(extr) }
  explvar = (variances[choice[1]]+variances[choice[2]])/sum(variances)
  # names
  #  clrnames=paste("clr.",colnames(x),sep="")
  clrnames=colnames(cxclr)
  if(choice[1] == 1){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste("first axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[1] == 2){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste("second axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[1] == 3){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste("third axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[1] >= 4){
    expl=100*variances[choice[1]]/sum(variances)
    xlabel=paste(paste(choice[1],"th axis",sep=""), " var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] == 1){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste("first axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] == 2){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste("second axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] == 3){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste("third axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
  if(choice[2] >= 4){
    expl=100*variances[choice[2]]/sum(variances)
    ylabel=paste(paste(choice[2],"th axis",sep=""), " var% ",format(expl,nsmall=1,digits=3),sep="")}
  
  if(punch==0){pun="n"}
  if(punch==1){pun="p"}
  if(punch==2){pun="n"}
  if(punch==3){pun="b"}
  
  # pdf output
  filenam=paste(fileT,".pdf",sep="")
  if(optpdf==1){
    if(is.null(filename)==FALSE){filenam=filename}
    pdf(filenam, width=5, height=5, fam="Times")}
  
  plot(sc[,choice],col=colist,type=pun,cex=0.8,asp=1,
       xlim=xylimit,ylim=xylimit,main=mainT,
       xlab=xlabel,ylab=ylabel,pch=pchpoint)
  # only form biplot: unit circle on variables
  if(circ==TRUE & biscale==0){
    theta = seq(from=0, to=(2*pi),length=150)
    xc=cos(theta)
    yc=sin(theta)
    lines(xc,yc,col="grey70")
  }
  
  #  this is for changing punch in place of color
  #  plot(sc[,choice],col="black",type=pun,cex=0.8,asp=1,
  #       pch=colpoint,
  #       xlim=xylimit,ylim=xylimit,main=mainT,
  #       xlab=xlabel,ylab=ylabel)
  
  if(punch==2){
    #    text(sc[,choice],labels=(1:nrow(sc)),col=colpoint,cex=0.8)       
    text(sc[,choice],labels=pchpoint,col=colpoint,cex=0.8)
  }
  for(i in 1:ncol(x)){
    xx=rbind(c(0,0),ld[i,choice])
    lines(xx, col=colray,lwd=lwdray)
    xtext = ld[i,choice[1]]
    ytext = ld[i,choice[2]]
    text(xtext,ytext,labels=clrnames[i],col=colray,
         pos=2,offset=0.3,cex=cextext)
  }
  if(!is.null(nextr)){
    for(iex in 1:nextr){
      nnrow = ncol(x)+iex
      xxetr = rbind(c(0,0),ld[nnrow,choice])
      lines(xxetr,col=colextr,lwd=lwdray)
      xtextr = ld[nnrow,choice[1]]
      ytextr = ld[nnrow,choice[2]]
      text(xtextr,ytextr,labels=clrnames[nnrow],col=colextr,
           pos=2,offset=0.3,cex=cextext)
    } }
  
  if(optpdf==1){
    dev.off()
  }
  lout = list("U"=U,"V"=V,"D"=D,"explvar"=explvar,"totvar"=totvar,
              "extrvar"=extrvar)
  return(lout)
}


