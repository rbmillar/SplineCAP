#' Does post-processing of canonical analysis object from vegan to extract spline gradients
#'
#' @param object The fitted rda/cca or capscale object
#' @param Bdf Integer vector giving the degrees of freedom (number of columns) of the bases used.
#' @param Gnames Gradient names corresponding to the bases.
#' @param predB List containing the new set of bases for prediction.
#' @param plot Draw plot of the gradient. predB must be supplied.
#' @param ... Arguments to be passed on to the plot command.
#'
#' @return Gradt Fitted gradient values
#' @return Env.scores Correlation scores for spline and linear covariates, for adding to a biplot
#' @return pred.Gradt List containing the fitted gradient values for predB
#'
#' @export
#'
#' @examples
#'data(PitTraps)
#'Locn=PitTraps$Locn
#'distance=PitTraps$distance
#'Counts=PitTraps[,-c(1,2)]

#'B=bSpline(distance,degree=2,df=5)   #Spline basis for distance
#'Bpred=predict(B,newx=min(distance):max(distance)) #For prediction
#'CAP.fit=capscale((Counts>0)~B+Locn,dist="bray") #Fit RDA on PCOs
#'Gradients=splineGradient(CAP.fit,Bdf=5,predB=list(Bpred),type="l")
#'
#'#Draw triplot
#'env=Gradients$Env.scores
#'ordiplot(CAP.fit,display = c("sp", "wa", "cn")[1:2],scaling=3,
#'             cex=0.6,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
#'arrows(0,0,env[1,],env[2,],length=0.05,col="blue",lwd=1.5)
#'text(env[1,],env[2,],colnames(env),cex=0.8,adj=c(0,0.5))
#'
#splineGradient function

splineGradient=function(object,Bdf,Gnames=NULL,predB=NULL,plot=T,...)
  {
  n.B=length(Bdf)
  B.index=list()
  B.index[[1]]=1:Bdf[1]
  cat("*** A basis is assumed to be in columns 1 to", Bdf[1])
  if(n.B>=2) {
    for(b in 2:n.B) B.index[[b]]=max(B.index[[b-1]])+c(1,Bdf[b])
    for(b in 2:n.B) cat(" and columns", min(B.index[[b]]),"to",max(B.index[[b]]))
  }
  cat("\n")

  #Propn of variability in Yhat explained by CAP axes
  CAP.varPropn=round( summary(object)$cont$importance[2,1:2] ,4)

  #Recover the X matrix of covariates (incl bases)
  X=qr.X(object$CCA$QR)
  BetaU=coef(object) #The covariate (incl bases) scores
  Z=X%*%BetaU #Fitted values in covariate space

  n.sites=nrow(X)
  Gradt=array(NA,c(n.sites,2,n.B));
  for(b in 1:n.B) {
    index=B.index[[b]]
    Gradt[,,b]=X[,index]%*%(BetaU[index,1:2]) #NB: Only using first two gradients
    colnames(Gradt[,,b])=paste0(Gnames[b],1:2)
  }

  #Calculate gradient and linear covariate correlation scores
  n.Bcols=sum(Bdf)
  X.lin=subset(X,select=-(1:n.Bcols)) #Linear environmental terms
  n.lincols=ncol(X.lin) #Number of linear covariate terms
  cor.env=matrix(NA,2,2*n.B+n.lincols)
  XG=X.lin #For holding spline and linear terms
  for(b in n.B:1) XG=cbind(Gradt[,,b],XG)
  cor.env=cor(Z[,1:2],XG)
  if(is.null(Gnames)) Gnames=LETTERS[1:n.B]
  colnames(cor.env)=c(paste0(rep(Gnames,rep(2,n.B)),rep(1:2,n.B)),colnames(X.lin))

  return.list=list(PropnVarExplained=CAP.varPropn,Gradt=Gradt,Env.scores=cor.env)

  if(!is.null(predB)) { #Prediction requested for new basis
    pred.Gradt=list()
    for(b in 1:n.B) {
      if(!is.null(predB[[b]])) {
        pred.Gradt[[b]]=predB[[b]]%*%(BetaU[B.index[[b]],1:2])
        x=attributes(predB[[b]])$x
        if(plot) for(pco in 1:2)
          plot(x,pred.Gradt[[b]][,pco],ylab=paste("Spline gradient",pco),...)
      }
      else pred.Gradt[[b]]=NULL
    }
    return.list=c(return.list,pred.Gradt=pred.Gradt)
  }
  return(return.list)
}


#' Calculates the leave-one-out prediction sums of squares of a multivariate linear regression
#'
#' @param f.object Model formula
#'
#' @return A list containing residual sums of squares and leave-one-out prediction sums of squares
#'
#' @export
#'
#' @examples
#'data(PitTraps)
#'Locn=PitTraps$Locn
#'distance=PitTraps$distance
#'Counts=PitTraps[,-c(1,2)]
#'
#'B=bSpline(distance,degree=2,df=5)   #Spline basis for distance
#'First need to calculate PCOs, as these are the Y variables
#'D=as.matrix(vegdist(Counts>0)) #Bray-Curtis is default in vegdist
#'Y.pco=cmdscale(D,k=nrow(D)-1,eig=T)$points #The PCOs
#'LOO.prednError(Y.pco~B+Locn)
#'
#LOO.prednError function

LOO.prednError=function(f.object) {
  #fobj=reformulate(attributes(terms(f.object))$term.labels,intercept=F)
  #X=model.matrix(fobj)
  Within.mlm=lm(f.object,x=T,y=T)
  RSS=sum(resid(Within.mlm)^2)
  cat("RSS is",RSS)
  X=Within.mlm$x
  if(class(X)[1]=="numeric") X=as.matrix(X,ncol=1)
  Y=Within.mlm$y
  if(class(Y)[1]=="numeric") Y=as.matrix(Y,ncol=1)
  colnames(X)=paste0("x",1:ncol(X))
  colnames(Y)=paste0("y",1:ncol(Y))
  n.sites=nrow(X)
  predSqErr=rep(NA,n.sites)
  #Df=as.data.frame(cbind(X,Y))
  for(i in 1:n.sites) {
    XLOO=X[-i,]
    YLOO=Y[-i,]
    #Remark: Don't need to centre XLOO, since lm fits intercept
    #predict doesn't work if use lm(YLOO~XLOO)
    if(ncol(X)>=1) multi.lm=lm(YLOO~.-1,data=data.frame(X=XLOO))
    YihatLOO=predict(multi.lm,newdata=data.frame(X=t(X[i,])))
    predSqErr[i]=sum((Y[i,]-YihatLOO)^2) }
  cat("\nLOO predn sum of squares is", sum(predSqErr),"\n")
  return(invisible(list(RSS=RSS,LOO=sum(predSqErr))))
}




