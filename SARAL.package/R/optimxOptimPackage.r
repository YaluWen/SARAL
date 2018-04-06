### This is the required function ####
#library(optimx);
#library(psych)
#library(lars);
#library(mvtnorm);
#library(expm)
#library(msgps);
#library(rrBLUP)
#library(matrixcalc)
#library(nloptr)
## Likelihood ##
neg.likelihood.SAR.MLE.sigma<-function(para,env)
{
  Y=get("Y",env)
  W=get("W",env)
  theta=para#[1:(length(para)-1)]
  n=length(Y)
  In=diag(length(Y))
  C_theta=theta[1]*W[[1]];if(length(W)>1) {for(i in 2:length(W)) C_theta=C_theta+theta[i]*W[[i]]}
  tmp=In-C_theta
  sigma=sqrt(t(Y-mean(Y)) %*% t(tmp) %*% tmp %*% (Y-mean(Y))/n)  ## At MLE
  eigen.value=eigen(t(tmp) %*% tmp)$value
  eigen.value[abs(eigen.value)<1e-15]=abs(eigen.value[abs(eigen.value)<1e-15]);
  likeli=-n/2*log(sigma^2)+1/2*sum(log(eigen.value))-n/2
  return(as.numeric(-likeli))
}

# first derivatives #
neg.firstderi.SAR.MLE.sigma<-function(para,env)
{
  Y=get("Y",env)
  W=get("W",env)
  result=NULL;
  theta=para#[1:(length(para)-1)]
  #sigma=para[length(para)]
  n=length(Y)
  In=diag(length(Y))
  C_theta=theta[1]*W[[1]];if(length(W)>1) {for(i in 2:length(W)) C_theta=C_theta+theta[i]*W[[i]]}
  tmp=In-C_theta
  sigma=sqrt(t(Y-mean(Y)) %*% t(tmp) %*% tmp %*% (Y-mean(Y))/n)  ## At MLE
  tmp.inv=solve(tmp)
  for(i in 1:length(W))
    result=rbind(result,-psych::tr(tmp.inv %*% W[[i]])+1/sigma^2*t(Y-mean(Y)) %*% tmp %*% W[[i]] %*% (Y-mean(Y)))
  #result=rbind(result,-n/sigma+1/sigma^3*t(Y-mean(Y))%*% tmp %*% tmp %*% (Y-mean(Y)))
  -result
}

## Second Derivatives ##
neg.secondderi.SAR.MLE.sigma<-function(para,env)
{
  Y=get("Y",env)
  W=get("W",env)
  theta=para
  n=length(Y)
  In=diag(length(Y))
  C_theta=theta[1]*W[[1]];if(length(W)>1) {for(i in 2:length(W)) C_theta=C_theta+theta[i]*W[[i]]}
  tmp=In-C_theta
  sigma=sqrt(t(Y-mean(Y)) %*% t(tmp) %*% tmp %*% (Y-mean(Y))/n)  ## At MLE
  tmp.inv=solve(tmp)
  result=matrix(NA,length(para),length(para));
  for(i in 1:length(theta))
    for(j in i:length(theta))
    {
      result[i,j]=result[j,i]=-psych::tr(tmp.inv %*% W[[j]] %*% tmp.inv %*% W[[i]])-1/sigma^2*t(Y-mean(Y)) %*% W[[j]] %*% W[[i]] %*% (Y-mean(Y))
    }
  -result
}

## MML scores ##
mml.dfgps<-function (dfgps_result, intercept = FALSE, stand.coef = FALSE)
{
  N <- dfgps_result$N;
  tau2<-dfgps_result$RSS/N;
  y <- dfgps_result$y;
  y <- y-mean(y)
  k <- dfgps_result$df
  mml_vec <- N * log(2 * pi * tau2) + N + (k+1)*log(sum(y^2)/2/tau2)+log(k+1)-2*gamma((k+3)/2);
  mml_index <- which.min(mml_vec)
  mml_index <- as.integer(mml_index)
  p = dfgps_result$p
  p = as.integer(p)
  delta_t = dfgps_result$delta_t
  delta_t = as.double(delta_t)
  betavec = as.integer(dfgps_result$coefficient_index)
  tuning <- dfgps_result$tuning[mml_index]
  if (stand.coef == FALSE)
    tuning <- dfgps_result$tuning_stand[mml_index]
  beta_mml = msgps::coef.dfgps(dfgps_result, tuning, intercept, stand.coef)$coefficient
  df <- dfgps_result$df[mml_index]
  return(list(result = mml_vec, step = mml_index, tuning = tuning,
              coef = beta_mml, df = df))
}


## optimization ##

SARpredassomain<-function(data,compare='rrblup',assess=TRUE,methods=c("Cp","AICc","GCV","BIC","MML"),max.step=20)
{

  Y=data$train$Y
  W=data$train$IBS
  train.index=data$train.index;

  maxvalue=NULL;
  for(i in 1:length(W))
  {
    tmp=W[[i]];
    maxvalue=c(maxvalue,max(eigen(tmp)$values));
  }
  maxvalue=max(maxvalue);

  ## rescale by sqrt ##
  maxvalue1=maxvalue;
  maxvalue=NULL;
  for(i in 1:length(W))
  {
    data$train$IBS[[i]]=data$train$IBS[[i]]/sqrt(maxvalue1);
    data$IBS[[i]]=data$IBS[[i]]/sqrt(maxvalue1);
    tmp=data$train$IBS[[i]];
    maxvalue=c(maxvalue,max(eigen(tmp)$values));
  }
  maxvalue=max(maxvalue);
  W=data$train$IBS

  env=new.env();
  assign("Y",Y,env);
  assign("W",W,env);


  ## initial estimates ##
  converge=TRUE;
  set.seed(10);
  int1=runif(length(W),-1/maxvalue,1/maxvalue);
  opt1 <- optimx::optimx(int1,neg.likelihood.SAR.MLE.sigma,gr=neg.firstderi.SAR.MLE.sigma, hessian = TRUE,env=env,method="L-BFGS-B",itnmax=5000)
  second<-attr(opt1,"details")[1,]$nhatend;

  if(!opt1$kkt2 | is.na(opt1$kkt2))
  {
    step=1;
    while((!(opt1$kkt2)| is.na(opt1$kkt2)) & step < max.step)
    {
      int1=runif(length(W),-1/maxvalue,1/maxvalue);
      #opt1 <- optimx::optimx(int1,neg.likelihood.SAR.MLE.sigma,hessian = TRUE,env=env,method="L-BFGS-B")
      opt1 <- optimx::optimx(int1,neg.likelihood.SAR.MLE.sigma,gr=neg.firstderi.SAR.MLE.sigma, hessian = TRUE,env=env,method="L-BFGS-B",itnmax=5000)

      second<-attr(opt1,"details")[1,]$nhatend;
      step=step+1;
    }
    if(!(opt1$kkt2) & step == max.step)
    {
      opt1[1,1:length(get("W",env))]=0; # no variable would be selected;
      second=neg.secondderi.SAR.MLE.sigma(unlist(opt1[1,1:length(get("W",env))]),env=env);
      converge=FALSE;
    }
  }
  B=expm::sqrtm(second);

  ## The below is finding the penalized solution ###
  #  if(length(get("W",env))==1) B=sqrt(second[1,1])
  if(converge)
  {
    MLE=matrix(unlist(opt1[1,1:length(get("W",env))]),ncol=1)
    Y.lla=B %*% MLE
    X.lla=B
    for(i in 1:ncol(X.lla))
      X.lla[,i]=X.lla[,i]*(MLE[i,1])
    fit4 <- msgps::msgps(X.lla,Y.lla[,1],penalty="enet",alpha=0,intercept=FALSE)

    ## Getting coefficients and tests##
    COEFS=TESTS=NULL;
    vars=solve(second);
    for(i in 1:length(methods))
    {
      if(methods[i]!="MML")
      {
        coefs=coef(fit4)[,methods[i]]*MLE ;
        tests=coefs/sqrt(diag(vars));
      }
      if(methods[i]=="MML")
      {
        tmp=mml.dfgps(fit4$dfgps_result, intercept = fit4$intercept, stand.coef = fit4$stand.coef)
        coefs=tmp$coef*MLE ;
        tests=coefs/sqrt(diag(vars));
      }
      COEFS=rbind(COEFS,matrix(coefs,nrow=1));
      TESTS=rbind(TESTS,matrix(tests,nrow=1));
    }
    rownames(TESTS)=rownames(COEFS)=methods;

    ## Getting Predictions ##
    Wall=data$IBS;
    train.index=data$train.index;
    In=diag(nrow(Wall[[1]]));
    PRED=NULL;
    ASSESS=NULL;
    for(me in 1:length(methods))
    {
      predcoef=COEFS[me,];
      Ctheta=predcoef[1]*Wall[[1]];if(length(Wall)>1){for(i in 2:length(Wall)) Ctheta=Ctheta+predcoef[i]*Wall[[i]]}

      Q=t(In-Ctheta) %*% (In-Ctheta)
      VARS=solve(Q)
      pred2=mean(Y)+VARS[-train.index,train.index] %*% solve(VARS[train.index,train.index]) %*% (Y-mean(Y))
      pred=matrix(pred2,ncol=1)
      colnames(pred)=paste(methods[me],c("pred"),sep="_");
      PRED=cbind(PRED,pred);
      if(assess)
      {
        COR=cor(data$Y[-data$train.index],pred2,use="na.or.complete");
        MSE=mean((data$Y[-data$train.index]-pred2)^2, na.rm = TRUE);
        tmp=matrix(c(COR,MSE),nrow=1);
        colnames(tmp)= c("COR","MSE");
        rownames(tmp)=methods[me];
        ASSESS=rbind(ASSESS,tmp);
      }
    }
    if(assess)
    {
      PRED=cbind("YTrue"=data$Y[-data$train.index],PRED)
    }

    if(tolower(compare)=='rrblup')
    {
      testid=(1:length(data$Y))[-train.index];
      M=data$Geno
      rownames(M) <- 1:nrow(data$Geno)
      A <- rrBLUP::A.mat(M)
      Y1=data$Y
      Y1[testid]=NA
      datablup <- data.frame(y=Y1,gid=1:nrow(data$Geno))
      #predict breeding values
      ans <- rrBLUP::kin.blup(data=datablup,geno="gid",pheno="y",K=A);
      PRED=cbind(PRED,"rrBLUP"=ans$g[testid]+mean(Y1,na.rm=T));

      if(assess)
      {
        CORrrblup <- cor(data$Y[testid],ans$g[testid],use="na.or.complete");
        MSErrblup= mean((data$Y[testid]-(ans$g[testid]+mean(Y1,na.rm=T)))^2, na.rm = TRUE);
        tmp=matrix(c(CORrrblup,MSErrblup),nrow=1);
        colnames(tmp)= colnames(ASSESS)
        rownames(tmp)="rrBLUP";
        ASSESS=rbind(ASSESS,tmp);
      }
    }
    Result=list();
    Result$Pred=PRED;
    Result$Assess=ASSESS;
    #Result$Coefs=COEFS;
    #Result$Tests=TESTS;
    #Result$opt1=opt1;
  }
  if(!converge)
  {
    Result=list();
    Result$Pred=list();
    Result$Assess=list();
    #Result$Coefs=list();
    #Result$Tests=list();
    #Result$opt1=opt1;
    Result$Errors="MLE does not converge! Try new initial values!"
  }
  Result;
}

#' @title Genetic Risk Prediction Using a Spatial Autoregressive Model with Adaptive Lasso
#'
#' @description
#' \code{SARpred} returns the predicted values of the phenotypes.
#' @details
#' This function is used to predict outcomes based on individual's genetic profiles.
#'
#' @param Y A vector of phenotypes. Missing values should be coded as \code{NA}.
#' @param Geno A matrix of genotypes, with each row representing individual and column representing single nucleotide variant.
#' @param TrainID A bool vector with \code{true} indicating the individual is in the training sample.
#' @param genelist A list with each element represents a group.
#' @param Similarity A list with \code{i}-th element measuring the similarity for a particular chuck of genome listed in the \code{genelist[[i]]}.
#' @param SaveSimilarity A bool indicating whether the similarity should be saved for futher use. By default, similarity is not saved.
#' @param compare A string that indicates which method the SARAL compares to. Currently it only supports GBLUP.
#' @param assess A bool. If true, the method will calculate Pearson correlation and mean square error from the testing data. This option is ignored if all the phenotypes of the testing sample is \code{NA} (i.e. \code{sum(!is.na(Y[!TrainID]))==0}).
#' @param methods A string vector that determines which model selection criterion is used to select model for SARAL. By default, it uses Mallows' Cp (Cp), bias corrected AIC (AICc), generalized cross validataion (GCV), BIC, and information theoretic minimum message length (MML).
#' @return A list that contains the predicted values, the pearson correlations and mean square errors (if \code{assess=TRUE}), and error if algorithm does not converge. The list also contains a list of similarity measure (if \code{SaveSimilarity=TRUE}).
#' @examples
#' # This may take 2 hour, don't run if not needed #
#' Y=SARAL::example_SARAL$Y;
#' Geno=SARAL::example_SARAL$Geno;
#' TrainID=SARAL::example_SARAL$TrainID;
#' genelist=SARAL::example_SARAL$Genelist;
#' Result<-SARALpred(Y,Geno,TrainID,genelist,compare='rrblup',assess=TRUE,methods=c("Cp","AICc","GCV","BIC","MML"))
#' @export
SARALpred<-function(Y, Geno, TrainID,genelist=list(),Similarity=NULL,SaveSimilarity=FALSE,YTest=NULL,compare='rrblup',assess=TRUE,methods=c("Cp","AICc","GCV","BIC","MML"))
{
  YTrain=Y[TrainID];
  YTest=Y[!TrainID];
  if(sum(!is.na(YTest))==0 & assess==TRUE)
  {
    warning("Can't assess the prediction performance due to the lack of true values!")
    assess=FALSE;
  }

  ## checking genotypes for both training and testing ##
  if(is.null(colnames(Geno)))
  {
    stop("Genotypes do not have their names!")
  }
  GenoTrain=Geno[TrainID,];
  GenoTest=Geno[!TrainID,];

  ## delete indiviudals in the training set without phenotypes measured ##
  tmp=is.na(YTrain);
  if(sum(tmp)!=0){warning("Individuals without phenotypes in the training set are not used for training the prediction model!")}
  YTrain=YTrain[!tmp];
  GenoTrain=GenoTrain[!tmp,]

  data=list();
  Y=c(YTrain,YTest);
  Geno=rbind(GenoTrain,GenoTest);
  train.index=1:length(YTrain);
  data$Y=Y;
  data$train.index=train.index;
  data$Geno=Geno;
  data$Gene=genelist;
  IBS=list();IBSTrain=list();
  if(length(genelist)==0 & is.null(Similarity)){
    warning("no annotation is provided, and all the SNPs will be treated as if they come from the same gene");
    IBS[[1]]=as.matrix(dist(as.matrix(data$Geno)))
    IBSTrain[[1]]=IBS[[1]][train.index,train.index];
  }
  if(length(genelist)>0 & is.null(Similarity))
  {
    for(i in 1:length(genelist))
    {
      cat("calculating similartiy for the",i,"th gene!","\n");
      geno=data$Geno[,colnames(data$Geno) %in% genelist[[i]]];
      IBS[[i]]=as.matrix(dist(as.matrix(geno)));
      IBSTrain[[i]]=IBS[[i]][train.index,train.index];
    }
  }
  if(!is.null(Similarity))
  {
    if(length(genelist)!=length(Similarity))
    {
      stop("Similarity should be calculated according the genelist. The number of chucks(i.e. the length of the genelist) and the number of similarities are different! Please check!");
    }
    # checking the dimension of similarities
    for(i in 1:length(Similarity))
    {
      if((length(data$Y)!=nrow(Similarity[[i]])) | (length(data$Y)!=ncol(Similarity[[i]])))
      {
        stop("Similarity should be calculated for each pair of individuals, and thus it should be square matrix with the dimension equals to the total sample size! Please check!");
      }
    }
    IBS=Similarity;
    for(i in 1:length(genelist))
    {
      IBSTrain[[i]]=IBS[[i]][train.index,train.index];
    }
  }
  data$IBS=IBS;
  data$train=list();
  data$train$Y=YTrain;
  data$train$Geno=GenoTrain;
  data$train$IBS=IBSTrain;
  cat("Building Risk Prediction Models...");
  Result<-try(SARpredassomain(data,compare=compare,assess=assess,methods=methods),silent=TRUE);
  if(SaveSimilarity) Result$IBS=IBS;
  Result;
}
