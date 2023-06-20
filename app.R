################
### Packages ###
################

library(fda)
library(shinythemes)
library(shiny)
library(shinyWidgets)
library(caret)
library(e1071)
library(fda.usc)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(shinydashboard)
library(DescTools)
library(knitr)
library(kableExtra)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)


#####################
### data function ###
#####################

# We create a function called data to select the dataset we want to analyze. There are
# 13 simulated data and 1 real data set (growth data)
data <- function(dataset){
  color_1 <- "deepskyblue2"
  color_2 <- "deeppink"
  color_3 <- "green"
  if(dataset == "growth"){
    attach(growth)
    X1 <- t(growth$hgtf)
    X2 <- t(growth$hgtm)
    t <- growth$age
    nf <- ncol(growth$hgtf)
    nm <- ncol(growth$hgtm)
    rangeval = c(1,18)
    truelabels <- c(rep(1,nf),rep(0,nm)) # 1 are girls and 0 are boys
    data <- rbind(X1,X2)
    xlab <- "Age (years)"
    ylab <- "Height (cm)"
    main="Growth of boys (blue) and girls (pink)"
    colors <- c(color_1, color_2)
    pch <- c(16, 17)
  }else{
    if(dataset>=2 & dataset<=9){
      n <- 50
      rangeval = c(0,1)
      t <- seq(0,1,length=30)
      truelabels <- c(rep(1,n), rep(0,n)) # 1 is pink and 0 blue
      xlab <- "t"
      ylab <- "Value"
      colors <- c(color_1,color_2)
      pch <- c(16, 17)
      X <- matrix(0,n,length(t))
      sigma <- matrix(0,length(t),length(t))
      
      # Building sigma
      for (i in 1:length(t)){
        for (j in 1:length(t)){
          sigma[i,j] <- 0.3*exp((-1/0.3)*abs(t[j]-t[i]));
        }
      }
      
      #centered Gaussian process
      egauss <- mvrnorm(n ,rep(0, length(t)), sigma) #a sample from the specified multivariate normal distribution
      
      sigma2 <- matrix(0,length(t),length(t))
      
      # Construimos sigma
      for (i in 1:length(t)){
        for (j in 1:length(t)){
          sigma2[i,j] <- 0.5*exp((-1/0.2)*abs(t[j]-t[i]));
        }
      }
      
      #centered Gaussian process
      hgauss <- mvrnorm(n ,rep(0, length(t)), sigma) #a sample from the specified multivariate normal distribution
      
      
      ## MODEL 1
      
      X1 <- X
      for (i in 1:n){
        for (j in 1:length(t)){
          X1[i,j] <- 30*t[j]^(3/2)*(1-t[j])+egauss[i,j]
        }
      }
      
      if(dataset==2){  ## MODEL 2
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+0.5+egauss[i,j]
          }
        }
      }
      else if(dataset==3){  ## MODEL 3
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+0.75+egauss[i,j]
          }
        }
      }
      else if(dataset==4){  ## MODEL 4
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+1.+egauss[i,j]
          }
        }
      }
      else if(dataset==5){  ## MODEL 5
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+2*egauss[i,j]
          }
        }
      }
      else if(dataset==6){  ## MODEL 6
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+0.25*egauss[i,j]
          }
        }
      }
      else if(dataset==7){  ## MODEL 7
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]^(3/2)*(1-t[j])+hgauss[i,j]
          }
        }
      }
      else if(dataset==8){  ## MODEL 8
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]*(1-t[j])+hgauss[i,j]
          }
        }
      }
      else if(dataset==9){  ## MODEL 9
        X2 <- X
        for (i in 1:n){
          for (j in 1:length(t)){
            X2[i,j] <- 30*t[j]*(1-t[j])+egauss[i,j]
          }
        }
      }
      data <- rbind(X1, X2)
      main=paste("Model 1 (pink) & Model ", dataset, " (blue)")
    }
    if(dataset == 11| dataset == 12){
      n <- 50
      rangeval = c(0,1)
      t <- seq(0,1,length.out=150)
      truelabels <- c(rep(1,n), rep(0,n)) # 1 is pink and 0 blue
      xlab <- "t"
      ylab <- "Value"
      colors <- c(color_1, color_2)
      pch <- c(16, 17)
      P <- 150 #equidistant points in the interval I=[0,1]
      K <- 100 #components
      m1 <- t*(1-t)
      
      rho <- rep(0, K)
      for(k in 1:K){
        if(k<4)
          rho[k] <- 1/(k+1)
        else
          rho[k] <- 1/(k+1)^2
      }
      
      theta <- matrix(0,K,P)
      for (k in 1:K) {
        if (k%%2 == 0)
          theta[k, ] <- sqrt(2) * sin(k*pi*t)
        else if (k%%2 != 0 && k != 1)
          theta[k, ] <- sqrt( 2 ) * cos(( k-1)*pi* t)
        else
          theta[k, ] <- rep(1, P)
      }
      
      s1 <- 0
      for (k in 1:4) {
        s1 <- s1 + sqrt(rho[k]) * theta[k, ]
      }
      
      s2 <- 0
      for (k in 4:K) {
        s2 <- s2 + sqrt(rho[k]) * theta[k, ]
      }
      
      m2_1 <- m1 + s1
      m2_2 <- m1 + s2
      
      zX <- matrix(0, n, K)
      zY <- matrix(0, n, K)
      
      uX <- matrix(0, n, P)
      uY <- matrix(0, n, P)
      for (i in 1:n) {
        zX[i, ] <- rnorm(K)
        zY[i, ] <- rnorm(K)
        for (k in 1:K) {
          uX[i, ] <- uX[i, ] + sqrt(rho[k]) * (zX[i, k] * theta[k, ])
          uY[i, ] <- uY[i, ] + sqrt(rho[k]) * (zY[i, k] * theta[k, ])
        }
      }
      
      X_ <- matrix(0, n, P) ## MODEL 10
      Y1 <- matrix(0, n, P)
      if(dataset == 11){  ## MODEL 11
        for (i in 1:n){
          X_[i, ] <- m1 + uX[i, ]
          Y1[i, ] <- m2_1 + uY[i, ]
        }
      } else if(dataset == 12){  ## MODEL 12
        for (i in 1:n){
          X_[i, ] <- m1 + uX[i, ]
          Y1[i, ] <- m2_2 + uY[i, ]
        }
      }
      main = paste("Model 10 & Model ", dataset, " (blue)")
      data <- rbind(X_, Y1)
    }
    if(dataset ==13){
      n <- 50
      rangeval = c(0,pi/3)
      t <- seq(0, pi/3, length = 100)
      truelabels <- c(rep(2,n), rep(1,n), rep(0,n)) # 1 is pink and 0 blue
      xlab <- "t"
      ylab <- "Value"
      colors <- c(color_1, color_2, color_3)
      pch <- c(16, 17, 18)
      #Scenario 1
      a1 <- runif(n,-0.25,0.25)
      b1 <- c(0.3,1,0.2)
      c1 <- c(1/1.3, 1/1.2, 1/4)
      X1 <- matrix(0,n,length(t))
      X2 <- matrix(0,n,length(t))
      X3 <- matrix(0,n,length(t))
      eps <- mvrnorm(n,mu=rep(2,length(t)), Sigma=(0.4^2)*diag(length(t)))
      for (i in 1:n){
        for (j in 1:length(t)){
          X1[i,j] <- a1[i]+b1[1]+c1[1]*sin(1.3*t[j])+(t[j])^3 
          X2[i,j] <- a1[i]+b1[2]+c1[2]*sin(1.3*t[j])+(t[j])^3
          X3[i,j] <- a1[i]+b1[3]+c1[3]*sin(1.3*t[j])+(t[j])^3
        }
      }
      
      Y13 <- X1+eps  ## MODEL 13
      Y14 <- X2+eps  ## MODEL 14
      Y15 <- X3+eps  ## MODEL 15
      data <-rbind(Y13, Y14, Y15)
      main=paste("Model 13 (pink), Model 14 (blue) & Model 15 (green) - Scenario 1")
    }
    if(dataset == 14){
      n <- 50
      rangeval = c(0,pi/3)
      t <- seq(0, pi/3, length = 100)
      truelabels <- c(rep(2,n), rep(1,n), rep(0,n)) # 1 is pink and 0 blue
      xlab <- "t"
      ylab <- "Value"
      colors <- c(color_1, color_2, color_3)
      pch <- c(16, 17, 18)
      #Scenario 2
      a1 <- runif(n,-0.5,0.5)
      b1 <- c(1.1, 1.5, 2.2)
      c1 <- c(1.5,1.7,1.9)
      X1 <- matrix(0,n,length(t))
      X2 <- matrix(0,n,length(t))
      X3 <- matrix(0,n,length(t))
      eps <- mvrnorm(n,mu=rep(2,length(t)), Sigma=(0.4^2)*diag(length(t)))
      for (i in 1:n){
        for (j in 1:length(t)){
          X1[i,j] <- a1[i]+b1[1]+sin(c1[1]*pi*t[j])+cos(pi*(t[j])^2)
          X2[i,j] <- a1[i]+b1[2]+sin(c1[2]*pi*t[j])+cos(pi*(t[j])^2)
          X3[i,j] <- a1[i]+b1[3]+sin(c1[3]*pi*t[j])+cos(pi*(t[j])^2)
        }
      }
      
      Y13 <- X1+eps  ## MODEL 16
      Y14 <- X2+eps  ## MODEL 17
      Y15 <- X3+eps  ## MODEL 18
      data <-rbind(Y13, Y14, Y15)
      main=paste("Model 16 (pink), Model 17 (blue) & Model 18 (green)")
    }
    if(dataset == 15){
      n <- 50
      rangeval = c(0,pi/3)
      t <- seq(0, pi/3, length = 100)
      truelabels <- c(rep(2,n), rep(1,n), rep(0,n)) # 1 is pink and 0 blue
      xlab <- "t"
      ylab <- "Value"
      colors <- c(color_1, color_2, color_3)
      pch <- c(16, 17, 18)
      #Scenario 3
      a1 <- runif(n,-0.25,0.25)
      b1 <- c(1/1.8,1/1.7,1/1.5)
      c1 <- c(1.1,1.4,1.5)
      X1 <- matrix(0,n,length(t))
      X2 <- matrix(0,n,length(t))
      X3 <- matrix(0,n,length(t))
      eps <- mvrnorm(n,mu=rep(2,length(t)), Sigma=(0.3^2)*diag(length(t)))
      for (i in 1:n){
        for (j in 1:length(t)){
          X1[i,j] <- a1[i]+b1[1]*exp(c1[1]*t[j])-(t[j])^3
          X2[i,j] <- a1[i]+b1[2]*exp(c1[2]*t[j])-(t[j])^3
          X3[i,j] <- a1[i]+b1[3]*exp(c1[3]*t[j])-(t[j])^3
        }
      }
      
      Y13 <- X1+eps  ## MODEL 19
      Y14 <- X2+eps  ## MODEL 20
      Y15 <- X3+eps  ## MODEL 21
      data <-rbind(Y13, Y14, Y15)
      main=paste("Model 19 (pink), Model 20 (blue) & Model 21 (green)")
    }
  }
  datalist <- list(data = data, truelabels = truelabels, rangeval =rangeval, t=t,
                   xlab = xlab, ylab = ylab, main = main, colors = colors, pch=pch )
  return(datalist)
}


##########################
### Indexes functions ####
##########################

# The following functions calculate the epigraph (EI), hypograph (HI), generalized
# epigraph (MEI) and generalized hypograph (MHI) indexes
EI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.min <- curves[i,]
    for (j in 1:B){
      inpoints <- sum((curves[j,] >= env.min))
      if (inpoints == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (1-index)
}

HI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.max <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,]<=env.max)
      if (inpoints == lengthcurves)
      {index[i] <- index[i]+1}
    }
  }
  index <- index/B
  return (index)
}

MEI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.min <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,] >= env.min)
      index[i] <- index[i]+inpoints
    }
  }
  index <- index/(B*lengthcurves)
  return (1-index)
}

MHI <- function(curves){
  index <- numeric()
  B <- dim(curves)[1]
  lengthcurves <- dim(curves)[2]
  for (i in 1:B){
    index[i] <- 0
    env.max <- curves[i,]
    for (j in 1:B){
      inpoints <- sum(curves[j,]<=env.max)
      index[i] <- index[i]+inpoints
    }
  }
  index <- index/(B*lengthcurves)
  return (index)
}

##########################
### funspline function ###
##########################

# This function fits a smooth B-spline basis for a bunch of curves X and
# calculates first and second derivatives for each curve.
funspline <- function(data, t, rangeval, nbasis, norder, basisobj){
  ys <- smooth.basis(argvals = t, y = t(data), fdParobj = basisobj) #data obtained when applying a smoothed B-spline basis
  smooth <- t(eval.fd(t,ys$fd,0)) #smoothed data
  deriv <- t(eval.fd(t,ys$fd,1)) #first derivatives
  deriv2 <- t(eval.fd(t,ys$fd,2)) #second derivatives
  
  res <- list("spline"=ys, "smooth"=smooth, "deriv"=deriv, "deriv2"=deriv2) #object containing the data and 1st and 2nd derivatives
  return(res)
}


########################
### indexes function ###
########################

# This function transform functional data into multivariate one using the indexes

indexes <- function(data, t, rangeval, nbasis, norder, truelabels){
  
  basisobj <- create.bspline.basis(rangeval = rangeval, norder = as.numeric(norder), nbasis = as.numeric(nbasis))
  dtaX <- funspline(data = data, t = t, rangeval = rangeval, nbasis = nbasis, norder = norder, basisobj = basisobj)
  
  dta <- data #original data
  # dta <- dtaX$smooth # smoothed data coefs
  ddta <- dtaX$deriv #first derivatives
  d2dta <- dtaX$deriv2 #second derivatives
  
  # Applying the indexes to the original data
  dtaEI <- EI(dta)
  dtaHI <- HI(dta)
  dtaMEI <- MEI(dta)
  dtaMHI <- MHI(dta)
  
  # Applying the indexes to the data first derivatives
  ddtaEI <- EI(ddta)
  ddtaHI <- HI(ddta)
  ddtaMEI <- MEI(ddta)
  ddtaMHI <- MHI(ddta)
  
  # Applying the indexes to the data second derivatives
  d2dtaEI <- EI(d2dta)
  d2dtaHI <- HI(d2dta)
  d2dtaMEI <- MEI(d2dta)
  d2dtaMHI <- MHI(d2dta)
  
  r <- as.factor(truelabels)
  
  # New multivariate data set
  indexes <- data.frame(dtaEI, dtaHI, dtaMEI, dtaMHI,
                        ddtaEI, ddtaHI, ddtaMEI, ddtaMHI,
                        d2dtaEI, d2dtaHI, d2dtaMEI, d2dtaMHI, r)
  res <- list("indexes" = indexes, "spline"=dtaX$ys, "smooth"=dtaX$smooth, "deriv"=dtaX$deriv, "deriv2"=dtaX$deriv2)
  return(res)
}

###################################
### multivariatemodels function ###
###################################

# This function apply multivariate classification techniques to the multivariate
# dataset obtained from the indexes

multivariatemodels <- function(indexes){
  sample <- sample(c(TRUE, FALSE), nrow(indexes), replace=TRUE, prob=c(0.7,0.3))
  train  <- indexes[sample, ]
  test   <- indexes[!sample, ]
  # original indexes
  VARS1 <- c("dtaEI","dtaHI")
  VARS2 <- c("ddtaEI","ddtaHI")
  VARS3 <- c("d2dtaEI","d2dtaHI")
  VARS4 <- c(VARS1,VARS2)
  VARS5 <- c(VARS1,VARS3)
  VARS6 <- c(VARS2,VARS3)
  VARS7 <- c(VARS4,VARS3)
  
  # generalized indexes
  VARS8 <- c("dtaMEI","ddtaMEI")
  VARS9 <- c("dtaMEI","d2dtaMEI")
  VARS10 <- c("ddtaMEI","d2dtaMEI")
  VARS11 <- c(VARS8,"d2dtaMEI")
  
  # combining both indexes types
  VARS12 <- c(VARS1,"dtaMEI")
  VARS13 <- c(VARS2,"ddtaMEI")
  VARS14 <- c(VARS3,"d2dtaMEI")
  VARS15 <- c(VARS4,VARS8)
  VARS16 <- c(VARS5,VARS9)
  VARS17 <- c(VARS6,VARS10)
  VARS18 <- c(VARS7,VARS11)
  
  VARS <- c("VARS1","VARS2","VARS3","VARS4","VARS5","VARS6","VARS7","VARS8","VARS9",
            "VARS10","VARS11","VARS12","VARS13","VARS14","VARS15","VARS16","VARS17","VARS18")
  technique <- c("glm","naive_bayes","svmLinear","knn","ranger","nnet")
  #technique <- c("glm","naive_bayes")
  times <- c()
  accuracy <-c()
  predictors <- c()
  model_name <- c()
  techname <- c("GLM","NB","SVM","KNN","RF","NNET")
  results <- data.frame()
  for(i in 1:length(technique)){
    for(j in 1:length(VARS)){
      predictors[j] <- paste(eval(parse(text = VARS[j])), collapse = "+")
      tryCatch({
        t0 <- Sys.time()
        model <- train(as.formula(paste0("r","~",predictors[j])),
                       method = technique[i],
                       data = train,
                       metric ="Accuracy")
        predictions <- predict(model,test)
        CM <- confusionMatrix(predictions,test$r)
        t1 <- Sys.time()
        accuracy <- CM$overall[1]
        results<-rbind(results,data.frame(Model = paste(techname[i], VARS[j], sep="_"), Technique = techname[i],
                                          Rfunction = technique[i], Variables = VARS[j], Accuracy = accuracy, Time = t1-t0))  
      },error = function(e){})
    }
  }
  row.names (results) <- 1: nrow (results)
  results<-results[order(results$Time, decreasing = F),]
  results<-results[order(results$Accuracy, decreasing = T),]
  
  bestmodels <- results %>%
    group_by(Technique) %>%
    filter(Accuracy == max(Accuracy, na.rm=TRUE))
  colnames(bestmodels) <- c("Model", "Technique", "Rfunction", "Variables", "Accuracy", "Time")
  bestmodels<-bestmodels[order(bestmodels$Accuracy, decreasing = T),]
  
  bestmodels <- bestmodels %>%
    group_by(Technique) %>%
    filter(Time == min(Time, na.rm=TRUE))
  
  res <- list("results" = results, "bestmodels" = bestmodels)
  return(results)
}


##############################
### DD_classifier function ###
##############################

# This function contains functional classifiers based on depths:
# maximum depth classifier and DD^G classifier.

DD_classifier <- function(r_train, X_train_func, test_fd, X_test_func, r_test_df){
  depth<-c("FM","mode","RP")
  classifier <- c("DD1","DD2","DD3","glm","gam","lda","qda","knn","knn")
  for (i in 1:length(depth)) {
    t0 <- Sys.time()
    fda_class_depth <- classif.depth(r_train,X_train_func,X_test_func)
    Accuracy_depth <- sum(r_test==fda_class_depth$group.pred)/nrow(X_test_func)
    t1 <- Sys.time()
    results<-rbind(results, data.frame(Technique = paste("MD", depth[i], sep="_"), Accuracy = Accuracy_depth, Time = t1-t0))
    for (j in 1:length(classifier)) {
      t0 <- Sys.time()
      print(i)
      fda_class_DD <- classif.DD(r_train,X_train_func,control=list(draw=FALSE), 
                                 depth = depth[i], classif = classifier[j])
      print(i)
      pred_test_DD <- predict(fda_class_DD,test_fd,type="class")
      print(i)
      Accuracy_DD <- confusionMatrix(pred_test_DD, as.factor(r_test_df))$overall[1]
      t1 <- Sys.time()
      results<-rbind(results, data.frame(Technique = paste(classifier[j], depth[i], sep="_"), Accuracy = Accuracy_DD, Time = t1-t0))
    }
  }
  return(results)
}

################################
### classicalmodels function ###
################################

# Thus function apply functional classifiers to the original data.

classicalmodels <- function(data, t, rangeval, r, nbasis, norder){
  results <- data.frame()
  sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7,0.3))
  
  X_train <- data[sample,]
  r_train <- r[sample]
  X_test <- data[-sample,]
  r_test <- r[-sample]
  
  X_train_func <- fdata(X_train,argvals=t,rangeval=rangeval)
  X_test_func <- fdata(X_test,argvals=t,rangeval=rangeval)
  r_train_df <- data.frame(r_train)
  r_test_df <- r_test
  
  
  train_pairs <- list("df"=r_train_df,"x"=X_train_func)
  test_fd <- list("x"=X_test_func)
  
  
  basis <- create.bspline.basis(rangeval=rangeval, norder=as.numeric(norder), nbasis = as.numeric(nbasis))
  
  results_DD <- DD_classifier(r_train = r_train, X_train_func = X_train_func, test_fd = test_fd, 
                              X_test_func = X_test_func, r_test_df = r_test_df)
  results <- rbind(results, results_DD)
  
  return(results)
}


################################
###            UI            ###
################################

ui <- fluidPage(
  theme = shinytheme("yeti"),
  tags$head(tags$style(
    HTML('
    @import url("https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap");
      body {
        background-color: #FFFAFA;
      }
      dataTableOutput{
        background-color: white;
      }
         ')
  )),
  navbarPage("Functional Data Analysis",
             ### TAB 1
             tabPanel("Data Selection",
                      textOutput("SelectedData"),
                      sidebarLayout(
                        sidebarPanel(id="DataSelection",
                                     uiOutput("DataSelection"),
                                     selectInput(inputId = "Data",
                                                 label = "Select a dataset:",
                                                 choices = c("Growth"="growth",
                                                             "Model 1 & Model 2" = 2,
                                                             "Model 1 & Model 3" = 3,
                                                             "Model 1 & Model 4" = 4,
                                                             "Model 1 & Model 5" = 5,
                                                             "Model 1 & Model 6" = 6,
                                                             "Model 1 & Model 7" = 7,
                                                             "Model 1 & Model 8" = 8,
                                                             "Model 1 & Model 9" = 9,
                                                             "Model 10 & Model 11" = 11,
                                                             "Model 10 & Model 12" = 12,
                                                             "Model 13, Model 14 & Model 15 - Scenario 1" = 13,
                                                             "Model 13, Model 14 & Model 15 - Scenario 2" = 14,
                                                             "Model 13, Model 14 & Model 15 - Scenario 3" = 15),
                                                 selected = NULL),
                                     actionButton("Select","Select")
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Dataset", DT::dataTableOutput("DataTable")),
                            tabPanel("Information", verbatimTextOutput("Information")),
                            tabPanel("Plot", plotOutput("DataPlot"))
                          )
                        )
                      )
             ),### END TAB 1
             ### TAB 2
             navbarMenu("Data smoothing",
                        ### TAB PANEL 1
                        tabPanel("Data smoothing",
                                 sidebarLayout(
                                   sidebarPanel(id="SmoothDataSelection",
                                                uiOutput("nbasis"),
                                                uiOutput("norder"),
                                                actionButton("Apply","Apply")
                                   ),
                                   mainPanel(
                                     #textOutput("SmoothParameters"),
                                     tabsetPanel(
                                       tabPanel("Smoothed Data", DT::dataTableOutput("DataTableSmooth")),
                                       tabPanel("First Derivatives", DT::dataTableOutput("DataTableDeriv")),
                                       tabPanel("Second Derivatives", DT::dataTableOutput("DataTableDeriv2"))
                                     )
                                   )
                                 )
                        ),
                        # TAB PANEL 2
                        tabPanel("Plots",
                                 plotOutput("SmoothDataPlot")
                        )
             ), ### END TAB 2
             ### TAB 3
             navbarMenu("Indexes Application",
                        # TAB PANEL 1
                        tabPanel("Indexes Data",
                                 DT::dataTableOutput("Indexes")
                        ),
                        # TAB PANEL 2
                        tabPanel("Plots",
                                 plotOutput("IndexesPlots")
                        )
             ), ### END TAB 3
             ### TAB 4
             navbarMenu("Multivariate Analysis",
                        tabPanel("Multivariate results",
                                 DT::dataTableOutput("Results")
                        ),
                        tabPanel("Plots",
                                 fluidRow(
                                   box(title = "Results mean by model", plotOutput("Meanplot"))
                                 ),
                                 
                                 fluidRow(
                                   box(title = "", DT::dataTableOutput("Meantable"))
                                 )
                                 
                        ),
                        tabPanel("Print",
                                 h5("Results:"),
                                 verbatimTextOutput("Print")
                        )
             ),
             ### TAB 5
             navbarMenu("Functional Data Analysis",
                        tabPanel("FDA results",
                                 DT::dataTableOutput("FDAResults")
                        ),
                        tabPanel("Plots",
                                 # plotOutput("FDAResultsPlots"),
                                 #DT::dataTableOutput("FDAmeans")
                                 fluidRow(
                                   box(title = "Functional Data Analysis", plotOutput("FDAResultsPlots"))
                                 ),
                                 fluidRow(
                                   box(title = "", DT::dataTableOutput("FDAmeans"))
                                 )
                        ),
                        tabPanel("Print",
                                 verbatimTextOutput("Print2")
                        )
             )
  )
)#UI

################################
###        SERVER            ###
################################

server <- function(input, output, session) {
  #############
  ### TAB 1 ###
  #############
  v = reactiveValues(path = NULL)
  observeEvent(input$Select, {
    req(input$Data)
    datalist <- data(input$Data)
    v$data <- datalist$data
    v$truelabels <- datalist$truelabels
    v$rangeval <- datalist$rangeval
    v$t <- datalist$t
    v$xlab <- datalist$xlab
    v$ylab <- datalist$ylab
    v$main <- datalist$main
    v$colors <- datalist$colors
    v$pch <- datalist$pch
  })
  
  output$SelectedData <- renderText({
    paste("You chose", input$Data)
  })
  # TABPANEL 1
  output$DataTable <- DT::renderDataTable({
    DT::datatable(data.frame(v$data))      
  })
  # TABPANEL 2
  output$Information <- renderPrint({
    if (is.null(v$data))
      return(NULL)
    str(data.frame(t(v$data)))
  })
  # TABPANEL 3
  output$DataPlot <- renderPlot({
    if (is.null(v$data))
      return(NULL)
    matplot(v$t,t(v$data),type="l",lty=1,
            col=v$colors[factor(v$truelabels)],
            xlab= v$xlab, ylab= v$ylab,
            main= v$main )
  }, height = 750)
  
  output$SelectedData <- renderText({
    paste("You chose", input$Data)
  })
  output$nbasis <- renderUI({
    req(v$data)
    selectInput(inputId = "nbasis",
                label = "Number of basis: ",
                choices = seq(min(v$t)+1, length(v$t)),
                selected = 30)
  })
  output$norder <- renderUI({
    req(v$data)
    selectInput(inputId = "norder",
                label = "Order of b-splines: ",
                choices = seq(2, min(input$nbasis,20)),
                selected = 4)
  })
  #############
  ### TAB 2 ###
  #############
  s = reactiveValues(path = NULL)
  observeEvent(input$Apply, {
    datasets <- indexes (data = v$data, t = v$t, rangeval = v$rangeval, nbasis = input$nbasis, norder = input$norder, truelabels = v$truelabels)
    s$ys <- datasets$ys
    s$smooth <- datasets$smooth
    s$deriv <- datasets$deriv
    s$deriv2 <- datasets$deriv2
    s$ind.data <- datasets$indexes
  })
  output$DataTableSmooth <- DT::renderDataTable({
    if (is.null(s$smooth))
      return(NULL)
    DT::datatable(data.frame(round(s$smooth,4)))
  })
  output$DataTableDeriv <- DT::renderDataTable({
    if (is.null(s$deriv))
      return(NULL)
    DT::datatable(data.frame(round(s$deriv,4)))
  })
  output$DataTableDeriv2 <- DT::renderDataTable({
    if (is.null(s$deriv2))
      return(NULL)
    DT::datatable(data.frame(round(s$deriv2,4)))
  })
  # TAB PANEL 2
  output$SmoothDataPlot <- renderPlot({
    if (is.null(s$deriv2))
      return(NULL)
    par(mfrow =c(2,2) )
    matplot(v$t,t(v$data),type="l",lty=1,
            col=v$colors[factor(v$truelabels)],
            xlab= v$xlab, ylab= v$ylab,
            main= v$main )
    matplot(v$t,t(s$smooth),type="l",lty=1,
            col=v$colors[factor(v$truelabels)],
            xlab=v$xlab, ylab= v$ylab,
            main="Smoothed data")
    matplot(v$t,t(s$deriv),type="l",lty=1,
            col=v$colors[factor(v$truelabels)],
            xlab=v$xlab, ylab= v$ylab,
            main="First derivatives")
    matplot(v$t,t(s$deriv2),type="l",lty=1,
            col=v$colors[factor(v$truelabels)],
            xlab=v$xlab, ylab= v$ylab,
            main="Second derivatives")
  }, height = 1000)
  
  #############
  ### TAB 3 ###
  #############
  
  # TAB PANEL 1
  output$Indexes <- DT::renderDataTable({
    if (is.null(s$ind.data))
      return(NULL)
    DT::datatable(s$ind.data)
  })
  # TAB PANEL 2
  output$IndexesPlots <- renderPlot({
    if (is.null(s$ind.data))
      return(NULL)
    par(mfrow =c(3,3) )
    par(cex=1.5, cex.main=1.5, cex.lab = 1.25, cex.sub=1.25)
    matplot(v$t,t(v$data),type="l",lty=1,
            col= v$colors[factor(v$truelabels)],
            xlab=v$xlab,ylab=v$ylab,
            main=v$main, cex=1.5)
    plot(s$ind.data[,1],s$ind.data[,2],
         pch = v$pch[as.factor(v$truelabels)],
         col = v$colors[factor(v$truelabels)],
         xlim = c(0,1), ylim = c(0,1), xlab = "EI", ylab = "HI",cex=1.5)
    plot(s$ind.data[,3],s$ind.data[,4],
         pch = v$pch[as.factor(v$truelabels)],
         col = v$colors[factor(v$truelabels)],
         xlim = c(0,1), ylim = c(0,1), xlab = "MEI", ylab = "MHI",cex=1.5)
    
    matplot(v$t,t(s$deriv),type="l",lty=1,
            col= v$colors[factor(v$truelabels)],
            xlab=v$xlab,ylab=v$ylab,
            main="First derivatives",cex=1.5)
    plot(s$ind.data[,5],s$ind.data[,6],
         pch = v$pch[as.factor(v$truelabels)],
         col = v$colors[factor(v$truelabels)],
         xlim = c(0,1), ylim = c(0,1), xlab = "EI", ylab = "HI",cex=1.5)
    plot(s$ind.data[,7],s$ind.data[,8],
         pch = v$pch[as.factor(v$truelabels)],
         col = v$colors[factor(v$truelabels)],
         xlim = c(0,1), ylim = c(0,1), xlab = "MEI", ylab = "MHI",cex=1.5)
    
    matplot(v$t,t(s$deriv2),type="l",lty=1,
            col= v$colors[factor(v$truelabels)],
            xlab=v$xlab,ylab=v$ylab,
            main="Second derivatives",cex=1.5)
    plot(s$ind.data[,9],s$ind.data[,10],
         pch = v$pch[as.factor(v$truelabels)],
         col = v$colors[factor(v$truelabels)],
         xlim = c(0,1), ylim = c(0,1), xlab = "EI", ylab = "HI",cex=1.5)
    plot(s$ind.data[,11],s$ind.data[,12],
         pch = v$pch[as.factor(v$truelabels)],
         col = v$colors[factor(v$truelabels)],
         xlim = c(0,1), ylim = c(0,1), xlab = "MEI", ylab = "MHI",cex=1.5)
  }, height = 1300)
  
  #############
  ### TAB 4 ###
  #############
  # TAB PANEL 1
  
  output$Results <- DT::renderDataTable({
    mult_models <- data.frame()
    clas_models <- data.frame()
    for (b in 1:50) {
      print("THIS IS MODEL")
      print(b)
      set.seed(b)
      dt <- data (input$Data)
      # Our approach
      ind <- indexes (data = dt$data, t = dt$t, rangeval = dt$rangeval, nbasis = input$nbasis, norder = input$norder, truelabels = dt$truelabels)
      mult <- multivariatemodels(ind$indexes)
      mult_models <- rbind(mult_models, mult)
      # # Other models
      clas <- classicalmodels (data = dt$data, t = dt$t, rangeval = dt$rangeval, nbasis = input$nbasis, norder = input$norder, r = factor(dt$truelabels))
      clas_models <- rbind(clas_models, clas)
    }
    s$mult_models <- mult_models
    mean_mult <- mult_models %>%
      group_by(Model, Technique, Variables) %>%
      summarise(SD = sd(Accuracy), Accuracy = mean(Accuracy), Time = mean(Time))
    mean_mult<-mean_mult[order(mean_mult$Accuracy, decreasing = T),]
    s$mean_mult <- mean_mult
    
    s$clas_models <- clas_models
    mean_clas <- clas_models %>%
      group_by(Technique) %>%
      summarise(SD = sd(Accuracy), Accuracy = mean(Accuracy), Time = mean(Time))
    mean_clas<-mean_clas[order(mean_clas$Accuracy, decreasing = T),]
    s$mean_clas <- mean_clas
    
    DT::datatable(s$clas_models)
  })
  # TAB PANEL 2
  output$Meanplot <- renderPlot({
    if (is.null(s$mult_models))
      return(NULL)
    par(mfrow =c(1,2) )
    s$mult_models %>%
      ggplot( aes(x=Model, y=Accuracy, fill=Model)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      theme_ipsum() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      ggtitle("Our approach") +
      xlab("Model")+
      theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.3, hjust=0.8))
  }, height = 450, width = 1000)
  
  output$Meantable <- DT::renderDataTable({
    options(digits=4)
    DT::datatable(s$mean_mult)
  })
  
  # TAB PANEL 3
  output$Print <- renderPrint({
    if (is.null(s$mean_mult))
      return(NULL)
    options(digits=4)
    
    GLM <- NB <- SVM <- KNN <- RF <- NNET <- c()
    
    for (i in 1:nrow(s$mean_mult)) {
      if (s$mean_mult$Technique[i] == "GLM")  {
        GLM[i] <- c(i)
      }
      if (s$mean_mult$Technique[i] == "NB")  {
        NB[i] <- c(i)
      }
      if (s$mean_mult$Technique[i] == "SVM")  {
        SVM[i] <- c(i)
      }
      if (s$mean_mult$Technique[i] == "KNN")  {
        KNN[i] <- c(i)
      }
      if (s$mean_mult$Technique[i] == "RF")  {
        RF[i] <- c(i)
      }
      if (s$mean_mult$Technique[i] == "NNET")  {
        NNET[i] <- c(i)
      }
    }
    GLM <- na.omit(GLM)
    NB<- na.omit(NB)
    SVM <- na.omit(SVM)
    KNN <- na.omit(KNN)
    RF<- na.omit(RF)
    NNET <- na.omit(NNET)
    
    s$mean_mult %>%
      kable(format = 'latex', booktabs = TRUE, row.names = F)%>%
      row_spec(GLM, bold = F, color = "black", background = "#F7CFD8")%>%
      row_spec(NB, bold = F, color = "black", background = "#E2CEF7")%>%
      row_spec(SVM, bold = F, color = "black", background = "#CFE3F7")%>%
      row_spec(KNN, bold = F, color = "black", background = "#CEF6E3")%>%
      row_spec(RF, bold = F, color = "black", background = "#F4F7CE")%>%
      row_spec(NNET, bold = F, color = "black", background = "#F6D8CE")
  })
  
  #############
  ### TAB 5 ###
  #############
  # TAB PANEL 1
  output$FDAResults <- DT::renderDataTable({
    if (is.null(s$clas_models))
      return(NULL)
    DT::datatable(s$clas_models)
  })
  # TAB PANEL 2
  
  output$FDAmeans <- DT::renderDataTable({
    DT::datatable(s$mean_clas)
  })
  output$FDAResultsPlots <- renderPlot({
    if (is.null(s$mean_clas))
      return(NULL)
    s$clas_models %>%
      ggplot( aes(x=Technique, y=Accuracy, fill=Technique)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      theme_ipsum() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      ggtitle("Other approaches") +
      xlab("Technique")+
      theme(axis.text.x = element_text(angle = 90, size = 8,vjust = 0.3, hjust=0.8))
  }, height = 450, width = 1000)
  # TAB PANEL 3
  output$Print2 <- renderPrint({
    if (is.null(s$mean_clas))
      return(NULL)
    options(digits=4)
    s$mean_clas %>%
      kable(format = 'latex', booktabs = TRUE, row.names = F)
  })
}



# Run the application
shinyApp(ui = ui, server = server)