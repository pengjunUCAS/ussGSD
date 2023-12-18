
  ###
  ###**************************************************************************************************************************
  ### R program ussGSD (v2.3).
  ###**************************************************************************************************************************
  ### Methods and open source R code for efficient and flexible unmixing of single-sample grain-size distributions.
  ### Author: Jun Peng, Hunan University of Science and Technology, China. 
  ### Last updated, 2023.12.18.
  ### 
  ### Note that the current version of the program is developed for unmixing of
  ### single-sample grain-size distributions measured using the Malvern Mastersizer-2000/3000
  ### laser particle-size analyzer. So it can not ensure that the program will work equally 
  ### well for datasets measured from other types of particle-size analyzers.
  ### 
  ### Please contact Jun Peng (pengjun10@mails.ucas.ac.cn, or 656971673@qq.com) 
  ### if you have any problem during the application of the program or 
  ### if you want to report any bug encountered during the use of the program.
  ###
  ### References:
  ###  Peng, J., Zhao, H., Dong, Z.B., Zhang, Z.C., Yang, H.Y., Wang, X.L., 2022. Numerical methodologies and tools for efficient 
  ###  and flexible unmixing of single-sample grain-size distributions: Application to late Quaternary aeolian sediments from the 
  ###  desert-loess transition zone of the Tengger Desert. Sedimentary Geology, 438, 106211.
  ###
  ###  Peng, J., Wang, X.L., Zhao, H., Dong, Z.B, 2023. Single-sample unmixing and parametric end-member modelling of grain-size 
  ###  distributions with transformed probability density functions and their performance comparison using aeolian sediments. 
  ###  Sedimentary Geology, 445, 106328. 
  ###**************************************************************************************************************************
  
  ###
  ###****************************************************************************************
  ### Install external packages from CRAN.
  ###****************************************************************************************
  PACKAGE <- c("pracma", "minpack.lm", "parallel", "foreach", "snow", "doSNOW")
  idx <- which(!PACKAGE %in% installed.packages()[,1])
  if (length(idx)>0)  install.packages(PACKAGE[idx])
  ###****************************************************************************************   

  ###
  ###*******************************************************
  ### Load the required external R packages.
  ### Please make sure that the following packages have  
  ### been installed in your R calculation environment.
  ###=======================================================
  library(pracma)          # install.packages("pracma")
  ###
  library(minpack.lm)      # install.packages("minpack.lm")
  ###=======================================================

  
  ###
  ###***********************************************************************************
  ### Users may also need to install and load the following external R package 
  ### if they want to perform a parallel unmixing of grain-size distributions in a batch 
  ### pattern using the following R functions update_ussGSDbatchp() or ussGSDbatchp().
  ###===================================================================================
  library(parallel)      # install.packages("parallel")
  ###
  library(foreach)       # install.packages("foreach")
  ###
  library(snow)          # install.packages("snow")
  ###
  library(doSNOW)        # install.packages("doSNOW")
  ###=========================================================================================


  ###
  ###**************************************************************************************************************************
  ### Function summary_ussGSDbatch() is used for generated a pdf file containing the 
  ### a summary of unmixing results of grain-size distributions obtained from a batch model. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ### obj_batchgsd: an S3 object of class "gsdbatch" generated using the function
  ###               ussGSDbatch(), ussGSDbatchp(), update_ussGSDbatch(), or update_ussGSDbatchp().
  ###
  ###           sf: a real number indicating the factor used for scaling the standard deviation (SD) of
  ###               a unmixed grain-size component. sf<1 indicates that the SD will be shrinked to sf*SD
  ###               sf=1 means that the SD will remain unchanged during the analysis. 
  ###
  ###     preclude: A vector consisting of integers representing the indexes of samples 
  ###               which will not participate the analysis. 
  ###
  ###          nkm: An integer indicating the number of clusters to be generated (from 2 to 13)
  ###               during the application of the k-means clustering algorithm in analyzing 
  ###               the unmixed grain-size mode and median.
  ###
  ###          nfm: An integer indicating the number of clusters to be generated (from 2 to 13)
  ###               during the application of the Bayesian clustering algorithm based on the 
  ###               variance-accounted finite mixture model in analysing the unmixed mean grain 
  ###               sizes of individual components.
  ###
  ###     addsigma: A real number indicating the additional (relative) error that will be added 
  ###               in quadratic to the relative standard error of mean grain sizes of individual 
  ###               components before the implementation of the Bayesian clustering algorithm.
  ###
  ###         logx: A logical value indicating whether the grain-size levels will be ploted in a logged scale.
  ###
  ###          mps: A real number indicating the upper limit for the precision of data points visualized in the radial plot.  
  ###
  ###      outfile: A character indicating the prefix name of the PDF/CSV file the modelling results will  
  ###               be written to. The results will be returned to files with a default prefix name if outfile=NULL. 
  ###==========================================================================================================================
  ### The function 
  ### Generates a PDF file with postfix of "_summary" in the current working directory, and returns  
  ### an invisible list containing the following elements (these elements will be saved to CSV files  
  ### with postfix of "_gsfit" and "_pgsp"):
  ###
  ### gsfit: A data.frame summarizing the fitting results of individual GSDs, with column names of        
  ###        (1) NO-sample number, (2) Model-fitting model, (3) ncomp-number of components, 
  ###        (4) FOM-figure of merit, (5) RSS-residual sum of squares, (6) R2-R2 statistic, 
  ###        and (7) RSE-residual standard error.
  ###
  ###  pgsp: A matrix showing the pooled grain-size parameters and associated clustering results, 
  ###        with column names of (1) NO-sample number, (2) Median-component median, 
  ###        (3) MedianKM-KM cluster the median belongs to, (4) Mode-component mode,
  ###        (5) ModeKM-KM cluster the mode belongs to, (6) Mean-component mean,
  ###        (7) Sd-component standard deviation, (8) MeanFMM-FMM cluster the mean belongs to, 
  ###        and (9) Abundance-component abundance. 
  ###
  ### Note that if obj_batchgsd=NULL, the user needs to ensure that the function load_ussGSDbatch() has been called
  ### to import a RData file from the current working directory which contains an object of S3 class "batchgsd".
  ###==========================================================================================================================
  summary_ussGSDbatch <- function(obj_batchgsd=NULL, preclude=NULL, sf=1, nkm=NULL, nfm=NULL, addsigma=0, logx=TRUE, mps=NULL, outfile=NULL) {

      stopifnot(is.null(preclude) || is.numeric(preclude),
                is.numeric(sf), length(sf)==1L, sf>0, sf<2, 
                is.null(nkm) || (is.numeric(nkm) && length(nkm)==1 && nkm %in% (2:13)),
                is.null(nfm) || (is.numeric(nfm) && length(nfm)==1 && nfm %in% (2:13)),
                length(addsigma)==1, is.numeric(addsigma), addsigma>=0,
                length(logx)==1, is.logical(logx), is.null(mps) || (is.numeric(mps) && length(mps)==1),
                is.null(outfile) || (length(outfile)==1 && is.character(outfile)))

      ###
      if (!is.null(obj_batchgsd)) {

          if(inherits(obj_batchgsd,what="batchgsd")==FALSE)  stop("Error: [obj_batchgsd] should be an S3 object of class 'batchgsd'!")
          GSDbatch <- obj_batchgsd

      } else {

          if (!exists("gsdbatch"))  stop("Error: function [load_ussGSDbatch] has not been called!")
          GSDbatch <- get("GSDbatch", envir=gsdbatch)

      } # end if.

      ###
      xoxoxo <- attr(GSDbatch, "xoxoxo")

      ###
      N <- length(GSDbatch)
      
      ###
      if (!is.null(preclude)) { 
          
          if (!all(preclude %in% seq(N))) stop("Error: invalid argument [preclude]!")

          ###
          for (i in 1:length(preclude)) GSDbatch[[preclude[i]]]$FOM <- Inf

      } # end if.


      ###
      ookk <- sapply(GSDbatch, function(x) (!is.null(x)) && (is.finite(x$FOM)))
      if (length(ookk)==0)  stop("Error: no data can be used to generate the summary file!")

      ###
      GSDbatchx <- GSDbatch[ookk]
      Nx <- length(GSDbatchx)
      
      ### 
      FOM <- sapply(GSDbatchx, function(x) x$FOM)
      RSS <- sapply(GSDbatchx, function(x) x$RSS)
      R2 <- sapply(GSDbatchx, function(x) x$R2)
      RSE <- sapply(GSDbatchx, function(x) x$RSE)

      ###
      compvec <- sapply(GSDbatchx, function(x) nrow(x$gs.pars))

      ###
      abundance <- c(sapply(GSDbatchx, function(x) x$gs.pars[,1]), recursive=TRUE)

      ###
      meanGZ <- c(sapply(GSDbatchx, function(x) x$gs.pars[,2]), recursive=TRUE)

      ###
      medianGZ <- c(sapply(GSDbatchx, function(x) x$gs.pars[,3]), recursive=TRUE)

      ###
      modeGZ <- c(sapply(GSDbatchx, function(x) x$gs.pars[,4]), recursive=TRUE)

      ###
      sdGZ <- c(sapply(GSDbatchx, function(x) x$gs.pars[,5]), recursive=TRUE)*sf

      ###
      threshold <- max(mean(sdGZ/meanGZ)-2*sd(sdGZ/meanGZ),0.05)
      for (i in 1:Nx) {

          rsdvmin <- GSDbatchx[[i]]$gs.pars[,5,drop=TRUE]/GSDbatchx[[i]]$gs.pars[,2,drop=TRUE]
          
          ###
          if (any(rsdvmin<threshold)) cat("NOTE: unmixing results of the ", 
             (seq(N))[ookk][i], "-th GSD contains very small RSD values which will be rescaled!\n",sep="")

      } # end if.

      ###
      oldpar <- par("mar", "mfrow", "bg", "new")
      on.exit(par(oldpar))

      ###
      matKM <- matrix(nrow=8, ncol=2)

      ###
      for (i in 2:9) {

          ithKM <- stats::kmeans(log(modeGZ), centers=i, nstart=200)

          ###
          matKM[i-1,1] <- ithKM$tot.withinss 
          matKM[i-1,2] <- ithKM$betweenss

      } # end for.

      ### 
      par("mar"=c(5.1,4.1,4.1,4.1))

      ###
      plot(2:9, matKM[,1], type="o", yaxt="n", col="skyblue", cex=3, lwd=2, xlab="Number of KM clusters [nkm]", 
           ylab="Total within-cluster sum of squares (TWSS)", cex.lab=1.2, main="Variation of K-Means statistics with [nkm]", cex.main=1.5)

      ###
      axis(side=2, col="skyblue", col.ticks="skyblue", lwd=3)

      ###
      par("new"=TRUE)
          
      ###
      plot(2:9, matKM[,2], type="o", xaxt="n", yaxt="n", col="red", xlab="", ylab="", cex=3, lwd=2)

      ###
      axis(side=4, col="red",col.ticks="red",lwd=3)

      ###
      mtext("Total between-cluster sum of squares (TBSS)", side=4, line=2, cex=1.2)

      ###
      if (is.null(nkm))  {

          nkm <- as.numeric(readline(paste("Please enter the number of KM clusters (i.e., [nkm]): ")))
          if (!is.finite(nkm))  stop("[nkm] has not been provided!")

          ###
          if (nkm>9)  stop("[nkm] should not exceed 9!")

      } # end if.

      ###
      abline(v=nkm, lwd=3, lty=2)

      ###
      if (is.null(nfm))  nfm <- 0

      ###
      sdGZ0 <- sdGZ
      toosmallidx <- sdGZ/meanGZ<threshold
      if (length(toosmallidx)>0L)  sdGZ0[toosmallidx] <- meanGZ[toosmallidx]*threshold

      ###
      ### Optimize the finite mixture age model, function 1.
      ###-----------------------------------------------------------------
      fmmED <- function(ed, sed, ncomp, addsigma=0, iflog=TRUE) {

          ###
          if (iflog==TRUE) {

              sed <- sqrt((sed/ed)^2+addsigma^2)
              ed <- log(ed)

          } else {

              sed <- sqrt(sed^2+addsigma^2)

          } # end if.

          ###
          w <- 1.0/sed^2

          ### 
          ndat <- length(ed)
            
          ###
          pf <- matrix(nrow=ndat, ncol=ncomp)
          muvec <- min(ed) +(max(ed)-min(ed))/(ncomp+4)*((1:(ncomp+5))-1)
          maxval <- -1.0e20

          ###
          info <- 1
          for (kk in 1:6) {

              inip <- rep(1.0/ncomp, ncomp)
              inimu <- muvec[kk:(kk+ncomp-1)]

              ###
              repeat {

                  for (j in 1:ncomp)  pf[,j] <- inip[j]*sqrt(w)*exp(-0.5*w*(ed-inimu[j])^2)

                  ###
                  if (any(!is.finite(pf)))  stop("Optimization of FMM failed!")

                  ###
                  pp <- pf/rowSums(pf)
                  wp <- w*pp
                  pv <- pp*w*ed

                  ###
                  p1 <- colSums(pp)/ndat
                  mu1 <- colSums(pv)/colSums(wp)

                  ###
                  if (sum(abs(inip-p1))+sum(abs(inimu-mu1))<=1.0e-8)  break

                  ###
                  inip <- p1
                  inimu <- mu1

              } # end repeat.

              ###
              maxlik <- sum(log(1.0/sqrt(2.0*pi)*rowSums(pf)))
              bic <- -2.0*maxlik + (2.0*ncomp-1.0)*log(ndat)
            
              ###
              if (maxlik>maxval) {

                  cp <- inip
                  cmu <- inimu
                  cmaxlik <- maxlik
                  cbic <- bic

                  ###              
                  maxval <- maxlik
                  info <- 0
                
              } # end if.

          } # end for.

          ###
          if (info==1)  stop("Optimization of FMM failed!")
            
          ###
          p <- cp
          mu <- cmu
          maxlik <- cmaxlik
          bic <- cbic
            
          ###
          if (iflog==TRUE)  mu <- exp(mu)

          ###
          idx <- order(mu)
          pars <- cbind(p[idx],mu[idx])
          colnames(pars) <- c("Prop","Mu") 
          rownames(pars)<-paste(rep("Comp", ncomp), 1:ncomp, sep="") 

          ###
          output <- list("pars"=pars, "maxlik"=maxlik, "bic"=bic)
          return(output)

      } # end function fmmED. 

      ###
      ### Optimize the finite mixture age model, function 2.
      ###-----------------------------------------------------------------
      optFMM <- function(ed, sed, ncomp, maxcomp=9, addsigma=0, iflog=TRUE) {

          if (ncomp==0) {

              minval <- 1.0e20
              LIST <- vector(length=maxcomp-1, mode="list")
 
              ###
              for (i in 2:maxcomp) {

                  SAM <- try(fmmED(ed=ed, sed=sed, ncomp=i, addsigma=addsigma, iflog=iflog), silent=TRUE)
                    
                  ###
                  if (inherits(SAM,what="try-error")==FALSE)  LIST[[i-1]] <- SAM

              } # end for.
                       
              ###
              idx <- sapply(LIST,function(x) !is.null(x))

              ###
              if (sum(idx)>0)  {

                  BICS <- sapply(LIST[idx], function(x) x$bic)
                  idx2 <- diff(BICS)>0

                  if (sum(idx2)>0) {

                      idx4 <- which(idx)[which(idx2)[1]]

                  } else {
 
                      idx4 <- rev(which(idx))[1]

                  } # end if.

                  ###
                  ###print(idx)
                  ###print(BICS)
                  return(LIST[[idx4]])
 
              } else {

                  stop("Optimization of FMM failed!")

              } # end if. 

          } else {

              fmmED(ed=ed, sed=sed, ncomp=ncomp, addsigma=addsigma, iflog=iflog)

          } # end if.

      } # end function optFMM.
      
      ###
      ###FMM <- numOSL::RadialPlotter(cbind(meanGZ,sdGZ0),ncomp=nfm, maxcomp=13, plot=FALSE, addsigma=addsigma)
      FMM <- optFMM(ed=meanGZ, sed=sdGZ0, ncomp=nfm, maxcomp=9, addsigma=addsigma)

      ###
      nfm <- nrow(FMM$pars)

      ###
      if (logx==TRUE) {

          myTK <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)
          myLB <- c("0.1","0.2","0.5","1","2","5","10","20","50","100","200","500","1000","2000","5000")

      } else {

          myTK <- seq(from=0, to=2000, by=400)
          myLB <- as.character(myTK)

      } # end if. 
     
      ###
      myTK2 <- seq(from=0, to=100, by=20)   

      ###
      if (!is.null(outfile)) {

          mypdfname <- paste(outfile,"_batch_summary.pdf",sep="")
          WPDF <- try(pdf(file=mypdfname, width=7, height=7), silent=TRUE)
          if (inherits(WPDF,what="try-error")==TRUE)  stop(paste("Failed in write to ", mypdfname, ", this may because the file is already opened!",sep=""))

      } else {

          mypdfname <- paste(xoxoxo,"_summary.pdf",sep="")
          WPDF <- try(pdf(file=mypdfname, width=7, height=7), silent=TRUE)
          if (inherits(WPDF,what="try-error")==TRUE)  stop(paste("Failed in write to ", mypdfname, ", this may because the file is already opened!",sep=""))
              
      } # end if. 
      
      ###
      ### Visualize the measured and fitted GSDs.
      ###----------------------------------------
      layout(cbind(c(1,1,2,2),c(1,1,2,2)), respect=TRUE)
      par(mar=c(4.1, 4.1, 2.1, 2.1))

      ###
      gsl <- GSDbatchx[[1]]$gs.comp[,1]
      GSD <- sapply(GSDbatchx, function(x) x$gs.comp[,2])
      fitGSD <- sapply(GSDbatchx, function(x) x$gs.comp[,3])
      
      ### 
      matplot(gsl, GSD, log="x", type="l", lty=1, xlab="Grain size (um)", ylab="Measured volume (%)", 
              col="red", lwd=1.25, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, xaxt="n")   
      axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)

      ### 
      matplot(gsl, fitGSD, log="x", type="l", lty=1, xlab="Grain size (um)", ylab="Fitted volume (%)",
              col="skyblue", lwd=1.25, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, xaxt="n")
      axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)

      ###
      ### THE FIRST PART. 
      ###================================================================================================
      ###************************************************************************************************
      par(mfrow=c(1,1))
      par(mar=c(2.1, 2.1, 4.1, 2.1))

      ###
      plot(1,1,type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main="Summary of unmixing results", cex.main=2)

      ###
      MDL <- sapply(GSDbatchx, function(x) x$model)
      levMDL <- levels(factor(MDL))

      ###
      ### Function lgy() is used for generating legends of models.
      ###---------------------------------------------------------
      lgy <- function(y,yf) {

          nyf <- length(yf)

          ###
          vy <- vector(length=nyf)

          ###
          for (i in 1:nyf) {

              ym <- if (yf[i]=="weibull") "Weibull" else if (yf[i]=="lognormal") "Lognormal" else if
                       (yf[i]=="weibull0")  "Weibull0" else if (yf[i]=="lognormal0") "Lognormal0" else if
                       (yf[i]=="skewnormal0") "SN0" else if (yf[i]=="skewgnormal0") "SGN0" 

              ###
              vy[i] <- paste("Number of GSDs unmixed with ", ym,": ", sum(y==yf[i]), sep="")

          } # end for.

          ###
          return(vy)

      } # end function lgy.
       
      ###
      mylegend <- c(paste("Total number of single-sample GSDs:  ", N, sep=""), "", lgy(MDL, levMDL), "",
                    paste("Number of GSDs with FOM>5%:  ", sum(FOM>5), sep=""),
                    paste("Number of GSDs with RSS>1:  ", sum(RSS>1), sep=""),
                    paste("Number of GSDs with R2<0.995:  ", sum(R2<0.996), sep=""),       
                    paste("Number of GSDs with RSE>0.1:  ", sum(RSE>0.1), sep=""), "",
                    paste("Number of GSDs with ncomp=1 :  ", sum(compvec==1), sep=""),
                    paste("Number of GSDs with ncomp=2 :  ", sum(compvec==2), sep=""),
                    paste("Number of GSDs with ncomp=3 :  ", sum(compvec==3), sep=""),
                    paste("Number of GSDs with ncomp=4 :  ", sum(compvec==4), sep=""),
                    paste("Number of GSDs with ncomp=5 :  ", sum(compvec==5), sep=""),
                    paste("Number of GSDs with ncomp=6 :  ", sum(compvec==6), sep=""),
                    paste("Number of GSDs with ncomp=7 :  ", sum(compvec==7), sep=""),
                    paste("Number of GSDs with ncomp=8 :  ", sum(compvec==8), sep=""),
                    paste("Number of GSDs with ncomp=9 :  ", sum(compvec==9), sep=""),
                    paste("Number of GSDs with ncomp=10:  ", sum(compvec==10), sep=""),
                    paste("Number of GSDs with ncomp=11:  ", sum(compvec==11), sep=""),
                    paste("Number of GSDs with ncomp=12:  ", sum(compvec==12), sep=""),
                    paste("Number of GSDs with ncomp=13:  ", sum(compvec==13), sep=""))

      ###
      legend("top",legend=mylegend[1:(max(compvec)+10)], bty="n", cex=1.25)

      ### 
      layout(matrix(c(1,2,3,4),nrow=2, byrow=TRUE), respect=TRUE)
      par(mar=c(4.1, 4.1, 2.1, 2.1))

      ###
      XY1 <- density(FOM)
      plot(XY1, main="", xlab="Figure of merit (FOM) (%)", ylab="Kernel density", 
           col="skyblue", lwd=3, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, xlim=range(XY1$x))
       
      ###
      par("new"=TRUE)
      plot(sort(FOM), seq(FOM), pch=21, col="black", cex=1.5,  
           xlim=range(XY1$x), xaxt="n", yaxt="n", xlab="", ylab="")

      ###
      qtfom <- round(quantile(FOM, probs=c(0,0.25,0.5,0.75,1)),3)
      legend("right", legend=c(paste("Min: ",qtfom[1],sep=""), 
             paste("25%: ",qtfom[2],sep=""),paste("50%: ",qtfom[3],sep=""),
             paste("75%: ",qtfom[4],sep=""),paste("Max: ",qtfom[5],sep="")), 
             text.col="skyblue", bty="n")

      ###
      XY2 <- density(RSS)
      plot(XY2, main="", xlab="Residual sum of squares (RSS)", ylab="Kernel density", 
           col="red", lwd=3, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, xlim=range(XY2$x))

      ###
      par("new"=TRUE)
      plot(sort(RSS), seq(RSS), pch=22, col="black", cex=1.5, xlim=range(XY2$x), 
           xaxt="n", yaxt="n", xlab="", ylab="")

      ###
      qtrss <- round(quantile(RSS, probs=c(0,0.25,0.5,0.75,1)),3)
      legend("right", legend=c(paste("Min: ",qtrss[1],sep=""), 
             paste("25%: ",qtrss[2],sep=""),paste("50%: ",qtrss[3],sep=""),
             paste("75%: ",qtrss[4],sep=""),paste("Max: ",qtrss[5],sep="")), 
             text.col="red", bty="n")
      
      ###
      XY3 <- density(R2)
      plot(XY3, main="", xlab="R2 statistic (R2)", ylab="Kernel density", 
           col="purple", lwd=3, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, 
           xlim=range(XY3$x))

      ###
      par("new"=TRUE)
      plot(sort(R2), seq(R2), pch=23, col="black", cex=1.5, xlim=range(XY3$x), 
           xaxt="n", yaxt="n", xlab="", ylab="")

      ###
      qtr2 <- round(quantile(R2, probs=c(0,0.25,0.5,0.75,1)),3)
      legend("left", legend=c(paste("Min: ",qtr2[1],sep=""), 
             paste("25%: ",qtr2[2],sep=""),paste("50%: ",qtr2[3],sep=""),
             paste("75%: ",qtr2[4],sep=""),paste("Max: ",qtr2[5],sep="")), 
             text.col="purple", bty="n")

      ###  
      XY4 <- density(RSE)
      plot(XY4, main="", xlab="Residual standard error (RSE)", ylab="Kernel density", 
           col="grey60", lwd=3, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, 
           xlim=range(XY4$x))
      
      ###
      par("new"=TRUE)
      plot(sort(RSE), seq(RSE), pch=24, col="black", cex=1.5,  
           xlim=range(XY4$x), xaxt="n", yaxt="n", xlab="", ylab="")

      ###
      qtrse <- round(quantile(RSE, probs=c(0,0.25,0.5,0.75,1)),3)
      legend("right", legend=c(paste("Min: ",qtrse[1],sep=""), 
             paste("25%: ",qtrse[2],sep=""),paste("50%: ",qtrse[3],sep=""),
             paste("75%: ",qtrse[4],sep=""),paste("Max: ",qtrse[5],sep="")), 
             text.col="grey60", bty="n")

      ###
      ### THE SECOND and THIRD PARTS.
      ###================================================================================================
      ###************************************************************************************************
      colvec <- c("skyblue", "red", "blue", "plum4", "green", "brown", "yellow2", "deeppink",  
                  "black", "wheat", "grey60", "turquoise", "magenta1")

      ###
      for (k in 1:2) {

          if (k==1) kDDD <- medianGZ
          if (k==2) kDDD <- modeGZ

          ###
          KM <- stats::kmeans(log(kDDD), centers=nkm, nstart=1000)

          ###
          kmCST <- round(sort(exp(KM$centers)),2)

          ###
          rangeXV <- meansdABD <- matrix(nrow=13, ncol=2)
          normABD <- nxx <- sdXV <- vector(length=13)

          ###
          for (i in 1:nkm) {

              iinndx <- KM$cluster==i
              nxx[i] <- sum(iinndx)
              meansdABD[i,] <- round(c(mean(abundance[iinndx]), sd(abundance[iinndx])/sqrt(nxx[i])),2)
              rangeXV[i,] <- range(kDDD[iinndx])
              sdXV[i] <- round(sd(kDDD[iinndx])/sqrt(nxx[i]),2)
              normABD[i] <- round(sum(abundance[iinndx])/Nx,2)

          } # end for.

          ###
          odrgXV <- order(rowMeans(rangeXV))
          nxx <- nxx[odrgXV]
          propnxx <- round(nxx/sum(nxx)*100,2)
          meansdABD <- meansdABD[odrgXV,,drop=FALSE]    
          sdXV <- sdXV[odrgXV]
          normABD <- normABD[odrgXV] 

          ###
          ###-------------------------------------------------------------------------------------
          par(mfrow=c(1,1))
          par(mar=c(2.1, 2.1, 4.1, 2.1))

          ###
          plot(1, 1, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", cex.main=2,
               main=paste("KM clustering of unmixed ", ifelse(k==1, "median","mode"),sep=""))

          ###
          mylegend <- c("Algorithm: K-Means clustering","",
                        paste("(1) Average_", ifelse(k==1, "Median","Mode"), "\u00B1", "SD_of_Average_",ifelse(k==1, "Median","Mode"),",",sep=""),
                        paste("(2) Average_Proportion", "\u00B1", "SD_of_Average_Proportion, ",sep=""),
                        "(3) Proportion_of_Data_Points (sum up to 100\u0025),",
                        "(4) Normalized_Abundance (sum up to 100\u0025),",
                        "for each cluster are summarized as follows:","",
                        paste("Cluster1: ", kmCST[1],"\u00B1",sdXV[1],"um, ", meansdABD[1,1],"\u00B1",meansdABD[1,2],"\u0025, ", propnxx[1],"\u0025, ", normABD[1],"\u0025", sep=""),
                        paste("Cluster2: ", kmCST[2],"\u00B1",sdXV[2],"um, ", meansdABD[2,1],"\u00B1",meansdABD[2,2],"\u0025, ", propnxx[2],"\u0025, ",normABD[2],"\u0025",sep=""),
                        paste("Cluster3: ", kmCST[3],"\u00B1",sdXV[3],"um, ", meansdABD[3,1],"\u00B1",meansdABD[3,2],"\u0025, ", propnxx[3],"\u0025, ",normABD[3],"\u0025",sep=""),
                        paste("Cluster4: ", kmCST[4],"\u00B1",sdXV[4],"um, ", meansdABD[4,1],"\u00B1",meansdABD[4,2],"\u0025, ", propnxx[4],"\u0025, ",normABD[4],"\u0025",sep=""),
                        paste("Cluster5: ", kmCST[5],"\u00B1",sdXV[5],"um, ", meansdABD[5,1],"\u00B1",meansdABD[5,2],"\u0025, ", propnxx[5],"\u0025, ", normABD[5],"\u0025",sep=""),
                        paste("Cluster6: ", kmCST[6],"\u00B1",sdXV[6],"um, ", meansdABD[6,1],"\u00B1",meansdABD[6,2],"\u0025, ", propnxx[6],"\u0025, ", normABD[6],"\u0025",sep=""),
                        paste("Cluster7: ", kmCST[7],"\u00B1",sdXV[7],"um, ", meansdABD[7,1],"\u00B1",meansdABD[7,2],"\u0025, ", propnxx[7],"\u0025, ", normABD[7],"\u0025",sep=""),
                        paste("Cluster8: ", kmCST[8],"\u00B1",sdXV[8],"um, ", meansdABD[8,1],"\u00B1",meansdABD[8,2],"\u0025, ", propnxx[8],"\u0025, ", normABD[8],"\u0025",sep=""),
                        paste("Cluster9: ", kmCST[9],"\u00B1",sdXV[9],"um, ", meansdABD[9,1],"\u00B1",meansdABD[9,2],"\u0025, ", propnxx[9],"\u0025, ", normABD[9],"\u0025",sep=""),
                        paste("Cluster10: ", kmCST[10],"\u00B1",sdXV[10],"um, ", meansdABD[10,1],"\u00B1",meansdABD[10,2],"\u0025, ", propnxx[10],"\u0025, ", normABD[10],"\u0025",sep=""),
                        paste("Cluster11: ", kmCST[11],"\u00B1",sdXV[11],"um, ", meansdABD[11,1],"\u00B1",meansdABD[11,2],"\u0025, ", propnxx[11],"\u0025, ", normABD[11],"\u0025",sep=""),
                        paste("Cluster12: ", kmCST[12],"\u00B1",sdXV[12],"um, ", meansdABD[12,1],"\u00B1",meansdABD[12,2],"\u0025, ", propnxx[12],"\u0025, ", normABD[12],"\u0025",sep=""),
                        paste("Cluster13: ", kmCST[13],"\u00B1",sdXV[13],"um, ", meansdABD[13,1],"\u00B1",meansdABD[13,2],"\u0025, ", propnxx[13],"\u0025, ", normABD[13],"\u0025",sep=""))

          ###
          legend("top",legend=mylegend[1:(nkm+8)], bty="n", cex=1.25)
      
          ###
          ###-------------------------------------------------------------------------------
          par(mfrow=c(1,1))
          par(mar=c(4.1, 4.1, 2.1, 2.1))

          ###
          plot(x=kDDD, y=abundance, xaxt="n", yaxt="n", log=ifelse(logx,"x",""), type="n", lwd=3, 
               xlab=paste("The ",ifelse(k==1,"median","mode"), " of grain-size components (um)",  
               sep=""), ylab="Proportion (%)", mgp=c(2.3,1,0), cex.lab=1.5, cex.axis=1.5, cex=2)
          grid(lwd=2, col="grey90", lty=1)
          box(lwd=2)

          ###
          axis(side=1, at=myTK, labels=myLB, cex.axis=1.5, lwd=1.8)
          axis(side=2, at=myTK2, labels=myTK2, cex.axis=1.5, lwd=1.8)

          ###
          for (i in 1:nkm) {

              XV <- kDDD[KM$cluster==odrgXV[i]]
              YV <- abundance[KM$cluster==odrgXV[i]]

              ###
              xxv <- range(XV)
              yyv <- range(YV)
    
              ###
              polygon(x=c(xxv[1],xxv[1],xxv[2],xxv[2]), y=c(yyv[1],yyv[2],yyv[2],yyv[1]), border="grey95", col="grey95")

          } # end for.

          ###
          spreadXV <- exp(seq(from=log(range(gsl)[1]), to=log(range(gsl)[2]), 
                          by=(log(range(gsl)[2])-log(range(gsl)[1]))/(2*length(gsl))))

          ###
          for (i in 1:nkm) {

              XV <- kDDD[KM$cluster==odrgXV[i]]
              YV <- abundance[KM$cluster==odrgXV[i]]

              ###
              points(x=XV, y=YV, col=colvec[i], bg=colvec[i], pch=21, cex=1)

              ###
              ###if (length(XV)>=2) {
                  ###spreadXV <- seq(from=0.1*min(kDDD), to=max(kDDD)*1.9, by=max(diff(sort(XV)))/50)
              ###} else {
                  ###spreadXV <- seq(from=0.1*min(kDDD), to=max(kDDD)*1.9, by=XV/50)
              ###} # end if.

              ###
              normPDF <- dnorm(x=log(spreadXV), mean=mean(log(XV)), sd=sd(log(XV)))

              ###
              points(x=spreadXV, y=max(YV)*normPDF/max(normPDF), col=colvec[i], lwd=3, type="l")

          } # end for.

          ###
          legend(ifelse(logx,"topleft","topright"), legend=paste("Cluster ",1:nkm,sep=""), 
                 pch=21, col=colvec[1:nkm], pt.bg=colvec[1:nkm], bty="n")

          ###
          if (k==1) medianKM <- KM$cluster
          if (k==2) modeKM <- KM$cluster

      } # end for.

      ###
      ### THE FOURTH PART.
      ###================================================================================================
      ###************************************************************************************************
      fmmMAT <- matrix(nrow=length(meanGZ), ncol=nfm)

      ###
      yv_fmm <- log(meanGZ)
      seyv_fmm <- sqrt((sdGZ0/meanGZ)^2 + addsigma^2)

      ###
      for (i in 1:nfm) {

          p_fmm <- FMM$pars[i,1]
          mu_fmm <- log(FMM$pars[i,2])

          ###
          fmmMAT[,i] <- p_fmm/sqrt(2*pi)/seyv_fmm*exp(-0.5*(yv_fmm-mu_fmm)^2/seyv_fmm^2)

      } # end for.

      ###
      cluster_fmm <- apply(fmmMAT, MARGIN=1, which.max) 
      if (!all((1:nfm) %in% cluster_fmm)) {

          dev.off()
          stop(paste("Application of FMM clustering with nfm=",nfm," failed!",sep=""))

      } # end if.

      ###
      kmCST <- round(FMM$pars[,2],2)

      ###
      meansdABD <- matrix(nrow=13, ncol=2)
      normABD <- nxx <- sdXV <- vector(length=13)
      
      ###
      for (i in 1:nfm) {

          iinndx <- cluster_fmm==i
          nxx[i] <- sum(iinndx)
          meansdABD[i,] <- round(c(mean(abundance[iinndx]), sd(abundance[iinndx])/sqrt(nxx[i])),2)
          sdXV[i] <- round(sd(meanGZ[iinndx])/sqrt(nxx[i]),2) 
          normABD[i] <- round(sum(abundance[iinndx])/Nx,2)

      } # end for.

      ###
      propnxx <- round(nxx/sum(nxx)*100,2)
      
      ###
      ###---------------------------------------------------------------------------------
      par(mfrow=c(1,1))
      par(mar=c(2.1, 2.1, 4.1, 2.1))
      plot(1, 1, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main="FMM clustering of unmixed mean", cex.main=2)

      ###
      mylegend <- c("Algorithm: Bayesian clustering (Finite mixture model)", 
                    paste("The maximum logged likelihood value: ", round(FMM$maxlik,1),sep=""),
                    paste("The Bayesian Information Criteria value: ", round(FMM$bic,1),sep=""),"",
                    paste("(1) Average_Mean", "\u00B1", "SD_of_Average_Mean,",sep=""),
                    paste("(2) Average_Proportion", "\u00B1", "SD_of_Average_Proportion",sep=""),
                    "(3) Proportion_of_Data_Points (sum up to 100\u0025),",
                    "(4) Normalized_Abundance (sum up to 100\u0025),",
                    "for each cluster are summarized as follows:","",
                    paste("Cluster1: ", kmCST[1],"\u00B1",sdXV[1],"um, ", meansdABD[1,1],"\u00B1",meansdABD[1,2],"\u0025, ",propnxx[1],"\u0025, ", normABD[1],"\u0025",sep=""),
                    paste("Cluster2: ", kmCST[2],"\u00B1",sdXV[2],"um, ", meansdABD[2,1],"\u00B1",meansdABD[2,2],"\u0025, ",propnxx[2],"\u0025, ", normABD[2],"\u0025",sep=""),
                    paste("Cluster3: ", kmCST[3],"\u00B1",sdXV[3],"um, ", meansdABD[3,1],"\u00B1",meansdABD[3,2],"\u0025, ",propnxx[3],"\u0025, ", normABD[3],"\u0025",sep=""),
                    paste("Cluster4: ", kmCST[4],"\u00B1",sdXV[4],"um, ", meansdABD[4,1],"\u00B1",meansdABD[4,2],"\u0025, ",propnxx[4],"\u0025, ", normABD[4],"\u0025",sep=""),
                    paste("Cluster5: ", kmCST[5],"\u00B1",sdXV[5],"um, ", meansdABD[5,1],"\u00B1",meansdABD[5,2],"\u0025, ",propnxx[5],"\u0025, ", normABD[5],"\u0025",sep=""),
                    paste("Cluster6: ", kmCST[6],"\u00B1",sdXV[6],"um, ", meansdABD[6,1],"\u00B1",meansdABD[6,2],"\u0025, ",propnxx[6],"\u0025, ", normABD[6],"\u0025",sep=""),
                    paste("Cluster7: ", kmCST[7],"\u00B1",sdXV[7],"um, ", meansdABD[7,1],"\u00B1",meansdABD[7,2],"\u0025, ",propnxx[7],"\u0025, ", normABD[7],"\u0025",sep=""),
                    paste("Cluster8: ", kmCST[8],"\u00B1",sdXV[8],"um, ", meansdABD[8,1],"\u00B1",meansdABD[8,2],"\u0025, ",propnxx[8],"\u0025, ", normABD[8],"\u0025",sep=""),
                    paste("Cluster9: ", kmCST[9],"\u00B1",sdXV[9],"um, ", meansdABD[9,1],"\u00B1",meansdABD[9,2],"\u0025, ",propnxx[9],"\u0025, ", normABD[9],"\u0025",sep=""),
                    paste("Cluster10: ", kmCST[10],"\u00B1",sdXV[10],"um, ", meansdABD[10,1],"\u00B1",meansdABD[10,2],"\u0025",propnxx[10],"\u0025, ", normABD[10],"\u0025",sep=""),
                    paste("Cluster11: ", kmCST[11],"\u00B1",sdXV[11],"um, ", meansdABD[11,1],"\u00B1",meansdABD[11,2],"\u0025",propnxx[11],"\u0025, ", normABD[11],"\u0025",sep=""),
                    paste("Cluster12: ", kmCST[12],"\u00B1",sdXV[12],"um, ", meansdABD[12,1],"\u00B1",meansdABD[12,2],"\u0025",propnxx[12],"\u0025, ", normABD[12],"\u0025",sep=""),
                    paste("Cluster13: ", kmCST[13],"\u00B1",sdXV[13],"um, ", meansdABD[13,1],"\u00B1",meansdABD[13,2],"\u0025",propnxx[13],"\u0025, ", normABD[13],"\u0025",sep=""))

      ###
      legend("top",legend=mylegend[1:(nfm+10)], bty="n", cex=1.25)

      ###
      ###-------------------------------------------------------------------------------
      x <- seyv_fmm
      y <- yv_fmm
      centralGZ <- sum(y/x^2)/sum(1/x^2)
   
      ###
      xx <- 1/x
      yy <- (y-centralGZ)/x

      ###
      minPrecision <- 0
      if (is.null(mps)) {

          maxPrecision <-  max(xx)*1.2

      } else {

          if (mps<quantile(xx,probs=0.8))  {

              dev.off()
              stop(paste("Failed in plot the radial plot, try with a larger [mps], for example, mps=",round(max(xx)*1.2,2),"!",sep=""))

          } # end if.

          ###
          if (mps>3*max(xx)) {

              dev.off()
              stop(paste("Failed in plot the radial plot, try with a smalller [mps], for example, mps=",round(max(xx)*1.2,2),"!",sep=""))

          } # end if.

          ###
          maxPrecision <- mps

      } # end if. 

      ###
      zmin <- min(meanGZ)*0.9 
      zmax <- max(meanGZ)*1.1
   
      ###
      miny <- (log(zmin)-centralGZ)*maxPrecision
      maxy <- (log(zmax)-centralGZ)*maxPrecision

      ###
      par("mar"=c(4.1, 4.1, 2.1, 5.1))

      ###
      plot(xx, yy, xaxt="n", yaxt="n", bty="n", xlim=c(minPrecision, maxPrecision), 
           xaxs="i", yaxs="i", ylim=c(miny, maxy), xlab="", ylab="", type="n")
      box(lwd=2)

      ###
      locx <- axTicks(side=1)
      locx <- seq(from=min(locx), to=max(locx), by=(max(locx)-min(locx))/5)
      locy <- seq(from=miny, to=maxy, by=(maxy-miny)/5)

      ###
      axis(side=1, at=locx[-1], lwd=2, labels=round(100/locx[-1],1), 
           line=0, tck=0.02, padj=-4, cex.axis=1.2)
      mtext(text="Relative standard error", side=1, cex=1.2, line=-3)

      ###
      axis(side=1, at=locx, lwd=2.0, labels=locx, line=0, padj=0, cex.axis=1.2)
      mtext(text="Precision", side=1, cex=1.2, line=2)

      ###
      axis(side=2, at=c(-2,0,2), lwd=2, labels=FALSE, cex.axis=1.2)
      mtext(text="-2", side=2, at=-2, cex=0.8, line=1)
      mtext(text="0",side=2, at=0, cex=0.8, line=1)
      mtext(text="2", side=2, at=2, cex=0.8, line=1)
      mtext(text="Standardised Estimate", side=2, at=0, cex=1.2, line=2)

      ###
      axis(side=4, at=locy, labels=round(exp(locy/maxPrecision+centralGZ),2), lwd=2, cex.axis=1.2)
      abline(v=maxPrecision, lwd=2)
      mtext(text="The mean of grain-size components (um)", side=4, cex=1.2, line=3) 
   
      ###
      for (i in 1:nfm) {    
  
          y1 <- (log(FMM$pars[i,2])-centralGZ)*maxPrecision

          ###
          polygon(x=rep(c(minPrecision, maxPrecision), each=2), y=c(-2, 2, y1+2, y1-2), col="grey95", lty="blank")

      } # end for.

      ###
      for (i in 1:nfm) {   

          y1 <- (log(FMM$pars[i,2])-centralGZ)*maxPrecision

          ###
          lines(c(minPrecision, maxPrecision), c(0,y1), lwd=2, lty=2, col=colvec[i]) 
          
          ###
          points(xx[cluster_fmm==i], yy[cluster_fmm==i], pch=21, col=colvec[i], bg=colvec[i], cex=1)

      } # end for.

      ###
      legend("topleft", legend=paste("Cluster ",1:nfm,sep=""), pch=21, col=colvec[1:nfm], pt.bg=colvec[1:nfm], bty="n")

      ###
      ###-------------------------------------------------------------------------------
      par(mfrow=c(1,1))
      par(mar=c(4.1, 4.1, 2.1, 2.1))

      ###
      plot(x=meanGZ, y=abundance, xaxt="n", yaxt="n", log=ifelse(logx,"x",""), type="n",  
           lwd=3, xlab="The mean of grain-size components (um)", ylab="Proportion (%)", 
           mgp=c(2.3,1,0), cex.lab=1.5, cex.axis=1.5, cex=2)
      grid(lwd=2, col="grey90", lty=1)
      box(lwd=2)

      ###
      axis(side=1, at=myTK, labels=myLB, cex.axis=1.5, lwd=1.8)
      axis(side=2, at=myTK2, labels=myTK2, cex.axis=1.5, lwd=1.8) 

      ###
      for (i in 1:nfm) {

          XV <- meanGZ[cluster_fmm==i]
          YV <- abundance[cluster_fmm==i]
          seXV <- sdGZ0[cluster_fmm==i]

          ### Scale "seXV" using "addsigma".
          if (addsigma!=0)  seXV <- XV*sqrt((seXV/XV)^2+addsigma^2)
  
          ###
          nXV <- length(XV)

          ###
          xxv <- range(XV)
          yyv <- range(YV)

          ###
          points(x=XV, y=YV, col=colvec[i], bg=colvec[i], pch=21, cex=1)

          ###
          arrows(x0=XV-seXV/2, y0=YV, x1=XV+seXV/2, y1=YV, code=3, lwd=2, angle=90, length=0.05, col=colvec[i])
             
          ###
          ###if (nXV>=2) {
              ###spreadXV <- seq(from=0.1*min(meanGZ), to=max(meanGZ)*1.9, by=max(diff(sort(XV)))/50)
          ###} else {
              ###spreadXV <- seq(from=0.1*min(meanGZ), to=max(meanGZ)*1.9, by=XV/50)
          ###} # end if.

          ###      
          pdfMat <- matrix(nrow=length(spreadXV), ncol=nXV)

          ###
          for(j in 1:nXV) {

              pdfMat[,j] <- dnorm(x=log(spreadXV), mean=log(XV[j]), sd=seXV[j]/XV[j], log=FALSE)

          } # end if.

          ###
          pdfXV <- rowSums(pdfMat)

          ###
          points(spreadXV,max(YV)*pdfXV/max(pdfXV),type="l", col=colvec[i], lwd=3)

      } # end for.

      ###
      legend(ifelse(logx,"topleft","topright"), legend=paste("Cluster ",1:nfm,sep=""), 
             pch=21, col=colvec[1:nfm], pt.bg=colvec[1:nfm], bty="n")

      ###
      dev.off()

      ###
      gsfit <- data.frame((seq(N))[ookk], MDL, compvec, FOM, RSS, R2, RSE)

      ###
      rownames(gsfit) <- NULL
      colnames(gsfit) <- c("NO","Model","ncomp","FOM", "RSS", "R2", "RSE")
    
      ###
      pgsp <- cbind(rep((seq(N))[ookk],times=compvec), abundance, meanGZ, medianGZ, modeGZ, sdGZ0, cluster_fmm, medianKM, modeKM)

      ###
      rownames(pgsp) <- NULL
      colnames(pgsp) <- c("NO", "Proportion", "Mean", "Median", "Mode", "Sd", "Mean_Cluster", "Median_Cluster","Mode_Cluster")

      ###
      output <- list("gsfit"=gsfit, "pgsp"=pgsp)
     
      ###
      if (!is.null(outfile)) {

          myfilename1 <- paste(outfile,"_batch_gsfit.csv",sep="")
          WF1 <- try(write.csv(gsfit, file=myfilename1), silent=TRUE)
          if (inherits(WF1,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename1, ", this may because the file is already opened!",sep=""))
          
          ###
          myfilename2 <- paste(outfile,"_batch_pgsp.csv",sep="")
          WF2 <- try(write.csv(pgsp, file=myfilename2), silent=TRUE)
          if (inherits(WF2,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename2, ", this may because the file is already opened!",sep=""))

      } else {

          myfilename1 <- paste(xoxoxo,"_gsfit.csv",sep="")
          WF1 <- try(write.csv(gsfit, file=myfilename1), silent=TRUE)
          if (inherits(WF1,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename1, ", this may because the file is already opened!",sep=""))

          ###
          myfilename2 <- paste(xoxoxo,"_pgsp.csv",sep="")
          WF2 <- try(write.csv(pgsp, file=myfilename2), silent=TRUE)
          if (inherits(WF2,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename2, ", this may because the file is already opened!",sep=""))

      } # end if.

      ###
      invisible(output)

  } # end function summary_ussGSDbatch.
  ###==========================================================================================================================
  ###
  
  ###**************************************************************************************************************************
  ### Function update_ussGSDbatch() is used to update (re-analyze) the   
  ### unmixing results of some specified grain-size distributions. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ### obj_batchgsd: an S3 object of class "gsdbatch" generated using the function 
  ###               ussGSDbatch(), ussGSDbatchp(), update_ussGSDbatch(), or update_ussGSDbatchp().
  ###
  ###     sampleNO: An integer vector indicating the ID numbers of the grain-size  
  ###              distributions to be updated (re-analyzed). 
  ###
  ###        ncomp: An integer (with values range from 1 to 13) indicating the number of 
  ###               components for the grain-size distributions to be re-analyzed. 
  ###
  ###         auto: A logical value indicating whether performing a automatic grain-size unmixing, 
  ###               in this case the user needs not to specify the initials used for optimization.
  ###
  ###        model: A character indicating the model to be fitted, "weibull", "lognormal",
  ###               "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###               if model=NULL, the program will automatically determine the optimal model 
  ###               (between "weibull" and "lognormal") which yields a smaller FOM or RSS value.
  ###
  ###          mpd: A real number indicating the allowed lower limit of the volume percentage 
  ###               of a GSD that can be used to identify the location of component peaks.
  ###
  ###          ctf: A numeric value (between -1 and 1) representing a critical threshold factor  
  ###               that controls the identification of peaks from the grain-size distribution.  
  ###               This argument can be used to prevent identifying a false peak characterised  
  ###               by a large positive second-order derivative. Specially, ctf=0 indicates
  ###               that peaks with second-order derivatives above zero will be precluded. 
  ###               The degree of precluding will be increasingly suppressed as ctf increases
  ###               from -1 to 1. 
  ###
  ###         ntry: An integer value indicating the number of trials in a trial-and-error protocol. 
  ###
  ###          kkf: A numeric value controlling the range of values from which random starting   
  ###               parameters will be generated during the "trial-and-error" protocol.
  ###
  ###    startPars: if model="weibull" or "lognormal", [startPars] should be a two-column matrix 
  ###               containing starting parameters used for unmixing, the first column contains 
  ###               the modes of individual components, and the second row contains the alpha 
  ###               ("weibull") or sigma ("lognormal") value of individual components.
  ###
  ###               if model="weibull0" or "lognormal0", [startPars] should be a three-column 
  ###               matrix containing starting parameters used for unmixing, the first column  
  ###               contains the modes of individual components, the second column contains the 
  ###               maximum volume percentages of individual components, and the thrid column 
  ###               contains the alpha ("weibull0") or sigma ("lognormal0") value of individual 
  ###               components.
  ###
  ###               if model="skewnormal0" or "skewgnormal0", [startPars] should be a four-column 
  ###               matrix containing starting parameters used for unmixing, the first column is
  ###               the modes of individual components, the second column is the maximum volume
  ###               percentages of individual components,the thrid and fourth columns are the
  ###               alpha and omega values ("skewnormal0") or sigma and q values ("skewgnormal0") 
  ###               of individual components.  
  ###
  ###   alphaRange: A two-element vector indicating the lower and upper limits on alpha values of 
  ###               the Weibull ("weibull" or "weibull0") or Skew Normal ("skewnormal0") distributions 
  ###               generated from a Uniform distribution.
  ###
  ###   sigmaRange: A two-element vector indicating the lower and upper limits on sigma values 
  ###               of the Lognormal ("lognormal" or "lognormal0") or Skewed Generalized Normal 
  ###              ("skewgnormal0") distributions generated from a Uniform distribution.  
  ###
  ###   omegaRange: A two-element vector indicating the lower and upper limits on omega values of 
  ###               a Skew Normal ("skewnormal0") distribution generated from a Uniform distribution.
  ###
  ###       qRange: A two-element vector indicating the lower and upper limits on q values of a Skewed 
  ###               Generalized Normal ("skewgnormal0") distribution generated from a Uniform distribution.           
  ###
  ###     useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###               used as the independent variable during the unmixing process.
  ###
  ###      minfunc: A character indicating the objective to be mimimized, either "fom" or "rss", 
  ###               for the figure-of-merit value or the residual sum of squares, respectively.
  ###
  ###         trim: A logical value indicating whether the unmixing results will be trimed  
  ###               using a recursive optimization protocol.
  ###
  ###         mrsl: A numeric value indicating the minimum resolution between two adjacent     
  ###               components used to trim the unmixing results, which ranges from 0 to 2.
  ###               Peaks with a second-order derivative below zero and a resolution smaller
  ###               than mrsl will be deleted.
  ###
  ###       rmZero: A logical value indicating whether zero volume percentages will be removed.
  ###
  ###       saveRD: A logical value indicating whether the updated results will be saved to
  ###               the RDdata file in the current directory.
  ###==========================================================================================================================
  ### The function
  ### (1) Returns a invisible list of S3 class of "batchgsd" containing the unmixed results of individual samples.
  ###
  ### (2) Re-write the loaded RData file containing the updated unmixing results 
  ###     in the current working directory if obj_batchgsd=NULL and saveRD=TRUE.
  ###
  ### Note that if obj_batchgsd=NULL, the user needs to ensure that the function load_ussGSDbatch() 
  ### has been called to import a RData file containing an object of S3 class "batchgsd" that is 
  ### available from the current working directory.
  ###==========================================================================================================================
  update_ussGSDbatch <- function(obj_batchgsd=NULL, sampleNO, ncomp=NULL, auto=FALSE, model="weibull0", 
                                 mpd=0, ctf=0.1, ntry=50, kkf=0.1, startPars=NULL, alphaRange=NULL,  
                                 sigmaRange=NULL, omegaRange=NULL, qRange=NULL, useIndex=TRUE, minfunc="fom",
                                 trim=FALSE, mrsl=0.6, rmZero=TRUE, saveRD=TRUE) {

      stopifnot(is.numeric(sampleNO), 
                is.null(ncomp) || (length(ncomp)==1 && is.numeric(ncomp) && ncomp %in% 1:13),
                is.logical(auto), length(auto)==1,
                is.null(model) || (length(model)==1 && model %in% c("weibull","lognormal","weibull0", "lognormal0", "skewnormal0", "skewgnormal0")),
                is.numeric(mpd), length(mpd)==1, mpd>=0, mpd<3,
                is.numeric(ctf), length(ctf)==1, ctf>=-1, ctf<=1,
                is.numeric(ntry), length(ntry)==1, ntry>=2, 
                is.numeric(kkf), length(kkf)==1, kkf>0, kkf<1,
                is.null(startPars) || (is.matrix(startPars) && ncol(startPars)==2), 
                is.null(alphaRange) || (is.numeric(alphaRange) && length(alphaRange)==2),
                is.null(sigmaRange) || (is.numeric(sigmaRange) && length(sigmaRange)==2),
                is.null(omegaRange) || (is.numeric(omegaRange) && length(omegaRange)==2),
                is.null(qRange)|| (is.numeric(qRange) && length(qRange)==2),
                is.logical(useIndex), length(useIndex)==1,
                is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "rss"),
                is.logical(trim), length(trim)==1,
                is.numeric(mrsl), length(mrsl)==1, mrsl>=0, mrsl<=2,
                is.logical(rmZero), length(rmZero)==1,  
                is.logical(saveRD), length(saveRD)==1)


      if (!is.null(obj_batchgsd)) {

          if(inherits(obj_batchgsd,what="batchgsd")==FALSE)  stop("Error: [obj_batchgsd] should be an S3 object of class 'batchgsd'!")
          GSDbatch <- obj_batchgsd

      } else {

          if (!exists("gsdbatch"))  stop("Error: function [load_ussGSDbatch] has not been called!")
          GSDbatch <- get("GSDbatch", envir=gsdbatch)

      } # end if.

      ###
      xoxoxo <- attr(GSDbatch, "xoxoxo")

      ###
      allmyNO <- 1:length(GSDbatch)

      ###
      if (!all(sampleNO %in% allmyNO)) {

          XNO <- sampleNO[!(sampleNO %in% allmyNO)]

          ###
          cat("The following sample numbers are invalid:\n")
          print(XNO)

          ###
          stop()

      } # end if.

      ###
      N <- length(sampleNO)
      plot <- ifelse(N==1, TRUE, FALSE)

      ###
      for (i in seq(N)) {

          ###
          cat("The ", sampleNO[i], "-th single-sample GSD.\n",sep="")

          ###
          gsl <- GSDbatch[[sampleNO[i]]]$gs.comp[,1]
          gsd <- GSDbatch[[sampleNO[i]]]$gs.comp[,2]

          ###
          if (!is.null(model)) {

              if (model %in% c("weibull", "lognormal")) {
     
                  res_ussGSD <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                    viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=TRUE, outfile=NULL, 
                    rmZero=rmZero, plot=plot, sampleName=paste("sampleNO=",sampleNO[i],sep="")), 
                    silent=TRUE)

              } else if (model %in% c("weibull0", "lognormal0", "skewnormal0", "skewgnormal0")) {

                  res_ussGSD <- try(ussGSD0(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, omegaRange=omegaRange, qRange=qRange, useIndex=useIndex, 
                    minfunc=minfunc, trim=trim, mrsl=mrsl, viewAutoInis=FALSE, viewLM=FALSE, 
                    viewFit=FALSE, viewPb=TRUE, outfile=NULL, rmZero=rmZero, plot=plot, 
                    sampleName=paste("sampleNO=",sampleNO[i],sep="")), silent=TRUE)
              
              } # end if.
  
              ###
              if (inherits(res_ussGSD,what="try-error")==FALSE) { 

                  GSDbatch[[sampleNO[i]]] <- res_ussGSD 

              } else {

                  print(res_ussGSD)

              } # end if. 

          } else {

              res_ussGSDwb <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model="weibull", 
                mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=TRUE, outfile=NULL, 
                rmZero=rmZero, plot=plot, sampleName=paste("sampleNO=",sampleNO[i],sep="")), 
                silent=TRUE)

              ###
              res_ussGSDlg <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model="lognormal", 
                mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=TRUE, outfile=NULL, 
                rmZero=rmZero, plot=plot, sampleName=paste("sampleNO=",sampleNO[i],sep="")), 
                silent=TRUE)

              ###
              if (inherits(res_ussGSDwb,what="try-error")==FALSE && inherits(res_ussGSDlg,what="try-error")==FALSE) {

                  if (minfunc=="fom") {

                      minfuncVAL1 <- res_ussGSDwb$FOM
                      minfuncVAL2 <- res_ussGSDlg$FOM

                  } else if (minfunc=="rss") {

                      minfuncVAL1 <- res_ussGSDwb$RSS
                      minfuncVAL2 <- res_ussGSDlg$RSS

                  } # end if.

                  ###
                  if (minfuncVAL1<minfuncVAL2) {

                      GSDbatch[[sampleNO[i]]] <- res_ussGSDwb

                  } else {

                      GSDbatch[[sampleNO[i]]] <- res_ussGSDlg

                  } # end if. 


              } else if (inherits(res_ussGSDwb,what="try-error")==FALSE) {

                  GSDbatch[[sampleNO[i]]] <- res_ussGSDwb

                  ###
                  cat("Failed in unmixing the ", sampleNO[i], "-th sample using Lognormal!\n", sep="")
                  print(paste("<",sampleNO[i],"> ",attr(res_ussGSDlg,"condition"),sep=""))

              } else if (inherits(res_ussGSDlg,what="try-error")==FALSE) {

                  GSDbatch[[sampleNO[i]]] <- res_ussGSDlg

                  ###
                  cat("Failed in unmixing the ", sampleNO[i], "-th sample using Weibull!\n", sep="")
                  print(paste("<",sampleNO[i],"> ",attr(res_ussGSDwb,"condition"),sep=""))

              } else {

                  cat("Failed in unmixing the ", sampleNO[i], "-th sample using Weibull and Lognormal!\n", sep="")
                  print(paste("<",sampleNO[i],"> ",attr(res_ussGSDwb,"condition"),sep=""))
                  print(paste("<",sampleNO[i],"> ",attr(res_ussGSDlg,"condition"),sep=""))

              } # end if.

          } # end if.

      } # end for.  

      ###
      assign("GSDbatch", GSDbatch, envir=gsdbatch)
      if (saveRD==TRUE) save(GSDbatch, file=paste(xoxoxo, ".RData", sep="")) 

      ###
      invisible(GSDbatch)

  } # end function update_ussGSDbatch.
  ###==========================================================================================================================
  ###


  ###
  ###**************************************************************************************************************************
  ### Function plot_ussGSDbatch() is used for visualizing the unmixing results
  ### of grain-size distributions obtained from a batch model using a PDF file. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ### obj_batchgsd: an S3 object of class "gsdbatch" generated using the function 
  ###               ussGSDbatch(), ussGSDbatchp(), update_ussGSDbatch(), or update_ussGSDbatchp().
  ###
  ###      pdfName: A character indicating the name of the PDF file to be generated.
  ###               A default PDF name with postfix "_plot" will be used if pdfName=NULL. 
  ###
  ###        addvl: A logical value indicating whether vertical lines will be added to  
  ###               the unmixed grain-size components visualized using a pdf file. 
  ###
  ###        logxy: A character indicating whether the x- or y-axis will be logged in the plot, 
  ###               one of "", "x", "y", "xy", or NULL.
  ###
  ###          lwd: A numeric value giving the widths of lines in the PDF file.
  ###       
  ###          pch: An integer giving the type of symbols in the PDF file.
  ### 
  ###          cex: A numeric value giving the size of symbols in the PDF file.
  ###==========================================================================================================================
  ### The function generates a PDF file in the current working directory.
  ###
  ### Note that if obj_batchgsd=NULL, the user needs to ensure that the function load_ussGSDbatch() 
  ### has been called to import a RData file containing an object of S3 class "batchgsd" that is 
  ### available from the current working directory.
  ###==========================================================================================================================
  plot_ussGSDbatch <- function(obj_batchgsd=NULL, pdfName=NULL, addvl=TRUE, logxy="x", lwd=3, pch=21, cex=2) {

      if (!is.null(obj_batchgsd)) {

          if(inherits(obj_batchgsd,what="batchgsd")==FALSE)  stop("Error: [obj_batchgsd] should be an S3 object of class 'batchgsd'!")
          GSDbatch <- obj_batchgsd

      } else {

          if (!exists("gsdbatch"))  stop("Error: function [load_ussGSDbatch] has not been called!")
          GSDbatch <- get("GSDbatch", envir=gsdbatch)

      } # end if.

      ###
      N <- length(GSDbatch)
      sampleName <- names(GSDbatch)

      ###
      if (is.null(pdfName)) {

          xoxoxo <- attr(GSDbatch, "xoxoxo")

          ###
          myfilename <- paste(xoxoxo,".pdf",sep="")
          WF <- try(pdf(myfilename), silent=TRUE)
          if (inherits(WF,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename, ", this may because the file is already opened!",sep=""))

      } else {

          myfilename <- paste(pdfName,".pdf",sep="")
          WF <- try(pdf(myfilename), silent=TRUE)
          if (inherits(WF,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename, ", this may because the file is already opened!",sep=""))

      } # end if.

      ###
      if (any(sapply(GSDbatch,is.null))==TRUE) {

          oldpar <- par("mar", "mfrow")
          on.exit(par(oldpar))

      } # end if.

      ###
      for (j in 1:N) {

          sampleNamej <- paste(sampleName[j]," [NO=",j,"]",sep="")

          if (is.null(GSDbatch[[j]])) { 

              ###
              layout(matrix(c(1,1,1,2,1,1,1,3,1,1,1,4),ncol=3), respect=TRUE)
              par(mar=c(4,5,5,1.5)+0.1)
              plot(1, 1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", main=sampleNamej, cex.main=2.2) 
              legend("center", legend="Unmixing failed", yjust=2, ncol=1, cex=1.5, bty="n")

              ###
              par(mar=c(4,4,1,1)+0.1)
              plot(1, 1, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
              plot(1, 1, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
              plot(1, 1, type="n", xaxt="n", yaxt="n", xlab="", ylab="")

          } else {

              plot_ussGSD(obj_ussgsd=GSDbatch[[j]], sampleName=sampleNamej, addvl=addvl, logxy=logxy, lwd=lwd, pch=pch, cex=cex)

          } # end if. 

      } # end for. 

      ###
      dev.off()

  } # end function plot_ussGSDbatch.
  ###==========================================================================================================================
  ###

  ###
  ###**************************************************************************************************************************
  ### Function load_ussGSDbatch() is used for loading a existing dataset and creating 
  ### a new environment called "gsdbatch" to hold the loaded object "GSDbatch". 
  ###==========================================================================================================================
  ### The function contains the following arugement.
  ###
  ### obj_batchgsd: an S3 object of class "gsdbatch" generated using the function
  ###               ussGSDbatch(), ussGSDbatchp(), update_ussGSDbatch(), or update_ussGSDbatchp(), 
  ###               or a character indicating the name of the RData file generated using the function
  ###               ussGSDbatch() or ussGSDbatchp() that is available from the current working directory.
  ###==========================================================================================================================
  ### The function creats a global environment called "gsdbatch" used  for holding the loaded dataset.
  ###==========================================================================================================================
  load_ussGSDbatch <- function(obj_batchgsd) {

      if (!is.character(obj_batchgsd)) {

          if(inherits(obj_batchgsd,what="batchgsd")==FALSE)  stop("Error: [obj_batchgsd] should be an S3 object of class 'batchgsd'!")
          GSDbatch <- obj_batchgsd

      } else {

          file_name <- list.files()

          ###
          myFILE <- obj_batchgsd

          ### 
          if (!myFILE %in% file_name)  stop(paste("Cannot find file [",myFILE,"] in the current working directory!",sep=""))

          ###
          load(myFILE)

          ###
          if (!exists("GSDbatch"))  stop("Error: cannot find object [GSDbatch] from the loaded RData!")

          ###
          if(inherits(GSDbatch,what="batchgsd")==FALSE)  stop("Error: the loaded RData should be an S3 object of class 'batchgsd'!")

      } # end if.
      ###

      ###
      cat("Note: a new environment called [gsdbatch] is created to hold the loaded RData!\n",sep="")
      gsdbatch <<- new.env()

      ###
      assign("GSDbatch", GSDbatch, envir=gsdbatch)

  } # end function load_ussGSDbatch.
  ###==========================================================================================================================
  ###

  ###
  ###**************************************************************************************************************************
  ### Function ussGSDbatch() is used for performing the unmixing of single-sample 
  ### grain-size distributions for a number of samples in a batch pattern. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###   GSDfile: A CSV file containing the grain-size data used for analysis that is available from
  ###            the current working directory, or a data.frame (or matrix) imported using function 
  ###            read.table(). The first row is the grain-size levels, and the remaining rows are the 
  ###            volume percentages of individual samples to be unmixed. If missing, a CSV template 
  ###            will be generated automatically to guide the user to prepare the input dataset.
  ###
  ###     ncomp: An integer (from 1 to 13) indicating the number of components to be fitted. 
  ###
  ###      auto: A logical value indicating whether performing a automatic grain-size unmixing, 
  ###            in this case the user needs not to specify the initials used for optimization.
  ###
  ###     model: A character indicating the model to be fitted, "weibull", "lognormal",
  ###            "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###            if model=NULL, the program will automatically determine the optimal model 
  ###            (between "weibull" and "lognormal") which yields a smaller FOM or RSS value.
  ###
  ###       mpd: A real number indicating the allowed lower limit of the volume percentage 
  ###            of a GSD that can be used to identify the location of component peaks.
  ###
  ###       ctf: A numeric value (between -1 and 1) representing a critical threshold factor  
  ###            that controls the identification of peaks from the grain-size distribution.  
  ###            This argument can be used to prevent identifying a false peak characterized  
  ###            by a large positive second-order derivative. Specially, ctf=0 indicates
  ###            that peaks with second-order derivatives above zero will be precluded. 
  ###            The degree of precluding will be increasingly suppressed as ctf increases 
  ###            from -1 to 1. 
  ###
  ###      ntry: An integer indicating the number of trials in a trial-and-error protocol. 
  ###
  ###       kkf: A numeric value controlling the range of values from which random starting   
  ###            parameters will be generated during the "trial-and-error" protocol.
  ###
  ### startPars: if model="weibull" or "lognormal", [startPars] should be a two-column matrix 
  ###            containing starting parameters used for unmixing, the first column contains 
  ###            the modes of individual components, and the second row contains the alpha 
  ###            ("weibull") or sigma ("lognormal") value of individual components.
  ###
  ###            if model="weibull0" or "lognormal0", [startPars] should be a three-column 
  ###            matrix containing starting parameters used for unmixing, the first column  
  ###            contains the modes of individual components, the second column contains the 
  ###            maximum volume percentages of individual components, and the thrid column 
  ###            contains the alpha ("weibull0") or sigma ("lognormal0") value of individual 
  ###            components.
  ###
  ###            if model="skewnormal0" or "skewgnormal0", [startPars] should be a four-column 
  ###            matrix containing starting parameters used for unmixing, the first column is
  ###            the modes of individual components, the second column is the maximum volume
  ###            percentages of individual components,the thrid and fourth columns are the
  ###            alpha and omega values ("skewnormal0") or sigma and q values ("skewgnormal0") 
  ###            of individual components.  
  ###
  ###alphaRange: A two-element vector indicating the lower and upper limits on alpha values of 
  ###            the Weibull ("weibull" or "weibull0") or Skew Normal ("skewnormal0") distributions 
  ###            generated from a Uniform distribution.
  ###
  ###sigmaRange: A two-element vector indicating the lower and upper limits on sigma values 
  ###            of the Lognormal ("lognormal" or "lognormal0") or Skewed Generalized Normal 
  ###           ("skewgnormal0") distributions generated from a Uniform distribution.  
  ###
  ###omegaRange: A two-element vector indicating the lower and upper limits on omega values of 
  ###            a Skew Normal ("skewnormal0") distribution generated from a Uniform distribution.
  ###
  ###    qRange: A two-element vector indicating the lower and upper limits on q values of a Skewed 
  ###            Generalized Normal ("skewgnormal0") distribution generated from a Uniform distribution.
  ###
  ###  useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###            used as the x-coordinate during the fitting process.
  ###
  ###   minfunc: A character indicating the objective to be mimimized, either "fom" or "rss",
  ###            for the figure-of-merit value or the residual sum of squares, respectively.
  ###
  ###      trim: A logical value indicating whether the unmixing results will be trimed  
  ###            using a recursive optimization protocol.
  ###
  ###      mrsl: A numeric value indicating the minimum resolution between two adjacent     
  ###            components used to trim the unmixing results, which ranges from 0 to 2.
  ###            Peaks with a second-order derivative below zero and a resolution smaller
  ###            than mrsl will be deleted.
  ###
  ###    rmZero: A logical value indicating whether removing zeros from the volume percentages.
  ### 
  ###   outfile: A character indicating the name of the PDF/RData file the unmixing results will be written to.  
  ###            The results will be returned to files with a default name if outfile=NULL. 
  ###
  ###     addvl: A logical value indicating whether vertical lines will be added to the unmixed 
  ###            grain-size components when visualizing the results using a PDF file. 
  ###
  ###     logxy: A character indicating whether the x- or y-axis will be logged in the PDF file, 
  ###            one of "", "x", "y", "xy", or NULL.
  ###
  ###       lwd: A numeric value giving the widths of lines in the PDF file.
  ###       
  ###       pch: An integer giving the type of symbols in the PDF file.
  ### 
  ###       cex: A numeric value giving the size of symbols in the PDF file.
  ###==========================================================================================================================
  ### The function
  ###
  ### (1) Returns a invisible list of S3 object of class "batchgsd" containing the unmixing results of individual samples.
  ###
  ### (2) Generates a PDF file showing the unmixed results in the current working directory.
  ###
  ### (3) Generates a RData file containing the unmixed results in the current working directory.
  ###==========================================================================================================================
  ussGSDbatch <- function(GSDdata="inputGSD.csv", ncomp=NULL, auto=TRUE, model="weibull0", mpd=0, 
                          ctf=0.1, ntry=50, kkf=0.1, startPars=NULL, alphaRange=NULL, sigmaRange=NULL, 
                          omegaRange=NULL, qRange=NULL, useIndex=TRUE, minfunc="fom", trim=TRUE, 
                          mrsl=0.6, rmZero=TRUE, outfile=NULL, addvl=TRUE, logxy="x", lwd=3, 
                          pch=21, cex=2) {
      ###
      stopifnot(is.null(ncomp) || (length(ncomp)==1 && is.numeric(ncomp) && ncomp %in% 1:13),
                is.logical(auto), length(auto)==1,
                is.null(model) || (length(model)==1 && model %in% c("weibull","lognormal","weibull0", "lognormal0", "skewnormal0", "skewgnormal0")),
                is.numeric(mpd), length(mpd)==1, mpd>=0, mpd<3,
                is.numeric(ctf), length(ctf)==1, ctf>=-1, ctf<=1,
                is.numeric(ntry), length(ntry)==1, ntry>=2, 
                is.numeric(kkf), length(kkf)==1, kkf>0, kkf<1,
                is.null(startPars) || (is.matrix(startPars) && ncol(startPars)==2), 
                is.null(alphaRange) || (is.numeric(alphaRange) && length(alphaRange)==2),
                is.null(sigmaRange) || (is.numeric(sigmaRange) && length(sigmaRange)==2),
                is.null(omegaRange) || (is.numeric(omegaRange) && length(omegaRange)==2),
                is.null(qRange)|| (is.numeric(qRange) && length(qRange)==2),
                is.logical(useIndex), length(useIndex)==1,
                is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "rss"),
                is.logical(trim), length(trim)==1,
                is.numeric(mrsl), length(mrsl)==1, mrsl>=0, mrsl<=2,
                is.logical(rmZero), length(rmZero)==1,  
                is.null(outfile) || (is.character(outfile) && length(outfile)==1),
                is.logical(addvl), length(addvl)==1, 
                is.character(logxy), length(logxy)==1,
                is.numeric(lwd), length(lwd)==1, 
                is.numeric(pch), length(pch)==1,
                is.numeric(cex), length(cex)==1)

      ###
      if (is.character(GSDdata) && file.exists(GSDdata)==FALSE) {

          TempTable <- data.frame(rbind(c("GSL",exp(log(0.02)+(log(2000)-log(0.02))/100*(0:100))),
          c("GSD1",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0.12,0.18,0.22,0.26,0.29,0.32,0.34,0.35,0.35,0.36,0.36,0.36,
            0.37,0.4,0.43,0.47,0.53,0.58,0.64,0.7,0.75,0.82,0.89,0.96,1.04,1.12,1.2,1.27,1.34,1.4,1.47,1.53,1.61,1.7,1.81,1.96,2.14,
            2.37,2.65,2.97,3.32,3.69,4.06,4.41,4.71,4.92,5.03,5.02,4.86,4.56,4.13,3.6,2.98,2.35,1.7,1.18,0.75,0.06,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD2",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.07,0.09,0.11,0.13,0.14,0.15,0.15,0.15,0.15,0.15,0.15,
            0.15,0.15,0.17,0.19,0.21,0.24,0.28,0.3,0.33,0.34,0.35,0.36,0.36,0.35,0.35,0.35,0.35,0.36,0.37,0.4,0.42,0.45,0.48,
            0.5,0.49,0.47,0.42,0.36,0.31,0.29,0.36,0.54,0.89,1.43,2.2,3.16,4.28,5.47,6.62,7.6,8.3,8.63,8.51,7.95,7.02,5.81,
            4.45,3.13,1.71,0.32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD3",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.06,0.08,0.1,0.12,0.13,0.14,0.14,0.14,0.14,0.13,0.13,
            0.13,0.13,0.14,0.16,0.18,0.2,0.23,0.25,0.27,0.29,0.3,0.31,0.31,0.32,0.32,0.33,0.35,0.37,0.4,0.43,0.47,0.5,0.51,
            0.51,0.47,0.41,0.33,0.25,0.21,0.24,0.4,0.73,1.28,2.04,3.04,4.18,5.43,6.63,7.68,8.43,8.8,8.71,8.18,7.24,6.03,
            4.66,3.31,1.97,0.62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD4",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.09,0.13,0.16,0.19,0.22,0.24,0.25,0.26,0.27,0.27,0.27,0.28,
            0.29,0.31,0.35,0.39,0.43,0.48,0.53,0.57,0.62,0.66,0.71,0.76,0.82,0.88,0.94,1.01,1.07,1.14,1.2,1.27,1.34,1.42,1.5,
            1.59,1.7,1.83,2,2.21,2.46,2.78,3.15,3.58,4.03,4.49,4.91,5.25,5.47,5.52,5.38,5.05,4.53,3.88,3.14,2.39,1.67,1.15,
            0.39,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD5",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.09,0.13,0.16,0.19,0.21,0.23,0.24,0.25,0.25,0.25,0.25,0.25,0.26,
            0.28,0.3,0.34,0.37,0.41,0.45,0.48,0.52,0.56,0.61,0.66,0.72,0.79,0.85,0.93,0.99,1.06,1.12,1.17,1.22,1.26,1.3,1.36,1.44,
            1.56,1.75,2,2.34,2.77,3.28,3.86,4.46,5.05,5.56,5.94,6.13,6.1,5.84,5.36,4.69,3.9,3.07,2.27,1.46,0.57,0.05,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD6",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.07,0.09,0.11,0.13,0.14,0.15,0.16,0.16,0.16,0.15,0.15,0.15,0.15,0.17,
            0.19,0.21,0.24,0.27,0.3,0.32,0.34,0.36,0.37,0.38,0.4,0.42,0.45,0.48,0.53,0.58,0.63,0.68,0.72,0.73,0.71,0.65,0.57,0.48,0.4,
            0.39,0.5,0.76,1.22,1.92,2.82,3.91,5.07,6.23,7.23,7.96,8.31,8.24,7.72,6.86,5.72,4.47,3.23,2.11,1.24,0.54,0.15,0.02,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD7",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.03,0.07,0.1,0.12,0.15,0.17,0.18,0.19,0.2,0.2,0.2,0.2,0.2,0.22,0.24,0.26,0.3,
            0.34,0.39,0.43,0.47,0.51,0.55,0.59,0.63,0.67,0.71,0.76,0.81,0.87,0.94,1.02,1.11,1.2,1.29,1.38,1.45,1.5,1.54,1.56,1.59,1.66,
            1.78,2,2.33,2.8,3.4,4.08,4.81,5.49,6.05,6.4,6.48,6.24,5.73,4.96,4.04,3.06,2.1,1.06,0.2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD8",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.1,0.14,0.17,0.2,0.23,0.25,0.26,0.27,0.28,0.28,0.28,0.29,0.31,0.33,0.37,0.42,
            0.48,0.54,0.6,0.66,0.72,0.78,0.84,0.9,0.97,1.03,1.1,1.16,1.23,1.29,1.37,1.44,1.53,1.61,1.71,1.8,1.9,2.01,2.12,2.25,2.41,2.61,2.85,
            3.14,3.48,3.87,4.26,4.63,4.94,5.12,5.14,4.98,4.61,4.08,3.41,2.68,1.93,1.31,0.31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD9",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.1,0.15,0.19,0.22,0.25,0.28,0.29,0.3,0.31,0.32,0.32,0.33,0.35,0.37,0.41,0.45,
            0.5,0.55,0.6,0.65,0.7,0.76,0.82,0.88,0.95,1.03,1.1,1.17,1.23,1.3,1.37,1.45,1.54,1.66,1.82,2.01,2.26,2.55,2.89,3.26,3.66,4.06,4.43,
            4.75,4.99,5.13,5.14,5.03,4.78,4.42,3.94,3.39,2.79,2.19,1.6,1.11,0.71,0.15,0.02,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)), 
          stringsAsFactors=FALSE)
          
          ###
          colnames(TempTable) <- NULL
          write.csv(TempTable, file=GSDdata, row.names=FALSE)
          cat(paste("An editable GSD input template [", GSDdata, "] has been created at: ", paste(getwd(),"/", GSDdata,"\n",sep=""),sep=""))
             
      } else  {

          if (!is.null(outfile) && !is.character(outfile))  stop("Error: [outfile] should be NULL or a character!")

          ###
          if (!is.null(outfile))   xoxoxo <- pdfName <- rdName <- paste(outfile,"_batch",sep="") 

          ###
          if (is.character(GSDdata)) {

              if (length(GSDdata)!=1)  stop("Error: file name [GSDdata] should be of length one!")

              ###
              if (is.null(outfile))  { 

                  nxo <- nchar(GSDdata)
                  xoxoxo <- pdfName <- rdName <- paste(substr(GSDdata, start=1, stop=nxo-4),"_batch",sep="")

              } # end if.

              ###
              GSDdata <- read.csv(GSDdata, header=FALSE)
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!")

              ###
              gsl <- as.numeric(GSDdata[1,-1])
              GSD <- GSDdata[-1,-1,drop=FALSE]
          
              ###
              sampleName <- as.character(GSDdata[-1,1])

         } else {

              if (!is.matrix(GSDdata) && !is.data.frame(GSDdata))  stop("Error: [GSDdata] should be matrix or data.frame!")
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!")

              ###
              gsl <- as.numeric(GSDdata[1,])
              GSD <- GSDdata[-1,,drop=FALSE] 

              ###
              if (is.null(outfile)) {

                  mthCALL <- (as.character(match.call()))[2L]
                  xoxoxo <- pdfName <- rdName <- paste(mthCALL,"_batch",sep="")
 
              } # end if.

              ###
              sampleName <- paste("GSD",1:nrow(GSD),sep="")
          
          } # end if.

          ###
          if(any(diff(gsl)<=0))  stop("Error: grain-size levels in the first row of [GSDdata] are not of increasing order!")

          ###
          N <- nrow(GSD)
          GSDbatch <- vector(length=N, mode="list")
          cat("Unmixing of N=", N, " GSDs in a batch pattern is in progress, please wait, ...\n", sep="")
          pb <- txtProgressBar(min=1, max=N, initial=1, style=3)

          ###
          for (i in 1:N) {

              gsd <- as.numeric(GSD[i,])

              ###
              if (!is.null(model)) {
 
                  if (model %in% c("weibull", "lognormal")) {

                      res_ussGSD <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                        mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                        sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                        viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                        rmZero=rmZero, plot=FALSE), silent=TRUE)

                  } else if (model %in% c("weibull0", "lognormal0", "skewnormal0", "skewgnormal0")) {

                      res_ussGSD <- try(ussGSD0(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                        mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                        sigmaRange=sigmaRange, omegaRange=omegaRange, qRange=qRange, useIndex=useIndex, 
                        minfunc=minfunc, trim=trim, mrsl=mrsl, viewAutoInis=FALSE, viewLM=FALSE, 
                        viewFit=FALSE, viewPb=FALSE, outfile=NULL, rmZero=rmZero, plot=FALSE), 
                        silent=TRUE)

                  } # end if.

                  ###
                  if (inherits(res_ussGSD,what="try-error")==FALSE) { 

                      GSDbatch[[i]] <- res_ussGSD 

                  } else {

                      cat("Failed in unmixing the ", i, "-th sample!\n", sep="")
                      print(paste("<",sampleName[i],"> ",attr(res_ussGSD,"condition"),sep=""))

                  } # end if. 

              } else {

                  res_ussGSDwb <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model="weibull", 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                    viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                    rmZero=rmZero, plot=FALSE), silent=TRUE)

                  ###
                  res_ussGSDlg <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model="lognormal", 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                    viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                    rmZero=rmZero, plot=FALSE), silent=TRUE)                                    

                  ###
                  if (inherits(res_ussGSDwb,what="try-error")==FALSE && inherits(res_ussGSDlg,what="try-error")==FALSE) {

                      if (minfunc=="fom") {

                          minfuncVAL1 <- res_ussGSDwb$FOM
                          minfuncVAL2 <- res_ussGSDlg$FOM

                      } else if (minfunc=="rss") {

                          minfuncVAL1 <- res_ussGSDwb$RSS
                          minfuncVAL2 <- res_ussGSDlg$RSS

                      } # end if.

                      ###
                      if (minfuncVAL1<minfuncVAL2) {

                          GSDbatch[[i]] <- res_ussGSDwb

                      } else {

                          GSDbatch[[i]] <- res_ussGSDlg

                      } # end if. 

                  } else if (inherits(res_ussGSDwb,what="try-error")==FALSE) {

                      GSDbatch[[i]] <- res_ussGSDwb
                      cat("Failed in unmixing the ", i, "-th sample using Lognormal!\n", sep="")
                      print(paste("<",sampleName[i],"> ",attr(res_ussGSDlg,"condition"),sep=""))

                  } else if (inherits(res_ussGSDlg,what="try-error")==FALSE) {

                      GSDbatch[[i]] <- res_ussGSDlg
                      cat("Failed in unmixing the ", i, "-th sample using Weibull!\n", sep="")
                      print(paste("<",sampleName[i],"> ",attr(res_ussGSDwb,"condition"),sep=""))

                  } else {

                      cat("Failed in unmixing the ", i, "-th sample using Weibull and Lognormal!\n", sep="")
                      print(paste("<",sampleName[i],"> ",attr(res_ussGSDwb,"condition"),sep=""))
                      print(paste("<",sampleName[i],"> ",attr(res_ussGSDlg,"condition"),sep=""))

                  } # end if.

              } # end if. 

              ###
              setTxtProgressBar(pb,i) 

          } # end for.

          ###
          cat("\n")
          close(pb)
      
          ###
          names(GSDbatch) <- sampleName
          class(GSDbatch) <- "batchgsd"
          attr(GSDbatch, "xoxoxo") <- xoxoxo

          ###
          cat("Note: a new environment called [gsdbatch] is created to hold the unmixing results!\n",sep="")
          gsdbatch <<- new.env()
 
          ###
          assign("GSDbatch", GSDbatch, envir=gsdbatch)

          ###
          plot_ussGSDbatch(obj_batchgsd=GSDbatch, pdfName=pdfName, addvl=addvl, logxy=logxy, lwd=lwd, pch=pch, cex=cex) 
          save(GSDbatch, file=paste(rdName,".RData",sep=""))

          ###
          invisible(GSDbatch)

      } # end if.

  } # end function ussGSDbatch.
  ###==========================================================================================================================
  ###
  

  ###
  ###**************************************************************************************************************************
  ### Function plot_ussGSD() is used for ploting the unmixing result of a single-sample grain-size distribution.
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###obj_ussgsd: A S3 objective of class "ussgsd" generated using function ussGSD().
  ###
  ###sampleName: A character indicating the name of the grain-size distribution. 
  ### 
  ###     addvl: A logical value indicating whether vertical lines will be added to 
  ###            the plot showing the unmixed grain-size components.
  ###
  ###     logxy: A character indicating whether the x- or y-axis will be logged in the plot, 
  ###            one of "", "x", "y", "xy", or NULL. 
  ###
  ###       lwd: A numeric value giving the widths of lines in the plot.
  ###       
  ###       pch: An integer giving the type of symbols in the plot.
  ### 
  ###       cex: A numeric value giving the size of symbols in the plot.
  ###==========================================================================================================================
  ### The function returns a plot showing the unmixing result of a single-samples grain-size distribution.
  ### The plot consists of four sub-plots containing the following results:
  ###
  ### (1) The measured and fitted grain-size volume percentage, and the unmixed individual grain-size components.
  ###
  ### (2) A histogram showing the residuals calculated using the best set of optimized parameters.
  ### 
  ### (3) A histogram showing all figure-of-merit (FOM) values obtained from the "trial-and-error" protocol.
  ###
  ### (4) A histogram showing all residual sum square (RSS) values obtained from the "trial-and-error" protocol.
  ###==========================================================================================================================
  plot_ussGSD <- function(obj_ussgsd, sampleName="", addvl=TRUE, logxy="x", lwd=3, pch=21, cex=2) {

      if (inherits(obj_ussgsd,what="ussgsd")==FALSE)  stop("Error: [obj_ussgsd] should be an object of class 'ussgsd'!")

      ###
      origin_xd <- obj_ussgsd$gs.comp[,1]
      yd <- obj_ussgsd$gs.comp[,2]

      ###
      if (ncol(obj_ussgsd$gs.comp)>3) {

          yd_fit <- obj_ussgsd$gs.comp[,3]
          mat <- obj_ussgsd$gs.comp[,-(1:3),drop=FALSE]

      } # end if.
       
      ###
      reserveidx <- obj_ussgsd$reserveidx
      tmdxv <- obj_ussgsd$mdxv
      tmdyv <- obj_ussgsd$mdyv
      model <- obj_ussgsd$model
      ncomp <- obj_ussgsd$ncomp
      residuals <- obj_ussgsd$residuals
      fom0vec <- obj_ussgsd$FOMs
      rss0vec <- obj_ussgsd$RSSs
      rsl <- obj_ussgsd$rsl
      fom <- obj_ussgsd$FOM
      rss <- obj_ussgsd$RSS
      R2 <- obj_ussgsd$R2
      RSE <- obj_ussgsd$RSE
      ntry <- obj_ussgsd$ntry
      acpt <- obj_ussgsd$acpt
      normTEST <- obj_ussgsd$normTest

      ###
      if (is.null(logxy)) { 

          boolean1x <- obj_ussgsd$expGSlev

          ###
          mdyv <- obj_ussgsd$mdyv

          if (!is.null(mdyv)) {

              boolean2x <- 16*min(mdyv)<max(mdyv)

          } else {

              boolean2x <- FALSE

          } # end if. 
            
          ###      
          if (boolean1x && boolean2x) {

              logxy <- "xy"

          } else if (boolean1x) {
                    
              logxy <- "x"

          } else if (boolean2x) {

              logxy <- "y"

          } else {

              logxy <- ""

          } # end if. 

      } # end if.
      ###

      ###
      oldpar <- par("mar", "mfrow")
      on.exit(par(oldpar))

      ###
      myTK <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)
      myLB <- c("0.1","0.2","0.5","1","2","5","10","20","50","100","200","500","1000","2000","5000")
      ###
      lineCol <- c("deepskyblue", "orangered", "purple",   "violetred", "yellowgreen", "lightblue", 
                   "goldenrod", "forestgreen", "blue",  "plum", "tan", "violet", "grey50") 

      ###
      layout(matrix(c(1,1,1,2,1,1,1,3,1,1,1,4),ncol=3), respect=TRUE)
      par(mar=c(4,5,5,1.5)+0.1)

      ###
      if (length(reserveidx)>0) {

          plot(origin_xd[reserveidx], yd[reserveidx], type="n", cex.lab=2, log=logxy, cex.axis=2, mgp=c(3.2,1,0), 
               xlab="Grain size (um)", ylab="Volume percentage (%)", xaxt="n", main=sampleName, cex.main=2.2) 

      } else {

          plot(origin_xd, yd, type="n", cex.lab=2, log=logxy, cex.axis=2, mgp=c(3.2,1,0), 
               xlab="Grain size (um)", ylab="Volume percentage (%)", xaxt="n", main=sampleName, cex.main=2.2) 

      } # end if. 

      ###
      XaxisCentral <- median(axTicks(side=1))

      ###
      box(lwd=2,lwd=2)

      ###
      axis(side=1, at=myTK, labels=myLB, cex.axis=2)
 
      ###
      points(origin_xd[reserveidx], yd[reserveidx], type="p", lwd=2, pch=pch, col="grey50", cex=cex) 

      ###
      mypos <- mean(yd[reserveidx][origin_xd[reserveidx]<=XaxisCentral])<mean(yd[reserveidx][origin_xd[reserveidx]>XaxisCentral]) 

      ###
      if (ncol(obj_ussgsd$gs.comp)>3) {

          points(origin_xd, yd_fit, col="grey30", type="l", lwd=lwd)
  
          ###
          for (i in seq(ncomp))  {

              points(origin_xd, mat[,i], col=lineCol[i], type="l", lwd=lwd)  

              ###
              if (addvl==TRUE) {               
                 
                  points(x=c(0,tmdxv[i],tmdxv[i]),y=c(tmdyv[i],tmdyv[i],0), col=lineCol[i], type="l", lwd=lwd)

              } # end if.

           } # end for.     

           ###
           legend(ifelse(mypos,"topleft","topright"),
                  legend=c("Measured","Fitted", paste("Comp.",seq(ncomp),sep=""),
                  if(model=="weibull") "Weibull" else if (model=="lognormal") "Lognormal" else if 
                  (model=="weibull0") "Weibull0" else if (model=="lognormal0") "Lognormal0" else if
                  (model=="skewnormal0") "SN0" else if (model=="skewgnormal0") "SGN0", 
                  paste("mrsl:",round(min(rsl),2),sep=""), paste("RSS:",round(rss,2),sep=""), 
                  paste("FOM:",round(fom,2),sep="")), col=c("grey50", "grey30", lineCol[seq(ncomp)],rep("grey30",4)), 
                  pch=c(pch, rep(NA,ncomp+1),rep(4,4)), lty=c(NA, rep("solid",ncomp+1),NA,NA,NA,NA), 
                  yjust=2, ncol=1, cex=1.5, pt.cex=c(rep(2,ncomp+2),rep(1,4)), bty="o", lwd=3.0, pt.bg="white")

           ###
           par(mar=c(4,4,1,1)+0.1)

           ###
           gxhtvcz1 <- normTEST$statistic
           gxhtvcz2 <- normTEST$p.value

           ###
           hist(residuals, main=NULL, col="skyblue", xlab="Residuals", border="white", mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=1.5)
           box(lwd=2)
           
           ###
           legend("topleft",legend=c(paste("SampleSize: ",length(residuals),sep=""),
                  paste("Mean: ",round(mean(residuals),3),sep=""), paste("Sd: ",round(sd(residuals),3),sep=""),
                  paste("Statistic: ",round(gxhtvcz1,2),sep=""), paste("p-value: ",round(gxhtvcz2,3),sep="")), 
                  yjust=2, ncol=1, cex=1, pt.cex=1, bty="n", text.col="red")
          
           ###
           hist(fom0vec, main=NULL, col="purple",border="white", xlab="FOM (%)", mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=1.5)
           box(lwd=2)

           ###
           legend("topright",legend=c(paste("ntry: ",acpt,"/",ntry,sep=""), paste("Min: ",round(min(fom0vec),2),sep=""),
                  paste("Mean ",round(mean(fom0vec),2),sep=""), paste("Max: ",round(max(fom0vec),2),sep="")), 
                  yjust=2, ncol=1, cex=1, pt.cex=1, bty="n", text.col="red")

           ###
           hist(rss0vec, main=NULL, col="red", border="white", xlab="RSS", mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=1.5)
           box(lwd=2)

           ###
           legend("topright",legend=c(paste("ntry: ",acpt,"/",ntry,sep=""), paste("Min: ",round(min(rss0vec),2),sep=""),
                  paste("Mean ",round(mean(rss0vec),2),sep=""), paste("Max: ",round(max(rss0vec),2),sep="")), 
                  yjust=2, ncol=1, cex=1, pt.cex=1, bty="n", text.col="purple")

      } else {

           legend(ifelse(mypos,"topleft","topright"), legend=c("Measured", model), col=c("grey50",NA), 
                  pch=c(pch,NA), lty=c(NA,NA), yjust=2, ncol=1, cex=1.5, pt.cex=2, bty="o", lwd=3.0, pt.bg="white")

           ###
           par(mar=c(4,4,1,1)+0.1)
           plot(1,1,type="n", xlab="Residuals", ylab="", mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=1.5)
           box(lwd=2)

           ###
           plot(1,1,type="n", xlab="FOM (%)", ylab="", mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=1.5)
           box(lwd=2)

           ###
           plot(1,1,type="n", xlab="RSS", ylab="", mgp=c(2.5,1,0), cex.axis=1.5, cex.lab=1.5)
           box(lwd=2)

      } # end if.

  } # end function plot_ussGSD.
  ###==========================================================================================================================
  ###

  ###
  ###**************************************************************************************************************************
  ### Function ussGSD() is used for unmixing a single-sample grain-size distribution. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###         gsl: A numeric vector containing the grain-size levels (in unit um). 
  ###
  ###         gsd: A nummeric vector cotaining the volume percentages (in unit %) of the grain-size distribution. 
  ###
  ###       ncomp: An integer (from 1 to 13) indicating the number of components to be unmixed. 
  ###
  ###        auto: A logical value indicating whether performing a automatic grain-size unmixing, 
  ###              in this case the user needs not to specify the initials used for optimization.
  ###
  ###       model: A character indicating the model to be used, "weibull", or "lognormal".
  ###
  ###         mpd: A real number indicating the allowed lower limit of the volume percentage 
  ###              of a GSD that can be used to identify the location of component peaks. 
  ###
  ###         ctf: A numeric value (between -1 and 1) representing a critical threshold factor  
  ###              that controls the identification of peaks from the grain-size distribution.  
  ###              This argument can be used to prevent identifying a false peak characterized  
  ###              by a large positive second-order derivative. Specially, ctf=0 indicates
  ###              that peaks with second-order derivatives above zero will be precluded. 
  ###              The precluding effect will be increasingly suppressed as ctf increases 
  ###              from -1 to 1. 
  ###
  ###        ntry: An integer indicating the number of trials in a trial-and-error protocol. 
  ###
  ###         kkf: A numeric value controlling the range of values (sigma, or beta) from which    
  ###              random starting parameters will be generated during the "trial-and-error" 
  ###              protocol.
  ###
  ###   startPars: A two-column matrix containing starting parameters used for unmixing, 
  ###              the first column contains the modes of individual components, and the 
  ###              second column contains the alpha ("weibull") or sigma ("lognormal") 
  ###              value of individual components.            
  ###
  ###  alphaRange: A two-element vector indicating the lower and upper limits on alpha values of 
  ###              a Weibull distribution generated from a Uniform distribution.
  ###
  ###  sigmaRange: A two-element vector indicating the lower and upper limits on sigma values 
  ###              of a Lognormal distribution generated from a Uniform distribution.
  ###
  ###    useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###              used as the x-coordinate during the unmixing process.
  ###
  ###     minfunc: A character indicating the objective to be mimimized, either "fom" or "rss",
  ###              for the figure-of-merit value or the residual sum of squares, respectively.
  ###
  ###        trim: A logical value indicating whether the unmixing results will be trimed  
  ###              using a recursive optimization protocol.
  ###
  ###        mrsl: A numeric value indicating the minimum resolution between two adjacent     
  ###              components used to trim the unmixing results, which ranges from 0 to 2.
  ###              Peaks with a second-order derivative below zero and a resolution smaller
  ###              than mrsl will be deleted. 
  ###
  ###viewAutoInis: A logical value indicating whether the automatically generated initial modal 
  ###              sizes will be visualized before the unmixing process if auto=TRUE.
  ###
  ###      viewLM: A logical value indicating whether the optimization using the Levenberg-Marquardt
  ###              nonlinear least-squares algorithm will be output onto the screen.
  ###
  ###     viewFit: A logical value indicating whether the trial-and-error unmixing results will be visualized in a plot.
  ###
  ###      viewPb: A logical value indicating whether the progress bar will be visualized during calculation.
  ###
  ###     outfile: A character indicating the name of the returned CSV file containing the unmixing results 
  ###              of individual components. The CSV file will be saved in the current working directory.
  ###
  ###      rmZero: A logical value indicating whether zero volume percentages of GSD will be removed.
  ###
  ###        plot: A logical value indicating whether the unmixing result will be visualized in a plot.
  ###
  ###  sampleName: A character indicating the name of the grain-size distribution.
  ### 
  ###       addvl: A logical value indicating whether vertical lines will be added      
  ###              to the plot showing the unmixed grain-size components. 
  ###
  ###       logxy: A character indicating whether the x- or y-axis will be logged in the plot,  
  ###              one of "", "x", "y", "xy", or NULL.  
  ###
  ###         lwd: A numeric value giving the widths of lines in the plot.
  ###       
  ###         pch: An integer giving the type of symbols in the plot.
  ### 
  ###         cex: A numeric value giving the size of symbols in the plot.
  ###==========================================================================================================================
  ### The function returns an invisible list of S3 object of class "ussgsd" containing the following elements.
  ###
  ###      model: A character indicating the model used for unmixing.
  ###
  ###      ncomp: An integer indicating the number of components used for unmixing.
  ###
  ###   expGSlev: A logical value indicating whether the grain-size levels are measured in exponential scale.  
  ###
  ### reserveidx: An integer vector indicating the indices of grain-size data points used for plotting.
  ###  
  ###   fit.data: A two-column matrix containing the data used for unmixing.
  ###
  ###   fit.pars: A matrix containing the optimized parameters of the distribution function. 
  ###
  ###   normTest: A list containing the results of normality test of residuals.
  ###
  ###       FOMs: A numeric vector containing the figure-of-merit values generated during the trial-and-error protocol.
  ###
  ###       RSSs: A numeric vector containing the residual sum of squares generated during the trial-and-error protocol.
  ###
  ###       ntry: The number of trials in a trial-and-error protocol.
  ###
  ###       acpt: An integer indicating the number of accpeted trials.
  ###
  ###    gs.comp: A matrix containing the grain-size levels, the measured volume percentages, the predicted 
  ###             volume percentages, and the unmixed volume percentages of individual components.
  ###    
  ###       mdxv: A numeric vector containing the peak locations of individual components.
  ###
  ###       mdyv: A numeric vector containing the peak volume percentages (i.e., the heights or  
  ###             the maximum amplitude) of individual components.
  ###
  ###  residuals: A numeric vector containing the unmixing residuals.
  ###
  ###    gs.pars: A matrix containing the optimized parameters of the grain-size distrubtion, including
  ###             the proportions, means, meadians, modes, and standard deviations of individual components.
  ###
  ###         sp: A matrix containing the shape parameters of unmixed individual grain-size components,
  ###             [x1] the indice corresponding to the half volume percentage of a component (left side),
  ###             [x2] the indice corresponding to the half volume percentage of a component (right side),
  ###             [xm] the indice corresponding to the max volume percentage of a component,
  ###             [d1] the half-width at the left side of a component,
  ###             [d2] the half-width at the right side of a component,
  ###             [thw] the total half-width of a component, 
  ###             [sf] the symmetry factors of a component.
  ###
  ###        rsl: A numeric vector containing the resolutions between adjacent components.
  ###
  ###     deltaH: A numberic vector containing the difference in entropy before and after mixing between adjacent components.
  ###
  ###        FOM: A numeric value indicating the calculated figure-of-merit value.
  ###
  ###        RSS: A numeric value indicating the calculated residual sum of squares.
  ###   
  ###         R2: A numeric value indicating the calculated R2 statistic.
  ###
  ###        RSE: A numeric value indicating the calculated residual standard error. 
  ###==========================================================================================================================
  ### In addition, the function automatically generates a plot showing the 
  ### unmixing results produced by function plot_ussGSD() if plot=TRUE.
  ###==========================================================================================================================
  ussGSD <- function(gsl, gsd, ncomp=NULL, auto=FALSE, model="weibull", mpd=0, ctf=0.1, ntry=50, kkf=0.1, 
                     startPars=NULL, alphaRange=NULL, sigmaRange=NULL, useIndex=TRUE, minfunc="fom", trim=FALSE,
                     mrsl=0.6, viewAutoInis=TRUE, viewLM=FALSE, viewFit=FALSE, viewPb=TRUE, outfile=NULL, rmZero=TRUE,    
                     plot=TRUE, sampleName="", addvl=TRUE, logxy="x", lwd=3, pch=21, cex=2) {

        ###
        stopifnot(length(gsl)==length(gsd),
                  is.null(ncomp) || (length(ncomp)==1 && is.numeric(ncomp) && ncomp %in% 1:13),
                  is.logical(auto), length(auto)==1,
                  is.character(model), length(model)==1, model %in% c("weibull","lognormal"),
                  is.numeric(mpd), length(mpd)==1, mpd>=0, mpd<3,
                  is.numeric(ctf), length(ctf)==1, ctf>=-1, ctf<=1,
                  is.numeric(ntry), length(ntry)==1, ntry>=2, 
                  is.numeric(kkf), length(kkf)==1, kkf>0, kkf<1,
                  is.null(startPars) || (is.matrix(startPars) && ncol(startPars)==2), 
                  is.null(alphaRange) || (is.numeric(alphaRange) && length(alphaRange)==2),
                  is.null(sigmaRange) || (is.numeric(sigmaRange) && length(sigmaRange)==2),
                  is.logical(useIndex), length(useIndex)==1,
                  is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "rss"),
                  is.logical(trim), length(trim)==1,
                  is.numeric(mrsl), length(mrsl)==1, mrsl>=0, mrsl<=2,
                  is.logical(viewAutoInis), length(viewAutoInis)==1, 
                  is.logical(viewLM), length(viewLM)==1, 
                  is.logical(viewFit), length(viewFit)==1,
                  is.logical(viewPb), length(viewPb)==1,
                  is.null(outfile) || (is.character(outfile) && length(outfile)==1),
                  is.logical(rmZero), length(rmZero)==1,  
                  is.logical(plot), length(plot)==1,
                  is.character(sampleName), length(sampleName)==1, 
                  is.logical(addvl), length(addvl)==1, 
                  is.character(logxy), length(logxy)==1,
                  is.numeric(lwd), length(lwd)==1, 
                  is.numeric(pch), length(pch)==1,
                  is.numeric(cex), length(cex)==1)

        ###
        origin_xd <- as.numeric(gsl)
        yd <- as.numeric(gsd)

        ###
        if (any(!is.finite(origin_xd)))  stop("Error: argument [gsl] contains non-finite value!")
        if (any(!is.finite(yd)))  stop("Error: argument [gsd] contains non-finite value!")

        ### Check if the grain size levels are of log-scale.
        xd0xd0 <- round(diff(log(as.numeric(origin_xd))),1)
        YESORNO <- all(xd0xd0==min(xd0xd0)) || all(xd0xd0 %% min(xd0xd0)==0)
        expGSlev <- ifelse(YESORNO, TRUE, FALSE)

        ###
        nd <- length(yd)     

        ###
        if (useIndex==TRUE) {

            xd <- seq(origin_xd)

            ###
            if(expGSlev==TRUE) {
            
                if (!all(xd0xd0==xd0xd0[1]))  warning("Logged grain-size levels do not increase with equal steps!")

                ###
                xd1 <- log(origin_xd)      

            } else {
            
                xd2xd2 <- round(diff(as.numeric(origin_xd)),1)
                if (!all(xd2xd2==xd2xd2[1]))  warning("Grain-size levels do not increase with equal steps!")
 
                ###
                xd1 <- origin_xd 

            } # end if.

            ### Build a linear relationship between xd and xd1.
            ###------------------------------------------------
            pab <- as.numeric(lm(xd1~xd)$coefficients)

        } else {

            xd <- origin_xd
            
        } # end if.

        ###
        MINX <- min(xd)
        MAXX <- max(xd)

        ###
        ### Function used for remove data points.
        ###--------------------------------------------
        removeData <- function(yd, nd, v0) {

            if (yd[1]>v0 && yd[nd]>v0) {

                reserveidx <- seq(nd)

            } else {

                negzero11 <- which(diff(as.numeric(yd<=v0))<0)  
                negzero22 <- if (length(negzero11)>0)  {1:negzero11[1]}  else  {nd*1000}

                ###
                negzero33 <- which(diff(as.numeric(rev(yd)<=v0))<0)
                negzero44 <-  if (length(negzero33)>0)  {seq(from=nd, to=1, by=-1)[negzero33[1]]:nd}  else  {nd*1000+1}

                ###
                if (yd[1]>v0)  {

                    reserveidx <- (1:nd)[-negzero44]

                } else if (yd[nd]>v0) {

                    reserveidx <- (1:nd)[-negzero22]

               } else {

                    reserveidx <- (1:nd)[-c(negzero22,negzero44)]

               } # end if.

           } # end if.

           return(reserveidx)

        } # end function removeData.

        ###
        if (rmZero==TRUE) {

            reserveidx <- removeData(yd=yd, nd=nd, v0=0)

        } else {

            reserveidx <- seq(nd)

        } # end if.

        ###
        if (is.null(alphaRange))  { 

            alphaRange <- if (useIndex==TRUE) c(6,30) else c(2,10)

        } # end if.

        if (is.null(sigmaRange))  {

            sigmaRange <- if (useIndex==TRUE) c(0.01,0.2) else c(0.1,1)

        } # end if.

        ###
        myTK <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)
        myLB <- c("0.1","0.2","0.5","1","2","5","10","20","50","100","200","500","1000","2000","5000")

        ###
        boolean1 <- is.null(startPars)
        boolean2 <- is.null(ncomp) && auto==FALSE
        boolean3 <- auto==TRUE && viewAutoInis==TRUE
        boolean4 <- auto==FALSE
        boolean5 <- viewFit==TRUE

        ###
        if ( (boolean1 && (boolean2 || (boolean3 || boolean4))) || boolean5) {

            oldpar <- par("mar", "mfrow", "bg", "new")
            on.exit(par(oldpar))

        } # end if. 

        ###
        if (is.null(startPars)) {

            ###
            ### Find the modes of peaks from the 
            ### second-order derivatives of the GSD.
            ###-----------------------------------------------
            ###d2y <- tgcd::savgol(yd, drv=2, hwd=4, pod=2)
            d2y <- pracma::savgol(yd, fl=9, forder=2, dorder=2)

            ###
            reserveidx0 <- removeData(yd=yd, nd=nd, v0=mpd)

            ###
            peakIDX <- pracma::findpeaks(x=-d2y[reserveidx0], minpeakdistance=6, sortstr=TRUE)

            ###
            if (inherits(peakIDX, what="matrix")==FALSE)  peakIDX <- matrix(peakIDX, nrow=1L)

            ###
            CTFV <- ctf*min(-d2y[reserveidx0])

            ###
            peakIDX <- peakIDX[peakIDX[,1]>CTFV,,drop=FALSE]

            ###
            if(nrow(peakIDX)==0) {

                mat <- cbind(origin_xd, yd, rep(NA,nd))
                rownames(mat) <- NULL
                colnames(mat) <- c("GSlev","Volume","Fit.Volume")

                ###
                output <- list("model"=model, "ncomp"=ncomp, "expGSlev"=expGSlev, 
                               "reserveidx"=reserveidx, "gs.comp"=mat, "FOM"=Inf, "RSS"=Inf)
                class(output) <- "ussgsd"

                ###
                if (plot==TRUE) {

                    plot_ussGSD(output, sampleName=sampleName, addvl=addvl, lwd=lwd, cex=cex, pch=pch)

                } # end if.

                ###
                cat("Note: no component can be identified from the GSD!\n")

                ###
                return(invisible(output))

            } # end if.
            ###

            ###
            mdgs0 <- origin_xd[reserveidx0][peakIDX[,2]]
            mdgsy0 <- d2y[reserveidx0][peakIDX[,2]]

            ###
            if (is.null(ncomp) && auto==FALSE) {

                plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                     xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0), 
                     log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5,
                     main=paste("Identified ", length(mdgs0), " potential components",sep=""))

                ###
                box(lwd=2)

                ###
                axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
            
                ###         
                abline(h=-CTFV, col="black", lty=2, lwd=4)
                abline(v=mdgs0, col="grey70", lty=1, lwd=6)

                ###
                par("new"=TRUE)
                plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                     col="red", lwd=4, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                ###
                ncomp <- as.numeric(readline(paste("Please enter the number of components (i.e., [ncomp]): ")))
                if (!is.finite(ncomp))  stop("[ncomp] has not been provided!")

            } # end if. 

            ###
            if (auto==TRUE) {

                ### An automatic pattern.
                mdgs <- mdgs0
                mdgsy <- mdgsy0

                ###
                Nmdgs <- length(mdgs)

                ###
                if (is.null(ncomp)) {

                    ncomp <- Nmdgs

                } else {

                    if (ncomp>Nmdgs) {

                        cat("Note: [ncomp=",ncomp,"] exceeds the number of ", 
                            "automatically identified peaks(k=", Nmdgs,")!\n",sep="")

                        ###
                        ncomp <- Nmdgs

                    } # end if. 

                } # end if. 

                ###
                if (viewAutoInis==TRUE) {

                    plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                         xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0),  
                         log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5, 
                         main=paste("Number of components: ",ncomp, " out of ",Nmdgs, sep=""))

                    ###
                    box(lwd=2)

                    ###
                    axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                    axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                    points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
  
                    ###
                    abline(h=-CTFV, col="black", lty=2, lwd=4)
                    abline(v=mdgs, col="grey70", lty=1, lwd=6)
    
                    ###
                    abline(v=mdgs[1:ncomp], lwd=4, col="skyblue")
                    points(x=mdgs[1:ncomp], y=mdgsy[1:ncomp], pch=23, bg="red", col="red", cex=1.6) 
                    text(x=mdgs[1:ncomp], y=mdgsy[1:ncomp], labels=paste("C",1:ncomp,sep=""), 
                         col="red", cex=1.2)

                    ###
                    par("new"=TRUE)
                    plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                         col="red", lwd=3, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                } # end if. 

            } else {

                ### An interactive pattern.
                ###------------------------
                mdgs <- mdgsy <- vector(length=ncomp)

                ###
                for (i in 1:ncomp) {

                    plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                         xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0), 
                         log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5, 
                         main=paste("Number of components: ",i-1, " out of ",ncomp, sep=""))

                    ###
                    box(lwd=2)

                    ###
                    axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                    axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                    points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
           
                    ###
                    abline(h=-CTFV, col="black", lty=2, lwd=4)
                    abline(v=mdgs0, col="grey70", lty=1, lwd=6)

                    ###
                    par("new"=TRUE)
                    plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                         col="red", lwd=4, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                    ###
                    if (i>1)  { 

                        abline(v=mdgs[1:(i-1)], lwd=3, col="skyblue") 
                        points(x=mdgs[1:(i-1)], y=mdgsy[1:(i-1)], pch=23, bg="red", col="red", cex=1.5)
                        text(x=mdgs[1:(i-1)], y=mdgsy[1:(i-1)], labels=paste("C",1:(i-1),sep=""), col="red", cex=1.2)

                    } # end if.
     
                    ###
                    LLL <- try(locator(n=1), silent=TRUE)
                    if (inherits(LLL,what="try-error")==TRUE)  stop("Failed in parameter initialisation by clicking!")

                    ###
                    mdgs[i] <- LLL$x
                    mdgsy[i] <- LLL$y
            
                } # end for.

                ###
                apxvy <- approx(x=origin_xd[reserveidx0],y=-d2y[reserveidx0],xout=mdgs)$y
                odapxvy <- order(apxvy, decreasing=TRUE)

                ###
                mdgs <- mdgs[odapxvy]
                mdgsy <- mdgsy[odapxvy]

                ###
                plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                     xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0),
                     log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5, 
                     main=paste("Number of components: ",ncomp, " out of ",ncomp, sep=""))

                ###
                box(lwd=2)

                ###
                axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
  
                ###
                abline(h=-CTFV, col="black", lty=2, lwd=4)
                abline(v=mdgs, col="grey70", lty=1, lwd=6)
                abline(v=mdgs, lwd=3, col="skyblue") 

                ###
                par("new"=TRUE)
                plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                     col="red", lwd=4, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                ###
                points(x=mdgs, y=mdgsy, pch=23, bg="red", col="red", cex=1.5)
                text(x=mdgs, y=mdgsy, labels=paste("C",(1:ncomp),sep=""), col="red", cex=1.2)

            } # end if.

        } else {

            if (is.null(ncomp)) { 

                ncomp <- nrow(startPars)

            } else {

                if (ncomp!=nrow(startPars))  stop("Error: [ncomp] must be equal to the number of rows of [startPars]!")

            } # end if.

            ###
            mdgs <- startPars[,1,drop=TRUE]    

        } # end if.  

        ###
        ### Change the scale of mdgs.
        ###--------------------------
        if (useIndex==TRUE) {

            if (expGSlev==TRUE) {

                mdgs <- (log(mdgs)-pab[1])/pab[2]

            } else {

                mdgs <- (mdgs-pab[1])/pab[2]

            } # end if.

        } # end if.

        ###
        mdgs <- mdgs[1:ncomp]

        ###
        ### The function to be minimised.
        ###---------------------------------------------------
        minfn <- function(p, xd, yd, mdgs, model, minfunc)  {

             nd <- length(xd)

             ###
             hg <- rep(.Machine$double.xmax, nd)
             if (any(!is.finite(p)))  return(hg)

             ###
             ncomp <- length(p)/2

             ###
             if (model=="weibull") {
 
                 beta <- abs(p[1:ncomp]) + 1

                 ###
                 Mode <- abs(p[(ncomp+1):(2*ncomp)])

                 ###
                 idxm1 <- which(Mode<(1-kkf)*mdgs)
                 Mode[idxm1] <- (1-kkf)*mdgs[idxm1]

                 ###
                 idxm2 <- which(Mode>(1+kkf)*mdgs)
                 Mode[idxm2] <- (1+kkf)*mdgs[idxm2]

                 ###
                 theta <- Mode/((beta-1.0)/beta)^(1.0/beta) 

                 ###
                 mat <- matrix(nrow=nd, ncol=ncomp)

                 ###
                 for (i in 1:ncomp)  {

                     mat[,i] <- (beta[i]/theta[i])*((xd/theta[i])^(beta[i]-1.0))*exp(-(xd/theta[i])^beta[i])

                 } # end for. 

             } else if (model=="lognormal") {

                 beta <- abs(p[1:ncomp]) + 0.001

                 ###
                 Mode <- abs(p[(ncomp+1):(2*ncomp)])

                 ###
                 idxm1 <- which(Mode<(1-kkf)*mdgs)
                 Mode[idxm1] <- (1-kkf)*mdgs[idxm1]

                 ###
                 idxm2 <- which(Mode>(1+kkf)*mdgs)
                 Mode[idxm2] <- (1+kkf)*mdgs[idxm2]

                 ###
                 theta <- log(Mode)+beta^2

                 ###
                 mat <- matrix(nrow=nd, ncol=ncomp)

                 ###
                 for (i in 1:ncomp)  {

                     mat[,i] <- 1/sqrt(2*pi)/beta[i]/xd*exp(-0.5*((log(xd)-theta[i])/beta[i])^2)

                 } # end for. 

             } # end if.
             ### 

             ### Calculate abundances using a linear algebra method.
             ###----------------------------------------------------
             coef <- try(solve(a=t(mat)%*%mat, b=t(mat)%*%as.matrix(yd)), silent=TRUE)

             ###
             if (inherits(coef,what="try-error")==TRUE || any(!is.finite(coef)) || any(coef<0)) return(hg) 
 
             ###
             coef <- as.numeric(coef)

             ###
             mat <- matrix(nrow=nd, ncol=ncomp)
 
             for (i in 1:ncomp)   { 

                 if (model=="weibull") {

                     mat[,i] <- coef[i]*(beta[i]/theta[i])*((xd/theta[i])^(beta[i]-1.0))*exp(-(xd/theta[i])^beta[i])

                 } else if (model=="lognormal") {

                     mat[,i] <- coef[i]/sqrt(2*pi)/beta[i]/xd*exp(-0.5*((log(xd)-theta[i])/beta[i])^2)

                 } # end if. 

             } # end for.

             ###
             if (all(is.finite(mat))) {

                 yd_fit <- rowSums(mat)
              
                 ### The objective to be minimised.
                 ###------------------------------ 
                 if (minfunc=="fom") {

                     rsdlv <- sqrt(abs(yd-yd_fit)/sum(yd_fit))

                 } else if (minfunc=="rss") {

                     rsdlv <- yd-yd_fit

                 } # end if. 

                 ###
                 return(rsdlv)
 
             } else {
 
                 return(hg) 

             } # end if. 

       } # end function minfn.
       ###
       ###--------------------------------------------------------------------------------
       ### Calculate shape parameters of grain-size components.
       ###-----------------------------------------------------
       calShape <- function(y, x)  {
           
           ny <- length(y)
           maxloc <- which.max(y)
           hmaxval <- max(y)/2
           Tm <- x[maxloc]
           
           ###
           T1 <- suppressWarnings(try(approx(x=y[1L:maxloc], y=x[1L:maxloc], xout=hmaxval)$y, silent=TRUE))
           T2 <- suppressWarnings(try(approx(x=y[maxloc:ny], y=x[maxloc:ny], xout=hmaxval)$y, silent=TRUE))

           ###
           if (inherits(T1,what="try-error")==FALSE) { d1 <- Tm-T1 } else { T1 <- d1 <- NA } # end if.
           if (inherits(T2,what="try-error")==FALSE) { d2 <- T2-Tm } else { T2 <- d2 <- NA } # end if. 
 
           ###          
           thw <- T2-T1
           sf <- d2/thw

           ###         
           return(c("x1"=T1, "x2"=T2, "xm"=Tm, "d1"=d1, "d2"=d2, "thw"=thw, "sf"=sf))

       } # end function calShape.
       ###--------------------------------------------------------------------------------

       ###
       lineCol <- c("deepskyblue", "orangered", "purple",  "violetred",  "yellowgreen", "lightblue",  
                    "goldenrod", "forestgreen", "blue",  "plum", "tan", "violet", "grey50") 

       ###
       gdcp <- 0
       minobj <- .Machine$double.xmax
       acpt <- 0

       ###
       if (viewLM==FALSE && viewPb==TRUE) {

           cat(paste("ncomp=",ncomp, ", model=", model, ", ntry=", ntry, ".\n", sep="")) 
           cat("Unmixing of single-sample GSD (ussGSD) is in progress, please wait, ...\n", sep="")
            
           ###
           pb <- txtProgressBar(min=1, max=ntry, initial=1, style=3)
           
       } # end if. 

       ###
       fom0vec <- rss0vec <- c()
       ### Implement a trial-and-error protocol.
       ###---------------------------------------------------------------------------------
       for (i in seq(ntry))  {

           p0 <- vector(length=2*ncomp)

           ###
           if (is.null(startPars)) {

               if (model=="weibull") {
               
                   ### Simulate random alpha values.
                   ###-----------------------------
                   p0[1:ncomp] <- runif(n=ncomp, min=alphaRange[1], max=alphaRange[2])

               } else if (model=="lognormal") {

                   ### Simulate random sigma values.
                   ###------------------------------
                   p0[1:ncomp] <- exp(runif(n=ncomp, min=log(sigmaRange[1]), max=log(sigmaRange[2])))

               } # end if.

           } else {

               parsCOL2 <- startPars[,2,drop=TRUE]

               ###
               p0[1:ncomp] <- runif(n=ncomp, min=(1-kkf)*parsCOL2, max=(1+kkf)*parsCOL2)

           } # end if.
 
           ###
           p0[(ncomp+1):(2*ncomp)] <- runif(n=ncomp, min=(1-kkf)*mdgs[1:ncomp], max=(1+kkf)*mdgs[1:ncomp])
     
           ### Parameter optimisation using the Levenberg-Marquardt algorithm.
           ###----------------------------------------------------------------
           optLM <- try(minpack.lm::nls.lm(par=p0, lower=NULL, upper=NULL, fn=minfn, jac=NULL, 
                        control=nls.lm.control(maxiter=1024), xd, yd, mdgs, model, minfunc), silent=TRUE)
           ###
           if (viewLM==TRUE) print(optLM)

           ###
           if (inherits(optLM,what="try-error")==FALSE) {

               p <- optLM$par

               ###
               if (model=="weibull") {

                   beta0 <- abs(p[1:ncomp]) + 1

                   ###
                   Mode0 <- abs(p[(ncomp+1):(2*ncomp)])

                   ###
                   idxm1 <- which(Mode0<(1-kkf)*mdgs)
                   Mode0[idxm1] <- (1-kkf)*mdgs[idxm1]

                   ###
                   idxm2 <- which(Mode0>(1+kkf)*mdgs)
                   Mode0[idxm2] <- (1+kkf)*mdgs[idxm2]

                   ###
                   theta0 <- Mode0/((beta0-1.0)/beta0)^(1.0/beta0) 
    
                   ###
                   mat0 <- matrix(nrow=nd, ncol=ncomp)

                   ###
                   for (j in 1:ncomp)  {

                       mat0[,j] <- (beta0[j]/theta0[j])*((xd/theta0[j])^(beta0[j]-1.0))*exp(-(xd/theta0[j])^beta0[j])

                   } # end for.

                   ###
                   Mean0 <- theta0*gamma(1+1/beta0)

                   ###
                   Median0 <- theta0*(log(2))^(1/beta0)  

               } else if (model=="lognormal") {

                   beta0 <- abs(p[1:ncomp]) + 0.001

                   ###
                   Mode0 <- abs(p[(ncomp+1):(2*ncomp)])

                   ###
                   idxm1 <- which(Mode0<(1-kkf)*mdgs)
                   Mode0[idxm1] <- (1-kkf)*mdgs[idxm1]

                   ###
                   idxm2 <- which(Mode0>(1+kkf)*mdgs)
                   Mode0[idxm2] <- (1+kkf)*mdgs[idxm2]

                   ###
                   theta0 <- log(Mode0)+beta0^2
    
                   ###
                   mat0 <- matrix(nrow=nd, ncol=ncomp)

                   ###
                   for (j in 1:ncomp)  {

                       mat0[,j] <- 1/sqrt(2*pi)/beta0[j]/xd*exp(-0.5*((log(xd)-theta0[j])/beta0[j])^2)

                   } # end for.

                   ###
                   Mean0 <- exp(theta0+0.5*beta0^2)

                   ###
                   Median0 <- exp(theta0)

               } # end if.

               ###
               coef0 <- try(solve(a=t(mat0)%*%mat0, b=t(mat0)%*%as.matrix(yd)), silent=TRUE)

               ### Accept the unmixing results of the "trial-and-error" protocol 
               ### if the following conditions are satisfied.
               ###------------------------------------------------------------------
               OOKK <- inherits(coef0,what="try-error")==FALSE && all(is.finite(coef0)) && all(coef0>0) && 
                       all(is.finite(beta0)) && all(is.finite(theta0)) && 
                       all(is.finite(Mean0))  && all(Mean0>MINX & Mean0<MAXX) &&
                       all(is.finite(Median0)) && all(Median0>MINX & Median0<MAXX) && 
                       all(is.finite(Mode0)) && all(Mode0>MINX & Mode0<MAXX)  

               ###
               if (OOKK==TRUE)  {

                   acpt <- acpt + 1

                   ###
                   coef0 <- as.numeric(coef0)

                   ###
                   for (j in 1:ncomp)  mat0[,j] <- coef0[j]*mat0[,j]

                   ###
                   yd_fit0 <- rowSums(mat0)
                   fom0 <- sum(abs(yd-yd_fit0))/sum(yd_fit0)*100
                   rss0 <- sum((yd-yd_fit0)^2)

                   ###
                   fom0vec <- c(fom0vec, fom0)
                   rss0vec <- c(rss0vec, rss0)

                   ### Visualise the unmixing results of the "trial-and-error" protocol.
                   ###------------------------------------------------------------------ 
                   if (viewFit==TRUE) {

                       mat_vfp <- mat0
                       mat_vfp <- mat_vfp[,order(Mode0),drop=FALSE] 

                       ###
                       par(mar=c(6,6,5,5)+0.1)

                       ###
                       plot(origin_xd[reserveidx], yd[reserveidx], type="n", cex.lab=1.5, 
                            log="x", cex.axis=1.5, mgp=c(2.5, 1, 0), xlab="Grain size (um)", 
                            ylab="Volume percentage (%)", xaxt="n")

                       ###
                       box(lwd=2)

                       ###
                       XaxisCentral <- median(axTicks(side=1))

                       ###
                       axis(side=1, at=myTK, labels=myLB, cex.axis=1.5)

                       ###
                       points(origin_xd[reserveidx], yd[reserveidx], type="p", lwd=2, pch=21, col="gray50", cex=1.5) 

                       ###
                       points(origin_xd, rowSums(mat_vfp), col="grey30", type="l", lwd=3.0)

                       ###
                       for (j in 1:ncomp) {

                           points(origin_xd, mat_vfp[,j], col=lineCol[j], type="l", lwd=3)
                       
                       } # end for.
                          
                       ###
                       mypos <- mean(yd[reserveidx][origin_xd[reserveidx]<=XaxisCentral])<
                                mean(yd[reserveidx][origin_xd[reserveidx]>XaxisCentral]) 

                       ###
                       legend(ifelse(mypos,"topleft","topright"),
                              legend=c("Measured","Fitted", paste("Comp.",seq(ncomp),sep=""),
                              paste("Count:",acpt,sep=""), paste("RSS:",round(rss0,2),sep=""), 
                              paste("FOM:",round(fom0,2),sep="")), 
                              col=c("grey50", "grey30", lineCol[seq(ncomp)],NA,NA,NA), 
                              pch=c(21, rep(NA,ncomp+1),NA,NA,NA), 
                              lty=c(NA, rep("solid",ncomp+1),NA,NA,NA), 
                              yjust=2, ncol=1, cex=1, pt.cex=1, bty="o", 
                              lwd=3.0, pt.bg="white")

                   } # end if. 

                   ###
                   vxfght <- ifelse(minfunc=="fom", fom0, rss0)            

                   ### Update the optimal results.
                   if (vxfght<minobj) {

                       beta <- beta0               
                       theta <- theta0
                       coef <- coef0

                       ###
                       Mean <- Mean0
                       Median <- Median0
                       Mode <- Mode0

                       ###
                       mat <- mat0
                       minobj <- vxfght 
                      
                       ###
                       gdcp <- 1 
                       
                   } # end if.

               } # end if.

           } # end if.

           ###
           if (viewLM==FALSE && viewPb==TRUE)  setTxtProgressBar(pb, i)

       } # end for.
       ###---------------------------------------------------------------------------------
    
       ###
       if (viewLM==FALSE && viewPb==TRUE)  close(pb) 

       ###
       if (gdcp==0) {

           mat <- cbind(origin_xd, yd, rep(NA,nd))
           rownames(mat) <- NULL
           colnames(mat) <- c("GSlev","Volume","Fit.Volume")

           ###
           output <- list("model"=model, "ncomp"=ncomp, "expGSlev"=expGSlev, 
                          "reserveidx"=reserveidx, "gs.comp"=mat, "FOM"=Inf, "RSS"=Inf)
           class(output) <- "ussgsd"

           ###
           if (plot==TRUE) {

               plot_ussGSD(output, sampleName=sampleName, addvl=addvl, lwd=lwd, cex=cex, pch=pch)

           } # end if.

           ###
           return(invisible(output))

       } # end if.

       ###
       Np <- 3*ncomp
       yd_fit <- rowSums(mat)
       residuals <- yd-yd_fit
       fom <- sum(abs(yd-yd_fit))/sum(yd_fit)*100
       rss <- sum((yd-yd_fit)^2)
       sst <- sum((yd-mean(yd))^2)
       R2 <- 1 - rss/sst
       RSE <- sqrt(rss/(nd-Np))
       
       ###
       prop <- coef/sum(coef)*100

       ###
       gs.pars <- cbind(prop, Mean, Median, Mode)
       oidx <- order(gs.pars[,4])
       gs.pars <- gs.pars[oidx,,drop=FALSE]

       ###
       mat <- mat[,oidx,drop=FALSE]
       sp <- t(apply(mat, MARGIN=2L, calShape, seq(xd)))
       rownames(sp) <- paste("Comp.", seq(ncomp), sep="")

       ###
       if (ncomp>1) {

           resolvec <- name_resolvec <- vector(length=ncomp-1) 
           deltaHvec <- name_deltaHvec <- vector(length=ncomp-1)

           ###
           for (i in seq(ncomp-1)) {

               resolvec[i] <- (sp[i+1,3]-sp[i,3])/(sp[i,5]+sp[i+1,4])
               name_resolvec[i] <- paste("Comp.",i,".",i+1L,sep="")

               ###
               px12idx <- mat[,i]>0 & mat[,i+1]>0
               px1 <- mat[px12idx,i]/sum(mat[px12idx,i])
               px2 <- mat[px12idx,i+1]/sum(mat[px12idx,i+1])
               px12 <- px1 + px2
       
               ###
               deltaHvec[i] <- (-sum(px12*log2(px12))) - (-sum(px1*log2(px1))) - (-sum(px2*log2(px2)))
               name_deltaHvec[i] <- paste("Comp.",i,".",i+1L,sep="")

           } # end for.

           ###
           names(resolvec) <- name_resolvec
           names(deltaHvec) <- name_deltaHvec

       } else {

           resolvec <- deltaHvec <- Inf

       } # end if.
 
       ###
       fit.data <- cbind(xd,yd)
       rownames(fit.data) <- colnames(fit.data) <- NULL
       ###

       ###
       if (model=="weibull") {

           SDp <- sqrt( theta^2*(gamma(1+2/beta)-(gamma(1+1/beta))^2) )

           ###
           fit.pars <- cbind(coef, beta, theta)
           colnames(fit.pars) <- c("Coef", "Alpha", "Beta")

       } else {
 
           SDp <- sqrt( (exp(beta^2)-1)*exp(2*theta+beta^2) )

           ###
           fit.pars <- cbind(coef, theta, beta)
           colnames(fit.pars) <- c("Coef", "Mu", "Sigma")

       } # end if.

       ###
       fit.pars <- fit.pars[oidx,,drop=FALSE]
       rownames(fit.pars) <- paste("Comp.", seq(ncomp), sep="")

       ###
       SDp <- SDp[oidx]

       ###
       normTEST <- unclass(shapiro.test(x=residuals))

       ###
       if (useIndex==TRUE) {

           if (expGSlev==TRUE) {

               gs.pars[,2] <- exp(pab[1]+pab[2]*gs.pars[,2])
               gs.pars[,3] <- exp(pab[1]+pab[2]*gs.pars[,3])
               gs.pars[,4] <- exp(pab[1]+pab[2]*gs.pars[,4]) 

               ###
               SDp <-  gs.pars[,2]*pab[2]*SDp      

           } else {

               gs.pars[,2] <- pab[1]+pab[2]*gs.pars[,2]
               gs.pars[,3] <- pab[1]+pab[2]*gs.pars[,3]
               gs.pars[,4] <- pab[1]+pab[2]*gs.pars[,4]

               ###
               SDp <- pab[2]*SDp 

           } # end if.

       } # end if.

       ###
       gs.pars <- cbind(gs.pars, SDp)

       ###
       colnames(gs.pars) <- c("Proportion", "Mean", "Median", "Mode", "SD")   
       rownames(gs.pars) <- paste("Comp.", seq(ncomp), sep="")

       ###
       ### Implement a recursive unmixing protocol.
       ###-----------------------------------------------------------------
       if (trim==TRUE) {
           
           if (is.null(startPars)) {

               weakPK <- which(peakIDX[,1]<=0)

           } else {

               weakPK <- NA

           } # end if. 

           ###
           if ((all(is.finite(resolvec)) && min(resolvec)<mrsl) && ncomp>1 && length(weakPK)>=1) {

               ncomp <- ncomp - 1
               ###
               overLAPidx <- which.min(resolvec)

               ###
               vecxxx <- c(overLAPidx,overLAPidx+1)

               ###
               leaveCOMPidx <- vecxxx[which.min(gs.pars[vecxxx,1])]

               ###
               if (model=="weibull") {

                   startPars <- cbind(gs.pars[,4], fit.pars[,2])

               } else {

                   startPars <- cbind(gs.pars[,4], fit.pars[,3])

               } # end if.

               ###
               startPars <- startPars[-leaveCOMPidx,,drop=FALSE]

               ### 
               recursiveFIT <- ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                 mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                 sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                 viewAutoInis=viewAutoInis, viewLM=viewLM,  viewFit=viewFit, viewPb=viewPb,  
                 outfile=outfile, plot=plot, rmZero=rmZero, sampleName=sampleName, addvl=addvl,  
                 logxy=logxy, lwd=lwd, pch=pch, cex=cex)

               ###
               return(invisible(recursiveFIT))
      
           } # end if.

       } # end if.

       ###
       tmdxv <- gs.pars[,4]
       names(tmdxv) <- NULL
       tmdyv <- vector(length=ncomp)

       ###     
       for (i in seq(ncomp))  {

           tmdyv[i] <- approx(x=origin_xd, y=mat[,i], xout=tmdxv[i])$y

       } # end for. 

       ###
       mat <- cbind(origin_xd, yd, yd_fit, mat)
       rownames(mat) <- NULL
       colnames(mat) <- c("GSlev","Volume","Fit.Volume",paste("Comp.", seq(ncomp), sep=""))

       ###
       if (!is.null(outfile)) { 

           write.csv(mat, file=paste(outfile, ".csv", sep = ""))

       } # end if.

       ###
       output <- list("model"=model, "ncomp"=ncomp, "expGSlev"=expGSlev, "reserveidx"=reserveidx,
                      "fit.data"=fit.data, "fit.pars"=fit.pars, "normTest"=normTEST, "FOMs"=fom0vec, 
                      "RSSs"=rss0vec, "ntry"=ntry, "acpt"=acpt, "gs.comp"=mat, "mdxv"=tmdxv, "mdyv"=tmdyv, 
                      "residuals"=residuals, "gs.pars"=gs.pars, "sp"=sp, "rsl"=resolvec, "deltaH"=deltaHvec, 
                      "FOM"=fom, "RSS"=rss, "R2"=R2, "RSE"=RSE)

       ###
       class(output) <- "ussgsd"

       ###
       if (plot==TRUE) {

           plot_ussGSD(obj_ussgsd=output, sampleName=sampleName, addvl=addvl, logxy=logxy, lwd=lwd, cex=cex, pch=pch)

       } # end if.

       ###
       invisible(output)

  } # end function ussGSD.
  ###==========================================================================================================================
  ###


  ###
  ###*******************************************************************************************************************************
  ### Function ussGSD0() is used for unmixing single-sample grain-size distribution using transformed probability density functions. 
  ###===============================================================================================================================
  ### The function contains the following arguments.
  ###
  ###         gsl: A numeric vector containing the grain-size levels (in unit um). 
  ###
  ###         gsd: A nummeric vector cotaining the volume percentages (in unit %) of the grain-size distribution. 
  ###
  ###       ncomp: An integer (from 1 to 13) indicating the number of components to be unmixed. 
  ###
  ###        auto: A logical value indicating whether performing a automatic grain-size unmixing, 
  ###              in this case the user needs not to specify the initials used for optimization.
  ###
  ###       model: A character indicating the model to be used, "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###
  ###         mpd: A real number indicating the allowed lower limit of the volume percentage 
  ###              of a GSD that can be used to identify the location of component peaks.
  ###
  ###         ctf: A numeric value (between -1 and 1) representing a critical threshold factor  
  ###              that controls the identification of peaks from the grain-size distribution.  
  ###              This argument can be used to prevent identifying a false peak characterized  
  ###              by a large positive second-order derivative. Specially, ctf=0 indicates
  ###              that peaks with second-order derivatives above zero will be precluded. 
  ###              The precluding effect will be increasingly suppressed as ctf increases 
  ###              from -1 to 1. 
  ###
  ###        ntry: An integer indicating the number of trials in a trial-and-error protocol. 
  ###
  ###         kkf: A numeric value controlling the range of values (sigma, or beta) from which    
  ###              random starting parameters will be generated during the "trial-and-error" 
  ###              protocol.
  ###
  ###   startPars: if model="weibull0" or "lognormal0", [startPars] should be a three-column matrix 
  ###              containing starting parameters used for unmixing, the first column contains 
  ###              the modes of individual components, the second column contains the maximum 
  ###              volume percentages of individual components, and the thrid column contains
  ###              the alpha ("weibull0") or sigma ("lognormal0") value of individual components.
  ###
  ###              if model="skewnormal0" or "skewgnormal0", [startPars] should be a four-column 
  ###              matrix containing starting parameters used for unmixing, the first column is
  ###              the modes of individual components, the second column is the maximum volume
  ###              percentages of individual components,the thrid and fourth columns are the
  ###              alpha and omega values ("skewnormal0") or sigma and q values ("skewgnormal0") 
  ###              of individual components.              
  ###
  ###  alphaRange: A two-element vector indicating the lower and upper limits on alpha values of 
  ###              the Weibull ("weibull0") or Skew Normal ("skewnormal0") distributions generated 
  ###              from a Uniform distribution.
  ###
  ###  sigmaRange: A two-element vector indicating the lower and upper limits on sigma values 
  ###              of the Lognormal ("lognormal0") or Skewed Generalized Normal ("skewgnormal0")  
  ###              distributions generated from a Uniform distribution.
  ###
  ###  omegaRange: A two-element vector indicating the lower and upper limits on omega values of 
  ###              a Skew Normal ("skewnormal0") distribution generated from a Uniform distribution.
  ###
  ###      qRange: A two-element vector indicating the lower and upper limits on q values of a Skewed 
  ###              Generalized Normal ("skewgnormal0") distribution generated from a Uniform distribution.
  ###
  ###    useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###              used as the independent variable during the unmixing process.
  ###
  ###     minfunc: A character indicating the objective to be mimimized, either "fom" or "rss",
  ###              for the figure-of-merit value or the residual sum of squares, respectively.
  ###
  ###        trim: A logical value indicating whether the unmixing results will be trimed  
  ###              using a recursive optimization protocol.
  ###
  ###        mrsl: A numeric value indicating the minimum resolution between two adjacent     
  ###              components used to trim the unmixing results, which ranges from 0 to 2.
  ###              Peaks with a second-order derivative below zero and a resolution smaller
  ###              than mrsl will be deleted.
  ###
  ###viewAutoInis: A logical value indicating whether the automatically generated initial modal 
  ###              sizes will be visualized before the unmixing process if auto=TRUE.
  ###
  ###      viewLM: A logical value indicating whether the optimization using the Levenberg-Marquardt
  ###              nonlinear least-squares algorithm will be output onto the screen.
  ###
  ###     viewFit: A logical value indicating whether the trial-and-error unmixing results will be visualized in a plot.
  ###
  ###      viewPb: A logical value indicating whether the progress bar will be visualized during calculation.
  ###
  ###     outfile: A character indicating the name of the returned CSV file containing the unmixing results 
  ###              of individual components. The CSV file will be saved in the current working directory.
  ###
  ###      rmZero: A logical value indicating whether zero volume percentages of GSD will be removed.
  ###
  ###        plot: A logical value indicating whether the unmixing result will be visualized in a plot.
  ###
  ###  sampleName: A character indicating the name of the grain-size distribution.
  ### 
  ###       addvl: A logical value indicating whether vertical lines representing unmixed modal sizes      
  ###              will be added to the plot showing the unmixed grain-size components. 
  ###
  ###       logxy: A character indicating whether the x- or y-axis will be logged in the plot,  
  ###              one of "", "x", "y", "xy", or NULL.  
  ###
  ###         lwd: A numeric value giving the widths of lines in the plot.
  ###       
  ###         pch: An integer giving the type of symbols in the plot.
  ### 
  ###         cex: A numeric value giving the size of symbols in the plot.
  ###==========================================================================================================================
  ### The function returns an invisible list of S3 object of class "ussgsd" containing the following elements.
  ###
  ###      model: A character indicating the model used for unmixing.
  ###
  ###      ncomp: An integer indicating the number of components used for unmixing.
  ###
  ###   expGSlev: A logical value indicating whether the grain-size levels are measured in exponential scale.  
  ###
  ### reserveidx: An integer vector indicating the indices of grain-size data points used for plotting.
  ###  
  ###   fit.data: A two-column matrix containing the data used for unmixing.
  ###
  ###   fit.pars: A matrix containing the optimized parameters of the distribution function. 
  ###
  ###   normTest: A list containing the results of normality test of residuals.
  ###
  ###       FOMs: A numeric vector containing the figure-of-merit values generated during the trial-and-error protocol.
  ###
  ###       RSSs: A numeric vector containing the residual sum of squares generated during the trial-and-error protocol.
  ###
  ###       ntry: The number of trials in a trial-and-error protocol.
  ###
  ###       acpt: An integer indicating the number of accpeted trials.
  ###
  ###    gs.comp: A matrix containing the grain-size levels, the measured volume percentages, the predicted 
  ###             volume percentages, and the unmixed volume percentages of individual components.
  ###    
  ###       mdxv: A numeric vector containing the peak locations of individual components.
  ###
  ###       mdyv: A numeric vector containing the peak volume percentages (i.e., the heights or  
  ###             the maximum amplitude) of individual components.
  ###
  ###  residuals: A numeric vector containing the unmixing residuals.
  ###
  ###    gs.pars: A matrix containing the optimized parameters of the grain-size distrubtion, including
  ###             the proportions, means, meadians, modes, and standard deviations of individual components.
  ###
  ###         sp: A matrix containing the shape parameters of unmixed individual grain-size components,
  ###             [x1] the indice corresponding to the half volume percentage of a component (left side),
  ###             [x2] the indice corresponding to the half volume percentage of a component (right side),
  ###             [xm] the indice corresponding to the max volume percentage of a component,
  ###             [d1] the half-width at the left side of a component,
  ###             [d2] the half-width at the right side of a component,
  ###             [thw] the total half-width of a component, 
  ###             [sf] the symmetry factors of a component.
  ###
  ###        rsl: A numeric vector containing the resolutions between adjacent components.
  ###
  ###     deltaH: A numberic vector containing the difference in entropy before and after mixing between two adjacent components.
  ###
  ###        FOM: A numeric value indicating the optimal figure-of-merit value.
  ###
  ###        RSS: A numeric value indicating the optimal residual sum of squares.
  ###   
  ###         R2: A numeric value indicating the optimal R2 statistic.
  ###
  ###        RSE: A numeric value indicating the optimal residual standard error. 
  ###==========================================================================================================================
  ### In addition, the function automatically generates a plot showing the 
  ### unmixing results produced by function plot_ussGSD() if plot=TRUE.
  ###==========================================================================================================================
  ussGSD0 <- function(gsl, gsd, ncomp=NULL, auto=FALSE, model="weibull0", mpd=0, ctf=0.1, ntry=50, kkf=0.1, 
                      startPars=NULL, alphaRange=NULL, sigmaRange=NULL, omegaRange=NULL, qRange=NULL, useIndex=TRUE,
                      minfunc="fom", trim=FALSE, mrsl=0.6, viewAutoInis=TRUE, viewLM=FALSE, viewFit=FALSE, viewPb=TRUE,  
                      outfile=NULL, rmZero=TRUE, plot=TRUE, sampleName="", addvl=TRUE, logxy="x", lwd=3, pch=21, cex=2) {

        ###
        stopifnot(length(gsl)==length(gsd),
                  is.null(ncomp) || (length(ncomp)==1 && is.numeric(ncomp) && ncomp %in% 1:13),
                  is.logical(auto), length(auto)==1,
                  is.character(model), length(model)==1, model %in% c("weibull0","lognormal0","skewnormal0","skewgnormal0"),
                  is.numeric(mpd), length(mpd)==1, mpd>=0, mpd<3,
                  is.numeric(ctf), length(ctf)==1, ctf>=-1, ctf<=1,
                  is.numeric(ntry), length(ntry)==1, ntry>=2, 
                  is.numeric(kkf), length(kkf)==1, kkf>0, kkf<1,
                  is.null(startPars)  || (is.matrix(startPars) && ncol(startPars) %in% c(3,4)), 
                  is.null(alphaRange) || (is.numeric(alphaRange) && length(alphaRange)==2),
                  is.null(sigmaRange) || (is.numeric(sigmaRange) && length(sigmaRange)==2),
                  is.null(omegaRange) || (is.numeric(omegaRange) && length(omegaRange)==2),
                  is.null(qRange) || (is.numeric(qRange) && length(qRange)==2),
                  is.logical(useIndex), length(useIndex)==1,
                  is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "rss"),
                  is.logical(trim), length(trim)==1,
                  is.numeric(mrsl), length(mrsl)==1, mrsl>=0, mrsl<=2, 
                  is.logical(viewAutoInis), length(viewAutoInis)==1, 
                  is.logical(viewLM), length(viewLM)==1, 
                  is.logical(viewFit), length(viewFit)==1,
                  is.logical(viewPb), length(viewPb)==1,
                  is.null(outfile) || (is.character(outfile) && length(outfile)==1),
                  is.logical(rmZero), length(rmZero)==1,  
                  is.logical(plot), length(plot)==1,
                  is.character(sampleName), length(sampleName)==1, 
                  is.logical(addvl), length(addvl)==1, 
                  is.character(logxy), length(logxy)==1,
                  is.numeric(lwd), length(lwd)==1, 
                  is.numeric(pch), length(pch)==1,
                  is.numeric(cex), length(cex)==1)

        ###
        origin_xd <- as.numeric(gsl)
        yd <- as.numeric(gsd)

        ###
        if (any(!is.finite(origin_xd)))  stop("Error: argument [gsl] contains non-finite value!")
        if (any(!is.finite(yd)))  stop("Error: argument [gsd] contains non-finite value!")

        ###
        if (model=="skewnormal0" && useIndex==FALSE)
        cat("Warning: apply [skewnormal0] with [useIndex=FALSE] may yield unreasonable result!\n\n")
        if (model=="skewgnormal0" && useIndex==FALSE)
        cat("Warning: apply [skewgnormal0] with [useIndex=FALSE] may yield unreasonable result!\n\n")

        ###
        ### Check if the grain size levels are of log-scale.
        xd0xd0 <- round(diff(log(as.numeric(origin_xd))),1)
        YESORNO <- all(xd0xd0==min(xd0xd0)) || all(xd0xd0 %% min(xd0xd0)==0)
        expGSlev <- ifelse(YESORNO, TRUE, FALSE)

        ###
        nd <- length(yd)     

        ###
        if (useIndex==TRUE) {

            xd <- seq(origin_xd)

            ###
            if(expGSlev==TRUE) {
            
                if (!all(xd0xd0==xd0xd0[1]))  warning("Logged grain-size levels do not increase with equal steps!")

                ###
                xd1 <- log(origin_xd)  

            } else {
            
                xd2xd2 <- round(diff(as.numeric(origin_xd)),1)
                if (!all(xd2xd2==xd2xd2[1]))  warning("Grain-size levels do not increase with equal steps!")
 
                ###
                xd1 <- origin_xd 

            } # end if.

            ### Build a linear relationship between xd and xd1.
            ###------------------------------------------------
            pab <- as.numeric(lm(xd1~xd)$coefficients)

        } else {

            xd <- origin_xd

        } # end if.

        ###
        MINX <- min(xd)
        MAXX <- max(xd)

        ###
        ### Function used for remove data points.
        ###--------------------------------------------
        removeData <- function(yd, nd, v0) {

            if (yd[1]>v0 && yd[nd]>v0) {

                reserveidx <- seq(nd)

            } else {

                negzero11 <- which(diff(as.numeric(yd<=v0))<0)  
                negzero22 <- if (length(negzero11)>0)  {1:negzero11[1]}  else  {nd*1000}

                ###
                negzero33 <- which(diff(as.numeric(rev(yd)<=v0))<0)
                negzero44 <-  if (length(negzero33)>0)  {seq(from=nd, to=1, by=-1)[negzero33[1]]:nd}  else  {nd*1000+1}

                ###
                if (yd[1]>v0)  {

                    reserveidx <- (1:nd)[-negzero44]

                } else if (yd[nd]>v0) {

                    reserveidx <- (1:nd)[-negzero22]

               } else {

                    reserveidx <- (1:nd)[-c(negzero22,negzero44)]

               } # end if.

           } # end if.

           return(reserveidx)

        } # end function removeData.

        ###
        if (rmZero==TRUE) {

            reserveidx <- removeData(yd=yd, nd=nd, v0=0)

        } else {

            reserveidx <- seq(nd)

        } # end if.

        ###
        if (is.null(alphaRange))  { 

            if (model=="weibull0") {

                alphaRange <- if (useIndex==TRUE) c(6,30) else c(2,10)

            } else if (model=="skewnormal0") {

                alphaRange <- if (useIndex==TRUE) c(-5,5) else c(2,10)

            } # end if.

        } # end if.

        ###
        if (is.null(sigmaRange))  {

            if (model=="lognormal0") {

                sigmaRange <- if (useIndex==TRUE) c(0.01,0.2) else c(0.1,1)

            } else if (model=="skewgnormal0") {

                sigmaRange <- if (useIndex==TRUE) c(3,9) else c(10,30)

            } # end if. 

        } # end if.

        ###
        if (is.null(omegaRange)) {

            omegaRange <- if (useIndex==TRUE) c(5,10) else c(5,20)

        } # end if.

        ###
        if (is.null(qRange)) {

            qRange <- if (useIndex==TRUE) c(0.6,0.9) else c(0.5,0.9) 

        } # end if.

        ###
        myTK <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)
        myLB <- c("0.1","0.2","0.5","1","2","5","10","20","50","100","200","500","1000","2000","5000")

        ###
        boolean1 <- is.null(startPars)
        boolean2 <- is.null(ncomp) && auto==FALSE
        boolean3 <- auto==TRUE && viewAutoInis==TRUE
        boolean4 <- auto==FALSE
        boolean5 <- viewFit==TRUE

        ###
        if ( (boolean1 && (boolean2 || (boolean3 || boolean4))) || boolean5) {

            oldpar <- par("mar", "mfrow", "bg", "new")
            on.exit(par(oldpar))

        } # end if. 

        ###
        if (is.null(startPars)) {

            ###
            ### Find the modes of peaks from the 
            ### second-order derivatives of the GSD.
            ###-----------------------------------------------
            ###d2y <- tgcd::savgol(yd, drv=2, hwd=4, pod=2)
            d2y <- pracma::savgol(yd, fl=9, forder=2, dorder=2)

            ###
            reserveidx0 <- removeData(yd=yd, nd=nd, v0=mpd)

            ###
            peakIDX <- pracma::findpeaks(x=-d2y[reserveidx0], minpeakdistance=6, sortstr=TRUE)

            ###
            if (inherits(peakIDX, what="matrix")==FALSE)  peakIDX <- matrix(peakIDX, nrow=1L)

            ###
            CTFV <- ctf*min(-d2y[reserveidx0])

            ###
            peakIDX <- peakIDX[peakIDX[,1]>CTFV,,drop=FALSE]

            ###
            if(nrow(peakIDX)==0) {

                mat <- cbind(origin_xd, yd, rep(NA,nd))
                rownames(mat) <- NULL
                colnames(mat) <- c("GSlev","Volume","Fit.Volume")

                ###
                output <- list("model"=model, "ncomp"=ncomp, "expGSlev"=expGSlev, 
                               "reserveidx"=reserveidx, "gs.comp"=mat, "FOM"=Inf, "RSS"=Inf)
                class(output) <- "ussgsd"

                ###
                if (plot==TRUE) {

                    plot_ussGSD(output, sampleName=sampleName, addvl=addvl, lwd=lwd, cex=cex, pch=pch)

                } # end if.

                ###
                cat("Note: no component can be identified from the GSD!\n")

                ###
                return(invisible(output))

            } # end if.
            ###

            ###
            mdgs0 <- origin_xd[reserveidx0][peakIDX[,2]]
            mdgsy0 <- d2y[reserveidx0][peakIDX[,2]]

            ###
            if (is.null(ncomp) && auto==FALSE) {

                plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                     xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0),
                     log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5,
                     main=paste("Identified ", length(mdgs0), " potential components",sep=""))

                ###
                box(lwd=2)

                ###
                axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
            
                ###         
                abline(h=-CTFV, col="black", lty=2, lwd=4)
                abline(v=mdgs0, col="grey70", lty=1, lwd=6)

                ###
                par("new"=TRUE)
                plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                     col="red", lwd=4, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                ###
                ncomp <- as.numeric(readline(paste("Please enter the number of components (i.e., [ncomp]): ")))

            } # end if. 

            ###
            if (auto==TRUE) {

                ### An automatic pattern.
                mdgs <- mdgs0
                mdgsy <- mdgsy0

                ###
                Nmdgs <- length(mdgs)

                ###
                if (is.null(ncomp)) {

                    ncomp <- Nmdgs

                } else {

                    if (ncomp>Nmdgs) {

                        cat("Note: [ncomp=",ncomp,"] exceeds the number of ", 
                            "automatically identified peaks(k=", Nmdgs,")!\n",sep="")

                        ###
                        ncomp <- Nmdgs

                    } # end if. 

                } # end if. 

                ###
                if (viewAutoInis==TRUE) {

                    plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                         xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0),  
                         log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5, 
                         main=paste("Number of components: ",ncomp, " out of ",Nmdgs, sep=""))

                    ###
                    box(lwd=2)

                    ###
                    axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                    axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                    points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
  
                    ###
                    abline(h=-CTFV, col="black", lty=2, lwd=4)
                    abline(v=mdgs, col="grey70", lty=1, lwd=6)
    
                    ###
                    abline(v=mdgs[1:ncomp], lwd=4, col="skyblue")
                    points(x=mdgs[1:ncomp], y=mdgsy[1:ncomp], pch=23, bg="red", col="red", cex=1.6) 
                    text(x=mdgs[1:ncomp], y=mdgsy[1:ncomp], labels=paste("C",1:ncomp,sep=""), 
                         col="red", cex=1.2)

                    ###
                    par("new"=TRUE)
                    plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                         col="red", lwd=3, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                } # end if. 

            } else {

                ### An interactive pattern.
                ###------------------------
                mdgs <- mdgsy <- vector(length=ncomp)

                ###
                for (i in 1:ncomp) {

                    plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                         xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0),
                         log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5,
                         main=paste("Number of components: ",i-1, " out of ",ncomp, sep=""))

                    ###
                    box(lwd=2)

                    ###
                    axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                    axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                    points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
           
                    ###
                    abline(h=-CTFV, col="black", lty=2, lwd=4)
                    abline(v=mdgs0, col="grey70", lty=1, lwd=6)

                    ###
                    par("new"=TRUE)
                    plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                         col="red", lwd=4, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                    ###
                    if (i>1)  { 

                        abline(v=mdgs[1:(i-1)], lwd=3, col="skyblue") 
                        points(x=mdgs[1:(i-1)], y=mdgsy[1:(i-1)], pch=23, bg="red", col="red", cex=1.5)
                        text(x=mdgs[1:(i-1)], y=mdgsy[1:(i-1)], labels=paste("C",1:(i-1),sep=""), col="red", cex=1.2)

                    } # end if.
     
                    ###
                    LLL <- try(locator(n=1), silent=TRUE)
                    if (inherits(LLL,what="try-error")==TRUE)  stop("Failed in parameter initialisation by clicking!")

                    ###
                    mdgs[i] <- LLL$x
                    mdgsy[i] <- LLL$y
            
                } # end for.

                ###
                apxvy <- approx(x=origin_xd[reserveidx0],y=-d2y[reserveidx0],xout=mdgs)$y
                odapxvy <- order(apxvy, decreasing=TRUE)

                ###
                mdgs <- mdgs[odapxvy]
                mdgsy <- mdgsy[odapxvy]

                ###
                plot(origin_xd[reserveidx0], d2y[reserveidx0], type="n", 
                     xlab="Grain size (um)", ylab="Second-order derivative of GSD", mgp=c(2.5,1,0),
                     log="x", cex.axis=1.2, cex.lab=1.5, xaxt="n", yaxt="n", cex.main=1.5, 
                     main=paste("Number of components: ",ncomp, " out of ",ncomp, sep=""))

                ###
                box(lwd=2)

                ###
                axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
                axis(side=2, col="purple", lwd=4, cex.axis=1.2)
                points(origin_xd[reserveidx0], d2y[reserveidx0], type="l", col="purple", lwd=4)
  
                ###
                abline(h=-CTFV, col="black", lty=2, lwd=4)
                abline(v=mdgs, col="grey70", lty=1, lwd=6)
                abline(v=mdgs, lwd=3, col="skyblue") 

                ###
                par("new"=TRUE)
                plot(origin_xd[reserveidx0], yd[reserveidx0], log="x", type="l", 
                     col="red", lwd=4, lty=3, xlab="", ylab="", xaxt="n", yaxt="n")

                ###
                points(x=mdgs, y=mdgsy, pch=23, bg="red", col="red", cex=1.5)
                text(x=mdgs, y=mdgsy, labels=paste("C",(1:ncomp),sep=""), col="red", cex=1.2)

            } # end if.

            ###
            mdvp <- approx(x=origin_xd[reserveidx0],y=yd[reserveidx0],xout=mdgs)$y

        } else {

            if (is.null(ncomp)) { 

                ncomp <- nrow(startPars)

            } else {

                if (ncomp!=nrow(startPars))  stop("Error: [ncomp] must be equal to the number of rows of [startPars]!")

            } # end if.

            ###
            mdgs <- startPars[,1,drop=TRUE] 
            mdvp <- startPars[,2,drop=TRUE]  

        } # end if.  

        ###
        ### Change the scale of mdgs.
        ###--------------------------
        if (useIndex==TRUE) {

            if (expGSlev==TRUE) {

                mdgs <- (log(mdgs)-pab[1])/pab[2]

            } else {

                mdgs <- (mdgs-pab[1])/pab[2]

            } # end if.

        } # end if.

        ###
        mdgs <- mdgs[1:ncomp]
        mdgsORD <- order(mdgs)

        ### The error function.
        erf <- function(x) { 2 * pnorm(x * sqrt(2)) - 1 }

        ###
        ### The function to be minimised.
        ###---------------------------------------------------
        minfn <- function(p, xd, yd, mdgs, model, minfunc)  {

             nd <- length(xd)

             ###
             hg <- rep(.Machine$double.xmax, nd)
             if (any(!is.finite(p)))  return(hg)

             ###
             if (model=="weibull0") {

                 ncomp <- length(p)/3
 
                 ###
                 xm <- abs(p[1:ncomp])

                 ###
                 vm <- abs(p[(ncomp+1):(2*ncomp)])

                 ###
                 alpha <- abs(p[(2*ncomp+1):(3*ncomp)])+1

                 ###
                 idxm1 <- which(xm<(1-kkf)*mdgs)
                 xm[idxm1] <- (1-kkf)*mdgs[idxm1]

                 ###
                 idxm2 <- which(xm>(1+kkf)*mdgs)    
                 xm[idxm2] <- (1+kkf)*mdgs[idxm2]

                 ###
                 mat <- matrix(nrow=nd, ncol=ncomp)

                 ###
                 for (i in 1:ncomp)  {

                     mat[,i] <- vm[i]*(xd/xm[i])^(alpha[i]-1)*exp((alpha[i]-1)/alpha[i]*(1-(xd/xm[i])^alpha[i]))

                 } # end for. 

             } else if (model=="lognormal0") {

                 ncomp <- length(p)/3

                 ###
                 xm <- abs(p[1:ncomp])

                 ###
                 vm <- abs(p[(ncomp+1):(2*ncomp)])

                 ###
                 sigma <- abs(p[(2*ncomp+1):(3*ncomp)])+0.001

                 ###
                 idxm1 <- which(xm<(1-kkf)*mdgs)
                 xm[idxm1] <- (1-kkf)*mdgs[idxm1]

                 ###
                 idxm2 <- which(xm>(1+kkf)*mdgs)    
                 xm[idxm2] <- (1+kkf)*mdgs[idxm2]

                 ###
                 mat <- matrix(nrow=nd, ncol=ncomp)

                 ###
                 for (i in 1:ncomp)  {

                     mat[,i] <- xm[i]/xd*vm[i]*exp(-log(xm[i]/xd)*(log(xm[i]/xd)+2*(sigma[i])^2)/2/(sigma[i])^2)

                 } # end for. 

             } else if (model=="skewnormal0") {

                 ncomp <- length(p)/4

                 ###
                 xm <- abs(p[1:ncomp])

                 ###
                 vm <- abs(p[(ncomp+1):(2*ncomp)])

                 ###
                 alpha <- p[(2*ncomp+1):(3*ncomp)]

                 ###
                 omega <- abs(p[(3*ncomp+1):(4*ncomp)])+0.001

                 ###
                 idxm1 <- which(xm<(1-kkf)*mdgs)
                 xm[idxm1] <- (1-kkf)*mdgs[idxm1]

                 ###
                 idxm2 <- which(xm>(1+kkf)*mdgs)    
                 xm[idxm2] <- (1+kkf)*mdgs[idxm2]

                 ###
                 mat <- matrix(nrow=nd, ncol=ncomp)

                 ###
                 for (i in 1:ncomp)  {

                     delta <- alpha[i]/sqrt(1+(alpha[i])^2)

                     ###
                     v1 <- sqrt(1-2*delta^2/pi)*(4-pi)/4*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^1.5

                     ###
                     v2 <- sign(alpha[i])/2*exp(-2*pi/abs(alpha[i]))

                     ###
                     D <- omega[i]*(sqrt(2/pi)*delta-v1-v2)

                     ###
                     mat[,i] <- vm[i]*exp(-(2*D+xd-xm[i])*(xd-xm[i])/2/(omega[i])^2)*
                     (1+erf(alpha[i]*(xd-xm[i]+D)/sqrt(2)/omega[i]))/(1+erf(alpha[i]*D/sqrt(2)/omega[i]))

                 } # end for. 
                 
             } else if (model=="skewgnormal0") {

                 ncomp <- length(p)/4

                 ###
                 xm <- abs(p[1:ncomp])

                 ###
                 vm <- abs(p[(ncomp+1):(2*ncomp)])

                 ###
                 sigma <- abs(p[(2*ncomp+1):(3*ncomp)])+0.001

                 ###
                 qv <- abs(p[(3*ncomp+1):(4*ncomp)])+0.001

                 ###
                 pv <- 2 + 6*(1-qv)^5

                 ###
                 idxm1 <- which(xm<(1-kkf)*mdgs)
                 xm[idxm1] <- (1-kkf)*mdgs[idxm1]

                 ###
                 idxm2 <- which(xm>(1+kkf)*mdgs)    
                 xm[idxm2] <- (1+kkf)*mdgs[idxm2]

                 ###
                 mat <- matrix(nrow=nd, ncol=ncomp)

                 ###
                 ###x0 <- seq(from=0,to=MAXX/3,by=MAXX/3/300)
                 x0 <- seq(from=0,to=MAXX/4,by=MAXX/4/200)
                 
                 ###
                 for (i in 1:ncomp)  {

                     yf1 <- exp(x0/sigma[i]*qv[i])
                     yf2 <- exp(x0/sigma[i]/qv[i])

                     ###
                     y0 <- abs(yf1*qv[i]+yf2/qv[i])/(yf1+yf2)*exp(-0.5*(abs(log((yf1+yf2)/2)))^pv[i])
       
                     ###
                     if (any(!is.finite(y0)))  return(hg)  

                     ###plot(x0,y0,main=paste(sigma[i],",",qv[i],sep=""))               
                     L <- x0[which.max(y0)]

                     ###
                     pf1 <- exp((xd-xm[i]+L)/sigma[i]*qv[i])
                     pf2 <- exp((xd-xm[i]+L)/sigma[i]/qv[i])
       
                     ###
                     pf3 <- exp(L/sigma[i]*qv[i])
                     pf4 <- exp(L/sigma[i]/qv[i])

                     ### 
                     mat[,i] <- vm[i]*abs(pf1*qv[i]+pf2/qv[i])*(pf3+pf4)/abs(pf3*qv[i]+pf4/qv[i])/(pf1+pf2)*
                                exp(-0.5*(abs(log((pf1+pf2)/2)))^pv[i]+0.5*(abs(log((pf3+pf4)/2)))^pv[i])
                 
                 } # end for.

             } # end if.
             ### 

             ###
             if (all(is.finite(mat))) {

                 yd_fit <- rowSums(mat)
              
                 ### The objective to be minimised.
                 ###------------------------------ 
                 if (minfunc=="fom") {

                     rsdlv <- sqrt(abs(yd-yd_fit)/sum(yd_fit))

                 } else if (minfunc=="rss") {

                     rsdlv <- yd-yd_fit

                 } # end if. 

                 ###
                 return(rsdlv)
 
             } else {
 
                 return(hg) 

             } # end if. 

       } # end function minfn.
       ###
       ###--------------------------------------------------------------------------------
       ### Calculate shape parameters of grain-size components.
       ###-----------------------------------------------------
       calShape <- function(y, x)  {
           
           ny <- length(y)
           maxloc <- which.max(y)
           hmaxval <- max(y)/2
           Tm <- x[maxloc]
           
           ###
           T1 <- suppressWarnings(try(approx(x=y[1L:maxloc], y=x[1L:maxloc], xout=hmaxval)$y, silent=TRUE))
           T2 <- suppressWarnings(try(approx(x=y[maxloc:ny], y=x[maxloc:ny], xout=hmaxval)$y, silent=TRUE))

           ###
           if (inherits(T1,what="try-error")==FALSE) { d1 <- Tm-T1 } else { T1 <- d1 <- NA } # end if.
           if (inherits(T2,what="try-error")==FALSE) { d2 <- T2-Tm } else { T2 <- d2 <- NA } # end if. 
 
           ###          
           thw <- T2-T1
           sf <- d2/thw

           ###         
           return(c("x1"=T1, "x2"=T2, "xm"=Tm, "d1"=d1, "d2"=d2, "thw"=thw, "sf"=sf))

       } # end function calShape.
       ###--------------------------------------------------------------------------------

       ###
       lineCol <- c("deepskyblue", "orangered", "purple",   "violetred",  "yellowgreen", "lightblue",   
                    "goldenrod", "forestgreen", "blue",  "plum", "tan", "violet", "grey50") 

       ###
       gdcp <- 0
       minobj <- .Machine$double.xmax
       acpt <- 0

       ###
       if (viewLM==FALSE && viewPb==TRUE) {

           cat(paste("ncomp=",ncomp, ", model=", model, ", ntry=", ntry, ".\n", sep="")) 
           cat("Unmixing of single-sample GSD (ussGSD)"," is in progress, please wait, ...\n", sep="") 

           ###
           pb <- txtProgressBar(min=1, max=ntry, initial=1, style=3)
       } # end if. 

       ###
       fom0vec <- rss0vec <- c()
       ### Implement a trial-and-error protocol.
       ###---------------------------------------------------------------------------------
       for (i in seq(ntry))  {

           if (model %in% c("weibull0","lognormal0")) {

               p0 <- vector(length=3*ncomp)

           } else if (model %in% c("skewnormal0","skewgnormal0")) {

               p0 <- vector(length=4*ncomp)

           } # end if.

           ###
           if (is.null(startPars)) {

               if (model %in% c("weibull0","skewnormal0")) {
               
                   ### Simulate random alpha values.
                   ###-----------------------------
                   p0[(2*ncomp+1):(3*ncomp)] <- runif(n=ncomp, min=alphaRange[1], max=alphaRange[2])

               } else if (model %in% c("lognormal0","skewgnormal0")) {

                   ### Simulate random sigma values.
                   ###------------------------------
                   p0[(2*ncomp+1):(3*ncomp)] <- exp(runif(n=ncomp, min=log(sigmaRange[1]), max=log(sigmaRange[2])))

                   if (model=="skewgnormal0") {

                       p0[(2*ncomp+1):(3*ncomp)] <- (p0[(2*ncomp+1):(3*ncomp)])[mdgsORD]

                   } # end if.

               } # end if.

               ###
               if (model=="skewnormal0") {

                   ### Simulate random omega values.
                   ###------------------------------
                   p0[(3*ncomp+1):(4*ncomp)] <- exp(runif(n=ncomp, min=log(omegaRange[1]), max=log(omegaRange[2])))[mdgsORD]

               } else if (model=="skewgnormal0") {

                   ### Simulate random q values.
                   ###------------------------------
                   p0[(3*ncomp+1):(4*ncomp)] <- runif(n=ncomp, min=qRange[1], max=qRange[2])

               } # end if.

           } else {

               parsCOL3 <- startPars[,3,drop=TRUE]

               ### Tackle minus values.
               minmaxv <- t(apply(cbind((1-kkf)*parsCOL3,(1+kkf)*parsCOL3),MARGIN=1,sort))
               p0[(2*ncomp+1):(3*ncomp)] <- runif(n=ncomp, min=minmaxv[,1], max=minmaxv[,2])         

               ###
               if (model %in% c("skewnormal0","skewgnormal0")) {

                   parsCOL4 <- startPars[,4,drop=TRUE]

                   ###
                   p0[(3*ncomp+1):(4*ncomp)] <- runif(n=ncomp, min=(1-kkf)*parsCOL4, max=(1+kkf)*parsCOL4)

               } # end if.

           } # end if.
 
           ###
           p0[1:ncomp] <- runif(n=ncomp, min=(1-kkf)*mdgs[1:ncomp], max=(1+kkf)*mdgs[1:ncomp])

           ###
           p0[(ncomp+1):(2*ncomp)] <- runif(n=ncomp, min=(1-kkf)*mdvp[1:ncomp], max=(1+kkf)*mdvp[1:ncomp])

           ###
           ### Parameter optimisation using the Levenberg-Marquardt algorithm.
           ###----------------------------------------------------------------
           optLM <- try(minpack.lm::nls.lm(par=p0, lower=NULL, upper=NULL, fn=minfn, jac=NULL, 
                        control=nls.lm.control(maxiter=1024), xd, yd, mdgs, model, minfunc), silent=TRUE)
           ###
           if (viewLM==TRUE) print(optLM)

           ###
           if (inherits(optLM,what="try-error")==FALSE) {

               p <- optLM$par

               ###
               if (model=="weibull0") {

                   xm0 <- abs(p[1:ncomp]) 

                   ###
                   idxm1 <- which(xm0<(1-kkf)*mdgs)
                   xm0[idxm1] <- (1-kkf)*mdgs[idxm1]

                   ###
                   idxm2 <- which(xm0>(1+kkf)*mdgs)    
                   xm0[idxm2] <- (1+kkf)*mdgs[idxm2]

                   ###
                   vm0 <- abs(p[(ncomp+1):(2*ncomp)])

                   ###
                   alpha0 <- abs(p[(2*ncomp+1):(3*ncomp)])+1

                   ###
                   beta0 <- xm0/((alpha0-1.0)/alpha0)^(1.0/alpha0) 

                   ###
                   Mean0 <- beta0*gamma(1+1/alpha0)

                   ###
                   Median0 <- beta0*(log(2))^(1/alpha0)  
    
                   ###
                   mat0 <- matrix(nrow=nd, ncol=ncomp)

                   ###
                   for (j in 1:ncomp)  {

                       mat0[,j] <- vm0[j]*(xd/xm0[j])^(alpha0[j]-1)*exp((alpha0[j]-1)/alpha0[j]*(1-(xd/xm0[j])^alpha0[j]))

                   } # end for.

               } else if (model=="lognormal0") {

                   xm0 <- abs(p[1:ncomp])

                   ###
                   idxm1 <- which(xm0<(1-kkf)*mdgs)
                   xm0[idxm1] <- (1-kkf)*mdgs[idxm1]

                   ###
                   idxm2 <- which(xm0>(1+kkf)*mdgs)    
                   xm0[idxm2] <- (1+kkf)*mdgs[idxm2]

                   ###
                   vm0 <- abs(p[(ncomp+1):(2*ncomp)])

                   ###
                   sigma0 <- abs(p[(2*ncomp+1):(3*ncomp)])+0.001

                   ###
                   theta0 <- log(xm0)+sigma0^2

                   ###
                   Mean0 <- exp(theta0+0.5*sigma0^2)

                   ###
                   Median0 <- exp(theta0)

                   ###
                   mat0 <- matrix(nrow=nd, ncol=ncomp)

                   ###
                   for (j in 1:ncomp)  {

                       mat0[,j] <- xm0[j]/xd*vm0[j]*exp(-log(xm0[j]/xd)*(log(xm0[j]/xd)+2*(sigma0[j])^2)/2/(sigma0[j])^2)

                   } # end for.

               } else if (model=="skewnormal0") {

                   xm0 <- abs(p[1:ncomp])

                   ###
                   idxm1 <- which(xm0<(1-kkf)*mdgs)
                   xm0[idxm1] <- (1-kkf)*mdgs[idxm1]

                   ###
                   idxm2 <- which(xm0>(1+kkf)*mdgs)    
                   xm0[idxm2] <- (1+kkf)*mdgs[idxm2]

                   ###
                   vm0 <- abs(p[(ncomp+1):(2*ncomp)])

                   ###
                   alpha0 <- p[(2*ncomp+1):(3*ncomp)]

                   ###
                   omega0 <- abs(p[(3*ncomp+1):(4*ncomp)])+0.001

                   ###
                   delta0 <- alpha0/sqrt(1+(alpha0)^2)

                   ###
                   v10 <- sqrt(1-2*delta0^2/pi)*(4-pi)/4*(delta0*sqrt(2/pi))^3/(1-2*delta0^2/pi)^1.5

                   ###
                   v20 <- sign(alpha0)/2*exp(-2*pi/abs(alpha0))

                   ###
                   D0 <- omega0*(sqrt(2/pi)*delta0-v10-v20)

                   ###
                   Mean0 <- xm0-D0+omega0*delta0*sqrt(2/pi)

                   ###
                   Median0 <- xm0

                   ###
                   mat0 <- matrix(nrow=nd, ncol=ncomp)

                   ###
                   for (j in 1:ncomp)  {

                       mat0[,j] <- vm0[j]*exp(-(2*D0[j]+xd-xm0[j])*(xd-xm0[j])/2/(omega0[j])^2)*
                      (1+erf(alpha0[j]*(xd-xm0[j]+D0[j])/sqrt(2)/omega0[j]))/(1+erf(alpha0[j]*D0[j]/sqrt(2)/omega0[j]))

                   } # end for. 

               } else if (model=="skewgnormal0") {

                   xm0 <- abs(p[1:ncomp])

                   ###
                   idxm1 <- which(xm0<(1-kkf)*mdgs)
                   xm0[idxm1] <- (1-kkf)*mdgs[idxm1]

                   ###
                   idxm2 <- which(xm0>(1+kkf)*mdgs)    
                   xm0[idxm2] <- (1+kkf)*mdgs[idxm2]

                   ###
                   vm0 <- abs(p[(ncomp+1):(2*ncomp)])

                   ###
                   sigma0 <- abs(p[(2*ncomp+1):(3*ncomp)])+0.001

                   ###
                   qv0 <- abs(p[(3*ncomp+1):(4*ncomp)])+0.001
                   
                   ###
                   pv0 <- 2 + 6*(1-qv0)^5

                   ###
                   ###x0 <- seq(from=0,to=MAXX/3,by=MAXX/3/300)
                   x0 <- seq(from=0,to=MAXX/4,by=MAXX/4/200)
                   
                   ###
                   L0 <- vector(length=ncomp)

                   ###
                   mat0 <- matrix(nrow=nd, ncol=ncomp)

                   ###
                   for (j in 1:ncomp)  {

                       yf1 <- exp(x0/sigma0[j]*qv0[j])
                       yf2 <- exp(x0/sigma0[j]/qv0[j])

                       ###
                       y0 <- abs(yf1*qv0[j]+yf2/qv0[j])/(yf1+yf2)*exp(-0.5*(abs(log((yf1+yf2)/2)))^pv0[j])

                       ###
                       if (any(!is.finite(y0)))  { 

                           mat0[,j] <- NA 

                       } else {

                           L0[j] <- x0[which.max(y0)]

                           ###
                           pf1 <- exp((xd-xm0[j]+L0[j])/sigma0[j]*qv0[j])
                           pf2 <- exp((xd-xm0[j]+L0[j])/sigma0[j]/qv0[j])

                           ###
                           pf3 <- exp(L0[j]/sigma0[j]*qv0[j])
                           pf4 <- exp(L0[j]/sigma0[j]/qv0[j])

                           ###
                           mat0[,j] <- vm0[j]*abs(pf1*qv0[j]+pf2/qv0[j])*(pf3+pf4)/abs(pf3*qv0[j]+pf4/qv0[j])/(pf1+pf2)*
                                       exp(-0.5*(abs(log((pf1+pf2)/2)))^pv0[j]+0.5*(abs(log((pf3+pf4)/2)))^pv0[j])
                       
                       } # end if.

                   } # end for.

                   ###
                   kurt0 <- 2-pv0

                   ###
                   skew0 <- -6*sign(qv0)*(1-qv0)^2*(1+1.856*kurt0)

                   ###
                   Mean0 <- xm0-L0+skew0/6*(1+0.856*kurt0)

                   ###
                   Median0 <- xm0        

               } # end if.

               ###
               ### Accept the unmixing results of the "trial-and-error" protocol 
               ### if the following conditions are satisfied.
               ###------------------------------------------------------------------
               OOKK <- all(is.finite(mat0)) && all(is.finite(Mean0)) && all(Mean0>MINX & Mean0<MAXX) &&
                       all(is.finite(Median0)) && all(Median0>MINX & Median0<MAXX)

               ###
               if (OOKK==TRUE)  {

                   acpt <- acpt + 1

                   ###
                   yd_fit0 <- rowSums(mat0)
                   fom0 <- sum(abs(yd-yd_fit0))/sum(yd_fit0)*100
                   rss0 <- sum((yd-yd_fit0)^2)

                   ###
                   fom0vec <- c(fom0vec, fom0)
                   rss0vec <- c(rss0vec, rss0)

                   ### Visualise the unmixing results of the "trial-and-error" protocol.
                   ###------------------------------------------------------------------ 
                   if (viewFit==TRUE) {

                       mat_vfp <- mat0
                       mat_vfp <- mat_vfp[,order(xm0),drop=FALSE] 

                       ###
                       par(mar=c(6,6,5,5)+0.1)

                       ###
                       plot(origin_xd[reserveidx], yd[reserveidx], type="n", cex.lab=1.5, 
                            log="x", cex.axis=1.5, mgp=c(2.5, 1, 0), xlab="Grain size (um)", 
                            ylab="Volume percentage (%)", xaxt="n")

                       ###
                       box(lwd=2)

                       ###
                       XaxisCentral <- median(axTicks(side=1))

                       ###
                       axis(side=1, at=myTK, labels=myLB, cex.axis=1.5)

                       ###
                       points(origin_xd[reserveidx], yd[reserveidx], type="p",  lwd=2, pch=21, col="gray50", cex=1.5) 

                       ###
                       points(origin_xd, rowSums(mat_vfp), col="grey30", type="l", lwd=3.0)

                       ###
                       for (j in 1:ncomp) {

                           points(origin_xd, mat_vfp[,j], col=lineCol[j], type="l", lwd=3)
                       
                       } # end for.
                          
                       ###
                       mypos <- mean(yd[reserveidx][origin_xd[reserveidx]<=XaxisCentral])<
                                mean(yd[reserveidx][origin_xd[reserveidx]>XaxisCentral]) 

                       ###
                       legend(ifelse(mypos,"topleft","topright"),
                              legend=c("Measured","Fitted", paste("Comp.",seq(ncomp),sep=""),
                              paste("Count:",acpt,sep=""), paste("RSS:",round(rss0,2),sep=""), 
                              paste("FOM:",round(fom0,2),sep="")), 
                              col=c("grey50", "grey30", lineCol[seq(ncomp)],NA,NA,NA), 
                              pch=c(21, rep(NA,ncomp+1),NA,NA,NA), 
                              lty=c(NA, rep("solid",ncomp+1),NA,NA,NA), 
                              yjust=2, ncol=1, cex=1, pt.cex=1, bty="o", 
                              lwd=3.0, pt.bg="white")

                   } # end if. 

                   ###
                   vxfght <- ifelse(minfunc=="fom", fom0, rss0)            

                   ### Update the optimal results.
                   if (vxfght<minobj) {

                       xm <- xm0              
                       vm <- vm0

                       ###
                       if (model=="weibull0") {

                           alpha <- alpha0
                           beta <- beta0

                       } else if (model=="lognormal0") {

                           sigma <- sigma0
                           theta <- theta0
                   
                       } else if (model=="skewnormal0") {

                           alpha <- alpha0
                           omega <- omega0
                           delta <- delta0
 
                           ###
                           D <- D0

                       } else if (model=="skewgnormal0") {

                           sigma <- sigma0
                           qv <- qv0
                           pv <- pv0
                           L <- L0

                           ###
                           kurt <- kurt0
                           skew <- skew0

                       } # end if.

                       ###
                       Mean <- Mean0
                       Median <- Median0
                       Mode <- xm0
                       
                       ###
                       mat <- mat0
                       minobj <- vxfght 
                      
                       ###
                       gdcp <- 1 
                       
                   } # end if.

               } # end if.

           } # end if.

           ###
           if (viewLM==FALSE && viewPb==TRUE)  setTxtProgressBar(pb, i)

       } # end for.
       ###---------------------------------------------------------------------------------
    
       ###
       if (viewLM==FALSE && viewPb==TRUE)  close(pb) 

       ###
       if (gdcp==0) {

           mat <- cbind(origin_xd, yd, rep(NA,nd))
           rownames(mat) <- NULL
           colnames(mat) <- c("GSlev","Volume","Fit.Volume")

           ###
           output <- list("model"=model, "ncomp"=ncomp, "expGSlev"=expGSlev, 
                          "reserveidx"=reserveidx, "gs.comp"=mat, "FOM"=Inf, "RSS"=Inf)
           class(output) <- "ussgsd"

           ###
           if (plot==TRUE) {

               plot_ussGSD(output, sampleName=sampleName, addvl=addvl, lwd=lwd, cex=cex, pch=pch)

           } # end if.

           ###
           return(invisible(output))

       } # end if.

       ###
       if (model %in% c("weibull0","lognormal0")) {

           Np <- 3*ncomp

       } else if (model %in% c("skewnormal0","skewgnormal0")) {

           Np <- 4*ncomp

       } # end if.

       ###
       yd_fit <- rowSums(mat)
       residuals <- yd-yd_fit
       fom <- sum(abs(yd-yd_fit))/sum(yd_fit)*100
       rss <- sum((yd-yd_fit)^2)
       sst <- sum((yd-mean(yd))^2)
       R2 <- 1 - rss/sst
       RSE <- sqrt(rss/(nd-Np))

       ###
       if (model=="weibull0") {

           coef <- vm*xm/(alpha-1)/exp(-(alpha-1)/alpha)

       } else if (model=="lognormal0") {

           coef <- vm*xm*sqrt(2*pi)*sigma/exp(-0.5*sigma^2)

       } else if (model=="skewnormal0") {

           coef <- vm*omega*sqrt(2*pi)/exp(-0.5*D^2/omega^2)/(1+erf(alpha*D/sqrt(2)/omega))

       } else if (model=="skewgnormal0") {

           coef <- vm*2^(1+1/pv)*sigma*gamma(1+1/pv)*(exp(qv*L/sigma)+exp((1/qv)*L/sigma))/    
           exp(-0.5*(abs(log(0.5*(exp(qv*L/sigma)+exp((1/qv)*L/sigma)))))^pv)/abs(qv*exp(qv*L/sigma)+1/qv*exp((1/qv)*L/sigma))

       } # end if.

       ###
       prop <- coef/sum(coef)*100

       ### Numerical approximation of the medians.
       ###-------------------------------------------------
       if (model %in% c("skewnormal0", "skewgnormal0")) {

           Median <- vector(length=ncomp)

           ###
           for (j in 1:ncomp) {

               Median[j] <- approx(x=cumsum(mat[,j]/sum(mat[,j])),y=xd, xout=0.5,ties="ordered")$y

           } # end for.

       } # end if.

       ###
       gs.pars <- cbind(prop, Mean, Median, Mode)
       oidx <- order(gs.pars[,4])
       gs.pars <- gs.pars[oidx,,drop=FALSE]

       ###
       mat <- mat[,oidx,drop=FALSE]
       sp <- t(apply(mat, MARGIN=2L, calShape, seq(xd)))
       rownames(sp) <- paste("Comp.", seq(ncomp), sep="")

       ###
       if (ncomp>1) {

           resolvec <- name_resolvec <- vector(length=ncomp-1) 
           deltaHvec <- name_deltaHvec <- vector(length=ncomp-1)

           ###
           for (i in seq(ncomp-1)) {

               resolvec[i] <- (sp[i+1,3]-sp[i,3])/(sp[i,5]+sp[i+1,4])
               name_resolvec[i] <- paste("Comp.",i,".",i+1L,sep="")

               ###
               px12idx <- mat[,i]>0 & mat[,i+1]>0
               px1 <- mat[px12idx,i]/sum(mat[px12idx,i])
               px2 <- mat[px12idx,i+1]/sum(mat[px12idx,i+1])
               px12 <- px1 + px2
       
               ###
               deltaHvec[i] <- (-sum(px12*log2(px12))) - (-sum(px1*log2(px1))) - (-sum(px2*log2(px2)))
               name_deltaHvec[i] <- paste("Comp.",i,".",i+1L,sep="")

           } # end for.

           ###
           names(resolvec) <- name_resolvec
           names(deltaHvec) <- name_deltaHvec

       } else {

           resolvec <- deltaHvec <- Inf

       } # end if.
 
       ###
       fit.data <- cbind(xd,yd)
       rownames(fit.data) <- colnames(fit.data) <- NULL
       ###

       ###
       if (model=="weibull0") {

           ###
           SDp <- sqrt( beta^2*(gamma(1+2/alpha)-(gamma(1+1/alpha))^2) )

           ###
           fit.pars <- cbind(xm, vm, alpha)
           colnames(fit.pars) <- c("Xm", "Vm", "Alpha")

       } else if (model=="lognormal0") {
          
           ###
           SDp <- sqrt( (exp(sigma^2)-1)*exp(2*theta+sigma^2) )

           ###
           fit.pars <- cbind(xm, vm, sigma)
           colnames(fit.pars) <- c("Xm", "Vm", "Sigma")

       } else if (model=="skewnormal0") {

           ###
           SDp <- sqrt(omega^2*(1-2*delta^2/pi))

           ###
           fit.pars <- cbind(xm, vm, alpha, omega)
           colnames(fit.pars) <- c("Xm", "Vm", "Alpha", "Omega")

       } else if (model=="skewgnormal0") {

           ###
           SDp <- sqrt(sigma^2*abs((1+1.856*kurt)*(1-abs(skew)/3)))

           ###
           fit.pars <- cbind(xm, vm, sigma, qv)
           colnames(fit.pars) <- c("Xm", "Vm", "Sigma", "q")

       } # end if.

       ###
       fit.pars <- fit.pars[oidx,,drop=FALSE]
       rownames(fit.pars) <- paste("Comp.", seq(ncomp), sep="")

       ###
       SDp <- SDp[oidx]

       ###
       normTEST <- unclass(shapiro.test(x=residuals))

       ###
       if (useIndex==TRUE) {

           if (expGSlev==TRUE) {

               gs.pars[,2] <- exp(pab[1]+pab[2]*gs.pars[,2])
               gs.pars[,3] <- exp(pab[1]+pab[2]*gs.pars[,3])
               gs.pars[,4] <- exp(pab[1]+pab[2]*gs.pars[,4]) 

               ###
               SDp <-  gs.pars[,2]*pab[2]*SDp      

           } else {

               gs.pars[,2] <- pab[1]+pab[2]*gs.pars[,2]
               gs.pars[,3] <- pab[1]+pab[2]*gs.pars[,3]
               gs.pars[,4] <- pab[1]+pab[2]*gs.pars[,4]

               ###
               SDp <- pab[2]*SDp 

           } # end if.

       } # end if.

       ###
       gs.pars <- cbind(gs.pars, SDp)

       ###
       colnames(gs.pars) <- c("Proportion", "Mean", "Median", "Mode", "SD")   
       rownames(gs.pars) <- paste("Comp.", seq(ncomp), sep="")

       ###
       ### Implement a recursive unmixing protocol.
       ###-----------------------------------------------------------------
       if (trim==TRUE) {

            if (is.null(startPars)) {

               weakPK <- which(peakIDX[,1]<=0)

           } else {

               weakPK <- NA

           } # end if. 

           ###
           if ((all(is.finite(resolvec)) && min(resolvec)<mrsl) && ncomp>1 && length(weakPK)>=1) {

               ncomp <- ncomp - 1

               ###
               overLAPidx <- which.min(resolvec)

               ###
               vecxxx <- c(overLAPidx,overLAPidx+1)

               ###
               leaveCOMPidx <- vecxxx[which.min(gs.pars[vecxxx,1])]

               ### 
               if (model %in% c("weibull0","lognormal0")) {
          
                   startPars <- cbind(gs.pars[,4], fit.pars[,2:3])

               } else if (model %in% c("skewnormal0","skewgnormal0")) {

                   startPars <- cbind(gs.pars[,4], fit.pars[,2:4])

               } # end if.

               ###
               startPars <- startPars[-leaveCOMPidx,,drop=FALSE]

               ### 
               recursiveFIT <- ussGSD0(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                 mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                 sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                 viewAutoInis=viewAutoInis, viewLM=viewLM, viewFit=viewFit, viewPb=viewPb,  
                 outfile=outfile, plot=plot, rmZero=rmZero, sampleName=sampleName, addvl=addvl, 
                 logxy=logxy, lwd=lwd, pch=pch, cex=cex)

               ###
               return(invisible(recursiveFIT))
      
           } # end if.

       } # end if.

       ###
       tmdxv <- gs.pars[,4]
       names(tmdxv) <- NULL
       tmdyv <- vector(length=ncomp)

       ###     
       for (i in seq(ncomp))  {

           tmdyv[i] <- approx(x=origin_xd, y=mat[,i], xout=tmdxv[i])$y

       } # end for. 

       ###
       mat <- cbind(origin_xd, yd, yd_fit, mat)
       rownames(mat) <- NULL
       colnames(mat) <- c("GSlev","Volume","Fit.Volume",paste("Comp.", seq(ncomp), sep=""))

       ###
       if (!is.null(outfile))  write.csv(mat, file=paste(outfile, ".csv", sep = ""))

       ###
       output <- list("model"=model, "ncomp"=ncomp, "expGSlev"=expGSlev, "reserveidx"=reserveidx,
                      "fit.data"=fit.data, "fit.pars"=fit.pars, "normTest"=normTEST, "FOMs"=fom0vec, 
                      "RSSs"=rss0vec, "ntry"=ntry, "acpt"=acpt, "gs.comp"=mat, "mdxv"=tmdxv, "mdyv"=tmdyv, 
                      "residuals"=residuals, "gs.pars"=gs.pars, "sp"=sp, "rsl"=resolvec, "deltaH"=deltaHvec, 
                      "FOM"=fom, "RSS"=rss, "R2"=R2, "RSE"=RSE)

       ###
       class(output) <- "ussgsd"

       ###
       if (plot==TRUE) {

           plot_ussGSD(obj_ussgsd=output, sampleName=sampleName, addvl=addvl, logxy=logxy, lwd=lwd, cex=cex, pch=pch)

       } # end if.

       ###
       invisible(output)

  } # end function ussGSD0.
  ###==========================================================================================================================
  ###
  
  
  ###
  ###**************************************************************************************************************************
  ### Function update_ussGSDbatchp() is used to parallelly update (re-analyze)   
  ### the unmixing results of some specified grain-size distributions. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ### obj_batchgsd: an S3 object of class "gsdbatch" generated using the function 
  ###               ussGSDbatch(), ussGSDbatchp(), update_ussGSDbatch(), or update_ussGSDbatchp().
  ###
  ###     sampleNO: An integer vector indicating the ID numbers of the grain-size distributions to be updated (re-analyzed). 
  ###
  ###        ncomp: An integer (with values range from 1 to 13) indicating the number of 
  ###               components for the grain-size distributions to be re-analysed. 
  ###
  ###         auto: A logical value indicating whether performing a automatic grain-size unmixing, 
  ###               in this case the user needs not to specify the initials used for optimisation.
  ###
  ###        model: A character indicating the model to be fitted, "weibull", "lognormal",
  ###               "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###               if model=NULL, the program will automatically determine the optimal model 
  ###               (between "weibull" and "lognormal") which yields a smaller FOM or RSS value.
  ###
  ###          mpd: A real number indicating the allowed lower limit of the volume percentage 
  ###               of a GSD that can be used to identify the location of component peaks.
  ###
  ###          ctf: A numeric value (between -1 and 1) representing a critical threshold factor  
  ###               that controls the identification of peaks from the grain-size distribution.  
  ###               This argument can be used to prevent identifying a false peak characterised  
  ###               by a large positive second-order derivative. Specially, ctf=0 indicates
  ###               that peaks with second-order derivatives above zero will be precluded. 
  ###               The precluding effect will be increasingly suppressed as ctf increases 
  ###               from -1 to 1.
  ###
  ###         ntry: An integer value indicating the number of trials in a trial-and-error protocol. 
  ###
  ###          kkf: A numeric value controlling the range of values from which random starting   
  ###               parameters will be generated during the "trial-and-error" protocol.
  ###
  ###    startPars: if model="weibull" or "lognormal", [startPars] should be a two-column matrix 
  ###               containing starting parameters used for unmixing, the first column contains 
  ###               the modes of individual components, and the second row contains the alpha 
  ###               ("weibull") or sigma ("lognormal") value of individual components.
  ###
  ###               if model="weibull0" or "lognormal0", [startPars] should be a three-column 
  ###               matrix containing starting parameters used for unmixing, the first column  
  ###               contains the modes of individual components, the second column contains the 
  ###               maximum volume percentages of individual components, and the thrid column 
  ###               contains the alpha ("weibull0") or sigma ("lognormal0") value of individual 
  ###               components.
  ###
  ###               if model="skewnormal0" or "skewgnormal0", [startPars] should be a four-column 
  ###               matrix containing starting parameters used for unmixing, the first column is
  ###               the modes of individual components, the second column is the maximum volume
  ###               percentages of individual components,the thrid and fourth columns are the
  ###               alpha and omega values ("skewnormal0") or sigma and q values ("skewgnormal0") 
  ###               of individual components.  
  ###
  ###   alphaRange: A two-element vector indicating the lower and upper limits on alpha values of 
  ###               the Weibull ("weibull" or "weibull0") or Skew Normal ("skewnormal0") distributions 
  ###               generated from a Uniform distribution.
  ###
  ###   sigmaRange: A two-element vector indicating the lower and upper limits on sigma values 
  ###               of the Lognormal ("lognormal" or "lognormal0") or Skewed Generalized Normal 
  ###              ("skewgnormal0") distributions generated from a Uniform distribution.  
  ###
  ###   omegaRange: A two-element vector indicating the lower and upper limits on omega values of 
  ###               a Skew Normal ("skewnormal0") distribution generated from a Uniform distribution.
  ###
  ###       qRange: A two-element vector indicating the lower and upper limits on q values of a Skewed 
  ###               Generalized Normal ("skewgnormal0") distribution generated from a Uniform distribution. 
  ###
  ###     useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###               used as the x-coordinate during the unmixing process.
  ###
  ###      minfunc: A character indicating the objective to be mimimized, either "fom" or "rss", 
  ###               for the figure-of-merit value or the residual sum of squares, respectively.
  ###
  ###         trim: A logical value indicating whether the unmixing results will be trimed  
  ###               using a recursive optimization protocol.
  ###
  ###         mrsl: A numeric value indicating the minimum resolution between two adjacent     
  ###               components used to trim the unmixing results, which ranges from 0 to 2.
  ###               Peaks with a second-order derivative below zero and a resolution smaller
  ###               than mrsl will be deleted.
  ###
  ###          ncr: The number of cores to be used during the parallel unmixing process. If
  ###               the total number of available cores (NumberOfCluster) is smaller than ncr, 
  ###               NumberOfCluster will be used instead of ncr.
  ###
  ###       rmZero: A logical value indicating whether zero volume percentages will be removed.
  ###
  ###       saveRD: A logical value indicating whether the updated results will be saved to
  ###               the loaded RDdata file in the current directory.
  ###==========================================================================================================================
  ### The function
  ### (1) Returns a invisible list of S3 class of "batchgsd" containing the unmixed results of individual samples.
  ###
  ### (2) Re-write the loaded RData file containing the updated unmixed results in the current working directory.
  ###
  ### Note that if obj_batchgsd=NULL, the user needs to ensure that the function load_ussGSDbatch() 
  ### has been called to import a RData file containing an object of S3 class "batchgsd" that is 
  ### available from the current working directory.
  ###==========================================================================================================================
  update_ussGSDbatchp <- function(obj_batchgsd=NULL, sampleNO, ncomp=NULL, auto=TRUE, model="weibull0", 
                                  mpd=0, ctf=0.1, ntry=50, kkf=0.1, startPars=NULL, alphaRange=NULL, 
                                  sigmaRange=NULL, omegaRange=NULL, qRange=NULL, useIndex=TRUE, minfunc="fom", 
                                  trim=FALSE, mrsl=0.6, ncr=NULL, rmZero=TRUE, saveRD=TRUE) {

      stopifnot(is.numeric(sampleNO), 
                is.null(ncomp) || (length(ncomp)==1 && is.numeric(ncomp) && ncomp %in% 1:13),
                is.logical(auto), length(auto)==1,
                is.null(model) || (length(model)==1 && model %in% c("weibull","lognormal","weibull0", "lognormal0", "skewnormal0", "skewgnormal0")),
                is.numeric(mpd), length(mpd)==1, mpd>=0, mpd<3,
                is.numeric(ctf), length(ctf)==1, ctf>=-1, ctf<=1,
                is.numeric(ntry), length(ntry)==1, ntry>=2, 
                is.numeric(kkf), length(kkf)==1, kkf>0, kkf<1,
                is.null(startPars) || (is.matrix(startPars) && ncol(startPars)==2), 
                is.null(alphaRange) || (is.numeric(alphaRange) && length(alphaRange)==2),
                is.null(sigmaRange) || (is.numeric(sigmaRange) && length(sigmaRange)==2),
                is.null(omegaRange) || (is.numeric(omegaRange) && length(omegaRange)==2),
                is.null(qRange) || (is.numeric(qRange) && length(qRange)==2),
                is.logical(useIndex), length(useIndex)==1,
                is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "rss"),
                is.logical(trim), length(trim)==1,
                is.numeric(mrsl), length(mrsl)==1, mrsl>=0, mrsl<=2,
                is.null(ncr) || (length(ncr)==1 && is.numeric(ncr)), 
                is.logical(rmZero), length(rmZero)==1,  
                is.logical(saveRD), length(saveRD)==1)

      if (!is.null(obj_batchgsd)) {

          if(inherits(obj_batchgsd,what="batchgsd")==FALSE)  stop("Error: [obj_batchgsd] should be an S3 object of class 'batchgsd'!")
          GSDbatch <- obj_batchgsd

      } else {

          if (!exists("gsdbatch"))  stop("Error: function [load_ussGSDbatch] has not been called!")
          GSDbatch <- get("GSDbatch", envir=gsdbatch)

      } # end if.

      ###
      xoxoxo <- attr(GSDbatch, "xoxoxo")

      ###
      allmyNO <- 1:length(GSDbatch)
      
      ###
      if (!all(sampleNO %in% allmyNO)) {

          XNO <- sampleNO[!(sampleNO %in% allmyNO)]

          ###
          cat("The following sample numbers are invalid:\n")
          print(XNO)

          ###
          stop()

      } # end if.

      ###
      N <- length(sampleNO)

      ###
      if (N<9)  stop("Error: Parallel optimisation is forbidden for a sample size of N<9!")

      ###
      NumberOfCluster <- parallel::detectCores()

      ###
      if (!is.null(ncr)) {

          if (!ncr %in% (2:1000))  stop("Error: [ncr] must be an integer ranging from 2 to 1000!")

          ###
          if (ncr<=NumberOfCluster) {

              NumberOfCluster <- ncr

          } else {

              cat("Note: [ncr=",ncr,"] exceeds the number of available cores (",NumberOfCluster,")!\n",sep="")

          } # end if.

      } # end if.
               
      ###
      cl <- snow::makeCluster(NumberOfCluster, outfile="")
      doSNOW::registerDoSNOW(cl)

      ###
      cat("Parallel unmixing of N=", N, " GSDs in a batch pattern is in progress, please wait, ...\n", sep="")

      ###
      pb <- txtProgressBar(max=N, style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)

      ### Parallel unmixing using the R function ussGSD(). 
      ###-------------------------------------------------------------------------------------
      updateGSDbatch <- foreach::foreach(kk=1:N, .inorder=TRUE, .options.snow=opts, .packages=c("pracma","minpack.lm"), 
                                         .export=c("ussGSD","ussGSD0"), .verbose=FALSE) %dopar% {

          ###
          gsl <- GSDbatch[[sampleNO[kk]]]$gs.comp[,1]
          gsd <- GSDbatch[[sampleNO[kk]]]$gs.comp[,2]

          ###
          if (!is.null(model)) {

              if (model %in% c("weibull","lognormal")) {

                  res_ussGSD <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                    viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                    rmZero=rmZero, plot=FALSE), silent=TRUE)

              } else if (model %in% c("weibull0", "lognormal0", "skewnormal0", "skewgnormal0")) {

                  res_ussGSD <- try(ussGSD0(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=auto, model=model, 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, omegaRange=omegaRange, qRange=qRange, useIndex=useIndex, 
                    minfunc=minfunc, trim=trim, mrsl=mrsl, viewAutoInis=FALSE, viewLM=FALSE, 
                    viewFit=FALSE, viewPb=FALSE, outfile=NULL, rmZero=rmZero, plot=FALSE), 
                    silent=TRUE)

              } # end if.
    
              ###
              if (inherits(res_ussGSD,what="try-error")==FALSE) { 

                  return(res_ussGSD)

              } else {

                  cat("Failed in unmixing the ", sampleNO[kk], "-th sample!\n", sep="")
                  print(paste("<",sampleNO[kk],"> ",attr(res_ussGSD,"condition"),sep=""))
                  return(NULL)

              } # end if. 

          } else {

              res_ussGSDwb <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp,  auto=auto, model="weibull", 
                mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                rmZero=rmZero, plot=FALSE), silent=TRUE)

              ###
              res_ussGSDlg <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp,  auto=auto, model="lognormal", 
                mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                rmZero=rmZero, plot=FALSE), silent=TRUE)

              ###
              if (inherits(res_ussGSDwb,what="try-error")==FALSE && inherits(res_ussGSDlg,what="try-error")==FALSE) {

                  if (minfunc=="fom") {

                      minfuncVAL1 <- res_ussGSDwb$FOM
                      minfuncVAL2 <- res_ussGSDlg$FOM

                  } else if (minfunc=="rss") {

                      minfuncVAL1 <- res_ussGSDwb$RSS
                      minfuncVAL2 <- res_ussGSDlg$RSS

                  } # end if.

                  ###
                  if (minfuncVAL1<minfuncVAL2) {

                      return(res_ussGSDwb)

                  } else {

                      return(res_ussGSDlg)

                  } # end if. 

              } else if (inherits(res_ussGSDwb,what="try-error")==FALSE) {

                  cat("Failed in unmixing the ", sampleNO[kk], "-th sample using Lognormal!\n", sep="")
                  print(paste("<",sampleNO[kk],"> ",attr(res_ussGSDlg,"condition"),sep=""))
                  return(res_ussGSDwb)

              } else if (inherits(res_ussGSDlg,what="try-error")==FALSE) {

                  cat("Failed in unmixing the ", sampleNO[kk], "-th sample using Weibull!\n", sep="")
                  print(paste("<",sampleNO[kk],"> ",attr(res_ussGSDwb,"condition"),sep=""))
                  return(res_ussGSDlg)

              } else {

                  cat("Failed in unmixing the ", sampleNO[kk], "-th sample using Weibull and Lognormal!\n", sep="")
                  print(paste("<",sampleNO[kk],"> ",attr(res_ussGSDwb,"condition"),sep=""))
                  print(paste("<",sampleNO[kk],"> ",attr(res_ussGSDlg,"condition"),sep=""))
                  return(NULL)

              } # end if.           

          } # end if.

      } # end foreach. 

      ###
      cat("\n")
      close(pb)

      ###
      foreach::registerDoSEQ()
      snow::stopCluster(cl)

      ###
      for (kk in 1:N) {

          GSDbatch[[sampleNO[kk]]] <- updateGSDbatch[[kk]] 

      } # end for.

      ###
      assign("GSDbatch", GSDbatch, envir=gsdbatch)
      if (saveRD==TRUE)  save(GSDbatch, file=paste(xoxoxo, ".RData", sep=""))

      ###
      invisible(GSDbatch)

  } # end function update_ussGSDbatchp.
  ###==========================================================================================================================
  ###


  ###
  ###**************************************************************************************************************************
  ### Function ussGSDbatchp() is used for parallel unmixing of single-sample  
  ### grain-size distributions for a number of samples in a batch pattern. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###   GSDfile: A CSV file containing the grain-size data used for analysis that is available from
  ###            the current working directory, or a data.frame (or matrix) imported using function 
  ###            read.table(). The first row is the grain-size levels, and the remaining rows are the 
  ###            volume percentages of individual samples to be unmixed. If missing, a CSV template 
  ###            will be generated automatically to guide the user to prepare the input dataset.
  ###
  ###     ncomp: An integer (from 1 to 13) indicating the number of components to be fitted. 
  ###
  ###     model: A character indicating the model to be fitted, "weibull", "lognormal",
  ###            "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###            if model=NULL, the program will automatically determine the optimal model 
  ###            (between "weibull" and "lognormal") which yields a smaller FOM or RSS value.
  ###
  ###       mpd: A real number indicating the allowed lower limit of the volume percentage 
  ###            of a GSD that can be used to identify the location of component peaks.
  ###
  ###       ctf: A numeric value (between -1 and 1) representing a critical threshold factor  
  ###            that controls the identification of peaks from the grain-size distribution.  
  ###            This argument can be used to prevent identifying a false peak characterised  
  ###            by a large positive second-order derivative. Specially, ctf=0 indicates
  ###            that peaks with second-order derivatives above zero will be precluded. 
  ###            The degree of precluding will be increasingly suppressed as ctf increases 
  ###            from -1 to 1.
  ###
  ###      ntry: An integer indicating the number of trials in a trial-and-error protocol. 
  ###
  ###       kkf: A numeric value controlling the range of values from which random starting   
  ###            parameters will be generated during the "trial-and-error" protocol.
  ###
  ### startPars: if model="weibull" or "lognormal", [startPars] should be a two-column matrix 
  ###            containing starting parameters used for unmixing, the first column contains 
  ###            the modes of individual components, and the second row contains the alpha 
  ###            ("weibull") or sigma ("lognormal") value of individual components.
  ###
  ###            if model="weibull0" or "lognormal0", [startPars] should be a three-column 
  ###            matrix containing starting parameters used for unmixing, the first column  
  ###            contains the modes of individual components, the second column contains the 
  ###            maximum volume percentages of individual components, and the thrid column 
  ###            contains the alpha ("weibull0") or sigma ("lognormal0") value of individual 
  ###            components.
  ###
  ###            if model="skewnormal0" or "skewgnormal0", [startPars] should be a four-column 
  ###            matrix containing starting parameters used for unmixing, the first column is
  ###            the modes of individual components, the second column is the maximum volume
  ###            percentages of individual components,the thrid and fourth columns are the
  ###            alpha and omega values ("skewnormal0") or sigma and q values ("skewgnormal0") 
  ###            of individual components.  
  ###
  ###alphaRange: A two-element vector indicating the lower and upper limits on alpha values of 
  ###            the Weibull ("weibull" or "weibull0") or Skew Normal ("skewnormal0") distributions 
  ###            generated from a Uniform distribution.
  ###
  ###sigmaRange: A two-element vector indicating the lower and upper limits on sigma values 
  ###            of the Lognormal ("lognormal" or "lognormal0") or Skewed Generalized Normal 
  ###           ("skewgnormal0") distributions generated from a Uniform distribution.  
  ###
  ###omegaRange: A two-element vector indicating the lower and upper limits on omega values of 
  ###            a Skew Normal ("skewnormal0") distribution generated from a Uniform distribution.
  ###
  ###    qRange: A two-element vector indicating the lower and upper limits on q values of a Skewed 
  ###            Generalized Normal ("skewgnormal0") distribution generated from a Uniform distribution.
  ###
  ###  useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###            used as the x-coordinate during the fitting process.
  ###
  ###   minfunc: A character indicating the objective to be mimimized, either "fom" or "rss",
  ###            for the figure-of-merit value or the residual sum of squares, respectively.
  ###
  ###      trim: A logical value indicating whether the unmixing results will be trimed  
  ###            using a recursive optimisation protocol.
  ###
  ###      mrsl: A numeric value indicating the minimum resolution between two adjacent     
  ###            components used to trim the unmixing results, which ranges from 0 to 2.
  ###            Peaks with a second-order derivative below zero and a resolution smaller
  ###            than mrsl will be deleted.
  ###
  ###       ncr: The number of cores to be used during the parallel unmixing process. If
  ###            the total number of available cores (NumberOfCluster) is smaller than ncr, 
  ###            NumberOfCluster will be used instead of ncr.
  ###
  ###    rmZero: A logical value indicating whether removing zeros from the volume percentages.
  ### 
  ###   outfile: A character indicating the name of the PDF/RData file the unmixing results will be written to.  
  ###            The results will be returned to files with a default name if outfile=NULL. 
  ###
  ###     addvl: A logical value indicating whether vertical lines will be added to the unmixed 
  ###            grain-size components when visualizing the results using a PDF file. 
  ###
  ###     logxy: A character indicating whether the x- or y-axis will be logged in the PDF file, 
  ###            one of "", "x", "y", "xy", or NULL.
  ###
  ###       lwd: A numeric value giving the widths of lines in the PDF file.
  ###       
  ###       pch: An integer giving the type of symbols in the PDF file.
  ### 
  ###       cex: A numeric value giving the size of symbols in the PDF file.
  ###==========================================================================================================================
  ### The function
  ###
  ### (1) Returns a invisible list of S3 object of class "batchgsd" containing the unmixing results of individual samples.
  ###
  ### (2) Generates a PDF file showing the unmixed results in the current working directory.
  ###
  ### (3) Generates a RData file containing the unmixed results in the current working directory.
  ###==========================================================================================================================
  ussGSDbatchp <- function(GSDdata="inputGSD.csv", ncomp=NULL, model="weibull0", mpd=0, ctf=0.1, ntry=50,  
                           kkf=0.1, startPars=NULL, alphaRange=NULL, sigmaRange=NULL, omegaRange=NULL, 
                           qRange=NULL, useIndex=TRUE, minfunc="fom", trim=TRUE, mrsl=0.6, ncr=NULL, 
                           rmZero=TRUE, outfile=NULL, addvl=TRUE, logxy="x", lwd=3, pch=21, cex=2) {

      stopifnot(is.null(ncomp) || (length(ncomp)==1 && is.numeric(ncomp) && ncomp %in% 1:13),
                is.null(model) || (length(model)==1 && model %in% c("weibull","lognormal","weibull0", "lognormal0", "skewnormal0", "skewgnormal0")),
                is.numeric(mpd), length(mpd)==1, mpd>=0, mpd<3,
                is.numeric(ctf), length(ctf)==1, ctf>=-1, ctf<=1,
                is.numeric(ntry), length(ntry)==1, ntry>=2, 
                is.numeric(kkf), length(kkf)==1, kkf>0, kkf<1,
                is.null(startPars) || (is.matrix(startPars) && ncol(startPars)==2), 
                is.null(alphaRange) || (is.numeric(alphaRange) && length(alphaRange)==2),
                is.null(sigmaRange) || (is.numeric(sigmaRange) && length(sigmaRange)==2),
                is.null(omegaRange) || (is.numeric(omegaRange) && length(omegaRange)==2),
                is.null(qRange) || (is.numeric(qRange) && length(qRange)==2),
                is.logical(useIndex), length(useIndex)==1,
                is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "rss"),
                is.logical(trim), length(trim)==1,
                is.numeric(mrsl), length(mrsl)==1, mrsl>=0, mrsl<=2,
                is.null(ncr) || (length(ncr)==1 && is.numeric(ncr)), 
                is.logical(rmZero), length(rmZero)==1,  
                is.null(outfile) || (is.character(outfile) && length(outfile)==1),
                is.logical(addvl), length(addvl)==1, 
                is.character(logxy), length(logxy)==1,
                is.numeric(lwd), length(lwd)==1, 
                is.numeric(pch), length(pch)==1,
                is.numeric(cex), length(cex)==1)

      ###
      if (is.character(GSDdata) && file.exists(GSDdata)==FALSE) {

          TempTable <- data.frame(rbind(c("GSL",exp(log(0.02)+(log(2000)-log(0.02))/100*(0:100))),
          c("GSD1",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0.12,0.18,0.22,0.26,0.29,0.32,0.34,0.35,0.35,0.36,0.36,0.36,
            0.37,0.4,0.43,0.47,0.53,0.58,0.64,0.7,0.75,0.82,0.89,0.96,1.04,1.12,1.2,1.27,1.34,1.4,1.47,1.53,1.61,1.7,1.81,1.96,2.14,
            2.37,2.65,2.97,3.32,3.69,4.06,4.41,4.71,4.92,5.03,5.02,4.86,4.56,4.13,3.6,2.98,2.35,1.7,1.18,0.75,0.06,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD2",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.07,0.09,0.11,0.13,0.14,0.15,0.15,0.15,0.15,0.15,0.15,
            0.15,0.15,0.17,0.19,0.21,0.24,0.28,0.3,0.33,0.34,0.35,0.36,0.36,0.35,0.35,0.35,0.35,0.36,0.37,0.4,0.42,0.45,0.48,
            0.5,0.49,0.47,0.42,0.36,0.31,0.29,0.36,0.54,0.89,1.43,2.2,3.16,4.28,5.47,6.62,7.6,8.3,8.63,8.51,7.95,7.02,5.81,
            4.45,3.13,1.71,0.32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD3",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.06,0.08,0.1,0.12,0.13,0.14,0.14,0.14,0.14,0.13,0.13,
            0.13,0.13,0.14,0.16,0.18,0.2,0.23,0.25,0.27,0.29,0.3,0.31,0.31,0.32,0.32,0.33,0.35,0.37,0.4,0.43,0.47,0.5,0.51,
            0.51,0.47,0.41,0.33,0.25,0.21,0.24,0.4,0.73,1.28,2.04,3.04,4.18,5.43,6.63,7.68,8.43,8.8,8.71,8.18,7.24,6.03,
            4.66,3.31,1.97,0.62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD4",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.09,0.13,0.16,0.19,0.22,0.24,0.25,0.26,0.27,0.27,0.27,0.28,
            0.29,0.31,0.35,0.39,0.43,0.48,0.53,0.57,0.62,0.66,0.71,0.76,0.82,0.88,0.94,1.01,1.07,1.14,1.2,1.27,1.34,1.42,1.5,
            1.59,1.7,1.83,2,2.21,2.46,2.78,3.15,3.58,4.03,4.49,4.91,5.25,5.47,5.52,5.38,5.05,4.53,3.88,3.14,2.39,1.67,1.15,
            0.39,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD5",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.09,0.13,0.16,0.19,0.21,0.23,0.24,0.25,0.25,0.25,0.25,0.25,0.26,
            0.28,0.3,0.34,0.37,0.41,0.45,0.48,0.52,0.56,0.61,0.66,0.72,0.79,0.85,0.93,0.99,1.06,1.12,1.17,1.22,1.26,1.3,1.36,1.44,
            1.56,1.75,2,2.34,2.77,3.28,3.86,4.46,5.05,5.56,5.94,6.13,6.1,5.84,5.36,4.69,3.9,3.07,2.27,1.46,0.57,0.05,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD6",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.07,0.09,0.11,0.13,0.14,0.15,0.16,0.16,0.16,0.15,0.15,0.15,0.15,0.17,
            0.19,0.21,0.24,0.27,0.3,0.32,0.34,0.36,0.37,0.38,0.4,0.42,0.45,0.48,0.53,0.58,0.63,0.68,0.72,0.73,0.71,0.65,0.57,0.48,0.4,
            0.39,0.5,0.76,1.22,1.92,2.82,3.91,5.07,6.23,7.23,7.96,8.31,8.24,7.72,6.86,5.72,4.47,3.23,2.11,1.24,0.54,0.15,0.02,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD7",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.03,0.07,0.1,0.12,0.15,0.17,0.18,0.19,0.2,0.2,0.2,0.2,0.2,0.22,0.24,0.26,0.3,
            0.34,0.39,0.43,0.47,0.51,0.55,0.59,0.63,0.67,0.71,0.76,0.81,0.87,0.94,1.02,1.11,1.2,1.29,1.38,1.45,1.5,1.54,1.56,1.59,1.66,
            1.78,2,2.33,2.8,3.4,4.08,4.81,5.49,6.05,6.4,6.48,6.24,5.73,4.96,4.04,3.06,2.1,1.06,0.2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD8",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.1,0.14,0.17,0.2,0.23,0.25,0.26,0.27,0.28,0.28,0.28,0.29,0.31,0.33,0.37,0.42,
            0.48,0.54,0.6,0.66,0.72,0.78,0.84,0.9,0.97,1.03,1.1,1.16,1.23,1.29,1.37,1.44,1.53,1.61,1.71,1.8,1.9,2.01,2.12,2.25,2.41,2.61,2.85,
            3.14,3.48,3.87,4.26,4.63,4.94,5.12,5.14,4.98,4.61,4.08,3.41,2.68,1.93,1.31,0.31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          c("GSD9",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04,0.1,0.15,0.19,0.22,0.25,0.28,0.29,0.3,0.31,0.32,0.32,0.33,0.35,0.37,0.41,0.45,
            0.5,0.55,0.6,0.65,0.7,0.76,0.82,0.88,0.95,1.03,1.1,1.17,1.23,1.3,1.37,1.45,1.54,1.66,1.82,2.01,2.26,2.55,2.89,3.26,3.66,4.06,4.43,
            4.75,4.99,5.13,5.14,5.03,4.78,4.42,3.94,3.39,2.79,2.19,1.6,1.11,0.71,0.15,0.02,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)), 
          stringsAsFactors=FALSE)
          
          ###
          colnames(TempTable) <- NULL
          write.csv(TempTable, file=GSDdata, row.names=FALSE)
          cat(paste("An editable GSD input template [", GSDdata, "] has been created at: ", paste(getwd(),"/", GSDdata,"\n",sep=""),sep=""))
             
      } else  {

          if (!is.null(outfile) && !is.character(outfile))  stop("Error: [outfile] should be NULL or a character!")

          ###
          if (!is.null(outfile))  xoxoxo <- pdfName <- rdName <- paste(outfile,"_batch",sep="")

          ###
          if (is.character(GSDdata)) {

              if (length(GSDdata)!=1)  stop("Error: file name [GSDdata] should be of length one!")

              ###
              if (is.null(outfile))  { 

                  nxo <- nchar(GSDdata)
                  xoxoxo <- pdfName <- rdName <- paste(substr(GSDdata, start=1, stop=nxo-4), "_batch",sep="")
                   
              } # end if.

              ###
              GSDdata <- read.csv(GSDdata, header=FALSE)
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!")

              ###
              gsl <- as.numeric(GSDdata[1,-1])
              GSD <- GSDdata[-1,-1,drop=FALSE]
          
              ###
              sampleName <- as.character(GSDdata[-1,1])

         } else {

              if (!is.matrix(GSDdata) && !is.data.frame(GSDdata))  stop("Error: [GSDdata] should be matrix or data.frame!")
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!")

              ###
              gsl <- as.numeric(GSDdata[1,])
              GSD <- GSDdata[-1,,drop=FALSE] 
      
              ###
              if (is.null(outfile))  { 

                  mthCALL <- (as.character(match.call()))[2L]
                  xoxoxo <- pdfName <- rdName <- paste(mthCALL,"_batch",sep="")

              } # end if.

              ###
              sampleName <- paste("GSD",1:nrow(GSD),sep="")
          
          } # end if.

          ###
          if(any(diff(gsl)<=0))  stop("Error: grain-size levels in the first row of [GSDdata] are not of increasing order!")
      
          ###
          N <- nrow(GSD)

          ###
          if (N<9)  stop("Error: Parallel optimisation is forbidden for a sample size of N<9!")
           
          ###
          NumberOfCluster <- parallel::detectCores()

          ###
          if (!is.null(ncr)) {

              if (!ncr %in% (2:1000))  stop("Error: [ncr] must be an integer ranging from 2 to 1000!")

              ###
              if (ncr<=NumberOfCluster) {

                  NumberOfCluster <- ncr

              } else {

                  cat("Note: [ncr=",ncr,"] exceeds the number of available cores (",NumberOfCluster,")!\n",sep="")

              } # end if.

          } # end if.
               
          ###
          cl <- snow::makeCluster(NumberOfCluster, outfile="")
          doSNOW::registerDoSNOW(cl)

          ###
          cat("Parallel unmixing of N=", N, " GSDs in a batch pattern is in progress, please wait, ...\n", sep="")

          ###
          pb <- txtProgressBar(max=N, style=3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress=progress)

          ### Parallel unmixing using the R function ussGSD(). 
          ###-------------------------------------------------------------------------------------
          GSDbatch <- foreach::foreach(kk=1:N, .inorder=TRUE, .options.snow=opts, .packages=c("pracma","minpack.lm"), 
                                       .export=c("ussGSD","ussGSD0"), .verbose=FALSE) %dopar% {   

              gsd <- as.numeric(GSD[kk,])

              ###
              if (!is.null(model)) {

                  if (model %in% c("weibull","lognormal")) {

                      res_ussGSD <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=TRUE, model=model, 
                        mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                        sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                        viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                        rmZero=rmZero, plot=FALSE), silent=TRUE)

                  } else if (model %in% c("weibull0", "lognormal0", "skewnormal0", "skewgnormal0")) {

                      res_ussGSD <- try(ussGSD0(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=TRUE, model=model, 
                        mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                        sigmaRange=sigmaRange, omegaRange=omegaRange, qRange=qRange, useIndex=useIndex, 
                        minfunc=minfunc, trim=trim, mrsl=mrsl, viewAutoInis=FALSE, viewLM=FALSE, 
                        viewFit=FALSE, viewPb=FALSE, outfile=NULL, rmZero=rmZero, plot=FALSE), 
                        silent=TRUE)

                  } # end if.

                  ###
                  if (inherits(res_ussGSD,what="try-error")==FALSE) { 

                      return(res_ussGSD)

                  } else {

                      cat("Failed in unmixing the ", kk, "-th sample!\n", sep="")
                      print(paste("<",sampleName[kk],"> ",attr(res_ussGSD,"condition"),sep=""))
                      return(NULL)

                  } # end if. 

              } else {

                  res_ussGSDwb <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=TRUE, model="weibull", 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                    viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                    rmZero=rmZero, plot=FALSE), silent=TRUE)

                  ###
                  res_ussGSDlg <- try(ussGSD(gsl=gsl, gsd=gsd, ncomp=ncomp, auto=TRUE, model="lognormal", 
                    mpd=mpd, ctf=ctf, ntry=ntry, kkf=kkf, startPars=startPars, alphaRange=alphaRange, 
                    sigmaRange=sigmaRange, useIndex=useIndex, minfunc=minfunc, trim=trim, mrsl=mrsl, 
                    viewAutoInis=FALSE, viewLM=FALSE, viewFit=FALSE, viewPb=FALSE, outfile=NULL, 
                    rmZero=rmZero, plot=FALSE), silent=TRUE)

                  ###
                  if (inherits(res_ussGSDwb,what="try-error")==FALSE && inherits(res_ussGSDlg,what="try-error")==FALSE) {

                      if (minfunc=="fom") {

                          minfuncVAL1 <- res_ussGSDwb$FOM
                          minfuncVAL2 <- res_ussGSDlg$FOM

                      } else if (minfunc=="rss") {

                          minfuncVAL1 <- res_ussGSDwb$RSS
                          minfuncVAL2 <- res_ussGSDlg$RSS

                      } # end if.

                      ###
                      if (minfuncVAL1<minfuncVAL2) {

                          return(res_ussGSDwb)

                      } else {

                          return(res_ussGSDlg)

                      } # end if. 

                  } else if (inherits(res_ussGSDwb,what="try-error")==FALSE) {

                      cat("Failed in unmixing the ", kk, "-th sample using Lognormal!\n", sep="")
                      print(paste("<",sampleName[kk],"> ",attr(res_ussGSDlg,"condition"),sep=""))
                      return(res_ussGSDwb)

                  } else if (inherits(res_ussGSDlg,what="try-error")==FALSE) {

                      cat("Failed in unmixing the ", kk, "-th sample using Weibull!\n", sep="")
                      print(paste("<",sampleName[kk],"> ",attr(res_ussGSDwb,"condition"),sep=""))
                      return(res_ussGSDlg)

                  } else {

                      cat("Failed in unmixing the ", kk, "-th sample using Weibull and Lognormal!\n", sep="")
                      print(paste("<",sampleName[kk],"> ",attr(res_ussGSDwb,"condition"),sep=""))
                      print(paste("<",sampleName[kk],"> ",attr(res_ussGSDlg,"condition"),sep=""))
                      return(NULL)

                  } # end if.
              
              } # end if. 

          } # end foreach.

          ###
          cat("\n")
          close(pb)

          ###
          foreach::registerDoSEQ()
          snow::stopCluster(cl)

          ###
          names(GSDbatch) <- sampleName
          class(GSDbatch) <- "batchgsd"
          attr(GSDbatch, "xoxoxo") <- xoxoxo

          ###
          cat("Note: a new environment called [gsdbatch] is created to hold the unmixing results!\n",sep="")
          gsdbatch <<- new.env()

          ###
          assign("GSDbatch", GSDbatch, envir=gsdbatch)

          ###
          plot_ussGSDbatch(obj_batchgsd=GSDbatch, pdfName=pdfName, addvl=addvl, logxy=logxy, lwd=lwd, pch=pch, cex=cex)
          save(GSDbatch, file=paste(rdName,".RData",sep=""))

          ###
          invisible(GSDbatch)

      } # end if.

  } # end function ussGSDbatchp.
  ###==========================================================================================================================
  ###
