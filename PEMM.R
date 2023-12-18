  ###
  ###**************************************************************************************************************************
  ### R functions PEMM() and PEMMp().
  ###**************************************************************************************************************************
  ### Parametric end-member modelling using transformed probability density functions.
  ### Author: Jun Peng, Hunan University of Science and Technology, China. 
  ### Last updated, 2023.12.17.
  ### 
  ### Note that the current version of the program is developed for parametric end-member modelling 
  ### of single-sample grain-size distributions measured using the Malvern Mastersizer-2000/3000
  ### laser particle-size analyzer. So it can not ensure that the program will work equally 
  ### well for datasets measured from other types of particle-size analyzers.
  ### 
  ### Please contact Jun Peng (pengjun10@mails.ucas.ac.cn, or 656971673@qq.com) 
  ### if you have any problem during the application of the program or 
  ### if you want to report any bug encountered during the use of the program.
  ###
  ### Reference: 
  ###  Peng, J., Wang, X.L., Zhao, H., Dong, Z.B, 2023. Single-sample unmixing and parametric end-member modelling of grain-size 
  ###  distributions with transformed probability density functions and their performance comparison using aeolian sediments. 
  ###  Sedimentary Geology, 445, 106328. 
  ###**************************************************************************************************************************

 ###
  ###****************************************************************************************
  ### Install external packages from CRAN.
  ###****************************************************************************************
  PACKAGE <- c("parallel", "foreach", "snow", "doSNOW")
  idx <- which(!PACKAGE %in% installed.packages()[,1])
  if (length(idx)>0)  install.packages(PACKAGE[idx])
  ###****************************************************************************************   

  ###
  ###***********************************************************************************
  ### Users may also need to install and load the following external R package if 
  ### they want to use a parallel version of the R function PEMM(), i.e., PEMMp().
  ###===================================================================================
  library(parallel)             # install.packages("parallel")
  ###
  library(foreach)              # install.packages("foreach")
  ###
  library(snow)                 # install.packages("snow")
  ###
  library(doSNOW)               # install.packages("doSNOW")
  ###=========================================================================================

  ###
  ###**************************************************************************************************************************
  ### Function getPEMM() is used for getting the results of parametric end-member modelling. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###     ncomp: an integer indicating the number of components.
  ###
  ###   outfile: A character indicating the prefix name of the PDF/CSV file the modelling results will  
  ###            be written to. The results will be returned to files with a default prefix name if outfile=NULL.
  ###
  ###        nr: An integer indicating the number of rows generated in the PDF file.
  ###
  ###        nc: An integer indicating the number of columns generated in the PDF file.
  ###==========================================================================================================================
  ### The function returns the following elements:
  ###
  ###      fitGSD: A matrix containing the fitted/predicted grain-size distributions.
  ###
  ###   abundance: A matrix containing the abundances of grain-size distributions.
  ###
  ###   endmember: A matrix containing the optimized parametric endmembers.
  ###
  ###       Angle: A vector containing the minimized angles calculated using the measured and predicted 
  ###              grain-size distributions.
  ###
  ###   meanAngle: A numerical value indicating the mean of the minimized angles.
  ###         
  ### sdmeanAngle: A numerical value indicating the standard deviation of the mean of the 
  ###              minimized angles.
  ###
  ###          R2: A vector containing the squared correlation coefficients calculated using the measured 
  ###              and predicted grain-size distributions.
  ###
  ###      meanR2: A numerical value indicating the mean of the squared correlation coefficients.
  ###
  ###    sdmeanR2: A numerical value indicating the standard deviation of the mean of the squared 
  ###              correlation coefficients.
  ###
  ###         FOM: A vector containing the figure-of-merit value calculated using the measured 
  ###              and predicted grain-size distributions.
  ###
  ###     meanFOM: A numerical value indicating the mean FOM values.
  ###
  ###   sdmeanFOM: A numerical value indicating the standard deviation of the mean of the FOM values.
  ###
  ###    fit.pars: A matrix containing the optimized parameters describing parametric endmembers.
  ###
  ###     gs.pars: A matrix containing the optimized parameters of the grain-size distrubtion, 
  ###              including the means, meadians, modes, and normalized abundance.
  ###==========================================================================================================================
  getPEMM <- function(obj_pemm, ncomp=NULL, outfile=NULL, nr=10, nc=9)  {

      ###
      stopifnot(inherits(obj_pemm, what="pemm")==TRUE,
                is.null(ncomp) || (length(ncomp)==1 && is.numeric(ncomp) && ncomp>=1), 
                is.null(outfile) || (is.character(outfile) && length(outfile)==1),
                is.numeric(nr), length(nr)==1, nr>0, is.numeric(nc), length(nc)==1, nc>0)
    
      ###
      maxcomp <- length(obj_pemm)
      idx <- sapply(obj_pemm, is.list)
      meanAngles <- sapply(obj_pemm[idx], function(x)  x$meanAngle)
      sdmeanAngles <- sapply(obj_pemm[idx], function(x)  x$sdmeanAngle)
      meanR2s <- sapply(obj_pemm[idx], function(x)  x$meanR2)
      sdmeanR2s <- sapply(obj_pemm[idx], function(x)  x$sdmeanR2)
      meanFOMs <- sapply(obj_pemm[idx], function(x)  x$meanFOM)
      sdmeanFOMs <- sapply(obj_pemm[idx], function(x)  x$sdmeanFOM)

      ### 
      if (is.null(ncomp))  {

          par("mar"=c(5.1,4.1,4.1,4.1))
          plot((1:maxcomp)[idx], meanAngles, type="o", yaxt="n", col="skyblue", cex=3, lwd=2, xlab="Number of end-members [k]", 
               ylab="Mean angle", cex.lab=1.2, main="Variation of PEMM statistics with [k]", cex.main=1.5,mgp=c(2,1,0))
          axis(side=2, col="skyblue", col.ticks="skyblue", lwd=4)
          suppressWarnings(arrows(x0=(1:maxcomp)[idx], y0=meanAngles[idx]-sdmeanAngles[idx]/2.0, 
              x1=(1:maxcomp)[idx], y1=meanAngles[idx]+sdmeanAngles[idx]/2.0,code=3, lwd=1.6, angle=90, length=0.05, col="skyblue"))

          ###
          par("new"=TRUE)
          plot((1:maxcomp)[idx], meanR2s, type="o", xaxt="n", yaxt="n", col="red", xlab="", ylab="", cex=3, lwd=2)
          axis(side=4, col="red",col.ticks="red",lwd=3)
          mtext("Mean R2", side=4, line=2, cex=1.2)
          suppressWarnings(arrows(x0=(1:maxcomp)[idx], y0=meanR2s[idx]-sdmeanR2s[idx]/2.0, 
              x1=(1:maxcomp)[idx], y1=meanR2s[idx]-sdmeanR2s[idx]/2.0, code=3, lwd=1.6, angle=90, length=0.05, col="red"))

          ###
          k <- as.numeric(readline(paste("Please enter the optimal number of end-member (i.e., [k]): ")))
          if (!is.finite(k))  stop("[ncomp] has not been provided!")
          abline(v=k, lwd=3, lty=2)

      } else {
 
          k <- ncomp

      } # end if.

      ###
      if (k>maxcomp)  stop(paste("[ncomp] should not exceed maxcomp=", maxcomp,"!",sep="")) 

      ###
      if (!is.list(obj_pemm[[k]]))  stop(paste("PEMM result with ", k, " end-members is not available!", sep=""))

      ###
      model <- attr(obj_pemm, "model")
      gsl <- attr(obj_pemm, "gsl")
      GSD <- attr(obj_pemm, "GSD")
      sampleName <- attr(obj_pemm, "sampleName")
      
      ###
      fitGSD <- obj_pemm[[k]]$fitGSD
      abundance <- obj_pemm[[k]]$abundance
      endmember <- obj_pemm[[k]]$endmember
      Angle <- obj_pemm[[k]]$Angle
      names(Angle) <- sampleName
      meanAngle <- obj_pemm[[k]]$meanAngle
      sdmeanAngle <- obj_pemm[[k]]$sdmeanAngle
      R2 <- obj_pemm[[k]]$R2
      names(R2) <- sampleName
      meanR2 <- obj_pemm[[k]]$meanR2
      sdmeanR2 <- obj_pemm[[k]]$sdmeanR2
      FOM <- obj_pemm[[k]]$FOM
      names(FOM) <- sampleName
      meanFOM <- obj_pemm[[k]]$meanFOM
      sdmeanFOM <- obj_pemm[[k]]$sdmeanFOM
      fit.pars <- obj_pemm[[k]]$fit.pars
      gs.pars <- obj_pemm[[k]]$gs.pars
      
      ###
      m <- nrow(GSD)
      tabundance <- t(abundance)
      compSig <- vector(length=m, mode="list")
      xm <- gs.pars[,3]
      abd <- gs.pars[,4]
     
      ###
      if (is.null(outfile))  outfile <- attr(obj_pemm, "xoxoxo") 

      ###
      for (i in 1:m)   compSig[[i]] <- t(tabundance[,i]*endmember)

      ###
      FL1 <- cbind(abundance, Angle, R2, FOM)
      colnames(FL1) <- c(paste("EM",1:k,sep=""),"Angle","R2","FOM")
      rownames(FL1) <- sampleName
      myfilename1 <- paste(outfile,"_",paste("abundance_",k,model,".csv",sep=""),sep="")
      WF1 <- try(write.csv(FL1, file=myfilename1), silent=TRUE)
      if (inherits(WF1,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename1, ", this may because the file is already opened!",sep=""))

      ###
      FL2 <- cbind(gsl,t(endmember),t(fitGSD))
      colnames(FL2) <- c("GSL",paste("EM",1:k,sep=""),sampleName)
      myfilename2 <- paste(outfile,"_",paste("endmember_fitGSD_",k,model,".csv",sep=""),sep="")
      WF2 <- try(write.csv(FL2, file=myfilename2), silent=TRUE)
      if (inherits(WF2,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename2, ", this may because the file is already opened!",sep=""))

      ###
      colvec <- c("skyblue", "red", "blue", "plum4", "green", "brown", "yellow2", 
                  "deeppink", "black", "wheat", "grey60", "turquoise", "magenta1")

      ###
      myTK <- c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000)
      myLB <- c("0.1","0.2","0.5","1","2","5","10","20","50","100","200","500","1000","2000","5000")

      ###
      myfilename3 <- paste(outfile,"_",k,model,".pdf",sep="")
      WF3 <- try(pdf(myfilename3), silent=TRUE)
      if (inherits(WF3,what="try-error")==TRUE)  stop(paste("Failed in write to ", myfilename3, ", this may because the file is already opened!",sep=""))

      ###
      ### The first plot. 
      par(mar=c(7.1, 6.1, 6.1, 4.1))
      plot((1:maxcomp)[idx], meanAngles, type="o", yaxt="n", col="skyblue", cex=3, lwd=2, xlab="Number of end-members [k]", 
           ylab="Mean angle", cex.lab=1.5, cex.axis=1.2, main="Variation of PEMM statistics with k", cex.main=1.5,mgp=c(2.3,1,0))
      axis(side=2, col="skyblue", col.ticks="skyblue", lwd=4)
      suppressWarnings(arrows(x0=(1:maxcomp)[idx], y0=meanAngles[idx]-sdmeanAngles[idx]/2.0, x1=(1:maxcomp)[idx], 
          y1=meanAngles[idx]+sdmeanAngles[idx]/2.0, code=3, lwd=1.6, angle=90, length=0.05, col="skyblue"))

      ###
      par("new"=TRUE)
      plot((1:maxcomp)[idx], meanR2s, type="o", xaxt="n", yaxt="n", col="red", xlab="", ylab="", cex=3, lwd=2)
      axis(side=4, col="red",col.ticks="red",lwd=3)
      mtext("Mean R2", side=4, line=2, cex=1.2)
      suppressWarnings(arrows(x0=(1:maxcomp)[idx], y0=meanR2s[idx]-sdmeanR2s[idx]/2.0, x1=(1:maxcomp)[idx],
          y1=meanR2s[idx]-sdmeanR2s[idx]/2.0, code=3, lwd=1.6, angle=90, length=0.05, col="red"))

      ### The second plot.
      ###------------------------------------
      par(mfrow=c(1,1))
      par(mar=c(2.1, 2.1, 4.1, 2.1))

      ###
      plot(1,1,type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main="Summary of results", cex.main=2)
   
      ###
      mylegend <- c(paste("Number of single-sample GSDs:  ", m, sep=""), "",paste("End-member model: ", model, sep=""), 
                    paste("Number of end-members: ", k, sep=""), 
                    paste("Mean angle:  ", round(meanAngle,3),"\u00B1",round(sdmeanAngle,3), " \u00B0", sep=""), 
                    paste("Mean R2   :  ", round(meanR2,3), "\u00B1",round(sdmeanR2,3), sep=""),
                    paste("Mean FOM  :  ", round(meanFOM,3), "\u00B1", round(sdmeanFOM,3), " \u0025", sep=""), "",
                    paste("Mode, abundance of EM1:  ", round(xm[1],3), " um, ", round(abd[1],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM2:  ", round(xm[2],3), " um, ", round(abd[2],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM3:  ", round(xm[3],3), " um, ", round(abd[3],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM4:  ", round(xm[4],3), " um, ", round(abd[4],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM5:  ", round(xm[5],3), " um, ", round(abd[5],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM6:  ", round(xm[6],3), " um, ", round(abd[6],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM7:  ", round(xm[7],3), " um, ", round(abd[7],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM8:  ", round(xm[8],3), " um, ", round(abd[8],3), " \u0025", sep=""),
                    paste("Mode, abundance of EM9:  ", round(xm[9],3), " um, ", round(abd[9],3), " \u0025", sep=""))

      ###
      legend("top",legend=mylegend[1:(8+k)], bty="n", cex=1.25)

      ###
      ### The third plot.
      ###------------------------------------
      par(mar=c(7.1, 6.1, 6.1, 4.1))

      ###
      matplot(gsl, t(endmember), type="l", lwd=3, lty=1, col=colvec[1:k], xlab="Grain size (um)", xaxt="n",
              ylab="Probability density", mgp=c(2.3,1,0), log="x", cex.axis=1.2, cex.lab=1.5, cex.main=1.5,
              main="Distribution of end-members")

      ###
      box(lwd=2)
      axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)
      legend("topleft", legend=paste("EM.",1:k,sep=""), col=colvec[1:k],lwd=2, bty="n")

      ###
      ### The fourth plot.
      ###------------------------------------
      layout(matrix(c(1,2,3,4),nrow=2, byrow=TRUE), respect=TRUE)
      par(mar=c(4.1, 4.1, 2.1, 2.1))

      ### 
      matplot(gsl, t(GSD), log="x", type="l", lty=1, xlab="Grain size (um)", ylab="Measured volume (%)", 
              col="red", lwd=1.25, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, xaxt="n")   
      axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)

      ### 
      matplot(gsl, t(fitGSD), log="x", type="l", lty=1, xlab="Grain size (um)", ylab="Fitted volume (%)",
              col="skyblue", lwd=1.25, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25, xaxt="n")
      axis(side=1, at=myTK, labels=myLB, cex.axis=1.2)

      ### 
      plot(1:m, Angle, type="p", lty=1, xlab="Sample index", ylab="Optimized angle (Degree)", 
           col="brown", lwd=3, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25)
      abline(h=meanAngle, col="red", lty=2, lwd=2)

      ### 
      plot(1:m, R2, type="p", lty=1, xlab="Sample index", ylab="Optimized R2", 
           col="purple", lwd=3, mgp=c(2.3,1,0), cex.axis=1.25, cex.lab=1.25)
      abline(h=meanR2, col="red", lty=2, lwd=2)

      ###
      ### The fifth plot.
      ###-------------------------------------
      layout(matrix(seq(nr*nc),nrow=nr, ncol=nc, byrow=TRUE), respect=TRUE)
      par(mar=c(0.3,0.3,0.3,0.3))

      ### 
      for (i in 1:m) {
 
          ###
          plot(gsl, GSD[i,], type="l", lty=1, lwd=2, log="x", xaxt="n", yaxt="n", bty="n", col="black")
          matpoints(gsl, compSig[[i]], type="l", lty=1, lwd=2, col=colvec[1:k])
          points(gsl, fitGSD[i,], type="l", lty=1, lwd=1.6, col="grey70")

          ###
          legend("topleft", legend=c(sampleName[i], paste("Angle=",round(Angle[i],2),sep=""), 
                 paste("R2=",round(R2[i],2),sep=""),paste("FOM=",round(FOM[i],2),sep="")), bty="n",cex=0.5)

          ###
          box(lwd=3)

      } # end for.
               
      ###  
      dev.off()
          
      ### 
      output <- list("Angle"=Angle, "meanAngle"=meanAngle, "sdmeanAngle"=sdmeanAngle, "R2"=R2, 
                     "meanR2"=meanR2, "sdmeanR2"=sdmeanR2, "FOM"=FOM, "meanFOM"=meanFOM, 
                     "sdmeanFOM"=sdmeanFOM, "fit.pars"=fit.pars, "gs.pars"=gs.pars)

      ###
      return(output)
      
  } # end function getPEMM.
  ###
  
  
  ###
  ###**************************************************************************************************************************
  ### Function PEMM() is used for parametric end-member modelling of 
  ### grain-size distributions using transformed probability density functions. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###    GSDmat: A matrix storing the grain-size data used for analysis. The first row is 
  ###            the grain-size levels, and the remaining rows are the volume percentages 
  ###            of individual samples to be modelled. Or a character indicating the name of
  ###            the CSV file containing the grain-size data used for analysis. A template
  ###            of the CSV file with name "inputGSD" can be generated using the command "PEMM()".
  ###
  ###   maxcomp: An integer (from 1 to 8) indicating the maximum number of endmembers to be modelled. 
  ###
  ###     model: A character indicating the model to be fitted,
  ###            "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###
  ###   minfunc: A character indicating the objective to be mimimized, either "fom" or "angle",
  ###            for the mean figure-of-merit values or the mean angles calculated using measured and fitted GSDs.
  ###
  ###  useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###            used as the independent variable during the modelling process.
  ###==========================================================================================================================
  ### The function returns an invisible list of S3 class object "pemm" containing 
  ### the results for different numbers of end-members.
  ###========================================================================================================================== 
  PEMM <- function(GSDdata="inputGSD.csv", maxcomp=6, model="weibull0", minfunc="angle", useIndex=TRUE) {

      ###
      stopifnot(length(maxcomp)==1, is.numeric(maxcomp), maxcomp %in% 1:8,
                is.character(model), length(model)==1, model %in% c("weibull0","lognormal0","skewnormal0","skewgnormal0"),
                is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "angle"),
                is.logical(useIndex), length(useIndex)==1)

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

          if (is.character(GSDdata)) {

              if (length(GSDdata)!=1)  stop("Error: file name [GSDdata] should be of length one!")

              ###
              nxo <- nchar(GSDdata)
              xoxoxo <- substr(GSDdata, start=1, stop=nxo-4)

              ###
              GSDdata0 <- read.csv(GSDdata, header=FALSE)
              GSDdata <- GSDdata0[,-1,drop=FALSE]
              sampleName <- as.character(GSDdata0[-1,1,drop=TRUE])
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!") 

          } else {

              if (!is.matrix(GSDdata) && !is.data.frame(GSDdata))  stop("Error: [GSDdata] should be matrix or data.frame!")
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!")
      
              ###
              mthCALL <- (as.character(match.call()))[2L]
              xoxoxo <- mthCALL
              sampleName <- paste("GSD",1:(nrow(GSDdata)-1),sep="")
          
          } # end if.

          ###
          LIST <- vector(length=maxcomp, mode="list")
          names(LIST) <- paste("k",1:maxcomp,sep="")
          class(LIST) <- "pemm"

          ###
          attr(LIST, "model") <- model
          attr(LIST, "gsl") <- as.numeric(GSDdata[1,])
          attr(LIST, "GSD") <- as.matrix(GSDdata[-1,,drop=FALSE])
          attr(LIST, "sampleName") <- sampleName
          attr(LIST, "xoxoxo") <- xoxoxo
          
          ###
          for (i in 1:maxcomp) {

              myPEMM <- try(pemm(GSDdata, ncomp=i, model=model, minfunc=minfunc, useIndex=useIndex, viewPb=TRUE),silent=TRUE)
   
              ###
              if (inherits(myPEMM, what="try-error")==FALSE) {
                  
                  LIST[[i]] <- myPEMM

              } else {

                  cat("Failed in PEMM optimization when ", i, " end-members were applied!\n", sep="")
                  print(attr(myPEMM,"condition"))
                  LIST[[i]] <- NA
 
              } # end if.

          } # end for.

          invisible(LIST)
          
      } # end if.

  } # end function PEMM.
  ### 
  
  ###
  ###**************************************************************************************************************************
  ### Function pemm() is used for parametric end-member modelling of grain-size distributions
  ### using transformed probability density functions with a specified number of components. 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###    GSDmat: A matrix storing the grain-size data used for analysis. The first row is 
  ###            the grain-size levels, and the remaining rows are the volume percentages 
  ###            of individual samples to be modelled. 
  ###
  ###     ncomp: An integer (from 1 to 8) indicating the number of endmembers to be modelled. 
  ###
  ###     model: A character indicating the model to be fitted,
  ###            "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###
  ###  useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###            used as the independent variable during the modelling process.
  ###
  ###    viewPb: A logical value indicating whether the progress bar will be visualized during calculation. 
  ###==========================================================================================================================
  ### The function returns an invisible list of S3 class object "pemm" containing the following elements.
  ###
  ###      fitGSD: A matrix containing the fitted/predicted grain-size distributions.
  ###
  ###   abundance: A matrix containing the abundances of grain-size distributions.
  ###
  ###   endmember: A matrix containing the optimized parametric endmembers.
  ###
  ###       Angle: A vector containing the minimized angles calculated using the measured and predicted 
  ###              grain-size distributions.
  ###
  ###   meanAngle: A numerical value indicating the mean of the minimized angles.
  ###         
  ### sdmeanAngle: A numerical value indicating the standard deviation of the mean of the 
  ###              minimized angles.
  ###
  ###          R2: A vector containing the squared correlation coefficients calculated using the measured 
  ###              and predicted grain-size distributions.
  ###
  ###      meanR2: A numerical value indicating the mean of the squared correlation coefficients.
  ###
  ###    sdmeanR2: A numerical value indicating the standard deviation of the mean of the squared 
  ###              correlation coefficients.
  ###
  ###         FOM: A vector containing the figure-of-merit value calculated using the measured 
  ###              and predicted grain-size distributions.
  ###
  ###     meanFOM: A numerical value indicating the mean FOM values.
  ###
  ###   sdmeanFOM: A numerical value indicating the standard deviation of the mean of the FOM values.
  ###
  ###    fit.pars: A matrix containing the optimized parameters describing parametric endmembers.
  ###
  ###     gs.pars: A matrix containing the optimized parameters of the grain-size distrubtion, 
  ###              including the means, meadians, and modes.
  ###==========================================================================================================================
  pemm <- function(GSDdata, ncomp, model="weibull0", minfunc="angle", useIndex=TRUE, viewPb=TRUE) {

      if (!is.matrix(GSDdata) && !is.data.frame(GSDdata)) stop("Error: [GSDdata] should be matrix or data.frame!")
      if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!")

      ###
      gsl <- as.numeric(GSDdata[1,])
      GSD <- as.matrix(GSDdata[-1,,drop=FALSE])

      ###
      X <- GSD
      m <- nrow(X)
      n <- ncol(X)

      ###
      xd <- gsl
      
      ### Check if the grain size levels are of log-scale.
      xd0xd0 <- round(diff(log(as.numeric(xd))),1)
      YESORNO <- all(xd0xd0==min(xd0xd0)) || all(xd0xd0 %% min(xd0xd0)==0)
      expGSlev <- ifelse(YESORNO, TRUE, FALSE)

      ###
      if (useIndex==TRUE) {

          xd1 <- seq(n)

          ###
          if(expGSlev==TRUE) {

              if (!all(xd0xd0==xd0xd0[1]))  warning("Logged grain-size levels do not increase with equal steps!")
              
              ###
              pab <- as.numeric(lm(xd1~log(xd))$coefficients)
              
          } else {
          
              xd2xd2 <- round(diff(as.numeric(xd)),2)
              if (!all(xd2xd2==xd2xd2[1]))  warning("Grain-size levels do not increase with equal steps!")
 
              ###
              pab <- as.numeric(lm(xd1~xd)$coefficients)
          
          } # end if.

      } else {

          xd1 <- xd

      } # end if.

      ### 
      ### The error function.
      ###---------------------
      erf <- function(x) { 2 * pnorm(x * sqrt(2)) - 1 }

      ###
      ### The function to be minimized.
      ###-------------------------------------
      minfn <- function(p,xd,X,model,minfunc) {

          X <- as.matrix(X)

          ###
          m <- nrow(X)
          n <- ncol(X)

          ###
          ###
          if (model %in% c("weibull0","lognormal0")) {

              ncomp <- length(p)/2

          } else if (model %in% c("skewnormal0","skewgnormal0")) {

              ncomp <- length(p)/3

          } # end if.

          ###
          B <- matrix(nrow=ncomp, ncol=n)

          ###
          Mode <- Mean <- Median <- vector(length=ncomp)

          ###
          if (model=="lognormal0") {

              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  sigma <- abs(p[i+ncomp])+0.001
                   
                  ###
                  B[i,] <- 1/sqrt(2*pi)/sigma/xd*exp(-0.5*((log(xd/xm)-sigma^2)/sigma)^2)
                  B[i,] <- B[i,]/sum(B[i,])

              } # end for.

          } else if (model=="weibull0") {

              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  alpha <- abs(p[i+ncomp])+1

                  ###
                  B[i,] <- 1/xm*(alpha-1)*(xd/xm)^(alpha-1)*exp(-(alpha-1)/alpha*(xd/xm)^alpha)
                  B[i,] <- B[i,]/sum(B[i,])

              } # end for.

          } else if (model=="skewnormal0") {

              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  alpha <- p[i+ncomp]
                  omega <- abs(p[i+2*ncomp])+0.001 
                  delta <- alpha/sqrt(1+(alpha)^2)

                  ###
                  v1 <- sqrt(1-2*delta^2/pi)*(4-pi)/4*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^1.5

                  ###
                  v2 <- sign(alpha)/2*exp(-2*pi/abs(alpha))

                  ###
                  D <- omega*(sqrt(2/pi)*delta-v1-v2)                

                  ###
                  B[i,] <-  1/sqrt(2*pi)/omega*exp(-0.5*(xd-xm+D)^2/omega^2)*(1+erf(alpha*(xd-xm+D)/sqrt(2)/omega))
                  B[i,] <- B[i,]/sum(B[i,])

              } # end for.

          } else if (model=="skewgnormal0") {

              x0 <- seq(from=0,to=max(xd)/4,by=max(xd)/4/200)
                
              ###
              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  sigma <- abs(p[i+ncomp])+0.001

                  ###
                  qv <- abs(p[i+2*ncomp])+0.001
                  if (qv>1)  qv <- 1

                  ###
                  pv <- 2 + 6*(1-qv)^5

                  ###
                  yf1 <- exp(x0/sigma*qv)
                  yf2 <- exp(x0/sigma/qv)

                  ###
                  y0 <- abs(yf1*qv+yf2/qv)/(yf1+yf2)*exp(-0.5*(abs(log((yf1+yf2)/2)))^pv)
                  if (any(!is.finite(y0)))  return(1e30)  
                        
                  ###            
                  L <- x0[which.max(y0)]

                  ###
                  pf1 <- exp((xd-xm+L)/sigma*qv)
                  pf2 <- exp((xd-xm+L)/sigma/qv)
       
                  ###
                  B[i,] <-  1/(2^(1+1/pv)*sigma*gamma(1+1/pv))*abs(pf1*qv+pf2/qv)/(pf1+pf2)*exp(-0.5*(abs(log(0.5*(pf1+pf2))))^pv)
                  B[i,] <- B[i,]/sum(B[i,])  

              } # end for.

          } # end if.

          ###
          BB <- try(solve(B%*%t(B)),silent=TRUE)
          if (inherits(BB,what="try-error")==TRUE)  return(1e30)

          ###
          M <- X %*% t(B) %*% BB 
          M[M<0] <- 0

          ###
          for (i in 1:m) {

              M[i,] <- M[i,]/sum(M[i,])*100

          } # end for.

          ###
          X1 <- M %*% B

          ###
          if (all(is.finite(X)) && all(is.finite(X1))) {

              ###
              if (minfunc=="angle") {

                  return(mean(abs(acos(rowSums(X*X1)/sqrt(rowSums(X^2))/sqrt(rowSums(X1^2))))*180/pi))

              } else if (minfunc=="fom") {

                  return(mean(rowSums(abs(X-X1))/rowSums(X1)*100))

              } # end if. 

          } else {

              return(1e30)

          } # end if.

      } # end function minfn.

      ### 
      ### The function used for calculating matrix M, B, and X1.
      ###-------------------------------------------------------
      calMBX <- function(p,xd,X,model) {

          X <- as.matrix(X)

          ###
          m <- nrow(X)
          n <- ncol(X)

          ###
          if (model %in% c("weibull0","lognormal0")) {

              ncomp <- length(p)/2

          } else if (model %in% c("skewnormal0","skewgnormal0")) {

              ncomp <- length(p)/3

          } # end if.

          ###
          B <- matrix(nrow=ncomp, ncol=n)

          ###
          Mode <- Mean <- Median <- vector(length=ncomp)

          ###
          if (model=="lognormal0") {

              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  sigma <- abs(p[i+ncomp])+0.001
                   
                  ###
                  B[i,] <- 1/sqrt(2*pi)/sigma/xd*exp(-0.5*((log(xd/xm)-sigma^2)/sigma)^2)
                  B[i,] <- B[i,]/sum(B[i,])

                  ###
                  theta <- log(xm)+sigma^2
                  Mean[i] <- exp(theta+0.5*sigma^2)
                  Median[i] <- exp(theta)
                  Mode[i] <- xm

              } # end for.
              
          } else if (model=="weibull0") {

              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  alpha <- abs(p[i+ncomp])+1

                  ###
                  B[i,] <- 1/xm*(alpha-1)*(xd/xm)^(alpha-1)*exp(-(alpha-1)/alpha*(xd/xm)^alpha)
                  B[i,] <- B[i,]/sum(B[i,])

                  ###
                  beta <- xm/((alpha-1)/alpha)^(1/alpha)          
                  Mean[i] <- beta*gamma(1+1/alpha)
                  Median[i] <- beta*(log(2))^(1/alpha) 
                  Mode[i] <- xm

              } # end for. 

          } else if (model=="skewnormal0") {

              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  alpha <- p[i+ncomp]
                  omega <- abs(p[i+2*ncomp])+0.001 
                  delta <- alpha/sqrt(1+(alpha)^2)

                  ###
                  v1 <- sqrt(1-2*delta^2/pi)*(4-pi)/4*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^1.5
                  v2 <- sign(alpha)/2*exp(-2*pi/abs(alpha))
                  D <- omega*(sqrt(2/pi)*delta-v1-v2)                

                  ###
                  B[i,] <-  1/sqrt(2*pi)/omega*exp(-0.5*(xd-xm+D)^2/omega^2)*(1+erf(alpha*(xd-xm+D)/sqrt(2)/omega))
                  B[i,] <- B[i,]/sum(B[i,])

                  ###
                  Mean[i] <- xm-D+omega*delta*sqrt(2/pi)
                  Median[i] <- xm
                  Mode[i] <- xm

              } # end for.

          } else if (model=="skewgnormal0") {

              ###
              x0 <- seq(from=0,to=max(xd)/4,by=max(xd)/4/200)
              
              ###
              for (i in 1:ncomp) {

                  xm <- abs(p[i])
                  sigma <- abs(p[i+ncomp])+0.001

                  ###
                  qv <- abs(p[i+2*ncomp])+0.001
                  if (qv>1)  qv <- 1

                  ###
                  pv <- 2 + 6*(1-qv)^5

                  ###
                  yf1 <- exp(x0/sigma*qv)
                  yf2 <- exp(x0/sigma/qv)

                  ###
                  y0 <- abs(yf1*qv+yf2/qv)/(yf1+yf2)*exp(-0.5*(abs(log((yf1+yf2)/2)))^pv)
                       
                  ###             
                  L <- x0[which.max(y0)]

                  ###
                  pf1 <- exp((xd-xm+L)/sigma*qv)
                  pf2 <- exp((xd-xm+L)/sigma/qv)
       
                  ###
                  B[i,] <-  1/(2^(1+1/pv)*sigma*gamma(1+1/pv))*abs(pf1*qv+pf2/qv)/(pf1+pf2)*exp(-0.5*(abs(log(0.5*(pf1+pf2))))^pv)
                  B[i,] <- B[i,]/sum(B[i,])

                  ###
                  kurt <- 2-pv
                  skew <- -6*sign(qv)*(1-qv)^2*(1+1.856*kurt)
                  Mean[i] <- xm-L+skew/6*(1+0.856*kurt)
                  Median[i] <- xm  
                  Mode[i] <- xm      

              } # end for.

          } # end if.

          ###
          M <- X %*% t(B) %*% solve(B%*%t(B)) 
          M[M<0] <- 0

          ###
          for (i in 1:m) {

              M[i,] <- M[i,]/sum(M[i,])*100

          } # end for.

          ###
          X1 <- M %*% B

          ###
          Angle <- R2 <- FOM <- vector(length=m)

          ###
          for (i in 1:m) {

              Angle[i] <- abs(acos(sum(X[i,]*X1[i,])/sqrt(sum(X[i,]^2))/sqrt(sum(X1[i,]^2))))*180/pi
              R2[i] <- (sum((X[i,]-mean(X[i,]))*(X1[i,]-mean(X1[i,]))))^2/sum((X[i,]-mean(X[i,]))^2)/sum((X1[i,]-mean(X1[i,]))^2)
              FOM[i] <- sum(abs(X[i,]-X1[i,]))/sum(X1[i,])*100

          } # end for.

          ###
          list("M"=M, "B"=B, "X1"=X1, "Angle"=Angle, "R2"=R2, "FOM"=FOM, "Mode"=Mode, "Mean"=Mean, "Median"=Median)

      } # end function calMBX.

      ###
      ### Function for generating combinations.
      comb_next <- function(n, k, a, done)  {

          if (done==TRUE)  {

              a <- seq(k)
              done <- FALSE

          } else {
          
              if (a[k] < n)  {

                  a[k] <- a[k] + 1L

                  ###
                  output <- list("done"=done, "a"=a) 
                  return(output)

              } # end if.
  
              ###
              for (i in seq(k, 2L, by=-1L))  {
            
                  if (a[i-1L] < n-k+i-1L)  {

                      a[i-1L] <- a[i-1L] + 1L

                      ###
                      for (j in i:k)  a[j] <- a[i-1L] + j - (i-1L)

                      ###
                      output <- list("done"=done, "a"=a) 
                      return(output)
                      
                  } # end if.

              } # end for.

              ###
              done <- TRUE

          } # end if.
       
          ###
          output <- list("done"=done, "a"=a) 
          return(output)

      } # end function comb_next.

      ###
      ncomb <- function(n, k)  prod((n-k+1L):n) / prod(seq(k))

      ###
      meanGSD <- colMeans(GSD)       
      reserveidx <- which(meanGSD>0.05*max(meanGSD)) 
      minamaxb <- range(gsl[reserveidx])

      ###
      if (expGSlev==TRUE) {

          mdgs <- exp(log(minamaxb[1])+(log(minamaxb[2])-log(minamaxb[1]))/7*(0:7))

      } else {

          mdgs <- minamaxb[1]+(minamaxb[2]-minamaxb[1])/7*(0:7)

      } # end if.
          
      ###
      if (useIndex==TRUE) {
      
          if (expGSlev==TRUE) {

              mdgs <- pab[1] + pab[2]*log(mdgs)
              
          } else {
          
              mdgs <- pab[1] + pab[2]*mdgs
          
          } # end if. 

      } # end if.

      ###
      if (model %in% c("lognormal0","weibull0")) {

          p0 <- vector(length=2*ncomp)

      } else if (model %in% c("skewnormal0","skewgnormal0")) {

          p0 <- vector(length=3*ncomp)

      } # end if.

      ###
      ntry <- ncomb(n=length(mdgs), k=ncomp)

      ###
      if (viewPb==TRUE) {

          cat(paste("model=", model, ", ncomp=", ncomp, ", ntry=", ntry, ".\n", sep="")) 
          cat("Parametric end-member modelling is in progress, please wait, ...\n", sep="") 

          ###
          pb <- txtProgressBar(min=1, max=ntry, initial=1, style=3)

      } # end if.

      ###
      minVAL <- 1e30

      ###
      done <- TRUE 
      a <- 1:ncomp 

      ###
      for (i in 1:ntry) {
               
          r_comb <- comb_next(n=length(mdgs), k=ncomp, a=a, done=done)
          a <- r_comb$a
          done <- r_comb$done

          ###
          p0[1:ncomp] <- mdgs[a]

          ###
          if (model=="lognormal0") {

              if (useIndex==TRUE) {

                  p0[(ncomp+1):(2*ncomp)] <- 0.05

              } else {

                  p0[(ncomp+1):(2*ncomp)] <- 0.5

              } # end if.

          } else if (model=="weibull0") {

              if (useIndex==TRUE) {

                  p0[(ncomp+1):(2*ncomp)] <- runif(n=ncomp,min=6,max=20)

              } else {

                  p0[(ncomp+1):(2*ncomp)] <- runif(n=ncomp,min=2,max=5)

              } # end if.

          } else if (model=="skewnormal0") {

              p0[(ncomp+1):(2*ncomp)] <- 0

              ###
              if (useIndex==TRUE) {

                  p0[(2*ncomp+1):(3*ncomp)] <- 6

              } else {

                  p0[(2*ncomp+1):(3*ncomp)] <- 50

              } # end if.

          } else if (model=="skewgnormal0") {

              if (useIndex==TRUE) {

                  p0[(ncomp+1):(2*ncomp)] <- 5

              } else {

                  p0[(ncomp+1):(2*ncomp)] <- 20

              } # end if.

              ###
              p0[(2*ncomp+1):(3*ncomp)] <- 0.7

          } # end if.

          ###
          OPT <- try(nlminb(start=p0, objective=minfn, xd=xd1, X=X, model=model, minfunc=minfunc),silent=TRUE)    
         
          ###
          if (inherits(OPT,what="try-error")==FALSE && OPT$objective<minVAL) {

              minVAL <- OPT$objective
                
              ###
              bestOPT <- OPT

          } # end if.

          ###
          if (viewPb==TRUE)  setTxtProgressBar(pb, i)

      } # end for.

      ###
      if (viewPb==TRUE)  close(pb)

      ###
      if (!exists("bestOPT"))  stop("PEMM optimization with ncomp=",ncomp," failed!")

      ###
      if (model %in% c("lognormal0","weibull0")) {

          fit.pars <- matrix(bestOPT$par,ncol=2)

      } else if (model %in% c("skewnormal0", "skewgnormal0")) {

          fit.pars <- matrix(bestOPT$par,ncol=3)

      } # end if.

      ###
      fit.pars[,1] <- abs(fit.pars[,1])

      ###
      MBX <- calMBX(p=c(fit.pars[order(fit.pars[,1]),,drop=FALSE]),xd1,X,model=model) 

      ###
      abundance <- MBX$M
      endmember <- MBX$B
      fitGSD <- MBX$X1 
      Angle <- MBX$Angle
      R2 <- MBX$R2
      FOM <- MBX$FOM
      Mode <- MBX$Mode
      Mean <- MBX$Mean
      Median <- MBX$Median

      ###
      meanAngle <- mean(Angle)
      sdmeanAngle <- sd(Angle)/sqrt(m)
      meanR2 <- mean(R2)
      sdmeanR2 <- sd(R2)/sqrt(m)
      meanFOM <- mean(FOM)
      sdmeanFOM <- sd(FOM)/sqrt(m)

      ###
      if (model=="lognormal0") {

          fit.pars[,2] <- abs(fit.pars[,2])+0.001
          colnames(fit.pars) <- c("xm0","sigma")

      } else if (model=="weibull0") {

          fit.pars[,2] <- abs(fit.pars[,2])+1
          colnames(fit.pars) <- c("xm0","alpha")

      } else if (model=="skewnormal0") {

          fit.pars[,3] <- abs(fit.pars[,3])+0.001 
          colnames(fit.pars) <- c("xm0", "alpha","omega")

      } else if (model=="skewgnormal0") {

          fit.pars[,2] <- abs(fit.pars[,2])+0.001
          fit.pars[,3] <- abs(fit.pars[,3])+0.001
          if (any(fit.pars[,3]>1)) fit.pars[fit.pars[,3]>1,3] <- 1
          colnames(fit.pars) <- c("xm0", "sigma","q")

      } # end if.

      ###
      fit.pars <- fit.pars[order(fit.pars[,1]),,drop=FALSE]
      rownames(fit.pars) <- paste("EM",1:ncomp,sep="")
     
      ###
      if (useIndex==TRUE) {
      
          if (expGSlev==TRUE) {

              Mode <- exp((Mode-pab[1])/pab[2])
              Mean <- exp((Mean-pab[1])/pab[2])
              Median <- exp((Median-pab[1])/pab[2])
              
          } else {
          
              Mode <- (Mode-pab[1])/pab[2]
              Mean <- (Mean-pab[1])/pab[2]
              Median <- (Median-pab[1])/pab[2]
              
          } # end if. 

      } # end if.

      ###
      gs.pars <- cbind(Mean, Median, Mode, colSums(abundance)/sum(abundance)*100)
      idx <- order(Mode)
      gs.pars <- gs.pars[idx,,drop=FALSE]
      colnames(gs.pars) <- c("Mean", "Median", "Mode", "Abundance")   
      rownames(gs.pars) <- paste("EM", seq(ncomp), sep="")

      ### 
      output <- list("fitGSD"=fitGSD, "abundance"=abundance[,idx,drop=FALSE], "endmember"=endmember[idx,,drop=FALSE], 
                     "Angle"=Angle, "meanAngle"=meanAngle, "sdmeanAngle"=sdmeanAngle, "R2"=R2, "meanR2"=meanR2, 
                     "sdmeanR2"=sdmeanR2, "FOM"=FOM, "meanFOM"=meanFOM, "sdmeanFOM"=sdmeanFOM, "fit.pars"=fit.pars, 
                     "gs.pars"=gs.pars)

      ###
      invisible(output)

  } # end function pemm.
  ###


  ###
  ###**************************************************************************************************************************
  ### Function PEMMp() is a parallel version of the R function PEMM(). 
  ###==========================================================================================================================
  ### The function contains the following arguments.
  ###
  ###    GSDmat: A matrix storing the grain-size data used for analysis. The first row is 
  ###            the grain-size levels, and the remaining rows are the volume percentages 
  ###            of individual samples to be modelled. Or a character indicating the name of
  ###            the CSV file containing the grain-size data used for analysis. A template
  ###            of the CSV file with name "inputGSD" can be generated using the command "PEMMp()".
  ###
  ###   maxcomp: An integer (from 1 to 8) indicating the maximum number of endmembers to be modelled. 
  ###
  ###     model: A character indicating the model to be fitted,
  ###            "weibull0", "lognormal0", "skewnormal0", or "skewgnormal0".
  ###
  ###   minfunc: A character indicating the objective to be mimimized, either "fom" or "angle",
  ###            for the mean figure-of-merit values or the mean angles calculated using measured and fitted GSDs.
  ###
  ###  useIndex: A logical value indicating whether the index of a grain-size level will be  
  ###            used as the independent variable during the modelling process.
  ###
  ###       ncr: The number of cores to be used during the parallel unmixing process. If
  ###            the total number of available cores (NumberOfCluster) is smaller than ncr, 
  ###            NumberOfCluster will be used instead of ncr.
  ###==========================================================================================================================
  ### The function returns an invisible list of S3 class object "pemm" containing 
  ### the results for different numbers of end-members.
  ###========================================================================================================================== 
  PEMMp <- function(GSDdata="inputGSD.csv", maxcomp=6, model="weibull0", minfunc="angle", useIndex=TRUE, ncr=NULL) {

      ###
      stopifnot(length(maxcomp)==1, is.numeric(maxcomp), maxcomp %in% 1:8,
                is.character(model), length(model)==1, model %in% c("weibull0","lognormal0","skewnormal0","skewgnormal0"),
                is.character(minfunc), length(minfunc)==1, minfunc %in% c("fom", "angle"),
                is.logical(useIndex), length(useIndex)==1)

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

          if (is.character(GSDdata)) {

              if (length(GSDdata)!=1)  stop("Error: file name [GSDdata] should be of length one!")

              ###
              nxo <- nchar(GSDdata)
              xoxoxo <- substr(GSDdata, start=1, stop=nxo-4)

              ###
              GSDdata0 <- read.csv(GSDdata, header=FALSE)
              GSDdata <- GSDdata0[,-1,drop=FALSE]
              sampleName <- as.character(GSDdata0[-1,1,drop=TRUE])
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!") 

          } else {

              if (!is.matrix(GSDdata) && !is.data.frame(GSDdata)) stop("Error: [GSDdata] should be matrix or data.frame!")
              if(nrow(GSDdata)<2)  stop("Error: [GSDdata] should contain at least two rows!")
      
              ###
              mthCALL <- (as.character(match.call()))[2L]
              xoxoxo <- mthCALL
              sampleName <- paste("GSD",1:(nrow(GSDdata)-1),sep="")
          
          } # end if.
          
          ###
          NumberOfCluster <- parallel::detectCores()

          ###
          if (!is.null(ncr)) {

              if (!ncr %in% (2:100)) 
              stop("Error: [ncr] must be an integer ranging from 2 to 100!")

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
          cat("Parametric end-member modelling (with 1 to ",maxcomp," end-members) is in progress, please wait, ...\n", sep="") 

          ###
          pb <- txtProgressBar(max=maxcomp, style=3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress=progress)

          ###
          info1 <- info2 <- c()

          ###
          optPEMM <- foreach::foreach(i=1:maxcomp, .inorder=TRUE,  
                     .options.snow=opts, .export="pemm", .verbose=FALSE) %dopar% { 

              myPEMM <- try(pemm(GSDdata, ncomp=i, model=model, minfunc=minfunc, useIndex=useIndex, viewPb=FALSE), silent=TRUE)
   
              ###
              if (inherits(myPEMM, what="try-error")==FALSE) {
                  
                  OPT <- myPEMM

              } else {

                  info1 <- c(info1, paste("Failed in PEMM optimization when ", i, " end-members were applied!", sep=""))
                  info2 <- c(info2, attr(myPEMM,"condition"))
                  OPT <- NA
 
              } # end if.

              ###
              return(OPT)

          } # end foreach.

          ###
          names(optPEMM) <- paste("k",1:maxcomp,sep="")
          class(optPEMM) <- "pemm"

          ###
          attr(optPEMM, "model") <- model
          attr(optPEMM, "gsl") <- as.numeric(GSDdata[1,])
          attr(optPEMM, "GSD") <- as.matrix(GSDdata[-1,,drop=FALSE])
          attr(optPEMM, "sampleName") <- sampleName
          attr(optPEMM, "xoxoxo") <- xoxoxo

          ###
          ninfo <- length(info1)
          if (ninfo>0L) {

              print("******")
              for (i in 1:ninfo) {

                  print(info1[i])
                  print(info2[i])
                  print("******")

              } # end for.

          } # end if.

          ###
          invisible(optPEMM)
          
      } # end if.

  } # end function PEMMp.
  ###
