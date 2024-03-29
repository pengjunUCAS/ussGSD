
    ######
    #######################################################################################################################################
    ################################ TEMPLATE-1.###########################################################################################
    #######################################################################################################################################
    rm(list=ls())
    ### Load the R program ussGSD.
    ###--------------------------------------------------------------------------------------------------
    source("ussGSD.R")

    ###
    ### See all available functions in program "ussGSD".
    ###--------------------------------------------------------------------------------------------------
    ls()[-1]
   
    ###
    ### Import the grain-size dataset used for unmixing.
    ###--------------------------------------------------------------------------------------------------
    DAT <- read.table("GrainSizeData.txt")

    ###
    ### Automatic unmixing of a single-sample grain-size 
    ### distribution using Lognormal distributions.
    ###--------------------------------------------------------------------------------------------------
    USS <- ussGSD(gsl=DAT[1,], gsd=DAT[2,], auto=TRUE, model="lognormal", viewFit=TRUE, viewLM=TRUE)

    ###
    ### Access elements in the R object "USS".
    ###--------------------------------------------------------------------------------------------------
    class(USS)

    ### Results contained in object "USS".
    names(USS)

    ### Optimized grain-size parameters.
    USS$gs.pars

    ### Calculated resolutions between components.
    USS$rsl

    ### Unmixed grain-size components.
    USS$gs.comp


    ### Generate a new plot based on the results in object "USS".
    ###--------------------------------------------------------------------------------------------------
    plot_ussGSD(USS, sampleName="Sample xx1", addvl=FALSE, logxy="xy", pch=25)

    ###
    ### Unmixing of a single-sample grain-size distribution with an interactive pattern,
    ### i.e., by clicking the mouse to select grain-size peaks.
    ###--------------------------------------------------------------------------------------------------
    USS1 <- ussGSD(gsl=DAT[1,], gsd=DAT[3,], ncomp=4, model="lognormal", auto=FALSE)
    
    ### Differences in entropy before and after mixing between adjacent components.
    USS1$deltaH

    ###
    USS2 <- ussGSD(gsl=DAT[1,], gsd=DAT[6,], ncomp=4, model="weibull", auto=FALSE)

    ### Locations of grain-size components. 
    USS2$mdxv
    USS2$mdyv

    ###
    ### Implementing a recursive unmixing protocol with four different GSDs.
    ###--------------------------------------------------------------------------------------------------
    ### (1).
    USSa <- ussGSD(gsl=DAT[1,], gsd=DAT[96,], auto=TRUE, model="weibull", trim=TRUE, mrsl=0.71) 

    ### (2).
    USSb <- ussGSD0(gsl=DAT[1,], gsd=DAT[97,], auto=TRUE, model="weibull0", trim=TRUE, mrsl=0.71) 

    ### (3).
    USSc <- ussGSD(gsl=DAT[1,], gsd=DAT[10,], auto=TRUE, model="weibull", trim=TRUE, mrsl=0.71) 

    ### (4).
    USSd <- ussGSD0(gsl=DAT[1,], gsd=DAT[42,], auto=TRUE, model="weibull0", trim=TRUE, mrsl=0.71)    
    #####################################################################################################
    ### END TEMPLATE-1.





    #####################################################################################################################################
    ################################ TEMPLATE-2.#########################################################################################
    #####################################################################################################################################
    rm(list=ls())
    source("ussGSD.R")

    ###
    ### Read the grain-size dataset used for unmixing.
    ###--------------------------------------------------------------------------------------------------
    GrainSizeData <- read.table("GrainSizeData.txt")

    ###
    ### Parallel unmixing of a number of single-sample grain-size distributions in a batch pattern.
    ###--------------------------------------------------------------------------------------------------
    ### Use function ussGSDbatch() instead if parallel calculation cannot be supported.  
    USSbatch1 <- ussGSDbatchp(GrainSizeData, model=NULL, trim=TRUE, mrsl=0.71) 

    ### See the sample names.
    names(USSbatch1)

    ###
    ### Access the unmixing results of the sixth sample.
    ###--------------------------------------------------------------------------------------------------
    USSbatch1$GSD6

    ###
    ### Access the shape parameters of the sixth sample.
    ###--------------------------------------------------------------------------------------------------
    USSbatch1$GSD6$sp

    ###
    ### Visualise the unmixing results of the 20th sample.
    ###--------------------------------------------------------------------------------------------------
    plot_ussGSD(USSbatch1$GSD20)

    ###
    ### Find the ID numbers of grain-size distributions that have FOM values larger than 4. 
    ###--------------------------------------------------------------------------------------------------
    largeFOMidx <- sapply(USSbatch1, function(x) !is.null(x) && is.finite(x$FOM) && x$FOM>4)
    seq(length(USSbatch1))[largeFOMidx]

    ###
    ### Successively update unmixing results stored in the file "ussGSDbatch1.RData".
    ###--------------------------------------------------------------------------------------------------
    ### (1).
    update_ussGSDbatch(sampleNO=c(12,18,29,30,31,72), ncomp=3, auto=TRUE, model="weibull", ctf=1) 

    ### (2).
    update_ussGSDbatch(sampleNO=c(11,33,51,57,83,89,94), ncomp=4, auto=TRUE, model="weibull", ctf=1)

    ### (3).
    update_ussGSDbatch(sampleNO=41, ncomp=5, auto=TRUE, model="weibull", ctf=1) 

    ### (4).clicking the mouse to select grain-size peaks
    update_ussGSDbatch(sampleNO=65, ncomp=4, auto=FALSE, model="weibull", ctf=1) 

    ### (5).
    update_ussGSDbatch(sampleNO=82, ncomp=4, auto=TRUE, model="lognormal", ctf=1) 

    ### (6).clicking the mouse to select grain-size peaks 
    update_ussGSDbatch(sampleNO=86, ncomp=4, auto=FALSE, model="lognormal", ctf=1) 

    
    ###
    ### Visualise the updated unmixing results.
    ###--------------------------------------------------------------------------------------------------
    plot_ussGSDbatch()

    ###
    ### Summary the final unmixing results.
    ###--------------------------------------------------------------------------------------------------
    summary_ussGSDbatch(nkm=6)
    #####################################################################################################
    ### END TEMPLATE-2.
    


    
    
    #####################################################################################################################################
    ################################ TEMPLATE-3.#########################################################################################
    #####################################################################################################################################
    rm(list=ls())
    source("ussGSD.R")

    ###
    ### Load the unmixed results in file "GrainSizeData_batch.RData".
    ###--------------------------------------------------------------------------------------------------
    load_ussGSDbatch("GrainSizeData_batch.RData")

    ### Update unmixing results of certain GSDs.
    update_ussGSDbatch(sampleNO=c(12,18,29,30,31,72), ncomp=3, auto=TRUE, model="weibull", ctf=1)
 
    ### Visualise the results in file "GrainSizeData_batch.RData".
    ###--------------------------------------------------------------------------------------------------
    plot_ussGSDbatch(logxy="xy", lwd=2, pch=23)

    ### Summarize the results.
    summary_ussGSDbatch(nkm=5)
    #####################################################################################################
    ### END TEMPLATE-3. 
    
    

    
    
    #######################################################################################################################################
    ################################ TEMPLATE-4.###########################################################################################
    #######################################################################################################################################
    rm(list=ls())
    ### Load the R program ussGSD.
    ###--------------------------------------------------------------------------------------------------
    source("ussGSD.R")
   
    ###
    ### Read the grain-size dataset used for unmixing.
    ###--------------------------------------------------------------------------------------------------
    DAT <- read.table("GrainSizeData.txt")

    ###
    ### Automatic unmixing of a single-sample grain-size 
    ### distribution using transformed Skew Normal distributions.
    ###--------------------------------------------------------------------------------------------------
    USS1 <- ussGSD0(gsl=DAT[1,], gsd=DAT[2,], auto=TRUE, model="skewnormal0", viewFit=TRUE)
    
    ### Check the fitting quality.
    USS1$FOM
    
    ###
    ### Automatic unmixing of a single-sample grain-size 
    ### distribution using transformed Skewed Generalized Normal distributions.
    ###--------------------------------------------------------------------------------------------------
    USS2 <- ussGSD0(gsl=DAT[1,], gsd=DAT[2,], auto=TRUE, model="skewgnormal0", viewFit=TRUE)
    
    ### Check the fitting quality.
    USS2$FOM
    
    ###
    ### Parallel unmixing of a number of single-sample grain-size distributions in a batch pattern
    ### using transformed Skew Normal distributions.
    ###--------------------------------------------------------------------------------------------------
    ### Use function ussGSDbatch() instead if parallel calculation cannot be supported.  
    ussGSDbatchp(DAT, model="skewnormal0", ntry=30, trim=TRUE,  mrsl=0.71, outfile="ussGSD")     

    ###Summary the above unmixing results.
    summary_ussGSDbatch(nkm=6)                        
    #####################################################################################################
    ### END TEMPLATE-4.




    
    #################################################################################################################################
    ################################ TEMPLATE-5.#####################################################################################
    #################################################################################################################################
    rm(list=ls())
    source("ussGSD.R")

    ###
    ### Load the unmixing results stored in RData "ussGSDbatch0".
    ### Please download the file and put it into the current R working directory.
    ###--------------------------------------------------------------------------------------------------
    load_ussGSDbatch("ussGSDbatch0.RData")
    
    ###
    ### Visualise the unmixing results of RData "ussGSDbatch0".
    ###--------------------------------------------------------------------------------------------------
    plot_ussGSDbatch()

    ###
    ### Summary the unmixing results of RData "ussGSDbatch0".
    ###--------------------------------------------------------------------------------------------------
    summary_ussGSDbatch(nkm=6)


    ###
    ### Classification of the unmixing results of RData "ussGSDbatch0" according to the 
    ### hierarchical clustering algorithm using the unmixed abundances and modes of grain sizes.
    ###--------------------------------------------------------------------------------------------------

    ###
    ### Argument settings.
    ###-----------------------------------------------------------------------------------
    ### [gst] is a character indicating the grain-size type ("mode", "median", or "mean") 
    ### to be analyzed during the application of the hierarchical clustering algorithm 
    ### in analyzing the unmixed grain-size distributions (components).
    gst <- "mode"

    ### [nhc] is an integer indicating the number of clusters to be generated 
    ### (from 2 to 13) during the application of the hierarchical clustering algorithm 
    ### in analyzing the unmixed grain-size distributions (components). 
    nhc <- 7

    ### [normthd] is a character indicating the method used for data normalization, 
    ###  one of "minmax", "meansd", or "sd".
    normthd <- "minmax"

    ### [ mdmc] is a character indicating the method used for distance matrix calculation
    ### in the function dist(), one of "euclidean", "maximum", "manhattan", "canberra",
    ### "binary" or "minkowski".
    mdmc <- "maximum"

    ### [agm] is a character indicating the agglomeration method to be used in  
    ### the function hclust(), one of "ward.D", "ward.D2", "single", "complete", 
    ### "average", "mcquitty", "median", or "centroid". 
    agm <- "complete"

    ### [comb] is an integer (from 1 to 10) indicating the type of the combination used 
    ### for generating the data used for hierarchical clustering.              
    comb <- 5
     
    ### [nr] is an integer indicating the number of rows generated in the PDF file for 
    ### the unmixing results.                            
    nr <- 11

    ### [nc] is an integer indicating the number of columns generated in the PDF file 
    ### for the unmixing results.
    nc <- 10
    ###-----------------------------------------------------------------------------------
   
    ###
    ### Get the data.
    GSDbatch <- get("GSDbatch", envir=gsdbatch)

    ### Add a label to the name of the data.
    xoxoxo <- attr(GSDbatch, which="xoxoxo")

    ### Return the data length.
    N <- length(GSDbatch)

    ### Find GSDs have successfully unmixed results. 
    ookk <- sapply(GSDbatch, function(x) (!is.null(x)) && (is.finite(x$FOM)))

    ###
    if (length(ookk)==0) stop("Error: no data can be used for analysis!")

    ### 
    GSDbatchx <- GSDbatch[ookk]
    Nx <- length(GSDbatchx)

    ### Return abundances of GSDs.
    abundance <- c(sapply(GSDbatchx, function(x) x$gs.pars[,1]), recursive=TRUE)

    ###
    if (gst=="mode") { idx <- 4 } else if (gst=="median") 
                     { idx <- 3 } else { idx <- 2}

    ### Set the upper limits of x and y axis for plotting.
    xupper <- 2.25*max(c(sapply(GSDbatchx, function(x) x$gs.pars[,idx]), recursive=TRUE))

    ###
    yupper1 <- 1.25*max(abundance)
    yupper2 <- 1.05*max(sapply(GSDbatchx, function(x) max(x$gs.comp[,2])))

    ###
    ### Prepare for the data matrix used for hierarchical clustering.
    ###--------------------------------------------------------------

    ###
    mat <- c()
   
    ### Loop through individual GSDs.
    for (i in 1:Nx) {

        abd <- GSDbatchx[[i]]$gs.pars[,1,drop=TRUE]
        md <- GSDbatchx[[i]]$gs.pars[,idx,drop=TRUE]

        ###
        abdmd <- sum(abd*log(md))

        ###
        domtabd <- if (max(abd)>90) 5 else 
                   if (max(abd)>60) 3 else 
                   if (max(abd)>30) 1 else 0

        ###
        x1 <- log(md[which.min(abd)])
        y1 <- min(abd)
 
        ###
        x2 <- log(md[which.max(abd)])
        y2 <- max(abd)

        ###
        x3 <- log(min(md))
        y3 <- abd[which.min(md)]

        ###
        x4 <- log(max(md))
        y4 <- abd[which.max(md)]

        ###
        if (length(abd)>1) {

            ###
            x11 <- log(md[(order(abd,decreasing=FALSE))[2]])
            y11 <- (sort(abd,decreasing=FALSE))[2]

            ###
            x22 <- log(md[(order(abd,decreasing=TRUE))[2]])
            y22 <- (sort(abd,decreasing=TRUE))[2]

            ###
            x33 <- log((sort(md,decreasing=FALSE))[2])
            y33 <- abd[(order(md,decreasing=FALSE))[2]]

            ###
            x44 <- log((sort(md,decreasing=TRUE))[2])
            y44 <- abd[(order(md,decreasing=TRUE))[2]]
              
            rsdx <- sd(log(md))/mean(log(md))

            ###
            rsdy <- sd(abd)/mean(abd)

        } else {

            ###
            x11 <- x1
            y11 <- y1

            ###
            x22 <- x2
            y22 <- y2

            ###
            x33 <- x3
            y33 <- y3
 
            ###
            x44 <- x4
            y44 <- y4

            ###
            rsdx <- rsdy <- 0
                  
        } # end if. 

        ###
        proxy <- c(length(abd), abdmd, domtabd, 
                   x1, y1, x2, y2, 
                   x3, y3, x4, y4, 
                   rsdx, rsdy, 
                   x11, y11, x22, y22, 
                   x33, y33, x44, y44)
  
        ###
        nameproxy <- c("ncomp","weightGZ","dominance",
                       "GZminABD1","minABD1", "GZmaxABD1","maxABD1", 
                       "minGZ1", "ABDminGZ1", "maxGZ1", "ABDmaxGZ1",
                       "rsdGZ", "rsdABD",
                       "GZminABD2","minABD2", "GZmaxABD2","maxABD2", 
                       "minGZ2", "ABDminGZ2", "maxGZ2", "ABDmaxGZ2")
   
        ###           
        ### Select the type of data combination.
        ###-------------------------------------
        combIDX <- if (comb==1) {

            seq(11)

        } else if (comb==2) {

            seq(21)

        } else if (comb==3) {

            c(1,2,3, seq(from=5,to=21,by=2))

        } else if (comb==4) {

            c(1,2,3, seq(from=5,to=11,by=2))

        } else if (comb==5) {

            seq(from=5,to=21,by=2)

        } else if (comb==6) {

            seq(from=5,to=11,by=2)       

        } else if (comb==7) {

            c(1,2,3, seq(from=4,to=20,by=2))      

        } else if (comb==8) {

            c(1,2,3, seq(from=4,to=10,by=2)) 

        } else if (comb==9) {

            seq(from=4,to=20,by=2)

        } else if (comb==10) {

            seq(from=4,to=10,by=2)

        } # end if.   

        mat <- rbind(mat, proxy[combIDX])

    } # end for.

    ###
    ### Select the method used for data normalization.
    ###-----------------------------------------------
    if (normthd=="minmax") {

        mat1 <- apply(mat, MARGIN=2, function(x) (x-min(x))/(max(x)-min(x)))

    } else if (normthd=="meansd") {
      
        mat1 <- apply(mat, MARGIN=2, function(x) (x-mean(x))/sd(x))

    } else if (normthd=="sd") {

        mat1 <- apply(mat, MARGIN=2, function(x) x/sd(x))

    } # end if.

    ###
    ### Perform the hierarchical clustering algorithm.
    ###-----------------------------------------------
    HC <- hclust(dist(mat1,method=mdmc)^2, method=agm)
    clutIDX <- as.numeric(cutree(tree=HC, k=nhc))

    ###
    ### Output the plots into a PDF file.
    pdf(paste(xoxoxo,"_gspattern.pdf",sep=""))

    ### Colors used for plotting.
    colvec <- c("skyblue", "red", "blue", "plum4", "green", "brown", 
                "yellow2", "deeppink", "black", "wheat", "grey60", 
                "turquoise", "magenta1")

    ###
    ### Ticks of x and y axies used for plotting. 
    myTK1 <- c(0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000)
    myTK2 <- seq(from=0, to=100, by=20) 

    ###
    ### The first plot.
    ###-----------------
    par(mfrow=c(1,1))
    par(mar=c(2.1, 2.1, 4.1, 2.1))
    plot(1,1,type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", 
         main="Pattern of unmixed GSDs", cex.main=2)

    ### Set the legend of the plot.
    mylegend <- c(paste("Total number of analysed single-sample GSDs:  ", N,sep=""), "",
                  paste("Method used for data normalisation:  ", normthd, sep=""),
                  paste("Method used for agglomeration:  ", agm, sep=""), 
                  paste("Method used for distance matrix calculation:  ", mdmc, sep=""), "",
                  paste("Number of GSDs belonging to type 1 :  ", sum(clutIDX==1), sep=""),
                  paste("Number of GSDs belonging to type 2 :  ", sum(clutIDX==2), sep=""),
                  paste("Number of GSDs belonging to type 3 :  ", sum(clutIDX==3), sep=""),
                  paste("Number of GSDs belonging to type 4 :  ", sum(clutIDX==4), sep=""),
                  paste("Number of GSDs belonging to type 5 :  ", sum(clutIDX==5), sep=""),
                  paste("Number of GSDs belonging to type 6 :  ", sum(clutIDX==6), sep=""),
                  paste("Number of GSDs belonging to type 7 :  ", sum(clutIDX==7), sep=""),
                  paste("Number of GSDs belonging to type 8 :  ", sum(clutIDX==8), sep=""),
                  paste("Number of GSDs belonging to type 9 :  ", sum(clutIDX==9), sep=""),
                  paste("Number of GSDs belonging to type 10:  ", sum(clutIDX==10), sep=""),
                  paste("Number of GSDs belonging to type 11:  ", sum(clutIDX==11), sep=""),
                  paste("Number of GSDs belonging to type 12:  ", sum(clutIDX==12), sep=""),
                  paste("Number of GSDs belonging to type 13:  ", sum(clutIDX==13), sep=""))
    legend("top",legend=mylegend[1:(nhc+6)], bty="n", cex=1.25)

    ###
    ### The second plot.
    ###-----------------
    layout(matrix(seq(nr*nc),nrow=nr, ncol=nc, byrow=TRUE), respect=TRUE)
    par(mar=c(0.3,0.3,0.3,0.3))

    ### Loop through each cluster.
    for (i in 1:nhc) {

        IDX0 <- which(clutIDX==i)
        ndx <- length(IDX0)

        ###
        for (j in 1:ndx) {

            jthidx <- IDX0[j]
 
            ###
            plot(GSDbatchx[[jthidx]]$gs.pars[,c(idx,1),drop=FALSE], type="h", lwd=5, 
                 log="x", xaxt="n", yaxt="n", xlim=c(0.2,xupper), ylim=c(0,yupper1), 
                 col=colvec[i])

            ###
            lines(GSDbatchx[[jthidx]]$gs.pars[,c(idx,1),drop=FALSE], lwd=3, col="purple")
            axis(side=1, at=myTK1, labels=FALSE, cex.axis=2)
            axis(side=2, at=myTK2, labels=FALSE, cex.axis=2)

            ###
            box(lwd=3, col=colvec[i])
            text(x=1, y=80, labels=jthidx, cex=1.25, col=colvec[i])

            ###
            par("new"=TRUE)
            plot(GSDbatchx[[jthidx]]$gs.comp[,c(1,2),drop=FALSE], type="l",  
                 lwd=1, log="x", xaxt="n", yaxt="n", bty="n",  
                 xlim=c(0.2,xupper), ylim=c(0,yupper2), col="grey30")

        } # end for.

    } # end for.
    
    ###
    dev.off()

    ###  
    ### Output the results using a CSV file. 
    ###-------------------------------------------------------------------------------------- 
    ### [gspattern] is a matrix showing the results and patterns of grain-size distributions, 
    ### with column names of NO-sample number, ncomp-number of components, weightGZ-weighted 
    ### grain size, dominance-degree of dominance, GZminABD1-GZ at min abundance, minABD1-min 
    ### abundance, GZmaxABD1-GZ at max abundance, maxABD1-max abundance, minGZ1-min GZ,  
    ### maxGZ1-max GZ, ABDminGZ1-abundance at min GZ, ABDmaxGZ1-abundance at max GZ,  
    ### rsdGZ-relative standard deviation of GZ, rsdABD-relative standard deviation of 
    ### abundance, GZminABD2-GZ at second min abundance, minABD2-second min abundance,     
    ### GZmaxABD2-GZ at second max abundance, maxABD2-second max abundance, minGZ2-second 
    ### min GZ, ABDminGZ2-abundance at second min GZ, maxGZ2-second max GZ, ABDmaxGZ2-abundance 
    ### at second max GZ. Where GZ denotes the median, mode, or mean of individual grain size           
    ### component, depending on argument [gst]. 

    ###
    colnames(mat) <- nameproxy[combIDX]
    gspattern <- cbind("NO"=(seq(N))[ookk], mat, "type"=clutIDX)

    ###                                           
    write.csv(gspattern, file=paste(xoxoxo,"_gspattern.csv",sep=""))
    #####################################################################################################
    ### END TEMPLATE-5.
