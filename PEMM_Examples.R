
    
    ### Load the functions.
    source("PEMM.R")

    ### Load the GSD dataset.
    GrainSizeData <- read.table("GrainSizeData.txt") 

    ### Perform parametric end-member modelling with the Skew Normal distribution,
    ### using the R function PEMM(). You may use function PEMMp() instead.
    RES_PEMM <- PEMMp(GrainSizeData, model="skewnormal0", useIndex=TRUE)

 
    ### Names of elements stored in object "RES_PEMM".
    names(RES_PEMM)

    ### Names of elements stored in "RES_PEMM$k1".
    names(RES_PEMM$k1)

    ### Mean angles for different numbers of end-members.
    sapply(RES_PEMM, function(x) x$meanAngle)

    ### Mean R2 values for different numbers of end-members.
    sapply(RES_PEMM, function(x) x$meanR2)

    ### Mean FOM values for different numbers of end-members.
    sapply(RES_PEMM, function(x) x$meanFOM)

    ### Extract the results of 5 end-members.
    getPEMM(RES_PEMM, ncomp=5)
