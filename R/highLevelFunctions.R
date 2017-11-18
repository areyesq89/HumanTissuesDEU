
fitAllExons <- function( countMatrix, modelFrame, mm, shrink, priorsd, testableVector, dispThis,
                          dispOthers, verbose=FALSE) {
    cat("preparing objects\n")
    sf <- modelFrame$sizeFactor
    rows <- rownames(countMatrix)
    cat("start fitting\n")
    allCoefs <- lapply( rows,
                       function(i, countMatrix, mm,
                                modelFrame, shrink, sf, dispThis,
                                dispOthers, testableVector,
                                priorsd ) {
                           x <- names( testableVector ) %in% i
                           y <- rownames( countMatrix ) %in% i
                           if( verbose ){
                               cat( sprintf( "fitting %s\n", i) ) }
                           if( which(x) %% 1000 == 0 ){
                               cat( sprintf( "fitting %s\n", i ) ) }
                           if( testableVector[x] ) {
                               cnts <- as.vector( countMatrix[y,] )
                               disps <- pmax( 3e-3,
                                             ifelse( modelFrame$exon=="this",
                                                    dispThis[x], dispOthers[x] ) )
                               beta0 <- rep( 2, ncol(mm) )
                               ans <- try( shrinkageFit( mm, cnts, disps,
                                                        sf, beta0, shrink, priorsd ) )
                               if( !inherits( ans, "try-error" ) ){
                                   ans
                               }else{
                                   rep( NA, ncol(mm) )
                               }
                           }else{
                               rep( NA, ncol(mm) )
                           }
                       },
                       countMatrix=countMatrix, mm=mm, modelFrame=modelFrame,
                       shrink=shrink, sf=sf, dispThis=dispThis, dispOthers=dispOthers,
                       testableVector=testableVector, priorsd=priorsd )
    allCoefs <- do.call( rbind, allCoefs )
    rownames(allCoefs) <- rownames(countMatrix)
    colnames(allCoefs) <- colnames(mm)
    allCoefs
}

#' @title Fit Relative Exon Usage Coefficients (REUCs)
#' @description
#' This function fits relative exon usage coefficients
#' for each exonic region of a \code{DEXSeqDataSet} object.
#'
#' @param dxd A \code{DEXSeqDataSet} object
#' @param dispThis A numeric vector of dispersion estimates
#' @param dispOthers A numeric vector of dispersion estimates
#' @param priorsd Prior estimates
#' @param bjp A BiocParam instance
#' 
#' @return A vector of relative exon usage coefficients
#' @export
#' @import DEXSeq
#' @import statmod
fitAllExonsParallel <- function(dxd, dispThis, dispOthers, priorsd, bjp){
    indexes1 <- getIndexList( dxd )
    mm <- modelMatrixREUC( dxd, indexes1 )  
    shrink <- 1:ncol(mm) %in% c( indexes1[["sexIdx"]], indexes1[["crossIdx"]]$col )  
    modelFrame <- colData(dxd)
    sizeFactors <- modelFrame$sizeFactor
    beta0 <- rep(2, ncol(mm))
    modelFrame <- colData(dxd)
    countMatrix <- counts(dxd)
    n <- bpnworkers(bjp)
    toSplit <- sort(rep(1:n, length.out=length(rownames(countMatrix))))
    matList <- lapply( split( rownames(countMatrix), toSplit ), function(x){
        countMatrix[x,]
    })
    cat(sprintf("estimating coefs with %s cores\n", n))
    allCoefsWeakShrinkage <-
        bplapply( matList,
            function(x, modelFrame, mm, shrink, priorsd,
                     testableVector, dispThis, dispOthers, verbose){
                fitAllExons( x, modelFrame, mm, shrink, priorsd,
                            testableVector, dispThis, dispOthers, verbose )
               }, modelFrame=modelFrame, mm=mm, shrink=shrink,
                 priorsd=priorsd, testableVector=testableVector,
                 dispThis=dispThis, dispOthers=dispOthers,
               verbose=TRUE, BPPARAM=bjp )
    do.call(rbind, allCoefsWeakShrinkage)
}

estimateDispersionsForMat <- function( countMatrix, coefsMat, testableVector,
                                      mm, sizeFactors, isThis){
    dispDf <- lapply( rownames(countMatrix), function(i){
        cat(sprintf("fitting %s\n", i))
        x <- rownames(countMatrix) %in% i
        xi <- rownames(coefsMat) %in% i
        stopifnot(all(rownames(coefsMat) == names(testableVector)))
        if( !testableVector[xi] ){
            return(c( NA, NA ) )
        }
        muhat <- as.vector(exp( mm %*% coefsMat[xi,] ) * sizeFactors)
        y <- as.vector(countMatrix[x,])
        disps <- try( estimateDispForExon( muhat, y, mm, isThis ), silent=TRUE )
        if( !inherits( disps, "try-error" ) ){
            disps
        }else{
            c( NA, NA )
        }
    })
    dispDf <- do.call(rbind, dispDf)
    rownames(dispDf) <- rownames(countMatrix)
    dispDf
}

#' @title Fit dispersion estimates to estimate REUCs
#' @description
#' Function to estimate dispersions for each exonic region
#' 
#'
#' @param dxd A \code{DEXSeqDataSet} object
#' @param coefsMat A matrix of relative exon usage coefficients
#' @param bjps A BiocParam instance
#' @param mm A model matrix
#' 
#' @return A vector of dispersion estimates
#' @export
#' @import DEXSeq
#' @import statmod
#' @importFrom abind abind
estimateDispersionsParallel <- function(dxd, coefsMat, bjp, mm=NULL){
    countMatrix <- counts(dxd)
    sizeFactors <- colData(dxd)$sizeFactor
    isThis <- colData(dxd)$exon == "this"
    indexes1 <- getIndexList(dxd)
    if( is.null(mm) ){
        mm <- modelMatrixREUC(dxd, indexes1)
    }
    n <- bpnworkers(bjp)
    toSplit <- sort(rep(1:n, length.out=length(rownames(countMatrix))))
    matList <- lapply( split( rownames(countMatrix), toSplit ), function(x){
        countMatrix[x,]
    })
    dispDf <- bplapply( matList,
             function(x, coefsMatL, testableVectorL, mmL, sizeFactorsL, isThisL ){
                 estimateDispersionsForMat( x, coefsMatL,
                                           testableVectorL, mmL, sizeFactorsL,
                                           isThisL )
             }, coefsMatL=coefsMat, testableVectorL=testableVector,
             mmL=mm, sizeFactorsL=sizeFactors, isThisL=isThis,
             BPPARAM=bjp)
    do.call(rbind, dispDf)
}

arrangeInto3DArray <- function(dxd, allCoefs){
    crossIdx <- getIndexList(dxd)[["crossIdx"]]
    crossCoefs <- do.call( abind, c( along=3,
                                    tapply( 1:nrow(crossIdx),
                                           crossIdx$tissue,
                                           function(i){
                                               x <- allCoefs[,crossIdx$col[i]];
                                               colnames(x) <- crossIdx$individual[i];
                                               x
                                           } ) ) )
    names(dimnames(crossCoefs)) <- c("exon", "individual", "tissue")
    crossCoefs
}
