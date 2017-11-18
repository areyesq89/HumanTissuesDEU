#' @title Plot REUCs
#' @description
#' Given a gene identifier and the vector of REUCs,
#' this function plots heatmaps for each exonic region
#' of the gene.
#'
#' @param crossCoefs A 3D matrix of relative exon usage coefficients 
#' @param geneName A gene identifier string.
#' @param exons A character vector with exon identifiers
#' @param colLim Range of values for colors
#' @param asZ Logical
#' 
#' @return A heatmap of relative exon usage coefficients
#' @export
#' @import GenomicFeatures
#' @import dplyr
#' @import Gviz
plotGeneREUCs <- function(crossCoefs, geneName, exons=NULL, colLim=NULL, asZ=TRUE){
    spID <- strsplit( dimnames(crossCoefs)[["exon"]], ":" )
    exonRows <- which(sapply(spID, "[[", 1) %in% geneName)
    coefsThisGene <- crossCoefs[exonRows,,]
    dimnames(coefsThisGene)$exon <- sapply( spID[exonRows], "[[", 2 )
    coefsThisGene <- coefsThisGene[ apply( coefsThisGene, 1, function(x){!any( is.na(x) )} ),,]
    if( asZ ){
        for( i in dimnames(coefsThisGene)$exon ){
            pr <- coefsThisGene[i,,]
            pr <- (pr - median(pr))/sd(pr)
            coefsThisGene[i,,] <- pr
        }
    }
    if( is.null(colLim) ){
        qtls <- quantile( coefsThisGene, seq(0, 1, 0.01))
        colLim <- max(abs(qtls[c("2%", "98%")]))
    }
    dfT <- expand.grid(
        tissue=dimnames(coefsThisGene)$tissue,
        individual=dimnames(coefsThisGene)$individual,
        exon=dimnames(coefsThisGene)$exon )
    dfT <- mutate( dfT, REUC=NA_real_ )
    for( i in seq_len(nrow(dfT)) ){
        dfT$REUC[i] <- coefsThisGene[dfT$exon[i], dfT$individual[i], dfT$tissue[i]]
    }
    bPal <- colorRampPalette(brewer.pal(9, name="RdBu"), interpolate="spline", bias=1)(3)
    dfT <- mutate( dfT, REUC=pmin(pmax(REUC, -colLim), colLim) )
    if( !is.null( exons ) ){
        dfT <- filter( dfT, exon %in% exons )
    }
    g1 <- ggplot( dfT, aes(individual, tissue)) + geom_tile(aes(fill=REUC), colour="black") + 
        facet_wrap(~exon) + 
        scale_fill_gradient2( low=bPal[3], high=bPal[1], mid=bPal[2], midpoint=0, limits=c(-colLim, colLim)) +
        theme( axis.text.x=element_blank(), strip.text.x=element_text(color="black", size=12),
               axis.text.y=element_text(color="black", size=11), axis.title=element_text(size=14),
               legend.text=element_text(size=12), legend.title=element_text(size=13) )
    g1
}

#' @title Plot sashimi plots
#' @description
#' This function receives bam files and plots sashimi plots
#' for the specified genomic coordinates.
#'
#' @param sampleVec Path to bam files.
#' @param nameVec Names of the labels for the shingles.
#' @param geneID A gene identifier.
#' @param transcriptDb A transcriptDb object.
#' @param coords Genomic coordinates to plot.
#' @param geneTrack A GeneTrack object from the Gviz package.
#' @param transcriptIntrons If specified, only junctions annotated in the transcript identifiers specified are considered.
#' @param offset Number of bases to increase the window size.
#' @param type Either 'coverage' or 'sashimi', specified what should be plotted.
#' @param plotTranscripts Logical. Whether to plot transcripts.
#' @param rnaYLim Limits for the y-axis for the RNA-seq tracks.
#' @param fantomObject A GenomicRanges object containing FANTOM data.
#' @param fantomTissues Which tissues from the FANTOM data to plot.
#' @param fantomLabs A vector of labels for the FANTOM tracks.
#' @param fantomYLim Limits for the y-axis for the FANTOM tracks.
#' @param highlight A GenomicRange to highlight.
#' @param highlightOffset Number of bases to increase the highlight area.
#' @param ... Further parameters for the plotTrack function.
#' 
#' @return A sashimi plot.
#' @export
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import dplyr
#' @import Gviz
plotSashimi <- function( sampleVec, nameVec, geneID=NULL, transcriptDb,
                         coords=NULL, geneTrack, transcriptIntrons=NULL, offset=0,
                         type=c("coverage", "sashimi"), plotTranscripts=FALSE, rnaYLim=NULL,
                         fantomObject=NULL, fantomTissues=NULL, fantomLabs=NULL, fantomYLim=NULL,
                         highlight=NULL, highlightOffset=0, ...){
    intronsByTx <- intronsByTranscript(transcriptDb, use.names=TRUE)
    introns <- unique(unlist(intronsByTx))
    if( !is( rnaYLim, "list") ){
        rnaYLim <- list( rnaYLim )
    }
    sampleTracks <-
        mapply( function(x, y, z){
                   AlignmentsTrack(
#                       file.path( path, "alignments", x,
#                                sprintf("%s_Aligned.sortedByCoord.out.bam", x) ),
                       x, isPaired=TRUE,
                       name=y, ylim=z, type=type) }, sampleVec, nameVec, rnaYLim )
    if( is.null(coords) ){
        coords <- range( geneTrack@range[geneTrack@range$transcript %in% geneID])
        geneTrack <- geneTrack[geneTrack@range$transcript %in% geneID,]
    }
    if( !is.null( fantomObject) ){
        fantomObject <- subsetByOverlaps( fantomObject, coords )
        fantomRanges <- rowRanges( fantomObject )
        mcols( fantomRanges ) <- NULL
        this <- colData(fantomObject)$exon == "this"
        spCols <- split(seq_len(ncol(fantomObject[,this])),
                      colData( fantomObject )$tissue[this] )
        cageExpr <- sapply( spCols[fantomTissues], function(x){
                        rowMeans(log2( featureCounts( fantomObject,
                            normalized=TRUE ) + 1 )[,x]) })
        names( fantomLabs ) <- fantomTissues
        fantomTracks <- lapply( colnames(cageExpr), function(x){
                            gapRanges <- gaps( fantomRanges )
                            rangesToPlot <- fantomRanges
                            mcols( gapRanges ) <- 0
                            colnames( mcols( gapRanges ) ) <- x
                            mcols( rangesToPlot ) <- cageExpr[,x]
                            colnames( mcols( rangesToPlot ) ) <- x
                            allR <- sort( c( rangesToPlot, gapRanges  ))
                            allR <- allR[strand( allR ) == unique( strand(coords) ),]
                            DataTrack( allR, name=fantomLabs[x], ylim=fantomYLim, type="histogram") } )
        sampleTracks <- c( sampleTracks, fantomTracks )
    }
    chr=unique(as.character(seqnames(coords)))
    start=start(coords) - offset
    end=end(coords) + offset
    if( !is.null(transcriptIntrons) ){
        thisIntrons <- unique(unlist(intronsByTx[transcriptIntrons]))
    }else{
        thisIntrons <- introns
    }
    geneTrackER <- geneTrack
    if( plotTranscripts ){
        geneTrack <- GeneRegionTrack( transcriptDb )
    }
    sampleTracks <- c( sampleTracks, geneTrack )
    if( !is.null( highlight ) ){
        if( is.character( highlight ) ){
            highlight <- range( geneTrackER@range[geneTrackER@range$transcript %in% geneID &
                                                    geneTrackER@range$exon %in% highlight,] )
            start( highlight ) <- start( highlight ) - highlightOffset
            end( highlight ) <- end( highlight ) + highlightOffset
        }
        sampleTracks <- HighlightTrack( trackList=sampleTracks,
                                       start=start(highlight),
                                       width=width(highlight),
                                       chromosome=as.character(seqnames(highlight) ) ) 
    }
    plotTracks( c( sampleTracks, GenomeAxisTrack()),
        from=start, to=end, chromosome=chr,
               sashimiFilter=thisIntrons, 
               lwd.sashimiMax=5, fontcolor="black", ...)
}

getSampleIdentifiers <- function( tissues, individual, dxd ){
    sapply( tissues, function(x){
               unique( as.character(
                   colData(dxd)[colData(dxd)$tissue %in% x &
                                    colData(dxd)$individual %in% individual,
                                "sample"] ) )
           })
}

#transcriptDb <- loadDb( file.path(path, "objects", "GRCh38.sqlite") )
#options( ucscChromosomeNames=FALSE ) 
#load(file.path(path, "objects", "geneTrack.RData"))
