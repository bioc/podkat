\name{plot}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,AssocTestResultRanges,missing-method}
\alias{plot,AssocTestResultRanges,character-method}
\alias{plot,AssocTestResultRanges,GRanges-method}
\alias{plot,GenotypeMatrix,missing-method}
\alias{plot,GenotypeMatrix,factor-method}
\alias{plot,GenotypeMatrix,numeric-method}
\alias{plot,GRanges,character-method}
\title{Plotting functions}
\description{
  Functions for visualizing association test results by means of a
  Manhattan plot and for visualizing genotypes
}
\usage{
\S4method{plot}{AssocTestResultRanges,missing}(x, y, cutoff=0.05,
     which=c("p.value", "p.value.adj", "p.value.resampled",
             "p.value.resampled.adj"), showEmpty=FALSE,
     as.dots=FALSE, pch=19, col=c("darkgrey", "grey"), scol="red",
     lcol="red", xlab=NULL, ylab=NULL, ylim=NULL, lwd=1, cex=1,
     cexXaxs=1, cexYaxs=1, srt=0, adj=c(0.5, 1), ...)
\S4method{plot}{AssocTestResultRanges,character}(x, y, cutoff=0.05,
     which=c("p.value", "p.value.adj", "p.value.resampled",
             "p.value.resampled.adj"), showEmpty=FALSE,
     as.dots=FALSE, pch=19, col=c("darkgrey", "grey"), scol="red",
     lcol="red", xlab=NULL, ylab=NULL, ylim=NULL, lwd=1, cex=1, 
     cexXaxs=1, cexYaxs=1, srt=0, adj=c(0.5, 1), ...)
\S4method{plot}{AssocTestResultRanges,GRanges}(x, y, cutoff=0.05,
     which=c("p.value", "p.value.adj", "p.value.resampled",
             "p.value.resampled.adj"), showEmpty=FALSE,
     as.dots=FALSE, pch=19, col="darkgrey", scol="red", lcol="red",
     xlab=NULL, ylab=NULL, ylim=NULL, lwd=1, cex=1,
     cexXaxs=1, cexYaxs=1, ...)
\S4method{plot}{GenotypeMatrix,missing}(x, y, col="black",
     labRow=NULL, labCol=NULL, cexXaxs=(0.2 + 1 / log10(ncol(x))),
     cexYaxs=(0.2 + 1 / log10(nrow(x))), srt=90, adj=c(1, 0.5))
\S4method{plot}{GenotypeMatrix,factor}(x, y, col=rainbow(length(levels(y))),
     labRow=NULL, labCol=NULL, cexXaxs=(0.2 + 1 / log10(ncol(x))),
     cexYaxs=(0.2 + 1 / log10(nrow(x))), srt=90, adj=c(1, 0.5))
\S4method{plot}{GenotypeMatrix,numeric}(x, y, col="black", ccol="red", lwd=2,
     labRow=NULL, labCol=NULL, cexXaxs=(0.2 + 1 / log10(ncol(x))),
     cexYaxs=(0.2 + 1 / log10(nrow(x))), srt=90, adj=c(1, 0.5))
\S4method{plot}{GRanges,character}(x, y, alongGenome=FALSE, 
     type=c("r", "s", "S", "l", "p", "b", "c", "h", "n"),
     xlab=NULL, ylab=NULL, col="red", lwd=2,
     cexXaxs=(0.2 + 1 / log10(length(x))), cexYaxs=1,
     frame.plot=TRUE, srt=90, adj=c(1, 0.5), ...)
}
\arguments{
  \item{x}{an object of class
    \code{\linkS4class{AssocTestResultRanges}},
    \code{\linkS4class{GenotypeMatrix}}, or \code{\linkS4class{GRanges}}}
  \item{y}{a character string, \code{GRanges} object, or factor}
  \item{cutoff}{significance threshold}
  \item{which}{a character string specifying which p-values to plot;
    if \dQuote{p.value} (default), raw p-values are plotted.
    Other options are \dQuote{p.value.adj} (adjusted p-values),
    \dQuote{p.value.resampled} (resampled p-values), and
    \dQuote{p.value.resampled.adj} (adjusted resampled p-values).
    If the requested column is not present in the input object \code{x},
    the function stops with an error message.}
  \item{showEmpty}{if \code{FALSE} (default), p-values of regions that
    did not contain any variants are omitted from the plot.}
  \item{as.dots}{if \code{TRUE}, p-values are
    plotted as dots/characters in the center of the genomic region.
    If \code{FALSE} (default), p-values are plotted as lines stretching
    from the starts to the ends of the corresponding genomic regions.}
  \item{pch}{plotting character used to plot a single p-value, ignored
    if \code{as.dots=FALSE}; see \code{\link{points}} for details.}
  \item{col}{plotting color(s); see details below}
  \item{scol}{color for plotting significant p-values (i.e. the ones
    passing the significance threshold)}
  \item{lcol}{color for plotting the significance threshold line}
  \item{xlab}{x axis label; if \code{NULL} (default) or \code{NA},
    \code{plot} makes an automatic choice}
  \item{ylab}{y axis label; if \code{NULL} (default) or \code{NA},
    \code{plot} makes an automatic choice}
  \item{ylim}{y axis limits; if \code{NULL} (default) or \code{NA},
     \code{plot} makes an automatic choice; if user-specified,
     \code{ylim} must be a two-element numeric vector with the first
     element being 0 and the second element being a positive value.}
  \item{lwd}{line thickness; in Manhattan plots, this parameter
    corresponds to the thickness of the significance threshold line.
    When plotting genotype matrices along with continuous traits, this
    is the thickness of the line that corresponds to the trait.}
  \item{cex}{scaling factor for plotting p-values; see
    \code{\link{points}} for details.}
  \item{labRow,labCol}{row and column labels; set to \code{NA} to switch
    labels off; if \code{NULL}, rows are labeled by sample names
    (\code{rownames(x)}) and columns are labeled by variant names
    (\code{names(variantInfo(x))}).}
  \item{cexXaxs,cexYaxs}{scaling factors for axes labels}
  \item{ccol}{color of the line that plots the continuous trait along
    with a genotype matrix}
  \item{srt}{rotation angle of text labels on horizontal axis
    (see \code{\link{text}} for details); ignored if standard numerical
    ticks and labels are used.}
  \item{adj}{adjustment of text labels on horizontal axis
    (see \code{\link{text}} for details); ignored if standard numerical
    ticks and labels are used.}
  \item{alongGenome}{plot along the genome or per region (default); see
    details below.}
  \item{type}{type of plot; see \code{\link{plot.default}} for
    details. Additionally, the type \dQuote{r} is available (default)
    which plots horizontal lines along the regions of \code{x}.}
  \item{frame.plot}{whether or not to frame the plotting area (see
    \code{\link[graphics:plot]{plot}}; default: \code{TRUE})}
  \item{...}{all other arguments are passed to
    \code{\link[graphics:plot]{plot}}.}
}
\details{
  If \code{plot} is called for an
  \code{\linkS4class{AssocTestResultRanges}}
  object without specifying the second argument \code{y},
  a so-called Manhattan plot is produced. The x axis corresponds to the
  genome on which the \code{\linkS4class{AssocTestResultRanges}}
  \code{x} is based and the y axis shows absolute values of
  log-transformed p-values. The \code{which} argument determines
  which p-value is plotted, i.e. raw p-values, adjusted p-values,
  resampled p-values, or adjusted resampled p-values. The \code{cutoff}
  argument allows for setting a significance threshold above which
  p-values are plotted in a different color (see above).

  The optional \code{y} argument can be used for two purposes: (1) if it is
  a character vector containing chromosome names (sequence names), it
  can be used for specifying a subset of one or more chromosomes to be
  plotted. (2) if \code{y} is a \code{\linkS4class{GRanges}} object of length
  1 (if longer, \code{plot} stops with an error), only the genomic
  region corresponding to \code{y} is plotted.

  The \code{col} argument serves for specifying the color for plotting
  insignificant p-values (i.e. the ones
  above the significance threshold); if the number of colors is
  smaller than the number of chromosomes, the vector is recycled.
  If \code{col} is a single color, all insignificant p-values are
  plotted in the same color. If \code{col} has two elements (like the
  default value), the insignificant p-values of different chromosomes
  are plotted with alternating colors. It is also possible to produce
  density plots of p-values by using semi-transparent colors (see,
  e.g., \code{\link{rgb}} or \code{\link{hsv}} for information on how
  to use the \code{alpha} channel).

  If \code{plot} is called for a \code{\linkS4class{GenotypeMatrix}}
  object \code{x} and no \code{y} argument, the matrix is visualized in
  a heatmap-like fashion, where two major alleles are displayed in white,
  two minor alleles are displayed in the color passed as \code{col}
  argument, and the heterozygotous case (one minor, one major) is
  displayed in the color passed as \code{col} argument, but with 50\%
  transparency. The arguments \code{cexYaxs} and \code{cexXaxs} can be
  used to change the scaling of the axis labels.

  If \code{plot} is called for a \code{\linkS4class{GenotypeMatrix}}
  object \code{x} and a factor \code{y}, then the factor \code{y} is
  interpreted as a binary trait. In this case, the rows of the genotype
  matrix \code{x} are reordered such that rows/samples with the same
  label are plotted next to each other. Each such group can be plotted in a
  different color. For this purpose, a vector of colors can be passed
  as \code{col} argument.

  If \code{plot} is called for a \code{\linkS4class{GenotypeMatrix}}
  object \code{x} and a numeric vector \code{y}, then the vector \code{y} is
  interpreted as a continuous trait. In this case, the rows of the genotype
  matrix \code{x} are reordered according to the trait vector \code{y}
  and the genotype matrix is plotted as described above. The trait
  \code{y} is superimposed in the plot in color \code{ccol} and with
  line width \code{lwd}. If the null model has been trained with
  covariates, it also makes sense to plot the genotype against the
  null model residuals, since these are exactly the values that the
  genotypes were tested against.
  
  If \code{plot} is called for a \code{\linkS4class{GRanges}}
  object \code{x} and a character string \code{y}, then \code{plot}
  checks whether \code{x} has a metadata column named \code{y}.
  If it exists, this column is plotted against the regions in
  \code{x}. If \code{alongGenome} is \code{FALSE} (which is the
  default), the regions
  in \code{x} are arranged along the horizontal plot axis
  with equal widths and in the same order as contained in \code{x}.
  If the regions in \code{x} are named, then the names are used as
  axis labels and the argument \code{cexXaxs} can be used to scale the
  font size of the names. If \code{alongGenome} is \code{TRUE},
  the metadata column is plotted against genomic positions. The knots
  of the curves are then positioned at the positions given in the 
  \code{\linkS4class{GRanges}} object \code{x}. For types \dQuote{s},
  \dQuote{S}, \dQuote{l}, \dQuote{p}, \dQuote{b}, \dQuote{c}, and
  \dQuote{h}, knots are placed in the middle of the genomic regions
  contained in \code{x} if they are longer than one base. For type
  \dQuote{r}, regions are plotted as lines exactly stretching between
  the start and end coordinates of each region in \code{x}.
}
\value{
  returns an invisible numeric vector of length 2 containing the
  y axis limits
}
\author{Ulrich Bodenhofer \email{bodenhofer@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/podkat}}
\seealso{\code{\linkS4class{AssocTestResultRanges}},
  \code{\linkS4class{GRanges}}}
\examples{
## load genome description
data(hgA)

## partition genome into overlapping windows
windows <- partitionRegions(hgA)

## load genotype data from VCF file
vcfFile <- system.file("examples/example1.vcf.gz", package="podkat")
Z <- readGenotypeMatrix(vcfFile)

## plot some fraction of the genotype matrix
plot(Z[, 1:25])

## read phenotype data from CSV file (continuous trait + covariates)
phenoFile <- system.file("examples/example1log.csv", package="podkat")
pheno <-read.table(phenoFile, header=TRUE, sep=",")

## train null model with all covariates in data frame 'pheno'
nm.log <- nullModel(y ~ ., pheno)

## perform association test
res <- assocTest(Z, nm.log, windows)
res.adj <- p.adjust(res)

## plot results
plot(res)
plot(res, cutoff=1.e-5, as.dots=TRUE)
plot(res.adj, which="p.value.adj")
plot(res.adj, reduce(windows[3:5]), which="p.value.adj")

## filter regions
res.adj.f <- filterResult(res.adj, filterBy="p.value.adj")

## plot genotype grouped by target
sel <- which(overlapsAny(variantInfo(Z), reduce(res.adj.f)))
plot(Z[, sel], factor(pheno$y))
plot(Z[, sel], residuals(nm.log), srt=45)

## compute contributions
contrib <- weights(res.adj.f, Z, nm.log)
contrib[[1]]

## plot contributions
plot(contrib[[1]], "weight.raw")
plot(contrib[[1]], "weight.contribution", type="b", alongGenome=TRUE)
}
\keyword{methods}
