#' ---
#' title: "Scatterplot code"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' ## Libraries
suppressPackageStartupMessages(library(genomeIntervals))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(GenomicRanges))

#' ## Create a color palette
pal=brewer.pal(8,"Dark2")

#' # Data
#' Read in the gff file
gff3 <- readGff3("/mnt/picea/projects/u2015027/intgen_w_exon_file.gff3",quiet = TRUE)

#' Plot the score distribution. As I do not know how you run your command, I can't be sure what the score is, but you can probably use it as a threshold to limit the number of hits.
plot(density(gff3$score),col=pal[1],lwd=2,main="Score distribution")

#' # Analysis
#' First we find the number of loci
message(paste("There are",nrow(gff3),"loci"))

#' Next, we look at the origin of the hit (i.e. the trinity transcript coordinates)
target <- getGffAttribute(gff3,"Target")

#' Here we create a genomic range object to store the transcript sequence name and coordinates
target.grng <- GRanges(seqnames = sub(" .*","",target),
                ranges = IRanges(start=as.integer(gsub(".*_i[0-9]+ | .*","",target)),
                                 end=as.integer(sub(".* ","",target))))

#' And we report the number of transcripts and their cumulative mapped size
message(paste("There are",length(seqlevels(target.grng)),"transcripts"))
plot(density(sum(split(width(target.grng),seqnames(target.grng)))),
     main="mapped bp coverage per transcript",lwd=2,col=pal[1])

#' Next we look at the genomic loci mapped by the transcripts
#' 
#' First creating a grange object, which we "reduce" to remove overlapped region (e.g. think alternative splicing)
genome.grng <- reduce(GRanges(seqnames = seqnames(gff3),
                       ranges = IRanges(start=gff3[,1],
                                        end=gff3[,2])))

#' And we report on those
message(paste("Mapping to",length(genome.grng),"genomic loci"))
message(paste("Located on",length(seqlevels(genome.grng)),"scaffolds"))
message(paste("Assuming a 1 gene per scaffold (as is the case for most spruce scaffolds), we have the following distribution of exons"))
barplot(table(table(seqnames(genome.grng))),las=2,main="Presumptive exon number distribution")
message(paste("Covering the following amount of sequence per scaffold"))
plot(density(sum(split(width(genome.grng),seqnames(genome.grng)))),
     main="exon bp coverage per scaffold",lwd=2,col=pal[1])

#' # SessionInfo
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
