#########################################
# Get all possible APA clusters via hatal
# v2-0.1.1
#########################################

#########################################
## haptal
library(hatpal2)
#Sys.setlocale("LC_ALL","English")

#########################################
## Load paras
# args=commandArgs(T)
# in_bamfile <- args[1]
# in_nSBfile <- args[2]
# out_dir <- args[3]
# gtf_file <- args[4]

in_bamfile <- "data/pbmc_F256_c1c2.bam"
in_nSBfile <- "data/pbmc_F256_c1c2.bam"
out_dir <- "tmp/"
fasta_file <- "data/hg19_chr1_chr2.fa"
gtf_file <- "anno_hg19/Homo_sapiens.GRCh37.87.chr.gtf"


#########################################
## Extracting the chromosomes information
chr_info <- ExtractChrinfo(in_bamfile)

#########################################
## Generate coordinate files
ExtractAPA(in_bamfile, in_nSBfile, chr_info, out_dir)

#########################################
## Negative strand
SA_prebed <- paste0(out_dir, "negative_strand/","sa.prebed")
UNI_prebed <- paste0(out_dir, "negative_strand/","rc_uni.prebed")
END_prebed <- paste0(out_dir, "negative_strand/","rc_end.prebed")

ClusterAPA(SA_prebed, UNI_prebed, paste0(out_dir, "negative_strand/"))
FilterAPAc(paste0(out_dir, "negative_strand/APA_clusters.out"),
           fasta_file, "-", paste0(out_dir, "negative_strand/"))
MapAPAc(END_prebed, SA_prebed,
        paste0(out_dir, "negative_strand/APA_clusters.filtered.out"), paste0(out_dir, "negative_strand/"))
AnnoAPAc(paste0(out_dir, "negative_strand/APA_map.out"), gtf_file, "-", TRUE) ### !!!change it to FALSE!

#########################################
## Positive strand
SA_prebed <- paste0(out_dir, "positive_strand/","sa.prebed")
UNI_prebed <- paste0(out_dir, "positive_strand/","rc_uni.prebed")
END_prebed <- paste0(out_dir, "positive_strand/","rc_end.prebed")

ClusterAPA(SA_prebed, UNI_prebed, paste0(out_dir, "positive_strand/"))
FilterAPAc(paste0(out_dir, "positive_strand/APA_clusters.out"),
           fasta_file, "-", paste0(out_dir, "positive_strand/"))
MapAPAc(END_prebed, SA_prebed,
        paste0(out_dir, "positive_strand/APA_clusters.filtered.out"), paste0(out_dir, "positive_strand/"))
AnnoAPAc(paste0(out_dir, "positive_strand/APA_map.out"), gtf_file, "+", TRUE) ### !!!change it to FALSE!

#########################################
## Creat a count matrix
CountAPAc(paste0(out_dir, "positive_strand/APA_map.out.anno"),
        paste0(out_dir, "negative_strand/APA_map.out.anno"),
        in_bamfile, paste0(out_dir, "_CBlist.txt"), chr_info, out_dir)

