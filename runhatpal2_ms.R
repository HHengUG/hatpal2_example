#########################################
# Get all possible APA clusters via hatal
# v2-0.1.1
#########################################

#########################################
## haptal
library(hatpal2)
library(data.table)
#Sys.setlocale("LC_ALL","English")

#########################################
## Load paras

in_bamfile1 <- "data/pbmc_F256_c1c2.bam"
in_nSBfile1 <- "data/pbmc_F256_c1c2.bam"
in_bamfile2 <- "data/pbmc_F256_c1c2.bam"
in_nSBfile2 <- "data/pbmc_F256_c1c2.bam"
in_bamfile3 <- "data/pbmc_F256_c1c2.bam"
in_nSBfile3 <- "data/pbmc_F256_c1c2.bam"

out_dir1 <- "tmp/tmp1"
out_dir2 <- "tmp/tmp2"
out_dir3 <- "tmp/tmp3"

fasta_file <- "data/hg19_chr1_chr2.fa"

gtf_file <- "anno_hg19/Homo_sapiens.GRCh37.87.chr.gtf"


#########################################
## Extracting the chromosomes information
chr_info1 <- ExtractChrinfo(in_bamfile1)
chr_info2 <- ExtractChrinfo(in_bamfile2)
chr_info3 <- ExtractChrinfo(in_bamfile3)


#########################################
## Generate coordinate files
ExtractAPA(in_bamfile1, in_nSBfile1, chr_info1, out_dir1)
ExtractAPA(in_bamfile2, in_nSBfile2, chr_info2, out_dir2)
ExtractAPA(in_bamfile3, in_nSBfile3, chr_info3, out_dir3)


#########################################
## Negative strand
SA_prebed1 <- paste0(out_dir1, "negative_strand/","sa.prebed")
SA_prebed2 <- paste0(out_dir2, "negative_strand/","sa.prebed")
SA_prebed3 <- paste0(out_dir3, "negative_strand/","sa.prebed")
SA_prebed <-paste0(out_dir1, "negative_strand/","sa_all.prebed")
allprebed <- rbind(fread(SA_prebed1),
                   fread(SA_prebed2),
                   fread(SA_prebed3))
allprebed[, V3 := sum(V3), by = .(V1, V2)]
allprebed <- unique(allprebed)
fwrite(allprebed, SA_prebed, quote = FALSE, sep = "\t", col.names = FALSE)
rm(allprebed)
gc()


UNI_prebed1 <- paste0(out_dir1, "negative_strand/","rc_uni.prebed")
UNI_prebed2 <- paste0(out_dir2, "negative_strand/","rc_uni.prebed")
UNI_prebed3 <- paste0(out_dir3, "negative_strand/","rc_uni.prebed")
UNI_prebed <-paste0(out_dir1, "negative_strand/","uni_all.prebed")
allprebed <- rbind(fread(UNI_prebed1),
                   fread(UNI_prebed2),
                   fread(UNI_prebed3))
allprebed[, V3 := sum(V3), by = .(V1, V2)]
allprebed <- unique(allprebed)
fwrite(allprebed, UNI_prebed, quote = FALSE, sep = "\t", col.names = FALSE)
rm(allprebed)
gc()


END_prebed1 <- paste0(out_dir1, "negative_strand/","rc_end.prebed")
END_prebed2 <- paste0(out_dir2, "negative_strand/","rc_end.prebed")
END_prebed3 <- paste0(out_dir3, "negative_strand/","rc_end.prebed")
END_prebed <-paste0(out_dir1, "negative_strand/","END_all.prebed")
allprebed <- rbind(fread(END_prebed1),
                   fread(END_prebed2),
                   fread(END_prebed3))
allprebed[, V3 := sum(V3), by = .(V1, V2)]
allprebed <- unique(allprebed)
fwrite(allprebed, END_prebed, quote = FALSE, sep = "\t", col.names = FALSE)
rm(allprebed)
gc()


ClusterAPA(SA_prebed, UNI_prebed, paste0(out_dir1, "negative_strand/"))
FilterAPAc(paste0(out_dir1, "negative_strand/APA_clusters.out"),
           fasta_file, "-", paste0(out_dir1, "negative_strand/"))
MapAPAc(END_prebed, SA_prebed,
        paste0(out_dir1, "negative_strand/APA_clusters.filtered.out"), paste0(out_dir1, "negative_strand/"))
AnnoAPAc(paste0(out_dir1, "negative_strand/APA_map.out"), gtf_file, "-", FALSE) ### !!!change it to FALSE!

#########################################
## Positive strand
SA_prebed1 <- paste0(out_dir1, "positive_strand/","sa.prebed")
SA_prebed2 <- paste0(out_dir2, "positive_strand/","sa.prebed")
SA_prebed3 <- paste0(out_dir3, "positive_strand/","sa.prebed")
SA_prebed <-paste0(out_dir1, "positive_strand/","sa_all.prebed")
allprebed <- rbind(fread(SA_prebed1),
                   fread(SA_prebed2),
                   fread(SA_prebed3))
allprebed[, V3 := sum(V3), by = .(V1, V2)]
allprebed <- unique(allprebed)
fwrite(allprebed, SA_prebed, quote = FALSE, sep = "\t", col.names = FALSE)
rm(allprebed)
gc()


UNI_prebed1 <- paste0(out_dir1, "positive_strand/","rc_uni.prebed")
UNI_prebed2 <- paste0(out_dir2, "positive_strand/","rc_uni.prebed")
UNI_prebed3 <- paste0(out_dir3, "positive_strand/","rc_uni.prebed")
UNI_prebed <-paste0(out_dir1, "positive_strand/","uni_all.prebed")
allprebed <- rbind(fread(UNI_prebed1),
                   fread(UNI_prebed2),
                   fread(UNI_prebed3))
allprebed[, V3 := sum(V3), by = .(V1, V2)]
allprebed <- unique(allprebed)
fwrite(allprebed, UNI_prebed, quote = FALSE, sep = "\t", col.names = FALSE)
rm(allprebed)
gc()


END_prebed1 <- paste0(out_dir1, "positive_strand/","rc_end.prebed")
END_prebed2 <- paste0(out_dir2, "positive_strand/","rc_end.prebed")
END_prebed3 <- paste0(out_dir3, "positive_strand/","rc_end.prebed")
END_prebed <-paste0(out_dir1, "positive_strand/","END_all.prebed")
allprebed <- rbind(fread(END_prebed1),
                   fread(END_prebed2),
                   fread(END_prebed3))
allprebed[, V3 := sum(V3), by = .(V1, V2)]
allprebed <- unique(allprebed)
fwrite(allprebed, END_prebed, quote = FALSE, sep = "\t", col.names = FALSE)
rm(allprebed)
gc()

ClusterAPA(SA_prebed, UNI_prebed, paste0(out_dir1, "positive_strand/"))
FilterAPAc(paste0(out_dir1, "positive_strand/APA_clusters.out"),
           fasta_file, "-", paste0(out_dir1, "positive_strand/"))
MapAPAc(END_prebed, SA_prebed,
        paste0(out_dir1, "positive_strand/APA_clusters.filtered.out"), paste0(out_dir1, "positive_strand/"))
AnnoAPAc(paste0(out_dir1, "positive_strand/APA_map.out"), gtf_file, "+", FALSE) ### !!!change it to FALSE!

#########################################
## Creat a count matrix
CountAPAc(paste0(out_dir1, "positive_strand/APA_map.out.anno"),
          paste0(out_dir1, "negative_strand/APA_map.out.anno"),
          in_bamfile1, paste0(out_dir1, "_CBlist.txt"), chr_info1, out_dir1)

CountAPAc(paste0(out_dir1, "positive_strand/APA_map.out.anno"),
          paste0(out_dir1, "negative_strand/APA_map.out.anno"),
          in_bamfile2, paste0(out_dir2, "_CBlist.txt"), chr_info2, out_dir2)

CountAPAc(paste0(out_dir1, "positive_strand/APA_map.out.anno"),
          paste0(out_dir1, "negative_strand/APA_map.out.anno"),
          in_bamfile3, paste0(out_dir3, "_CBlist.txt"), chr_info3, out_dir3)
