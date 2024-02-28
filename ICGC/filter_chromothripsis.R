library(tidyverse)

setwd('/path/to/project/ICGC_data/fuse')


# PCAWG_chromothripsisOverlap.txt.gz
# 
# File containing annotated chromothripsis event calls for all PCAWG samples witht he following columns:
#   
#   eventID: unique ID concatenating "chromothripsis", PCAWG tumour wgs aliquot ID, and chromosome+arm
# samplename: CAWG tumour wgs aliquot ID
# chr: chromosome+arm involved
# start: start position of the cluster of breakpoints
# end: end position of the cluster of breakpoints
# size: end-start
# pClonal: probability of being clonal=normalised product of probabilities of being clonal across SVs within the region (normalised so that pClonal+pSubclonal=1)
# pSubclonal: probability of being subclonal=normalised product of probabilities of being clonal across SVs within the region (normalised so that pClonal+pSubclonal=1)
# nbSVs: number of SVs with at least one mate within chromosome:start-end
# countsClonalNA: sum of probabilities of being clonal[NA] across SVs
# countsClonalEarly: sum of probabilities of being clonal[early] across SVs
# countsClonalLate: sum of probabilities of being clonal[late] across SVs
# averageSegmentSize: mean of the distance between consecutive breakpoints
# pvalue.exponentialTest: P-value of the exponential test= KS-test that segment sizes follow an exponential with rate=1/mean(averageSegmentSize)
# cov3StateMode: fraction covered by the 3 copy number states with the highest number of segments
# cov3StateMax: max fraction of region covered by 3 copy number states
# cov3stateMinExpected: minimum fraction expected given number of segments to call chromothripsis
# N.breakpoints.SV.BAF: total number of breakpoints (SV-related or BAF/logr-segmentation related)
# CNA.sizes.states: total size of segments per state
# CNA.number.states: number of segments per state
# chromosome: chromosome
# p.RandomOrder: P-Value of a spearman correlation: H0=the orders of the mate positions of the intrachromosomal SVs are random
# p.RandomJoins: P-Value of a Chi^2 test for the randomness of the type of SVs: H0=probabilities for H2H, T2T, Inv, and Del are equal, i.e. (1/4,1/4,1/4,1/4)
# fractionInterChromosomalSV: fraction of inter-chromosomal SVs
# nSV.per.chromosomes: number of SVs/mates falling in the region per chromosome (do not use!)
# overlapEvents: how many SVs link region to other cluster of SVs
# chrInvolved: which additional chromosome is involved
# NchrInvolved: how many additional chromosome involved
# isAmplified: is the region amplified (see methods)
# histology: histology of the sample
# distantTeloCentro: is it distant from the telomere (see methods)
# Pass_CNA3StatesCoverage: pass criteria for limited copy number states (see methods)
# Pass_SegmentSizesBelow0.05: pass exponential test (p-value<0.05)
# Pass_MinNbSegmentsAbove30: pass number of breakpoints above 30
# Pass_RandomJoinsAbove0.05: pass random SV types
# Pass_RandomMatePositionAbove0.05: pass random position
# Pass_MinNbSVsAbove10: pass min number of SV=10
# Pass_ALL: does it pass all criteria for chromothripsis (see methods)
# RetrievedFromSVConnection: was it retrieved through a connection with other chromothripsis region
# FinalCalls: final calls ("Chromothripsis" or "Cluster")
# LipoSarclike: is it liposarc-like (see methods)
# OverlapCTBetween2Callers: this is a boolean value and is true for final calls, i.e. overlapping with the set of chromothripsis calls reported by Cortés-Ciriano et al 1


c <- read_tsv("chromothripsis.tsv")
p <- read_tsv("purity_ploidy.tsv")

f <- read_tsv("pcawg_files.tsv")
f <- rename(f, "icgc_donor_id"="ICGC Donor")
f <- f %>% select(c("samplename", "icgc_donor_id"))
f
d <- merge(c, p, all=T)
d <- merge(d, f, by="samplename", all=T)

clin <- read_tsv("clinical.tsv")
clin
d <- merge(d, clin, by="icgc_donor_id", all=T)

d <- d %>% arrange(samplename)

breastadenoca <- d %>% filter(histology=="Breast-AdenoCA")
write.table(breastadenoca, file="breast.csv")

write.table(d, file="merged.csv", row.names = FALSE)

chrom_sample <- d %>% filter(FinalCalls=="Chromothripsis") %>% select(samplename) %>% unique() %>% as.list()
chrom_sample
clin %>% filter(submitted_donor_id %in% chrom_sample)
