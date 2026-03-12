getwd()
set.seed(0)
#Manu Singh 5 April 2020
#Run regionR to get sinificant values for family. The idea is to check if the peaks are enriched in a family.

#BiocManager::install("regioneR")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library("regioneR")

args = commandArgs(trailingOnly=TRUE)

Human_class2_TEs <- read.delim(args[1], stringsAsFactors = F, header=F)

All_Peaks <- read.delim(args[2], stringsAsFactors = F, header=F)
All_Peaks_GR <- toGRanges(All_Peaks)

human.genome <- getGenomeAndMask("hg38", mask=NA)$genome

human.chrs <- filterChromosomes(human.genome, keep.chr=unlist(unique(All_Peaks[,1])))

Total_TEs_family <- (unique(Human_class2_TEs$V4))
Total_TEs_family <- Total_TEs_family[!grepl("^\\(.+\\)n$", Total_TEs_family)]
Total_TEs_familyDF <- as.data.frame(Total_TEs_family)
Total_TEs_familyDF$Pvalue <- 2 

Save_Model_Also <- list()
#for(numfam in 1:length(Total_TEs_family)){
for(numfam in 1:10){
  
  cat("Family Number", numfam, Total_TEs_family[numfam],"\n")
  Take_one_family <- Human_class2_TEs[(Human_class2_TEs$V4%in%Total_TEs_family[numfam]),]
  
  Take_one_familyV2 <- Take_one_family[(((Take_one_family[,1]))%in%unlist(unique(All_Peaks[,1]))),]
  
  
  Take_one_familyV2_GR <- toGRanges(Take_one_familyV2)
  
  
  #Human_class2_TEsV2 <- Human_class2_TEs[(((Human_class2_TEs[,1]))%in%unlist(unique(All_Peaks[,1]))),]
  #Human_class2_TEs_GR <- toGRanges(Human_class2_TEsV2)
  #numOverlaps(Take_one_familyV2_GR, All_Peaks_GR)
  #random.RS <- resampleRegions(Take_one_familyV2_GR, universe=Human_class2_TEs_GR)
  
  pt <- permTest(A=Take_one_familyV2_GR, ntimes=1000, randomize.function=randomizeRegions, genome=human.chrs,
                 evaluate.function=numOverlaps, B=All_Peaks_GR, verbose=TRUE, alternative = "greater", force.parallel = TRUE)
  Total_TEs_familyDF[numfam,2] <- pt$numOverlaps$pval
  Save_Model_Also[[numfam]] <- pt
  
}

#saveRDS(Total_TEs_familyDF,file="Total_TEs_familyDF_TakingGenomehg19ForRandomization.rds")
#saveRDS(Save_Model_Also,file="Save_Model_Also_TakingGenomehg38ForRandomization.rds")

write.csv(Total_TEs_familyDF, paste0(args[3]), row.names = FALSE)


























