## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gwid)

## ----out.width="100%",include=TRUE, fig.align="center", fig.cap=c("gwid pipeline"), echo=FALSE----
knitr::include_graphics("../man/figures/final-copy-arrow.png")

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("gwid")

## ----eval = FALSE-------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("soroushmdg/gwid")

## -----------------------------------------------------------------------------
# install.packages("piggyback")
piggyback::pb_download(repo = "soroushmdg/gwid",
            tag = "v0.0.1",
            dest = tempdir())
ibd_data_file <- paste0(tempdir(), "//chr3.ibd")
genome_data_file <- paste0(tempdir(), "//chr3.gds")
phase_data_file <- paste0(tempdir(), "//chr3.vcf")
case_control_data_file <- paste0(tempdir(), "//case-cont-RA.withmap.Rda")

## ----example------------------------------------------------------------------
library(gwid)

# case-control data
case_control <- gwid::case_control(case_control_rda = case_control_data_file)
names(case_control) #cases and controls group
summary(case_control) # in here, we only consider cases,cont1,cont2,cont3
#groups in the study
case_control$cases[1:3] # first three subject names of cases group

# read SNP data (use SNPRelate to convert it to gds) and count number of
#minor alleles  
snp_data_gds <- gwid::build_gwas(gds_data = genome_data_file,
                                 caco = case_control,
                                 gwas_generator = TRUE)
class(snp_data_gds)
names(snp_data_gds)
head(snp_data_gds$snps) # it has information about counts of minor alleles 
#in each location.

# read haplotype data (output of beagle)
haplotype_data <- gwid::build_phase(phased_vcf = phase_data_file,
                                    caco = case_control)
class(haplotype_data)
names(haplotype_data)
dim(haplotype_data$Hap.1) #22302 SNP and 1911 subjects

# read IBD data (output of Refined-IBD)
ibd_data <- gwid::build_gwid(ibd_data = ibd_data_file,
                             gwas = snp_data_gds)
class(ibd_data)
ibd_data$ibd # refined IBD output
ibd_data$res # count number of IBD for each SNP location 


## ----fig.width=7--------------------------------------------------------------
# plot count of IBD in chromosome 3
plot(ibd_data,
     y = c("cases","cont1"),
     ly=FALSE) 

# Further investigate location between 117M and 122M
# significant number of IBD's in group cases, compare to cont1, cont2 and cont3.
plot(ibd_data,
     y = c("cases","cont1"),
     snp_start = 117026294,
     snp_end = 122613594,
     ly=FALSE) 


## ----, fig.width=7------------------------------------------------------------
model_fisher <- gwid::fisher_test(ibd_data,
                                  case_control,
                                  reference = "cases",           
                                  snp_start = 117026294,
                                  snp_end = 122613594)

class(model_fisher)

plot(model_fisher, 
     y = c("cont1","cont2"),
     ly=FALSE,
     log_transformation = TRUE)

plot(model_fisher, 
     y = c("cont1","cont2"),
     QQplot = TRUE)


## ----eval=FALSE---------------------------------------------------------------
#  model_permutation <- gwid::permutation_test(ibd_data,gwas = snp_data_gds,
#                                              reference = "cases",
#                                              snp_start = 117026294,
#                                              snp_end = 122613594,
#                                              nperm = 100)
#  plot(model_permutation,
#       y = c("cont1","cont2","cont3"),
#       log_transformation=TRUE)

## -----------------------------------------------------------------------------
hap_str <- gwid::haplotype_structure(ibd_data,
                                     phase = haplotype_data,
                                     w = 10,
                                     snp_start = 117026294,
                                     snp_end = 122613594)
class(hap_str)

hap_str[sample(1:nrow(hap_str),size = 5),] # structures column 
#have haplotype of length w=10 

## -----------------------------------------------------------------------------
haplo_freq <- gwid::haplotype_frequency(hap_str)

## ----fig.width=7--------------------------------------------------------------

# plot haplotype counts in first window (nwin=1).
 plot(haplo_freq,
   y = c("cases", "cont1"),
   plot_type = "haplotype_structure_frequency",
   nwin = 1, type = "version1",
   ly=FALSE
 )


