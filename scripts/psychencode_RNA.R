
############# Script for testing enrichment of bipolar risk variants in C4A coexpression network ###########

library(openxlsx)
library(dplyr)
library(ggplot2)
library(tidyverse)  
library(data.table)
library(ggstatsplot)
library(ggrepel)
library(ggpubr)

### Read in Input RNA-seq data from http://resource.psychencode.org/
meta<- read.csv("/data/PEC_capstone_data_map_clinical.csv") #metadata for all samples
datExpr.reg<- as.data.frame(data.table::fread("/data/DER-01_PEC_Gene_expression_matrix_normalized.txt",header = TRUE)) # normalized counts matrix with controls + disorders (SCZ,BP,etc)

colnames(meta)[2] <- "ID"
#scz_names<- meta$ID[which(meta$diagnosis == "Schizophrenia")] #558 individuals with SCZ
bd_names <- meta$ID[which(meta$diagnosis == "Bipolar Disorder")] #216 individuals
hc_names<- meta$ID[which(meta$diagnosis == "Control")] #986 indiviudals
rownames(datExpr.reg) <- datExpr.reg$gene_id
datExpr.reg <- datExpr.reg[,-1]

#check if C4A is present in the data
grep("ENSG00000244731", rownames(datExpr.reg))

dat <- datExpr.reg[grep("ENSG00000244731", rownames(datExpr.reg)),] #subset data to only contain C4A
dat[2,] <- meta$diagnosis[match(colnames(dat), meta$ID)]
dat[3,] <- meta$sex[match(colnames(dat), meta$ID)]
dat[4,] <- meta$ethnicity[match(colnames(dat), meta$ID)]
dat[5,] <- meta$ageDeath[match(colnames(dat), meta$ID)]
dat[6,] <- meta$study[match(colnames(dat), meta$ID)]
rownames(dat) <- c("C4A", "Diagnosis", "sex", "ethnicity", "ageDeath", "study")

dat1 <- as.data.frame(t(dat))
dat2 <- dplyr::filter(dat1, Diagnosis %in% c("Control", "Bipolar Disorder"))
dat2$C4A <- as.numeric(dat2$C4A)
dat2[, 2:6] <- as.factor(dat2[, 2:6])

## Check for confounding factors and correct for it using combat if needed
model <- glm(Diagnosis ~ C4A + ageDeath + sex + ethnicity + study,
data = dat2, family = binomial(link = 'logit'))
summary(model)


##### Plot C4A mRNA levels across groups ###
ggbetweenstats(
data  = dat2,
x     = Diagnosis,
y     = C4A,
title = "Distribution of C4A expression"
)

res.aov <- aov(C4A ~ Diagnosis, data = dat2)
summary(res.aov)

TukeyHSD(res.aov)

########## Construct seeded co-expression network of C4A and test enrichment of GWAS hits for Bipolar disorder

## load functions required
OR <- function(q,k,m,t) {
# 2 x 2 table:
#         inTest   !inTest
# inRef     q        k
# !inRef    m        t
fisher.out <- fisher.test(matrix(c(q, k - q, m - q, t - m - k + q), 2, 2),
conf.int = TRUE)
OR <- fisher.out$estimate
pval <- fisher.out$p.value
upCI <- fisher.out$conf.int[1]
downCI <- fisher.out$conf.int[2]
output <- c(OR, pval, upCI, downCI)
names(output) <- c("OR", "Fisher p", "-95%CI", "+95%CI")
return(output)
}
# Count overlaps and run the analysis
ORA <- function(testpath, refpath, testbackground, refbackground) {
testpath = testpath[testpath %in% testbackground]
refpath = refpath[refpath %in% refbackground]
q <- length(intersect(testpath, refpath)) # overlapped pathway size
k <- length(intersect(refpath, testbackground)) # input gene set
m <- length(intersect(testpath, refbackground)) # input module
t <- length(intersect(testbackground, refbackground)) # Total assessed background (intersect reference and test backgrounds)
empvals <- OR(q, k, m, t)
tmpnames <- names(empvals)
empvals <- as.character(c(empvals, q, k, m, t, 100*signif(q/k, 3)))
names(empvals) <- c(tmpnames, "Overlap", "Reference List", "Input List", "Background", "% List Overlap")
return(empvals)
}

# Load pre-computed C4A network files from Kim et al., 2021
load("C4A-network.rdata")
c4a<- network

### Netowrk enrichment of genetic risk variants associated to disorders

#load risk genes from Mullins et al., 2021 for bipolar disorder

bp_risk <- openxlsx::read.xlsx("Bipolar_risk_genes_GWAS.xlsx")
bd1_risk <- openxlsx::read.xlsx("BD1_only_temp.xlsx")

# test bipolar disorder (type1) enrichment in C4A networks
ORA(refpath = c4a$Gene[c4a$`34.FDR` < .05 & c4a$`34.R` < 0], testpath = unique(bd1_risk), testbackground = c4a$Gene, refbackground = c4a$Gene)
 #                OR            Fisher p              -95%CI              +95%CI             Overlap      
#"0.499006473137796" "0.015018261815649" "0.254014290453813" "0.888951593509675"                "12"   
    
 ORA(refpath = c4a$Gene[c4a$`34.FDR` < .05 & c4a$`34.R` > 0], testpath = unique(bd1_risk), testbackground = c4a$Gene, refbackground = c4a$Gene)
 #                  OR               Fisher p                 -95%CI                 +95%CI                Overlap         
 #  "1.97030306421926" "0.000264810463952911"     "1.36106113723605"     "2.78647234755508"                   "39"                      


#load risk genes from Trubetskoy et al., 2022 for schizophrenia

scz_risk <- openxlsx::read.xlsx("Trubetskoy_extended_scz_GWAS.xlsx")

core_scz<- scz_risk$Symbol.ID[which(scz_risk$Extended.GWAS == "NO")]
all_scz<- scz_risk$Symbol.ID[which(scz_risk$Extended.GWAS == "YES")]
all_scz <- unique(all_scz)
prior_scz<- scz_risk$Symbol.ID[which(scz_risk$Prioritised == 1)] #prioritised gene list curated by PGC authors

# test scz enrichment in C4A networks
ORA(refpath = c4a$Gene[c4a$`34.FDR` < .05 & c4a$`34.R` < 0], testpath = unique(all_scz), testbackground = c4a$Gene, refbackground = c4a$Gene)
#                    OR               Fisher p                 -95%CI                 +95%CI                Overlap 
#    "1.78314230133616" "2.73443423567088e-05"     "1.36434398189782"      "2.3025822079709"                   "72" 


ORA(refpath = c4a$Gene[c4a$`34.FDR` < .05 & c4a$`34.R` > 0], testpath = unique(all_scz), testbackground = c4a$Gene, refbackground = c4a$Gene)
#                  OR             Fisher p               -95%CI               +95%CI              Overlap 
# "0.609323339806113" "0.0132355092663909"  "0.392940266731395"   "0.90675963901901"                 "26" 

ORA(refpath = c4a$Gene[c4a$`34.FDR` < .05 & c4a$`34.R` > 0], testpath = unique(prior_scz), testbackground = c4a$Gene, refbackground = c4a$Gene)
#                 OR            Fisher p              -95%CI              +95%CI             Overlap 
#"0.515341818234545"   "0.1694343233889" "0.163881510002303"  "1.24146510374463"                 "5" 


 ORA(refpath = c4a$Gene[c4a$`34.FDR` < .05 & c4a$`34.R` < 0], testpath = unique(prior_scz), testbackground = c4a$Gene, refbackground = c4a$Gene)
#                    OR               Fisher p                 -95%CI                 +95%CI                Overlap 
#    "2.97220752678201" "1.25721354426271e-05"     "1.82353725073903"     "4.68420315181716"                   "25" 


#### Plot the enrichment results
plot <- data.frame(row.names = c("Schizophrenia", "Bipolar Disorder"))
 #                p-val_positive   p-val_negative  log_positive  log_negative
#Schizophrenia    0.0132355092663909 2.734434e-05 1.878259 4.563133
#Bipolar Disorder 0.000264810463952911 0.015018261815649 3.577065 1.82338

#add pvalues from above to the data-frame and plot log10 values for each enrichment test
#add threshold line for bonferroni correction = (0.05/4)

p1<- ggplot(plot, mapping = aes(x=log_positive, y =rownames(plot))) + geom_bar(stat = "identity", fill="darkred") + theme_classic() + xlim(limits=c(0,7)) + geom_vline(xintercept = -log10(0.05/4),linetype = 'dashed', color = 'black')

p2<- ggplot(plot, mapping = aes(x=log_negative, y =rownames(plot))) + geom_bar(stat = "identity", fill="darkblue") + theme_classic() + xlim(limits=c(0,7)) + geom_vline(xintercept = -log10(0.05/4),linetype = 'dashed', color = 'black')

ggarrange(p2,p1)
