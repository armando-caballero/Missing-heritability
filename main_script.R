# Load libraries
set.seed(1234)
require("plyr")
require("magrittr")
require ("ggplot2")
require("ppcor")
require("tidyr")
require("tibble")
require("purrr")
install.packages("~/GWASmining_1.0.tar.gz",repos=NULL,type="source")
require("GWASmining")

# Load pre-processed catalog
load ("catalog_C50.2.rda")
catalog$trait.type = factor(catalog$trait.type, levels=sort(levels(catalog$trait.type)))
catalog$pmid = factor(catalog$pmid)

# Review PMIDs with effect.type NA
catalog %<>%
  dplyr::filter(pmid != 18193043) %>%
  dplyr::filter(pmid != 19448622) %>%
  dplyr::filter(pmid != 20686565) %>%
  dplyr::filter(pmid != 21102462) %>%
  dplyr::filter(pmid != 22139419) %>%
  dplyr::filter(pmid != 23599027)

catalog %<>%
  dplyr::filter(pmid != 22037903) # Too large h2

catalog[which(catalog$pmid==23263863),]$effect.type = "BETA"
catalog[which(catalog$pmid==24026423),]$effect.type = "BETA"
catalog[which(catalog$pmid==25705162),]$effect.type = "BETA"


# All SNPs in a trait should have the same effect.type
catalog = dplyr::filter(catalog, pmid != 28828242) # T2D
catalog$cluster = as.character(catalog$cluster)
catalog[which(catalog$trait=="Obesity"),]$ cluster = "Obesity"
catalog$cluster = as.factor(catalog$cluster)

# Filter by frequency
catalog %<>% dplyr::filter (!is.na(q))

# Remove PMIDs with excessive high BETA values
catalog = dplyr::mutate(catalog, va = 2*effect*effect*q*(1-q))
va_per_pmid_trait_beta = ddply(dplyr::filter(catalog,effect.type=="BETA"), .(pmid,cluster), function(x) sum(x$va, na.rm = TRUE))
invalid_pmids = dplyr::filter(va_per_pmid_trait_beta, V1>1)$pmid
catalog = catalog[-which(catalog$pmid %in% invalid_pmids),]
rm(va_per_pmid_trait_beta, invalid_pmids)

# Remove traits without minimum 3 PMIDs and 30 genes
tmp = ddply(catalog, .(cluster), summarise, N = length(unique(gene)), works = length(levels(factor(pmid))))
tmp = dplyr::filter(tmp, N >= 30 & works >= 3)
catalog = catalog[catalog$cluster %in% tmp$cluster,]
catalog = droplevels(catalog)
rm(tmp)
catalog$va = NULL


# Sort PMIDs by trait and N and update frequency and effect sizes
acum = ddply ( catalog, .(cluster), function(x) minimum_dataset (x, criteria="N", all.snps = TRUE, cumulative=TRUE))
acum = ddply(acum, .(cluster), function(x) update_parameters_trait(x))
acum = ddply(acum, .(cluster), function(x) update_parameters_trait(x,"q"))
acum = subset(acum, cumulative_genes >= 30)
acum = droplevels(acum)
acum %<>% dplyr::mutate (Va = 2*(effect^2)*(q)*(1-q), maf = ifelse(q>0.5, 1-q,q))

######################
eRR = function (OR, faa) {
  return (OR/(1+faa*(OR-1)))
}

efaa = function (k, q, RR1, RR2) {
  p = 1 - q
  p_aa = p*p
  p_Aa = 2*p*q
  p_AA = q*q
  return (k/(p_aa + p_Aa*0.5*RR1 + p_AA*RR2)) # additive effect on RR1
}

iterative_faa = function (k,q,OR) {
  faa = 0.5
  values = vector()
  count = 1
  repeat {
    faa = efaa(k,q,eRR(OR,faa),eRR(OR,faa))
    values[count] = faa
    if (length(values) == 10 & sd(values)==0) {
      break
    } else {
      count = count + 1
      if (count > 10) {
        count = 1
      }
    }
  }
  return (mean(values))
}

get_Vg = function (k,q,OR) {
  
  p = 1 -q
  p_aa = p*p
  p_Aa = 2*p*q
  p_AA = q*q
  
  faa = iterative_faa(k,q,OR)
  RR = eRR (OR, faa)
  faa = k/(p_aa + p_Aa*RR + p_AA*RR)
  fAa= RR*faa
  fAA= RR*faa
  
  mu_aa=0
  T = qnorm(1-faa)
  mu_Aa = T-qnorm(1-fAa)
  mu_AA = T-qnorm(1-fAA)
  mu_all= p_Aa*mu_Aa + p_AA*mu_AA
  
  Vg = p_aa * (mu_aa-mu_all)^2 + p_Aa*(mu_Aa-mu_all)^2 + p_AA*(mu_AA-mu_all)^2
  return (Vg/(1+Vg)) 
}

mean.tnorm3 = function(x) {
  dnorm(x)/pnorm(x,lower.tail=FALSE)
}

# OR effect are converted into BETA
acum[which(acum$cluster =="Basal cell carcinoma"),]$Va = data.frame(acum[which(acum$cluster =="Basal cell carcinoma"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.0034264,q,effect)))$Vg
acum[which(acum$cluster =="Colorectal cancer"),]$Va = data.frame(acum[which(acum$cluster =="Colorectal cancer"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.0041,q,effect)))$Vg
acum[which(acum$cluster =="Coronary heart disease"),]$Va = data.frame(acum[which(acum$cluster =="Coronary heart disease"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.05468181,q,effect)))$Vg
acum[which(acum$cluster =="Digestive disease"),]$Va = data.frame(acum[which(acum$cluster =="Digestive disease"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.003205,q,effect)))$Vg
acum[which(acum$cluster =="Obesity"),]$Va = data.frame(acum[which(acum$cluster =="Obesity"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.2,q,effect)))$Vg
acum[which(acum$cluster =="Parkinson's disease"),]$Va = data.frame(acum[which(acum$cluster =="Parkinson's disease"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.002,q,effect)))$Vg
acum[which(acum$cluster =="Primary biliary cholangitis"),]$Va = data.frame(acum[which(acum$cluster =="Primary biliary cholangitis"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.00025555,q,effect)))$Vg
acum[which(acum$cluster =="Prostate cancer"),]$Va = data.frame(acum[which(acum$cluster =="Prostate cancer"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.14,q,effect)))$Vg
acum[which(acum$cluster =="Psoriasis"),]$Va = data.frame(acum[which(acum$cluster =="Psoriasis"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.011,q,effect)))$Vg
acum[which(acum$cluster =="Rheumatoid arthritis"),]$Va = data.frame(acum[which(acum$cluster =="Rheumatoid arthritis"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.01,q,effect)))$Vg
acum[which(acum$cluster =="Schizophrenia"),]$Va = data.frame(acum[which(acum$cluster =="Schizophrenia"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.0033,q,effect)))$Vg
acum[which(acum$cluster =="Systemic lupus erythematosus"),]$Va = data.frame(acum[which(acum$cluster =="Systemic lupus erythematosus"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.0004,q,effect)))$Vg
acum[which(acum$cluster =="Testicular germ cell tumor"),]$Va = data.frame(acum[which(acum$cluster =="Testicular germ cell tumor"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.005,q,effect)))$Vg
acum[which(acum$cluster =="Type 1 diabetes"),]$Va = data.frame(acum[which(acum$cluster =="Type 1 diabetes"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.0033,q,effect)))$Vg
acum[which(acum$cluster =="Type 2 diabetes"),]$Va = data.frame(acum[which(acum$cluster =="Type 2 diabetes"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.08,q,effect)))$Vg
acum[which(acum$cluster =="Ulcerative colitis"),]$Va = data.frame(acum[which(acum$cluster =="Ulcerative colitis"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.00377,q,effect)))$Vg
acum[which(acum$cluster =="Vitiligo"),]$Va = data.frame(acum[which(acum$cluster =="Vitiligo"),] %>% dplyr::rowwise() %>% dplyr::mutate (Vg = get_Vg(0.012,q,effect)))$Vg
acum[which(acum$effect.type=="OR"),]$effect = sqrt(acum[which(acum$effect.type=="OR"),]$Va / (2.0 *(1-acum[which(acum$effect.type=="OR"),]$q)*acum[which(acum$effect.type=="OR"),]$q ))

# Summarise cumulated data
lnorm = ddply (acum, .(cluster,count), function(x)
  if (x$effect.type =="BETA") {
    distribution(x$effect, mode = "analytic", "lnorm")
  } else if (x$effect.type =="OR") {
    distribution(x$effect, mode = "analytic", "lnorm")
  })
lnorm$parameter1 = as.numeric(as.character(lnorm$parameter1))
lnorm$parameter2 = as.numeric(as.character(lnorm$parameter2))
x = ddply (acum, .(cluster,count), summarise, effect.type = levels(droplevels(effect.type)))
x2 = ddply (acum, .(cluster,count), summarise, effect.mean =mean(effect))
xq = ddply (acum, .(cluster,count), summarise, q.mean =mean(q,na.rm = TRUE))
x3 = ddply (acum, .(cluster,count), summarise, trait.type = levels(droplevels(trait.type)))
x4 = ddply (acum, .(cluster,count), summarise, ngenes =mean(cumulative_genes))
x5 = ddply (acum, .(cluster,count), summarise, N =mean(cumulative_N))
x6 = ddply (acum, .(cluster,count), summarise, Va =sum(Va,na.rm = TRUE))
x7 = ddply (acum, .(cluster,count), summarise, effect.median =median(effect))
x8 = ddply (acum, .(cluster,count), summarise, effect.sd =sd(effect))

lnorm = merge(lnorm,x) %>% merge(x2) %>% merge (xq) %>% merge(x3) %>% merge(x4) %>% merge(x5) %>% merge(x7) %>% merge(x8) %>% merge(x6)

# Final dataset
final = ddply(acum,.(cluster), function(x) x[which(x$count==max(x$count)),]  )