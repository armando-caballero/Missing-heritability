library("magrittr")
library("plyr")
install.packages("~/GWASmining_1.0.tar.gz",repos=NULL,type="source")
require("GWASmining")
library("ggplot2")

# To use with "acum" dataset

sim_h2_obs = function (TRAIT, data_acum, add_genes = 0, total_genes = 0, sample_genes = 0, dist = "lnorm") {
  
  NBOOT = 10000
  data_acum %<>% dplyr::filter ( cluster==TRAIT)

  # Acum trait y modelo "q"
  acum_trait = dplyr::filter(data_acum, count == max(count))

  # data3
  lnorm = ddply (data_acum, .(cluster,count), function(x) distribution(x$effect, mode = "analytic", dist))
  lnorm$parameter1 = as.numeric(as.character(lnorm$parameter1))
  lnorm$parameter2 = as.numeric(as.character(lnorm$parameter2))
  x2 = ddply (data_acum, .(cluster,count), summarise, effect.mean =mean(effect))
  xq = ddply (data_acum, .(cluster,count), summarise, q.mean =mean(q,na.rm = TRUE))
  x4 = ddply (data_acum, .(cluster,count), summarise, ngenes =mean(cumulative_genes))

  lnorm = merge(lnorm,x2) %>% merge (xq) %>% merge(x4)
  rm(x2,x4,xq)
  
  lnorm_maf = ddply(acum, .(cluster, count), function(x) distribution(x$maf, mode="analytic", dist="norm"))
  lnorm_maf = dplyr::select (lnorm_maf, -dist, -kstest)
  colnames(lnorm_maf) = c("cluster","count","mean_maf","sd_maf")
  lnorm = merge(lnorm, lnorm_maf); rm (lnorm_maf)
  lnorm$mean_maf = as.numeric(as.character(lnorm$mean_maf))
  lnorm$sd_maf = as.numeric(as.character(lnorm$sd_maf))

  MIN = 3
  tbl = table(lnorm$cluster)
  df3 = droplevels(lnorm[lnorm$cluster %in% names(tbl)[tbl >= MIN],,drop=FALSE])

  print("par1")
  par1 = dlply(df3, .(cluster), function(x) nls (data=x, parameter1 ~ a*ngenes^b, start = list(a=0.15,b=0.23))) %>% ldply(coef) # beta -0.15 0.23 # gamma -0.15,0.23
  print("par2")
  if (dist!="exp") par2 = dlply(df3, .(cluster), function(x) nls (data=x, parameter2 ~ a*ngenes^b, start = list(a=0.15,b=-0.023))) %>% ldply(coef)
  print("maf1")
  par1_maf = dlply(df3, .(cluster), function(x) nls (data=x, mean_maf ~ a*ngenes^b, start = list(a=-0.15,b=0.23))) %>% ldply(coef)
  print("maf2")
  par2_maf = dlply(df3, .(cluster), function(x) nls (data=x, sd_maf ~ a*ngenes^b, start = list(a=0.15,b=-0.23))) %>% ldply(coef)
  
  mexp = function (a, b, x) {
    a * x ^b
  }

  if (total_genes == 0) {
    N_genes = nrow(acum_trait) + add_genes
  } else {
    N_genes = total_genes
  }
  
  if (sample_genes == 0) {
    sample_genes = N_genes
  } else {
    print("par1")
    print(mexp(par1$a, par1$b, N_genes))
    if (dist != "exp") {
      print("par2")
      print(mexp (par2$a, par2$b, N_genes))
    }
  }

  par1 = mexp(par1$a, par1$b, N_genes)
  if (dist != "exp") par2 = mexp (par2$a, par2$b, N_genes)
  par1_maf = mexp(par1_maf$a, par1_maf$b, N_genes)
  par2_maf = mexp (par2_maf$a, par2_maf$b, N_genes)

  result = list()
  result$N_genes = N_genes
  result$sample_genes = sample_genes
  result$h2 = vector(mode = "double",length = NBOOT)
  result$par1 = par1
  if (dist != "exp") result$par2 = par2
  result$mean_maf = par1_maf
  result$sd_maf = par2_maf
  
  for (i in 1:NBOOT){
    #effects = rlnorm (n = N_genes, meanlog = lnorm2$parameter1, sdlog = lnorm2$parameter2 )
    if (dist=="lnorm") effects = rlnorm (n = sample_genes, meanlog = par1, sdlog = par2)
    else if (dist=="beta") effects = rbeta (n = sample_genes, shape1 = par1, shape2 = par2)
    else if (dist=="gamma")  effects = rgamma (n = sample_genes, shape = par1, rate = par2)
    else if (dist=="exp")  effects = rexp (n = sample_genes, rate = par1)
    else print("ERORRRRR")

    mafs = rnorm (n = sample_genes, mean = par1_maf, sd = par2_maf)
    #q = mexp (fit_q[1], fit_q[2], effects)
    #h2 = sum(effects * effects* 2 * q * (1-q))
    mafs = ifelse(mafs > 0.5, 0.5, mafs)
    mafs = ifelse(mafs < 0.0, 0.0, mafs)
    h2 = sum(effects * effects* 2 * mafs * (1- mafs))
    result$h2[i] = h2
    
  }

  return(result)

}

# Running example
set.seed(1234)
TRAIT = "Prostate cancer"
df = data.frame()
MIN = 106
MAX = 400
min_x = MIN
max_x = MAX
step = 1
all = list()
count = 1
DIST = "gamma"

for (i in seq(min_x,max_x,step)) {
  print(paste(i," / ", max_x, sep=""))
  #all[[count]] = sim_h2_obs(TRAIT, acum, total_genes = i)
  h2 = sim_h2_obs(TRAIT, acum, total_genes = i, dist = DIST)
  h2 = cbind (h2$h2, i)
  df = rbind(df, h2)
  count = count + 1
}

colnames (df) = c("h2","i")
matrix = ddply(df, .(i), summarise, median_h2 = median(h2), h2max = quantile(h2,probs = 0.95), h2min=quantile(h2, probs=0.05))
load("df_h2.rda")
df_h2 = dplyr::filter(df_h2, cluster==TRAIT)
df_h2 %<>% dplyr::mutate (meanh2 = (minh2+maxh2)/2)
test = dplyr::mutate (matrix, best_mean = abs(median_h2 - df_h2$meanh2))
test[which(test$best_mean == min(test$best_mean)),]

PLOT = ggplot(data=df, aes(y=h2, x=factor(i)))
PLOT = PLOT + geom_boxplot(fill="lightblue", outlier.alpha = 0.01, outlier.colour = "lightblue")
PLOT = PLOT + geom_hline(yintercept = df_h2$h2C, colour="blue", linetype = "dashed") 
PLOT = PLOT + geom_hline(yintercept = df_h2$minh2, colour="red", size=0.5)
PLOT = PLOT + geom_hline(yintercept = mean(c(df_h2$minh2,df_h2$maxh2)), colour="red", size=0.5, linetype="dashed")
PLOT = PLOT + geom_hline(yintercept = df_h2$maxh2, colour="red", size=0.5)
PLOT = PLOT + theme_bw() + ylim(0,2) + xlab ("N genes")
PLOT = PLOT + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle=90), axis.title.y = element_text(angle = 0, vjust = 0.5))
PLOT = PLOT + scale_x_discrete(breaks = c(min_x, seq(max_x/10,max_x,max_x/10)))
PLOT

best_result = test[which(test$best_mean == min(test$best_mean)),]
mejor_h2 = best_result$median_h2

asciende = TRUE
if (mejor_h2 > test[nrow(test),]$median_h2) asciende = FALSE
max_genes = MAX
if (asciende) {
  if (length(test[which(test$h2min <= mejor_h2),]$i)>0) max_genes = max(test[which(test$h2min <= mejor_h2),]$i)
}
if (!asciende) {
  if (length(test[which(test$h2max >= mejor_h2),]$i)>0) max_genes = max(test[which(test$h2max >= mejor_h2),]$i)
}

asciende = TRUE
if (mejor_h2 < test[1,]$median_h2) asciende = FALSE
min_genes = MIN
if (asciende) {
  if (length(test[which(test$h2max >= mejor_h2),]$i)>0) min_genes = min(test[which(test$h2max >= mejor_h2),]$i)
}
if (!asciende) {
  if (length(test[which(test$h2min <= mejor_h2),]$i)>0) min_genes = min(test[which(test$h2min <= mejor_h2),]$i)
}
min_genes
max_genes

best_N = matrix %>% dplyr::mutate(diff = median_h2 - df_h2$minh2)
best_N = best_N %>%  dplyr::arrange(abs(diff)) 
best_N = head(best_N,1)$i
best_N
quantile(sim_h2_obs(TRAIT, acum, total_genes = best_N, sample_genes = df_h2$known_genes)$h2, probs = c(0.025, 0.5, 0.975))
matrix[which(matrix$medianh2 < 0.21),] %>% dplyr::filter(i == max(i))
