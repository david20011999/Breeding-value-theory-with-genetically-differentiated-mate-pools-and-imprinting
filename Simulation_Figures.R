# Run from the project root with Rscript Figures/Simulation_Figures.R.
# Requires David's AlphaSimR extension at the path passed to load_all below.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
figure_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else normalizePath("Figures")

# Figures 5 and 6

# ============================================================
# 1. LIBRARIES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ragg)

devtools::load_all("~/Estancia- gametic models/AlphaSimR/")

# This seed applies to future runs;
#the seed behind David's existing saved plotting data is unknown.
simulation_seed <- 20261108L
set.seed(simulation_seed)

nind <- 2000
nchr <- 10
segsites <- 1000
nqtlperchr <- 1000

founderPop = quickHaplo(nInd=nind, nChr=nchr, segSites=segsites)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = nqtlperchr, var = 1, gamma=TRUE,shape= 0.001)
trt = SP$traits[[1]]

qtlEffects = data.frame(
  marker = paste(rep(1:nchr, each = nqtlperchr), trt@lociLoc, sep = "_"),
  additiveEffect = trt@addEff,
  imprintingEffect = c(trt@addEff * sample(
    c(1, -1), size = trt@nLoci, replace = T
  ))[sample(
    1:length(trt@addEff),
    size = length(trt@addEff),
    replace = F
  )]
)

# Import in SimParam
SP$importTrait(markerNames = qtlEffects$marker,
               addEff = qtlEffects$additiveEffect,
               #domEff = qtlEffects$dominanceEffect,
               impEff = qtlEffects$imprintingEffect,
               name = "New_Trait")

SP$setSexes(sexes = "yes_sys")

a <- SP$traits[[2]]@addEff
i <- SP$traits[[2]]@impEff

#Create population
Pop = newPop(founderPop, simParam=SP)


Pop@ebv<- bvSS(Pop) #selectInd is selection using bv without impriting. bv() function use the impritning one.

Sel_Males_1 <- selectInd(Pop, 100,sex = "M", simParam=SP,use = "gv",trait = 2, selectTop = TRUE)
Sel_Males_2 <- selectInd(Pop, 100,sex = "M", simParam=SP,use = "ebv",trait = 2, selectTop = TRUE)

GenoMales_1 <- pullQtlGeno(Sel_Males_1,trait = 2, simParam=SP)
pM_1 <- apply(GenoMales_1 * 0.5, 2, mean, na.rm = T)

GenoMales_2 <- pullQtlGeno(Sel_Males_2,trait = 2, simParam=SP)
pM_2 <- apply(GenoMales_2 * 0.5, 2, mean, na.rm = T)

Sel_Females_1 <- selectInd(Pop, 100,sex = "F", simParam=SP,use = "gv",trait = 2, selectTop = FALSE)
Sel_Females_2 <- selectInd(Pop, 100,sex = "F", simParam=SP,use = "ebv",trait = 2, selectTop = TRUE)

GenoFemales_1 <- pullQtlGeno(Sel_Females_1,trait = 2, simParam=SP)
pF_1 <- apply(GenoFemales_1 * 0.5, 2, mean, na.rm = T)

GenoFemales_2 <- pullQtlGeno(Sel_Females_2,trait = 2, simParam=SP)
pF_2 <- apply(GenoFemales_2 * 0.5, 2, mean, na.rm = T)

df_1 <- data.frame(Dif=abs(pM_1 - pF_1),Value = abs(a),Chr=rep(1:nchr, each=nqtlperchr),Strategy="Antagonistic selection",Effect="Additive")
df_2 <- data.frame(Dif=abs(pM_1 - pF_1),Value = abs(i),Chr=rep(1:nchr, each=nqtlperchr),Strategy="Antagonistic selection",Effect="Imprinting")
df_3 <- data.frame(Dif=abs(pM_2 - pF_2),Value = abs(a),Chr=rep(1:nchr, each=nqtlperchr),Strategy="Synergistic selection",Effect="Additive")
df_4 <- data.frame(Dif=abs(pM_2 - pF_2),Value = abs(i),Chr=rep(1:nchr, each=nqtlperchr),Strategy="Synergistic selection",Effect="Imprinting")

df <- rbind(df_1,df_2,df_3,df_4)
df$Chr<-as.factor(df$Chr)

source(file.path(figure_dir, "figure_style.R"))
nQTL <- length(a)
stopifnot(length(i) == nQTL, nQTL == nchr * nqtlperchr)
effplot <- rbind(
  data.frame(QTL_position=seq_len(nQTL), Value=a,
             Chr=rep(seq_len(nchr),each=nqtlperchr), Effect="Additive"),
  data.frame(QTL_position=seq_len(nQTL), Value=i,
             Chr=rep(seq_len(nchr),each=nqtlperchr), Effect="Imprinting"))
f5_data <- list(effplot=effplot, nchr=nchr, nQTL=nQTL,
                chr_boundaries=seq_len(nchr-1)*nqtlperchr+.5,
                chr_midpoints=(seq_len(nchr)-.5)*nqtlperchr+.5)
# Keep both saved plotting datasets tied to the same run and code version.
simulation_metadata <- list(
  seed=simulation_seed, rng_kind=RNGkind(), generated_at=Sys.time(),
  script_md5=unname(tools::md5sum(file.path(figure_dir,"Simulation_Figures.R"))),
  session_info=sessionInfo())
f5_data$simulation_metadata <- simulation_metadata
f6_data <- list(df=df, simulation_metadata=simulation_metadata)
saveRDS(f5_data, file.path(figure_dir,"plot_F5_objects.rds"))
saveRDS(f6_data, file.path(figure_dir,"plot_F6_objects.rds"))
p1 <- plot_simulation_effects(f5_data)
p2 <- plot_selection_signatures(f6_data)
if (interactive()) { print(p1); print(p2) }
save_paper_plot(p1, "simulation_effects", paper_column_mm, 76)
save_paper_plot(p2, "simulation_results", paper_column_mm, 92)
