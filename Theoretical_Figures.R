# Regenerate the manuscript genotype-colour Figures 2-3 and Figure 3 supplement:
# Rscript Figures/Theoretical_Figures.R --genotype-colours
# Figure 3 only: Rscript Figures/Figure3_Variants.R
# Blue genotype-shading backups: Rscript Figures/Theoretical_Figures.R --alpha-only
# Contains Figures 2-4; simulation and cattle plots have their own scripts.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
figure_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else normalizePath("Figures")

# Regenerate only Figure 4 without touching the accepted Figures 2-3.
if ("--delta-only" %in% commandArgs(trailingOnly = TRUE)) {
  source_lines <- readLines(file.path(figure_dir, "Theoretical_Figures.R"))
  first <- match("# Figure 4", source_lines)
  stopifnot(!is.na(first))
  eval(parse(text = source_lines[(first + 1L):length(source_lines)]))
  quit(save = "no")
}

# Retain the blue genotype-shading backups as the default script output.
genotype_colours <- "--genotype-colours" %in% commandArgs(trailingOnly = TRUE)
alpha_figure_height <- 120

alpha_stem <- function(stem) {
  if (genotype_colours) paste0(stem, "_genotype_colours") else stem
}

# Five discrete Set1 colours, arranged in hue order across both figures:
# purple -> pink -> red -> orange -> yellow. Pink/red/orange link the hets.
# Use a fixed subset of the NINE-colour palette; independently requesting
# three/four colours would change the shared mapping.
# https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
# This is an ordered categorical key, not a monotonic lightness scale.
# Avoid the blue/purple and red/green pairs that are difficult under common
# red-green colour-vision deficiencies. Simulations favour this subset over
# the blue/purple alternatives; pink/orange remain closer under tritanopia.
# Positions: homozygote 0/0, het d-i, pooled het d, het d+i, homozygote 1/1.
# Figure 2 uses positions 1/3/5; Figure 3 uses 1/2/4/5. These are category
# positions, not a continuous mapping of genetic value (d+i can exceed a).
# Genotypes use the manuscript's maternal/paternal allele order.
# Colour follows the individual whose Value is plotted: the parent above,
# the progeny below. Genotype remains the original x-axis predictor.
genotype_palette <- setNames(
  RColorBrewer::brewer.pal(9, "Set1")[c(4, 8, 1, 5, 6)],
  c("0/0", "1/0", "het", "0/1", "1/1")
)
# A fixed blue shade identifies each genotype in the backup figures:
# ref/ref (0/0) is darkest, alt/alt (1/1) is lightest in every panel.
backup_genotype_palette <- setNames(
  RColorBrewer::brewer.pal(9, "Blues")[c(9, 7, 5, 3, 1)],
  names(genotype_palette)
)
alpha_fill_scale <- function(imprinting = FALSE) {
  keys <- if (imprinting) c("0/0", "1/0", "0/1", "1/1") else c("0/0", "het", "1/1")
  labels <- if (imprinting) keys else c("0/0", "0/1, 1/0", "1/1")
  palette <- if (genotype_colours) genotype_palette else backup_genotype_palette
  scale_fill_manual(
    values = palette[keys], breaks = keys, labels = labels,
    name = if (imprinting) "Ordered genotype\n(mat/pat)" else "Genotype",
    guide = guide_legend(order = 1,
                        theme = if (imprinting) theme(legend.title = element_text(hjust = 0.5)) else NULL,
                        override.aes = list(shape = 21, size = 4, alpha = 1))
  )
}
alpha_size_guide <- function() {
  guide_legend(order = 2,
               override.aes = list(shape = 21, fill = "grey75", alpha = 1))
}

parent_fit_label <- expression(mu[G]+A[j])
breeding_fit_label <- expression(mu[G[plain(progeny)]]+B[j*","*s])
# Show only the fitted-value diamond and breeding-value square in the legend.
# Keep all line types in the plots, without a line legend.
alpha_model_scales <- function(imprinting = FALSE) {
  within <- if (imprinting) c("mat", "pat") else "within"
  line_keys <- c(within, "B")
  line_types <- c(within="solid", mat="solid", pat="dashed",
                  B="solid")
  list(
    scale_shape_manual(values=c(A=23, B=22), limits=c("A", "B"), breaks=c("A", "B"),
      labels=c(parent_fit_label, breeding_fit_label), name=NULL,
      guide=guide_legend(order=3)),
    scale_linetype_manual(values=line_types[line_keys], limits=line_keys, guide="none")
  )
}

# Render PNGs from the vector master: this also avoids lost nested-facet
# headings seen with the local ragg device. Ghostscript must be on PATH.
render_alpha_png <- function(stem) {
  gs <- Sys.which("gs")
  if (!nzchar(gs)) stop("Ghostscript (gs) is required to render the alpha-plot PNGs.")
  status <- system2(gs, c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
                         "-sDEVICE=png16m", "-r300",
                         paste0("-sOutputFile=", shQuote(file.path(figure_dir, paste0(stem, ".png")))),
                         shQuote(file.path(figure_dir, paste0(stem, ".pdf")))))
  if (status != 0L) stop("PNG rendering failed for ", stem)
}

#Figure 2


# ============================================================
# 1. LIBRARIES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(scales)
library(scico)
library(ragg)


# ============================================================
# 2. PARAMETERS
# ============================================================

pPop <- 0.6
pMales <- 0.9
pFemales <- 0.3

qPop <- 1 - pPop
qMales <- 1 - pMales
qFemales <- 1 - pFemales

a <- 1
d <- 0.75
i <- 0


# ============================================================
# 3. FUNCTIONS TO GENERATE PROBABILITY DISTRIBUTIONS
# ============================================================

# Within-generation regression

make_within <- function(p, a, d) {
  
  q <- 1 - p
  
  data.frame(
    Genotype = c(2, 1, 0),
    ValueGenotype = c("1/1", "het", "0/0"),
    Value = c(a, d, -a),
    probability = c(
      p^2,
      2 * p * q,
      q^2
    )
  )
}


# Between-generation regression

make_between <- function(sex, p_male, p_female, a, d, i) {
  
  qm <- 1 - p_male
  qf <- 1 - p_female
  
  if (sex == "Male") {
    
    data.frame(
      Genotype = c(2, 2, 1, 1, 1, 1, 0, 0),
      ValueGenotype = c("1/1", "het", "1/1", "het",
                        "het", "0/0", "0/0", "het"),
      
      Value = c(
        a,
        d + i,
        a,
        d + i,
        d - i,
        -a,
        -a,
        d - i
      ),
      
      probability = c(
        p_male^2 * p_female,
        p_male^2 * qf,
        p_male * qm * p_female,
        p_male * qm * qf,
        p_male * qm * p_female,
        p_male * qm * qf,
        qm^2 * qf,
        qm^2 * p_female
      )
    )
    
  } else {
    
    data.frame(
      Genotype = c(2, 2, 1, 1, 1, 1, 0, 0),
      ValueGenotype = c("1/1", "het", "1/1", "het",
                        "het", "0/0", "0/0", "het"),
      
      Value = c(
        a,
        d - i,
        a,
        d - i,
        d + i,
        -a,
        -a,
        d + i
      ),
      
      probability = c(
        p_female^2 * p_male,
        p_female^2 * qm,
        p_female * qf * p_male,
        p_female * qf * qm,
        p_female * qf * p_male,
        p_female * qf * qm,
        qf^2 * qm,
        qf^2 * p_male
      )
    )
  }
}


# ============================================================
# 4. GENERATE ALL EXPERIMENTAL SCENARIOS
# ============================================================

scenarios <- list()

for (sex in c("Male", "Female")) {
  
  for (mating in c("Classical benchmark", "Differentiated mate pools")) {
    
    for (method in c(
      "Within-generation (w)",
      "Between-generation (b)"
    )) {
      
      # Select allele frequencies
      
      if (mating == "Classical benchmark") {
        
        pm <- pPop
        pf <- pPop
        
      } else {
        
        pm <- pMales
        pf <- pFemales
        
      }
      
      # Generate distributions
      
      if (method == "Within-generation (w)") {
        
        p <- if (sex == "Male") pm else pf
        
        tmp <- make_within(
          p = p,
          a = a,
          d = d
        )
        
      } else {
        
        tmp <- make_between(
          sex = sex,
          p_male = pm,
          p_female = pf,
          a = a,
          d = d,
          i = i
        )
      }
      
      # Add scenario information
      
      tmp$Sex <- sex
      tmp$Method <- method
      tmp$Mating <- mating
      
      scenarios[[paste(sex, method, mating, sep = "_")]] <- tmp
      
    }
  }
}

def_lm <- bind_rows(scenarios)


# ============================================================
# 5. VALIDATE PROBABILITY DISTRIBUTIONS
# ============================================================

probability_check <- def_lm %>%
  group_by(Sex, Method, Mating) %>%
  summarise(
    total_probability = sum(probability),
    .groups = "drop"
  )

print(probability_check)

stopifnot(
  all(abs(probability_check$total_probability - 1) < 1e-10)
)


# ============================================================
# 6. CALCULATE WEIGHTED REGRESSION COEFFICIENTS
# ============================================================

alpha_results <- lapply(
  split(
    def_lm,
    interaction(
      def_lm$Sex,
      def_lm$Method,
      def_lm$Mating,
      drop = TRUE
    )
  ),
  
  function(dat) {
    
    # Between-generation regression uses 2 * Value
    
    response <- if (
      dat$Method[1] == "Between-generation (b)"
    ) {
      2 * dat$Value
    } else {
      dat$Value
    }
    
    # Weighted regression
    
    fit <- lm(
      response ~ Genotype,
      data = dat,
      weights = probability
    )
    
    data.frame(
      Sex = dat$Sex[1],
      Method = dat$Method[1],
      Mating = dat$Mating[1],
      
      regression = unname(
        coef(fit)["Genotype"]
      ),
      
      intercept = unname(
        coef(fit)["(Intercept)"]
      ),
      
      response_mean = weighted.mean(
        response,
        dat$probability
      )
    )
  }
)

alphas <- bind_rows(alpha_results)


# ============================================================
# 7. GENERATE REGRESSION LINES ON THE GENETIC-VALUE SCALE
# ============================================================

line_df <- alphas %>%
  crossing(
    Genotype = seq(0, 2, length.out = 100)
  ) %>%
  mutate(
    fitted_value = (intercept + regression * Genotype) /
      if_else(Method == "Between-generation (b)", 2, 1)
  )


fitted_df <- alphas %>%
  filter(Method == "Within-generation (w)") %>%
  crossing(Genotype = 0:2) %>%
  mutate(fitted_value = intercept + regression * Genotype)

# Shift B by the progeny mean: mu_Gprogeny+B = 2*E(G_progeny)-mu_Gprogeny.
breeding_df <- line_df %>%
  filter(Method == "Between-generation (b)") %>%
  mutate(fitted_value = 2*fitted_value-response_mean/2)
square_df <- alphas %>%
  filter(Method == "Between-generation (b)") %>%
  crossing(Genotype=0:2) %>%
  mutate(fitted_value=intercept+regression*Genotype-response_mean/2)

# ============================================================
# 8. AGGREGATE POINT PROBABILITIES
# ============================================================

# Combine observations sharing the same genotype and value

plot_df <- def_lm %>%
  group_by(
    Sex,
    Method,
    Mating,
    Genotype,
    Value,
    ValueGenotype
  ) %>%
  summarise(
    probability = sum(probability),
    .groups = "drop"
  )


# ============================================================
# 9. FORMAT REGRESSION COEFFICIENT ANNOTATIONS
# ============================================================

alphas <- alphas %>%
  mutate(
    alpha_subscript = if_else(Sex == "Female", "f", "m"),
    alpha_superscript = if_else(Method == "Within-generation (w)", "w", "b"),
    alpha_label = if_else(Method == "Between-generation (b)",
      sprintf("paste(alpha[%s]^b, ' = ', 2*b[%s], ' = %.2f')",
              alpha_subscript, alpha_subscript, regression),
      sprintf("alpha[%s]^w == '%.2f'", alpha_subscript, regression))
  )


# ============================================================
# 10. FACTOR ORDER AND FACET LABELS
# ============================================================

old_methods <- c("Within-generation (w)", "Between-generation (b)")
new_methods <- c("Within-generation", "Between-generation")
for (frame_name in c("plot_df", "line_df", "fitted_df", "breeding_df", "square_df", "alphas")) {
  dat <- get(frame_name)
  dat$Method <- factor(dat$Method, levels = old_methods, labels = new_methods)
  dat$Mating <- factor(dat$Mating,
    levels = c("Classical benchmark", "Differentiated mate pools"),
    labels = c("Undifferentiated mate pools", "Differentiated mate pools"))
  dat$Sex <- factor(dat$Sex, levels = c("Female", "Male"),
                    labels = c("Female", "Male"))
  assign(frame_name, dat)
}

# ============================================================
# 11. PUBLICATION-QUALITY FIGURE
# ============================================================

p <- ggplot() +
  
  # ----------------------------------------------------------
# Probability distribution
# ----------------------------------------------------------

geom_point(
  data = plot_df,
  aes(
    x = Genotype,
    y = Value,
    fill = ValueGenotype,
    size = probability
  ),
  shape = 21,
  colour = "grey20",
  stroke = 0.25,
  alpha = 1
) +
  
  # ----------------------------------------------------------
# Weighted regression lines
# ----------------------------------------------------------

geom_line(data=subset(line_df,Method=="Within-generation"),
  aes(Genotype,fitted_value,linetype="within"),colour="black",linewidth=.65) +
geom_line(data=breeding_df,aes(Genotype,fitted_value,linetype="B"),colour="black",linewidth=.65) +

# Fixed-size white diamonds distinguish fitted values from genotype bubbles.
geom_point(
  data = fitted_df,
  aes(x = Genotype, y = fitted_value, shape = "A"),
  size = 2.2,
  fill = "white",
  colour = "black",
  stroke = 0.55
) +
geom_point(data=square_df,aes(Genotype,fitted_value,shape="B"),
  size=2.2,fill="white",colour="black",stroke=.55) +
  
# ----------------------------------------------------------
# Regression coefficients
# ----------------------------------------------------------

geom_text(
  data = subset(alphas, !(Method=="Between-generation" &
                          Mating=="Differentiated mate pools" & Sex=="Male")),
  aes(x=1, y=-1.65, label=alpha_label),
  parse=TRUE, vjust=.5, size=3.2, colour="black"
) +
geom_label(
  data = subset(alphas, Method=="Between-generation" &
                         Mating=="Differentiated mate pools" & Sex=="Male"),
  aes(x=2.3, y=-1.65, label=alpha_label),
  parse=TRUE, hjust=1, vjust=.5, size=3.2, colour="black",
  fill=scales::alpha("white", .75), linewidth=0,
  label.padding=grid::unit(.03,"lines"), label.r=grid::unit(0,"lines"),
  show.legend=FALSE
) +
  
  # ----------------------------------------------------------
# Nested facets
# ----------------------------------------------------------

ggh4x::facet_nested(
  Method ~ Mating + Sex
) +
  
  # ----------------------------------------------------------
# Fill identifies genotype; bubble area encodes frequency.
alpha_fill_scale() +
alpha_model_scales() +
  
  # Bubble size: area proportional to probability
# ----------------------------------------------------------

scale_size_area(
  max_size = 8,
  limits = c(0, 0.9),
  breaks = c(0.10, 0.30, 0.50, 0.70, 0.90),
  labels = scales::label_percent(accuracy = 1),
  name = "Frequency",
  guide = alpha_size_guide()
) +
  
  # ----------------------------------------------------------
# X-axis
# ----------------------------------------------------------

scale_x_continuous(
  breaks = c(0, 1, 2),
  limits = c(-0.4, 2.4),
  expand = expansion(mult = 0)
) +
  
  # ----------------------------------------------------------
# Y-axis
# ----------------------------------------------------------

scale_y_continuous(
  limits = c(-2.45, 2.05),
  breaks = -2:2,
  expand = expansion(mult = 0)
) +
  
  labs(
    x = expression("Allele dosage, " * N[1]),
    y = "Genetic value"
  ) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 11,
  base_family = "sans"
) +
  
  theme(
    
    # Background grid
    
    panel.grid.major.x = element_line(
      colour = "#E5E5E5",
      linewidth = 0.3
    ),
    
    panel.grid.major.y = element_line(
      colour = "#E5E5E5",
      linewidth = 0.3
    ),
    
    panel.grid.minor = element_blank(),
    
    # Axes
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    axis.ticks.length = grid::unit(2, "mm"),
    
    axis.text = element_text(
      colour = "black",
      size = 9
    ),
    
    axis.title = element_text(
      colour = "black",
      size = 10
    ),
    
    # Facets
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      colour = "black",
      face = "bold",
      size = 9.5
    ),
    
    panel.spacing.x = grid::unit(3, "mm"),
    panel.spacing.y = grid::unit(4, "mm"),
    
    # Legends
    
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8.5
    ),
    
    legend.key.width = grid::unit(8, "mm"),
    legend.spacing.x = grid::unit(5, "mm"),
    legend.spacing.y = grid::unit(0, "mm"),
    
    # General
    
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    plot.margin = margin(
      3, 8, 3, 5,
      unit = "mm"
    )
  )


# ============================================================
# 12. DISPLAY
# ============================================================

if (interactive()) print(p)


# ============================================================
# 13. EXPORT
# ============================================================

# Double-column PNG

# PNG is rendered from the vector PDF below.


# Vector PDF

ggsave(
  filename = file.path(figure_dir, paste0(alpha_stem("alpha_AD_def"), ".pdf")),
  plot = p,
  width = 180,
  height = alpha_figure_height,
  units = "mm",
  device = grDevices::pdf,
  useDingbats = FALSE,
  bg = "white"
)
render_alpha_png(alpha_stem("alpha_AD_def"))


#Figure 3

# ============================================================
# 1. LIBRARIES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(viridis)
library(scales)
library(ragg)

# ============================================================
# 2. PARAMETERS
# ============================================================

pPop <- 0.6
pMales <- 0.9
pFemales <- 0.3

a <- 1
d <- 0.75
i <- 0.5

qPop <- 1 - pPop
qMales <- 1 - pMales
qFemales <- 1 - pFemales

# ============================================================
# 3. FUNCTIONS TO GENERATE THE PROBABILITY DISTRIBUTIONS
# ============================================================

# Within-generation regression

make_within <- function(origin, p, a, d, i) {
  origin <- match.arg(origin, c("mat", "pat"))
  mat <- c(0, 0, 1, 1)
  pat <- c(0, 1, 0, 1)
  q <- 1 - p
  data.frame(
    Genotype = if (origin == "mat") mat else pat,
    ValueGenotype = c("0/0", "0/1", "1/0", "1/1"),
    Value = c(-a, d + i, d - i, a),
    probability = c(q^2, q * p, p * q, p^2)
  )
}


# Between-generation regression

make_between <- function(sex, p_male, p_female, a, d, i) {
  
  qm <- 1 - p_male
  qf <- 1 - p_female
  
  if (sex == "Male") {
    
    data.frame(
      Genotype = c(2, 2, 1, 1, 1, 1, 0, 0),
      ValueGenotype = c("1/1", "0/1", "1/1", "0/1",
                        "1/0", "0/0", "0/0", "1/0"),
      
      Value = c(
        a,
        d + i,
        a,
        d + i,
        d - i,
        -a,
        -a,
        d - i
      ),
      
      probability = c(
        p_male^2 * p_female,
        p_male^2 * qf,
        p_male * qm * p_female,
        p_male * qm * qf,
        p_male * qm * p_female,
        p_male * qm * qf,
        qm^2 * qf,
        qm^2 * p_female
      )
    )
    
  } else {
    
    data.frame(
      Genotype = c(2, 2, 1, 1, 1, 1, 0, 0),
      ValueGenotype = c("1/1", "1/0", "1/1", "1/0",
                        "0/1", "0/0", "0/0", "0/1"),
      
      Value = c(
        a,
        d - i,
        a,
        d - i,
        d + i,
        -a,
        -a,
        d + i
      ),
      
      probability = c(
        p_female^2 * p_male,
        p_female^2 * qm,
        p_female * qf * p_male,
        p_female * qf * qm,
        p_female * qf * p_male,
        p_female * qf * qm,
        qf^2 * qm,
        qf^2 * p_male
      )
    )
  }
}


# ============================================================
# 4. GENERATE ALL EXPERIMENTAL SCENARIOS
# ============================================================

# This base-R function is also checked independently by check_demo_figures.R.
make_adi_scenarios <- function(p_common, p_female, p_male, a, d, i) {
  scenarios <- list()
  for (sex in c("Female", "Male")) {
    for (mating in c("Common reference and mate population", "Differentiated mate pools")) {
      pf <- if (mating == "Common reference and mate population") p_common else p_female
      pm <- if (mating == "Common reference and mate population") p_common else p_male
      p_ref <- if (sex == "Female") pf else pm
      for (origin in c("mat", "pat", "none")) {
        method <- if (origin == "none") "Between-generation (b)" else "Within-generation (w)"
        tmp <- if (origin == "none") {
          make_between(sex, pm, pf, a, d, i)
        } else {
          make_within(origin, p_ref, a, d, i)
        }
        tmp$Sex <- sex
        tmp$Method <- method
        tmp$Mating <- mating
        tmp$Origin <- origin
        tmp$ReferenceP <- p_ref
        scenarios[[paste(sex, method, mating, origin, sep = "_")]] <- tmp
      }
    }
  }
  do.call(rbind, scenarios)
}

def_lm <- make_adi_scenarios(pPop, pFemales, pMales, a, d, i)



# This base-R builder is also evaluated by check_demo_figures.R.
make_adi_display <- function(scenarios) {
  bubbles <- fitted <- segments <- curves <- squares <- annotations <- list()
  for (pool in unique(scenarios$Mating)) for (sex in c("Female", "Male")) {
    key <- paste(pool, sex)
    tag <- function(z, method) {
      z$Mating <- pool; z$Sex <- sex; z$Method <- method
      z
    }
    z <- subset(scenarios, Mating == pool & Sex == sex & Origin == "mat")
    z$mat <- as.numeric(substr(z$ValueGenotype, 1, 1))
    z$pat <- as.numeric(substr(z$ValueGenotype, 3, 3))
    z$Genotype <- z$mat + z$pat
    fit <- lm(Value ~ mat + pat, data=z, weights=probability)
    prediction <- as.numeric(predict(fit))
    mu <- weighted.mean(z$Value, z$probability)
    bubbles[[paste0(key,"w")]] <- tag(z[c("Genotype","Value","ValueGenotype","probability")], "Within-generation")
    fitted[[key]] <- tag(data.frame(Genotype=z$Genotype, ValueGenotype=z$ValueGenotype,
      fitted_value=prediction, reference_mean=mu, A=prediction-mu), "Within-generation")
    edge <- data.frame(from=c("0/0","0/1","0/0","1/0"),
                       to=c("1/0","1/1","0/1","1/1"), Origin=c("mat","mat","pat","pat"))
    from <- match(edge$from,z$ValueGenotype); to <- match(edge$to,z$ValueGenotype)
    edge$Genotype <- z$Genotype[from]; edge$GenotypeEnd <- z$Genotype[to]
    edge$Value <- prediction[from]; edge$ValueEnd <- prediction[to]
    edge$alpha <- as.numeric(coef(fit)[edge$Origin])
    stopifnot(all(abs(edge$ValueEnd-edge$Value-edge$alpha)<1e-10))
    segments[[key]] <- tag(edge,"Within-generation")
    annotations[[paste0(key,"w")]] <- tag(data.frame(Genotype=1,
      Value=c(-1.4,-2.05), regression=as.numeric(coef(fit)[c("mat","pat")]),
      alpha_label=sprintf("alpha[plain(%s)]^w == '%.2f'",c("mat","pat"),coef(fit)[c("mat","pat")])),
      "Within-generation")

    z <- subset(scenarios,Mating==pool & Sex==sex & Origin=="none")
    z <- aggregate(probability ~ Genotype+Value+ValueGenotype,data=z,FUN=sum)
    fit <- lm(Value ~ Genotype,data=z,weights=probability)
    mu <- weighted.mean(z$Value,z$probability)
    bubbles[[paste0(key,"b")]] <- tag(z,"Between-generation")
    x <- seq(0,2,length.out=101)
    expected <- as.numeric(predict(fit,data.frame(Genotype=x)))
    curves[[key]] <- tag(data.frame(Genotype=x, expected=expected,
      breeding=2*expected-mu, progeny_mean=mu),"Between-generation")
    expected <- as.numeric(predict(fit,data.frame(Genotype=0:2)))
    squares[[key]] <- tag(data.frame(Genotype=0:2, fitted_value=2*expected-mu,
      expected=expected, progeny_mean=mu, B=2*(expected-mu)),"Between-generation")
    role <- if (sex=="Female") "f" else "m"
    annotations[[paste0(key,"b")]] <- tag(data.frame(
      Genotype=if (pool=="Differentiated mate pools" && sex=="Male") 2.3 else 1,
      Value=-1.65,
      regression=2*unname(coef(fit)[2]),
      alpha_label=sprintf("paste(alpha[%s]^b, ' = ', 2*b[%s], ' = %.2f')",role,role,2*coef(fit)[2])),
      "Between-generation")
  }
  bind <- function(parts) {
    z <- do.call(rbind,parts)
    z$Mating <- factor(z$Mating,c("Common reference and mate population","Differentiated mate pools"),
                      labels=c("Undifferentiated mate pools","Differentiated mate pools"))
    z$Sex <- factor(z$Sex,c("Female","Male"))
    z$Method <- factor(z$Method,c("Within-generation","Between-generation"))
    z
  }
  lapply(list(bubbles=bubbles,fitted=fitted,segments=segments,curves=curves,
              squares=squares,annotations=annotations),bind)
}
adi_display <- make_adi_display(def_lm)

p <- ggplot() +
  geom_segment(data=adi_display$segments,
    aes(x=Genotype,xend=GenotypeEnd,y=Value,yend=ValueEnd,linetype=Origin),
    colour="black",linewidth=.55) +
  geom_point(data=adi_display$bubbles,
    aes(Genotype,Value,fill=ValueGenotype,size=probability),
    shape=21,colour="grey20",stroke=.25) +
  geom_line(data=adi_display$curves,aes(Genotype,breeding,linetype="B"),colour="black",linewidth=.65) +
  geom_point(data=adi_display$fitted,
    aes(Genotype,fitted_value,shape="A"),size=2.2,fill="white",colour="black",stroke=.55) +
  geom_point(data=adi_display$squares,
    aes(Genotype,fitted_value,shape="B"),size=2.2,fill="white",colour="black",stroke=.55) +
  geom_text(data=subset(adi_display$annotations,
      !(Method=="Between-generation" & Mating=="Differentiated mate pools" & Sex=="Male")),
    aes(Genotype,Value,label=alpha_label),parse=TRUE,size=3) +
  geom_label(data=subset(adi_display$annotations,
      Method=="Between-generation" & Mating=="Differentiated mate pools" & Sex=="Male"),
    aes(Genotype,Value,label=alpha_label),parse=TRUE,hjust=1,vjust=.5,size=3,
    fill=scales::alpha("white",.88),colour="black",linewidth=0,
    label.padding=grid::unit(.08,"lines"),label.r=grid::unit(0,"lines"),
    show.legend=FALSE) +
  facet_nested(Method ~ Mating+Sex,axes="x") +
  alpha_fill_scale(imprinting=TRUE) +
  scale_size_area(max_size=8,limits=c(0,.9),breaks=c(.1,.3,.5,.7,.9),
    labels=scales::label_percent(accuracy=1),name="Frequency",guide=alpha_size_guide()) +
  alpha_model_scales(imprinting=TRUE) +
  scale_x_continuous(breaks=0:2,limits=c(-.4,2.4),expand=c(0,0)) +
  scale_y_continuous(limits=c(-2.45,2.05),breaks=-2:2,expand=c(0,0)) +
  labs(x=expression("Allele dosage, "*N[1]),y="Genetic value") +
  theme_classic(base_size=11,base_family="sans") +
  theme(panel.grid.major=element_line(colour="#E5E5E5",linewidth=.3),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black",linewidth=.4),
    axis.ticks=element_line(colour="black",linewidth=.4),axis.ticks.length=grid::unit(2,"mm"),
    axis.text=element_text(colour="black",size=9),axis.title=element_text(size=10),
    strip.background=element_blank(),strip.text=element_text(face="bold",size=9.5),
    panel.spacing.x=grid::unit(3,"mm"),panel.spacing.y=grid::unit(4,"mm"),
    legend.title=element_text(size=9,face="bold"),legend.text=element_text(size=8.5),
    legend.key.width=grid::unit(8,"mm"),legend.spacing.y=grid::unit(0,"mm"),
    plot.margin=margin(3,8,3,5,unit="mm"))
if (interactive()) print(p)
if (genotype_colours) {
  source(file.path(figure_dir, "Figure3_Variants.R"))
} else {
  ggsave(filename=file.path(figure_dir,"alpha_ADI_def.pdf"),
         plot=p,width=180,height=alpha_figure_height,units="mm",device=grDevices::pdf,
         useDingbats=FALSE,bg="white")
  render_alpha_png("alpha_ADI_def")
}

# These options regenerate Figures 2-3 only; stop before the Figure 4 heat map.
if (genotype_colours || "--alpha-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no")

# Figure 4
library(dplyr)
library(scico)
source(file.path(figure_dir, "figure_style.R"))
a <- 1
d_values <- seq(0, 1.5, 0.5)
i_values <- seq(-0.5, 0.5, 0.5)
p_values <- seq(0, 1, 0.05)

delta.df <- expand.grid(p_female=p_values, p_male=p_values,
                        d=d_values, i=i_values) %>%
  mutate(q_female=1-p_female, q_male=1-p_male,
    alpha_female = a + d * (q_male - p_male) - i,
    alpha_male   = a + d * (q_female - p_female) + i,
    Delta = alpha_male - alpha_female
  )

heat_factors <- function(z) {
  z$d <- factor(z$d, levels=d_values, labels=paste0("d == ",d_values))
  z$i <- factor(z$i, levels=i_values, labels=paste0("i == ",i_values))
  z
}
zero.df <- heat_factors(subset(delta.df, d > abs(i)))
# Signs and the zero contour make the direction readable without colour.
signs <- do.call(rbind,lapply(seq_len(length(d_values)*length(i_values)), function(k) {
  z <- expand.grid(d=d_values,i=i_values)[k,]
  lower <- -2*z$d+2*z$i; upper <- 2*z$d+2*z$i
  if (lower < 0 && upper > 0) {
    negative <- c(.10,.84)
    positive <- c(.90,.16)
    data.frame(d=z$d,i=z$i,p_male=c(negative[1],positive[1]),
               p_female=c(negative[2],positive[2]),sign=c("-","+"))
  } else {
    data.frame(z,p_male=.5,p_female=.5,
               sign=if (lower==0 && upper==0) "0" else if (upper<=0) "-" else "+")
  }
}))
signs <- heat_factors(signs)
delta.df <- heat_factors(delta.df)

p <- ggplot(delta.df,aes(p_male,p_female,fill=Delta)) +
  geom_tile(width=.05,height=.05) +
  geom_contour(data=zero.df,aes(p_male,p_female,z=Delta),inherit.aes=FALSE,
               breaks=0,colour="white",linewidth=.65) +
  geom_contour(data=zero.df,aes(p_male,p_female,z=Delta),inherit.aes=FALSE,
               breaks=0,colour="#404040",linewidth=.25) +
  geom_label(data=signs,aes(p_male,p_female,label=sign),inherit.aes=FALSE,
             size=9/(72.27/25.4),fill=scales::alpha("white",.75),linewidth=0,
             label.padding=grid::unit(.03,"lines")) +
  facet_grid(i ~ d,labeller=label_parsed) +
  scale_fill_scico(palette="vik",limits=c(-4,4),midpoint=0,
    breaks=c(-4,-2,0,2,4),name=expression(Delta * alpha^b),
    guide=guide_colorbar(title.position="top",title.hjust=.5,
                        barwidth=grid::unit(48,"mm"),barheight=grid::unit(3,"mm"),
                        frame.colour="black",ticks.colour="black")) +
  scale_x_continuous(breaks=c(0,1),expand=c(0,0)) +
  scale_y_continuous(breaks=c(0,1),expand=c(0,0)) +
  coord_fixed(xlim=c(0,1),ylim=c(0,1),expand=FALSE) +
  labs(x=expression("Male allele frequency, "*p[m]),
       y=expression("Female allele frequency, "*p[f])) +
  paper_theme() +
  theme(panel.grid=element_blank(),panel.spacing=grid::unit(3.5,"mm"),
        strip.text=element_text(size=8.5),legend.position="bottom",
        legend.box.spacing=grid::unit(1,"mm"),legend.margin=margin(0,0,0,0),
        axis.title=element_text(size=9),plot.margin=margin(2,2,2,2,unit="mm"))
if (interactive()) print(p)
save_paper_plot(p,"delta_ADI",paper_column_mm,94)

#Figure3_alt
# ============================================================
# 1. LIBRARIES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(viridis)
library(scales)
library(ragg)

# ============================================================
# 2. PARAMETERS & DISTRIBUTIONS
# ============================================================

pPop <- 0.6
pMales <- 0.9
pFemales <- 0.3

a <- 1
d <- 0.75
i <- 0.5

qPop <- 1 - pPop
qMales <- 1 - pMales
qFemales <- 1 - pFemales

make_within <- function(origin, p, a, d, i) {
  origin <- match.arg(origin, c("mat", "pat"))
  mat <- c(0, 0, 1, 1)
  pat <- c(0, 1, 0, 1)
  q <- 1 - p
  data.frame(
    Genotype = if (origin == "mat") mat else pat,
    ValueGenotype = c("0/0", "0/1", "1/0", "1/1"),
    Value = c(-a, d + i, d - i, a),
    probability = c(q^2, q * p, p * q, p^2)
  )
}

make_between <- function(sex, p_male, p_female, a, d, i) {
  qm <- 1 - p_male
  qf <- 1 - p_female
  
  if (sex == "Male") {
    data.frame(
      Genotype = c(2, 2, 1, 1, 1, 1, 0, 0),
      ValueGenotype = c("1/1", "0/1", "1/1", "0/1", "1/0", "0/0", "0/0", "1/0"),
      Value = c(a, d + i, a, d + i, d - i, -a, -a, d - i),
      probability = c(p_male^2 * p_female, p_male^2 * qf, p_male * qm * p_female, 
                      p_male * qm * qf, p_male * qm * p_female, p_male * qm * qf, 
                      qm^2 * qf, qm^2 * p_female)
    )
  } else {
    data.frame(
      Genotype = c(2, 2, 1, 1, 1, 1, 0, 0),
      ValueGenotype = c("1/1", "1/0", "1/1", "1/0", "0/1", "0/0", "0/0", "0/1"),
      Value = c(a, d - i, a, d - i, d + i, -a, -a, d + i),
      probability = c(p_female^2 * p_male, p_female^2 * qm, p_female * qf * p_male, 
                      p_female * qf * qm, p_female * qf * p_male, p_female * qf * qm, 
                      qf^2 * qm, qf^2 * p_male)
    )
  }
}

make_adi_scenarios <- function(p_common, p_female, p_male, a, d, i) {
  scenarios <- list()
  for (sex in c("Female", "Male")) {
    for (mating in c("Common reference and mate population", "Differentiated mate pools")) {
      pf <- if (mating == "Common reference and mate population") p_common else p_female
      pm <- if (mating == "Common reference and mate population") p_common else p_male
      p_ref <- if (sex == "Female") pf else pm
      for (origin in c("mat", "pat", "none")) {
        method <- if (origin == "none") "Between-generation (b)" else "Within-generation (w)"
        tmp <- if (origin == "none") {
          make_between(sex, pm, pf, a, d, i)
        } else {
          make_within(origin, p_ref, a, d, i)
        }
        tmp$Sex <- sex; tmp$Method <- method; tmp$Mating <- mating; tmp$Origin <- origin
        scenarios[[paste(sex, method, mating, origin, sep = "_")]] <- tmp
      }
    }
  }
  do.call(rbind, scenarios)
}

def_lm <- make_adi_scenarios(pPop, pFemales, pMales, a, d, i)

# ============================================================
# 3. BUILDER SIMPLIFICADO: Jerarquía para ggh4x
# ============================================================

make_adi_display <- function(scenarios) {
  bubbles <- lines <- annotations <- list()
  
  for (pool in unique(scenarios$Mating)) for (sex in c("Female", "Male")) {
    key <- paste(pool, sex)
    
    # Adaptado para aceptar un sub-método
    tag <- function(z, method, sub_method) {
      z$Mating <- pool; z$Sex <- sex; 
      z$Method <- method; z$SubMethod <- sub_method
      z
    }
    
    # -------------------------------------------------------------
    # WITHIN-GENERATION
    # -------------------------------------------------------------
    z_base <- subset(scenarios, Mating == pool & Sex == sex & Origin == "mat")
    z_base$mat <- as.numeric(substr(z_base$ValueGenotype, 1, 1))
    z_base$pat <- as.numeric(substr(z_base$ValueGenotype, 3, 3))
    
    fit <- lm(Value ~ mat + pat, data=z_base, weights=probability)
    mu <- weighted.mean(z_base$Value, z_base$probability)
    
    # -- Maternal --
    z_mat <- z_base
    z_mat$Genotype <- z_mat$mat
    bubbles[[paste0(key,"w_mat")]] <- tag(z_mat[c("Genotype","Value","ValueGenotype","probability")], 
                                          "Within-generation", "Maternal")
    
    mean_X_mat <- weighted.mean(z_mat$Genotype, z_mat$probability)
    int_mat <- mu - coef(fit)["mat"] * mean_X_mat
    lines[[paste0(key,"w_mat")]] <- tag(data.frame(Genotype=c(0, 1), Value=int_mat + coef(fit)["mat"]*c(0, 1)), 
                                        "Within-generation", "Maternal")
    
    annotations[[paste0(key,"w_mat")]] <- tag(data.frame(Genotype=0.5, Value=-2.05, alpha_label=sprintf("alpha[plain(mat)]^w == '%.2f'",coef(fit)["mat"])), 
                                              "Within-generation", "Maternal")
    
    # -- Paternal --
    z_pat <- z_base
    z_pat$Genotype <- z_pat$pat
    bubbles[[paste0(key,"w_pat")]] <- tag(z_pat[c("Genotype","Value","ValueGenotype","probability")], 
                                          "Within-generation", "Paternal")
    
    mean_X_pat <- weighted.mean(z_pat$Genotype, z_pat$probability)
    int_pat <- mu - coef(fit)["pat"] * mean_X_pat
    lines[[paste0(key,"w_pat")]] <- tag(data.frame(Genotype=c(0, 1), Value=int_pat + coef(fit)["pat"]*c(0, 1)), 
                                        "Within-generation", "Paternal")
    
    annotations[[paste0(key,"w_pat")]] <- tag(data.frame(Genotype=0.5, Value=-2.05, alpha_label=sprintf("alpha[plain(pat)]^w == '%.2f'",coef(fit)["pat"])), 
                                              "Within-generation", "Paternal")
    
    # -------------------------------------------------------------
    # BETWEEN-GENERATION
    # -------------------------------------------------------------
    z <- subset(scenarios, Mating==pool & Sex==sex & Origin=="none")
    z <- aggregate(probability ~ Genotype+Value+ValueGenotype, data=z, FUN=sum)
    fit_b <- lm(Value ~ Genotype, data=z, weights=probability)
    
    bubbles[[paste0(key,"b")]] <- tag(z, "Between-generation", " ") # Espacio en blanco para no separar
    lines[[paste0(key,"b")]] <- tag(data.frame(Genotype=c(0, 2), Value=predict(fit_b, data.frame(Genotype=c(0, 2)))), 
                                    "Between-generation", " ")
    
    role <- if (sex=="Female") "f" else "m"
    annotations[[paste0(key,"b")]] <- tag(data.frame(
      Genotype=if (pool=="Differentiated mate pools" && sex=="Male") 2.3 else 1,
      Value=-1.65,
      alpha_label=sprintf("paste(alpha[%s]^b, ' = ', 2*b[%s], ' = %.2f')",role,role,2*coef(fit_b)[2])),
      "Between-generation", " ")
  }
  
  bind <- function(parts) {
    z <- do.call(rbind, parts)
    z$Mating <- factor(z$Mating, c("Common reference and mate population","Differentiated mate pools"),
                       labels=c("Undifferentiated mate pools","Differentiated mate pools"))
    z$Sex <- factor(z$Sex, c("Female","Male"))
    z$Method <- factor(z$Method, c("Within-generation", "Between-generation"))
    z$SubMethod <- factor(z$SubMethod, c("maternal", "paternal", " ")) # Orden de la fila anidada
    z
  }
  
  lapply(list(bubbles=bubbles, lines=lines, annotations=annotations), bind)
}

adi_display <- make_adi_display(def_lm)

# ============================================================
# 4. PLOT FINAL (Limpio)
# ============================================================

p <- ggplot() +
  # 1. Recta de regresión
  geom_line(data=adi_display$lines, aes(x=Genotype, y=Value), colour="black", linewidth=0.65) +
  # 2. Burbujas
  geom_point(data=adi_display$bubbles,
             aes(Genotype, Value, fill=ValueGenotype, size=probability),
             shape=21, colour="grey20", stroke=0.25) +
  # 3. Etiquetas
  geom_text(data=subset(adi_display$annotations,
                        !(Method=="Between-generation" & Mating=="Differentiated mate pools" & Sex=="Male")),
            aes(Genotype, Value, label=alpha_label), parse=TRUE, size=3) +
  geom_label(data=subset(adi_display$annotations,
                         Method=="Between-generation" & Mating=="Differentiated mate pools" & Sex=="Male"),
             aes(Genotype, Value, label=alpha_label), parse=TRUE, hjust=1, vjust=0.5, size=3,
             fill=scales::alpha("white", 0.88), colour="black", linewidth=0,
             label.padding=grid::unit(0.08,"lines"), label.r=grid::unit(0,"lines"),
             show.legend=FALSE) +
  facet_nested(Method + SubMethod ~ Mating + Sex, axes="x") +
  # alpha_fill_scale(imprinting=TRUE) + # (Usa tus funciones personalizadas de colores aquí)
  scale_size_area(max_size=8, limits=c(0, 0.9), breaks=c(0.1, 0.3, 0.5, 0.7, 0.9),
                  labels=scales::label_percent(accuracy=1), name="Frequency") +
  # alpha_model_scales(imprinting=TRUE) + # (Usa tus escalas personalizadas aquí)
  scale_x_continuous(breaks=0:2, limits=c(-0.4, 2.4), expand=c(0,0)) +
  scale_y_continuous(limits=c(-2.45, 2.05), breaks=-2:2, expand=c(0,0)) +
  labs(x="Allele dosage (gametic for Within, total for Between)", y="Genetic value") +
  theme_classic(base_size=11, base_family="sans") +
  theme(panel.grid.major=element_line(colour="#E5E5E5", linewidth=0.3),
        panel.grid.minor=element_blank(), axis.line=element_line(colour="black", linewidth=0.4),
        axis.ticks=element_line(colour="black", linewidth=0.4), axis.ticks.length=grid::unit(2,"mm"),
        axis.text=element_text(colour="black", size=9), axis.title=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(face="bold", size=9.5),
        panel.spacing.x=grid::unit(3,"mm"), panel.spacing.y=grid::unit(4,"mm"),
        legend.title=element_text(size=9, face="bold"), legend.text=element_text(size=8.5),
        legend.key.width=grid::unit(8,"mm"), legend.spacing.y=grid::unit(0,"mm"),
        plot.margin=margin(3, 8, 3, 5, unit="mm"))
if (interactive()) print(p)
save_paper_plot(p,"delta_ADI_alt",paper_column_mm,94)
