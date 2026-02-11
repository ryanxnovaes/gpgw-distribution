# ============================================================
# Project: Gamma Power Generalized Weibull Distribution
# Created by: Arioane P. Soares & Fernando A. Peña-Ramírez
# Edited, improved, and revised by: Ryan Novaes Pereira
# Date: February 2026
# ============================================================

# ============================================================
# 00. Load Required Libraries
# ============================================================

library(AdequacyModel)
library(readr)

# ============================================================
# 01. Load and Prepare the Dataset
# ============================================================

source("codes/GPGWfunc.R")

scimago <- read_delim("data/database.csv", delim = ";",
                      escape_double = FALSE, trim_ws = TRUE)

x <- as.numeric(gsub(",", ".", trimws(scimago$SJR)))
x <- x[is.finite(x)]

# ============================================================
# 02. Initial values
# ============================================================

init_from_data <- function(x, k) {
  as.numeric(quantile(
    x,
    probs = seq(1/(k+1), k/(k+1), length.out = k)
  ))
}

init <- list(
  p1 = mean(x),
  p2 = init_from_data(x,2),
  p3 = init_from_data(x,3),
  p4 = init_from_data(x,4)
)

# ============================================================
# 03. Storage objects
# ============================================================

results <- matrix(NA, nrow = 14, ncol = 7)
colnames(results) <- c("AIC","CAIC","BIC","HQIC","W*","A*","KS")
rownames(results) <- c(
  "EGPW","Kw-W","GPW","EW","EXP","ENH","NH",
  "GEXP","GRay","Weibull","GWeibull","Rayleigh","GPGW","GNH"
)

estim <- matrix(NA, nrow = 28, ncol = 4)
colnames(estim) <- c("alpha","lambda","nu","a")
rownames(estim) <- c(
  "EGPW","", "Kw-W","", "GPW","", "EW","", "EXP","",
  "ENH","", "NH","", "GEXP","", "GRay","",
  "Weibull","", "GWeibull","", "Rayleigh","",
  "GPGW","", "GNH",""
)

# ============================================================
# 04. Model fitting
# ============================================================

nc <- 0
mods <- vector("list", 14)

fit_and_store <- function(i, pdf, cdf, starts, domain, npar) {
  
  mod <- try(
    goodness.fit(
      pdf    = pdf,
      cdf    = cdf,
      starts = starts,
      data   = x,
      method = "CG",
      domain = domain,
      mle    = NULL
    ),
    silent = TRUE
  )
  
  if (inherits(mod, "try-error")) return(NULL)
  
  results[i, ] <<- c(
    mod$AIC, mod$CAIC, mod$BIC,
    mod$HQIC, mod$W, mod$A,
    mod$KS$statistic
  )
  
  r <- 2 * i - 1
  
  tryCatch({
    
    k <- length(mod$mle)  
    
    estim[r,     1:k] <<- format(mod$mle,  digits = 4)
    estim[r + 1, 1:k] <<- format(mod$Erro, digits = 4)
    
  }, error = function(e) {
    warning(
      paste("Estimation storage failed for model", rownames(results)[i])
    )
  })
  
  return(mod)
}

models <- list(
  list(f_EGPW,       F_EGPW,       init$p4, c(0,1e15), 4),
  list(f_kww,        F_kww,        init$p4, c(0,1e15), 4),
  list(f_GPW,        F_GPW,        init$p3, c(0,1e15), 3),
  list(f_expweibull, F_expweibull, init$p3, c(0,1e15), 3),
  list(fdp_EXP,      Fda_EXP,      init$p1, c(0,1e15), 1),
  list(fdp_enh,      fda_enh,      init$p3, c(0,Inf),  2),
  list(fdp_nh,       fda_nh,       init$p2, c(0,Inf),  2),
  list(fdp_GEXP,     Fda_GEXP,     init$p2, c(0,1e15), 2),
  list(fdp_GR,       Fda_GR,       init$p2, c(0,1e15), 2),
  list(fdp_we,       fda_we,       init$p2, c(0,Inf),  2),
  list(fdp_GW,       Fda_GW,       init$p3, c(0,1e15), 3),
  list(f_Ray,        F_Ray,        init$p1, c(0,Inf),  1),
  list(fdp_gam,      Fda_gam,      init$p4, c(0,1e15), 4),
  list(fdp_gam1,     Fda_gam1,     init$p3, c(0,1e15), 3)
)


mods <- vector("list", length(models))

for (i in seq_along(models)) {
  mods[[i]] <- fit_and_store(
    i      = i,
    pdf    = models[[i]][[1]],
    cdf    = models[[i]][[2]],
    starts = models[[i]][[3]],
    domain = models[[i]][[4]],
    npar   = models[[i]][[5]]
  )
}

print(results, quote = FALSE)

# ============================================================
# Final Figure: Histogram + densities (top) | Boxplot (bottom)
# ============================================================

x <- x[is.finite(x)]
x_sub <- x[x < 4]

# pdf(
#   file      = "hist_App.pdf",
#   width     = 6,
#   height    = 5.6,
#   family    = "Times",
#   pointsize = 12
# )

layout(
  matrix(c(1, 2), nrow = 2),
  heights = c(3, 1.55)
)

# ------------------------------------------------------------
# Panel 1: Histogram + densities
# ------------------------------------------------------------

par(mar = c(2.2, 4, 1, 1))  

hist(
  x_sub,
  probability = TRUE,
  breaks      = "FD",
  col         = "gray80",
  border      = "gray80",
  xlab        = "",
  ylab        = "Density",
  main        = ""
)

models_plot <- list(
  list(name="GPGW", fun=fdp_gam,      mod=mods[[13]], col=1, lty=1),
  list(name="EGPW", fun=f_EGPW,       mod=mods[[1]],  col=2, lty=2),
  list(name="KWW",  fun=f_kww,        mod=mods[[2]],  col=3, lty=3),
  list(name="EW",   fun=f_expweibull, mod=mods[[4]],  col=4, lty=4)
)

for (m in models_plot) {
  if (!is.null(m$mod)) {
    curve(
      m$fun(m$mod$mle, x),
      from = min(x_sub),
      to   = max(x_sub),
      add  = TRUE,
      lwd  = 2.5,
      col  = m$col,
      lty  = m$lty
    )
  }
}

legend(
  "topright",
  legend = sapply(models_plot, `[[`, "name"),
  col    = sapply(models_plot, `[[`, "col"),
  lty    = sapply(models_plot, `[[`, "lty"),
  lwd    = 2.5,
  bty    = "n",
  cex    = 1.05
)

# ------------------------------------------------------------
# Panel 2: Boxplot (closer to histogram)
# ------------------------------------------------------------

par(mar = c(3.8, 4, 0.1, 1))  

boxplot(
  x_sub,
  horizontal = TRUE,
  col        = "gray80",
  outline    = TRUE,
  frame      = FALSE,
  axes       = FALSE,
  pch        = 1,
  cex        = 0.6
)

axis(1)
mtext("SJR", side = 1, line = 2.4)

# dev.off()
