## specialists_analyses
# follows script "data_simulation.R"
library(phytools)
library(geiger)
library(caper)
source("R/function_traitplot2.R")

## data_preparation
inv_logit <- function(x) exp(x)/(1 + exp(x))
pm_tree$node.label <- NULL
phy_fe <- keep.tip(pm_tree, fe$splist$pm_names)
phy_oc <- keep.tip(pm_tree, oc$splist$pm_names)
phy_mt <- keep.tip(pm_tree, mt$splist$pm_names)
calc_phydist <- function(x, dat, phy, splist, R = 9999, pool = "spec"){
  phylodist_funs <- function(sel_spec, phy){
    phy <- keep.tip(phy, sel_spec)
    dmat <- cophenetic(phy)
    dmat_lt <- dmat[lower.tri(dmat)]
    c(phydivers = sum(phy$edge.length), phydist = sum(dmat_lt),  
      meannn = mean(apply(dmat, 1, function(x) min(x[x != 0])))) #mean nearest neighbour
  }
  if (pool == "spec") sppool <- splist$pm_names[splist$Specialist > 0.5]
  if (pool == "gene") sppool <- splist$pm_names[splist$Specialist < 0.5]
  if (pool == "all") sppool <- splist$pm_names
  phy_all <- keep.tip(phy, sppool)
  sel_spec <- sppool[sppool %in% colnames(dat)[x > 0.5]]
  if (pool == "all") sel_spec <- splist$pm_names[splist$Specialist > 0.5]
  selections <- sapply(1:R, function(x) sample(phy_all$tip.label, length(sel_spec)))  
  null_mod <- apply(selections, 2, phylodist_funs, phy_all)
  obs <- phylodist_funs(sel_spec, phy_all)
  c(SES = (obs - rowMeans(null_mod)) / apply(null_mod, 1, sd), 
    pval = rowSums(null_mod < obs) / (R + 1))
}

## analyses
# H1
phydist_full <- function(dat, phy){
  calc_phydist(1, dat$pres[, -1], phy, dat$splist, pool = "all")
}

d_fe <- phylo.d(fe$splist, phy_fe, pm_names, Specialist, 10000)
d_oc <- phylo.d(oc$splist, phy_oc, pm_names, Specialist, 10000)
d_mt <- phylo.d(mt$splist, phy_mt, pm_names, Specialist, 10000)

tabS1 <- round(rbind(fe = phydist_full(fe, phy_fe),
  oc = phydist_full(oc, phy_oc),
  mt = phydist_full(mt, phy_mt)), 3)

fit_mkdelta <- function(splist, phy){
  aux_spec <- setNames(splist$Specialist, splist$pm_names)
  set.seed(10)
  mod_mkd <- fitDiscrete(multi2di(phy), factor(aux_spec), model = "ARD", 
    transform = "delta", bounds = list(delta = c(0,20)), 
    control = list(niter = 200))
  mod_mk <- fitDiscrete(multi2di(phy), factor(aux_spec), model = "ARD",  
    control = list(niter = 200))
  list(mod_mkd = mod_mkd, mod_mk = mod_mk, dAIC = AIC(mod_mk) - AIC(mod_mkd))
}

mkdelta_fe <- fit_mkdelta(fe$splist, phy_fe)
# mkdelta_oc <- fit_mkdelta(oc$splist, phy_oc)
# mkdelta_mt <- fit_mkdelta(mt$splist, phy_mt)

# Figure 1
# requires function "trait.plot2"
cols <- data.frame(cols = c("white","red"), row.names = 0:1)
aux_plot <- function(phy, dat, cols, leg = FALSE, title = ""){
  par(mai = c(0,0,0,0))
  t_spe <- data.frame(Specialist = dat$splist$Specialist, row.names = dat$splist$pm_names)
  cla <- dat$splist$family[match(phy$tip.label, dat$splist$pm_names)]
  cla[!cla %in% names(sort(table(cla), decreasing = TRUE)[1:15])] <- NA
  trait.plot2(phy, t_spe, cols, legend = FALSE, 
    class = cla, margin = 0.45, cex.lab = 0.7, g_lwd = 3, w = 1/30)
  if (leg) legend("bottomleft", fill = cols$cols, 
    legend = c("Nonspecialist","Specialist"), bty = "n")
  legend("topleft", bty = "n", legend = "", title = title, cex = 2.2, title.col = cols[2,])
}

cols2 <- c("blue", "green", "black")
cols2 <- c(4,3,1)
# png("figures/Fig1_spec.png", width = 480*10, height = 480*10*3, res = 72*15)
par(mfrow = c(3,1))
cols <- data.frame(cols = c("white",cols2[1]), row.names = 0:1)
aux_plot(phy_fe, fe, cols, title = "Fens")
cols <- data.frame(cols = c("white",cols2[2]), row.names = 0:1)
aux_plot(phy_oc, oc, cols, title = "Outcrops")
cols <- data.frame(cols = c("white",cols2[3]), row.names = 0:1)
aux_plot(phy_mt, mt, cols, title = "Mountaintops")
#legend("bottomleft", fill = cols2, 
#  legend = c("Fen specialist", "Outcrop specialist", "Mountaintop specialist"), bty = "n")

# H2
aux_calc_phydist <- function(dat, phy, pool = "spec"){
  data.frame(t(apply(dat$pres[, -1], 1, calc_phydist, dat$pres[, -1], 
    phy, dat$splist, pool = pool)))
}
pd_fe <- aux_calc_phydist(fe, phy_fe)
pd_oc <- aux_calc_phydist(oc, phy_oc)
pd_mt <- aux_calc_phydist(mt, phy_mt)
pdg_fe <- aux_calc_phydist(fe, phy_fe, "gene")
pdg_oc <- aux_calc_phydist(oc, phy_oc, "gene")
pdg_mt <- aux_calc_phydist(mt, phy_mt, "gene")


phyd_ins <- function(phyd, ins, col, xt = NULL, yt){
  ins_aux <- ins
  mod <- lm(phyd ~ ins_aux)
  new_x <- seq(min(ins), max(ins), length.out = 100)
  pred <- predict(mod, list(ins_aux = new_x), interval = "confidence")
  plot(range(c(fe$ins$TE, oc$ins$TE, mt$ins$TE)), range(phyd), type = "n", axes = FALSE, 
    xlab = "",
    ylab = expression('PD'[SES]))
  abline(h = c(-1,1) * 1.96, col = "grey", lwd = 2, lty = 2)
  points(phyd ~ ins, col = col)
  lines(new_x, pred[, 1], col = col, lwd = 2)
  lines(new_x, pred[, 2], col = col, lty = 2)
  lines(new_x, pred[, 3], col = col, lty = 2)
  box(bty = "l")
  axis(2, las = 2, labels = -3:3, at = -3:3)
  aux <- round(c(summary(mod)$coefficients[2, 4], summary(mod)$adj.r.squared), 2)
  if (!is.null(xt)) text(xt, yt, col = col,
    labels = bquote(atop(P:~.(aux[1]),R^2:~.(aux[2]))), adj = c(0,1))
  invisible(summary(mod))
}

# Table S2
extr_fun <- function(x, ins){
  mod <- summary(lm(x ~ ins))
  setNames(c(mod$coefficients[c(2,8)], mod$adj.r.squared), c("Slope", "P", "Adj-Rsq"))
}
round(rbind(apply(pd_fe, 2, extr_fun, fe$ins$TE),
  apply(pd_oc, 2, extr_fun, oc$ins$TE),
  apply(pd_mt, 2, extr_fun, mt$ins$TE)), 3)[, 1:3]

# Figure 2
cols <- c("blue", "green", "black")
cols <- c(4,3,1)
# png("figures/Fig2_phydivers.png", height = 480*10, width = 480*10, res = 72*12.6)
par(mai = c(0.1,0.55,0.1,0.1), mfrow = c(4,1))
phyd_ins(pd_fe$SES.phydivers, fe$ins$TE, col = cols[1], -1.4, par("usr")[4])
text(0.3, -1, "Fens", col = cols[1], cex = 2.2)
phyd_ins(pd_oc$SES.phydivers, oc$ins$TE, col = cols[2], -1.4, par("usr")[4])
text(0.3, -1, "Outcrops", col = cols[2], cex = 2.2)
phyd_ins(pd_mt$SES.phydivers, mt$ins$TE, col = cols[3], 6, par("usr")[4])
text(7, 0, "Mountaintops", col = cols[3], cex = 2.2)
plot(range(c(fe$ins$TE, oc$ins$TE, mt$ins$TE)), 0:1, type = "n", ann = FALSE, axes = FALSE)
axis(1, line = -10, lwd = 0, lwd.tick = 1)
text(mean(par("usr")[1:2]), 0.85, "Insularity (target effect)")
#par(new = TRUE, mfrow = c(1,1))
#plot(range(c(fe$ins$TE, oc$ins$TE, mt$ins$TE)), 0:1, type = "n", ann = FALSE, axes = FALSE)
#legend(mean(par("usr")[1:2])-1, 0.1, pch = 1, lwd = 1, col = cols, 
#  legend = c("Fens", "Outcrops", "Mountaintops"), bty = "n", cex = 0.7)

# Figure S1
# png("figures/FigS1_phydivers.png", height = 480*10, width = 480*10, res = 72*12.6)
par(mai = c(0.1,0.55,0.1,0.1), mfrow = c(4,1))
phyd_ins(pdg_fe$SES.phydivers, fe$ins$TE, col = cols[1], -1.4, par("usr")[4])
text(0.3, 0, "Fens", col = cols[1], cex = 2.2)
phyd_ins(pdg_oc$SES.phydivers, oc$ins$TE, col = cols[2], -1.4, par("usr")[4])
text(0.3, 0, "Outcrops", col = cols[2], cex = 2.2)
phyd_ins(pdg_mt$SES.phydivers, mt$ins$TE, col = cols[3], 6, par("usr")[4])
text(7, -1, "Mountaintops", col = cols[3], cex = 2.2)
plot(range(c(fe$ins$TE, oc$ins$TE, mt$ins$TE)), 0:1, type = "n", ann = FALSE, axes = FALSE)
axis(1, line = -10, lwd = 0, lwd.tick = 1)
text(mean(par("usr")[1:2]), 0.85, "Insularity (target effect)")
#par(new = TRUE, mfrow = c(1,1))
#plot(range(c(fe$ins$TE, oc$ins$TE, mt$ins$TE)), 0:1, type = "n", ann = FALSE, axes = FALSE)
#legend(mean(par("usr")[1:2])-1, 0.1, pch = 1, lwd = 1, col = cols, 
#  legend = c("Fens", "Outcrops", "Mountaintops"), bty = "n", cex = 0.7)

# Q1
dclon_fe <- phylo.d(fe$traits[, c("pm_names", "clon")], phy_fe, pm_names, clon, 10000)
dclon_oc <- phylo.d(oc$traits[, c("pm_names", "clon")], phy_oc, pm_names, clon, 10000)
dclon_mt <- phylo.d(mt$traits[, c("pm_names", "clon")], phy_mt, pm_names, clon, 10000)

vars <- c("Size_of_the_belowground_bud_bank_with_root_buds_included",
  "Depth_of_the_belowground_bud_bank_with_root_buds_included",
  "Lateral_spreading_distance_by_clonal_growth")
calc_physig <- function(tr, dat, phy){
  cdat <- comparative.data(phy, dat, "pm_names", na.omit = FALSE)
  form <- formula(paste(tr, "~ 1"))
  mod <- pgls(form, cdat, lambda = "ML") 
  c(lambda = mod$param.CI$lambda$opt, P = mod$param.CI$lambda$bounds.p[1])
}
tr_fe <- sapply(vars, calc_physig, fe$traits, phy_fe)
tr_oc <- sapply(vars, calc_physig, oc$traits, phy_oc)

# Table S3
round(cbind(c(dclon_fe$DEstimate, dclon_fe$Pval1, as.vector(tr_fe)), 
  c(dclon_oc$DEstimate, dclon_oc$Pval1, as.vector(tr_oc))), 3) 
round(c(dclon_mt$DEstimate, dclon_mt$Pval1), 3)

# H3
plr <- function(dat, phy, vars, onlyClon = FALSE){
  preds <- "clon"
  if (!onlyClon) preds <- c(preds, vars[1:2])
  form <- paste("Specialist ~ ", paste(preds, collapse = " + "))
  dat$traits$Specialist <- dat$splist$Specialist[match(dat$traits$pm_names, dat$splist$pm_names)]
  rownames(dat$traits) <- dat$traits$pm_names
  dat$traits <- dat$traits[rowSums(is.na(dat$traits[, c("Specialist", preds)])) == 0, ]
  phy <- keep.tip(phy, dat$traits$pm_names)
  binaryPGLMM(form, data = dat$traits, phy = phy)
}
sp_fe <- plr(fe, phy_fe, vars)
sp_oc <- plr(oc, phy_oc, vars)
sp_mt <- plr(mt, phy_mt, vars, onlyClon = TRUE)

calc_com <- function(tr, dat, phy){
  pres <- dat$pres
  splist <- dat$splist
  traits <- dat$traits
  if (!tr %in% colnames(traits)) return(NULL)
  bin <- tr == "clon"
  sp_dat <- traits[!is.na(traits[, tr]), ]
  rownames(sp_dat) <- sp_dat$pm_names
  sp_phy <- keep.tip(phy, sp_dat$pm_names)
  form <- as.formula(paste0(tr, "~ 1"))
  if (!bin) {
    cdat <- comparative.data(sp_phy, sp_dat, "pm_names", na.omit = FALSE)
    sp_mod <- pgls(form, cdat, lambda = "ML")
    mod_out <- c(sp_mod$param.CI$lambda$opt, sp_mod$param.CI$lambda$bounds.p[1])
  } else {
    sp_mod <- binaryPGLMM(form, data = sp_dat, phy = sp_phy)
    mod_out <- c(s2 = sp_mod$s2, Pr = sp_mod$P.H0.s2)
  }
  means_com <- function(pres_com){
    com_sp <- intersect(names(pres_com)[pres_com > 0], 
      splist$pm_names[splist$Specialist > 0.5])
    com_dat <- traits[traits$pm_names %in% com_sp, ]
    com_dat <- com_dat[!is.na(com_dat[, tr]), ]
    rownames(com_dat) <- com_dat$pm_names
    com_phy <- keep.tip(phy, com_dat$pm_names)
    com_null <- mean(com_dat[, tr])
    if (com_null %in% 0:1) com_est <- com_null else {
      if (!bin) {
      cdat <- comparative.data(com_phy, com_dat, "pm_names")
      com_est <- coef(pgls(form, cdat, lambda = mod_out[1]))
      } else
        com_est <- inv_logit(binaryPGLMM(form, data = com_dat, phy = com_phy, 
          maxit.reml = 0, s2.init = mod_out[1])$B)
    }
    c(com_null, com_est)
  }
  means <- apply(pres[, -1], 1, means_com)
  rownames(means) <- c("NoPhysig", "EstPhysig")
  list(means = means, physig = mod_out)
}

vars <- c("Size_of_the_belowground_bud_bank_with_root_buds_included",
  "Depth_of_the_belowground_bud_bank_with_root_buds_included",
  "Lateral_spreading_distance_by_clonal_growth")
est_fe <- lapply(c("clon", vars), calc_com, fe, phy_fe)
est_oc <- lapply(c("clon", vars), calc_com, oc, phy_oc)
est_fe[[4]]$means <- est_fe[[4]]$means * 100
est_oc[[4]]$means <- est_oc[[4]]$means * 100
est_mt <- lapply(c("clon", vars), calc_com, mt, phy_mt)

# Figure S2
aux_lm_plot <- function(x, ys, col){
  y <- ys[1, ]
  y_phy <- ys[2, ]
  apply(rbind(x, ys), 2, function(xx) lines(rep(xx[1], nrow(ys)), sort(xx[1:nrow(ys) + 1]), col = col))
  points(x, y, col = col, pch = 16)
  points(x, y_phy, col = col, pch = 4)
  new_x <- seq(min(x), max(x), length.out = 100)
  mod <- lm(y ~ x)
  mod_phy <- lm(y_phy ~ x)
  lines(new_x, coef(mod)[1] + coef(mod)[2] * new_x, col = col, lwd = 2)
  lines(new_x, coef(mod_phy)[1] + coef(mod_phy)[2] * new_x, col = col, lwd = 2, lty = 2)
  mod_pf <- NULL
  if (nrow(ys) == 3) {
    y_phyfull <- ys[3, ]
    points(x, y_phyfull, col = col, pch = 5)
    mod_pf <- lm(y_phyfull ~ x)
    lines(new_x, coef(mod_pf)[1] + coef(mod_pf)[2] * new_x, col = col, lwd = 2, lty = 3)
  }
  invisible(list(NoPhysig = summary(mod), EstPhysig = summary(mod_phy)))
}
aux_traitplot <- function(pars_fe, pars_oc, pars_mt, ylab){
  aux_mt <- mt$ins$TE
  if (is.null(pars_mt)) aux_mt <- NULL
  plot(range(c(fe$ins$TE, oc$ins$TE, aux_mt)), range(c(pars_fe$means, pars_oc$means, pars_mt$means)), type = "n", axes = FALSE, xlab = "Insularity (target effect)", ylab = ylab)
  out_fe <- aux_lm_plot(fe$ins$TE, pars_fe$means, col = 4)
  out_oc <- aux_lm_plot(oc$ins$TE, pars_oc$means, col = 3)
  out_mt <- if (!is.null(pars_mt)) aux_lm_plot(mt$ins$TE, pars_mt$means, col = 1) else NULL
  box(bty = "l")
  axis(1)
  axis(2)
  list(FE = out_fe, OC = out_oc, MT = out_mt)
}

# png("figures/FigS2_traits.png", height = 480*10, width = 480*10, res = 72*12)
layout(matrix(c(1,1,3,3,5,2,2,4,4,5), nrow = 5))
par(mai = c(0.7,0.7,0.1,0.1))
mod_out <- Map(aux_traitplot, est_fe, est_oc, est_mt,
  c("Proportion of clonal plants", "Number of buds", "Depth of buds [cm]", "Lateral spread [cm]"))
par(mai = c(0,0,0,0))
plot(0:1, 0:1, type = "n", axes = FALSE, ann = FALSE)
legend(0.3,1, pch = c(16, 4), legend = c("None", "Estimated"), title = "Phylogenetic effect", bty = "n", lty = 1:3)
legend(0.55,1, fill = c(4,3,1), legend = c("Fens", "Outcrops", "Mountaintops"), title = "Archipelago", bty = "n")

ext <- lapply(mod_out, lapply, sapply, function(x) round(c(coef(x)["x", c("Estimate", "Pr(>|t|)")], adjRsq = x$adj.r.squared), 3))
do.call(cbind, lapply(ext, function(x) x$FE))
do.call(cbind, lapply(ext, function(x) x$OC))
do.call(cbind, lapply(ext, function(x) x$MT))

#_
