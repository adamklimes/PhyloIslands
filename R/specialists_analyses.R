## specialists_analyses
# follows script "data_simulation.R"
library(phytools)
library(geiger)
library(caper)

## data_preparation
inv_logit <- function(x) exp(x)/(1 + exp(x))
phy_fe <- keep.tip(pm_tree, fe$splist$pm_names)
phy_oc <- keep.tip(pm_tree, oc$splist$pm_names)
phy_mt <- keep.tip(pm_tree, mt$splist$pm_names)
phy_fe$node.label <- NULL
phy_oc$node.label <- NULL
phy_mt$node.label <- NULL

## analyses
# H1
aux_PGLMM <- function(splist, phy){
  dat <- data.frame(Specialist = splist$Specialist, 
    row.names = splist$pm_names)
  mod <- binaryPGLMM(Specialist ~ 1, data = dat, phy = phy)
  list(mod = mod, p_Spec = round(c(nonPhyl = mean(dat$Specialist), 
    Phyl = inv_logit(mod$B)), 3))
}

mod_PGLMM_fe <- aux_PGLMM(fe$splist, phy_fe)
mod_PGLMM_oc <- aux_PGLMM(oc$splist, phy_oc)
mod_PGLMM_mt <- aux_PGLMM(mt$splist, phy_mt)

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
mkdelta_oc <- fit_mkdelta(oc$splist, phy_oc)
mkdelta_mt <- fit_mkdelta(mt$splist, phy_mt)

# Figure 1
# requires function "trait.plot2"
cols <- data.frame(cols = c("white","red"), row.names = 0:1)
aux_plot <- function(phy, dat, cols, leg = FALSE){
  par(mai = c(0,0,0,0))
  t_spe <- data.frame(Specialist = dat$splist$Specialist, row.names = dat$splist$pm_names)
  cla <- dat$splist$family[match(phy$tip.label, dat$splist$pm_names)]
  cla[!cla %in% names(sort(table(cla), decreasing = TRUE)[1:15])] <- NA
  trait.plot2(phy, t_spe, cols, legend = FALSE, 
    class = cla, margin = 0.45, cex.lab = 0.7, g_lwd = 3, w = 1/30)
  if (leg) legend("bottomleft", fill = cols$cols, 
    legend = c("Nonspecialist","Specialist"), bty = "n")
}

cols2 <- c("blue", "green", "black")
# png("figures/Fig1_spec.png", width = 480*10, height = 480*10*3, res = 72*15)
par(mfrow = c(3,1))
cols <- data.frame(cols = c("white",cols2[1]), row.names = 0:1)
aux_plot(phy_fe, fe, cols)
cols <- data.frame(cols = c("white",cols2[2]), row.names = 0:1)
aux_plot(phy_oc, oc, cols)
cols <- data.frame(cols = c("white",cols2[3]), row.names = 0:1)
aux_plot(phy_mt, mt, cols)
legend("bottomleft", fill = cols2, 
  legend = c("Fen specialist", "Outcrop specialist", "Mountaintop specialist"), bty = "n")

# H2
calc_phydist <- function(x, dat, phy, splist, R = 9999){
  phylodist_funs <- function(sel_spec, phy){
    phy <- keep.tip(phy, sel_spec)
    dmat <- cophenetic(phy)
    dmat_lt <- dmat[lower.tri(dmat)]
    c(phydivers = sum(phy$edge.length), phydist = sum(dmat_lt), 
      phydistm = mean(dmat_lt),  
      meannn = mean(apply(cophenetic(phy), 1, function(x) min(x[x != 0])))) #mean nearest neighbour
  }
  phy_all <- keep.tip(phy, colnames(dat)[x > 0.5])
  sel_spec <- phy_all$tip.label[phy_all$tip.label %in% 
    splist$pm_names[splist$Specialist > 0.5]]
  selections <- sapply(1:R, function(x) sample(phy_all$tip.label, length(sel_spec)))  
  rowSums(apply(selections, 2, phylodist_funs, phy) < phylodist_funs(sel_spec, phy)) / (R + 1)
}
aux_calc_phydist <- function(dat, phy){
  data.frame(t(apply(dat$pres[, -1], 1, calc_phydist, dat$pres[, -1], 
    phy, dat$splist)))
}
pd_fe <- aux_calc_phydist(fe, phy_fe)
pd_oc <- aux_calc_phydist(fe, phy_oc)
pd_mt <- aux_calc_phydist(fe, phy_mt)

phyd_ins <- function(phyd, ins, col, xt = NULL, yt){
  points(phyd ~ ins, col = col)
  ins_aux <- ins
  mod <- lm(phyd ~ ins_aux)
  new_x <- seq(min(ins), max(ins), length.out = 100)
  pred <- predict(mod, list(ins_aux = new_x), interval = "confidence")
  lines(new_x, pred[, 1], col = col, lwd = 2)
  lines(new_x, pred[, 2], col = col, lty = 2)
  lines(new_x, pred[, 3], col = col, lty = 2)
  aux <- round(c(summary(mod)$coefficients[2, 4], summary(mod)$adj.r.squared), 2)
  if (!is.null(xt)) text(xt, yt, col = col,
    labels = bquote(atop(P-value:~.(aux[1]),Adj-R^2:~.(aux[2]))))
  invisible(summary(mod))
}

# Table S1
extr_fun <- function(x, ins){
  mod <- summary(lm(x ~ ins))
  c(mod$coefficients[c(2,8)], mod$adj.r.squared)
}
round(rbind(apply(pd_fe, 2, extr_fun, fe$ins$TE),
  apply(pd_oc, 2, extr_fun, oc$ins$TE),
  apply(pd_mt, 2, extr_fun, mt$ins$TE)), 3)

# Figure 2
cols <- c("blue", "green", "black")
# png("figures/Fig2_phydivers.png", height = 480*10, width = 480*10, res = 72*10)
par(mai = c(0.85,0.85,0.1,0.1))
plot(c(fe$ins$TE, oc$ins$TE, mt$ins$TE), 
  c(pd_fe$phydivers, pd_oc$phydivers, 
  pd_mt$phydivers), type = "n", axes = FALSE, 
  xlab = "Insularity (target effect)", 
  ylab = "Phylogentic diversity (percentile)")
phyd_ins(pd_fe$phydivers, fe$ins$TE, col = cols[1], 7.5, 0.5)
phyd_ins(pd_oc$phydivers, oc$ins$TE, col = cols[2], 6, 0.7)
phyd_ins(pd_mt$phydivers, mt$ins$TE, col = cols[3], 3, 0.9)
box(bty = "l")
axis(1)
axis(2)
legend("topright", pch = 1, lwd = 1, col = cols, 
  legend = c("Fens", "Outcrops", "Mountaintops"), bty = "n")

# Number of species, specialists and effects of insularity
calc_nspec <- function(pres, sp_list){
  t(apply(pres[, -1], 1, function(x, sp_list, x_names) 
    c(sum(x == 1), sum(x == 1 & 
    sp_list$Specialist[match(x_names, sp_list$pm_names)])),
  sp_list, colnames(pres)[-1]))
}
aux_lm <- function(x, y, ...){
  mod <- lm(y ~ x)
  new_x <- seq(min(x), max(x), length.out = 100)
  plot(x, y, axes = FALSE, ...)
  lines(new_x, coef(mod)[1] + new_x * coef(mod)[2])
  box(bty = "l")
  axis(1)
  axis(2)
  summary(mod)
}
nspec_fe <- calc_nspec(fe$pres, fe$splist)
nspec_oc <- calc_nspec(oc$pres, oc$splist)
nspec_mt <- calc_nspec(mt$pres, mt$splist)

aux_lm(fe$ins$TE, nspec_fe[, 1], xlab = "Target effect", ylab = "Number of species")
aux_lm(fe$ins$TE, nspec_fe[, 2], xlab = "Target effect", ylab = "Number of specialists")
aux_lm(fe$ins$TE, nspec_fe[, 2]/nspec_fe[, 1], xlab = "Target effect", ylab = "Proportion of specialists")

aux_lm(oc$ins$TE, nspec_oc[, 1], xlab = "Target effect", ylab = "Number of species")
aux_lm(oc$ins$TE, nspec_oc[, 2], xlab = "Target effect", ylab = "Number of specialists")
aux_lm(oc$ins$TE, nspec_oc[, 2]/nspec_oc[, 1], xlab = "Target effect", ylab = "Number of specialists")

aux_lm(mt$ins$TE, nspec_mt[, 1], xlab = "Target effect", ylab = "Number of species")
aux_lm(mt$ins$TE, nspec_mt[, 2], xlab = "Target effect", ylab = "Number of specialists")
aux_lm(mt$ins$TE, nspec_mt[, 2]/nspec_mt[, 1], xlab = "Target effect", ylab = "Proportion of specialists")

# H3 and H4
calc_com <- function(tr, dat, phy){
  pres <- dat$pres
  splist <- dat$splist
  traits <- dat$traits
  bin <- tr == "clon"
  onlyclon <- tr == "Lateral_spreading_distance_by_clonal_growth"
  sp_spec <- splist$pm_names[splist$Specialist > 0.5]
  if (onlyclon) sp_spec <- sp_spec[which(traits$clon[match(sp_spec, traits$pm_names)] > 0.5)]
  sp_dat <- traits[traits$pm_names %in% sp_spec, ]
  sp_dat <- sp_dat[!is.na(sp_dat[,tr]), ]
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
    com_sp <- intersect(names(pres_com)[pres_com > 0], sp_spec)
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

clon_fe <- calc_com("clon", fe, phy_fe)
clon_oc <- calc_com("clon", oc, phy_oc)
clon_mt <- calc_com("clon", mt, phy_mt)
vars <- c("Size_of_the_belowground_bud_bank_with_root_buds_included",
  "Depth_of_the_belowground_bud_bank_with_root_buds_included",
  "Lateral_spreading_distance_by_clonal_growth")
est_fe <- lapply(vars, calc_com, fe, phy_fe)
est_oc <- lapply(vars, calc_com, oc, phy_oc)
est_fe[[3]]$means <- est_fe[[3]]$means * 100
est_oc[[3]]$means <- est_oc[[3]]$means * 100

# Figures 3 and S1
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
aux_traitplot <- function(pars_fe, pars_oc, ylab){
  plot(range(c(fe$ins$TE, oc$ins$TE)), range(c(pars_fe$means, pars_oc$means)), type = "n", axes = FALSE, xlab = "Insularity (target effect)", ylab = ylab)
  out_fe <- aux_lm_plot(fe$ins$TE, pars_fe$means, col = "blue")
  out_oc <- aux_lm_plot(oc$ins$TE, pars_oc$means, col = "green")
  box(bty = "l")
  axis(1)
  axis(2)
  list(FE = out_fe, OC = out_oc)
}

# png("figures/Fig3_clon.png", height = 480*10, width = 480*10, res = 72*10)
par(mai = c(0.85,0.85,0.1,0.1))
plot(range(c(fe$ins$TE, oc$ins$TE, mt$ins$TE)), 
  range(c(clon_fe$means, clon_oc$means, clon_mt$means)*100), 
  type = "n", axes = FALSE, xlab = "Insularity (target effect)", 
  ylab = "Clonality [%]")
aux_lm_plot(fe$ins$TE, clon_fe$means * 100, col = "blue")
aux_lm_plot(mt$ins$TE, clon_mt$means * 100, col = "black")
aux_lm_plot(oc$ins$TE, clon_oc$means * 100, col = "green")
box(bty = "l")
axis(1)
axis(2)
legend("bottomright", pch = c(16, 4), legend = c("None", "Estimated"), title = "Phylogenetic effect", bty = "n", lty = 1:2)
legend(6.5,35, fill = c("blue", "green", "black"), legend = c("Fens", "Outcrops", "Mountaintops"), title = "Archipelago", bty = "n")

# png("figures/FigS1_traits.png", height = 480*10, width = 480*10, res = 72*13)
par(mai = c(0.7,0.7,0.1,0.1), mfrow = c(2,2))
mod_out <- Map(aux_traitplot, est_fe, est_oc, 
  c("Number of buds", "Depth of buds [cm]", "Lateral spread [cm]"))
par(new = TRUE, mfrow = c(1,1))
plot(0:1, 0:1, type = "n", axes = FALSE, ann = FALSE)
legend(0.65,0.15, pch = c(16, 4), legend = c("None", "Estimated"), title = "Phylogenetic effect", bty = "n", lty = 1:3, cex = 0.7)
legend(0.65,0.3, fill = c("blue", "green"), legend = c("Fens", "Outcrops"), title = "Archipelago", bty = "n", cex = 0.7)

ext <- lapply(mod_out, lapply, sapply, function(x) round(c(coef(x)["x", c("Estimate", "Pr(>|t|)")], x$adj.r.squared), 3))
do.call(cbind, lapply(ext, function(x) x$FE))
do.call(cbind, lapply(ext, function(x) x$OC))

#_
