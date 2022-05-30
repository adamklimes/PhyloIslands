## data_simulation
# very simple simulation of required data structure (no phylogenetic signal)
library(ape)

simul_data <- function(nspec = 200, nisl = 20){
  spec <- paste0("Spec_", 1:nspec)
  fam <- sample(LETTERS[1:10], nspec, replace = TRUE)
  t_spe <- rbinom(nspec, 1, 0.2)
  t_clo <- rbinom(nspec, 1, 0.5)
  t_bud_d <- rnorm(nspec, 4)
  t_bud_n <- rnorm(nspec, 22, 3)
  t_las <- rnorm(nspec, 0.07, 0.02)
  t_las[t_clo == 0] <- NA
  splist <- data.frame(pm_names = spec, Specialist = t_spe, family = fam)
  traits <- data.frame(pm_names = spec, clon = t_clo, 
    Depth_of_the_belowground_bud_bank_with_root_buds_included = t_bud_d,
    Size_of_the_belowground_bud_bank_with_root_buds_included = t_bud_n,
    Lateral_spreading_distance_by_clonal_growth = t_las)
  pres <- data.frame(lapply(spec, function(x, nisl) rbinom(nisl, 1, 0.3), nisl))
  colnames(pres) <- spec 
  ins <- data.frame(TE = rnorm(nisl, 5))
  list(splist = splist, traits = traits, pres = pres, ins = ins)
}
set.seed(21)
fe <- simul_data()
oc <- simul_data()
mt <- simul_data()
pm_tree <- rcoal(max(length(fe$splist$pm_names), length(oc$splist$pm_names), 
  length(mt$splist$pm_names)))
pm_tree$tip.label <- paste0("Spec_", seq_along(pm_tree$tip.label))
#_
