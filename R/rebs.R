################################################################################
# GOA rougheye and blackspotted rockfish: reproduce figures and results for the
# Tier 4 Tech Memo
################################################################################

# setup ---

# devtools::install_github("noaa-afsc/tier4tools")
library(tier4tools)
library(ggplot2)
library(patchwork)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# SPR inputs ----

ages <- 3:52
use_plus_group <- TRUE

# Recruitment proportions, base case is 50/50
rec_prop <- c(RE = 0.5, BS = 0.5)

# Natural mortality (assumed identical for both species)
M_val <- 0.042

# Functional maturity-at-age vectors (ages 3-52)
maturity_RE <- c(
  0.004, 0.003, 0.003, 0.003, 0.003, 0.004, 0.005, 0.009, 0.017, 0.037,
  0.085, 0.175, 0.297, 0.412, 0.503, 0.569, 0.616, 0.649, 0.672, 0.688,
  0.699, 0.705, 0.708, 0.708, 0.707, 0.705, 0.702, 0.698, 0.693, 0.688,
  0.683, 0.677, 0.671, 0.664, 0.658, 0.651, 0.644, 0.637, 0.630, 0.622,
  0.615, 0.607, 0.599, 0.591, 0.583, 0.575, 0.567, 0.559, 0.551, 0.543
)

maturity_BS <- c(
  0.013, 0.013, 0.013, 0.013, 0.014, 0.014, 0.016, 0.018, 0.021, 0.024,
  0.030, 0.037, 0.045, 0.054, 0.066, 0.079, 0.093, 0.108, 0.124, 0.139,
  0.154, 0.169, 0.182, 0.195, 0.208, 0.220, 0.232, 0.243, 0.255, 0.266,
  0.276, 0.286, 0.294, 0.302, 0.308, 0.314, 0.319, 0.323, 0.327, 0.331,
  0.334, 0.338, 0.341, 0.345, 0.349, 0.353, 0.358, 0.364, 0.370, 0.370
)

# Shared selectivity-at-age vector (ages 3-52)
selex_shared <- c(
  0.001, 0.001, 0.003, 0.007, 0.014, 0.028, 0.048, 0.066, 0.064,
  0.057, 0.052, 0.076, 0.203,
  rep(1, length(ages) - 13)
)

# Growth (von Bertalanffy parameters in cm)
vb_RE <- list(type = "vb", Linf = 54.4, t0 = -0.36, k = 0.102)
vb_BS <- list(type = "vb", Linf = 52.6, t0 =  0.20, k = 0.062)

# Weight-length parameters (W = alpha * L^beta)
wl_RE <- list(type = "wl", alpha = 1.14e-5, beta = 3.098)
wl_BS <- list(type = "wl", alpha = 8.13e-6, beta = 3.177)

inp <- spr_input(
  ages = ages,
  species = list(
    RE = list(
      len_at_age = vb_RE,
      wt_at_age = wl_RE,
      maturity = list(type = "vector", m = maturity_RE),
      selectivity = list(type = "vector", s = selex_shared),
      M = M_val
    ),
    BS = list(
      len_at_age = vb_BS,
      wt_at_age = wl_BS,
      maturity = list(type = "vector", m = maturity_BS),
      selectivity = list(type = "vector", s = selex_shared),
      M = M_val
    )
  ),
  rec_prop = c(RE = 0.5, BS = 0.5),
  use_plus_group = TRUE
)

# Check implied schedules (single combined panel)
inp_plots <- plot_spr_inputs(inp)

# select a group of plots using optional helper function that relies on
# the "patchwork" R library
grid_spr_inputs(inp_plots, order = c("mat_selex", "wt"), ncol = 1, guides = "keep") +
  patchwork::plot_annotation(tag_levels = 'A')
ggsave(filename = file.path(out_dir, "rebs_inputs.png"),
       width = 6,
       height = 8,
       units = "in",
       dpi = 300)

# run spr ----

spr_uncon <- run_spr(
  inp,
  spr_targets = c(0.40, 0.35),
  multispecies_constraint = "none",
  diagnostics = TRUE
)

spr_con <- run_spr(
  inp,
  spr_targets = c(0.40, 0.35),
  multispecies_constraint = "all_species",
  diagnostics = TRUE
)

spr_uncon$F_spr_total
# F40_total F35_total
# 0.048     0.059
spr_con$F_spr_total_constrained
# F40_total_constrained F35_total_constrained
# 0.031                 0.036

plot_spr_curves(spr_con, which = "total", digits = 3) +
plot_spr_curves(spr_uncon, which = "species", digits = 3) +
  patchwork::plot_layout(ncol = 1, guides = "keep") +
  patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(out_dir, "rebs_spr.png"),
       width = 6,
       height = 8,
       units = "in",
       dpi = 300)

# decomposition of SBPR by age ---
plots <- plot_spr_decomp(
  spr_con,
  which_panels = c("contrib", "removed", "removed_prop", "survivorship"),
  drop_plus_from_plot = TRUE, # default
  include_plus_caption = TRUE # default
)

## individual plots
# plots$contrib
# plots$removed
# plots$removed_prop
# plots$survivorship

# select a group of plots using optional helper function that relies on the
# "patchwork" R library
grid_spr_decomp(plots[c("contrib", "survivorship")], ncol = 1) +
  patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(out_dir, "rebs_decomp_F40_contrib_surv.png"),
       width = 6,
       height = 8,
       units = "in",
       dpi = 300)

grid_spr_decomp(plots[c("removed", "removed_prop")], ncol = 1) +
  patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(out_dir, "rebs_decomp_F40_removals.png"),
       width = 6,
       height = 8,
       units = "in",
       dpi = 300)

# sensitivity to R0 proportions ----

prop_grid <- seq(0.2, 0.8, by = 0.05)

# unconstrained
sens_uncon <- spr_sensitivity(
  x_base = inp,
  prop = prop_grid,
  spr_targets = c(0.40, 0.35),
  multispecies_constraint = "none"
)

sens_con <- spr_sensitivity(
  x_base = inp,
  prop = prop_grid,
  spr_targets = c(0.40, 0.35),
  multispecies_constraint = "all_species"
)

range(sens_con$F40_total)
range(sens_con$F35_total)

p_recprop <- plot_recprop_sensitivity(
  sens_uncon = sens_uncon,
  sens_con = sens_con,
  targets = c(0.40, 0.35),
  include_yield_share = TRUE
)
sens_con |>
  dplyr::filter(prop == 0.5) |>
  dplyr::select(share40_RE, share35_RE,
                share40_RE_constrained, share35_RE_constrained)
p_recprop$main + p_recprop$yield_share +
  patchwork::plot_layout(ncol = 1, guides = "collect") +
  patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(out_dir, "rebs_sensitivity.png"),
       width = 7,
       height = 8,
       units = "in",
       dpi = 300)
