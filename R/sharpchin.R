################################################################################
# GOA sharpchin rockfish: reproduce figures and results for the Tier 4 Tech Memo
################################################################################

# setup ---

# devtools::install_github("noaa-afsc/tier4tools")
library(tier4tools)
library(ggplot2)
library(patchwork)

out_dir <- "results"   # change as desired
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# SPR inputs ----

ages <- 0:100

inp <- spr_input(
  ages = ages,
  species = list(
    Sharpchin = list(
      len_at_age = list(type = "vb", Linf = 35.02, k = 0.122, t0 = -0.75),
      wt_at_age = list(type = "wl", alpha = 1.19e-5, beta = 3.056),
      maturity = list(type = "logistic", a50 = 10.06, delta = -1.96),
      selectivity = list(type = "knife_edge", a_c = 10),
      M = 0.059
    )
  ),
  use_plus_group = FALSE
)

# Check implied schedules (single combined panel)
inp_plots <- plot_spr_inputs(inp)

# select a group of plots using optional helper function that relies on
# the "patchwork" R library
grid_spr_inputs(inp_plots, order = c("mat_selex", "wt"), ncol = 1, guides = "keep") +
  patchwork::plot_annotation(tag_levels = 'A')
ggsave(filename = file.path(out_dir, "sharpchin_inputs.png"),
       width = 6,
       height = 8,
       units = "in",
       dpi = 300)

# run spr ----
spr_out <- run_spr(
  inp,
  spr_targets = c(0.40, 0.35),
  multispecies_constraint = "none",
  diagnostics = TRUE
)

spr_out$F_spr_total
# F40_total F35_total
# 0.066     0.080

plot_spr_curves(spr_out)

ggsave(filename = file.path(out_dir, "sharpchin_spr.png"),
       width = 6,
       height = 4,
       units = "in",
       dpi = 300)


# decomposition of SBPR by age ---
plots <- plot_spr_decomp(spr_out,
                         which_panels = c("contrib", "removed", "survivorship"),
                         drop_plus_from_plot = TRUE)

plots
# select a group of plots using optional helper function that relies on
# the "patchwork" R library
grid_spr_decomp(plots, guides = "collect") +
  patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(out_dir, "sharpchin_decomp_F40.png"),
       width = 6,
       height = 9,
       units = "in",
       dpi = 300)

# sensitivity (M x knife-edge ac) and heatmaps ----

sens <- spr_sensitivity(
  inp,
  M = seq(0.05, 0.12, by = 0.01),
  selex_ac = seq(5, 15, by = 1),
  spr_targets = c(0.40),
  multispecies_constraint = "none"
)

range(sens$F40_total)

p1 <- plot_sensitivity_heatmap(
  sens,
  x = "selex_ac",
  y = "M",
  fill = "F40_total",
  title = "Sensitivity of F40 to M and selectivity assumptions"
)

p2 <- plot_sensitivity_heatmap(
  sens,
  x = "selex_ac",
  y = "M",
  fill = "YPR_at_F40",
  title = "Sensitivity of YPR evaluated at F40"
)

p1 + p2 + patchwork::plot_layout(ncol = 1) + patchwork::plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(out_dir, "sharpchin_sensitivity.png"),
       width = 6,
       height = 8,
       units = "in",
       dpi = 300)
