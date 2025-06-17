library(mixmustr)
library(ggplot2)
mus <- tracer_parameters$mus
out <- compare_mixing_proportions(
  synthetic_df_divergent, synthetic_df_convergent, mus
)
ggsave("man/figures/comp_fig.png", out, dev = "png")
