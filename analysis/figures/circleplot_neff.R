#### 0. set-up ####
if (!require("packcircles")) { install.packages("packcircles") }
if (!require("data.table")) { install.packages("packcircles") }
if (!require("dplyr")) { install.packages("packcircles") }
if (!require("openxlsx")) { install.packages("openxlsx") }
if (!require("ggplot2")) { install.packages("ggplot2") }
if (!require("here")) { install.packages("ggplot2") }

#### 1. circle plot ####
main <- readWorkbook(here("analysis/cohort_summary/cohort_summary.xlsx"), sheet = "main")

packing <- circleProgressiveLayout(main$neff, sizetype='area')
data <- cbind(main, packing)
dat.gg <- circleLayoutVertices(packing, npoints=50)
dat.gg$ancestry <- rep(main$ancestry, each = 51)

p <- dat.gg %>%
  mutate(
    ancestry = toupper(ancestry),
    ancestry = factor(ancestry, levels = c("EUR", "AFR", "EAS", "AMR", "SAS")),
  ) %>%
  ggplot() + 
  # Make the bubbles
  geom_polygon(aes(x, y, group = id, fill = ancestry), colour = "black", alpha = 0.6) +
  # Add text in the center of each bubble + control its size
  geom_text(data = data, aes(x, y, size = neff, label = toupper(ancestry)), show.legend = F) +
  scale_size_continuous(range = c(1,4)) +
  scale_fill_manual(values = c("#a6cee3", "#AB5C79", "#b2df8a", "#33a02c", "#1f78b4")) +
  # General theme:
  theme_void() + 
  coord_equal() +
  guides(fill = guide_legend(title="Ancestry"))

height <- 12
width <- height * 0.8224638
ggsave(
  plot = p,
  filename = paste0(here(
    "analysis/figures/circleplot_neff.png"
  )),
  width = width,
  height = height,
  units = "cm",
  dpi = 400,
  bg = "transparent"
)








