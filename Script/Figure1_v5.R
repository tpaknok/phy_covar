# Figure 1 (conceptual): triangle -> sympatry -> pairwise -> community
# Requires: ggplot2, patchwork, ape
library(ggplot2)
library(patchwork)
library(ape)

set.seed(11)
HI  <- "#1f4e9c"   # high native distribution overlap
LO  <- "#e07b00"   # low native distribution overlap
NEU <- "grey25"
NOTE_COL <- "grey40"

# ---------------------------------------------------------------- helpers ----
# Species pairs on a triangle: close relatives always similar,
# distant relatives can be similar or very different. `top` = upper-edge slope.
triangle <- function(n = NULL, top = 0.95, x = NULL, skew = 3) {
  if (is.null(x)) x <- runif(n, 0.03, 1)
  y <- x * top * runif(length(x))^skew + rnorm(length(x), 0, 0.012)
  data.frame(x = x, y = pmax(y, 0.01))
}

tri_poly <- function(top) data.frame(x = c(0, 1, 1), y = c(0, 0, top))

theme_concept <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(linewidth = 0.6,
                                   arrow = arrow(length = unit(2, "mm"), type = "closed")),
          plot.title = element_text(hjust = 0.5, size = base + 0.5),
          axis.title.x = element_text(margin = margin(t = 7)),
          axis.title.y = element_text(margin = margin(r = 7)),
          plot.margin = margin(4, 8, 4, 4))
}

lims <- list(scale_x_continuous(limits = c(0, 1.08), expand = c(0, 0)),
             scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)),
             coord_cartesian(clip = "off"))

row_title <- function(txt) {
  wrap_elements(grid::textGrob(txt, x = 0.01, hjust = 0,
                               gp = grid::gpar(fontface = "bold", fontsize = 12)))
}

# Rectangular, root-at-top drawing of an ultrametric tree as segments.
# Tips are aligned at `ybot`; a younger tree (smaller root depth) is shorter.
tree_segments <- function(tr, x0, w, ybot, hscale) {
  tr <- ladderize(tr)
  n <- Ntip(tr); depth <- node.depth.edgelength(tr); H <- max(depth)
  xpos <- numeric(n + tr$Nnode)
  tips_in_order <- tr$edge[tr$edge[, 2] <= n, 2]
  xpos[tips_in_order] <- seq_len(n)
  for (nd in rev(sort(unique(tr$edge[, 1])))) {          # children before parents
    ch <- tr$edge[tr$edge[, 1] == nd, 2]
    xpos[nd] <- mean(range(xpos[ch]))
  }
  X <- function(v) x0 + (v - 1) / (n - 1) * w
  Y <- function(d) ybot + (H - d) * hscale
  vert <- data.frame(x = X(xpos[tr$edge[, 2]]), xend = X(xpos[tr$edge[, 2]]),
                     y = Y(depth[tr$edge[, 1]]), yend = Y(depth[tr$edge[, 2]]))
  horiz <- do.call(rbind, lapply(unique(tr$edge[, 1]), function(nd) {
    ch <- tr$edge[tr$edge[, 1] == nd, 2]
    data.frame(x = X(min(xpos[ch])), xend = X(max(xpos[ch])),
               y = Y(depth[nd]), yend = Y(depth[nd]))
  }))
  root <- n + 1
  stem <- data.frame(x = X(xpos[root]), xend = X(xpos[root]),
                     y = Y(depth[root]), yend = Y(depth[root]) + 0.03)
  rbind(vert, horiz, stem)
}

# Two example communities (branch lengths in units where the species pool's root = 1).
# Low PD variance: near-star (very rapid radiation) -> all species almost equally related.
tree_low <- read.tree(text = paste0(
  "(((A:0.68,B:0.68):0.01,(C:0.68,(D:0.67,E:0.67):0.01):0.01):0.01,",
  "((F:0.68,G:0.68):0.01,H:0.69):0.01);"))
# High PD variance: deep split plus recent radiations -> close AND distant relatives.
tree_high <- read.tree(text = paste0(
  "((((A:0.08,B:0.08):0.04,(C:0.06,D:0.06):0.06):0.50,(E:0.10,F:0.10):0.52):0.38,",
  "((G:0.07,H:0.07):0.08,I:0.15):0.85);"))
pair_pd <- function(tr) { d <- cophenetic(tr); d[upper.tri(d)] / 2 }  # pool max = 1

# ======================================================= (A) the triangle ====
dA <- triangle(500, 0.95)
pA1 <- ggplot(dA, aes(x, y)) +
  geom_polygon(data = tri_poly(0.95), fill = NEU, alpha = 0.08) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 0.95, linetype = "dashed", colour = NEU) +
  geom_point(colour = NEU, alpha = 0.7, size = 1.4) +
  annotate("segment", x = 0.12, y = 0.48, xend = 0.08, yend = 0.06, colour = "grey60", linewidth = 0.3) +
  annotate("text", x = 0.05, y = 0.52, label = "close relatives:\nalways similar",
           hjust = 0, vjust = 0, size = 3.5, colour = "grey35") +
  annotate("segment", x = 0.7, y = 0.82, xend = 0.92, yend = 0.5, colour = "grey60", linewidth = 0.3) +
  annotate("text", x = 0.3, y = 0.86, label = "distant relatives:\nsimilar or very different",
           hjust = 0, vjust = 0, size = 3.5, colour = "grey35") +
  lims + labs(x = "Phylogenetic distance", y = "Ecological difference") + theme_concept()

pA2 <- ggplot() +
  annotate("segment", x = 0.05, xend = 0.95, y = 0.5, yend = 0.5, linewidth = 2.2,
           arrow = arrow(length = unit(3.5, "mm"), type = "closed")) +
  annotate("text", x = 0.5, y = 0.68, label = "Statistical\nmodel", size = 3.5, vjust = 0) +
  annotate("text", x = 0.5, y = 0.3, label = "(points =\nspecies pairs)", size = 3.5,
           colour = NOTE_COL, vjust = 1) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

mini_line <- function(slope, ttl, ylab, tag) {
  ggplot() +
    annotate("segment", x = 0, xend = 1, y = 0.08, yend = 0.08 + slope,
             colour = NEU, linewidth = 1.2) +
    annotate("text", x = 0.95, y = 0.12 + slope, label = tag, hjust = 1, vjust = 0,
             size = 3.5, colour = NOTE_COL) +
    lims + labs(title = ttl, x = "Phylogenetic\ndistance", y = ylab) + theme_concept()
}
pA3 <- mini_line(0.08, "Mean", "Mean ecological\ndifference", "weak")
pA4 <- mini_line(0.85, "Dispersion", "Variation among\npairs", "strong")
rowA <- pA1 + pA2 + pA3 + pA4 + plot_layout(widths = c(1.35, 0.45, 0.8, 0.8))

# ================================================= (B) role of sympatry =====
dBl <- triangle(500, 0.5)
dBh <- triangle(500, 1.0)
pB1 <- ggplot(dBl, aes(x, y)) +
  geom_polygon(data = tri_poly(0.5), fill = LO, alpha = 0.12) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 0.5, colour = LO) +
  geom_point(colour = LO, alpha = 0.75, size = 1.4) +
  annotate("text", x = 0.04, y = 1.03, label = "Limited past competition:\nlimited divergence",
           hjust = 0, vjust = 1, size = 3.5, colour = "grey30") +
  lims + labs(title = "Low native distribution overlap", x = "Phylogenetic distance",
              y = "Ecological difference") +
  theme_concept() + theme(plot.title = element_text(colour = LO))

pB2 <- ggplot(dBh, aes(x, y)) +
  geom_polygon(data = tri_poly(1.0), fill = HI, alpha = 0.12) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1.0, colour = HI) +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0.5, colour = LO, linetype = "dotted",
           linewidth = 0.7) +
  annotate("text", x = 1.01, y = 0.5, label = "low-overlap\nedge", hjust = 0, size = 3.5,
           colour = LO) +
  geom_point(colour = HI, alpha = 0.75, size = 1.4) +
  annotate("text", x = 0.04, y = 1.03,
           label = "Character displacement:\ndivergence amplified,\nespecially among\ndistant relatives",
           hjust = 0, vjust = 1, size = 3.5, colour = "grey30") +
  lims + labs(title = "High native distribution overlap", x = "Phylogenetic distance",
              y = "Ecological difference") +
  theme_concept() + theme(plot.title = element_text(colour = HI),
                          plot.margin = margin(4, 40, 4, 4))
rowB <- pB1 + pB2

# ============================================ (C) pairwise analyses =========
# Example population dynamics: a shared environmental driver; close relatives
# respond alike (synchronous), distant relatives respond differently (asynchronous).
tt  <- 1:30
env <- as.numeric(stats::filter(rnorm(40), rep(1 / 4, 4), sides = 1))[11:40]
env <- (env - mean(env)) / sd(env)
ts_df <- rbind(
  data.frame(t = tt, sp = "Species 1", pair = "Close relatives",
             N = 6 + env + rnorm(30, 0, 0.25)),
  data.frame(t = tt, sp = "Species 2", pair = "Close relatives",
             N = 2.5 + 0.9 * env + rnorm(30, 0, 0.25)),
  data.frame(t = tt, sp = "Species 1", pair = "Distant relatives",
             N = 6 + env + rnorm(30, 0, 0.25)),
  data.frame(t = tt, sp = "Species 2", pair = "Distant relatives",
             N = 2.5 - 0.8 * env + rnorm(30, 0, 0.35)))
ts_df$pair <- factor(ts_df$pair, levels = c("Close relatives", "Distant relatives"),
                     labels = c("Ecologically similar pair: synchronous",
                                "Ecologically dissimilar pair: asynchronous"))

pC1 <- ggplot(ts_df, aes(t, N, linetype = sp)) +
  geom_line(colour = NEU, linewidth = 0.55) +
  facet_wrap(~pair, ncol = 1) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  labs(x = "Time", y = "Abundance") +
  theme_concept() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10, hjust = 0, colour = "grey20"),
        legend.position = "bottom", legend.text = element_text(size = 10),
        legend.key.width = unit(6, "mm"), legend.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(4, "mm"))

pair_pred <- function(hi, lo, ttl, ylab) {
  d <- data.frame(x = c(0, 1, 0, 1), y = c(0.1, 0.1 + hi, 0.1, 0.1 + lo),
                  grp = factor(rep(c("High overlap", "Low overlap"), each = 2),
                               levels = c("High overlap", "Low overlap")))
  ggplot(d, aes(x, y, colour = grp)) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = c(HI, LO), name = NULL) +
    lims + labs(title = ttl, x = "Phylogenetic\ndistance", y = ylab) +
    theme_concept() +
    theme(legend.position = c(0.04, 0.99), legend.justification = c(0, 1),
          legend.text = element_text(size = 10), legend.key.width = unit(4, "mm"),
          legend.key.height = unit(3, "mm"), legend.background = element_blank())
}
pC2 <- pair_pred(0.12, 0.04, "Mean", "Interspecific asynchrony /\nstability difference")
pC3 <- pair_pred(0.88, 0.25, "Dispersion", "Variation among pairs") +
  theme(legend.position = "none")
rowC <- pC1 + pC2 + pC3 + plot_layout(widths = c(1.15, 0.85, 0.85))

# =================================== (D) PD variance within communities =====
bg <- triangle(500, 0.95)
community_plot <- function(tr, ttl) {
  pd <- pair_pd(tr)
  dc <- triangle(top = 0.95, x = pd + rnorm(length(pd), 0, 0.008))  # slight jitter so tied distances stay visible
  segs <- tree_segments(tr, x0 = 0.05, w = 0.36, ybot = 0.6, hscale = 0.38)
  ggplot(bg, aes(x, y)) +
    geom_polygon(data = tri_poly(0.95), fill = NEU, alpha = 0.05) +
    annotate("segment", x = 0, y = 0, xend = 1, yend = 0.95, linetype = "dashed", colour = NEU) +
    geom_point(colour = "grey82", size = 0.9) +
    geom_point(data = dc, colour = NEU, size = 1.9, alpha = 0.9) +
    geom_segment(data = segs, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = "grey15", linewidth = 0.45, lineend = "square") +
    lims + labs(title = ttl, x = "Phylogenetic distance", y = "Ecological difference") +
    theme_concept()
}
pD1 <- community_plot(tree_low, "Low phylogenetic\ndistance variance")
pD2 <- community_plot(tree_high, "High phylogenetic\ndistance variance")

comm_pred <- function(hi, lo, ttl, ylab, legend = TRUE) {
  d <- data.frame(x = c(0, 1, 0, 1), y = c(0.1, 0.1 + hi, 0.1, 0.1 + lo),
                  grp = factor(rep(c("High overlap", "Low overlap"), each = 2),
                               levels = c("High overlap", "Low overlap")))
  ggplot(d, aes(x, y, colour = grp)) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = c(HI, LO), name = NULL) +
    lims + labs(title = ttl, x = "Phylogenetic distance\nvariance", y = ylab) +
    theme_concept() +
    theme(legend.position = if (legend) c(0.04, 0.99) else "none",
          legend.justification = c(0, 1), legend.text = element_text(size = 7),
          legend.key.width = unit(4, "mm"), legend.key.height = unit(3, "mm"),
          legend.background = element_blank(),
          axis.title.x = element_text(size = 8, lineheight = 0.95, margin = margin(t = 7)),
          plot.margin = margin(4, 10, 4, 4))
}
pD3 <- comm_pred(0.12, 0.04, "Mean", "Mean phylogenetic signal\nin population stability")
pD4 <- comm_pred(0.88, 0.25, "Dispersion", "Variation in phylogenetic\nsignal among communities",
                 legend = FALSE)
rowD <- pD1 + pD2 + pD3 + pD4 + plot_layout(widths = c(1, 1, 0.78, 0.78))

# ================================================================ assemble ===
fig <- row_title("(A) Phylogenetic distance and ecological divergence") / rowA /
  row_title("(B) Role of sympatry (native distribution overlap as a proxy)") / rowB /
  row_title("(C) Pairwise analyses: interspecific asynchrony and population stability") / rowC /
  row_title("(D) Community analyses: phylogenetic distance variance") / rowD +
  plot_layout(heights = c(0.09, 1, 0.09, 1.05, 0.09, 1.05, 0.09, 1.05))

ggsave("Figure/Figure1_v5_R.tif", fig, width = 9, height = 13, dpi = 600,compression = "lzw")

