#!/usr/bin/env Rscript
# ============================================================================
# theme_graphpad.R  —  shared GraphPad Prism-like theme, palettes & helpers
# ----------------------------------------------------------------------------
# Sourced by plot_figures.R so every figure uses ONE reusable theme and ONE
# colour system (no duplication).
#
# GraphPad Prism conventions implemented here:
#   * solid, bold BLACK axis lines on x and y (open L-shape, no panel box)
#   * tick marks pointing OUTWARD, no grid lines, white background
#   * bold, legible sans-serif axis titles / tick labels sized for print
#   * colourblind-safe (Okabe-Ito) high-contrast palette shared everywhere
#   * helpers for mean/median +/- error, significance brackets, Wilson CIs
# ============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ----------------------------------------------------------------------------
# COLOUR SYSTEM  (Okabe-Ito: safe for deuteran/protan/tritan colour vision)
# ----------------------------------------------------------------------------
gp_platform_cols <- c(
  "ONT"      = "#0072B2",  # blue
  "Illumina" = "#D55E00",  # vermillion
  "PacBio"   = "#009E73"   # bluish green
)

# tint a colour toward white by fraction p (0 = base, 1 = white); grDevices only
.gp_tint <- function(col, p) {
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb((1 - p) * m[1] + p, (1 - p) * m[2] + p, (1 - p) * m[3] + p)
}
# n distinguishable shades within one platform hue (light -> base)
.gp_shades <- function(base, n) {
  if (n <= 1) return(unname(base))
  grDevices::colorRampPalette(c(.gp_tint(base, 0.55), base))(n)
}

# condition-level palette: shades grouped within each platform's hue so the
# platform identity stays readable even across many library-prep conditions.
.gp_ont <- c("ONT Barcode", "ONT Adapter+Barcode", "ONT Rapid Barcode",
             "ONT Ligation Adapter", "ONT Rapid Adapter", "ONT PCR Adapter",
             "ONT cDNA Adapter")
.gp_ill <- c("Illumina Nextera", "Illumina TruSeq", "Illumina SmallRNA",
             "Illumina BGI/DNBSEQ", "Illumina i7 index+adapter",
             "Illumina i5 index+adapter")
.gp_pb  <- c("PacBio barcode+adapter", "PacBio SMRTbell adapter",
             "PacBio C2 primer")
gp_condition_cols <- c(
  setNames(.gp_shades(gp_platform_cols[["ONT"]],      length(.gp_ont)), .gp_ont),
  setNames(.gp_shades(gp_platform_cols[["Illumina"]], length(.gp_ill)), .gp_ill),
  setNames(.gp_shades(gp_platform_cols[["PacBio"]],   length(.gp_pb)),  .gp_pb)
)

# semantic palettes (shared so every panel speaks the same visual language)
gp_source_cols   <- c("Real" = "#333333", "Decoy" = "#BDBDBD")
gp_match_cols    <- c("Perfect (100%)" = "#333333", "Partial (\u226595% id)" = "#BDBDBD",
                      "Perfect" = "#333333", "Partial" = "#BDBDBD")
gp_position_cols <- c("Terminal\n(\u2264500 bp from end)" = "#D55E00",
                      "Internal" = "#7F7F7F")
gp_strand_cols   <- c("forward" = "#0072B2", "reverse" = "#E69F00")

# ----------------------------------------------------------------------------
# THE THEME
# ----------------------------------------------------------------------------
theme_graphpad <- function(base_size = 11, base_family = "sans") {
  ln <- 0.6  # axis / tick line weight (pt) — reads crisp at journal column width
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # --- open, bold axis lines; explicitly NO panel box / grid ---
      axis.line        = element_line(colour = "black", linewidth = ln,
                                      lineend = "square"),
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      # --- ticks pointing OUTWARD ---
      axis.ticks        = element_line(colour = "black", linewidth = ln),
      axis.ticks.length = unit(3.2, "pt"),
      # --- bold, legible sans-serif type ---
      axis.title    = element_text(face = "bold", size = base_size, colour = "black"),
      axis.title.x  = element_text(margin = margin(t = 4)),
      axis.title.y  = element_text(margin = margin(r = 4), angle = 90),
      axis.text     = element_text(size = base_size - 1, colour = "black"),
      # --- panel labels / titles: bold, top-left, consistent ---
      plot.title    = element_text(face = "bold", size = base_size + 1, hjust = 0,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 2.5, colour = "grey35", hjust = 0,
                                   margin = margin(b = 6)),
      plot.tag      = element_text(face = "bold", size = base_size + 4, hjust = 0),
      plot.tag.position = c(0, 1),
      # --- facets: clean strip, no box ---
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size - 1, colour = "black"),
      # --- legend: no key box, compact ---
      legend.position  = "bottom",
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      plot.margin      = margin(8, 12, 8, 8)
    )
}

# ----------------------------------------------------------------------------
# STAT / ANNOTATION HELPERS
# ----------------------------------------------------------------------------
# Wilson score interval for a binomial proportion (no package dependency).
# Returns a data.frame with p, lo, hi (all as fractions in [0, 1]).
gp_wilson <- function(x, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- ifelse(n > 0, x / n, NA_real_)
  denom  <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  data.frame(p = p, lo = pmax(0, centre - half), hi = pmin(1, centre + half))
}

# GraphPad-style significance code from a p-value.
gp_p_stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 1e-4, "****",
      ifelse(p < 1e-3, "***",
        ifelse(p < 1e-2, "**",
          ifelse(p < 0.05, "*", "ns")))))
}

# GraphPad-style significance bracket. `drop` is a multiplicative tick length
# so it renders correctly on BOTH linear and log10 y-axes (y assumed > 0).
add_sig_bracket <- function(x1, x2, y, label, drop = 0.04, size = 3,
                            family = "sans", lwd = 0.5) {
  y0 <- y * (1 - drop)
  list(
    annotate("segment", x = x1, xend = x2, y = y,  yend = y,  linewidth = lwd),
    annotate("segment", x = x1, xend = x1, y = y,  yend = y0, linewidth = lwd),
    annotate("segment", x = x2, xend = x2, y = y,  yend = y0, linewidth = lwd),
    annotate("text", x = (x1 + x2) / 2, y = y, label = label,
             vjust = -0.25, size = size, family = family)
  )
}

# Median + IQR overlay (GraphPad "median with interquartile range").
gp_median_iqr <- function(colour = "black", size = 0.45, width = 0.0) {
  stat_summary(fun = median,
               fun.min = function(z) stats::quantile(z, 0.25),
               fun.max = function(z) stats::quantile(z, 0.75),
               geom = "pointrange", colour = colour, linewidth = 0.6,
               size = size)
}

# ----------------------------------------------------------------------------
# NESTING COLLAPSE (mirrors scripts/nesting.py :: collapse_nested_events)
# ----------------------------------------------------------------------------
# Within each group defined by `group_cols`, drop hit loci that are STRICTLY
# contained inside a LONGER construct's subject footprint (a bare tag sitting
# inside its adapter-embedded form), crediting the physical site to the longest
# construct. Equal-footprint / distinct-site hits are preserved, so legitimate
# cross-barcode multiplicity is untouched. `df` must carry sstart, send and the
# grouping columns. Returns the surviving rows (a subset of `df`). This makes the
# platform-pooled per-hit figures show the same nesting sensitivity as the
# Python three-way review scripts (report_addition18-20).
gp_collapse_nested <- function(df, group_cols = c("platform", "saccver"),
                               sstart = "sstart", send = "send") {
  n <- nrow(df)
  if (n < 2) return(df)
  lo <- pmin(df[[sstart]], df[[send]])
  hi <- pmax(df[[sstart]], df[[send]])
  gkey <- do.call(paste, c(df[group_cols], sep = "\r"))

  # order by group, locus start asc, locus end desc
  o  <- order(gkey, lo, -hi)
  g  <- gkey[o]; L <- lo[o]; H <- hi[o]
  new_g <- c(TRUE, g[-1] != g[-n])            # group starts
  gid   <- cumsum(new_g)

  # running max of hi within each group, and its lag (max over PREVIOUS rows)
  cmax_incl <- stats::ave(H, gid, FUN = cummax)
  cmax_prev <- c(-Inf, cmax_incl[-n])
  cmax_prev[new_g] <- -Inf

  # lo-blocks within a group (a block = one shared start coordinate)
  new_blk <- new_g | (L != c(L[1], L[-n]))
  blk     <- cumsum(new_blk)
  blk_first <- which(new_blk)                  # first row index of each block
  prevmax   <- cmax_prev[blk_first][blk]       # max hi of any strictly-smaller-lo hit
  blkmax    <- stats::ave(H, blk, FUN = max)   # max hi among same-lo hits

  # drop iff strictly contained: a smaller-lo hit reaches past hi, OR a same-lo
  # hit has a strictly larger hi. Equal-footprint duplicates are preserved.
  drop_sorted <- (prevmax >= H) | (blkmax > H)
  drop <- logical(n); drop[o] <- drop_sorted
  df[!drop, , drop = FALSE]
}

# GraphPad-style pairwise Mann-Whitney (Wilcoxon rank-sum) significance brackets
# for a numeric `value` split by a factor `group`. Compares the given pairs
# (list of length-2 character vectors) and returns a list of add_sig_bracket()
# annotations stacked at geometrically increasing heights (works on log10 y).
# `x_index` maps group level -> x position (defaults to factor level order).
gp_mwu_brackets <- function(value, group, pairs, y_start, step = 1.9,
                            x_index = NULL, size = 3.4, show_ns = TRUE) {
  group <- as.character(group)
  if (is.null(x_index)) {
    lv <- levels(factor(group)); x_index <- setNames(seq_along(lv), lv)
  }
  anns <- list(); y <- y_start
  for (pr in pairs) {
    a <- value[group == pr[1]]; b <- value[group == pr[2]]
    if (length(a) < 2 || length(b) < 2) next
    p <- tryCatch(suppressWarnings(stats::wilcox.test(a, b, exact = FALSE)$p.value),
                  error = function(e) NA_real_)
    lab <- gp_p_stars(p)
    if (lab == "ns" && !show_ns) next
    anns <- c(anns, add_sig_bracket(x_index[[pr[1]]], x_index[[pr[2]]], y,
                                    lab, size = size))
    y <- y * step
  }
  anns
}


# ----------------------------------------------------------------------------
# PUBLICATION-QUALITY OUTPUT DEVICES
# ----------------------------------------------------------------------------
# Shadow ggsave so every figure in both scripts is written with a device that
# embeds fonts and renders full Unicode (x, >=, en/em dashes, >>) instead of the
# base pdf()/png() devices that silently transliterate them:
#   * .pdf  -> cairo_pdf  (vector, UTF-8, font-embedded)
#   * .png  -> ragg::agg_png  (crisp raster text)
# Falls back to the ggplot2 default device when cairo/ragg are unavailable.
ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  device <- NULL
  if (grepl("\\.pdf$", filename, ignore.case = TRUE) && capabilities("cairo")) {
    device <- grDevices::cairo_pdf
  } else if (grepl("\\.png$", filename, ignore.case = TRUE) &&
             requireNamespace("ragg", quietly = TRUE)) {
    device <- ragg::agg_png
  }
  ggplot2::ggsave(filename = filename, plot = plot, ..., device = device)
}
