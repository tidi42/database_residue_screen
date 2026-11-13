suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(jsonlite)
  library(digest)
  library(grid)
})
options(warn = 2)
# Grobs are built on a null PDF device (same metrics as the default device), so no Rplots.pdf is written to the working directory.
grDevices::pdf(NULL)

# Paper 3 publication figures from one figure-source directory (contract publication_figure_sources_v1, written by the
# figure_sources stage). Main figures: figure1 (framework and evidence states), figure2 (matched reference and control
# recovery), figure3 (positional context), figure4 (a priori testability and record cases); supplementary figures: figureS1
# (query properties and alignment quality), figureS2 (Monte Carlo precision), figureS3 (positional robustness).
# Result panels show testable queries only: untestable and chance-limited parents and lengths are left out.
# The record-case figure (figureS4) is drawn by plot_record_cases.R. No panel contains a table: lookup-style content is
# exported as CSV (tables/ inside the output directory) and every such panel is drawn from the exported CSV.
# Artwork conventions: 180 mm wide, at most 240 mm tall, 300 dpi PNG and PDF; platform colours (ONT blue, Illumina vermilion,
# BGI/DNBSEQ pink, PacBio green); Unicode glyphs for >=, <= and the multiplication sign; square unit brackets on axis titles;
# US spelling; serial commas; en dashes in numeric ranges; thousands separators in counts. Union rows cover every condition.
# R object names stay figure1-figure4; the manuscript numbers them Figures 2-5 (the graphical abstract is Figure 1).
# Exported CSV names carry the manuscript numbers (e.g. fig2B_* for figure1 panel B, figS4B_* for figureS2, tableS6_* for Table S6).
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% c(3, 4))
input <- normalizePath(arguments[1])
theme_file <- normalizePath(arguments[2])
output <- arguments[3]
selection <- if (length(arguments) == 4) arguments[4] else "all"
manifest <- fromJSON(file.path(input, "manifest.json"))
stopifnot(manifest$contract == "publication_figure_sources_v1")
for (filename in names(manifest$files)) {
  stopifnot(!grepl("(^/|(^|/)\\.\\.(/|$))", filename),
            digest(file.path(input, filename), algo = "sha256", file = TRUE) == manifest$files[[filename]])
}
source(theme_file)
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
tables_directory <- file.path(output, "tables")
if (!dir.exists(tables_directory)) dir.create(tables_directory, recursive = TRUE)
figure_width_mm <- 180
figure_dpi <- 300
groups <- c("ONT", "Illumina", "BGI/DNBSEQ", "PacBio")
colours <- c(gp_platform_cols, "BGI/DNBSEQ" = "#CC79A7")
draws_per_parent <- manifest$counts$draws_per_parent
draw_prefixes <- as.numeric(unlist(manifest$draw_prefixes))
final_prefix <- max(draw_prefixes)
target_limit <- manifest$target_limit
parent_label <- function(id) gsub("_", " ", sub("_forward$", "", id))
# Count labels for breaks that survive the data limits (1, 100, 10k, 1M).
short_count <- function(x) ifelse(is.na(x), NA_character_, ifelse(x >= 1e6, paste0(x / 1e6, "M"), ifelse(x >= 1000, paste0(x / 1000, "k"),
  format(x, trim = TRUE, scientific = FALSE))))
# Platform colours (ONT blue, Illumina vermillion, BGI/DNBSEQ pink, PacBio green, as in Figure 2) are used for every
# element that is split by platform, including per-platform heatmaps (colour intensity = count) and the Illumina match markers.
# Categories that are not platforms use a separate family with no blue, orange, pink or green:
# burgundy (E = 10, rejected, mismatch, failure), olive (accepted, positive, included), slate (E = 1e12, non-positive,
# query and scaffold) and brown (excluded). Sequential heatmaps that are not platform-coded use a sand-to-burgundy ramp (S2B) or an olive ramp (S3C).
ink <- "#1A1A1A"
grey_mid <- "#7F7F7F"
grey_light <- "#D9D9D9"
burgundy <- "#8E2043"
burgundy_light <- "#F2DCE2"
olive <- "#6E7B12"
olive_light <- "#E6EAC2"
olive_pale <- "#F3F5E1"
slate <- "#3E4E5C"
slate_mid <- "#8795A1"
slate_light <- "#DDE3E8"
brown <- "#8C510A"
sand_light <- "#F3E6C4"
role_colours <- c(barcode = burgundy, scaffold = slate, unknown = grey_mid, linker = olive, flank = "#8C8C8C")
role_fills <- c(barcode = burgundy_light, scaffold = slate_light, unknown = "#F2F2F2", linker = olive_light, flank = "white")
sand_burgundy_ramp <- c("#FBF7EA", "#EEDDB0", "#D2B27A", "#A5764F", "#7A3B3A", "#420F1E")
olive_ramp <- c("#FBF9EC", "#EEE8B8", "#CDC677", "#949637", "#5B6420", "#2E360F")
theme_set(theme_graphpad(10) + theme(
  axis.text = element_text(size = 8.5, face = "plain"),
  axis.title = element_text(size = 9),
  plot.title = element_text(size = 11, face = "bold"),
  plot.title.position = "plot",
  plot.subtitle = element_text(size = 8.5, colour = "black"),
  strip.text = element_text(size = 9),
  legend.text = element_text(size = 8.5),
  legend.title = element_text(size = 8.5),
  legend.position = "bottom", panel.spacing = unit(3, "mm"),
  plot.margin = margin(5, 5, 5, 5)))
load_table <- function(name) read_tsv(file.path(input, paste0(name, ".tsv")),
                                     show_col_types = FALSE, na = c("", "not_estimated"))
group_order <- function(data) mutate(data, group = factor(group, levels = groups))
pt <- function(points) points / ggplot2::.pt
serial <- function(values) {
  values <- as.character(values)
  if (length(values) <= 2) paste(values, collapse = " and ") else paste0(paste(values[-length(values)], collapse = ", "), ", and ", values[length(values)])
}
independent_panel <- function(plot) wrap_elements(full = ggplotGrob(plot))
with_margin <- function(plot, top = 5, bottom = 5) plot + theme(plot.margin = margin(top, 5, bottom, 5))
panel_separator <- function() wrap_elements(full = linesGrob(x = unit(c(0, 1), "npc"), y = unit(0.5, "npc"),
  gp = gpar(col = "grey80", lwd = 0.6)))
separator_height <- 1.5
# patchwork gives every stacked element fixed top and bottom plot margins (and adds them once more around the figure); the
# relative heights share the rest. The reference layout is 240 mm tall with heights summing to 240 for a given number of stacked
# elements; shorter stacks keep the same size per height unit, so a panel is drawn at the same size in every layout.
stack_margin_mm <- sum(as.numeric(convertHeight(theme_get()$plot.margin[c(1, 3)], "mm")))
stack_height <- function(heights, reference_elements) {
  unit_mm <- (240 - (reference_elements + 1) * stack_margin_mm) / 240
  round((length(heights) + 1) * stack_margin_mm + sum(heights) * unit_mm, 2)
}

# Write a table to tables/ (never overwriting different content) and return the values read back,
# so that the panels are built from exactly what is exported.
export_csv <- function(data, name) {
  path <- file.path(tables_directory, name)
  temporary <- tempfile(fileext = ".csv")
  write_csv(data, temporary, na = "NA")
  if (file.exists(path)) {
    stopifnot(digest(temporary, algo = "sha256", file = TRUE) == digest(path, algo = "sha256", file = TRUE))
  } else {
    stopifnot(file.copy(temporary, path))
  }
  unlink(temporary)
  read_csv(path, show_col_types = FALSE, na = "NA")
}

# Panels drawn in millimetre coordinates (x 0-180, y downwards from the top of the panel).
canvas <- function(height, layers) {
  wrap_elements(full = ggplotGrob(ggplot() + layers +
    scale_x_continuous(limits = c(0, figure_width_mm), expand = c(0, 0)) +
    scale_y_reverse(limits = c(height, 0), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA))))
}
mm_text <- function(x, y, label, size = 8.5, face = "plain", hjust = 0, colour = "black", family = "sans", vjust = 0.5) {
  annotate("text", x = x, y = y, label = label, size = pt(size), fontface = face, hjust = hjust,
           vjust = vjust, colour = colour, family = family, lineheight = 1)
}
# Species names (Latin binomials, flag from the record-case verification) are set in italics with plotmath. A leading comma of the upright remainder is
# set with the italic name, so the plotmath italic correction falls inside the following space instead of before the comma.
mm_text_species <- function(x, y, organism, italic, rest, size = 8.5) {
  if (!isTRUE(italic)) return(mm_text(x, y, paste0(organism, rest), size))
  quote <- function(text) paste0('"', gsub('"', '\\\\"', text), '"')
  comma <- startsWith(rest, ",")
  label <- paste0("italic(", quote(paste0(organism, if (comma) "," else "")), ")*", quote(if (comma) substring(rest, 2) else rest))
  annotate("text", x = x, y = y, label = label, parse = TRUE, size = pt(size), hjust = 0, vjust = 0.5, colour = "black", family = "sans")
}
mm_rect <- function(xmin, xmax, ymin, ymax, fill = "white", colour = NA, linewidth = 0.3, linetype = "solid") {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, colour = colour,
           linewidth = linewidth, linetype = linetype)
}
mm_segment <- function(x, y, xend, yend, colour = "grey30", linewidth = 0.4, head = FALSE, linetype = "solid", lineend = "butt") {
  arguments <- list("segment", x = x, y = y, xend = xend, yend = yend, colour = colour, linewidth = linewidth,
                    linetype = linetype, lineend = lineend)
  if (head) arguments$arrow <- arrow(length = unit(1.4, "mm"), type = "closed", angle = 22)
  do.call(annotate, arguments)
}
panel_title <- function(label, size = 11) mm_text(0.5, 3.6, label, size = size, face = "bold")

saved <- list()
# Manuscript numbering: the graphical abstract is Figure 1, so the R figures figure1-figure4 are manuscript Figures 2-5.
# Supplementary figures are numbered in the order of first citation: R figureS2 is Figure S4 and the record cases
# (R figureS4, plot_record_cases.R) are Figure S2. File names keep the R names.
manuscript_numbers <- c(figure1 = "Figure 2", figure2 = "Figure 3", figure3 = "Figure 4", figure4 = "Figure 5",
                        figureS1 = "Figure S1", figureS2 = "Figure S4", figureS3 = "Figure S3")
save_figure <- function(name, plot, height, panels, source_tables, note, csv_tables = character(), legend_notes = character()) {
  if (selection != "all" && selection != name) return(invisible(NULL))
  stopifnot(height <= 240 + 1e-6, name %in% names(manuscript_numbers))
  for (extension in c("pdf", "png")) {
    filename <- file.path(output, paste0("paper3_", name, ".", extension))
    stopifnot(!file.exists(filename))
    ggplot2::ggsave(filename, plot, width = figure_width_mm, height = height, units = "mm", bg = "white", dpi = figure_dpi,
      device = if (extension == "pdf") grDevices::cairo_pdf else ragg::agg_png)
  }
  saved[[name]] <<- list(manuscript_figure = unname(manuscript_numbers[name]), panels = panels, source_tables = source_tables,
                         exported_csv = csv_tables, note = note, legend_notes = legend_notes,
                         width_mm = figure_width_mm, height_mm = height, png_dpi = figure_dpi)
  message("Rendered ", name)
}

# ---------------------------------------------------------------------------------------------
# Figure 1. Framework and evidence states
# ---------------------------------------------------------------------------------------------
counts <- manifest$counts
testability <- load_table("length_testability")
production_evalue <- unique(testability$production_evalue)
stopifnot(length(production_evalue) == 1)
box_bottom <- 2.45
box_top <- 4.9
box_centre <- (box_bottom + box_top) / 2
workflow_nodes <- data.frame(
  xpos = c(1.4, 4.5, 7.6, 10.7), ypos = box_centre,
  label = c(paste0("Query panel\n", scales::comma(counts$physical_references), " queries\n", scales::comma(counts$canonical_parents),
                   " parents\n", draws_per_parent, " draws / parent"),
            paste0("Matched search\nE = ", production_evalue, ", cap ", scales::comma(target_limit), "\nTestability by length"),
            "Single-HSP support\nStrict / tolerant\nOne HSP per match", "Record context\nExact version + length\nExplicit eligibility"))
workflow <- ggplot() +
  geom_segment(aes(x = c(2.85, 5.95, 9.05), xend = c(3.05, 6.15, 9.25), y = box_centre, yend = box_centre),
               arrow = arrow(length = unit(1.4, "mm")), linewidth = 0.4) +
  geom_rect(data = workflow_nodes, aes(xmin = xpos - 1.4, xmax = xpos + 1.4, ymin = box_bottom, ymax = box_top),
            fill = "#F7F8F9", colour = slate_mid, linewidth = 0.35) +
  geom_text(data = workflow_nodes, aes(xpos, ypos, label = label), size = 2.95, lineheight = 1.15) +
  geom_segment(aes(x = c(4.5, 10.7), xend = c(4.5, 10.7), y = box_bottom, yend = 1.72),
               arrow = arrow(length = unit(1.4, "mm")), linewidth = 0.4) +
  annotate("text", x = c(4.5, 10.7), y = 1.05, colour = c(burgundy, brown), fontface = "bold",
    label = c("Untestable or incomplete\nExcluded from counts", "Excluded or unresolved\nOutside linear-genomic set"), size = 2.9) +
  annotate("text", x = 0, y = 0, hjust = 0, size = 2.9,
    label = paste0(scales::comma(counts$physical_references), " queries, ", scales::comma(counts$canonical_parents),
                   " canonical parents, ", scales::comma(counts$distinct_searches), " search sequences.")) +
  coord_cartesian(xlim = c(0, 12.2), ylim = c(-0.1, box_top + 0.25), clip = "off") + theme_void() +
  labs(title = "A  From queries to qualified evidence") +
  theme(plot.title = element_text(size = 11, face = "bold"), plot.margin = margin(5, 5, 0, 5))

# B: matching rules as a drawn schematic (one 24-nt query bar above six subject alignments)
examples <- load_table("matching_examples")
stopifnot(nrow(examples) == 6, length(unique(examples$reference)) == 1, all(nchar(examples$reference) == 24))
matching_rules <- export_csv(data.frame(
  example = examples$example, query = examples$reference, aligned_subject = examples$subject,
  strict = examples$strict, tolerant = examples$tolerant), "fig2B_matching_rules.csv")
query_sequence <- matching_rules$query[1]
query_length <- nchar(query_sequence)
parse_alignment <- function(reference, subject) {
  if (grepl(" / ", subject, fixed = TRUE)) {
    parts <- strsplit(subject, " / ", fixed = TRUE)[[1]]
    ends <- cumsum(nchar(parts))
    stopifnot(length(parts) == 2, ends[2] == nchar(reference))
    return(data.frame(kind = "aligned", start = c(0, ends[1]), end = ends, lane = 1:2))
  }
  subject_chars <- strsplit(subject, "")[[1]]
  reference_chars <- strsplit(reference, "")[[1]]
  stopifnot(length(subject_chars) == length(reference_chars))
  kind <- ifelse(subject_chars == ".", "unaligned", ifelse(subject_chars == "-", "gap", "aligned"))
  runs <- rle(kind)
  bounds <- cumsum(runs$lengths)
  parts <- data.frame(kind = runs$values, start = bounds - runs$lengths, end = bounds, lane = 0)
  mismatch <- which(kind == "aligned" & subject_chars != reference_chars)
  if (length(mismatch)) parts <- rbind(parts, data.frame(kind = "mismatch", start = mismatch - 1, end = mismatch, lane = 0))
  parts
}
rule_x0 <- 66
rule_scale <- 3
rule_rows_y <- 31 + 7.6 * (seq_len(nrow(matching_rules)) - 1)
alignment_parts <- bind_rows(lapply(seq_len(nrow(matching_rules)), function(row) {
  mutate(parse_alignment(matching_rules$query[row], matching_rules$aligned_subject[row]), row = row)
})) %>% mutate(
  ycentre = rule_rows_y[row] + ifelse(lane == 0, 0, (lane - 1.5) * 3.3),
  half = ifelse(lane == 0, 1.8, 1.35),
  xmin = rule_x0 + start * rule_scale, xmax = rule_x0 + end * rule_scale)
bars <- filter(alignment_parts, kind %in% c("aligned", "unaligned"))
bars <- bars[order(bars$kind != "unaligned"), ]
mismatches <- filter(alignment_parts, kind == "mismatch")
gaps <- filter(alignment_parts, kind == "gap")
symbol_columns <- c(strict = 152, tolerant = 168)
draw_outcome <- function(x, y, accepted) {
  if (accepted) {
    mm_segment(c(x - 1.7, x - 0.5), c(y + 0.1, y + 1.5), c(x - 0.5, x + 1.9), c(y + 1.5, y - 1.7), colour = olive, linewidth = 0.9, lineend = "round")
  } else {
    mm_segment(c(x - 1.4, x - 1.4), c(y - 1.4, y + 1.4), c(x + 1.4, x + 1.4), c(y + 1.4, y - 1.4), colour = burgundy, linewidth = 0.9, lineend = "round")
  }
}
rules_height <- 85.8
rules_panel <- canvas(rules_height, c(
  list(panel_title("B  Matching rules: conceptual 24-nt examples"),
    mm_text(1.5, 13, "Example", 8.5, "bold"), mm_text(rule_x0, 13, "Alignment to the query", 8.5, "bold"),
    mm_text(symbol_columns[["strict"]], 13, "Strict", 8.5, "bold", hjust = 0.5),
    mm_text(symbol_columns[["tolerant"]], 13, "Tolerant", 8.5, "bold", hjust = 0.5),
    mm_text(1.5, 21, paste0("Query (", query_length, " nt)"), 8.5),
    mm_rect(rule_x0, rule_x0 + query_length * rule_scale, 21 - 1.8, 21 + 1.8, fill = brown),
    mm_text(1.5, rule_rows_y, matching_rules$example, 8.5),
    mm_rect(bars$xmin, bars$xmax, bars$ycentre - bars$half, bars$ycentre + bars$half,
      fill = ifelse(bars$kind == "aligned", sand_light, "white"), colour = brown, linewidth = 0.35,
      linetype = ifelse(bars$kind == "aligned", "solid", "dashed")),
    mm_rect(mismatches$xmin, mismatches$xmax, mismatches$ycentre - mismatches$half, mismatches$ycentre + mismatches$half, fill = burgundy),
    mm_segment(rule_x0 + gaps$start * rule_scale, gaps$ycentre, rule_x0 + gaps$end * rule_scale, gaps$ycentre, colour = brown, linewidth = 0.5)),
  lapply(seq_len(nrow(matching_rules)), function(row) list(
    draw_outcome(symbol_columns[["strict"]], rule_rows_y[row], matching_rules$strict[row]),
    draw_outcome(symbol_columns[["tolerant"]], rule_rows_y[row], matching_rules$tolerant[row]))),
  list(
    mm_rect(1.5, 5, 76.6, 79.4, fill = sand_light, colour = brown, linewidth = 0.35), mm_text(6.5, 78, "Aligned", 8),
    mm_rect(28, 31.5, 76.6, 79.4, fill = "white", colour = brown, linewidth = 0.35, linetype = "dashed"), mm_text(33, 78, "Unaligned query bases", 8),
    mm_rect(75, 78.5, 76.6, 79.4, fill = burgundy), mm_text(80, 78, "Mismatch", 8),
    mm_segment(102, 78, 106.5, 78, colour = brown, linewidth = 0.5), mm_text(108, 78, "Gap within one HSP", 8),
    # One footnote line; the conceptual query sequence is not printed.
    mm_text(1.5, 83, "Two separate bars are two HSPs, each assessed alone. Strict is nested within tolerant (\u226595% identity, \u226590% query span, no gaps).",
      8, colour = "grey25"))))

# C: evidence states as a decision tree. Two a-priori testability states precede completion; both are excluded from every
# count. Thresholds are read from the exported length_testability table.
informative_limit <- unique(testability$informativeness_limit)
testable_limit <- unique(testability$testable_chance_limit)
stopifnot(length(informative_limit) == 1, length(testable_limit) == 1)
evidence_states <- export_csv(data.frame(
  evidence_state = c("Untestable (excluded)", "Chance-limited (excluded)", "Raw zero", "Endpoint non-positive", "Cap-limited non-positive",
                     "Identified positive", "Positive, unresolved breadth", "Failed or incomplete (excluded)"),
  testability = c("censored_untestable", "chance_limited", rep("testable", 6)),
  counted_as = c("excluded", "excluded", "0", "0", "0", "1", "1", "excluded"),
  search_complete = c("n/a", "n/a", "Yes", "Yes", "Yes", "Yes", "Yes", "No"),
  raw_hsps = c("n/a", "n/a", "0", ">0", ">0", ">0", ">0", "Unknown"),
  accepted_support = c("n/a", "n/a", "No", "No", "No", "Yes", "Yes, additional breadth unresolved", "Unknown"),
  target_cap_reached = c("n/a", "n/a", "No", "No", "Yes", "Either", "Either", "Unknown")), "fig2C_evidence_states.csv")
state_of <- function(name) evidence_states$evidence_state[match(name, evidence_states$evidence_state)]
tree_nodes <- data.frame(
  id = c("T1", "T2", "Q1", "Q2", "Q3", "Q4", "Q5", "U1", "U2", "L1", "L2", "L3", "L4", "L5", "L6"),
  label = c(paste0("Score retainable at E = ", production_evalue, "\nand chance bound \u2264 ", informative_limit, "?"),
            paste0("Chance bound \u2264 ", testable_limit, "?"),
            "Search complete?", "Raw HSPs > 0?", "Accepted support?", "Reached target cap?", "Additional breadth resolved?",
            state_of(c("Untestable (excluded)", "Chance-limited (excluded)", "Failed or incomplete (excluded)", "Raw zero", "Endpoint non-positive",
                       "Cap-limited non-positive", "Positive, unresolved breadth", "Identified positive"))),
  kind = c(rep("question", 7), rep("state", 8)),
  x = c(30, 30, 30, 30, 30, 84, 30, 150, 150, 150, 150, 150, 150, 150, 30),
  y = c(11.5, 24, 36, 48, 63, 63, 79, 11.5, 24, 36, 48, 57, 69, 79, 91),
  w = c(52, 52, 52, 52, 52, 40, 52, rep(56, 7), 52),
  h = c(9.5, rep(8, 14)))
stopifnot(!anyNA(tree_nodes$label))
tree_edges <- export_csv(data.frame(
  from = c("T1", "T1", "T2", "T2", "Q1", "Q1", "Q2", "Q2", "Q3", "Q3", "Q4", "Q4", "Q5", "Q5"),
  to = c("U1", "T2", "U2", "Q1", "L1", "Q2", "L2", "Q3", "Q4", "Q5", "L3", "L4", "L6", "L5"),
  answer = c("no", "yes", "no", "yes", "no", "yes", "no", "yes", "no", "yes", "no", "yes", "yes", "no")), "fig2C_decision_tree_edges.csv")
excluded_states <- c("U1", "U2", "L1")
tree_nodes$linetype <- ifelse(tree_nodes$id %in% excluded_states, "dashed", "solid")
tree_nodes$fill <- ifelse(tree_nodes$kind == "question", "white", c(
  U1 = sand_light, U2 = sand_light, L1 = "white", L2 = slate_light, L3 = slate_light, L4 = slate_light, L5 = olive_pale, L6 = olive_light)[tree_nodes$id])
tree_nodes$border <- ifelse(tree_nodes$kind == "question", grey_mid, c(
  U1 = brown, U2 = brown, L1 = burgundy, L2 = slate, L3 = slate, L4 = slate, L5 = olive, L6 = olive)[tree_nodes$id])
node_at <- function(id) tree_nodes[match(id, tree_nodes$id), ]
tree_edge_layers <- lapply(seq_len(nrow(tree_edges)), function(index) {
  edge <- tree_edges[index, ]
  from <- node_at(edge$from)
  to <- node_at(edge$to)
  if (edge$from == "Q4") {
    x_branch <- from$x + from$w / 2 + 5
    return(list(
      mm_segment(from$x + from$w / 2, from$y, x_branch, from$y, linewidth = 0.4),
      mm_segment(x_branch, from$y, x_branch, to$y, linewidth = 0.4),
      mm_segment(x_branch, to$y, to$x - to$w / 2, to$y, linewidth = 0.4, head = TRUE),
      mm_text(x_branch + 2.3, to$y + ifelse(edge$answer == "no", -2.6, 2.6), edge$answer, 8, hjust = 0.5, colour = "grey25")))
  }
  if (from$x == to$x) {
    list(mm_segment(from$x, from$y + from$h / 2, to$x, to$y - to$h / 2, linewidth = 0.4, head = TRUE),
         mm_text(from$x + 1.6, (from$y + from$h / 2 + to$y - to$h / 2) / 2, edge$answer, 8, colour = "grey25"))
  } else {
    x_start <- from$x + from$w / 2
    x_end <- to$x - to$w / 2
    list(mm_segment(x_start, from$y, x_end, to$y, linewidth = 0.4, head = TRUE),
         mm_text((x_start + x_end) / 2, from$y - 2.6, edge$answer, 8, hjust = 0.5, colour = "grey25"))
  }
})
tree_height <- 97.5
# Notes that belong in the figure legend are recorded in the artwork manifest. The panel ends 2.5 mm below the last box.
figure1_legend_notes <- c(
  "Testability is fixed a priori per query length and database size. Dashed states are excluded from every count.")
stopifnot(max(tree_nodes$y + tree_nodes$h / 2) + 2.5 == tree_height)
state_panel <- canvas(tree_height, c(
  list(panel_title("C  Testability, completion, match support, and ambiguity are different states"),
    mm_rect(tree_nodes$x - tree_nodes$w / 2, tree_nodes$x + tree_nodes$w / 2, tree_nodes$y - tree_nodes$h / 2,
            tree_nodes$y + tree_nodes$h / 2, fill = tree_nodes$fill, colour = tree_nodes$border, linewidth = 0.4, linetype = tree_nodes$linetype),
    mm_text(tree_nodes$x, tree_nodes$y, tree_nodes$label, 8.5, hjust = 0.5)),
  tree_edge_layers))

figure1_heights <- c(49.5, separator_height, rules_height, separator_height, tree_height)
figure1 <- (independent_panel(workflow) / panel_separator() / rules_panel / panel_separator() / state_panel) +
  plot_layout(heights = figure1_heights)
save_figure("figure1", figure1, stack_height(figure1_heights, 5), c("Workflow", "Matching rules schematic", "Evidence-state decision tree"),
            c("matching_examples.tsv", "manifest.json", "length_testability.tsv"),
            "Conceptual method illustration",
            c("fig2B_matching_rules.csv", "fig2C_evidence_states.csv", "fig2C_decision_tree_edges.csv"),
            legend_notes = figure1_legend_notes)

# ---------------------------------------------------------------------------------------------
# Figure 2. Matched reference and control recovery
# ---------------------------------------------------------------------------------------------
# Native testable cohorts per endpoint (testable, complete and resolved parents). Excluded testability states stay out of every count.
# Platform and global unions cover every condition.
detection <- load_table("detection") %>% filter(population == "all_records")
endpoint_values <- c("eligible_parents", "reference_positive", "reference_rate", "control_draws", "control_positive", "control_rate",
                     "censored", "chance_limited", "partially_censored", "unresolved", "incomplete", "other_inapplicable")
native <- detection %>% filter(cohort == "native")
stopifnot(all(native %>% group_by(order) %>% summarise(n = n_distinct(planned_parents), .groups = "drop") %>% pull(n) == 1))
native <- native %>% select(scope, order, block, group, label, min_length, max_length, planned_parents, endpoint, all_of(endpoint_values)) %>%
  pivot_wider(names_from = endpoint, values_from = all_of(endpoint_values))
common <- detection %>% filter(cohort == "common_strict_tolerant") %>%
  select(order, endpoint, eligible_parents, reference_positive, control_positive, control_draws) %>%
  pivot_wider(names_from = endpoint, values_from = c(eligible_parents, reference_positive, control_positive, control_draws), names_prefix = "common_")
stopifnot(identical(common$eligible_parents_common_strict, common$eligible_parents_common_tolerant))
block_levels <- c("Reference classes", "Platform unions")
classes <- native %>% left_join(common, by = "order") %>%
  mutate(block = factor(block, levels = block_levels)) %>%
  arrange(block, match(group, c(groups, "Global union")), order)
stopifnot(!anyNA(classes$block))
stopifnot(all(classes$eligible_parents_strict + classes$censored_strict + classes$chance_limited_strict + classes$partially_censored_strict +
                classes$unresolved_strict + classes$incomplete_strict + classes$other_inapplicable_strict == classes$planned_parents))
pct <- function(positive, n) ifelse(n > 0, 100 * positive / n, NA_real_)
table1 <- export_csv(data.frame(
  block = as.character(classes$block), platform = classes$group, reference_class = classes$label,
  length_nt = ifelse(classes$min_length == classes$max_length, as.character(classes$min_length), paste0(classes$min_length, "-", classes$max_length)),
  planned_parents_n = classes$planned_parents,
  strict_testable_n = classes$eligible_parents_strict, strict_censored_n = classes$censored_strict, strict_chance_limited_n = classes$chance_limited_strict,
  strict_unresolved_n = classes$unresolved_strict, strict_incomplete_n = classes$incomplete_strict,
  ref_pos_strict = classes$reference_positive_strict, ref_pct_strict = pct(classes$reference_positive_strict, classes$eligible_parents_strict),
  ctrl_draws_strict = classes$control_draws_strict, ctrl_pos_strict = classes$control_positive_strict,
  ctrl_pct_strict = pct(classes$control_positive_strict, classes$control_draws_strict),
  tolerant_testable_n = classes$eligible_parents_tolerant, tolerant_censored_n = classes$censored_tolerant,
  tolerant_chance_limited_n = classes$chance_limited_tolerant, tolerant_partially_censored_n = classes$partially_censored_tolerant,
  tolerant_unresolved_n = classes$unresolved_tolerant, tolerant_incomplete_n = classes$incomplete_tolerant,
  ref_pos_tolerant = classes$reference_positive_tolerant, ref_pct_tolerant = pct(classes$reference_positive_tolerant, classes$eligible_parents_tolerant),
  ctrl_draws_tolerant = classes$control_draws_tolerant, ctrl_pos_tolerant = classes$control_positive_tolerant,
  ctrl_pct_tolerant = pct(classes$control_positive_tolerant, classes$control_draws_tolerant),
  common_parents_n = classes$eligible_parents_common_strict,
  common_ref_pos_strict = classes$reference_positive_common_strict, common_ref_pos_tolerant = classes$reference_positive_common_tolerant,
  common_ctrl_draws = classes$control_draws_common_strict,
  common_ctrl_pos_strict = classes$control_positive_common_strict, common_ctrl_pos_tolerant = classes$control_positive_common_tolerant),
  "table1_reference_classes.csv")
stopifnot(isTRUE(all.equal(table1$ref_pct_strict / 100, ifelse(table1$strict_testable_n > 0, classes$reference_rate_strict, NA_real_))),
          isTRUE(all.equal(table1$ctrl_pct_tolerant / 100, ifelse(table1$tolerant_testable_n > 0, classes$control_rate_tolerant, NA_real_))))
n_text <- function(n) scales::comma(n)
# Panel A shows the cohorts with at least one testable parent under either rule. Table 1 keeps every class.
arrow_rows <- table1 %>% filter(strict_testable_n > 0 | tolerant_testable_n > 0) %>% mutate(order = row_number(),
  row_name = case_when(platform == "Global union" ~ "All platforms", block == "Platform unions" ~ paste0(platform, ", all classes"),
                       TRUE ~ paste(platform, reference_class, sep = " / ")),
  n_label = paste0(n_text(strict_testable_n), " / ", n_text(tolerant_testable_n)))
stopifnot(nrow(arrow_rows) > 0)
arrow_data <- export_csv(bind_rows(
  transmute(arrow_rows, order, block, platform, reference_class, row_name, n_label, series = "Reference", rule = "strict", pct = ref_pct_strict, n = strict_testable_n),
  transmute(arrow_rows, order, block, platform, reference_class, row_name, n_label, series = "Reference", rule = "tolerant", pct = ref_pct_tolerant, n = tolerant_testable_n),
  transmute(arrow_rows, order, block, platform, reference_class, row_name, n_label, series = "Control mean", rule = "strict", pct = ctrl_pct_strict, n = strict_testable_n),
  transmute(arrow_rows, order, block, platform, reference_class, row_name, n_label, series = "Control mean", rule = "tolerant", pct = ctrl_pct_tolerant, n = tolerant_testable_n)
) %>% arrange(order, rule, series), "fig3A_detection_points.csv")
# Each block after the first gets an empty slot that carries its title.
slots <- arrow_rows %>% transmute(order, block, slot = order + cumsum(block != lag(block, default = first(block))))
row_count <- max(slots$slot)
slot_y <- setNames(row_count - slots$slot + 1, slots$order)
panel_colours <- c(colours, "Global union" = ink)
same_parents <- with(arrow_rows, strict_testable_n > 0 & strict_testable_n == tolerant_testable_n & strict_testable_n == common_parents_n)
arrow_data <- arrow_data %>% mutate(
  colour_key = platform, ypos = slot_y[as.character(order)] + ifelse(series == "Reference", 0.17, -0.17),
  fill = ifelse(series == "Reference", unname(panel_colours[platform]), "white"),
  shape = ifelse(rule == "strict", 21, 24))
row_labels <- arrow_rows %>% transmute(y = unname(slot_y[as.character(order)]), row_name, n_label)
block_breaks <- arrow_rows %>% mutate(next_block = lead(block), next_platform = lead(platform)) %>%
  filter(!is.na(next_block), next_block != block | (next_platform != platform & block == "Reference classes")) %>%
  transmute(y = unname(slot_y[as.character(order)]) - 0.5, strong = next_block != block)
block_titles <- arrow_rows %>% group_by(block) %>% summarise(order = min(order), .groups = "drop") %>% filter(block != "Reference classes") %>%
  mutate(y = unname(slot_y[as.character(order)]) + 1)
arrow_plot <- ggplot(filter(arrow_data, !is.na(pct)), aes(pct, ypos, colour = colour_key)) +
  geom_hline(data = block_breaks, aes(yintercept = y, linewidth = strong), colour = "grey80", inherit.aes = FALSE) +
  scale_linewidth_manual(values = c(`TRUE` = 0.6, `FALSE` = 0.3), guide = "none") +
  geom_vline(xintercept = 50, linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_segment(data = arrow_data %>% filter(order %in% arrow_rows$order[same_parents]) %>% select(order, series, rule, pct, ypos, colour_key) %>%
                 pivot_wider(names_from = rule, values_from = pct),
               aes(x = strict, xend = tolerant, y = ypos, yend = ypos), linewidth = 0.35, colour = "grey60", inherit.aes = FALSE) +
  geom_point(aes(fill = I(fill), shape = I(shape)), size = 1.9, stroke = 0.6) +
  geom_text(data = row_labels, aes(x = 104, y = y, label = n_label), inherit.aes = FALSE, hjust = 0, size = 2.75) +
  geom_text(data = block_titles, aes(x = 1, y = y, label = block), inherit.aes = FALSE, hjust = 0, size = 2.75, fontface = "italic", colour = "grey30") +
  annotate("text", x = 104, y = row_count + 1.25, label = "n strict / tolerant", hjust = 0, size = 2.75, fontface = "bold", lineheight = 0.95) +
  scale_colour_manual(values = panel_colours, guide = "none") +
  scale_x_continuous(limits = c(-3, 150), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.5, row_count + 1.9), breaks = row_labels$y, labels = row_labels$row_name, expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "A  Strict and tolerant detection in testable cohorts",
    subtitle = paste0("Circles: strict. Triangles: tolerant. Filled: references. Open: their draws (", draws_per_parent, " per parent).\n",
      "n: testable, complete, and resolved parents per rule. Gray line: same parents under both rules. Dotted line: 50%."),
    x = "Detection [%]", y = NULL) +
  theme(axis.text.y = element_text(size = 7.6), legend.position = "none", axis.ticks.y = element_blank(),
        axis.line.x = element_line(), plot.subtitle = element_text(size = 8))

breadth <- load_table("paired_breadth") %>% filter(group != "Global union") %>% group_order()
breadth_counts <- breadth %>% count(group, reference, control_mean, name = "parents")
breadth_summary <- breadth %>% group_by(group) %>% summarise(parents = n(), zero_pairs = sum(reference == 0 & control_mean == 0), .groups = "drop")
breadth_limit <- max(log10(breadth$reference + 1)) * 1.05
control_limit <- max(breadth$control_mean) * 1.1
breadth_labels <- setNames(paste0(breadth_summary$group, "\nN=", breadth_summary$parents,
  ", both zero: ", breadth_summary$zero_pairs), as.character(breadth_summary$group))
breadth_breaks <- c(0, 1, 10, 100, 1000)
control_breaks <- scales::breaks_pretty(n = 3)(c(0, max(breadth$control_mean)))
breadth_plot <- ggplot(breadth_counts, aes(control_mean, log10(reference + 1), colour = group)) +
  geom_function(fun = function(x) log10(x + 1), xlim = c(0, control_limit), inherit.aes = FALSE,
    linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_point(aes(size = parents), alpha = 0.7) +
  facet_wrap(~group, nrow = 1, drop = FALSE, labeller = labeller(group = breadth_labels)) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_size_continuous(range = c(1.2, 4), breaks = c(1, 100, 500), name = "Parents per point") +
  scale_x_continuous(limits = c(-0.08 * control_limit, control_limit), breaks = control_breaks[control_breaks <= control_limit]) +
  scale_y_continuous(limits = c(-0.06, breadth_limit), breaks = log10(breadth_breaks + 1), labels = breadth_breaks) +
  labs(title = "B  Paired retained-accession breadth",
       subtitle = paste0("Strict, testable parents of every condition. Dashed curve: reference count equals control mean.\nY-axis: log10(count + 1). X-axis: linear, limited to the observed control range (0 to ",
         format(max(breadth$control_mean), digits = 2), ")."),
      x = "Mean control accession count", y = "Reference\naccession count") +
  theme(strip.text = element_text(size = 8), legend.margin = margin(0, 0, 0, 0))
figure2_heights <- c(min(240 - 88 - separator_height, 40 + 4.6 * row_count), separator_height, 88)
figure2 <- (independent_panel(with_margin(arrow_plot, bottom = 0)) / panel_separator() /
    independent_panel(with_margin(breadth_plot, top = 2))) + plot_layout(heights = figure2_heights)
save_figure("figure2", figure2, sum(figure2_heights), c("Strict and tolerant detection in testable cohorts", "Paired breadth"),
  c("detection.tsv", "paired_breadth.tsv"),
  "Descriptive construct-level support in rule-specific testable cohorts. Untestable and chance-limited parents are excluded",
  c("table1_reference_classes.csv", "fig3A_detection_points.csv"))

# ---------------------------------------------------------------------------------------------
# Figure 3. Positional context (strict, a-priori testable events of every condition; unit-level sensitivities)
# ---------------------------------------------------------------------------------------------
reason_levels <- c(
  linear_genomic_dna_known_exclusions_removed = "Included: linear-genomic DNA",
  record_type_transcript = "Excluded: transcript record", record_type_amplicon = "Excluded: amplicon record",
  record_type_synthetic = "Excluded: synthetic record", record_type_vector = "Excluded: vector record",
  origin_synthetic = "Excluded: synthetic origin", origin_vector = "Excluded: vector origin",
  origin_artificial_other = "Excluded: other artificial origin", origin_engineered = "Excluded: engineered origin",
  targeted_record = "Excluded: targeted record", rna_record = "Excluded: RNA", circular_record = "Excluded: circular",
  molecule_not_established = "Unresolved: molecule", record_type_not_established = "Unresolved: record type",
  topology_not_established = "Unresolved: topology", molecule_conflicting = "Unresolved: conflicting molecule",
  record_type_conflicting = "Unresolved: conflicting record type", topology_conflicting = "Unresolved: conflicting topology",
  origin_conflicting = "Unresolved: conflicting origin", targeted_conflicting = "Unresolved: conflicting targeting",
  version_mismatch = "Unresolved: version mismatch", length_conflict = "Unresolved: length conflict",
  conflicting_record_evidence = "Unresolved: conflicting evidence", unresolved_id = "Unresolved: identifier",
  identifier_not_normalised = "Unresolved: unnormalized identifier", unresolved_identifier = "Unresolved: identifier")
flow <- load_table("eligibility_flow")
stopifnot(all(flow$reason %in% names(reason_levels)))
flow <- flow %>% mutate(reason_label = factor(unname(reason_levels[reason]), levels = rev(unique(unname(reason_levels)))))
flow_totals <- flow %>% group_by(state) %>% summarise(accessions = sum(accessions), .groups = "drop")
flow_count <- function(name) scales::comma(sum(flow_totals$accessions[flow_totals$state == name]))
stopifnot(sum(flow$accessions[!(flow$state == "unresolved" & flow$reason == "unresolved_identifier")]) == counts$testable_strict_accessions)
flow_plot <- ggplot(flow, aes(accessions, reason_label, colour = state, shape = state)) +
  geom_segment(aes(x = 0, xend = accessions, yend = reason_label), linewidth = 0.6) + geom_point(size = 2) +
  geom_text(aes(label = scales::comma(accessions)), hjust = -0.25, size = 3, colour = "black") +
  scale_colour_manual(values = c(included = olive, excluded = brown, unresolved = slate), guide = "none") +
  scale_shape_manual(values = c(included = 16, excluded = 16, unresolved = 1), guide = "none") +
  scale_y_discrete(drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "A  Eligibility for linear-genomic positional analysis",
      subtitle = paste0(scales::comma(sum(flow$accessions)), " accessions with strict matches of testable queries: ", flow_count("included"),
        " included, ", flow_count("excluded"), " excluded, ", flow_count("unresolved"), " unresolved.\n",
        "Reasons are applied in order (record type, origin, targeting, RNA, circular)."),
       x = "Distinct accession.version records", y = NULL)
position_policy <- manifest$position_policy
main_window <- position_policy$primary_window
main_minimum <- position_policy$primary_minimum_length
position <- load_table("position_events") %>% group_order()
main_position <- filter(position, record_length >= main_minimum)
position_totals <- main_position %>% group_by(group) %>%
  summarise(events = n(), accessions = n_distinct(accession_version), .groups = "drop") %>%
  mutate(label = paste0(as.character(group), " (", scales::comma(events), ifelse(events == 1, " event, ", " events, "), scales::comma(accessions),
    ifelse(accessions == 1, " accession)", " accessions)")))
legend_labels <- setNames(position_totals$label, as.character(position_totals$group))
sparse_position <- main_position %>% group_by(group) %>% filter(n() <= 5) %>%
  count(group, distance, name = "events") %>% arrange(group, distance) %>%
  mutate(cumulative = cumsum(events) / sum(events)) %>% ungroup()
position_ecdf <- ggplot(main_position, aes(log10(distance + 1), colour = group)) +
  stat_ecdf(geom = "step", pad = FALSE, linewidth = 0.55) +
  geom_point(data = sparse_position, aes(y = cumulative), size = 1.8) +
  geom_vline(xintercept = log10(main_window + 1), linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  scale_colour_manual(values = colours, labels = legend_labels, name = NULL) +
  scale_x_continuous(breaks = log10(c(0, 100, 500, 10000, 1000000) + 1), labels = c("0", "100", "500", "10k", "1M")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "B  Distance to the nearest stored record end",
       subtitle = paste0("Strict matches of testable queries in eligible linear-genomic records \u2265", main_minimum / 1000, " kb. Dashed line: ",
         main_window, " nt.\nMarkers show the observed steps in groups with \u22645 events."),
       x = "Distance [nt, log10(distance + 1)]", y = "Cumulative\nevent fraction") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
sensitivity <- load_table("position_sensitivity") %>% filter(group != "Global union") %>% group_order()
variant_labels <- c(all_events = "All events", one_per_accession = "One per accession", one_per_bioproject = "One per BioProject",
                    excluding_targeted_insdc = "Without INSDC-targeted records")
variant_labels <- variant_labels[names(variant_labels) %in% unique(load_table("position_sensitivity")$variant)]
variant_shapes <- c("All events" = 16, "One per accession" = 1, "One per BioProject" = 2, "Without INSDC-targeted records" = 5)[unname(variant_labels)]
stopifnot(all(sensitivity$variant %in% names(variant_labels)))
terminal <- export_csv(sensitivity %>% filter(window == main_window, minimum_length == main_minimum) %>%
  arrange(match(as.character(group), groups), match(variant, names(variant_labels))) %>%
  transmute(group = as.character(group), variant, unit_rule, events, terminal, expected, ratio, units_without_project), "fig4C_terminal_ratios.csv")
terminal <- terminal %>% mutate(group = factor(group, levels = rev(groups)),
    variant = factor(unname(variant_labels[variant]), levels = unname(variant_labels)),
    ypos = as.numeric(group) + ((length(variant_labels) + 1) / 2 - as.numeric(variant)) * 0.8 / length(variant_labels),
    label = ifelse(is.na(ratio), ifelse(events == 0, "no eligible events", paste0("ratio undefined (", terminal, "/", events, ")")),
                   paste0(sprintf("%.1f", ratio), "x (", scales::comma(terminal), "/", scales::comma(events), "), Exp ", signif(expected, 3))))
terminal_plot <- ggplot(terminal, aes(log10(ratio + 1), ypos, colour = group, shape = variant)) +
  geom_vline(xintercept = log10(2), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_point(size = 2.1, stroke = 0.6) + geom_text(aes(x = ifelse(is.na(ratio), 0, log10(ratio + 1)), label = label), hjust = -0.12, colour = "black", size = 2.4) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_shape_manual(values = variant_shapes, name = NULL) +
  scale_x_continuous(limits = c(-0.05, 3), breaks = log10(c(0, 1, 10, 100) + 1), labels = c("0", "1", "10", "100")) +
  scale_y_continuous(breaks = seq_along(levels(terminal$group)), labels = levels(terminal$group), limits = c(0.5, length(levels(terminal$group)) + 0.5)) +
  labs(title = "C  Terminal observed / expected ratios: unit and targeted-record sensitivity",
       subtitle = paste0("W = ", main_window, " nt, records \u2265", main_minimum / 1000, " kb. One event per unit is chosen by a position-independent hash.\n",
         "Labels: ratio (terminal/all events) and Exp, the expected terminal count under uniform placement."),
      x = "Observed / expected [log10(ratio + 1)]", y = NULL) +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
figure3 <- (independent_panel(with_margin(flow_plot, bottom = 0)) / panel_separator() /
    independent_panel(with_margin(position_ecdf, top = 2, bottom = 0)) / panel_separator() /
    independent_panel(with_margin(terminal_plot, top = 2))) + plot_layout(heights = c(0.95, 0.05, 1.0, 0.05, 1.45))
save_figure("figure3", figure3, 238, c("Accession eligibility", "End-distance ECDF", "Terminal ratios by unit"),
  c("eligibility_flow.tsv", "position_events.tsv", "position_sensitivity.tsv"),
  "Strict, a priori testable events of every condition. Ordered eligibility reasons. Descriptive geometry with unit-level sensitivities",
  "fig4C_terminal_ratios.csv")


# ---------------------------------------------------------------------------------------------
# Figure 4. Short-query calibration and record context
# ---------------------------------------------------------------------------------------------
# A: a-priori testability from the production calibration (length_testability) and benchmark planted exact-spike recovery at E = 10
# for three database sizes. E-value at the minimum qualifying score decides score retainability; the iid chance bound decides
# informativeness. Both are expected chance counts on one axis; the limits are read from the table.
class_labels <- c(testable = "Testable", chance_limited = "Chance-limited", censored_untestable = "Untestable")
curves <- export_csv(testability %>% transmute(rule = endpoint, length_nt = query_length, testability, testability_reason,
    score_retainable, partially_score_censored, minimum_qualifying_score, word_size, effective_search_space,
    evalue_at_minimum_score, expected_chance_bound, production_evalue, informativeness_limit, testable_chance_limit) %>%
  arrange(rule, length_nt), "fig5A_testability_curves.csv")
stopifnot(all(curves$testability %in% names(class_labels)))
spikes <- load_table("short_query_detectability")
database_labels <- c(core_nt_div10 = "core_nt / 10", core_nt = "core_nt", core_nt_x10 = "core_nt \u00d7 10")
stopifnot(all(spikes$database %in% names(database_labels)), all(spikes$permissive_diagnostic_recovered == 1))
spike_summary <- export_csv(spikes %>% group_by(database, length_nt = length) %>%
  summarise(parents = n_distinct(parent_id), planned_spikes = n(), calibrated_recovered = sum(production_calibrated_recovered),
    diagnostic_recovered = sum(permissive_diagnostic_recovered), exact_evalue = unique(calibrated_exact_evalue),
    effective_search_space = unique(effective_search_space), .groups = "drop") %>%
  arrange(match(database, names(database_labels)), length_nt), "fig5A_spike_recovery.csv")
length_breaks <- c(8, 10, 16, 24, 32, 40, 61, 90)
strict_label <- "Strict (exact full length)"
tolerant_label <- "Tolerant (minimum qualifying score)"
# The benchmark is named "benchmark" in visible text, with a multiplication sign for "x 10".
benchmark_series <- "Exact E-value at core_nt / 10 and \u00d7 10 (benchmark)"
missed_marker <- "Benchmark exact spikes of testable lengths missed at E = 10"
testable_label <- "Testable panel length"
y_floor <- 1e-12
y_ceiling <- 1e9
curve_data <- curves %>% mutate(rule = factor(recode(rule, strict = strict_label, tolerant = tolerant_label), levels = c(strict_label, tolerant_label)),
  class = factor(unname(class_labels[testability]), levels = unname(class_labels)),
  below = evalue_at_minimum_score < y_floor, y_plot = pmax(evalue_at_minimum_score, y_floor))
# Points mark the lengths that are testable under the rule of the facet. The curves cover every panel length.
testable_points <- curve_data %>% filter(testability == "testable") %>% mutate(point_label = factor(testable_label))
stopifnot(all(c(strict_label, tolerant_label) %in% testable_points$rule))
curve_lines <- bind_rows(
  curve_data %>% transmute(rule, length_nt, series = "E-value at minimum score (core_nt)", expectation = evalue_at_minimum_score),
  curve_data %>% transmute(rule, length_nt, series = "Chance bound (core_nt)", expectation = expected_chance_bound),
  spike_summary %>% filter(database != "core_nt") %>%
    transmute(rule = factor(strict_label, levels = levels(curve_data$rule)), length_nt, database,
              series = benchmark_series, expectation = exact_evalue)) %>%
  mutate(database = coalesce(database, "core_nt"))
core_spikes <- spike_summary %>% filter(database == "core_nt") %>% inner_join(curves %>% filter(rule == "strict"), by = "length_nt")
stopifnot(nrow(core_spikes) == n_distinct(spike_summary$length_nt), isTRUE(all.equal(core_spikes$exact_evalue, core_spikes$evalue_at_minimum_score)))
strict_testable_lengths <- curves$length_nt[curves$rule == "strict" & curves$testability == "testable"]
all_missed <- spike_summary %>% filter(calibrated_recovered < planned_spikes)
missed <- all_missed %>% filter(length_nt %in% strict_testable_lengths) %>%
  mutate(rule = factor(strict_label, levels = levels(curve_data$rule)), y_plot = pmin(pmax(exact_evalue, y_floor), y_ceiling),
         marker = missed_marker)
production_e <- unique(curves$production_evalue)
stopifnot(all(all_missed$exact_evalue > production_e),
          all(with(spike_summary, calibrated_recovered[exact_evalue <= production_e] == planned_spikes[exact_evalue <= production_e])))
spike_text <- if (nrow(missed)) {
  paste0("Crosses: benchmark exact spikes of testable lengths missed at E = 10 (", serial(sort(unique(missed$length_nt))), " nt).")
} else "All benchmark exact spikes of testable lengths were recovered at E = 10."
size_labels <- spike_summary %>% filter(database != "core_nt") %>% group_by(database) %>%
  filter(exact_evalue <= y_ceiling) %>% slice_min(length_nt, n = 1, with_ties = TRUE) %>% ungroup() %>%
  mutate(rule = factor(strict_label, levels = levels(curve_data$rule)), label = recode(database, core_nt_div10 = "/10", core_nt_x10 = "\u00d710"))
limit_lines <- data.frame(y = c(production_e, unique(curves$informativeness_limit), unique(curves$testable_chance_limit)),
  label = c(paste0("E = ", production_e), paste0("chance ", unique(curves$informativeness_limit)),
            paste0("chance ", unique(curves$testable_chance_limit))),
  colour = c(burgundy, slate, slate), linetype = c("solid", "dotted", "dashed"), vjust = c(-0.35, -0.35, 1.35))
evalue_plot <- ggplot(curve_data, aes(length_nt)) +
  geom_hline(data = limit_lines, aes(yintercept = y, colour = I(colour), linetype = I(linetype)), linewidth = 0.4) +
  geom_line(data = curve_lines, aes(y = expectation, linetype = series, group = interaction(series, database)), colour = "grey50", linewidth = 0.4) +
  geom_point(data = testable_points, aes(y = y_plot, fill = point_label, shape = below), colour = ink, size = 1.6, stroke = 0.35) +
  (if (nrow(missed)) geom_point(data = missed, aes(y = y_plot, size = marker), shape = 4, colour = burgundy, stroke = 0.7)) +
  geom_text(data = size_labels, aes(x = length_nt * 0.93, y = exact_evalue, label = label), hjust = 1, size = 2.4, colour = "grey30") +
  geom_text(data = mutate(limit_lines, rule = factor(tolerant_label, levels = levels(curve_data$rule))),
    aes(x = 125, y = y, label = label, colour = I(colour), vjust = vjust), hjust = 1, size = 2.4) +
  facet_wrap(~rule, nrow = 1) +
  scale_fill_manual(values = setNames(olive, testable_label), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 25), guide = "none") +
  (if (nrow(missed)) scale_size_manual(values = setNames(1.7, missed_marker), name = NULL)) +
  scale_linetype_manual(values = c("E-value at minimum score (core_nt)" = "solid", "Chance bound (core_nt)" = "longdash",
    setNames("dotted", benchmark_series)), name = NULL) +
  scale_x_log10(breaks = length_breaks, limits = c(6, 125), expand = expansion(mult = 0)) +
  scale_y_log10(breaks = 10^seq(-12, 9, 3), labels = scales::label_math(10^.x, format = log10)) +
  coord_cartesian(ylim = c(y_floor, y_ceiling), expand = FALSE, clip = "on") +
  labs(title = "A  A priori testability against query length",
    subtitle = paste0("Accepted core_nt dimensions, E = 10. Points: testable panel lengths at the E-value of the minimum qualifying score.\n",
      "Downward triangles: E-value below the axis. ", spike_text),
    x = "Query length [nt]", y = "Expected chance count") +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, colour = ink)), linetype = guide_legend(order = 2, nrow = 1),
         size = guide_legend(order = 3)) +
  theme(plot.subtitle = element_text(size = 8), legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"),
    legend.key.height = unit(2.6, "mm"), legend.key.width = unit(6, "mm"), legend.text = element_text(size = 7.5),
    legend.margin = margin(0, 0, 0, 0), strip.text = element_text(size = 8), panel.spacing = unit(5, "mm"))

case_headings <- c(terminal = "Terminal match", synthetic_or_vector = "Synthetic or vector record", internal_or_ambiguous = "Internal match")
case_flank_shown <- 25
# Record cases come from the record_cases stage (deterministic selection, 0-3 cases); every case is verified at sequence level.
cases <- load_table("record_cases")
verification <- load_table("record_case_verification")
cases <- cases[!is.na(cases$accession_version), ]
case_columns <- c("panel", "category", "accession_version", "organism", "organism_italic", "topology", "molecule", "record_type", "record_context",
  "record_length_nt", "match_start", "match_end", "match_length_nt", "distance_to_nearest_end_nt", "group", "query_id", "query_label",
  "display_orientation", "query_sequence", "subject_sequence", "display_left_flank", "display_right_flank", "segment_roles")
if (nrow(cases)) {
  stopifnot(nrow(cases) <= 3, identical(verification$accession_version, cases$accession_version), all(cases$category %in% names(case_headings)),
    !anyDuplicated(cases$category), all(verification$query_parent_matches), all(verification$searched_sequence_found),
    all(verification$hsp_full_length_exact), all(verification$record_window_equals_hsp), all(verification$displayed_window_equals_record),
    all(verification$packet_identical | is.na(verification$packet_identical)), all(cases$qseq == cases$sseq),
    all(cases$upper - cases$lower + 1 == nchar(cases$qseq)),
    all(lengths(strsplit(cases$display_segment_roles, ";", fixed = TRUE)) == nchar(cases$qseq)))
  cases <- cases[order(match(cases$category, names(case_headings))), ]
  verification <- verification[match(cases$accession_version, verification$accession_version), ]
  cases$display_left_flank[is.na(cases$display_left_flank)] <- ""
  cases$display_right_flank[is.na(cases$display_right_flank)] <- ""
  case_data <- export_csv(data.frame(
    panel = LETTERS[1 + seq_len(nrow(cases))], category = cases$category, accession_version = cases$accession_version,
    organism = cases$organism, organism_italic = verification$organism_italic, topology = cases$topology, molecule = cases$molecule,
    record_type = cases$record_type, record_context = cases$context, record_length_nt = cases$record_length, match_start = cases$lower,
    match_end = cases$upper, match_length_nt = nchar(cases$qseq),
    distance_to_nearest_end_nt = ifelse(cases$topology == "circular", NA, pmin(cases$lower - 1, cases$record_length - cases$upper)),
    group = cases$group, query_id = cases$query_id, query_label = parent_label(cases$query_id),
    display_orientation = cases$display_orientation, query_sequence = cases$qseq, subject_sequence = cases$sseq,
    display_left_flank = cases$display_left_flank, display_right_flank = cases$display_right_flank,
    segment_roles = cases$display_segment_roles), "fig5BD_record_cases.csv")
  stopifnot(all(case_data$distance_to_nearest_end_nt[case_data$topology == "linear"] == cases$distance_to_end[cases$topology == "linear"]))
  case_data$display_left_flank[is.na(case_data$display_left_flank)] <- ""
  case_data$display_right_flank[is.na(case_data$display_right_flank)] <- ""
} else {
  case_data <- export_csv(setNames(data.frame(matrix(character(), nrow = 0, ncol = length(case_columns))), case_columns), "fig5BD_record_cases.csv")
}
invisible(export_csv(data.frame(
  panel = case_data$panel, accession_version = case_data$accession_version, organism = case_data$organism, topology = case_data$topology,
  record_length_nt = case_data$record_length_nt, match_start = case_data$match_start, match_end = case_data$match_end,
  distance_to_nearest_end_nt = case_data$distance_to_nearest_end_nt, query = case_data$query_label,
  match_length_nt = case_data$match_length_nt, case_category = unname(case_headings[case_data$category])), "tableS8_record_case_summary.csv"))

track_x0 <- 16
track_x1 <- 170
alignment_x0 <- 14
alignment_dx <- 1.75
minimum_marker_mm <- 0.9
record_case_panel <- function(row) {
  circular <- row$topology == "circular"
  marker_colour <- unname(colours[row$group])
  stopifnot(!is.na(marker_colour))
  length_nt <- row$record_length_nt
  position_x <- function(position) track_x0 + (track_x1 - track_x0) * position / length_nt
  type_parts <- c(row$topology, if (row$record_type == "genomic") "genomic" else if (row$record_context == "vector") "vector", row$molecule)
  type_label <- paste(type_parts, collapse = " ")
  match_range <- paste0(format(row$match_start, big.mark = ",", trim = TRUE), "\u2013", format(row$match_end, big.mark = ",", trim = TRUE))
  length_label <- format(length_nt, big.mark = ",", trim = TRUE)
  layers <- list(
    mm_text(0.5, 3.6, paste0(row$panel, "  ", case_headings[[row$category]], ": ", row$accession_version), 10, "bold"),
    mm_text_species(0.5, 8.4, row$organism, row$organism_italic, paste0(", ", type_label, ". ", row$query_label, ", ", row$match_length_nt, "-nt exact match, ",
      ifelse(row$display_orientation == "-", "minus", "plus"), " orientation"), 8.3))
  if (circular) {
    centre_x <- 27; centre_y <- 26.5; radius <- 9.5
    angle <- seq(0, 2 * pi, length.out = 361)
    match_angle <- 2 * pi * ((row$match_start - 1 + row$match_end) / 2) / length_nt
    marker <- c(sin(match_angle), -cos(match_angle))
    # The ring starts at the top (position 1, dark tick). The match range is written in the text column, so it never collides with the ring.
    origin_text <- if (row$match_start == 1) " at the stored origin" else ""
    layers <- c(layers, list(
      annotate("path", x = centre_x + radius * sin(angle), y = centre_y - radius * cos(angle), colour = "grey70", linewidth = 1.6),
      mm_segment(centre_x, centre_y - radius + 1.3, centre_x, centre_y - radius - 1.3, colour = ink, linewidth = 0.5),
      mm_segment(centre_x + (radius - 2.6) * marker[1], centre_y + (radius - 2.6) * marker[2],
                 centre_x + (radius + 2.6) * marker[1], centre_y + (radius + 2.6) * marker[2], colour = marker_colour, linewidth = 1.8),
      mm_text(centre_x, centre_y, length_label, 8, hjust = 0.5),
      mm_text(centre_x + radius + 9, centre_y - 5, paste0("Match ", match_range, origin_text, "."), 8),
      mm_text(centre_x + radius + 9, centre_y, paste0("Circular record of ", length_label, " nt, drawn to scale."), 8),
      mm_text(centre_x + radius + 9, centre_y + 5, "Position 1 (top tick) is the stored origin of the circular record.", 8)))
    alignment_top <- 42.5
    height <- 51
  } else {
    marker_centre <- position_x((row$match_start - 1 + row$match_end) / 2)
    marker_width <- max(position_x(row$match_end) - position_x(row$match_start - 1), minimum_marker_mm)
    nearest_left <- (row$match_start - 1) <= (length_nt - row$match_end)
    end_x <- if (nearest_left) track_x0 else track_x1
    edge_x <- if (nearest_left) marker_centre - marker_width / 2 else marker_centre + marker_width / 2
    label_side <- if (marker_centre > 120) 1 else if (marker_centre < 40) 0 else 0.5
    below_resolution <- abs(edge_x - end_x) < 8
    layers <- c(layers, list(
      mm_rect(track_x0, track_x1, 16.2, 19.8, fill = "#EEEEEE", colour = "grey55", linewidth = 0.35),
      mm_rect(marker_centre - marker_width / 2, marker_centre + marker_width / 2, 14.6, 21.4, fill = marker_colour),
      mm_text(marker_centre + ifelse(label_side == 1, marker_width / 2, ifelse(label_side == 0, -marker_width / 2, 0)),
        12.6, paste0("Match ", match_range), 7.5, hjust = label_side),
      # End labels keep clear of the tick that marks a match end below the scale resolution.
      if (below_resolution && nearest_left) mm_text(track_x0 - 1.4, 23.4, "1", 7, hjust = 1) else mm_text(track_x0, 23.4, "1", 7, hjust = 0.5),
      mm_text(if (below_resolution && !nearest_left) min(track_x1 - 1.2, edge_x - 1.6) else track_x1 - 1.2, 23.4, length_label, 7, hjust = 1)))
    distance_label <- paste0(format(row$distance_to_nearest_end_nt, big.mark = ",", trim = TRUE), " nt to the nearest record end")
    if (!below_resolution) {
      layers <- c(layers, list(
        annotate("segment", x = end_x, xend = edge_x, y = 26.6, yend = 26.6, colour = ink, linewidth = 0.4,
          arrow = arrow(length = unit(1.3, "mm"), ends = "both", type = "closed", angle = 22)),
        mm_text((end_x + edge_x) / 2, 29.6, distance_label, 7.5, hjust = 0.5)))
    } else {
      layers <- c(layers, list(
        mm_segment(edge_x, 21.4, edge_x, 26, colour = ink, linewidth = 0.4),
        mm_text(edge_x, 28.4, paste0(distance_label, " (below the scale resolution)"), 7.5, hjust = if (nearest_left) 0 else 1)))
    }
    alignment_top <- 35
    height <- 43
  }
  left_flank <- row$display_left_flank
  right_flank <- row$display_right_flank
  # Up to 25 nt of flank per side; flanks shrink (and, for very long queries, the base pitch) so the alignment fits the width.
  available <- floor((figure_width_mm - alignment_x0 - 2) / alignment_dx)
  flank_shown <- max(0, min(case_flank_shown, floor((available - nchar(row$query_sequence)) / 2)))
  left_shown <- if (flank_shown) substring(left_flank, max(1, nchar(left_flank) - flank_shown + 1)) else ""
  right_shown <- if (flank_shown) substring(right_flank, 1, flank_shown) else ""
  subject <- strsplit(paste0(left_shown, row$subject_sequence, right_shown), "")[[1]]
  query <- c(rep(NA, nchar(left_shown)), strsplit(row$query_sequence, "")[[1]], rep(NA, nchar(right_shown)))
  roles <- c(rep("flank", nchar(left_shown)), strsplit(row$segment_roles, ";", fixed = TRUE)[[1]], rep("flank", nchar(right_shown)))
  roles[!roles %in% names(role_colours)] <- "unknown"
  dx <- min(alignment_dx, (figure_width_mm - alignment_x0 - 2) / length(subject))
  positions <- seq_along(subject)
  x <- alignment_x0 + (positions - 0.5) * dx
  in_match <- roles != "flank"
  runs <- rle(roles)
  run_end <- cumsum(runs$lengths)
  run_start <- run_end - runs$lengths + 1
  match_runs <- which(runs$values != "flank")
  layers <- c(layers, list(
    mm_rect(alignment_x0 + (run_start[match_runs] - 1) * dx, alignment_x0 + run_end[match_runs] * dx,
            alignment_top - 2.1, alignment_top + 6.6, fill = unname(role_fills[runs$values[match_runs]])),
    mm_text(6, alignment_top, "S", 8, "bold", family = "mono"), mm_text(6, alignment_top + 4.5, "Q", 8, "bold", family = "mono"),
    mm_text(x, alignment_top, subject, 7.5, hjust = 0.5, family = "mono", colour = unname(role_colours[roles])),
    mm_text(x[in_match], alignment_top + 4.5, query[in_match], 7.5, hjust = 0.5, family = "mono", colour = unname(role_colours[roles[in_match]]))))
  canvas(height, layers)
}
case_panels <- lapply(seq_len(nrow(case_data)), function(index) record_case_panel(case_data[index, ]))
case_heights <- ifelse(case_data$topology == "circular", 51, 43)
if (!length(case_panels)) {
  case_panels <- list(canvas(14, list(panel_title("B  Record cases"),
    mm_text(0.5, 10, "No record case was selected. The record-case selection tables list the category status.", 8.5))))
  case_heights <- 14
}
figure4_legend_notes <- c(
  "Record cases are deterministically selected by the record_cases stage and verified at sequence level for each case.",
  paste0("Positions are scaled within each record (separate scales). Markers in the platform color of the query, at least ", minimum_marker_mm, " mm wide."),
  paste0("Q and S: aligned query and subject. Slate: scaffold. Burgundy: barcode or index. Up to ", case_flank_shown,
         " nt of gray flank, clipped at record boundaries. Details: Figure S2."))
# Panel A keeps its reference height unit; each record diagram is 43 mm (linear) or 51 mm (circular).
evalue_height <- 74.6
figure4_heights <- c(evalue_height, as.vector(rbind(separator_height, case_heights)))
figure4 <- Reduce(`/`, c(list(independent_panel(with_margin(evalue_plot, bottom = 0))),
  unlist(lapply(case_panels, function(panel) list(panel_separator(), panel)), recursive = FALSE))) + plot_layout(heights = figure4_heights)
save_figure("figure4", figure4, stack_height(figure4_heights, 9),
  c("A priori testability and exact-spike recovery against query length", paste("Record diagram:", unname(case_headings[case_data$category]))),
  c("length_testability.tsv", "short_query_detectability.tsv", "record_cases.tsv", "record_case_verification.tsv"),
  "Testability classes from the production calibration. Benchmark exact spikes at three database sizes. Record cases from the record_cases stage",
  c("fig5A_testability_curves.csv", "fig5A_spike_recovery.csv", "fig5BD_record_cases.csv", "tableS8_record_case_summary.csv"),
  legend_notes = figure4_legend_notes)

# ---------------------------------------------------------------------------------------------
# Figure S1. Query properties and matching diagnostics
# ---------------------------------------------------------------------------------------------
condition_names <- c(
  ont_barcode = "ONT: Native barcode", ont_native_barcode_context_k40 = "ONT: K40 context",
  ont_native_barcode_context_k46 = "ONT: K46 context", ont_rapid_barcode_context_90 = "ONT: Rapid barcode context",
  ont_ligation_adapter = "ONT: Ligation adapter", ont_rapid_adapter = "ONT: Rapid adapter", ont_pcr_adapter = "ONT: PCR adapter",
  ont_cdna_adapter = "ONT: cDNA adapter",
  illumina_i7_index = "Illumina: bare i7", illumina_i5_index = "Illumina: bare i5",
  illumina_i7_index_adapter = "Illumina: embedded i7", illumina_i5_index_adapter = "Illumina: embedded i5",
  illumina_truseq = "Illumina: TruSeq", illumina_nextera = "Illumina: Nextera", illumina_smallrna = "Illumina: Small RNA",
  illumina_bgi = "BGI/DNBSEQ: adapters", pacbio_barcode = "PacBio: barcode", pacbio_barcode_adapter = "PacBio: barcode-hairpin",
  pacbio_smrtbell_adapter = "PacBio: SMRTbell / C2", pacbio_c2_primer = "PacBio: C2 primer")
# Panel A shows the strict-testable parents only.
properties <- load_table("query_properties") %>% filter(strict_testable) %>% group_order()
stopifnot(nrow(properties) > 0)
# Conditions without a declared label are shown under their condition name, after the declared ones.
extra_conditions <- sort(setdiff(unique(properties$condition), names(condition_names)))
condition_names <- c(condition_names, setNames(extra_conditions, extra_conditions))
properties <- properties %>%
  mutate(condition_label = factor(unname(condition_names[condition]), levels = rev(unname(condition_names[names(condition_names) %in% condition]))))
length_counts <- properties %>% count(group, condition_label, length, name = "parents")
composition_counts <- properties %>% count(group, condition_label, gc_fraction, name = "parents")
property_count_limit <- max(c(length_counts$parents, composition_counts$parents))
property_size <- scale_size_continuous(limits = c(1, property_count_limit), range = c(1, 3),
  breaks = c(1, 100, 200), name = "Parents at value")
first_testable <- function(rule) min(testability$query_length[testability$endpoint == rule & testability$testability == "testable"])
strict_first <- first_testable("strict")
tolerant_first <- first_testable("tolerant")
stopifnot(all(testability$testability[testability$endpoint == "strict" & testability$query_length < strict_first] != "testable"),
          all(testability$testability[testability$endpoint == "tolerant" & testability$query_length < tolerant_first] != "testable"))
stopifnot(min(properties$length) >= strict_first)
# The dotted line marks the first tolerant-testable panel length when strict-testable parents lie below it.
tolerant_line <- any(properties$length < tolerant_first)
length_plot <- ggplot(length_counts, aes(length, condition_label, colour = group)) +
  (if (tolerant_line) geom_vline(xintercept = tolerant_first - 0.5, colour = olive, linetype = "dotted", linewidth = 0.45)) +
  geom_point(aes(size = parents), alpha = 0.6) + property_size +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_continuous(limits = c(0, 120), breaks = sort(unique(c(0, strict_first, if (tolerant_line) tolerant_first, 60, 90, 120)))) +
  labs(title = "A  Query properties of strict-testable parents",
    x = paste0("Length [nt]. Strict testable from ", strict_first, " nt.",
               if (tolerant_line) paste0("\nDotted: tolerant testable from ", tolerant_first, " nt (first testable panel length).") else ""), y = NULL) +
  theme(axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 7.5))
composition_plot <- ggplot(composition_counts, aes(gc_fraction * 100, condition_label, colour = group)) +
  geom_point(aes(size = parents), alpha = 0.6) + property_size +
  scale_colour_manual(values = colours, guide = "none") + scale_x_continuous(limits = c(0, 100)) +
  labs(title = "", x = "GC [%]", y = NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
properties_panel <- ((length_plot | composition_plot) + plot_layout(widths = c(1.5, 0.75), guides = "collect")) &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.key.height = unit(3, "mm"),
    plot.margin = margin(5, 5, 0, 5))

# Independent benchmark. The by-length pair recall is reported in the text and exported for Table S9 (no panel).
benchmark_all <- load_table("benchmark_metrics")
benchmark <- benchmark_all %>% filter(database == "core_nt", population == "all_lengths") %>%
  mutate(arm_label = recode(search_space, local_fixture = "Local fixture", production_calibrated = "Calibrated E=10",
    permissive_diagnostic = "Diagnostic E=1e12")) %>%
  arrange(match(search_space, c("local_fixture", "production_calibrated", "permissive_diagnostic")), endpoint)
stopifnot(nrow(benchmark) == 6, all(benchmark_all$false_positive == 0),
  all(benchmark$positive_pairs + benchmark$negative_pairs == manifest$benchmarks$core_nt$pair_universe))
negative_counts <- benchmark %>% distinct(endpoint, negative_pairs)
stopifnot(nrow(negative_counts) == 2)
negative_note <- paste0("No truth-negative pairs accepted: ", paste0(negative_counts$endpoint, " 0/",
  scales::comma(negative_counts$negative_pairs, accuracy = 1), collapse = ", "), " (each arm).")
benchmark_length <- benchmark_all %>%
  filter(population == "length", search_space %in% c("production_calibrated", "permissive_diagnostic")) %>%
  mutate(arm_label = recode(search_space, production_calibrated = "E = 10", permissive_diagnostic = "E = 1e12"))
stopifnot(!anyNA(benchmark_length$recall), all(benchmark_length$false_positive == 0))
pair_recall <- export_csv(benchmark_length %>%
  arrange(match(database, c("core_nt", "core_nt_div10", "core_nt_x10")), match(search_space, c("production_calibrated", "permissive_diagnostic")), endpoint, length) %>%
  transmute(database, arm = search_space, arm_label, rule = endpoint, length_nt = length, truth_positive_pairs = positive_pairs,
    recovered_pairs = true_positive, missed_pairs = false_negative, pair_recall = recall), "tableS9_pair_recall_by_length.csv")
# At production scale an exact spike is recovered exactly when its exact E-value passes the production threshold.
strict_calibrated <- filter(pair_recall, database == "core_nt", arm == "production_calibrated", rule == "strict") %>%
  inner_join(testability %>% filter(endpoint == "strict") %>% select(length_nt = query_length, evalue_at_minimum_score), by = "length_nt")
stopifnot(all(strict_calibrated$pair_recall[strict_calibrated$evalue_at_minimum_score > production_evalue] == 0),
  all(strict_calibrated$pair_recall[strict_calibrated$evalue_at_minimum_score <= production_evalue] == 1),
  all(pair_recall$pair_recall[pair_recall$arm == "permissive_diagnostic" & pair_recall$rule == "strict"] == 1))
# Returned HSPs of the strict-testable reference sequences (figure_sources.alignment_quality).
quality <- load_table("quality_marginal") %>% group_order()
quality_hsps <- manifest$counts$strict_testable_reference_hsps
stopifnot(sum(load_table("quality_acceptance")$hsp_rows) >= quality_hsps)
quality_guides <- data.frame(measure = c("Identity", "Query span"), threshold = c(95, 90))
quality_plot <- ggplot(quality, aes(value, cumulative, colour = group)) +
  geom_step(linewidth = 0.5) + geom_point(data = quality %>% group_by(group, measure) %>% filter(n() == 1) %>% ungroup(), size = 1.5) +
  geom_vline(data = quality_guides, aes(xintercept = threshold), inherit.aes = FALSE,
    linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  facet_wrap(~measure, nrow = 1, scales = "free_x") + scale_colour_manual(values = colours, name = NULL) +
  scale_x_continuous(limits = function(limits) c(if (min(limits) >= 80) 80 else 0, 100),
    breaks = function(limits) if (min(limits) >= 75) seq(80, 100, 5) else seq(0, 100, 25)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title = "B  Cumulative returned-alignment distributions", subtitle = paste0("All ", scales::comma(quality_hsps), " returned HSPs of strict-testable reference queries. Rejected HSPs are included.\n",
        "Guides: 95% identity and 90% query span."),
      x = "Identity / query span [%]", y = "HSP fraction") +
    theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0), legend.key.height = unit(3, "mm"))
joint <- load_table("quality_joint") %>% group_order() %>%
  group_by(group, identity_bin, span_bin) %>% summarise(hsp_rows = sum(hsp_rows), .groups = "drop")
joint_identity_start <- min(80, floor(min(joint$identity_bin) / 10) * 10)
joint_plot <- ggplot(joint, aes(identity_bin + 0.5, span_bin + 0.5, fill = group, alpha = hsp_rows)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = 95, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 90, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~group, nrow = 1) + scale_fill_manual(values = colours, guide = "none") +
  scale_alpha_continuous(trans = "log10", range = c(0.4, 1), name = "HSPs per bin (color intensity)",
    breaks = c(1, 100, 10000, 100000), labels = short_count,
    guide = guide_legend(override.aes = list(fill = ink), nrow = 1)) +
  scale_x_continuous(breaks = seq(joint_identity_start, 100, 5) + 0.5, labels = seq(joint_identity_start, 100, 5)) +
  scale_y_continuous(breaks = seq(0, 100, 25) + 0.5, labels = seq(0, 100, 25)) +
  coord_cartesian(xlim = c(joint_identity_start, 101), ylim = c(0, 101), expand = FALSE) +
  labs(title = "C  Joint identity and span", subtitle = "HSPs of strict-testable reference queries in 1%-wide bins labeled by lower bound. 100 denotes exact 100%.\nColor intensity: HSPs per bin (log10). Dashed lines: 95% identity and 90% query span.",
      x = "Identity bin [%]", y = "Span bin [%]")
figureS1_heights <- c(90, separator_height, 60, separator_height, 84)
figure_s1 <- (wrap_elements(full = patchworkGrob(properties_panel)) / panel_separator() /
    independent_panel(with_margin(quality_plot, top = 2, bottom = 0)) / panel_separator() /
    independent_panel(with_margin(joint_plot, top = 2))) + plot_layout(heights = figureS1_heights)
save_figure("figureS1", figure_s1, sum(figureS1_heights), c("Query properties", "Marginal quality", "Joint quality"),
  c("query_properties.tsv", "benchmark_metrics.tsv", "quality_marginal.tsv", "quality_joint.tsv", "quality_acceptance.tsv"),
  "Independent current-version fixture QA",
  "tableS9_pair_recall_by_length.csv")

# ---------------------------------------------------------------------------------------------
# Target-limit comparison (500 versus 2,000 targets and production) and cap-reaching searches. No panel: the values are
# exported for Table S10.
# ---------------------------------------------------------------------------------------------
parent_names <- load_table("logical_queries") %>% filter(role == "real") %>% distinct(parent_id, parent_logical_id = logical_id)
parent_display <- parent_names %>% group_by(parent_id) %>%
  summarise(reference_parent = paste(parent_label(sort(parent_logical_id)), collapse = " / "), .groups = "drop")
target_summary <- load_table("target_limit_summary")
if (nrow(target_summary)) target_summary <- target_summary %>% arrange(match(group, c(groups, "Global union")), endpoint, population, gate)
invisible(export_csv(target_summary, "tableS10_target_limit_summary.csv"))
cap_membership <- load_table("cap_memberships")
caps <- load_table("cap_sequences") %>%
  mutate(display_label = ifelse(is.na(aliases), paste0("Draw ", substr(sequence_id, 5, 12)), parent_label(aliases)),
    role_label = recode(roles, real = "Reference", control = "Control", .default = "Reference and control")) %>%
  left_join(select(cap_membership, sequence_id, membership_parent_id = parent_id, draw_index, membership_role = role), by = "sequence_id") %>%
  left_join(parent_display, by = c("membership_parent_id" = "parent_id")) %>%
  mutate(row_label = paste0(groups, " ", reference_parent, ", ", ifelse(membership_role == "real", "reference", paste0("draw ", draw_index)),
                            ", ", length, " nt"))
stopifnot(!anyNA(caps$reference_parent), all(caps$targets >= target_limit))
invisible(export_csv(data.frame(
  display_label = caps$display_label, group = caps$groups, role = caps$role_label, length_nt = caps$length,
  raw_hsp_rows = caps$hsp_rows, accepted_strict_hsps = caps$strict_retained, accepted_tolerant_hsps = caps$tolerant_retained,
  search_sequence_id = caps$sequence_id, full_aliases = caps$aliases, logical_memberships = caps$logical_memberships,
  reference_parent = caps$reference_parent, control_draw_index = ifelse(caps$membership_role == "control", caps$draw_index, NA),
  search_label = caps$row_label, target_cap = rep(target_limit, nrow(caps))), "tableS10_cap_reaching_searches.csv"))
target_parents <- load_table("target_limit_parents")
if (nrow(target_parents)) {
  target_parents <- target_parents %>% left_join(parent_display, by = "parent_id") %>% relocate(reference_parent) %>%
    arrange(match(sub(";.*$", "", groups), groups), endpoint, population, reference_parent)
  stopifnot(!anyNA(target_parents$reference_parent))
}
invisible(export_csv(target_parents, "tableS10_target_limit_parents.csv"))

# ---------------------------------------------------------------------------------------------
# Figure S2. Control quality (Table S6) and Monte Carlo precision
# ---------------------------------------------------------------------------------------------
qc <- load_table("control_qc_summary")
qc_parents <- load_table("control_parent_qc") %>% group_order() %>%
  mutate(testable_identical = parent_identical_canonical * strict_testable, testable_repeated = repeated_canonical * strict_testable)
draw_qc <- load_table("control_draw_qc") %>% select(parent_id, invariants_passed) %>%
  left_join(select(qc_parents, parent_id, group), by = "parent_id")
stopifnot(!anyNA(draw_qc$group))
quality_by_group <- bind_rows(lapply(groups[groups %in% qc_parents$group], function(name) {
  parents <- filter(qc_parents, group == name)
  draw_rows <- filter(draw_qc, group == name)
  cbind(group = name, summarise(parents, parents = n(), distinct_sequences_median = median(distinct_canonical),
    distinct_sequences_min = min(distinct_canonical), distinct_sequences_max = max(distinct_canonical),
    parent_identical_canonical = sum(parent_identical_canonical), repeated_canonical = sum(repeated_canonical),
    parent_identical_exact = sum(parent_identical_exact), repeated_exact = sum(repeated_exact),
    draws_passing_invariants = sum(draw_rows$invariants_passed), draws_total = nrow(draw_rows),
    testable_parents = sum(strict_testable), testable_parent_identical_canonical = sum(testable_identical),
    testable_repeated_canonical = sum(testable_repeated)))
}))
quality_total <- cbind(group = "All groups", summarise(qc_parents, parents = n(), distinct_sequences_median = median(distinct_canonical),
  distinct_sequences_min = min(distinct_canonical), distinct_sequences_max = max(distinct_canonical),
  parent_identical_canonical = sum(parent_identical_canonical), repeated_canonical = sum(repeated_canonical),
  parent_identical_exact = sum(parent_identical_exact), repeated_exact = sum(repeated_exact),
  draws_passing_invariants = sum(draw_qc$invariants_passed), draws_total = nrow(draw_qc),
  testable_parents = sum(strict_testable), testable_parent_identical_canonical = sum(testable_identical),
  testable_repeated_canonical = sum(testable_repeated)))
control_quality <- export_csv(rbind(quality_by_group, quality_total), "tableS6_control_quality.csv")
total_row <- nrow(control_quality)
stopifnot(control_quality$group[total_row] == "All groups",
  control_quality$parent_identical_canonical[total_row] == qc$parent_identical_canonical, control_quality$repeated_canonical[total_row] == qc$repeated_canonical,
  control_quality$parent_identical_exact[total_row] == qc$parent_identical_exact, control_quality$repeated_exact[total_row] == qc$repeated_exact,
  control_quality$draws_passing_invariants[total_row] == qc$selected_draws, control_quality$draws_total[total_row] == qc$selected_draws,
  control_quality$testable_parents[total_row] == qc$testable_parents,
  control_quality$testable_parent_identical_canonical[total_row] == qc$testable_parent_identical_canonical,
  control_quality$testable_repeated_canonical[total_row] == qc$testable_repeated_canonical)

# Control draws of the strict-testable parents (legend note; all parents: Table S6).
testable_control_note <- paste0(scales::comma(qc$testable_parents), " strict-testable parents have ", scales::comma(qc$testable_draws),
  " control draws. ", if (qc$testable_parent_identical_canonical == 0 && qc$testable_repeated_canonical == 0)
    "None of these draws repeats within its parent or equals its parent." else
    paste0(scales::comma(qc$testable_parent_identical_canonical), " of these draws equal their parent and ",
           scales::comma(qc$testable_repeated_canonical), " repeat within their parent."),
  " Reverse complements count as identical.")

precision <- load_table("precision")
precision_selected <- precision %>% filter(condition == "all_conditions", population == "all_records", group != "Global union") %>%
  group_order() %>% mutate(endpoint = factor(endpoint, levels = c("strict", "tolerant")),
    measure = factor(measure, levels = c("detection", "breadth")),
    group_measure = paste(as.character(group), ifelse(measure == "detection", "Detection [%]", "log10(breadth + 1)"), sep = "\n"))
precision_levels <- as.vector(outer(groups, c("Detection [%]", "log10(breadth + 1)"), paste, sep = "\n"))
precision_selected <- precision_selected %>% mutate(group_measure = factor(group_measure, levels = precision_levels))
precision_selected <- precision_selected %>% mutate(
  display_estimate = ifelse(measure == "detection", estimate * 100, log10(estimate + 1)),
  display_lower = ifelse(measure == "detection", lower * 100, log10(lower + 1)),
  display_upper = ifelse(measure == "detection", upper * 100, log10(upper + 1)))
fixed_precision <- load_table("precision_fixed_cohort") %>% filter(population == "all_records", group != "Global union") %>%
  mutate(group_measure = factor(paste(group, ifelse(measure == "detection", "Detection [%]", "log10(breadth + 1)"), sep = "\n"), levels = precision_levels),
    display_estimate = ifelse(measure == "detection", estimate * 100, log10(estimate + 1)))
precision_selected <- precision_selected %>% mutate(
  marker_fill = ifelse(status == "conditional_bound", unname(colours[as.character(group)]), "white"),
  rule_label = factor(recode(as.character(endpoint), strict = "Strict", tolerant = "Tolerant"), levels = c("Strict", "Tolerant")))
fixed_precision <- fixed_precision %>% mutate(group = factor(group, levels = groups),
  rule_label = factor(recode(endpoint, strict = "Strict", tolerant = "Tolerant"), levels = c("Strict", "Tolerant")))
precision_plot <- ggplot(precision_selected, aes(draws_per_parent, display_estimate, colour = group, shape = rule_label, group = rule_label)) +
  geom_linerange(aes(ymin = display_lower, ymax = display_upper), linewidth = 0.45, position = position_dodge(width = 1.4)) +
  geom_point(aes(fill = I(marker_fill)), position = position_dodge(width = 1.4), size = 1.7, stroke = 0.6) +
  geom_line(data = fixed_precision, aes(draws_per_parent, display_estimate, colour = group, linetype = rule_label, group = rule_label),
    inherit.aes = FALSE, linewidth = 0.4) +
  facet_wrap(~group_measure, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = draw_prefixes) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_shape_manual(values = c(Strict = 21, Tolerant = 24), name = NULL) +
  scale_linetype_manual(values = c(Strict = "dotted", Tolerant = "dashed"), name = NULL) +
  guides(shape = guide_legend(override.aes = list(colour = ink, fill = ink)), linetype = guide_legend()) +
  labs(title = "A  Fixed-prefix estimates and Monte Carlo intervals",
      subtitle = "Filled markers: conditional bound. Open markers: resolved-subset diagnostic. Y scales differ between panels.\nLines: final-parent cohort at each draw prefix (dotted: strict, dashed: tolerant). Bars: interval limits.",
       x = "Draws per parent", y = "Control mean\nand interval") +
  theme(legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))
precision_final <- precision %>% filter(draws_per_parent == final_prefix, condition == "all_conditions") %>%
  mutate(group = factor(group, levels = rev(c(groups, "Global union"))),
    endpoint = factor(endpoint, levels = c("strict", "tolerant")),
    population = recode(population, all_records = "All records", linear_genomic = "Linear-genomic"),
    cell = factor(paste(endpoint, measure, sep = " / "), levels = c("strict / detection", "strict / breadth", "tolerant / detection", "tolerant / breadth")),
    status_label = case_when(precision_met %in% "yes" ~ "Met", precision_met %in% "no" ~ "Missed",
      status == "diagnostic_resolved_subset_no_coverage_claim" ~ "Subset", TRUE ~ "No cohort"),
    ratio = half_width / target_half_width)
precision_counts <- precision %>% filter(draws_per_parent == final_prefix)
precision_count_label <- function(data) {
  counts <- c(met = sum(data$precision_met %in% "yes"), missed = sum(data$precision_met %in% "no"),
    subset = sum(is.na(data$precision_met) & data$status == "diagnostic_resolved_subset_no_coverage_claim"),
    empty = sum(is.na(data$precision_met) & data$status == "no_resolved_matched_population"))
  stopifnot(sum(counts) == nrow(data))
  parts <- c(paste(counts[["met"]], "targets met"), paste(counts[["missed"]], "missed"), paste(counts[["subset"]], "resolved-subset diagnostics"),
    if (counts[["empty"]]) paste(counts[["empty"]], "without a resolved cohort"))
  paste0(nrow(data), " rows: ", serial(parts), ".")
}
precision_tiles <- export_csv(precision_final %>% transmute(population, group = as.character(group), rule = as.character(endpoint), measure,
  half_width, target_half_width, half_width_over_target = ratio, status = status_label) %>%
  arrange(population, match(group, c(groups, "Global union")), rule, measure), "figS4B_precision_tiles.csv")
precision_final <- precision_final %>% mutate(text_colour = ifelse(log10(ratio) > 1.2, "white", "black"),
  tile_label = paste0(status_label, "\n", vapply(signif(ratio, 3), format, character(1), scientific = FALSE, trim = TRUE)))
precision_status <- ggplot(precision_final, aes(cell, group, fill = ratio)) +
  geom_tile(colour = "white", linewidth = 0.7) + geom_text(aes(label = tile_label, colour = I(text_colour)), size = 2.9, lineheight = 0.95) +
  facet_wrap(~population, nrow = 1) +
  scale_fill_gradientn(colours = sand_burgundy_ramp, trans = "log10", name = "Half-width / target",
    breaks = c(1, 3, 10, 30, 100, 300), labels = c("1", "3", "10", "30", "100", "300"),
    guide = guide_colourbar(barwidth = unit(55, "mm"), barheight = unit(3, "mm"))) +
  scale_x_discrete(labels = c("strict / breadth" = "Strict\nbreadth", "strict / detection" = "Strict\ndetect.",
    "tolerant / breadth" = "Tolerant\nbreadth", "tolerant / detection" = "Tolerant\ndetect.")) +
  labs(title = "B  Final precision targets",
    subtitle = paste0(precision_count_label(filter(precision_counts, condition == "all_conditions")), " Met: half-width \u2264 target.",
      "\nTile text: status and half-width / target. Fill: the same ratio on a log10 scale. Values: Table S6."), x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
figureS2_heights <- c(88, separator_height, 88)
figure_s3 <- (independent_panel(with_margin(precision_plot, bottom = 0)) / panel_separator() /
    independent_panel(with_margin(precision_status, top = 2))) + plot_layout(heights = figureS2_heights)
save_figure("figureS2", figure_s3, sum(figureS2_heights), c("Prefix precision", "Precision targets"),
  c("control_draw_qc.tsv", "control_parent_qc.tsv", "control_qc_summary.tsv", "precision.tsv", "precision_fixed_cohort.tsv"),
  paste0("All final-stage ", precision_count_label(precision_counts), " Intervals describe conditional Monte Carlo error."),
  c("tableS6_control_quality.csv", "figS4B_precision_tiles.csv"), legend_notes = testable_control_note)

# ---------------------------------------------------------------------------------------------
# Figure S3. Positional robustness
# ---------------------------------------------------------------------------------------------
dense_position <- filter(position, group %in% c("ONT", "Illumina"))
sparse_position <- filter(position, group %in% c("BGI/DNBSEQ", "PacBio"))
all_position_totals <- position %>% group_by(group) %>% summarise(events = n(), accessions = n_distinct(accession_version), .groups = "drop")
all_position_labels <- setNames(paste0(all_position_totals$group, "\n", scales::comma(all_position_totals$events),
  " events\n", scales::comma(all_position_totals$accessions), ifelse(all_position_totals$accessions == 1, " accession", " accessions")),
  as.character(all_position_totals$group))
position_scatter <- ggplot(position, aes(log10(distance + 1), log10(record_length))) +
  geom_hex(data = dense_position, aes(fill = group, alpha = after_stat(count)), bins = 22) +
  geom_point(data = sparse_position, aes(colour = group), size = 1.8, alpha = 0.8) +
  geom_vline(xintercept = log10(main_window + 1), linetype = "dashed", linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = log10(main_minimum), linetype = "dashed", linewidth = 0.3, colour = "grey55") +
  facet_wrap(~group, nrow = 1, labeller = labeller(group = all_position_labels)) +
  scale_fill_manual(values = colours, guide = "none") +
  scale_alpha_continuous(trans = "log10", range = c(0.15, 1), breaks = c(1, 10, 100, 1000), labels = short_count,
    name = "Events / bin\n(intensity)", guide = guide_legend(override.aes = list(fill = ink))) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_continuous(breaks = log10(c(0, 1000, 1000000) + 1), labels = c("0", "1k", "1M")) +
  scale_y_continuous(breaks = 2:7, labels = c("100", "1k", "10k", "100k", "1M", "10M")) +
    labs(title = "A  Position and record length", subtitle = paste0("Strict matches of testable queries in eligible linear-genomic records of all lengths.\nHexagons: events per bin (log10). Points: sparse groups. Dashed: records \u2265", main_minimum / 1000, " kb, distance \u2264", main_window, " nt."),
      x = "Distance [nt, log10(distance + 1)]", y = "Length [nt, log10]") +
    theme(strip.text = element_text(size = 8), legend.position = "right")
density <- load_table("density") %>% group_order()
density_summary <- density %>% group_by(group) %>% summarise(accessions = n(), mean = mean(strict_orientation_deduplicated_events),
  median = median(strict_orientation_deduplicated_events), .groups = "drop")
single_event_share <- density %>% group_by(group) %>%
  summarise(share = mean(strict_orientation_deduplicated_events == 1), .groups = "drop")
# Nested queries (Illumina embedded adapters; ONT 24-nt barcode inside K40/K46) count once per query (nesting report), so the
# single-event share can fall below 50%; the y-axis floor follows the data.
ecdf_floor <- floor(10 * min(single_event_share$share)) / 10
density_plot <- ggplot(density, aes(strict_orientation_deduplicated_events, colour = group)) +
  stat_ecdf(geom = "step", pad = FALSE, linewidth = 0.6) + scale_x_log10() +
  scale_y_continuous(breaks = seq(ecdf_floor, 1, 0.1), labels = scales::percent) +
  coord_cartesian(ylim = c(ecdf_floor, 1)) +
  scale_colour_manual(values = colours, labels = setNames(paste0(as.character(density_summary$group), " (N=", format(density_summary$accessions, big.mark = ",", trim = TRUE), ")"),
    as.character(density_summary$group)), name = NULL) +
    labs(title = "B  Match-event counts in hit-positive accessions", subtitle = paste0("Strict matches of testable queries in all records. Single-event share: ",
        floor(100 * min(single_event_share$share)), "\u2013", floor(100 * max(single_event_share$share)), "% per group. Y-axis from ", round(100 * ecdf_floor), "%.\n",
        "Events are orientation-deduplicated per query. Nested queries (e.g. K40/K46 contexts) add events at one site."),
      x = "Events per accession [log10]", y = "Cumulative\nfraction")
# C: every window x minimum-length cell for every unit variant (no cell selection); the 25-nt window is diagnostic.
grid_export <- export_csv(load_table("position_sensitivity") %>%
  arrange(match(variant, names(variant_labels)), match(group, c(groups, "Global union")), minimum_length, window) %>%
  transmute(variant, unit_rule, group, minimum_length_nt = minimum_length, window_nt = window, events, terminal, expected, ratio,
            diagnostic_window = window == 25, units_without_project), "figS3C_positional_sensitivity.csv")
grid <- grid_export %>% filter(group != "Global union") %>%
  mutate(group = factor(group, levels = groups), variant = factor(unname(variant_labels[variant]), levels = unname(variant_labels)),
    cell = paste0(minimum_length_nt / 1000, " kb / ", window_nt, " nt", ifelse(diagnostic_window, " (diag.)", "")),
    x = ifelse(is.na(ratio), 0, log10(ratio + 1)), state = ifelse(is.na(ratio), "No eligible events", "Ratio"))
cell_levels <- grid %>% distinct(minimum_length_nt, window_nt, cell) %>% arrange(desc(minimum_length_nt), desc(window_nt)) %>% pull(cell)
grid <- grid %>% mutate(cell = factor(cell, levels = cell_levels))
sensitivity_plot <- ggplot(grid, aes(x, cell, colour = group, shape = variant)) +
  geom_vline(xintercept = log10(2), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(data = filter(grid, state == "Ratio"), size = 1.4, stroke = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(data = filter(grid, state != "Ratio"), shape = 4, colour = grey_mid, size = 1.4, position = position_dodge(width = 0.7)) +
  facet_wrap(~group, nrow = 1) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_shape_manual(values = variant_shapes, name = NULL) +
  scale_x_continuous(limits = c(-0.05, 3.1), breaks = log10(c(0, 1, 10, 100, 1000) + 1), labels = c("0", "1", "10", "100", "1k")) +
  labs(title = "C  Complete positional sensitivity",
       subtitle = paste0("All minimum-length / window cells and all unit variants. Values: exported CSV. Gray crosses: no eligible events.\n",
         "diag.: 25-nt diagnostic window. Dashed line: ratio 1. A record without a BioProject is its own unit."),
       x = "Observed / expected [log10(ratio + 1)]", y = "Minimum length / window") +
  theme(axis.text.y = element_text(size = 7.5), legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), strip.text = element_text(size = 8))
figureS3_heights <- c(75, separator_height, 60, separator_height, 102)
figure_s4 <- (independent_panel(with_margin(position_scatter, bottom = 0)) / panel_separator() /
    independent_panel(with_margin(density_plot, top = 2, bottom = 0)) / panel_separator() /
    independent_panel(with_margin(sensitivity_plot, top = 2))) + plot_layout(heights = figureS3_heights)
save_figure("figureS3", figure_s4, sum(figureS3_heights), c("Distance versus length", "Conditional HSP density", "Positional sensitivity by unit"),
  c("position_events.tsv", "density.tsv", "position_sensitivity.tsv"),
  "Strict, a priori testable events of every condition. Sparse events are drawn as exact points",
  "figS3C_positional_sensitivity.csv")

# ---------------------------------------------------------------------------------------------
# Supplementary Tables re-exported from the same source (Word tables are regenerated from these)
# ---------------------------------------------------------------------------------------------
precision_export <- function(measure_name) {
  precision %>% filter(draws_per_parent == final_prefix, condition == "all_conditions", measure == measure_name) %>%
    mutate(group = factor(group, levels = c(groups, "Global union"))) %>% arrange(group, endpoint, population) %>%
    transmute(group = as.character(group), rule = endpoint, population, planned_parents, matched_parents,
      mean = estimate, range_lower = lower, range_upper = upper, half_width, target_half_width,
      control_draws, unique_seeds, target_limited_parents, precision_met, status)
}
invisible(export_csv(precision_export("detection"), "tableS6_detection_precision.csv"))
invisible(export_csv(precision_export("breadth"), "tableS6_breadth_precision.csv"))
invisible(export_csv(data.frame(arm = benchmark$search_space, arm_label = benchmark$arm_label, rule = benchmark$endpoint,
  truth_positive = benchmark$positive_pairs, recovered = benchmark$true_positive, missed = benchmark$false_negative,
  truth_negative = benchmark$negative_pairs, accepted_negative = benchmark$false_positive, recall = benchmark$recall),
  "tableS9_independent_benchmark.csv"))
invisible(export_csv(spike_summary %>% transmute(database, length_nt, canonical_parents = parents, exact_spikes = planned_spikes,
    calibrated_e10_recovered = calibrated_recovered, diagnostic_e1e12_recovered = diagnostic_recovered,
    exact_evalue_at_database_scale = exact_evalue, effective_search_space), "tableS9_exact_spike_detectability.csv"))

if (length(saved)) {
  suffix <- if (selection == "all") "" else paste0("_", selection)
  # PNG files are byte-deterministic; cairo writes the current time into every PDF, so PDFs are listed by name and their
  # content digest (creation time excluded) is added by build_figures.sh.
  filenames <- paste0("paper3_", names(saved), ".png")
  csv_files <- sort(list.files(tables_directory, pattern = "\\.csv$"))
  result <- list(contract = "publication_figure_artwork_v1",
    source_manifest_sha256 = digest(file.path(input, "manifest.json"), algo = "sha256", file = TRUE),
    plot_script_sha256 = digest(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1]), algo = "sha256", file = TRUE),
    theme_sha256 = digest(theme_file, algo = "sha256", file = TRUE), figures = saved,
    R_version = R.version.string,
    package_versions = setNames(lapply(c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "jsonlite", "digest", "ragg", "hexbin"),
      function(package) as.character(packageVersion(package))), c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "jsonlite", "digest", "ragg", "hexbin")),
    files = setNames(lapply(filenames, function(name) digest(file.path(output, name), algo = "sha256", file = TRUE)), filenames),
    pdf_files = paste0("paper3_", names(saved), ".pdf"),
    tables = setNames(lapply(csv_files, function(name) digest(file.path(tables_directory, name), algo = "sha256", file = TRUE)), csv_files))
  write_json(result, file.path(output, paste0("paper3_figure_manifest", suffix, ".json")), auto_unbox = TRUE, pretty = TRUE)
}
