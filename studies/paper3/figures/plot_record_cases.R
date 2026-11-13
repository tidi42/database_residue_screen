suppressPackageStartupMessages({
  library(grid)
  library(readr)
  library(jsonlite)
  library(digest)
})
options(warn = 2)
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% c(2, 3))
input <- normalizePath(arguments[1])
output <- arguments[2]
manifest <- fromJSON(file.path(input, "manifest.json"))
# Record cases (figureS4, manuscript Figure S2) from one figure-source directory (contract publication_figure_sources_v1): the
# cases selected by the record_cases stage with their sequence-level verification. US spelling (gray, color); coordinate ranges
# use an en dash and name the display orientation ("minus orientation", as in Figure 5B-D).
stopifnot(manifest$contract == "publication_figure_sources_v1")
for (name in names(manifest$files)) {
  stopifnot(!grepl("(^/|(^|/)\\.\\.(/|$))", name),
    digest(file.path(input, name), algo = "sha256", file = TRUE) == manifest$files[[name]])
}
cases <- read_tsv(file.path(input, "record_cases.tsv"), show_col_types = FALSE, na = c("", "not_estimated"))
cases <- cases[!is.na(cases$accession_version), ]
# Sequence-level verification of every case; species names (Latin binomials) are set in italics.
verification <- read_tsv(file.path(input, "record_case_verification.tsv"), show_col_types = FALSE, na = c("", "not_estimated"))
stopifnot(nrow(cases) <= 3, !anyDuplicated(cases$accession_version))
if (nrow(cases)) {
  stopifnot(identical(verification$accession_version, cases$accession_version), all(verification$displayed_window_equals_record),
    all(verification$record_window_equals_hsp), all(verification$hsp_full_length_exact),
    all(verification$packet_identical | is.na(verification$packet_identical)))
  cases$organism_italic <- verification$organism_italic
  cases$packet_identical <- verification$packet_identical
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
width_mm <- 180
# An optional third argument selects the compact layout (journal limit of 240 mm). The default reproduces the 260-mm release.
compact <- length(arguments) == 3
stopifnot(!compact || arguments[3] == "compact")
height_mm <- if (compact) 240 else 260
case_pitch <- if (compact) 67 else 75
case_start <- if (compact) 19 else 21
alignment_offset <- if (compact) 36 else 37
alignment_pitch <- if (compact) 10 else 12
track_offset <- if (compact) 26.5 else 26
interpretation_offset <- if (compact) 56.5 else 61
separator_offset <- if (compact) 64.5 else 70
# Colour is reserved for platforms, so query segments use neutral greys (as in scripts/plot_paper3_figures.R).
role_colours <- c(barcode = "#8E2043", scaffold = "#3E4E5C", unknown = "#7F7F7F", linker = "#6E7B12", flank = "#8C8C8C")
role_fills <- c(barcode = "#F2DCE2", scaffold = "#DDE3E8", unknown = "#F2F2F2", linker = "#E6EAC2", flank = "white")
# Match markers use the platform colour of the query (as in Figures 3 and 5). Segment roles use the non-platform family.
platform_colour <- function(group) {
  colour <- c(ONT = "#0072B2", Illumina = "#D55E00", PacBio = "#009E73", "BGI/DNBSEQ" = "#CC79A7")[group]
  stopifnot(!anyNA(colour))
  unname(colour)
}

write_text <- function(text, xpos, top, size = 8, bold = FALSE, colour = "black", family = "sans", just = "left") {
  grid.text(text, x = unit(xpos, "mm"), y = unit(height_mm - top, "mm"), just = c(just, "centre"),
    gp = gpar(fontsize = size, fontface = if (bold) "bold" else "plain", col = colour, fontfamily = family))
}

write_species <- function(organism, rest, xpos, top, size, italic) {
  if (!isTRUE(italic)) return(write_text(paste0(organism, rest), xpos, top, size))
  gp <- gpar(fontsize = size, fontface = "italic", fontfamily = "sans")
  width <- convertWidth(grobWidth(textGrob(organism, gp = gp)), "mm", valueOnly = TRUE)
  grid.text(organism, x = unit(xpos, "mm"), y = unit(height_mm - top, "mm"), just = c("left", "centre"), gp = gp)
  write_text(rest, xpos + width, top, size)
}

draw_case <- function(row, index, top) {
  category_label <- switch(row$category, terminal = "Terminal match", synthetic_or_vector = "Synthetic or vector record",
    internal_or_ambiguous = "Internal match", row$category)
  write_text(paste0(LETTERS[index], "  ", category_label, ": ", row$accession_version), 6, top, 10, TRUE)
  roles <- strsplit(row$display_segment_roles, ";", fixed = TRUE)[[1]]
  stopifnot(length(roles) == nchar(row$qseq), row$qseq == row$sseq)
  # Panel segment kinds: barcodes and indexes are drawn as barcode, all other designed segments as scaffold.
  roles[roles == "index"] <- "barcode"
  roles[roles %in% c("front", "rear", "hairpin", "oligo", "a_tail")] <- "scaffold"
  roles[!roles %in% names(role_colours)] <- "unknown"
  query_label <- paste0("Query: ", row$query_id, "  |  ", nchar(row$qseq), "-nt strict exact match")
  write_text(query_label, 6, top + 5, 7.7)
  title_lines <- strwrap(row$title, width = 100)
  for (line_index in seq_along(title_lines)) {
    line <- title_lines[line_index]
    if (line_index == 1 && isTRUE(row$organism_italic) && startsWith(line, row$organism)) {
      write_species(row$organism, substring(line, nchar(row$organism) + 1), 6, top + 9.5, 7.3, TRUE)
    } else {
      write_text(line, 6, top + 9.5 + 3.2 * (line_index - 1), 7.3)
    }
  }
  metadata <- paste(row$molecule, row$topology, row$record_type, paste("origin:", row$origin),
    paste("linear-genomic:", row$linear_eligibility), paste("taxid", row$taxids), sep = " | ")
  write_text(metadata, 6, top + 16.3, 7.1)
  coordinate_label <- paste0("Subject ", format(row$lower, big.mark = ",", scientific = FALSE), "\u2013",
    format(row$upper, big.mark = ",", scientific = FALSE), " (", ifelse(row$display_orientation == "-", "minus", "plus"), " orientation) / ",
    format(row$record_length, big.mark = ",", scientific = FALSE), " nt, ",
    if (row$topology == "circular") "circular record, position 1 is its stored origin" else
      paste(format(row$distance_to_end, big.mark = ",", scientific = FALSE), "nt from nearest record end"))
  write_text(coordinate_label, 6, top + 21, 7.1)
  grid.lines(x = unit(c(16, 168), "mm"), y = unit(rep(height_mm - top - track_offset, 2), "mm"), gp = gpar(col = "grey60", lwd = 1))
  match_x <- 16 + 152 * ((row$lower + row$upper) / 2 - 1) / max(1, row$record_length - 1)
  grid.lines(x = unit(rep(match_x, 2), "mm"), y = unit(height_mm - top - track_offset - c(-2, 2), "mm"), gp = gpar(col = platform_colour(row$group), lwd = 2.4))
  write_text("1", 16, top + track_offset + 4, 6.5)
  write_text(format(row$record_length, big.mark = ",", scientific = FALSE), 168, top + track_offset + 4, 6.5, just = "right")
  left <- if (is.na(row$display_left_flank)) "" else row$display_left_flank
  right <- if (is.na(row$display_right_flank)) "" else row$display_right_flank
  # Two alignment lines of 60 bases: longer queries show correspondingly shorter flanks.
  flank <- max(0, floor((120 - nchar(row$qseq)) / 2))
  left <- if (flank) substring(left, max(1, nchar(left) - flank + 1)) else ""
  right <- if (flank) substring(right, 1, flank) else ""
  subject <- strsplit(paste0(left, row$sseq, right), "", fixed = TRUE)[[1]]
  query <- c(rep(".", nchar(left)), strsplit(row$qseq, "", fixed = TRUE)[[1]], rep(".", nchar(right)))
  base_roles <- c(rep("flank", nchar(left)), roles, rep("flank", nchar(right)))
  stopifnot(length(subject) == length(query), length(query) == length(base_roles), length(subject) <= 120)
  for (offset in seq(1, length(subject), by = 60)) {
    indices <- offset:min(length(subject), offset + 59)
    line_index <- (offset - 1) / 60
    subject_top <- top + alignment_offset + alignment_pitch * line_index
    write_text("S", 6, subject_top, 8, TRUE, family = "mono")
    write_text("Q", 6, subject_top + 4.5, 8, TRUE, family = "mono")
    for (position in seq_along(indices)) {
      base <- indices[position]
      xpos <- 14 + 2.6 * (position - 1)
      if (base_roles[base] != "flank") {
        grid.rect(x = unit(xpos, "mm"), y = unit(height_mm - subject_top - 2.1, "mm"),
          width = unit(2.6, "mm"), height = unit(9, "mm"), gp = gpar(fill = role_fills[base_roles[base]], col = NA))
      }
      write_text(subject[base], xpos, subject_top, 8, colour = role_colours[base_roles[base]], family = "mono", just = "centre")
      write_text(query[base], xpos, subject_top + 4.5, 8,
        colour = if (base_roles[base] == "flank") "grey80" else role_colours[base_roles[base]], family = "mono", just = "centre")
    }
  }
  # Verification result of this case (every check is asserted above from record_case_verification.tsv), one sentence per line.
  lines <- c("Verified: exact full-length HSP of the searched query. The record sequence at these coordinates equals the HSP.",
    if (isTRUE(row$packet_identical)) "The archived record sequence is identical.")
  for (line_index in seq_along(lines)) write_text(lines[line_index], 6, top + interpretation_offset + 3.4 * (line_index - 1), 7.3)
  if (index < nrow(cases)) grid.lines(x = unit(c(6, 174), "mm"), y = unit(rep(height_mm - top - separator_offset, 2), "mm"), gp = gpar(col = "grey80", lwd = 0.5))
}

draw_figure <- function() {
  grid.newpage()
  write_text("Record cases", 6, 6, 11, TRUE)
  write_text("Deterministically selected strict matches with exact-version local sequence verification.", 6, 12, 7.8)
  if (nrow(cases)) for (index in seq_len(nrow(cases))) draw_case(cases[index, ], index, case_start + case_pitch * (index - 1))
  footer_top <- if (compact) case_start + 3 * case_pitch + 3 else 246
  write_text("S: subject. Q: query in the searched (canonical) orientation.", 6, footer_top, 7.1)
  write_text("Gray bases: record flanks (up to 50 nt), clipped at stored record boundaries. Slate: scaffold. Burgundy: barcode or index. Marker: platform color.", 6, footer_top + 5, 7.1)
  write_text("Selection rules and complete case evidence: Table S8.", 6, footer_top + 10, 6.8)
}

for (extension in c("pdf", "png")) {
  path <- file.path(output, paste0("paper3_figureS4.", extension))
  stopifnot(!file.exists(path))
  if (extension == "pdf") grDevices::cairo_pdf(path, width = width_mm / 25.4, height = height_mm / 25.4)
  else ragg::agg_png(path, width = width_mm, height = height_mm, units = "mm", res = if (compact) 300 else 400)
  draw_figure()
  grDevices::dev.off()
}
files <- "paper3_figureS4.png"
report <- list(contract = "publication_record_case_artwork_v1",
  source_manifest_sha256 = digest(file.path(input, "manifest.json"), algo = "sha256", file = TRUE),
  figures = list(figureS4 = list(width_mm = width_mm, height_mm = height_mm, selected_cases = nrow(cases),
    source_tables = c("record_cases.tsv", "record_case_verification.tsv"))),
  code_sha256 = digest(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1]), algo = "sha256", file = TRUE),
  files = setNames(lapply(files, function(name) digest(file.path(output, name), algo = "sha256", file = TRUE)), files),
  pdf_files = "paper3_figureS4.pdf")
write_json(report, file.path(output, "paper3_case_figure_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
message("Rendered figureS4")