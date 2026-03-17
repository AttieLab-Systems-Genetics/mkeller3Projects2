# =========================================================
# LOCKED Manhattan Plot Script (GRCm39 fixed chr widths)
# - Assumes input CSV has EXACT columns: phenotype, Chr, Mbp, LOD
# - Chr must be: 1..19 or X (character or numeric)
# - Mbp must be MEGABASES (0–~200), not bp
# - Plots to the current plot window by default
# - Only thing you need to edit each time: input_file
# =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =========================================================
# 1) EDIT ONLY THIS LINE
# =========================================================
input_file <- "sexdietdiff.csv"   # e.g., "dietdiff.csv", "sexdiff.csv", "sexdietdiff.csv"

# =========================================================
# 2) OPTIONAL PLOT CONTROLS (safe to ignore)
# =========================================================
#plot_title <- paste0("Manhattan: ", tools::file_path_sans_ext(basename(input_file)))
x_label <- "Chromosome"
y_label <- "LOD"

# Y-axis limits (set to NULL for automatic)
y_min <- 9
y_max <- 50

# Point style
point_size    <- 6.0
outline_color <- "black"
outline_width <- 0.25
alpha_points  <- 1.0

# Alternating chromosome colors
chr_colors <- c("royalblue3", "orange2")

# Optional horizontal threshold line
threshold_y <- NULL  # e.g., 7

# Optional save (OFF by default)
save_file   <- TRUE
output_file <- "manhattan.png"  # ".png" or ".pdf"
plot_width  <- 16
plot_height <- 4
plot_dpi    <- 600

# =========================================================
# 3) GRCm39 chromosome lengths (bp) — fixed widths
# =========================================================
chr_levels <- c(as.character(1:19), "X")

grcm39_chr_len_bp <- c(
  "1"  = 195154279, "2"  = 181755017, "3"  = 159745316, "4"  = 156860686, "5"  = 151758149,
  "6"  = 149588044, "7"  = 144995196, "8"  = 130127694, "9"  = 124359700, "10" = 130530862,
  "11" = 121973369, "12" = 120092757, "13" = 120883175, "14" = 125139656, "15" = 104073951,
  "16" = 98008968,  "17" = 95294699,  "18" = 90720763,  "19" = 61420004,  "X"  = 169476592
)

chr_sizes <- tibble(
  Chr = chr_levels,                         # character join key
  chr_len = as.numeric(grcm39_chr_len_bp)   # bp
) %>%
  mutate(
    chr_offset = cumsum(chr_len) - chr_len,
    center     = chr_offset + chr_len / 2
  )

total_genome_len <- sum(chr_sizes$chr_len)

# =========================================================
# 4) READ + STRICT VALIDATION
# =========================================================
df <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c("phenotype", "Chr", "Mbp", "LOD")
missing <- setdiff(required_cols, names(df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "),
       "\nFound columns: ", paste(names(df), collapse = ", "))
}

# Normalize Chr safely (handles numeric 1.0 -> "1", chr1 -> "1", 20 -> X)
normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x <- gsub("^chromosome\\s*", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x <- ifelse(x == "20", "X", x)
  # If numeric-like, coerce to integer string (e.g., "1.0" -> "1")
  num <- suppressWarnings(readr::parse_number(x))
  is_num <- !is.na(num) & !(x %in% c("X", "Y"))
  x[is_num] <- as.character(as.integer(round(num[is_num])))
  x
}

df <- df %>%
  mutate(
    phenotype = as.character(phenotype),
    Chr = normalize_chr(Chr),
    Mbp = as.numeric(Mbp),
    LOD = as.numeric(LOD)
  ) %>%
  filter(!is.na(Chr), !is.na(Mbp), !is.na(LOD))

# Hard check Chr values
bad_chr <- setdiff(sort(unique(df$Chr)), chr_levels)
if (length(bad_chr) > 0) {
  stop("Unexpected Chr values found: ", paste(bad_chr, collapse = ", "),
       "\nExpected only: ", paste(chr_levels, collapse = ", "))
}

# =========================================================
# 5) MAP TO FIXED GENOME COORDINATES (x-axis constant)
# =========================================================
df <- df %>%
  mutate(
    BP = Mbp * 1e6,                 # Mbp -> bp
    Chr = as.character(Chr)
  ) %>%
  left_join(chr_sizes %>% select(Chr, chr_offset), by = "Chr") %>%
  mutate(BP_cum = BP + chr_offset,
         Chr = factor(Chr, levels = chr_levels))

if (any(is.na(df$BP_cum))) stop("Join to chromosome offsets failed; check Chr formatting.")

# Alternating fill per chromosome
df <- df %>%
  mutate(
    chr_index = as.integer(Chr),
    fill_col  = chr_colors[(chr_index %% length(chr_colors)) + 1]
  )

# Y limits
ymin <- if (is.null(y_min)) min(df$LOD, na.rm = TRUE) else y_min
ymax <- if (is.null(y_max)) max(df$LOD, na.rm = TRUE) else y_max
if (is.null(y_min)) ymin <- min(0, ymin)

# =========================================================
# 6) PLOT (plot window; optional save)
# =========================================================
p <- ggplot(df, aes(x = BP_cum, y = LOD)) +
  geom_point(
    aes(fill = fill_col),
    shape = 21,
    color = outline_color,
    stroke = outline_width,
    size = point_size,
    alpha = alpha_points
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = chr_sizes$center,
    labels = chr_sizes$Chr,
    limits = c(0, total_genome_len),
    expand = expansion(mult = c(0.002, 0.01))  # small left, slightly larger right
  ) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  theme_classic(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(vjust = 0.5)
  ) +
  #labs(x = x_label, y = y_label, title = plot_title)
   labs(x = x_label, y = y_label)
if (!is.null(threshold_y)) {
  p <- p + geom_hline(yintercept = threshold_y, linetype = 2, linewidth = 0.6, color = "grey30")
}

print(p)






if (isTRUE(save_file)) {
  if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    ggsave(output_file, plot = p, width = plot_width, height = plot_height, units = "in")
  } else {
    ggsave(output_file, plot = p, width = plot_width, height = plot_height, units = "in", dpi = plot_dpi)
  }
  message("Saved: ", output_file)
}
