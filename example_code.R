# ============================================================
# LIBRERIE
# ============================================================
library(synthdid)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# ============================================================
# DATI (modifica solo il path)
# ============================================================
datasdid <- read.csv(
  "data/datasdid.csv",
  stringsAsFactors = FALSE
)
# filtro come nel tuo codice
datasdid <- subset(
  datasdid,
  db_selected == 1 & tolower(trimws(as.character(territorio_enc))) != "caltanissetta"
)

# ============================================================
# NOMI VARIABILI (come nel tuo codice)
# ============================================================
UNIT <- "territorio_enc"
YVAR <- "index_pc_land"
TIME <- "anno"
DVAR <- "did"

# ============================================================
# OPZIONI (qui fai i “robustness” senza toccare tutto)
# ============================================================
DROP_YEAR <- NULL     # es: 2016 per eliminare il 2016; oppure NULL
TREAT_START_TRUE <- 2017
DO_TIME_PLACEBO <- FALSE   # tu non lo vuoi: lascia FALSE

K_LOO <- 50
B_SPACE <- 50
SEED <- 123

SAVE_INTERMEDIATE <- TRUE
CACHE_DIR <- "cache_sdid"
TAG <- "MAIN"

# ============================================================
# PREP DATA
# ============================================================
if (!is.null(DROP_YEAR)) {
  datasdid <- datasdid %>% filter(.data[[TIME]] != DROP_YEAR)
}

datasdid[[UNIT]] <- droplevels(as.factor(datasdid[[UNIT]]))
datasdid[[TIME]] <- as.integer(datasdid[[TIME]])
datasdid[[YVAR]] <- as.numeric(datasdid[[YVAR]])
datasdid[[DVAR]] <- as.integer(datasdid[[DVAR]])

# start trattamento (robusto)
treat_start <- min(datasdid[[TIME]][datasdid[[DVAR]] == 1], na.rm = TRUE)

# ============================================================
# HELPERS
# ============================================================
sdid_paths <- function(df, treat_var = DVAR) {
  df <- as.data.frame(df)
  df[[UNIT]] <- droplevels(as.factor(df[[UNIT]]))
  df[[TIME]] <- as.integer(df[[TIME]])
  df[[YVAR]] <- as.numeric(df[[YVAR]])
  df[[treat_var]] <- as.integer(df[[treat_var]])
  
  s <- synthdid::panel.matrices(
    df,
    unit = UNIT,
    time = TIME,
    outcome = YVAR,
    treatment = treat_var
  )
  
  est <- synthdid::synthdid_estimate(s$Y, s$N0, s$T0)
  
  Y  <- s$Y
  N0 <- s$N0
  
  years <- suppressWarnings(as.integer(colnames(Y)))
  if (all(is.na(years))) years <- sort(unique(df[[TIME]]))
  
  treated_ts <- colMeans(Y[(N0 + 1):nrow(Y), , drop = FALSE])
  
  wtab <- synthdid::synthdid_controls(est, weight.type = "omega", mass = 1)
  omega <- rep(0, N0); names(omega) <- rownames(Y)[1:N0]
  omega[rownames(wtab)] <- wtab[, 1]
  if (sum(omega) > 0) omega <- omega / sum(omega)
  
  synthetic_ts <- as.numeric(t(omega) %*% Y[1:N0, , drop = FALSE])
  
  out <- data.frame(
    year = years,
    treated = as.numeric(treated_ts),
    synthetic = as.numeric(synthetic_ts)
  )
  out$gap <- out$treated - out$synthetic
  
  list(series = out, est = est)
}

center_gap_refyear <- function(df, ref_year = 2013) {
  ref <- df$gap[df$year == ref_year]
  if (length(ref) == 0 || is.na(ref[1])) stop(paste("Ref year", ref_year, "not found in series"))
  df$gap_c <- df$gap - ref[1]
  df
}

theme_paper <- theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill = scales::alpha("white", 0.75), color = NA),
    legend.key = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

# =========================
# A) BASELINE
# =========================
db <- datasdid  # <-- usa il tuo dataset qui

treat_start <- min(db[[TIME]][db[[DVAR]] == 1], na.rm = TRUE)

base_out <- sdid_paths(db, DVAR)
base     <- base_out$series
tau_hat  <- base_out$est
se       <- sqrt(vcov(tau_hat, method = "placebo"))

pA <- ggplot() +
  geom_line(data = base, aes(year, synthetic, color = "Synthetic"),
            linetype = "dashed", linewidth = 1) +
  geom_point(data = base, aes(year, synthetic, color = "Synthetic"), size = 2) +
  geom_line(data = base, aes(year, treated, color = "Treated"), linewidth = 1) +
  geom_point(data = base, aes(year, treated, color = "Treated"), size = 2) +
  geom_vline(xintercept = treat_start, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Treated" = "#2297E6", "Synthetic" = "salmon")) +
  scale_x_continuous(breaks = seq(min(base$year, na.rm=TRUE), max(base$year, na.rm=TRUE), by = 2)) +
  labs(x = "Year", y = "LCI") +
  theme_paper

treated_units <- names(which(tapply(db[[DVAR]], db[[UNIT]], max, na.rm = TRUE) == 1))
control_units <- setdiff(levels(db[[UNIT]]), treated_units)

# =========================
# B) TIME PLACEBO (-1 anno) (se lo vuoi)
# =========================
placebo_year <- TREAT_START_TRUE - 1

if (DO_TIME_PLACEBO && placebo_year %in% unique(db[[TIME]])) {
  df_ypl <- db
  df_ypl$did_year_pl <- as.integer(df_ypl[[UNIT]] %in% treated_units & df_ypl[[TIME]] >= placebo_year)
  
  ypl <- sdid_paths(df_ypl, "did_year_pl")$series
  
  pB <- ggplot() +
    geom_line(data = ypl, aes(year, synthetic, color = "Synthetic"),
              linetype = "dashed", linewidth = 1) +
    geom_point(data = ypl, aes(year, synthetic, color = "Synthetic"), size = 2) +
    geom_line(data = ypl, aes(year, treated, color = "Treated"), linewidth = 1) +
    geom_point(data = ypl, aes(year, treated, color = "Treated"), size = 2) +
    geom_vline(xintercept = placebo_year, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("Treated" = "#2297E6", "Synthetic" = "salmon")) +
    scale_x_continuous(breaks = seq(min(ypl$year, na.rm=TRUE), max(ypl$year, na.rm=TRUE), by = 2)) +
    labs(x = "Year", y = "LCI") +
    theme_paper
} else {
  pB <- ggplot() + theme_void()
}

pB

# =========================
# C) LOO (grigio)
# =========================
tab_omega <- synthdid::synthdid_controls(tau_hat, weight.type = "omega", mass = 0.99)
donors_plot <- rownames(tab_omega)[1:min(K_LOO, nrow(tab_omega))]

loo_list <- vector("list", length(donors_plot))
names(loo_list) <- donors_plot

for (j in seq_along(donors_plot)) {
  d <- donors_plot[j]
  df_sub <- db %>% dplyr::filter(.data[[UNIT]] != d)
  
  res <- tryCatch(sdid_paths(df_sub, DVAR)$series, error = function(e) NULL)
  if (!is.null(res)) {
    loo_list[[j]] <- data.frame(year = res$year, synthetic = res$synthetic, loo = d)
  }
}

loo_syn <- dplyr::bind_rows(loo_list)
treat_start <- 2016

pC <- ggplot() +
  geom_line(data = loo_syn,
            aes(year, synthetic, group = loo, color = "LOO synthetic"),
            linewidth = 0.9, alpha = 0.35) +
  geom_line(data = base, aes(year, synthetic, color = "Synthetic"), linewidth = 1) +
  geom_line(data = base, aes(year, treated,   color = "Treated"),   linewidth = 1) +
  geom_vline(xintercept = treat_start, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Treated"="#2297E6", "Synthetic"="salmon", "LOO synthetic"="grey50")) +
  scale_x_continuous(breaks = seq(min(base$year, na.rm=TRUE), max(base$year, na.rm=TRUE), by = 2)) +
  labs(x = "", y = "Land Consuption Index") +
  theme_paper  +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1.2, "lines"),
    panel.border = element_blank(),                 # rimuove il box esterno
    axis.line.x.bottom = element_line(),            # mantiene solo assi bottom/left
    axis.line.y.left   = element_line()
  )


pC
# =========================
# D) SPACE PLACEBO (gap centrato su 2013)
# =========================
set.seed(SEED)
n_tr <- length(treated_units)
if (length(control_units) < n_tr) stop("Controlli < trattati: impossibile campionare placebo.")

# gap “vero” centrato su 2013
true_gap <- base %>%
  center_gap_refyear(ref_year = 2013) %>%
  dplyr::transmute(year, gap_c)

space_list <- vector("list", B_SPACE)
for (b in 1:B_SPACE) {
  fake_group <- sample(control_units, n_tr, replace = FALSE)
  
  df_sp <- db
  df_sp$did_space_pl <- as.integer(df_sp[[UNIT]] %in% fake_group & df_sp[[TIME]] >= treat_start)
  
  sp <- sdid_paths(df_sp, "did_space_pl")$series %>%
    center_gap_refyear(ref_year = 2013)
  
  space_list[[b]] <- data.frame(year = sp$year, gap_c = sp$gap_c, id = b)
}
space_gap <- dplyr::bind_rows(space_list)


pD <- ggplot() +
  geom_line(data = space_gap,
            aes(year, gap_c, group = id, color = "Synthetic"),
            linetype = "dashed", alpha = 0.35, linewidth = 0.7) +
  geom_line(data = true_gap,
            aes(year, gap_c, color = "Treated"),
            linewidth = 1.1) +
  geom_vline(xintercept = treat_start, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Treated"="#2297E6", "Synthetic"="black")) +
  scale_x_continuous(breaks = seq(min(true_gap$year, na.rm=TRUE), max(true_gap$year, na.rm=TRUE), by = 2)) +
  labs(x = "", y = "Gap") +
  theme_paper +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1.2, "lines"),
    panel.border = element_blank(),                 # rimuove il box esterno
    axis.line.x.bottom = element_line(),            # mantiene solo assi bottom/left
    axis.line.y.left   = element_line()
  )

pD

grid_out <- gridExtra::arrangeGrob(pC, pD, ncol = 1)
grid::grid.newpage(); grid::grid.draw(grid_out)



p <- synthdid_plot(tau_hat)




legend_font_size <- 18

p_final <- p +
  coord_cartesian(ylim = c(-.01, .03)) +
  scale_y_continuous(breaks = seq(-.01, .03, by = 0.010)) +
  geom_vline(xintercept = 2016, color = "black", size = 1) +
  scale_x_continuous(
    breaks = function(x) pretty(x, n = 6),
    labels = function(x) as.integer(round(x))
  ) +
  labs(x = "", y = "Land Consuption Index") +
  theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 20),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = legend_font_size),
    legend.title = element_blank()
  )


p_final


# ============================================================
# DONOR WEIGHTS (omega) - tabella semplice per referee
# ============================================================


tab_omega_all <- synthdid::synthdid_controls(tau_hat, weight.type = "omega", mass = 1)

# tab_omega_all: righe = donors con peso positivo, colonna 1 = omega
donor_weights <- data.frame(
  donor = rownames(tab_omega_all),
  omega = as.numeric(tab_omega_all[, 1]),
  row.names = NULL
)


donor_weights$omega <- donor_weights$omega / sum(donor_weights$omega)


donor_weights <- donor_weights %>% arrange(desc(omega))


cat("Numero di donors con peso > 0:", nrow(donor_weights), "\n")
cat("Somma dei pesi omega:", sum(donor_weights$omega), "\n")
cat("Top 10 donors:\n")
print(head(donor_weights, 10))




ggplot(donor_weights, aes(x = omega)) +
  geom_histogram(bins = 50, fill = "#2297E6", color = "white") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 0.02),
    breaks = seq(0, 0.001, by = 0.0002) 
  ) +
  scale_y_continuous(
    limits = c(0, 1000),
    breaks = seq(0, 1000, by = 200)      
  ) +
  coord_cartesian(xlim = c(0, 0.001), ylim = c(0, 1000)) +
  labs(
    x = "Donor weight (ω)",
    y = "Number of municipalities"
  ) +
  theme_paper +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1.2, "lines"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left   = element_line()