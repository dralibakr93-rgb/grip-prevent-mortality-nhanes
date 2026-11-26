# Stage 6. Main and supplementary figures, drawn from the result files of stages 1 to 5d.
# Each figure is saved as PDF, SVG, 600-dpi TIFF, 300-dpi PNG, and a ggplot2 object (.rds).

invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
source("code/00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(survival)
})

# design tokens
FONT <- "Arial"
INK <- "#1A1A1A"
MUTED <- "#595959"
REFC <- "#8C8C8C"
RULE <- "#BDBDBD"
COL_S <- "#0B5A8E"
TINT_S <- "#DCE9F4" # measured grip strength (blue)
COL_S2 <- "#7FB2DA" # higher strength / comparison (light blue)
COL_M <- "#C2410C"
TINT_M <- "#F6E0D5" # reported mobility limitation (vermillion)
COL_B <- "#3D3D3D" # both measures
NEUTRAL <- "#F2F2F2"
BASE <- 8 # pt at final printed size
pt <- function(x) x / ggplot2::.pt # text size in pt -> ggplot size units
theme_pub <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = INK),
      axis.ticks = element_line(linewidth = 0.3, colour = INK),
      axis.ticks.length = unit(1.1, "mm"),
      axis.text = element_text(colour = INK, size = base),
      axis.title = element_text(colour = INK, size = base + 0.5),
      plot.title = element_text(
        colour = INK, size = base + 1, face = "bold", hjust = 0,
        margin = margin(0, 0, 2, 0, "mm")
      ),
      plot.subtitle = element_text(colour = MUTED, size = base, hjust = 0),
      plot.margin = margin(1.5, 2, 1.5, 1.5, "mm"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "none"
    )
}
fmt_n <- function(x) formatC(round(x), big.mark = ",", format = "d")
f2 <- function(x) formatC(x, format = "f", digits = 2)
f1 <- function(x) formatC(x, format = "f", digits = 1)
fxm <- function(x, d = 3) {
  v <- formatC(x, format = "f", digits = d)
  v <- sub("^-(0\\.0+)$", "\\1", v)
  gsub("-", "\u2212", v)
}
ci_txt <- function(h, l, u) sprintf("%s (%s\u2013%s)", f2(h), f2(l), f2(u))
ci_to <- function(e, l, u, d = 3) sprintf("%s (%s to %s)", fxm(e, d), fxm(l, d), fxm(u, d))
pv <- function(p) if (p < 0.001) "P<0.001" else if (p < 0.01) sprintf(
  "P=%s",
  formatC(p, format = "f", digits = 3)
) else sprintf(
  "P=%s",
  formatC(p, format = "f", digits = 2)
)
a2 <- function(name) read.csv(file.path(P$csv, "analysis_outputs", paste0(name, ".csv")),
  stringsAsFactors = FALSE
)
ao <- a2

save_fig <- function(p, name, w_mm, h_mm, dir) {
  w <- w_mm / 25.4
  h <- h_mm / 25.4
  quartz(
    type = "pdf",
    file = file.path(dir, paste0(name, ".pdf")),
    width = w,
    height = h,
    family = FONT
  )
  print(p)
  invisible(dev.off())
  svglite::svglite(file.path(dir, paste0(name, ".svg")), width = w, height = h)
  print(p)
  invisible(dev.off())
  ragg::agg_tiff(file.path(dir, paste0(name, ".tiff")),
    width = w, height = h, units = "in", res = 600,
    compression = "lzw", background = "white"
  )
  print(p)
  invisible(dev.off())
  ragg::agg_png(file.path(dir, paste0(name, ".png")),
    width = w, height = h, units = "in", res = 300,
    background = "white"
  )
  print(p)
  invisible(dev.off())
  saveRDS(p, file.path(P$fig_src, paste0(name, ".ggplot.rds")))
  message("Saved ", name, sprintf(" (%.0f x %.0f mm)", w_mm, h_mm))
}

# Generic forest panel: label column | optional n/deaths column | estimates | value column.
# d: label, header (logical), est, lo, hi, col, and optionally n_txt; header rows carry only a label.
forest <- function(
  d,
  xlim,
  breaks,
  xlab,
  ref = 1,
  log = TRUE,
  title = NULL,
  val_head = "HR (95% CI)",
  val_fun = function(e, l, u) ci_txt(e, l, u), n_head = "n / deaths", band = NULL,
  widths = c(0.40, 0.13, 0.29, 0.18), xlabels = NULL
) {
  d$y <- rev(seq_len(nrow(d)))
  e <- d[!d$header, ]
  yl <- c(0.45, nrow(d) + 1.05)
  hy <- nrow(d) + 0.85
  has_n <- "n_txt" %in% names(d)
  lab <- ggplot(d) +
    geom_text(aes(x = ifelse(header, 0, 0.04), y = y, label = label),
      hjust = 0, size = pt(BASE), family = FONT,
      fontface = ifelse(d$header, "bold", "plain"), colour = ifelse(d$header, d$col, INK)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = yl, clip = "off") +
    labs(title = title) +
    theme_void(base_family = FONT) +
    theme(
      plot.title = element_text(face = "bold", size = BASE + 1, colour = INK, hjust = 0, margin = margin(0, 0, 1.5, 0, "mm")),
      plot.margin = margin(1.5, 0, 1.5, 1.5, "mm")
    )
  nc <- ggplot(e) +
    geom_text(aes(x = 1, y = y, label = if (has_n) n_txt else ""),
      hjust = 1,
      size = pt(BASE - 0.5),
      family = FONT,
      colour = MUTED
    ) +
    annotate("text",
      x = 1,
      y = hy,
      label = if (has_n) n_head else "",
      hjust = 1,
      size = pt(BASE - 0.5),
      family = FONT,
      fontface = "bold",
      colour = MUTED
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = yl, clip = "off") +
    labs(title = if (is.null(title)) NULL else " ") +
    theme_void(base_family = FONT) +
    theme(
      plot.title = element_text(size = BASE + 1, margin = margin(0, 0, 1.5, 0, "mm")),
      plot.margin = margin(1.5, 2, 1.5, 0, "mm")
    )
  f <- ggplot(e, aes(est, y))
  if (!is.null(band)) f <- f + annotate("rect",
    xmin = band[1],
    xmax = band[2],
    ymin = -Inf,
    ymax = nrow(d) + 0.5,
    fill = TINT_S,
    alpha = 0.6
  )
  f <- f + geom_vline(xintercept = ref, linetype = "22", colour = REFC, linewidth = 0.35) +
    geom_errorbar(aes(xmin = lo, xmax = hi),
      width = 0.28,
      orientation = "y",
      colour = e$col,
      linewidth = 0.45
    ) +
    geom_point(shape = 21, size = 2.1, colour = e$col, fill = e$col, stroke = 0.5)
  f <- f + (if (log) scale_x_log10(
    breaks = breaks,
    labels = if (is.null(xlabels)) format(breaks, drop0trailing = TRUE, trim = TRUE) else xlabels,
    expand = c(0, 0)
  )
  else scale_x_continuous(
    breaks = breaks,
    labels = if (is.null(xlabels)) waiver() else xlabels,
    expand = c(0, 0)
  )) +
    scale_y_continuous(expand = c(0, 0)) + coord_cartesian(xlim = xlim, ylim = yl) +
    labs(x = xlab, y = NULL, title = if (is.null(title)) NULL else " ") + theme_pub() +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
      plot.margin = margin(1.5, 3, 1.5, 0, "mm")
    )
  v <- ggplot(e) +
    geom_text(aes(x = 0, y = y, label = val_fun(est, lo, hi)),
      hjust = 0,
      size = pt(BASE - 0.5),
      family = FONT,
      colour = INK
    ) +
    annotate("text",
      x = 0,
      y = hy,
      label = val_head,
      hjust = 0,
      size = pt(BASE - 0.5),
      family = FONT,
      fontface = "bold",
      colour = MUTED
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = yl, clip = "off") +
    labs(title = if (is.null(title)) NULL else " ") +
    theme_void(base_family = FONT) +
    theme(
      plot.title = element_text(size = BASE + 1, margin = margin(0, 0, 1.5, 0, "mm")),
      plot.margin = margin(1.5, 1.5, 1.5, 0, "mm")
    )
  lab + nc + f + v + plot_layout(widths = widths)
}
hdr <- function(
  label,
  col = INK
) data.frame(
  label = label,
  header = TRUE,
  est = NA,
  lo = NA,
  hi = NA,
  col = col,
  n_txt = ""
)
row_ <- function(label, r, col, est = "HR", lo = "lower", hi = "upper")
  data.frame(
    label = label, header = FALSE, est = r[[est]], lo = r[[lo]], hi = r[[hi]], col = col,
    n_txt = if (!is.null(r$N)) sprintf("%s / %s", fmt_n(r$N), fmt_n(r$deaths)) else ""
  )

# FIGURE 1 (cohort flow)
flow <- res("cohort_flow")
gs <- res("grip_selection")
mr <- res("mobility_routing")
fig1v <- data.frame(
  key = c(flow$stage, paste0("grip_", gsub(" ", "_", gs$grip_status)), paste0("routing_", seq_len(nrow(mr)))),
  N = c(flow$N, gs$N, mr$N), deaths = c(flow$deaths, gs$deaths, rep(NA, nrow(mr))),
  label = c(flow$stage, gs$grip_status, mr$group)
)
write.csv(fig1v, file.path(P$csv, "figure1_flow_values.csv"), row.names = FALSE)
N <- setNames(flow$N, flow$stage)
Dth <- setNames(flow$deaths, flow$stage)
G <- setNames(gs$N, gs$grip_status)
fs <- res("cohort_summary")
mp <- res("mobility_prevalence")
n_younger <- mr$N[grepl("younger than 60", mr$group)]
n_uncl <- mr$N[grepl("unclassified", mr$group)]
stopifnot(
  N[["core"]] - N[["mobility_age60_79"]] == n_younger + n_uncl,
  N[["linked_height"]] - N[["core"]] == sum(G[names(G) != "bilateral valid maximal effort"])
)
bx <- list()
tx <- list()
ar <- list()
box <- function(x, y, w, h, fill = "white", col = INK, lwd = 0.3)
  bx[[length(bx) + 1]] <<- data.frame(
    xmin = x - w / 2,
    xmax = x + w / 2,
    ymin = y - h / 2,
    ymax = y + h / 2,
    fill = fill,
    col = col,
    lwd = lwd
  )
txt <- function(x, y, label, size = 8, face = "plain", col = INK, hjust = 0.5)
  tx[[length(tx) + 1]] <<- data.frame(
    x = x,
    y = y,
    label = label,
    size = size,
    face = face,
    col = col,
    hjust = hjust
  )
arrow_seg <- function(
  x,
  y,
  xend,
  yend
) ar[[length(ar) + 1]] <<- data.frame(
  x = x,
  y = y,
  xend = xend,
  yend = yend
)
MX <- 50
MW <- 86
EX <- 135
EW <- 66
stage <- function(y, l1, n) {
  box(MX, y, MW, 10.5)
  txt(MX, y + 2.1, l1, 8)
  txt(MX, y - 2.4, sprintf("n = %s", fmt_n(n)), 8, "bold")
}
excl <- function(y, n, why, label = "Excluded") {
  box(EX, y, EW, 10, fill = NEUTRAL, col = REFC, lwd = 0.25)
  txt(EX - EW / 2 + 2.5, y + 2.0, sprintf("%s, n = %s", label, fmt_n(n)), 7.5, "bold", INK, 0)
  txt(EX - EW / 2 + 2.5, y - 2.2, why, 7.5, "plain", INK, 0)
  arrow_seg(MX, y, EX - EW / 2, y)
}
Y <- c(130, 111, 92, 73, 54)
E <- (Y[-1] + Y[-5]) / 2
YC <- 33.5
YM <- 9
stage(Y[1], "NHANES 2011\u20132014 participants", N[["all_respondents"]])
stage(Y[2], "Aged 50\u201379 years with examination data", N[["age_50_79"]])
stage(Y[3], "No self-reported cardiovascular disease", N[["no_CVD_not_pregnant"]])
stage(Y[4], "Complete, in-range PREVENT inputs", N[["score_range"]])
stage(Y[5], "Mortality linkage and height available", N[["linked_height"]])
excl(E[1], N[["all_respondents"]] - N[["age_50_79"]], "Other ages or no examination data")
excl(
  E[2],
  N[["age_50_79"]] - N[["no_CVD_not_pregnant"]],
  "Cardiovascular disease, unknown, or pregnancy"
)
excl(
  E[3],
  N[["no_CVD_not_pregnant"]] - N[["score_range"]],
  "Missing or out-of-range PREVENT inputs"
)
excl(E[4], N[["score_range"]] - N[["linked_height"]], "No mortality linkage or height")
excl((Y[5] + YC) / 2 + 1, N[["linked_height"]] - N[["core"]], "No valid grip test in both hands")
for (k in 1:4) arrow_seg(MX, Y[k] - 5.25, MX, Y[k + 1] + 5.25)
arrow_seg(MX, Y[5] - 5.25, MX, YC + 7)
box(MX, YC, MW, 14, fill = TINT_S, col = COL_S, lwd = 0.7)
txt(MX, YC + 3.6, "Primary cohort, age 50\u201379 years", 8.5, "bold")
txt(MX, YC - 0.4, sprintf("n = %s; %s deaths", fmt_n(N[["core"]]), fmt_n(Dth[["core"]])), 8)
txt(
  MX,
  YC - 4.2,
  sprintf("median follow-up %s years", formatC(fs$median_followup[fs$cohort == "core"], format = "f", digits = 2)),
  7.5,
  "plain",
  MUTED
)
arrow_seg(MX, YC - 7, MX, YM + 7)
excl((YC + YM) / 2,
  n_younger + n_uncl,
  label = "Not in older cohort",
  why = sprintf(
    "Age 50\u201359 years (%s) or mobility unknown (%s)",
    fmt_n(n_younger),
    fmt_n(n_uncl)
  )
)
box(MX, YM, MW, 14, fill = TINT_M, col = COL_M, lwd = 0.7)
txt(MX, YM + 3.6, "Older cohort, age 60\u201379 years", 8.5, "bold")
txt(
  MX,
  YM - 0.4,
  sprintf("n = %s; %s deaths", fmt_n(N[["mobility_age60_79"]]), fmt_n(Dth[["mobility_age60_79"]])),
  8
)
txt(
  MX,
  YM - 4.2,
  sprintf("%s (%s%%) with mobility limitation", fmt_n(mp$with_limitation), formatC(100 * mp$prevalence, format = "f", digits = 1)),
  7.5,
  "plain",
  MUTED
)
B <- do.call(rbind, bx)
T <- do.call(rbind, tx)
A <- do.call(rbind, ar)
fig1 <- ggplot() +
  geom_rect(
    data = B,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = B$fill,
    colour = B$col,
    linewidth = B$lwd
  ) +
  geom_segment(
    data = A, aes(x = x, y = y, xend = xend, yend = yend), colour = INK, linewidth = 0.3,
    arrow = arrow(length = unit(1.4, "mm"), type = "closed")
  ) +
  geom_text(
    data = T,
    aes(x = x, y = y, label = label),
    size = pt(T$size),
    fontface = T$face,
    colour = T$col,
    hjust = T$hjust, family = FONT
  ) +
  coord_fixed(xlim = c(0, 174), ylim = c(0, 138), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_fig(fig1, "Figure_1", 174, 138, P$fig_main)

# FIGURE 2 (strength and mobility: association and added discrimination)
fe <- res("focal_estimates")
get_est <- function(model, term) {
  r <- fe[fe$model == model & fe$term == term, ]
  stopifnot(nrow(r) == 1)
  r
}
p_prim <- get_est("Primary: continuous strength", "lower_z")
p_c <- get_est("Mobility C: add measured strength", "lower_z")
p_ds <- get_est("Mobility D: strength and limitation", "lower_z")
p_b <- get_est("Mobility B: add reported limitation", "mobility")
p_dm <- get_est("Mobility D: strength and limitation", "mobility")
stopifnot(p_c$N == p_ds$N, p_b$N == p_ds$N, p_c$deaths == p_dm$deaths)
d2 <- rbind(
  hdr("Lower grip strength (per SD)", COL_S),
  row_("All adults aged 50\u201379 years", p_prim, COL_S),
  row_("Adults aged 60\u201379 years", p_c, COL_S),
  row_("  + adjusted for mobility limitation", p_ds, COL_S),
  hdr("Mobility limitation (yes vs no)", COL_M),
  row_("Adults aged 60\u201379 years", p_b, COL_M),
  row_("  + adjusted for grip strength", p_dm, COL_M)
)
p2a <- forest(d2,
  c(0.8, 5),
  c(1, 1.5, 2, 3, 4),
  "Hazard ratio for death (log scale)",
  title = "A  Association with all-cause mortality"
)
pm <- a2("prediction_metrics")
pm_ <- function(co, me, mo) {
  r <- pm[pm$cohort == co & pm$metric == me & pm$model == mo, ]
  stopifnot(nrow(r) == 1)
  r
}
d4 <- rbind(
  hdr("Adults aged 50\u201379 years"),
  data.frame(
    label = "Add grip strength",
    header = FALSE,
    est = pm_("primary", "dC", "M1")$apparent,
    lo = pm_("primary", "dC", "M1")$lower,
    hi = pm_("primary", "dC", "M1")$upper, col = COL_S, n_txt = ""
  ),
  hdr("Adults aged 60\u201379 years"),
  do.call(rbind, lapply(
    list(c("Mobility", "Add mobility limitation", COL_M), c("Grip", "Add grip strength", COL_S), c("Both", "Add both", COL_B)),
    function(z) {
      r <- pm_("mobility", "dC", z[1])
      data.frame(
        label = z[2],
        header = FALSE,
        est = r$apparent,
        lo = r$lower,
        hi = r$upper,
        col = z[3],
        n_txt = ""
      )
    }
  ))
)
d4$n_txt <- NULL
p2b <- forest(d4,
  c(-0.02, 0.06),
  seq(-0.02, 0.06, 0.02),
  "Change in C statistic when the measure is added to the base model",
  ref = 0, log = FALSE,
  title = "B  Added ability to rank risk of death",
  val_head = "Change in C (95% CI)",
  val_fun = function(e, l, u) ci_to(e, l, u),
  widths = c(0.40, 0.01, 0.38, 0.21), xlabels = fxm(seq(-0.02, 0.06, 0.02), 2)
)
fig2 <- wrap_elements(p2a) / wrap_elements(p2b) + plot_layout(heights = c(
  nrow(d2) + 1.6,
  nrow(d4) + 2.4
))
save_fig(fig2, "Figure_2", 174, 122, P$fig_main)

# FIGURE 3 (weaker vs stronger grip on the absolute scale: death and daily function)
rr <- res("standardized_5yr_risks")
qq <- res("risk_strength_contrasts")
R <- split(rr, rr$estimand)
stopifnot(
  abs(R$q25_risk$estimate - R$q75_risk$estimate - R$difference$estimate) < 1e-12,
  R$q25_risk$estimate > R$q75_risk$estimate
)
stepify <- function(t, v, lo = v, hi = v) {
  n <- length(t)
  ti <- c(1, rep(2:n, each = 2))
  vi <- c(rep(1:(n - 1), each = 2), n)
  data.frame(time = t[ti], estimate = v[vi], lower = lo[vi], upper = hi[vi])
}
cu <- ao("standardized_mortality_curves")
TMAX <- 7
cu <- cu[cu$time <= TMAX + 1e-9, ]
chk5 <- cu[abs(cu$time - 5) < 1e-9, ]
stopifnot(
  abs(chk5$estimate[chk5$series == "lower_strength_q25"] - R$q25_risk$estimate) < 1e-10,
  abs(chk5$estimate[chk5$series == "higher_strength_q75"] - R$q75_risk$estimate) < 1e-10,
  abs(chk5$estimate[chk5$series == "difference"] - R$difference$estimate) < 1e-10
)
ser <- function(s) {
  z <- cu[cu$series == s, ]
  z <- z[order(z$time), ]
  st <- stepify(z$time, 100 * z$estimate, 100 * z$lower, 100 * z$upper)
  st$series <- s
  st
}
c3 <- rbind(ser("lower_strength_q25"), ser("higher_strength_q75"))
nar <- unique(cu[
  cu$series == "difference" & abs(cu$time - round(cu$time)) < 1e-9,
  c("time", "at_risk")
])
ytop <- ceiling(max(c3$estimate) + 1)
p5 <- data.frame(
  x = 5, y = 100 * c(R$q25_risk$estimate, R$q75_risk$estimate), col = c(COL_S, COL_S2),
  lab = sprintf("%s%%", f2(100 * c(R$q25_risk$estimate, R$q75_risk$estimate)))
)
ends <- data.frame(
  x = TMAX + 0.12,
  y = c(tail(c3$estimate[c3$series == "lower_strength_q25"], 1), tail(c3$estimate[c3$series == "higher_strength_q75"], 1)),
  col = c(COL_S, COL_S2),
  lab = c("Weaker grip\n(25th percentile)", "Stronger grip\n(75th percentile)")
)
p3a <- ggplot(c3, aes(time, estimate, group = series)) +
  geom_vline(xintercept = 5, linetype = "22", colour = REFC, linewidth = 0.3) +
  geom_line(aes(colour = series), linewidth = 0.9) +
  geom_point(
    data = p5,
    aes(x, y),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    colour = p5$col,
    fill = p5$col
  ) +
  geom_text(
    data = p5,
    aes(x = x + c(-0.12, 0.12), y = y, label = lab),
    inherit.aes = FALSE,
    hjust = c(1, 0),
    vjust = c(-0.6, 1.7),
    size = pt(BASE), family = FONT, fontface = "bold", colour = p5$col
  ) +
  geom_text(
    data = ends,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 0,
    size = pt(BASE - 0.5),
    family = FONT,
    colour = ends$col, lineheight = 0.9
  ) +
  scale_colour_manual(values = c(lower_strength_q25 = COL_S, higher_strength_q75 = COL_S2)) +
  scale_x_continuous(breaks = 0:TMAX, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, ytop, 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, TMAX), ylim = c(0, ytop), clip = "off") +
  annotate("text",
    x = nar$time,
    y = -ytop * 0.25,
    label = fmt_n(nar$at_risk),
    size = pt(BASE - 1),
    family = FONT,
    colour = MUTED
  ) +
  annotate("text",
    x = -0.6,
    y = -ytop * 0.25,
    label = "At risk",
    hjust = 1,
    size = pt(BASE - 1),
    family = FONT,
    colour = MUTED,
    fontface = "bold"
  ) +
  labs(
    x = "Years since examination",
    y = "Cumulative mortality (%)",
    title = "A  Standardized cumulative mortality"
  ) +
  theme_pub() +
  theme(plot.margin = margin(1.5, 22, 8, 11, "mm"))
dfc <- cu[cu$series == "difference", ]
dfc <- dfc[order(dfc$time), ]
d3 <- stepify(dfc$time, 100 * dfc$estimate, 100 * dfc$lower, 100 * dfc$upper)
d5 <- dfc[abs(dfc$time - 5) < 1e-9, ]
p3b <- ggplot(d3, aes(time, estimate)) +
  geom_hline(yintercept = 0, linetype = "22", colour = REFC, linewidth = 0.35) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = TINT_S) +
  geom_line(colour = COL_S, linewidth = 0.9) +
  annotate("point",
    x = 5,
    y = 100 * d5$estimate,
    shape = 21,
    size = 2.2,
    colour = COL_S,
    fill = COL_S
  ) +
  annotate("text",
    x = 0.25,
    y = ceiling(max(d3$upper) + 0.5) - 0.2,
    hjust = 0,
    vjust = 1,
    size = pt(BASE - 0.5),
    family = FONT,
    colour = INK,
    label = sprintf(
      "At 5 years: %s points\n(95%% CI, %s to %s)",
      f2(100 * d5$estimate),
      f2(100 * d5$lower),
      f2(100 * d5$upper)
    )
  ) +
  scale_x_continuous(breaks = 0:TMAX, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, TMAX), ylim = c(-0.5, ceiling(max(d3$upper) + 0.5)), clip = "off") +
  labs(
    x = "Years since examination",
    y = "Difference (percentage points)",
    title = "B  Weaker minus stronger grip"
  ) +
  theme_pub() +
  theme(plot.margin = margin(1.5, 2, 8, 4, "mm"))
fde <- a2("functional_disability_estimates")
fsd <- a2("functional_disability_standardized")
fsd_ <- function(o, e) {
  r <- fsd[fsd$outcome == o & fsd$estimand == e, ]
  stopifnot(nrow(r) == 1)
  r
}
fde_ <- function(o) {
  r <- fde[fde$outcome == o, ]
  stopifnot(nrow(r) == 1)
  r
}
outs <- c(adl = "Basic daily activities (ADL)", iadl = "Household activities (IADL)")
bars <- do.call(rbind, lapply(names(outs), function(o) {
  lo <- fsd_(o, "q25_prevalence")
  hi <- fsd_(o, "q75_prevalence")
  data.frame(
    outcome = factor(outs[[o]], levels = outs),
    grp = factor(c("Weaker grip (25th percentile)", "Stronger grip (75th percentile)"),
      levels = c("Weaker grip (25th percentile)", "Stronger grip (75th percentile)")
    ),
    est = 100 * c(lo$estimate, hi$estimate),
    lo = 100 * c(lo$lower, hi$lower),
    hi = 100 * c(lo$upper, hi$upper)
  )
}))
stopifnot(all(bars$est[bars$grp == levels(bars$grp)[1]] > bars$est[bars$grp == levels(bars$grp)[2]]))
orl <- data.frame(
  outcome = factor(outs, levels = outs), y = 22.6,
  lab = vapply(names(outs), function(o) {
    e <- fde_(o)
    d <- fsd_(o, "difference")
    sprintf(
      "Difference %s points (95%% CI, %s to %s)\nOR per SD lower grip %s",
      f1(100 * d$estimate),
      f1(100 * d$lower),
      f1(100 * d$upper),
      sprintf("%s (%s to %s)", f2(e$OR), f2(e$lower), f2(e$upper))
    )
  }, "")
)
dodge <- position_dodge(width = 0.72)
p3c <- ggplot(bars, aes(outcome, est, fill = grp)) +
  geom_col(position = dodge, width = 0.68, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
    position = dodge,
    width = 0.16,
    linewidth = 0.4,
    colour = INK
  ) +
  geom_text(aes(y = hi + 0.6, label = sprintf("%s%%", f1(est))),
    position = dodge,
    vjust = 0,
    size = pt(BASE - 0.5),
    family = FONT,
    colour = INK
  ) +
  geom_text(
    data = orl,
    aes(outcome, y, label = lab),
    inherit.aes = FALSE,
    vjust = 1,
    size = pt(BASE - 1),
    family = FONT,
    colour = MUTED,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = c(COL_S, COL_S2), name = NULL) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Adjusted prevalence of\nany difficulty (%)",
    title = "C  Limitation in daily activities, adults aged 60\u201379 years"
  ) +
  theme_pub() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = BASE - 0.5),
    legend.key.size = unit(3.2, "mm"), legend.margin = margin(0, 0, 0, 0, "mm")
  )
fig3 <- wrap_elements(p3a + p3b + plot_layout(widths = c(
  1.35,
  1
))) / wrap_elements(p3c) + plot_layout(heights = c(
  1,
  0.95
))
save_fig(fig3, "Figure_3", 174, 160, P$fig_main)

# FIGURE 4 (physical function within guideline PREVENT categories)
gb <- a2("guideline_burden")
gdc <- a2("guideline_discordance")
cats5 <- c("<3%", "3% to <5%", "5% to <10%", ">=10%")
catl <- c("<3%", "3% to <5%", "5% to <10%", "\u226510%")
m5 <- c(
  low_grip = "Low grip strength (age 50\u201379 y)",
  anyfx = "Limitation in mobility or daily activities (age 60\u201379 y)"
)
pb <- do.call(rbind, lapply(names(m5), function(k) {
  r <- gb[gb$measure == k & gb$category_var == "asc_cat", ]
  r <- r[match(cats5, r$category), ]
  stopifnot(nrow(r) == 4, !anyNA(r$prevalence))
  data.frame(
    cat = factor(catl, levels = catl),
    grp = m5[[k]],
    est = 100 * r$prevalence,
    lo = 100 * r$prev_lower,
    hi = 100 * r$prev_upper
  )
}))
pb$grp <- factor(pb$grp, levels = unique(pb$grp))
p4a <- ggplot(pb, aes(cat, est, fill = grp)) +
  scale_x_discrete() +
  annotate("rect", xmin = 2.5, xmax = 4.6, ymin = 0, ymax = 62, fill = NEUTRAL) +
  annotate("text",
    x = 3.55,
    y = 61,
    label = "Statin recommended (\u22655%)",
    vjust = 1,
    size = pt(BASE - 0.5),
    family = FONT,
    colour = MUTED
  ) +
  geom_col(position = position_dodge(width = 0.74), width = 0.7, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
    position = position_dodge(width = 0.74),
    width = 0.16,
    linewidth = 0.4,
    colour = INK
  ) +
  geom_text(aes(y = hi + 1, label = sprintf("%s%%", formatC(est, format = "f", digits = 0))),
    position = position_dodge(width = 0.74), vjust = 0,
    size = pt(BASE - 0.5), family = FONT, colour = INK
  ) +
  scale_fill_manual(values = c(COL_S, COL_M), name = NULL) +
  scale_y_continuous(limits = c(0, 62), breaks = seq(0, 60, 10), expand = c(0, 0)) +
  labs(
    x = "PREVENT 10-year ASCVD risk",
    y = "Prevalence (%)",
    title = "A  Reduced physical function at every level of PREVENT risk"
  ) +
  theme_pub() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = BASE - 0.5),
    legend.key.size = unit(3.2, "mm"), legend.margin = margin(0, 0, 0, 0, "mm")
  )
g5 <- gdc[gdc$analysis == "grip, ASCVD 5%", ]
g5_ <- function(e) {
  r <- g5[g5$estimand == e, ]
  stopifnot(nrow(r) == 1, r$reportable)
  r
}
pm5 <- do.call(rbind, lapply(list(
  c("below_preserved", "Below 5%", "Preserved grip"), c("below_low", "Below 5%", "Low grip"),
  c("above_preserved", "5% or higher", "Preserved grip"), c("above_low", "5% or higher", "Low grip")
), function(z) {
  r <- g5_(z[1])
  data.frame(
    band = z[2],
    grp = z[3],
    est = 100 * r$estimate,
    lo = 100 * r$lower,
    hi = 100 * r$upper,
    pop = r$pop_millions,
    d = r$deaths
  )
}))
pm5$band <- factor(pm5$band, levels = c("Below 5%", "5% or higher"))
pm5$grp <- factor(pm5$grp, levels = c("Preserved grip", "Low grip"))
dsc <- g5_("discordance")
xb <- c(1 + 0.185, 2 - 0.185)
yb <- max(
  pm5$hi[pm5$band == "Below 5%" & pm5$grp == "Low grip"],
  pm5$hi[pm5$band == "5% or higher" & pm5$grp == "Preserved grip"]
) + 5.2
p4b <- ggplot(pm5, aes(band, est, fill = grp)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.7, colour = NA) +
  geom_errorbar(aes(ymin = pmax(lo, 0), ymax = hi),
    position = position_dodge(width = 0.74),
    width = 0.16,
    linewidth = 0.4,
    colour = INK
  ) +
  geom_text(aes(y = hi + 0.35, label = sprintf("%s%%\n%s deaths\n%s million", formatC(est, format = "f", digits = 1), fmt_n(d), formatC(pop, format = "f", digits = 1))),
    position = position_dodge(width = 0.74),
    vjust = 0,
    size = pt(BASE - 1),
    family = FONT,
    colour = INK,
    lineheight = 0.9
  ) +
  annotate("segment", x = xb[1], xend = xb[2], y = yb, yend = yb, linewidth = 0.35, colour = INK) +
  annotate("segment", x = xb, xend = xb, y = yb, yend = yb - 0.6, linewidth = 0.35, colour = INK) +
  annotate("text",
    x = mean(xb), y = yb + 0.4, vjust = 0, size = pt(BASE - 0.5), family = FONT, colour = INK,
    label = sprintf(
      "Difference %s percentage points (95%% CI, %s to %s)",
      fxm(100 * dsc$estimate, 1),
      fxm(100 * dsc$lower, 1),
      fxm(100 * dsc$upper, 1)
    )
  ) +
  scale_fill_manual(values = c(COL_S2, COL_S), name = NULL) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 15, 5), expand = c(0, 0)) +
  labs(
    x = "PREVENT 10-year ASCVD risk",
    y = "5-year mortality (%)",
    title = "B  Five-year mortality by grip strength below and above the statin threshold"
  ) +
  theme_pub() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = BASE - 0.5),
    legend.key.size = unit(3.2, "mm"), legend.margin = margin(0, 0, 0, 0, "mm")
  )
stopifnot(abs(pm5$est[2] - 100 * g5_("below_low")$estimate) < 1e-10)
fig4 <- wrap_elements(p4a) / wrap_elements(p4b) + plot_layout(heights = c(1, 1))
save_fig(fig4, "Figure_4", 174, 160, P$fig_main)

# FIGURE S1 (selection: input availability and the excluded groups)
re <- res("range_exclusions")
mi0 <- res("missingness")
inp <- c(
  tc = "Total cholesterol (130\u2013320 mg/dL)",
  hdl = "HDL cholesterol (20\u2013100 mg/dL)",
  sbp = "Systolic BP (90\u2013200 mm Hg)",
  egfr = "eGFR from creatinine (15\u2013140)",
  bmi = "Body-mass index (18.5\u201339.9 kg/m\u00b2)",
  height = "Height",
  diabetes = "Diabetes status",
  smoking = "Smoking status",
  bp_tx = "Antihypertensive treatment",
  statin = "Statin use",
  grip = "Valid maximal grip in both hands",
  mobility = "Mobility status (asked at age \u226560 y)"
)
mi0$variable[mi0$variable == "scr"] <- "egfr"
s1l <- do.call(rbind, lapply(names(inp), function(v) {
  m <- mi0$missing[mi0$variable == v]
  stopifnot(length(m) == 1)
  b <- if (v %in% re$variable) re$below[re$variable == v] else 0
  a <- if (v %in% re$variable) re$above[re$variable == v] else 0
  if (v %in% re$variable) stopifnot(re$missing[re$variable == v] == m)
  data.frame(
    var = inp[[v]],
    kind = c("Missing", "Below study range", "Above study range"),
    n = c(m, b, a)
  )
}))
s1l$var <- factor(s1l$var, levels = rev(inp))
s1l$kind <- factor(s1l$kind, levels = c("Above study range", "Below study range", "Missing"))
tot1 <- aggregate(n ~ var, s1l, sum)
lab1 <- do.call(rbind, lapply(levels(s1l$var), function(v) {
  z <- s1l[s1l$var == v, ]
  data.frame(var = v, tot = sum(z$n), lab = if (sum(z$n[z$kind != "Missing"]) > 0)
    sprintf(
      "%s missing; %s below; %s above",
      fmt_n(z$n[z$kind == "Missing"]),
      fmt_n(z$n[z$kind == "Below study range"]),
      fmt_n(z$n[z$kind == "Above study range"])
    )
  else sprintf("%s missing", fmt_n(z$n[z$kind == "Missing"])))
}))
lab1$var <- factor(lab1$var, levels = levels(s1l$var))
s1a <- ggplot(s1l, aes(n, var, fill = kind)) +
  geom_col(width = 0.68, colour = NA) +
  geom_text(
    data = lab1,
    aes(x = tot + 25, y = var, label = lab),
    inherit.aes = FALSE,
    hjust = 0,
    size = pt(BASE - 1),
    family = FONT,
    colour = INK
  ) +
  scale_fill_manual(
    values = c(Missing = "#9A9A9A", "Below study range" = COL_S2, "Above study range" = COL_S),
    breaks = c("Missing", "Below study range", "Above study range"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, 1650), breaks = seq(0, 1200, 300), expand = c(0, 0)) +
  labs(
    x = "Clinical candidates, n (of 3,853)",
    y = NULL,
    title = "A  Missing and out-of-range inputs"
  ) +
  theme_pub() +
  theme(
    plot.title.position = "plot",
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = BASE - 0.5),
    legend.key.size = unit(3, "mm"), axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )
sg <- a2("selection_by_gate")
glev <- c(
  "Incomplete PREVENT inputs",
  "PREVENT inputs outside equation ranges",
  "No eligible mortality linkage or no height",
  "Grip not tested",
  "Valid maximal-effort grip from 1 hand only",
  "No qualifying maximal-effort trial",
  "All excluded candidates",
  "Included (primary cohort)"
)
stopifnot(setequal(glev, sg$group))
sg <- sg[match(glev, sg$group), ]
gshort <- c(
  "Incomplete PREVENT inputs",
  "PREVENT inputs out of range",
  "No mortality linkage or height",
  "Grip not tested",
  "Valid grip in 1 hand only", "No maximal-effort trial", "All excluded", "Included cohort"
)
GREY <- "#9A9A9A"
# B: one measure (crude mortality), one axis; the 2 summary rows sit below a gap.
s1b <- data.frame(
  y = c(9.2, 8.2, 7.2, 6.2, 5.2, 4.2, 2.6, 1.6), lab = sprintf("%s (n = %s)", gshort, fmt_n(sg$N)),
  pct = sg$crude_mortality_pct, deaths = sg$deaths, linked = sg$linked,
  col = c(rep(GREY, 6), INK, COL_S)
)
s1b$txt <- ifelse(is.na(s1b$pct), "No linked follow-up",
  sprintf(
    "%s%%  (%s of %s)",
    formatC(s1b$pct, format = "f", digits = 1),
    fmt_n(s1b$deaths),
    fmt_n(s1b$linked)
  )
)
s1bp <- ggplot(s1b, aes(y = y)) +
  geom_rect(aes(xmin = 0, xmax = ifelse(is.na(pct), 0, pct), ymin = y - 0.34, ymax = y + 0.34),
    fill = s1b$col,
    colour = NA
  ) +
  geom_text(aes(x = ifelse(is.na(pct), 0, pct) + 0.4, label = txt),
    hjust = 0, size = pt(BASE - 1), family = FONT,
    colour = ifelse(s1b$col == GREY, INK, s1b$col)
  ) +
  scale_y_continuous(breaks = s1b$y, labels = s1b$lab, expand = c(0, 0.5)) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  labs(
    x = "Died during follow-up, % (crude)",
    y = NULL,
    title = "B  Mortality in each excluded group, in the order the criteria were applied"
  ) +
  theme_pub() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.title.position = "plot")
# C: excluded vs included, one small panel per characteristic, 2 bars each.
cvar <- c(
  wt_mean_bmi = "Mean BMI, kg/m\u00b2",
  wt_diabetes_pct = "Diabetes, %",
  wt_mean_age = "Mean age, y",
  wt_female_pct = "Women, %"
)
cg <- sg[sg$group %in% c("All excluded candidates", "Included (primary cohort)"), ]
clab <- sprintf(c("All excluded (n = %s)", "Included cohort (n = %s)"), fmt_n(cg$N))
s1c <- do.call(rbind, lapply(names(cvar), function(v) data.frame(
  var = factor(cvar[[v]], levels = cvar), grp = factor(clab, levels = rev(clab)),
  val = cg[[v]], col = c(INK, COL_S)
)))
s1c$txt <- formatC(s1c$val, format = "f", digits = 1)
s1cp <- ggplot(s1c, aes(val, grp)) +
  geom_col(fill = s1c$col, width = 0.62) +
  geom_text(aes(label = txt), hjust = -0.15, size = pt(BASE - 1), family = FONT, colour = s1c$col) +
  facet_wrap(~var, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.38)), n.breaks = 3) +
  labs(x = NULL, y = NULL, title = "C  All excluded candidates compared with the included cohort") +
  theme_pub() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = BASE, hjust = 0),
    panel.spacing = unit(4, "mm"), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = BASE - 1), plot.title.position = "plot"
  )
figS1 <- wrap_elements(s1a) / wrap_elements(s1bp) / wrap_elements(s1cp) + plot_layout(heights = c(
  1.15,
  0.78,
  0.34
))
save_fig(figS1, "Figure_S1", 174, 215, P$fig_supp)

# FIGURE S2 (unadjusted cumulative mortality by grip third)
km <- ao("km_by_grip_tertile")
nr <- ao("km_numbers_at_risk")
TK <- 8
levs <- c("Lowest third", "Middle third", "Highest third")
kcol <- c("Lowest third" = COL_S, "Middle third" = COL_S2, "Highest third" = "#9A9A9A")
ks <- do.call(rbind, lapply(levs, function(k) {
  z <- km[km$tertile == k & km$time <= TK, ]
  z <- z[order(z$time), ]
  z <- rbind(
    data.frame(time = 0, cum_mortality = 0, lower = 0, upper = 0),
    z[c("time", "cum_mortality", "lower", "upper")]
  )
  z <- rbind(
    z,
    data.frame(
      time = TK,
      cum_mortality = tail(z$cum_mortality, 1),
      lower = tail(z$lower, 1),
      upper = tail(z$upper, 1)
    )
  )
  st <- stepify(z$time, 100 * z$cum_mortality, 100 * z$lower, 100 * z$upper)
  st$tertile <- k
  st
}))
ks$tertile <- factor(ks$tertile, levels = levs)
ytk <- ceiling(max(ks$estimate) + 1)
nd <- unique(nr[c("tertile", "n", "deaths")])
lab_k <- data.frame(
  tertile = factor(levs, levels = levs), x = TK + 0.15,
  y = vapply(levs, function(k) tail(ks$estimate[ks$tertile == k], 1), 0),
  label = sprintf(
    "%s of grip strength\n%s deaths / %s",
    levs,
    fmt_n(nd$deaths[match(levs, nd$tertile)]),
    fmt_n(nd$n[match(levs, nd$tertile)])
  )
)
lab_k$y <- lab_k$y + c(0, 0.3, -0.5)
gk <- ggplot(ks, aes(time, estimate, group = tertile)) +
  geom_line(aes(colour = tertile), linewidth = 0.9) +
  geom_text(
    data = lab_k,
    aes(x = x, y = y, label = label, colour = tertile),
    hjust = 0,
    size = pt(BASE - 0.5),
    family = FONT,
    inherit.aes = FALSE, lineheight = 0.9
  ) +
  scale_colour_manual(values = kcol) +
  scale_x_continuous(breaks = 0:TK, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, ytk, 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, TK), ylim = c(0, ytk), clip = "off") +
  labs(x = "Years since examination", y = "Cumulative mortality (%), unadjusted") +
  theme_pub() +
  theme(plot.margin = margin(2, 40, 27, 22, "mm"))
for (i in seq_along(levs)) {
  z <- nr[nr$tertile == levs[i] & nr$time <= TK, ]
  gk <- gk + annotate("text",
    x = z$time,
    y = -ytk * (0.19 + 0.07 * i),
    label = fmt_n(z$at_risk),
    size = pt(BASE - 1),
    family = FONT,
    colour = kcol[[levs[i]]]
  ) +
    annotate("text",
      x = -0.45,
      y = -ytk * (0.19 + 0.07 * i),
      label = levs[i],
      hjust = 1,
      size = pt(BASE - 1),
      family = FONT,
      colour = kcol[[levs[i]]]
    )
}
gk <- gk + annotate("text",
  x = -0.45,
  y = -ytk * 0.19,
  label = "At risk",
  hjust = 1,
  size = pt(BASE - 1),
  family = FONT,
  colour = MUTED,
  fontface = "bold"
)
save_fig(gk, "Figure_S2", 174, 108, P$fig_supp)

# FIGURE S3 (grip distribution and shape of the association, by sex)
cv <- res("grip_spline_curve")
nl <- res("nonlinearity_tests")
sc <- res("strength_scaling")
ad <- readRDS(file.path(P$analytic, "study4_analytic_design_sample.rds"))
qq2 <- res("risk_strength_contrasts")
sx_panels <- lapply(1:2, function(sx) {
  nm <- c("Men", "Women")[sx]
  mu <- sc$mean_kg[sc$sex == sx]
  s <- sc$sd_kg[sc$sex == sx]
  cvs <- transform(cv, kg = mu + strength_z * s)
  xr <- range(cvs$kg)
  lin <- data.frame(kg = seq(xr[1], xr[2], length.out = 100))
  lin$HR <- exp(-p_prim$beta * (lin$kg - mu) / s)
  z <- ad[ad$core & ad$sex == sx, c("grip", "w")]
  dd <- density(z$grip, weights = z$w / sum(z$w), n = 512, bw = bw.nrd0(z$grip))
  dens <- data.frame(x = dd$x, y = dd$y)
  q <- c(qq2$q25[qq2$sex == sx], qq2$q75[qq2$sex == sx])
  top <- ggplot(dens, aes(x, y)) +
    geom_area(fill = TINT_S, colour = COL_S, linewidth = 0.4) +
    geom_vline(xintercept = q, linetype = "22", colour = INK, linewidth = 0.3) +
    annotate("text",
      x = q,
      y = Inf,
      label = sprintf("%s\n%s kg", c("25th", "75th"), f1(q)),
      vjust = 1.1,
      hjust = c(1.08, -0.08),
      size = pt(BASE - 1), family = FONT, colour = INK, lineheight = 0.9
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = xr) +
    labs(
      x = NULL,
      y = if (sx == 1) "Distribution" else NULL,
      title = sprintf("%s  %s", LETTERS[sx], nm)
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank()
    )
  bot <- ggplot(cvs, aes(kg, HR)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = TINT_S) +
    geom_hline(yintercept = 1, linetype = "22", colour = REFC, linewidth = 0.35) +
    geom_line(data = lin, aes(kg, HR), linetype = "42", colour = MUTED, linewidth = 0.45) +
    geom_line(colour = COL_S, linewidth = 0.9) +
    scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = xr, ylim = c(0.12, 8)) +
    labs(
      x = "Grip strength (kg)",
      y = if (sx == 1) "Hazard ratio for death (log scale)" else NULL
    ) +
    theme_pub()
  if (sx == 2) bot <- bot + annotate("text",
    x = xr[2], y = 7.5, hjust = 1, vjust = 1, size = pt(BASE - 0.5), family = FONT, colour = INK,
    label = sprintf(
      "Overall association: %s\nNonlinearity: %s",
      pv(nl$P[nl$test == "overall"]),
      pv(nl$P[nl$test == "nonlinear"])
    )
  )
  if (sx == 1) bot <- bot + annotate("text",
    x = xr[1] + 0.03 * diff(xr), y = 0.15, hjust = 0, vjust = 0, size = pt(BASE - 1), family = FONT,
    colour = MUTED, label = "Solid: flexible spline (band, 95% CI)\nDashed: main linear model"
  )
  top / bot + plot_layout(heights = c(0.38, 1))
})
figS3 <- wrap_elements(sx_panels[[1]]) | wrap_elements(sx_panels[[2]])
save_fig(figS3, "Figure_S3", 174, 110, P$fig_supp)

# FIGURE S6 (incremental prediction: change in each metric)
mrow <- data.frame(
  cohort = c("primary", "mobility", "mobility", "mobility"),
  model = c("M1", "Mobility", "Grip", "Both"),
  label = c("Add grip strength", "Add mobility limitation", "Add grip strength", "Add both measures"),
  y = c(5, 3, 2, 1), col = c(COL_S, COL_M, COL_S, COL_B)
)
s6h <- data.frame(y = c(6, 4), label = c("Adults aged 50\u201379 y", "Adults aged 60\u201379 y"))
mets <- c(dC = "C statistic", dAUC5 = "5-year AUC", NRI = "Continuous NRI", IDI = "IDI")
YL <- c(0.45, 6.6)
s6lab <- ggplot() +
  geom_text(
    data = s6h,
    aes(x = 0, y = y, label = label),
    hjust = 0,
    size = pt(BASE),
    family = FONT,
    fontface = "bold",
    colour = INK
  ) +
  geom_text(
    data = mrow,
    aes(x = 0.06, y = y, label = label),
    hjust = 0,
    size = pt(BASE),
    family = FONT,
    colour = mrow$col
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = YL, clip = "off") +
  labs(title = " ") +
  theme_void(base_family = FONT) +
  theme(
    plot.title = element_text(size = BASE + 0.5, margin = margin(0, 0, 1.5, 0, "mm")),
    plot.margin = margin(1.5, 0, 1.5, 1.5, "mm")
  )
mnum <- function(x) gsub("-", "\u2212", format(x, drop0trailing = TRUE, trim = TRUE))
s6p <- lapply(names(mets), function(me) {
  d <- do.call(rbind, lapply(seq_len(nrow(mrow)), function(i) {
    r <- pm[pm$cohort == mrow$cohort[i] & pm$metric == me & pm$model == mrow$model[i], ]
    stopifnot(nrow(r) == 1)
    data.frame(y = mrow$y[i], col = mrow$col[i], est = r$apparent, lo = r$lower, hi = r$upper)
  }))
  ggplot(d, aes(est, y)) +
    geom_hline(yintercept = 4, colour = RULE, linewidth = 0.25) +
    geom_vline(xintercept = 0, linetype = "22", colour = REFC, linewidth = 0.35) +
    geom_errorbar(aes(xmin = lo, xmax = hi),
      width = 0.28,
      orientation = "y",
      colour = d$col,
      linewidth = 0.5
    ) +
    geom_point(shape = 21, size = 2.1, colour = d$col, fill = d$col) +
    scale_x_continuous(n.breaks = 3, labels = mnum) +
    expand_limits(x = 0) +
    coord_cartesian(ylim = YL) +
    labs(x = NULL, y = NULL, title = mets[[me]]) +
    theme_pub() +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
      plot.title = element_text(size = BASE + 0.5, face = "bold", hjust = 0.5, margin = margin(0, 0, 1.5, 0, "mm")),
      plot.margin = margin(1.5, 3, 1.5, 1, "mm")
    )
})
figS6 <- (s6lab | s6p[[1]] | s6p[[2]] | s6p[[3]] | s6p[[4]]) + plot_layout(widths = c(
  1.15,
  1,
  1,
  1,
  1
)) +
  plot_annotation(
    caption = "Change compared with the base model (95% CI)",
    theme = theme(plot.caption = element_text(
      family = FONT,
      size = BASE + 0.5,
      colour = INK,
      hjust = 0.62
    ))
  )
save_fig(figS6, "Figure_S6", 174, 62, P$fig_supp)

# FIGURE S7 (physical function within guideline PREVENT categories: all measures)
msr <- c(
  low_grip = "Low grip strength",
  mobility = "Mobility limitation",
  adl = "ADL limitation",
  iadl = "IADL limitation",
  anyfx = "Any of the 3 limitations"
)
mcol <- c(
  "Low grip strength" = "#0072B2", "Mobility limitation" = "#D55E00", "ADL limitation" = "#E69F00",
  "IADL limitation" = "#CC79A7", "Any of the 3 limitations" = "#1A1A1A"
)
prev_rows <- function(cv, lv, lab) do.call(rbind, lapply(names(msr), function(m) {
  r <- gb[gb$measure == m & gb$category_var == cv & gb$category %in% lv, ]
  r <- r[match(lv, r$category), ]
  stopifnot(nrow(r) == length(lv))
  data.frame(
    measure = factor(msr[[m]], levels = msr),
    cat = factor(lab, levels = lab),
    est = 100 * r$prevalence,
    lo = 100 * r$prev_lower,
    hi = 100 * r$prev_upper
  )
}))
pa <- prev_rows(
  "asc_cat",
  c("<3%", "3% to <5%", "5% to <10%", ">=10%"),
  c("<3%", "3% to\n<5%", "5% to\n<10%", "\u226510%")
)
pb7 <- prev_rows("tot75", c("0", "1"), c("<7.5%", "\u22657.5%"))
dg <- position_dodge(width = 0.3)
prev_panel <- function(d, title, xlab, end_labels) {
  p <- ggplot(d, aes(cat, est, colour = measure, group = measure)) +
    geom_line(position = dg, linewidth = 0.55) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = dg, width = 0, linewidth = 0.4) +
    geom_point(position = dg, size = 1.9) +
    scale_colour_manual(values = mcol, name = NULL) +
    scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 20), expand = c(0, 0)) +
    labs(x = xlab, y = "Prevalence (%)", title = title) +
    theme_pub()
  if (end_labels) {
    e <- d[d$cat == levels(d$cat)[nlevels(d$cat)], ]
    e$lab <- sprintf("%s (%s%%)", e$measure, formatC(e$est, format = "f", digits = 0))
    e <- e[order(e$est), ]
    gap <- 4.2
    y <- e$est
    for (k in 2:length(y)) y[k] <- max(y[k], y[k - 1] + gap)
    e$ylab <- y
    p <- p + geom_text(
      data = e,
      aes(x = nlevels(d$cat) + 0.3, y = ylab, label = lab),
      hjust = 0,
      size = pt(BASE - 0.5),
      family = FONT,
      show.legend = FALSE
    ) +
      coord_cartesian(clip = "off") + theme(plot.margin = margin(1.5, 46, 1.5, 1.5, "mm"))
  }
  p
}
s7a <- prev_panel(pa, "A  By PREVENT-ASCVD risk", "PREVENT 10-year ASCVD risk", TRUE)
s7b <- prev_panel(
  pb7,
  "B  By total CVD risk",
  "PREVENT 10-year\ntotal CVD risk",
  FALSE
) + theme(plot.margin = margin(1.5, 4, 1.5, 1.5, "mm"))
shv <- list(
  c("share_below_asc5", "Below the statin threshold\n(PREVENT-ASCVD <5%)"),
  c("share_below_tot75", "Below the blood-pressure threshold\n(PREVENT total CVD <7.5%)")
)
s7s <- do.call(rbind, lapply(names(msr), function(m) do.call(rbind, lapply(shv, function(sv) {
  r <- gb[gb$measure == m & gb$category_var == sv[1], ]
  stopifnot(nrow(r) == 1)
  data.frame(
    measure = factor(msr[[m]], levels = rev(msr)),
    thr = factor(sv[2], levels = vapply(shv, `[`, "", 2)),
    est = 100 * r$prevalence
  )
}))))
s7c <- ggplot(s7s, aes(est, measure, fill = measure)) +
  geom_col(width = 0.68) +
  geom_text(aes(x = est + 1.5, label = sprintf("%s%%", formatC(est, format = "f", digits = 0))),
    hjust = 0,
    size = pt(BASE - 0.5),
    family = FONT,
    colour = INK
  ) +
  facet_wrap(~thr, nrow = 1) +
  scale_fill_manual(values = mcol) +
  scale_x_continuous(limits = c(0, 55), breaks = seq(0, 50, 25), expand = c(0, 0)) +
  labs(
    x = "Share of adults with the measure (%)",
    y = NULL,
    title = "C  Adults with each measure whose PREVENT estimate is below the treatment threshold"
  ) +
  theme_pub() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = BASE - 0.5, hjust = 0),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing.x = unit(8, "mm")
  )
figS7 <- wrap_elements(s7a + s7b + plot_layout(widths = c(
  2.2,
  1
))) / wrap_elements(s7c) + plot_layout(heights = c(
  1.15,
  1
))
save_fig(figS7, "Figure_S7", 174, 165, P$fig_supp)

# FIGURE S8 (5-year mortality within the treatment thresholds: all comparisons)
gax <- a2("guideline_discordance_adjusted_hr")
an <- list(
  c("grip, ASCVD 5%", "Grip strength, ASCVD threshold 5%", "preserved grip", "low grip"),
  c("grip, total CVD 7.5%", "Grip strength, total CVD threshold 7.5%", "preserved grip", "low grip"),
  c("reported function, ASCVD 5%", "Reported function, ASCVD threshold 5%", "no limitation", "any limitation"),
  c(
    "reported function, total CVD 7.5%",
    "Reported function, total CVD threshold 7.5%",
    "no limitation",
    "any limitation"
  )
)
stopifnot(all(vapply(an, `[`, "", 1) %in% gdc$analysis))
na_val <- function(e, l, u) ifelse(is.na(e), "<10 deaths", ci_to(e, l, u, 1))
acol <- function(a) if (grepl("grip", a[1])) COL_S else COL_M
s8m <- do.call(rbind, lapply(an, function(a) rbind(
  hdr(a[2], acol(a)),
  do.call(
    rbind,
    lapply(
      list(c("below_preserved", "Below", a[3]), c("below_low", "Below", a[4]), c("above_preserved", "At or above", a[3]), c("above_low", "At or above", a[4])),
      function(e) {
        r <- gdc[gdc$analysis == a[1] & gdc$estimand == e[1], ]
        stopifnot(nrow(r) == 1)
        data.frame(
          label = sprintf("%s, %s", e[2], e[3]), header = FALSE,
          est = if (r$reportable) 100 * r$estimate else NA,
          lo = if (r$reportable) 100 * r$lower else NA,
          hi = if (r$reportable) 100 * r$upper else NA,
          col = acol(a), n_txt = sprintf("%s / %s", fmt_n(r$deaths), fmt_n(r$N))
        )
      }
    )
  )
)))
s8a <- forest(s8m, c(0, 16), seq(0, 15, 5), "5-year all-cause mortality (%)",
  ref = 0, log = FALSE,
  title = "A  Five-year mortality in each group", val_head = "% (95% CI)", val_fun = na_val,
  n_head = "Deaths / n", widths = c(0.36, 0.16, 0.28, 0.20)
)
dlab <- c(
  rd_below = "Below threshold",
  rd_above = "At or above threshold",
  additive_interaction = "Difference between the 2",
  discordance = "Below and reduced vs above and preserved"
)
s8f <- do.call(rbind, lapply(an, function(a) {
  rows <- do.call(rbind, lapply(names(dlab), function(e) {
    r <- gdc[gdc$analysis == a[1] & gdc$estimand == e, ]
    stopifnot(nrow(r) == 1)
    if (!r$reportable) return(NULL)
    data.frame(
      label = dlab[[e]],
      header = FALSE,
      est = 100 * r$estimate,
      lo = 100 * r$lower,
      hi = 100 * r$upper,
      col = acol(a),
      n_txt = ""
    )
  }))
  rbind(hdr(a[2], acol(a)), rows)
}))
s8f$n_txt <- NULL
s8b <- forest(s8f, c(-8, 13), seq(-5, 10, 5), "Difference in 5-year mortality (percentage points)",
  ref = 0, log = FALSE,
  title = "B  Absolute differences, reduced minus preserved function", val_head = "Points (95% CI)",
  val_fun = function(e, l, u) ci_to(e, l, u, 1), widths = c(0.40, 0.01, 0.37, 0.22)
)
s8h <- do.call(
  rbind,
  lapply(
    list(c("asc5", "ASCVD threshold 5%"), c("tot75", "Total CVD threshold 7.5%")),
    function(t) {
      b <- gax[gax$threshold == t[1] & gax$term == "low_grip_below", ]
      u <- gax[gax$threshold == t[1] & gax$term == "low_grip_above", ]
      i <- gax[gax$threshold == t[1] & gax$term == "interaction", ]
      stopifnot(nrow(b) == 1, nrow(u) == 1, nrow(i) == 1)
      rbind(
        hdr(sprintf("%s (product term %s)", t[2], pv(i$P)), COL_S),
        data.frame(label = "Below threshold", header = FALSE, est = b$HR, lo = b$lower, hi = b$upper, col = COL_S, n_txt = ""),
        data.frame(
          label = "At or above threshold",
          header = FALSE,
          est = u$HR,
          lo = u$lower,
          hi = u$upper,
          col = COL_S,
          n_txt = ""
        )
      )
    }
  )
)
s8h$n_txt <- NULL
s8c <- forest(s8h,
  c(0.5, 25),
  c(0.5, 1, 2, 5, 10, 20),
  "Adjusted hazard ratio, low vs preserved grip (log scale)",
  title = "C  Relative scale", widths = c(0.40, 0.01, 0.37, 0.22)
)
figS8 <- wrap_elements(s8a) / wrap_elements(s8b) / wrap_elements(s8c) + plot_layout(heights = c(
  nrow(s8m) + 2,
  nrow(s8f) + 2.4,
  nrow(s8h) * 1.5 + 3.4
))
save_fig(figS8, "Figure_S8", 174, 225, P$fig_supp)

# FIGURE S4 (sensitivity analyses)
pd <- res("post_diagnostic_estimates")
rv <- res("additional_sensitivity_estimates")
mi <- res("imputation_result")
hs <- a2("prior_hospitalization_sensitivity")
pick <- function(tab, model, term = "lower_z") {
  r <- tab[tab$model == model & tab$term == term, ]
  stopifnot(nrow(r) == 1)
  r
}
alt3 <- a2("alternative_adjustment_sensitivity") # one row per adjustment basis, model, and term
alt3$model <- paste(alt3$adjustment, alt3$model, sep = "|")
alt3_ <- function(a, m, t = "lower_z") pick(alt3, paste(a, m, sep = "|"), t)
dA <- rbind(
  hdr("Main analysis"),
  row_("All adults aged 50\u201379 years", p_prim, COL_S),
  hdr("Follow-up and time scale"),
  row_("Age as the time scale", pick(fe, "Sensitivity: attained age"), COL_S),
  row_("First 2 years of follow-up excluded", pick(fe, "Sensitivity: two-year landmark"), COL_S),
  row_("Baseline risk stratified by PREVENT quartile", pick(pd, "Post-diagnostic: risk-stratified core"), COL_S),
  hdr("Grip measurement"),
  row_("One valid hand allowed", pick(fe, "Sensitivity: any valid hand"), COL_S),
  row_("Three maximal trials per hand required", pick(rv, "Three maximal-effort trials per hand"), COL_S),
  row_("No recent hand pain or stiffness", pick(fe, "Sensitivity: no reported hand pain"), COL_S),
  hdr("Adjustment"),
  row_("Diabetes including HbA1c \u22656.5%", pick(fe, "Sensitivity: HbA1c-inclusive diabetes"), COL_S),
  row_("Risk factors instead of PREVENT", pick(fe, "Sensitivity: individual clinical factors"), COL_S),
  row_("Risk factors with flexible age", pick(rv, "Clinical factors with matched flexible age"), COL_S),
  row_("Full PREVENT (with UACR and HbA1c)", pick(rv, "Full PREVENT with UACR and HbA1c"), COL_S),
  row_("+ education, income, cancer history", pick(rv, "Additional education income and cancer"), COL_S),
  row_("+ same, missing income imputed", mi, COL_S),
  row_("+ prior-year hospitalization", pick(hs, "Primary: adjusted for prior-year hospitalization"), COL_S),
  hdr("PREVENT version used for adjustment"),
  row_("PREVENT without its kidney term", alt3_("Kidney-neutral PREVENT (eGFR set to 90)", "Primary"), COL_S),
  row_("+ eGFR modeled separately", alt3_("Kidney-neutral PREVENT plus eGFR spline", "Primary"), COL_S),
  row_("PREVENT-ASCVD instead of total CVD", alt3_("PREVENT-ASCVD in place of total CVD", "Primary"), COL_S),
  row_("PREVENT heart failure instead of total CVD", alt3_("PREVENT-HF in place of total CVD", "Primary"), COL_S),
  hdr("Sample"),
  row_("No BMI range restriction", pick(rv, "BMI range relaxed"), COL_S),
  row_(
    "Hospitalized in prior year excluded",
    pick(hs, "Primary: excluding prior-year hospitalization"),
    COL_S
  )
)
panSA <- forest(dA,
  c(0.8, 3.2),
  c(1, 1.5, 2, 3),
  "Hazard ratio per SD lower grip strength (log scale)",
  title = "A  Grip strength, all adults aged 50\u201379 years", band = c(p_prim$lower, p_prim$upper)
)
pdm <- function(m, t) pick(pd, m, t)
dB <- rbind(
  hdr("Lower grip strength (per SD)", COL_S),
  row_("Main analysis", p_ds, COL_S),
  row_("Baseline risk stratified by PREVENT quartile", pdm("Post-diagnostic: risk-stratified mobility_cohort", "lower_z"), COL_S),
  row_("Spline knots set in older adults", pdm("Post-diagnostic: mobility-specific spline knots", "lower_z"), COL_S),
  row_("Hospitalized in prior year excluded", pick(hs, "Mobility D: excluding prior-year hospitalization"), COL_S),
  row_("PREVENT without its kidney term", alt3_("Kidney-neutral PREVENT (eGFR set to 90)", "Mobility D"), COL_S),
  row_("PREVENT-ASCVD instead of total CVD", alt3_("PREVENT-ASCVD in place of total CVD", "Mobility D"), COL_S),
  hdr("Mobility limitation (yes vs no)", COL_M),
  row_("Main analysis", p_dm, COL_M),
  row_("Baseline risk stratified by PREVENT quartile", pdm("Post-diagnostic: risk-stratified mobility_cohort", "mobility"), COL_M),
  row_("Spline knots set in older adults", pdm("Post-diagnostic: mobility-specific spline knots", "mobility"), COL_M),
  row_("Hospitalized in prior year excluded", pick(hs, "Mobility D: excluding prior-year hospitalization", "mobility"), COL_M),
  row_("PREVENT without its kidney term", alt3_("Kidney-neutral PREVENT (eGFR set to 90)", "Mobility D", "mobility"), COL_M),
  row_(
    "PREVENT-ASCVD instead of total CVD",
    alt3_("PREVENT-ASCVD in place of total CVD", "Mobility D", "mobility"),
    COL_M
  )
)
panSB <- forest(dB, c(0.8, 5), c(1, 1.5, 2, 3, 4), "Hazard ratio (log scale)",
  title = "B  Both measures together, adults aged 60\u201379 years"
)
figS7 <- wrap_elements(panSA) / wrap_elements(panSB) + plot_layout(heights = c(
  nrow(dA) + 2.5,
  nrow(dB) + 2.5
))
save_fig(figS7, "Figure_S4", 174, 240, P$fig_supp)

# FIGURE S9 (proportional-hazards diagnostics)
zz <- readRDS(file.path(P$derived_add, "ph_joint.rds"))
phj <- res("PH_joint_terms")
ph_panel <- function(coh, var, title, col, beta, term_key) {
  z <- zz[[coh]]
  sm <- survival:::plot.cox.zph(z,
    var = var,
    resid = FALSE,
    se = TRUE,
    df = 4,
    nsmo = 80,
    plot = FALSE
  )
  keep <- !duplicated(z$x)
  tt <- approx(z$x[keep], z$time[keep], sm$x, rule = 2)$y
  dfp <- data.frame(t = tt, y = sm$y[, 1], up = sm$y[, 2], lo = sm$y[, 3])
  pval <- phj$p[phj$cohort == coh & phj$term == term_key]
  ggplot(dfp, aes(t, y)) +
    geom_ribbon(aes(ymin = lo, ymax = up), fill = if (col == COL_M) TINT_M else TINT_S) +
    geom_hline(yintercept = 0, linetype = "22", colour = REFC, linewidth = 0.3) +
    geom_hline(yintercept = beta, linetype = "42", colour = MUTED, linewidth = 0.4) +
    geom_line(colour = col, linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 9), breaks = seq(0, 8, 2), expand = c(0, 0)) +
    labs(
      x = "Years since examination", y = "Log hazard ratio over time", title = title,
      subtitle = sprintf("Test for change over time: %s", pv(pval))
    ) +
    theme_pub()
}
s4a <- ph_panel("core", 1, "A  Grip strength, age 50\u201379 y", COL_S, p_prim$beta, "lower_z")
s4b <- ph_panel(
  "mobility_cohort",
  1,
  "B  Grip strength, age 60\u201379 y",
  COL_S,
  p_ds$beta,
  "lower_z"
) + labs(y = NULL)
s4c <- ph_panel(
  "mobility_cohort",
  2,
  "C  Mobility limitation",
  COL_M,
  p_dm$beta,
  "mobility"
) + labs(y = NULL)
figS8 <- s4a | s4b | s4c
save_fig(figS8, "Figure_S9", 174, 68, P$fig_supp)

# FIGURE S5 (5-year mortality difference under alternative models)
ars <- res("absolute_risk_sensitivities")
stt <- res("support_target_comparison")
mlab <- c(
  original = "Main model",
  risk_stratified = "Baseline risk stratified by PREVENT quartile",
  nonlinear_grip = "Flexible spline for grip strength"
)
tlab <- c(
  original_support = sprintf("Supported population (n = %s)", fmt_n(stt$N[stt$target == "original_support"])),
  full_target = sprintf("All participants (n = %s)", fmt_n(stt$N[stt$target == "full_target"]))
)
d5r <- do.call(
  rbind,
  lapply(
    names(mlab),
    function(m) rbind(hdr(mlab[[m]]), do.call(rbind, lapply(names(tlab), function(t) {
      e <- ars[ars$model == m & ars$target == t & ars$estimand == "difference", ]
      stopifnot(nrow(e) == 1)
      data.frame(
        label = tlab[[t]],
        header = FALSE,
        est = 100 * e$estimate,
        lo = 100 * e$lower,
        hi = 100 * e$upper,
        col = COL_S,
        n_txt = ""
      )
    })))
  )
)
stopifnot(abs(d5r$est[2] - 100 * R$difference$estimate) < 1e-10)
d5r$n_txt <- NULL
figS9 <- forest(d5r,
  c(-0.3, 3.5),
  0:3,
  "Difference in 5-year mortality, weaker minus stronger grip (percentage points)",
  ref = 0, log = FALSE,
  val_head = "Difference (95% CI)", widths = c(0.42, 0.01, 0.37, 0.20)
)
save_fig(figS9, "Figure_S5", 174, 66, P$fig_supp)

writeLines(capture.output(sessionInfo()), file.path(P$logs, "sessionInfo_figures.txt"))
