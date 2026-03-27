# =============================================================================
# Immortal Time Bias: Three Separate Figures
# Author: Ólafur B Davíðsson, PhD — Frameshift
# Date: 2026-03-26
#
# Reference: Brøndum-Jacobsen et al., Int J Epidemiol, 42(5):1486–1496, 2013
# DOI: 10.1093/ije/dyt168
# =============================================================================
#
# Produces three figures:
#   figures/fig1_nmsc.png       — NMSC study design, biased (ages 60–80)
#   figures/fig2_prize.png      — Prize analogy, biased (ages 70–80)
#   figures/fig3_correct.png    — Correct analysis from prize date
#
# Mortality: Gompertz model calibrated to Danish adult life tables
#   h(t) = a * exp(b * t), t = years since age 40
#   a = 0.00195, b = 0.0685
#   → age 40: 0.002/yr, age 70: 0.015/yr, age 80: 0.030/yr
# =============================================================================

# --- Packages ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(survival)
library(broom)

# --- Parameters --------------------------------------------------------------
set.seed(42)

n_pairs    <- 20000   # matched case-control pairs
max_age    <- 90      # administrative censoring age

gompertz_a <- 0.00195
gompertz_b <- 0.0685

# NMSC: typically diagnosed ages 60–80
nmsc_dx_min <- 60
nmsc_dx_max <- 80

# Prize: awarded ages 70–80 (shifted to maximise immortal time window)
prize_min   <- 70
prize_max   <- 80

# --- Shared colour palette (colorblind-friendly) ----------------------------
col_exposed   <- "#E69F00"   # orange — cases / prize winners
col_unexposed <- "#56B4E9"   # blue   — controls
col_immortal  <- "#CCCCCC"   # grey   — immortal time shading

# --- Helper: simulate age at death via Gompertz inverse CDF -----------------
sim_death_age <- function(n, a = gompertz_a, b = gompertz_b, entry_age = 40) {
  u <- runif(n)
  t <- (1 / b) * log(1 - (b / a) * log(1 - u))
  entry_age + pmax(t, 0)
}

# --- Shared ggplot theme -----------------------------------------------------
theme_fs <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 11, color = "grey40"),
      plot.caption     = element_text(size = 9,  color = "grey55", hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 10)
    )
}

# =============================================================================
# Simulate cohort
# =============================================================================

death_age_case    <- sim_death_age(n_pairs)
death_age_control <- sim_death_age(n_pairs)

obs_age_case    <- pmin(death_age_case,    max_age)
obs_age_control <- pmin(death_age_control, max_age)

event_case    <- as.integer(death_age_case    <= max_age)
event_control <- as.integer(death_age_control <= max_age)

# NMSC diagnosis age: uniform [60, 80]
nmsc_dx_age <- runif(n_pairs, nmsc_dx_min, nmsc_dx_max)
got_nmsc    <- as.integer(nmsc_dx_age <= obs_age_case)

# Prize age: uniform [70, 80]
prize_age   <- runif(n_pairs, prize_min, prize_max)
got_prize   <- as.integer(prize_age <= obs_age_case)

cat(sprintf("Cohort: %s matched pairs\n", format(n_pairs, big.mark = ",")))
cat(sprintf("NMSC diagnoses received: %s (%.1f%%)\n",
            sum(got_nmsc), 100 * mean(got_nmsc)))
cat(sprintf("Prizes received:         %s (%.1f%%)\n\n",
            sum(got_prize), 100 * mean(got_prize)))

# =============================================================================
# Figure 1: NMSC study design (biased)
# =============================================================================

nmsc_cases <- data.frame(
  time  = obs_age_case[got_nmsc == 1] - 40,
  event = event_case[got_nmsc == 1],
  group = "NMSC"
)
nmsc_controls <- data.frame(
  time  = obs_age_control[got_nmsc == 1] - 40,
  event = event_control[got_nmsc == 1],
  group = "No NMSC"
)
nmsc_data <- bind_rows(nmsc_cases, nmsc_controls) |>
  mutate(group = factor(group, levels = c("No NMSC", "NMSC")))

nmsc_fit <- survfit(Surv(time, event) ~ group, data = nmsc_data)
nmsc_cox <- coxph(Surv(time, event) ~ group, data = nmsc_data)
nmsc_hr  <- exp(coef(nmsc_cox))
nmsc_lo  <- exp(confint(nmsc_cox))[1]
nmsc_hi  <- exp(confint(nmsc_cox))[2]

cat(sprintf("=== Figure 1: NMSC (biased) ===\n"))
cat(sprintf("  HR = %.3f  [95%% CI: %.3f–%.3f]  (paper: 0.52)\n\n",
            nmsc_hr, nmsc_lo, nmsc_hi))

nmsc_km        <- tidy(nmsc_fit)
immortal_nmsc  <- mean(c(nmsc_dx_min, nmsc_dx_max)) - 40  # ~30 yrs

p1 <- ggplot(nmsc_km, aes(x = time, y = estimate, color = strata)) +
  annotate("rect",
           xmin = 0, xmax = immortal_nmsc,
           ymin = -Inf, ymax = Inf,
           fill = col_immortal, alpha = 0.30) +
  geom_step(linewidth = 1.0) +
  scale_color_manual(
    values = c("group=No NMSC" = col_unexposed, "group=NMSC" = col_exposed),
    labels = c("group=No NMSC" = "No NMSC", "group=NMSC" = "NMSC"),
    name   = NULL
  ) +
  scale_x_continuous(limits = c(0, 50),
                     breaks = seq(0, 50, 10),
                     labels = function(x) x + 40,
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0.30, 1.0),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Skin cancer patients appear to live longer",
    subtitle = "NMSC cases vs. matched controls, follow-up from age 40",
    x        = "Age (years)",
    y        = "Survival probability",
    caption  = NULL
  ) +
  theme_fs() +
  theme(legend.position = "bottom",
        legend.key = element_blank())

dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig1_nmsc.png", p1, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("figures/fig1_nmsc.pdf", p1, width = 8, height = 5, bg = "white")
cat("Saved figures/fig1_nmsc.png + .pdf\n\n")

# =============================================================================
# Figure 2: Prize analogy (biased)
# =============================================================================

prize_cases <- data.frame(
  time  = obs_age_case[got_prize == 1] - 40,
  event = event_case[got_prize == 1],
  group = "Prize winner"
)
prize_controls <- data.frame(
  time  = obs_age_control[got_prize == 1] - 40,
  event = event_control[got_prize == 1],
  group = "No prize"
)
prize_data <- bind_rows(prize_cases, prize_controls) |>
  mutate(group = factor(group, levels = c("No prize", "Prize winner")))

prize_fit <- survfit(Surv(time, event) ~ group, data = prize_data)
prize_cox <- coxph(Surv(time, event) ~ group, data = prize_data)
prize_hr  <- exp(coef(prize_cox))
prize_lo  <- exp(confint(prize_cox))[1]
prize_hi  <- exp(confint(prize_cox))[2]

cat(sprintf("=== Figure 2: Prize (biased) ===\n"))
cat(sprintf("  HR = %.3f  [95%% CI: %.3f–%.3f]  (no biology)\n\n",
            prize_hr, prize_lo, prize_hi))

prize_km       <- tidy(prize_fit)
immortal_prize <- mean(c(prize_min, prize_max)) - 40  # ~35 yrs

p2 <- ggplot(prize_km, aes(x = time, y = estimate, color = strata)) +
  annotate("rect",
           xmin = 0, xmax = immortal_prize,
           ymin = -Inf, ymax = Inf,
           fill = col_immortal, alpha = 0.30) +
  geom_step(linewidth = 1.0) +
  scale_color_manual(
    values = c("group=No prize" = col_unexposed, "group=Prize winner" = col_exposed),
    labels = c("group=No prize" = "No prize", "group=Prize winner" = "Prize winner"),
    name   = NULL
  ) +
  scale_x_continuous(limits = c(0, 50),
                     breaks = seq(0, 50, 10),
                     labels = function(x) x + 40,
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0.30, 1.0),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Lottery winners also appear to live longer",
    subtitle = "Prize awarded at random (ages 70-80), follow-up from age 40",
    x        = "Age (years)",
    y        = "Survival probability",
    caption  = NULL
  ) +
  theme_fs() +
  theme(legend.position = "bottom",
        legend.key = element_blank())

ggsave("figures/fig2_prize.png", p2, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("figures/fig2_prize.pdf", p2, width = 8, height = 5, bg = "white")
cat("Saved figures/fig2_prize.png + .pdf\n\n")

# =============================================================================
# Figure 3: Correct analysis — follow-up from prize date
# =============================================================================

idx <- which(got_prize == 1)

# Cases: follow-up from prize age
correct_cases <- data.frame(
  pair  = idx,
  time  = obs_age_case[idx] - prize_age[idx],
  event = event_case[idx],
  group = "Prize winner"
)

# Controls: pseudo-prize at the same age as their matched case
# Exclude controls who died before their pseudo-prize date
alive_at_pseudo <- death_age_control[idx] > prize_age[idx]

correct_controls <- data.frame(
  pair  = idx[alive_at_pseudo],
  time  = obs_age_control[idx[alive_at_pseudo]] - prize_age[idx[alive_at_pseudo]],
  event = event_control[idx[alive_at_pseudo]],
  group = "No prize"
)
correct_cases <- correct_cases[alive_at_pseudo, ]

correct_data <- bind_rows(correct_cases, correct_controls) |>
  mutate(group = factor(group, levels = c("No prize", "Prize winner")))

correct_fit <- survfit(Surv(time, event) ~ group, data = correct_data)
correct_cox <- coxph(Surv(time, event) ~ group, data = correct_data)
correct_hr  <- exp(coef(correct_cox))
correct_lo  <- exp(confint(correct_cox))[1]
correct_hi  <- exp(confint(correct_cox))[2]
correct_p   <- summary(correct_cox)$coefficients[, "Pr(>|z|)"]

cat(sprintf("=== Figure 3: Correct analysis ===\n"))
cat(sprintf("  HR = %.3f  [95%% CI: %.3f–%.3f]  p = %.3f\n\n",
            correct_hr, correct_lo, correct_hi, correct_p))

correct_km  <- tidy(correct_fit)
max_follow  <- max_age - mean(c(prize_min, prize_max))  # ~15 yrs

p3 <- ggplot(correct_km, aes(x = time, y = estimate, color = strata)) +
  geom_step(linewidth = 1.0) +
  scale_color_manual(
    values = c("group=No prize" = col_unexposed, "group=Prize winner" = col_exposed),
    labels = c("group=No prize" = "No prize", "group=Prize winner" = "Prize winner"),
    name   = NULL
  ) +
  scale_x_continuous(limits = c(0, max_follow),
                     breaks = seq(0, max_follow, 5),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0.30, 1.0),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Starting follow-up at the prize date: no effect",
    subtitle = "Both groups begin at prize age, controls assigned pseudo-prize date",
    x        = "Years since prize",
    y        = "Survival probability",
    caption  = NULL
  ) +
  theme_fs() +
  theme(legend.position = "bottom",
        legend.key = element_blank())

ggsave("figures/fig3_correct.png", p3, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("figures/fig3_correct.pdf", p3, width = 8, height = 5, bg = "white")
cat("Saved figures/fig3_correct.png + .pdf\n\n")

# =============================================================================
# Summary
# =============================================================================
cat("=== Summary ===\n")
cat(sprintf("  NMSC study (biased):      HR = %.3f  [paper: 0.52]\n", nmsc_hr))
cat(sprintf("  Prize analogy (biased):   HR = %.3f  [no biology]\n",  prize_hr))
cat(sprintf("  Correct analysis:         HR = %.3f  [true = 1.0]\n",  correct_hr))
