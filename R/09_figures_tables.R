CONTROL_COLOUR <- "#67A9CF"
MDD_COLOUR <- "#EF8A62"
CONTRAST_COLOUR <- "#7B4AB5"
CHARCOAL <- "#333333"

paper_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.tag = element_text(face = "bold", size = base_size + 4),
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_figure <- function(plot, stem, width, height) {
  figure_dir <- file.path(OUTPUT_DIR, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf, bg = "white")
}

slab_rows <- function(draws, label, y, colour, height = 0.34) {
  den <- density(draws, n = 512)
  tibble(x = den$x, y = y, density = den$y / max(den$y) * height,
         label = label, colour = colour)
}

slab_plot <- function(data, labels, x_label, zero = FALSE) {
  ggplot(data, aes(x = x, y = y, group = label, fill = colour, colour = colour)) +
    geom_ribbon(aes(ymin = y - density, ymax = y + density), alpha = 0.25,
                linewidth = 0) +
    geom_line(aes(y = y + density), linewidth = 0.55) +
    geom_segment(data = data %>% group_by(label, y, colour) %>%
                   summarise(x = quantile(x, 0.025), xend = quantile(x, 0.975),
                             .groups = "drop"),
                 aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.8) +
    geom_vline(xintercept = if (zero) 0 else NULL, colour = "grey70",
               linewidth = 0.5) +
    scale_colour_identity() + scale_fill_identity() +
    scale_y_continuous(breaks = seq_along(labels), labels = labels,
                       expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = x_label, y = NULL) + paper_theme()
}

build_figure2 <- function(inputs, ep_primary) {
  d <- ep_primary$prepared$trials %>%
    left_join(ep_primary$prepared$participants %>%
                select(participant_id, diagnosis), by = "participant_id")
  subject <- d %>%
    group_by(participant_id, diagnosis, required_effort_pct) %>%
    summarise(rating = mean(rating), .groups = "drop")
  group <- subject %>%
    group_by(diagnosis, required_effort_pct) %>%
    summarise(mean = mean(rating), se = sd(rating) / sqrt(n()), .groups = "drop")
  panel_a <- ggplot() +
    geom_line(data = subject,
              aes(required_effort_pct, rating, group = participant_id, colour = diagnosis),
              alpha = 0.10, linewidth = 0.35) +
    geom_ribbon(data = group,
                aes(required_effort_pct, ymin = mean - se, ymax = mean + se, fill = diagnosis),
                alpha = 0.18, colour = NA) +
    geom_line(data = group,
              aes(required_effort_pct, mean, colour = diagnosis), linewidth = 1.2) +
    geom_point(data = group,
               aes(required_effort_pct, mean, colour = diagnosis), size = 2.2) +
    scale_colour_manual(values = c(control = CONTROL_COLOUR, MDD = MDD_COLOUR),
                        labels = c(control = "Controls", MDD = "MDD")) +
    scale_fill_manual(values = c(control = CONTROL_COLOUR, MDD = MDD_COLOUR)) +
    scale_x_continuous(breaks = c(10, 30, 50, 70)) +
    labs(x = "Objective effort (%MVC)", y = "Subjective effort rating") +
    paper_theme() + theme(legend.position = "top")

  f <- ep_primary$summary$finite_draws
  group_slabs <- bind_rows(
    slab_rows(f$control_intercept, "Controls intercept", 4, CONTROL_COLOUR),
    slab_rows(f$MDD_intercept, "MDD intercept", 3, MDD_COLOUR),
    slab_rows(f$control_discriminability, "Controls discriminability", 2, CONTROL_COLOUR),
    slab_rows(f$MDD_discriminability, "MDD discriminability", 1, MDD_COLOUR)
  )
  panel_b <- slab_plot(
    group_slabs,
    c("MDD discriminability", "Controls discriminability", "MDD intercept", "Controls intercept"),
    "Posterior standardised estimate"
  )

  contrast_slabs <- bind_rows(
    slab_rows(f$MDD_minus_control_intercept, "Intercept", 2, CONTRAST_COLOUR),
    slab_rows(f$MDD_minus_control_discriminability, "Discriminability", 1, CONTRAST_COLOUR)
  )
  panel_c <- slab_plot(contrast_slabs, c("Discriminability", "Intercept"),
                       "Posterior MDD - Controls contrast", zero = TRUE)
  figure <- panel_a + (panel_b / panel_c) +
    patchwork::plot_layout(widths = c(1.55, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_figure(figure, "Figure2_primary_results", 11.5, 6.3)
  figure
}

build_figure3 <- function(inputs, ep_primary, dm, associations) {
  ep_means <- ep_primary$summary$participant_means
  dm_means <- dm$summary$participant_means
  clinical <- ep_means %>%
    left_join(inputs$participants %>% select(participant_id, phq9_total),
              by = "participant_id")
  cross <- ep_means %>% inner_join(dm_means, by = c("participant_id", "diagnosis", "source_group"))
  panel_a <- ggplot(cross %>% filter(diagnosis == "MDD"),
                    aes(discriminability_posterior_mean,
                        effort_sensitivity_posterior_mean)) +
    geom_point(colour = MDD_COLOUR, alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, colour = MDD_COLOUR, fill = MDD_COLOUR,
                alpha = 0.16, linewidth = 0.9) +
    labs(x = "Effort discriminability", y = "Decision-making effort sensitivity") +
    paper_theme()
  clinical_long <- clinical %>%
    select(participant_id, diagnosis, phq9_total,
           `Rating intercept` = intercept_posterior_mean,
           `Effort discriminability` = discriminability_posterior_mean) %>%
    pivot_longer(c(`Rating intercept`, `Effort discriminability`),
                 names_to = "measure", values_to = "estimate")
  panel_b <- ggplot(clinical_long %>% filter(diagnosis == "MDD"),
                    aes(phq9_total, estimate, colour = measure)) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    scale_colour_manual(values = c("Effort discriminability" = MDD_COLOUR,
                                   "Rating intercept" = "grey35")) +
    labs(x = "PHQ-9", y = "Posterior participant estimate") +
    paper_theme() + theme(legend.position = "top")
  post <- associations$cross_task %>% filter(diagnosis == "MDD") %>%
    mutate(term = factor(term, levels = rev(c(
      "EP_discriminability_vs_DM_effort", "EP_intercept_vs_DM_effort",
      "discriminability_minus_intercept_specificity"
    ))))
  panel_c <- ggplot(post, aes(mean, term)) +
    geom_vline(xintercept = 0, colour = "grey75") +
    geom_errorbarh(aes(xmin = q025, xmax = q975, colour = term), height = 0) +
    geom_point(aes(colour = term), size = 2.5) +
    scale_colour_manual(values = c(MDD_COLOUR, MDD_COLOUR, CHARCOAL), guide = "none") +
    scale_y_discrete(labels = c(
      "discriminability_minus_intercept_specificity" = "Difference",
      "EP_intercept_vs_DM_effort" = "Rating intercept",
      "EP_discriminability_vs_DM_effort" = "Effort discriminability"
    )) +
    labs(x = "Posterior Spearman association", y = NULL) + paper_theme()
  figure <- panel_a + panel_b + panel_c + patchwork::plot_annotation(tag_levels = "A")
  save_figure(figure, "Figure3_clinical_cross_task", 12, 4.2)
  figure
}

build_supplementary_figures <- function(ordinal_piecewise, observer, dm) {
  ordinal <- ordinal_piecewise$ordinal_subject$all
  p_s2a <- ggplot(ordinal, aes(diagnosis, ordinal_accuracy, colour = diagnosis)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, alpha = 0.35) +
    scale_colour_manual(values = c(control = CONTROL_COLOUR, MDD = MDD_COLOUR),
                        guide = "none") +
    labs(x = NULL, y = "Ordinal accuracy") + paper_theme()
  p_s2b <- ggplot(ordinal_piecewise$subject_means,
                  aes(required_effort_pct, rating, group = participant_id,
                      colour = diagnosis)) +
    geom_line(alpha = 0.10) +
    stat_summary(aes(group = diagnosis), fun = mean, geom = "line", linewidth = 1.1) +
    scale_colour_manual(values = c(control = CONTROL_COLOUR, MDD = MDD_COLOUR)) +
    scale_x_continuous(breaks = c(10, 30, 50, 70)) +
    labs(x = "Objective effort (%MVC)", y = "Subjective effort rating") +
    paper_theme() + theme(legend.position = "top")
  save_figure(p_s2a + p_s2b + patchwork::plot_annotation(tag_levels = "A"),
              "Figure_S2_ordinal_piecewise", 9.5, 4.2)

  obs <- observer$summary$summary
  group_terms <- obs %>%
    filter(term %in% c("control_prior_uncertainty", "MDD_prior_uncertainty",
                       "control_sensory_uncertainty", "MDD_sensory_uncertainty"))
  p_s3a <- ggplot(group_terms, aes(mean, reorder(term, mean), colour = term)) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = c(CONTROL_COLOUR, MDD_COLOUR,
                                   CONTROL_COLOUR, MDD_COLOUR), guide = "none") +
    labs(x = "Posterior uncertainty estimate", y = NULL) + paper_theme()
  selected <- obs %>%
    filter(term %in% c("MDD_minus_control_prior_uncertainty",
                       "MDD_minus_control_sensory_uncertainty",
                       "MDD_minus_control_effective_prior_weight"))
  p_s3b <- ggplot(selected, aes(mean, reorder(term, mean))) +
    geom_vline(xintercept = 0, colour = "grey75") +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0,
                   colour = CONTRAST_COLOUR) +
    geom_point(colour = CONTRAST_COLOUR, size = 2.4) +
    labs(x = "Posterior MDD - Controls contrast", y = NULL) + paper_theme()
  save_figure(p_s3a + p_s3b + patchwork::plot_annotation(tag_levels = "A"),
              "Figure_S3_dynamic_observer", 10.5, 4.2)

  dm_p <- dm$summary$participant_means %>%
    pivot_longer(c(reward_sensitivity_posterior_mean,
                   effort_sensitivity_posterior_mean),
                 names_to = "parameter", values_to = "value")
  p_s4 <- ggplot(dm_p, aes(diagnosis, value, colour = diagnosis)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, alpha = 0.35) +
    facet_wrap(~ parameter, scales = "free_y") +
    scale_colour_manual(values = c(control = CONTROL_COLOUR, MDD = MDD_COLOUR),
                        guide = "none") +
    labs(x = NULL, y = "Posterior participant estimate") + paper_theme()
  save_figure(p_s4, "Figure_S4_decision_making", 8.5, 4.2)

  if (!is.null(dm$recovery)) {
    recovery <- dm$recovery$participant
    p_a <- ggplot(recovery, aes(true_reward_sensitivity,
                                recovered_reward_sensitivity, colour = diagnosis)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey45") + geom_point(alpha = 0.6) +
      scale_x_log10() + scale_y_log10() +
      scale_colour_manual(values = c(control = CONTROL_COLOUR, MDD = MDD_COLOUR)) +
      labs(x = "Generating reward sensitivity", y = "Recovered reward sensitivity") +
      paper_theme() + theme(legend.position = "none")
    p_b <- ggplot(recovery, aes(true_effort_sensitivity,
                                recovered_effort_sensitivity, colour = diagnosis)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey45") + geom_point(alpha = 0.6) +
      scale_x_log10() + scale_y_log10() +
      scale_colour_manual(values = c(control = CONTROL_COLOUR, MDD = MDD_COLOUR)) +
      labs(x = "Generating effort sensitivity", y = "Recovered effort sensitivity") +
      paper_theme() + theme(legend.position = "none")
    ppc <- tibble(option2_rate = dm$summary$ppc_draws)
    p_c <- ggplot(ppc, aes(option2_rate)) +
      geom_density(fill = "grey70", colour = CHARCOAL, alpha = 0.45) +
      geom_vline(xintercept = dm$summary$ppc$observed_option2_rate,
                 colour = MDD_COLOUR) +
      labs(x = "Model-coded option 2 choice rate", y = "Density") +
      paper_theme() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    save_figure(p_a + p_b + p_c + patchwork::plot_annotation(tag_levels = "A"),
                "Figure_S5_recovery_PPC", 12, 4.1)
  }
}

write_prior_table <- function() {
  priors <- tribble(
    ~analysis, ~parameter, ~prior, ~constraint,
    "Primary EP", "Group intercept means", "Normal(0, 1)", "unbounded",
    "Primary EP", "Group effort-slope means", "Normal(0.75, 0.20)", "unbounded",
    "Primary EP", "Random-intercept SD", "Exponential(1)", "> 0",
    "Primary EP", "Random-slope SD", "Exponential(8)", "> 0",
    "Primary EP", "Random-effect correlation", "LKJ(2)", "correlation matrix",
    "Observer", "Initial prior mean", "Normal(50, 20)", "unbounded",
    "Observer", "Prior uncertainty", "log Normal(log(20), 0.75)", "1 to 80 after exp transform",
    "Observer", "Sensory uncertainty", "log Normal(log(20), 0.75)", "1 to 80 after exp transform",
    "Observer", "Learning rate", "Normal(logit(0.25), 0.8)", "0 to 1 after inverse-logit",
    "Early force", "Required-effort coefficient", "Normal(100, 40)", "unbounded",
    "Early force", "Early-estimate coefficient", "Normal(0, 40)", "unbounded",
    "Early force", "Early-estimate by effort", "Normal(0, 80)", "unbounded",
    "Decision making", "Group log reward sensitivity", "Normal(0, 1)", "unbounded; exponentiated",
    "Decision making", "Group log effort sensitivity", "Normal(0, 1)", "unbounded; exponentiated",
    "Decision making", "Between-participant SDs", "Exponential(1)", "> 0"
  )
  dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
  write_csv(priors, file.path(OUTPUT_DIR, "tables", "supplementary_priors_constraints.csv"))
}

build_manuscript_outputs <- function(inputs, table1, ep_primary, ep_sensitivity,
                                     ordinal_piecewise, observer, effort_prior,
                                     dm, associations) {
  figure1 <- file.path(REPO_ROOT, "assets", "Figure1_intro_schematic.png")
  if (file.exists(figure1)) {
    dir.create(file.path(OUTPUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
    file.copy(figure1, file.path(OUTPUT_DIR, "figures", "Figure1_intro_schematic.png"),
              overwrite = TRUE)
  }
  build_figure2(inputs, ep_primary)
  build_figure3(inputs, ep_primary, dm, associations)
  build_supplementary_figures(ordinal_piecewise, observer, dm)
  write_prior_table()
  invisible(TRUE)
}
