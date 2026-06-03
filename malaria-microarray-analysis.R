###############################################################################
# Malaria exposure and heterologous antibody responses
#
# Purpose:
#   Process protein microarray data from the Junju and Ngerenya childhood cohorts,
#   derive antigen-specific antibody summaries, generate manuscript figures, and
#   fit age-adjusted mixed-effects models comparing antibody responses by cohort
#   and malaria episode burden.
#
# Notes for reuse:
#   1. This script preserves the analytical logic of the original working script,
#      but removes exploratory clutter, duplicate package calls, hard-coded
#      setwd() calls, and repeated plotting blocks.
#   2. Update the paths in the USER SETTINGS section before running.
#   3. The script assumes GenePix-style output files with foreground-background
#      columns named `F635 Mean - B635` and `F532 Mean - B532`.
###############################################################################

## ============================================================================
## 0. Packages
## ============================================================================

required_packages <- c(
  "broom.mixed",
  "data.table",
  "dplyr",
  "ggbeeswarm",
  "ggplot2",
  "ggpubr",
  "hrbrthemes",
  "lme4",
  "lmerTest",
  "lubridate",
  "purrr",
  "readxl",
  "reshape2",
  "splines",
  "stringr",
  "tibble",
  "tidyr",
  "wesanderson",
  "openxlsx"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

## ============================================================================
## 1. User settings
## ============================================================================

# Replace these paths with paths on your system or, preferably, with paths inside
# a project directory such as data/raw, data/metadata, results/tables, and
# results/figures.

base_dir <- "Malaria immunology microarray/Analysis"

paths <- list(
  junju_template      = file.path(base_dir, "Malaria cohort5/Slide template/JUNJU_IMMUNOLOGY_SLIDES.xlsx"),
  ngerenya_template  = file.path(base_dir, "Malaria cohort5/Slide template/NGERENYA_IMMUNOLOGY_SLIDES.xlsx"),
  junju_results_dir  = file.path(base_dir, "Malaria cohort5/Results_txt/Junju"),
  ngerenya_results_dir = file.path(base_dir, "Malaria cohort5/Results_txt/Ngerenya"),
  misprints          = file.path(base_dir, "misprints.xlsx"),
  date_of_births     = file.path(base_dir, "Malaria cohort/malaria immunology cohort dobs.xlsx"),
  vaccination_data   = "Measles paper/mal.imm.cohort.full.vacccination.data.xlsx",
  output_tables_dir  = "results/tables",
  output_figures_dir = "results/figures"
)

# Create output directories if they do not already exist.
dir.create(paths$output_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$output_figures_dir, recursive = TRUE, showWarnings = FALSE)

# Main antigens used in the cross-sectional heterologous antibody analysis.
target_antigens <- c(
  "BPERT_V2",
  "COXSACB1",
  "CMVPP150",
  "EBV",
  "HSV1GDP",
  "MEASL_V2",
  "RUBV_V2",
  "H1N1HACALIF2009"
)

# Antigen used for single-antigen longitudinal plots.
primary_antigen <- "PFALAMA1_V2"

# Participants excluded because they had too few Ngerenya observations.
low_observation_ids <- c("N1060", "N1040", "N1027", "N1028")

## ============================================================================
## 2. Helper functions
## ============================================================================

read_all_excel_sheets <- function(file_path) {
  sheets <- readxl::excel_sheets(file_path)
  sheet_data <- lapply(sheets, function(sheet) {
    as.data.frame(readxl::read_excel(file_path, sheet = sheet))
  })
  names(sheet_data) <- sheets
  sheet_data
}

extract_slide_sample_ids <- function(template_sheet) {
  # The original template stores 24 miniarrays as 8 columns after the first
  # three header rows. The ordering below preserves the original logic.
  as.character(unlist(template_sheet[-1:-3, 1:8], use.names = FALSE))
}

clean_signal <- function(x) {
  ifelse(is.na(x), NA_real_, ifelse(x < 0, 0, x))
}

process_cohort_microarrays <- function(template_file, results_dir, cohort_name) {
  templates <- read_all_excel_sheets(template_file)
  sheets <- readxl::excel_sheets(template_file)
  raw_files <- list.files(results_dir, pattern = "\\.txt$", full.names = FALSE)

  # Junju raw files contain a known spelling error in some filenames.
  searchable_files <- gsub("Sliide", "Slide", raw_files)

  slide_data <- vector("list", length(sheets))
  names(slide_data) <- sheets

  for (sheet in sheets) {
    slide_name <- gsub("\\s+", "", sheet)
    file_index <- grep(
      paste0(tolower(slide_name), "_"),
      tolower(searchable_files)
    )

    if (length(file_index) == 0) {
      warning("No raw file found for ", cohort_name, ", ", slide_name)
      next
    }

    if (length(file_index) > 1) {
      warning("Multiple raw files found for ", cohort_name, ", ", slide_name, "; using the first match")
      file_index <- file_index[1]
    }

    raw_file <- file.path(results_dir, raw_files[file_index])
    sample_ids <- tibble(
      miniarray = 1:24,
      SampleID = extract_slide_sample_ids(templates[[sheet]])
    )

    dat <- data.table::fread(raw_file, skip = "Flags") %>%
      as_tibble() %>%
      mutate(miniarray = ceiling(Block / 1)) %>%
      left_join(sample_ids, by = "miniarray") %>%
      mutate(
        personID = stringr::str_split_fixed(SampleID, "_", 2)[, 1],
        year = as.numeric(stringr::str_split_fixed(SampleID, "_", 2)[, 2]),
        slide = as.numeric(gsub("slide", "", tolower(slide_name))),
        `F635 Mean - B635` = clean_signal(`F635 Mean - B635`),
        `F532 Mean - B532` = clean_signal(`F532 Mean - B532`),
        measurement = `F635 Mean - B635`,
        measurement_iga = `F532 Mean - B532`,
        cohort = cohort_name
      )

    slide_data[[sheet]] <- dat
  }

  full_data <- bind_rows(slide_data)

  summary_data <- full_data %>%
    filter(personID != "PC") %>%
    group_by(personID, year, Name, slide, cohort) %>%
    summarise(
      means = median(measurement, na.rm = TRUE),
      means_iga = median(measurement_iga, na.rm = TRUE),
      .groups = "drop"
    )

  list(raw = full_data, summary = summary_data)
}

remove_ngerenya_misprints <- function(raw_data, summary_data, misprint_file) {
  misprints <- readxl::read_xlsx(misprint_file) %>%
    reshape2::melt(id.vars = "slide") %>%
    as_tibble() %>%
    select(slide, Name = value) %>%
    filter(!is.na(Name)) %>%
    group_by(slide, Name) %>%
    mutate(n_printed = n()) %>%
    filter(n_printed == 1) %>%
    ungroup() %>%
    select(slide, Name) %>%
    mutate(misprint_flag = TRUE)

  raw_clean <- raw_data %>%
    left_join(misprints, by = c("slide", "Name")) %>%
    filter(is.na(misprint_flag)) %>%
    select(-misprint_flag)

  summary_clean <- summary_data %>%
    left_join(misprints, by = c("slide", "Name")) %>%
    filter(is.na(misprint_flag)) %>%
    select(-misprint_flag)

  list(raw = raw_clean, summary = summary_clean)
}

assign_cohort_from_id <- function(person_id) {
  case_when(
    grepl("n", tolower(person_id)) ~ "ngerenya",
    grepl("c", tolower(person_id)) ~ NA_character_,
    TRUE ~ "junju"
  )
}

add_age_metadata <- function(raw_data, summary_data, dob_file) {
  dob_data <- readxl::read_xlsx(dob_file) %>%
    distinct(studyno, .keep_all = TRUE) %>%
    mutate(personID = studyno) %>%
    select(-1:-3)

  summary_age <- summary_data %>%
    left_join(dob_data, by = "personID") %>%
    mutate(age.y = year - lubridate::year(dob))

  raw_age <- raw_data %>%
    left_join(dob_data, by = "personID") %>%
    mutate(age.y = year - lubridate::year(dob))

  list(raw = raw_age, summary = summary_age)
}

identify_slide_spillovers <- function(antigen_raw_data) {
  antigen_raw_data %>%
    filter(!personID %in% c("NC", "PC")) %>%
    group_by(personID, slide) %>%
    summarise(n_records = n(), .groups = "drop") %>%
    group_by(personID) %>%
    mutate(n_slides = n()) %>%
    filter(n_slides != 1) %>%
    group_by(personID) %>%
    mutate(spillover_flag = n_records == min(n_records)) %>%
    filter(spillover_flag) %>%
    ungroup() %>%
    select(personID, slide, spillover_flag)
}

remove_slide_spillovers <- function(data, spillovers) {
  data %>%
    left_join(spillovers, by = c("personID", "slide")) %>%
    filter(is.na(spillover_flag)) %>%
    select(-spillover_flag)
}

make_antigen_labels <- function(data) {
  data %>%
    mutate(
      antigen_label = case_when(
        Name == "BPERT_V2"        ~ "B. pertussis",
        Name == "COXSACB1"        ~ "Coxsackievirus B1",
        Name == "CMVPP150"        ~ "Cytomegalovirus",
        Name == "EBV"             ~ "Epstein-Barr virus",
        Name == "HSV1GDP"         ~ "Herpes simplex virus 1",
        Name == "MEASL_V2"        ~ "Measles",
        Name == "RUBV_V2"         ~ "Rubella",
        Name == "H1N1HACALIF2009" ~ "H1N1 influenza, 2009 pandemic strain",
        Name == "PFALAMA1_V2"     ~ "P. falciparum AMA1",
        TRUE                      ~ Name
      )
    )
}

summarise_mean_ci <- function(data, group_vars, value_col) {
  value_col <- rlang::ensym(value_col)

  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      mean_value = mean(!!value_col, na.rm = TRUE),
      median_value = median(!!value_col, na.rm = TRUE),
      sd_value = sd(!!value_col, na.rm = TRUE),
      n_value = sum(!is.na(!!value_col)),
      .groups = "drop"
    ) %>%
    mutate(
      se_value = sd_value / sqrt(n_value),
      lower_ci = mean_value - qt(1 - 0.05 / 2, n_value - 1) * se_value,
      upper_ci = mean_value + qt(1 - 0.05 / 2, n_value - 1) * se_value,
      lower_ci = pmax(lower_ci, 0),
      upper_ci = pmin(upper_ci, 65000)
    )
}

paper_theme <- function() {
  hrbrthemes::theme_ipsum() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

## ============================================================================
## 3. Load and clean microarray data
## ============================================================================

junju <- process_cohort_microarrays(
  template_file = paths$junju_template,
  results_dir = paths$junju_results_dir,
  cohort_name = "junju"
)

ngerenya <- process_cohort_microarrays(
  template_file = paths$ngerenya_template,
  results_dir = paths$ngerenya_results_dir,
  cohort_name = "ngerenya"
)

ngerenya <- remove_ngerenya_misprints(
  raw_data = ngerenya$raw,
  summary_data = ngerenya$summary,
  misprint_file = paths$misprints
)

raw_data <- bind_rows(junju$raw, ngerenya$raw) %>%
  mutate(cohort = assign_cohort_from_id(personID))

summary_data <- bind_rows(junju$summary, ngerenya$summary) %>%
  mutate(cohort = assign_cohort_from_id(personID))

# Remove sparse late Ngerenya years and children with fewer than four samples.
raw_data <- raw_data %>%
  filter(!(cohort == "ngerenya" & year %in% c(2016, 2017))) %>%
  filter(!personID %in% low_observation_ids) %>%
  mutate(cohort = factor(cohort, levels = c("ngerenya", "junju")))

summary_data <- summary_data %>%
  filter(!(cohort == "ngerenya" & year %in% c(2016, 2017))) %>%
  filter(!personID %in% low_observation_ids) %>%
  mutate(cohort = factor(cohort, levels = c("ngerenya", "junju")))

age_added <- add_age_metadata(raw_data, summary_data, paths$date_of_births)
raw_data <- age_added$raw
summary_data <- age_added$summary

openxlsx::write.xlsx(summary_data, file.path(paths$output_tables_dir, "microarray_antibody_summary_data.xlsx"))

## ============================================================================
## 4. Cohort sampling structure plot
## ============================================================================

cohort_plot_data <- summary_data %>%
  filter(Name == "ASTROVT1", !is.na(year)) %>%
  group_by(personID) %>%
  mutate(first_year = min(year, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(first_year, personID) %>%
  mutate(personID = factor(personID, levels = unique(personID)))

p_cohort_structure <- ggplot(
  cohort_plot_data,
  aes(x = year, y = personID, group = personID)
) +
  geom_line(linewidth = 0.2) +
  geom_point(size = 0.8) +
  facet_wrap(~cohort, scales = "free") +
  scale_x_continuous(breaks = sort(unique(cohort_plot_data$year))) +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "Sampling year", y = "Participant")

print(p_cohort_structure)
ggsave(
  file.path(paths$output_figures_dir, "cohort_sampling_structure.svg"),
  p_cohort_structure,
  width = 18,
  height = 4
)

## ============================================================================
## 5. Single-antigen longitudinal analysis
## ============================================================================

antigen_raw <- raw_data %>%
  filter(Name == primary_antigen, age.y >= 0) %>%
  mutate(
    measurement = clean_signal(measurement),
    measurement_iga = clean_signal(measurement_iga)
  )

antigen_summary <- summary_data %>%
  filter(Name == primary_antigen, age.y >= 0)

spillovers <- identify_slide_spillovers(antigen_raw)
antigen_raw <- remove_slide_spillovers(antigen_raw, spillovers)
antigen_summary <- remove_slide_spillovers(antigen_summary, spillovers)

age_summary <- summarise_mean_ci(
  antigen_summary,
  group_vars = c("cohort", "age.y", "Name"),
  value_col = means
)

p_primary_age <- ggplot(antigen_raw, aes(x = age.y, y = measurement)) +
  geom_quasirandom(alpha = 0.3, size = 0.3) +
  geom_point(
    data = age_summary,
    aes(y = mean_value, colour = cohort),
    size = 1.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    data = filter(age_summary, age.y > 0, age.y < 17),
    aes(y = mean_value, ymin = lower_ci, ymax = upper_ci, colour = cohort),
    width = 0.25,
    linewidth = 0.4,
    position = position_dodge(width = 0.5)
  ) +
  geom_smooth(
    data = age_summary,
    aes(y = mean_value, colour = cohort),
    method = "loess",
    se = FALSE,
    span = 0.4,
    linewidth = 0.5,
    linetype = "2121"
  ) +
  facet_wrap(~Name) +
  scale_x_continuous(breaks = sort(unique(antigen_raw$age.y))) +
  paper_theme() +
  labs(x = "Age, years", y = "Median fluorescence intensity", colour = "Cohort")

print(p_primary_age)
ggsave(
  file.path(paths$output_figures_dir, paste0(primary_antigen, "_age_trajectory.svg")),
  p_primary_age,
  width = 8,
  height = 6
)

## ============================================================================
## 6. Cross-sectional heterologous antibody comparison at age 10
## ============================================================================

heterologous_age10 <- summary_data %>%
  filter(Name %in% target_antigens, age.y == 10, !is.na(personID)) %>%
  make_antigen_labels() %>%
  remove_slide_spillovers(spillovers) %>%
  mutate(cohort = as.character(cohort))

p_age10 <- ggplot(
  heterologous_age10,
  aes(x = cohort, y = means, fill = cohort)
) +
  geom_quasirandom(alpha = 0.7, size = 0.9) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_y_continuous(trans = "log2") +
  ggpubr::stat_compare_means(
    comparisons = list(c("junju", "ngerenya")),
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    colour = "red",
    vjust = 0.65
  ) +
  facet_wrap(~antigen_label, nrow = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    strip.text.x = element_text(size = 10, angle = 90),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Antibody response, log2 MFI")

print(p_age10)
ggsave(
  file.path(paths$output_figures_dir, "age10_heterologous_antibody_comparison.svg"),
  p_age10,
  width = 10,
  height = 9
)

## ============================================================================
## 7. Age-adjusted mixed-effects models by antigen
## ============================================================================

model_data <- summary_data %>%
  filter(Name %in% target_antigens) %>%
  filter(!is.na(personID), !is.na(age.y), !is.na(means)) %>%
  filter(cohort %in% c("junju", "ngerenya")) %>%
  filter(is.finite(means), means > 0) %>%
  make_antigen_labels() %>%
  mutate(
    cohort = factor(as.character(cohort), levels = c("junju", "ngerenya")),
    age.y = as.numeric(age.y),
    log2_response = log2(means)
  )

fit_antigen_cohort_model <- function(df) {
  antigen_name <- unique(df$Name)
  antigen_label <- unique(df$antigen_label)

  if (nrow(df) < 20 || length(unique(df$personID)) < 10) {
    return(tibble(
      antigen = antigen_name,
      antigen_label = antigen_label,
      n_obs = nrow(df),
      n_ids = length(unique(df$personID)),
      beta = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      fold_change = NA_real_,
      note = "Insufficient data"
    ))
  }

  fit <- lmer(
    log2_response ~ splines::ns(age.y, df = 4) + cohort + (1 | personID),
    data = df,
    REML = FALSE
  )

  coef_tbl <- broom.mixed::tidy(
    fit,
    effects = "fixed",
    conf.int = TRUE,
    conf.method = "Wald"
  ) %>%
    filter(term == "cohortngerenya")

  tibble(
    antigen = antigen_name,
    antigen_label = antigen_label,
    n_obs = nrow(df),
    n_ids = length(unique(df$personID)),
    beta = coef_tbl$estimate[1],
    conf_low = coef_tbl$conf.low[1],
    conf_high = coef_tbl$conf.high[1],
    p_value = coef_tbl$p.value[1],
    fold_change = 2 ^ coef_tbl$estimate[1],
    note = NA_character_
  )
}

cohort_model_results <- model_data %>%
  group_split(Name) %>%
  map_dfr(fit_antigen_cohort_model) %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)

print(cohort_model_results)
write.csv(
  cohort_model_results,
  file.path(paths$output_tables_dir, "age_adjusted_cohort_model_results.csv"),
  row.names = FALSE
)

p_cohort_forest <- cohort_model_results %>%
  filter(!is.na(beta)) %>%
  mutate(antigen_label = factor(antigen_label, levels = rev(antigen_label))) %>%
  ggplot(aes(x = beta, y = antigen_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2) +
  geom_point(size = 2.5) +
  theme_classic() +
  labs(
    x = "Adjusted Ngerenya vs Junju difference, beta on log2 scale",
    y = NULL
  )

print(p_cohort_forest)
ggsave(
  file.path(paths$output_figures_dir, "age_adjusted_cohort_forest_plot.svg"),
  p_cohort_forest,
  width = 7,
  height = 4.5
)

## ============================================================================
## 8. Malaria episode burden and heterologous antibody responses
## ============================================================================

# This section requires an object or file containing malaria episode records.
# Expected minimum columns: studyno and pfpos.
# Optional age columns support early/late timing analyses.
# Load pf.ep before running this section, for example:
# pf.ep <- readRDS("data/processed/malaria_episode_records.rds")

if (exists("pf.ep")) {

  pf_episode_counts <- pf.ep %>%
    as_tibble() %>%
    filter(!is.na(studyno)) %>%
    mutate(personID = studyno) %>%
    group_by(personID) %>%
    summarise(
      malaria_episode_count = sum(pfpos == 1, na.rm = TRUE),
      .groups = "drop"
    )

  heterologous_data <- summary_data %>%
    filter(Name %in% target_antigens) %>%
    filter(!is.na(personID), !is.na(age.y), !is.na(means), means > 0) %>%
    make_antigen_labels() %>%
    transmute(
      personID,
      age.y,
      year,
      antigen = Name,
      antigen_label,
      heterologous_log2 = log2(means)
    ) %>%
    left_join(pf_episode_counts, by = "personID") %>%
    filter(!is.na(malaria_episode_count)) %>%
    group_by(antigen) %>%
    mutate(heterologous_z = as.numeric(scale(heterologous_log2))) %>%
    ungroup()

  model_main <- lmer(
    heterologous_z ~ malaria_episode_count + splines::ns(age.y, df = 4) +
      (1 | personID) + (1 | antigen),
    data = heterologous_data,
    REML = FALSE
  )

  main_term <- broom.mixed::tidy(
    model_main,
    effects = "fixed",
    conf.int = TRUE,
    conf.method = "Wald"
  ) %>%
    filter(term == "malaria_episode_count")

  print(summary(model_main))
  print(main_term)

  model_interaction <- lmer(
    heterologous_z ~ malaria_episode_count * splines::ns(age.y, df = 4) +
      (1 | personID) + (1 | antigen),
    data = heterologous_data,
    REML = FALSE
  )

  print(anova(model_main, model_interaction))

  p_episode_scatter <- ggplot(
    heterologous_data,
    aes(x = malaria_episode_count, y = heterologous_z)
  ) +
    geom_quasirandom(pch = 21, size = 0.8, stroke = 0.2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    theme_classic() +
    labs(
      x = "Febrile malaria episode count",
      y = "Heterologous antibody response, z score"
    )

  print(p_episode_scatter)
  ggsave(
    file.path(paths$output_figures_dir, "episode_count_vs_heterologous_response.svg"),
    p_episode_scatter,
    width = 5.5,
    height = 4.5
  )

  fit_single_antigen_episode_model <- function(df_antigen) {
    n_obs <- nrow(df_antigen)
    n_ids <- length(unique(df_antigen$personID))
    n_ages <- length(unique(df_antigen$age.y))

    if (n_obs < 20 || n_ids < 10 || n_ages < 4) {
      return(tibble(
        antigen = unique(df_antigen$antigen),
        antigen_label = unique(df_antigen$antigen_label),
        n_obs = n_obs,
        n_ids = n_ids,
        beta = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        note = "Insufficient data"
      ))
    }

    fit <- lmer(
      heterologous_log2 ~ malaria_episode_count + splines::ns(age.y, df = 4) +
        (1 | personID),
      data = df_antigen,
      REML = FALSE
    )

    coef_tbl <- broom.mixed::tidy(
      fit,
      effects = "fixed",
      conf.int = TRUE,
      conf.method = "Wald"
    ) %>%
      filter(term == "malaria_episode_count")

    tibble(
      antigen = unique(df_antigen$antigen),
      antigen_label = unique(df_antigen$antigen_label),
      n_obs = n_obs,
      n_ids = n_ids,
      beta = coef_tbl$estimate[1],
      conf_low = coef_tbl$conf.low[1],
      conf_high = coef_tbl$conf.high[1],
      p_value = coef_tbl$p.value[1],
      note = NA_character_
    )
  }

  single_antigen_episode_results <- heterologous_data %>%
    group_split(antigen) %>%
    map_dfr(fit_single_antigen_episode_model) %>%
    mutate(p_fdr = p.adjust(p_value, method = "fdr")) %>%
    arrange(p_value)

  print(single_antigen_episode_results)
  write.csv(
    single_antigen_episode_results,
    file.path(paths$output_tables_dir, "episode_count_single_antigen_results.csv"),
    row.names = FALSE
  )

  p_episode_forest <- single_antigen_episode_results %>%
    filter(!is.na(beta)) %>%
    mutate(antigen_label = factor(antigen_label, levels = rev(antigen_label))) %>%
    ggplot(aes(x = beta, y = antigen_label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2) +
    geom_point(size = 2.5) +
    theme_classic() +
    labs(
      x = "Association with malaria episode count, beta on log2 scale",
      y = NULL
    )

  print(p_episode_forest)
  ggsave(
    file.path(paths$output_figures_dir, "episode_count_single_antigen_forest_plot.svg"),
    p_episode_forest,
    width = 7,
    height = 4.5
  )
}

###############################################################################
# End of script
###############################################################################
