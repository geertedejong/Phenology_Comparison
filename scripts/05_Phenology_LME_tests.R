# Function to run linear mixed-effects model and generate boxplots with stats included
lmm_boxplot <- function(df, phase_id, species = NULL, title = "") {
  # Filter data by phenophase and year range
  filtered_data <- df %>% 
    filter(phase_ID == phase_id, Year >= 2016 & Year <= 2019) 
  
  # Further filter by species if provided
  if (!is.null(species)) {
    filtered_data <- filtered_data %>% filter(Spp == species)
  }
  
  # Check if 'obs' and 'phase_ID' have at least two levels
  if (nlevels(as.factor(filtered_data$obs)) < 2) {
    message("Skipping analysis for phase ", phase_id, " and species ", species, 
            " due to insufficient levels in 'obs'")
    return(NULL)
  }
  
  # Fit the linear mixed model (LMM)
  lmm_result <- lme4::lmer(phase_DATE ~ phase_ID * obs + (1|Year), data = filtered_data)
  lmm_summary <- broom.mixed::tidy(lmm_result)
  
  # Extract p-value for interaction term
  p_value <- lmm_summary %>% 
    filter(term == "phase_ID:obs") %>%
    pull(p.value)
  
  annotation_y <- max(filtered_data$phase_DATE, na.rm = TRUE) + 10
  
  # Generate boxplot with annotations
  plot <- filtered_data %>%
    ggplot(aes(x = obs, y = phase_DATE, fill = obs, col = obs)) +
    geom_boxplot(alpha = 0.8, outlier.colour = NA) + # Add transparency to the boxes
    geom_jitter(width = 0.2, size = 2, alpha = 0.2) + 
    hrbrthemes::scale_fill_ipsum() +
    hrbrthemes::scale_colour_ipsum() +
    labs(x = "Observation type", y = "DOY (2016 - 2019)", 
         title = title, fill = "Observation type") +
    theme_classic() +
    theme(legend.position = "none") +
    annotate("text", x = Inf, y = annotation_y, 
             label = paste("p =", format.pval(p_value)), 
             hjust = 1.1, vjust = 1.5, size = 4, color = "black", fontface = "italic")
  
  return(list(plot = plot, summary = lmm_summary))
}

# List of all phases
phases <- list(
  list(phase_id = "P1", species = NULL, title = "First Day 100% Snow Free"),
  list(phase_id = "P2", species = "ERIVAG", title = "First E. vaginatum Bud Appearance"),
  list(phase_id = "P2", species = "DRYINT", title = "First D. integrifolia Bud Appearance"),
  list(phase_id = "P3", species = "DRYINT", title = "First D. integrifolia Open Flower"),
  list(phase_id = "P4", species = "DRYINT", title = "First D. integrifolia Petal Shed"),
  list(phase_id = "P5", species = "DRYINT", title = "First D. integrifolia Twisting of Filament"),
  list(phase_id = "P2", species = "SALARC", title = "S. arctica First Leaf Bud Burst"),
  list(phase_id = "P5", species = "SALARC", title = "S. arctica First Leaf Turns Yellow"),
  list(phase_id = "P6", species = "SALARC", title = "S. arctica Last Leaf Turns Yellow")
)

# Use the function to run analysis for each phenophase and save plots/summaries
results <- map(phases, function(phases) {
  lmm_boxplot(pheno, phases$phase_id, phases$species, phases$title)
})

# Collect the plots and summaries
plots <- map(results, "plot")
lmm_summaries <- map_dfr(results, "summary")

# Combine and save the plots
lmm_pheno <- ggarrange(plotlist = plots, ncol = 3, nrow = 3)

ggsave(lmm_pheno, filename = "figures/lmm_phenocam_box_2024.png", height = 10, width = 12)

# Save the LMM summary table to a CSV file
write.csv(lmm_summaries, file = "figures/lmm_summary.csv", row.names = FALSE)

