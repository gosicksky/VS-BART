# Plot of the performance of 3 useful and 3 noise of Z
library(ggplot2)
library(dplyr)
library(readr)

df <- read_csv("~/Desktop/final1_with_alpha.csv")
df_filtered <- df %>%
  filter(alpha == 0.050) %>% filter(ntree == 20) %>%
  mutate(simulation = as.character(simulation),
         alpha = as.character(alpha))

simulation_mapping <- data.frame(
  simulation = c("sim1", "sim2", "sim3", "sim4"),
  num_noise_Z = c(3, 4, 5, 6)  # Number of noise Z for each scenario
)

df_summary <- df_filtered %>%
  group_by(simulation, alpha) %>%
  summarise(across(
    ends_with("_mean"),
    list(mean = ~ mean(.x, na.rm = TRUE), 
         sd = ~ ifelse(n() > 1, sd(.x, na.rm = TRUE), 0)),  # Compute SD from mean values
    .names = "{.col}_{.fn}"
  ))%>%
  left_join(simulation_mapping, by = "simulation") %>%
  arrange(num_noise_Z)

print(df_summary)
write.csv(df_summary, "summary.csv", row.names = FALSE)

df <- df %>%
  filter(alpha == 0.050) %>% filter(ntree == 20) %>%
  filter(simulation == "sim1") %>%
  mutate(simulation = as.character(simulation),
         alpha = as.character(alpha))
df_long <- df %>%
  select(precision_mi_mean, recall_mi_mean, f1_mi_mean) %>%
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value")

df_long$Metric <- factor(df_long$Metric,
                         levels = c("precision_mi_mean_mean", "recall_mi_mean_mean", "f1_mi_mean_mean"),
                         labels = c("Precision", "Recall", "F1"))

ggplot(df_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot() +
  labs(title = "Boxplots of Precision, Recall, and F1", y = "Value") +
  theme_minimal(base_size = 16)

df_plot <- data.frame(
  num_noise_Z = rep(df_summary$num_noise_Z, 3),
  Value = c(df_summary$precision_mi_mean_mean,
            df_summary$recall_mi_mean_mean,
            df_summary$f1_mi_mean_mean),
  Metric = factor(rep(c("Precision", "Recall", "F1-score"),
                      each = nrow(df_summary)))
)


ggplot(df_plot, aes(x = num_noise_Z, y = Value, group = Metric)) +
  geom_line(aes(linetype = Metric), color = "black", size = 1) +
  geom_point(aes(shape = Metric), color = "black", size = 3) +
  scale_linetype_manual(values = c("Precision" = "solid",
                                   "Recall" = "dashed",
                                   "F1-score" = "dotted")) +
  scale_shape_manual(values = c("Precision" = 16,
                                "Recall" = 17,
                                "F1-score" = 15)) +
  labs(
    x = "Number of Noise Cluster-Level Covariates (Z)",
    y = "Performance Metrics"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = c(0.95, 0.45),
    legend.justification = c("right", "bottom"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length.y = unit(-0.25, "cm"),        
    axis.text.y = element_text(color = "black"),       
    axis.ticks.y = element_line(color = "black"),      
    axis.line = element_line(size = 0.8, colour = "black")
  )


ggsave("performance_plot.eps", width = 8, height = 6, dpi = 400)

# Plot if Z1, Z2 and Z3 counts
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

df <- read_csv("final_Z123_counts_alpha_0.05.csv")
df$simulation <- factor(
  paste0("sim", seq_len(nrow(df))),
  levels = paste0("sim", 1:4),
  labels = c("3U,3N", "3U,4N", "3U,5N", "3U,6N")
)

df_long <- df %>%
  pivot_longer(
    cols = -simulation,
    names_to = c("Variable", "Method"),
    names_sep = "_",
    values_to = "Count"
  ) %>%
  mutate(
    Frequency = Count / 250,
    Variable  = factor(Variable, levels = c("Z1","Z2","Z3")),
    Method    = factor(recode(Method,
                              vip = "VIP",
                              mi = "MI",
                              within = "VIP Type"),
                       levels = c("MI","VIP","VIP Type"))
  )

ggplot(df_long, aes(x = simulation, y = Frequency, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black") +
  facet_wrap(~ Variable, ncol = 3) +
  scale_fill_manual(values = c("MI" = "#1b9e77",      # Teal green
                               "VIP" = "#d95f02",     # Orange
                               "VIP Type" = "#7570b3")) + # Purple
  labs(x = "Simulation Scenario", y = "Selection Frequency", fill = "Method") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(size = 0.8, colour = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave("variable_selection_plot.eps", width = 10, height = 6, dpi = 400)

# Plot of the performance of each scenerio with respect to different alpha values
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)

data <- read_csv("~/Desktop/final1_with_alpha.csv")
data$alpha <- as.numeric(data$alpha)

aggregated_results_vip <- data %>%
  group_by(simulation, alpha) %>%
  summarise(
    mean_f1 = mean(f1_vip_mean, na.rm = TRUE),
    mean_recall = mean(recall_vip_mean, na.rm = TRUE),
    mean_precision = mean(precision_vip_mean, na.rm = TRUE),
    mean_error_rate = mean(error_rate_vip_mean, na.rm = TRUE),
    mean_type_II = mean(type_II_vip_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(Method = "VIP")

aggregated_results_mi <- data %>%
  group_by(simulation, alpha) %>%
  summarise(
    mean_f1 = mean(f1_mi_mean, na.rm = TRUE),
    mean_recall = mean(recall_mi_mean, na.rm = TRUE),
    mean_precision = mean(precision_mi_mean, na.rm = TRUE),
    mean_error_rate = mean(error_rate_mi_mean, na.rm = TRUE),
    mean_type_II = mean(type_II_mi_mean, na.rm = TRUE),
    mean_type_I = mean(type_I_mi_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(Method = "MI")

aggregated_results_within <- data %>%
  group_by(simulation, alpha) %>%
  summarise(
    mean_f1 = mean(f1_within_mean, na.rm = TRUE),
    mean_recall = mean(recall_within_mean, na.rm = TRUE),
    mean_precision = mean(precision_within_mean, na.rm = TRUE),
    mean_error_rate = mean(error_rate_within_mean, na.rm = TRUE),
    mean_type_II = mean(type_II_within_mean, na.rm = TRUE),
    mean_type_I = mean(type_I_within_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(Method = "Within VIP")

aggregated_results <- bind_rows(aggregated_results_mi)

plot_data <- aggregated_results %>%
  pivot_longer(cols = c(mean_f1, mean_recall, mean_precision, mean_error_rate, mean_type_II), 
               names_to = "Metric", values_to = "Value")

plot_data$alpha <- as.factor(plot_data$alpha)
plot_data$simulation <- factor(plot_data$simulation, 
                               levels = c("sim1", "sim2", "sim3", "sim4"), 
                               labels = c("3U,3N", "3U,4N", "3U,5N", "3U,6N"))


library(latex2exp)

custom_theme <- theme_minimal(base_size = 18) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = c(0.95,0.85),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length.y = unit(-0.25, "cm"),         # Make y ticks point inward
    axis.text.y = element_text(color = "black"),       # Optional: Red y-axis labels
    axis.ticks.y = element_line(color = "black"),      # Optional: Red y-axis ticks
    axis.line = element_line(size = 0.8, colour = "black")
  )
aggregated_results$simulation <- factor(aggregated_results$simulation, 
                               levels = c("sim1", "sim2", "sim3", "sim4"), 
                               labels = c("3U,3N", "3U,4N", "3U,5N", "3U,6N"))

ggplot(aggregated_results, aes(x = factor(alpha), y = mean_f1, group = simulation)) +
  geom_line(aes(linetype = simulation), color = "black", size = 1) +
  geom_point(aes(shape = simulation), color = "black", size = 3) +
  labs(
    y = TeX(""),
    x = TeX("$\\alpha$ Level")
  ) +
  scale_linetype_manual(values = c(
    "3U,3N" = "solid",
    "3U,4N" = "dashed",
    "3U,5N" = "dotted",
    "3U,6N" = "dotdash"
  )) +
  scale_shape_manual(values = c(
    "3U,3N" = 16,
    "3U,4N" = 17,
    "3U,5N" = 15,
    "3U,6N" = 18
  )) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    panel.grid = element_blank(),
    axis.ticks.length.y = unit(-0.25, "cm"),         # Make y ticks point inward
    axis.text.y = element_text(color = "black"),       # Optional: Red y-axis labels
    axis.ticks.y = element_line(color = "black"),
    axis.line = element_line(size = 0.8, colour = "black")
  )

 ggsave("f1.eps", width = 8, height = 6, dpi = 400)

 ggplot(aggregated_results, aes(x = factor(alpha), y = mean_precision, group = simulation)) +
   geom_line(aes(linetype = simulation), color = "black", size = 1) +
   geom_point(aes(shape = simulation), color = "black", size = 3) +
   labs(
     x = TeX("$\\alpha$ Level"),
     y = ""
   ) +
   scale_linetype_manual(values = c(
     "3U,3N" = "solid",
     "3U,4N" = "dashed",
     "3U,5N" = "dotted",
     "3U,6N" = "dotdash"
   )) +
   scale_shape_manual(values = c(
     "3U,3N" = 16,
     "3U,4N" = 17,
     "3U,5N" = 15,
     "3U,6N" = 18
   )) +
   theme_minimal(base_size = 16) +
   theme(
     plot.title = element_text(size = 18, face = "bold"),
     axis.title = element_text(size = 20),
     axis.text = element_text(size = 22),
     legend.title = element_blank(),
     legend.text = element_text(size = 18),
     legend.position = c(1, 0.03),
     legend.justification = c("right", "bottom"),
     panel.grid = element_blank(),
     axis.ticks.length.y = unit(-0.25, "cm"),         # Make y ticks point inward
     axis.text.y = element_text(color = "black"),       # Optional: Red y-axis labels
     axis.ticks.y = element_line(color = "black"),
     axis.line = element_line(size = 0.8, colour = "black")
   )
 ggsave("precision.eps", width = 8, height = 6, dpi = 400)
 
 ggplot(aggregated_results, aes(x = factor(alpha), y = mean_recall, group = simulation)) +
   geom_line(aes(linetype = simulation), color = "black", size = 1) +
   geom_point(aes(shape = simulation), color = "black", size = 3) +
   labs(
     x = TeX("$\\alpha$ Level"),
     y = ""
   ) +
   scale_linetype_manual(values = c(
     "3U,3N" = "solid",
     "3U,4N" = "dashed",
     "3U,5N" = "dotted",
     "3U,6N" = "dotdash"
   )) +
   scale_shape_manual(values = c(
     "3U,3N" = 16,
     "3U,4N" = 17,
     "3U,5N" = 15,
     "3U,6N" = 18
   )) +
   theme_minimal(base_size = 16) +
   theme(
     plot.title = element_text(size = 18, face = "bold"),
     axis.title = element_text(size = 20),
     axis.text = element_text(size = 22),
     legend.title = element_blank(),
     legend.text = element_text(size = 18),
     legend.position = c(0.95, 0.05),
     legend.justification = c("right", "bottom"),
     panel.grid = element_blank(),
     axis.ticks.length.y = unit(-0.25, "cm"),         # Make y ticks point inward
     axis.text.y = element_text(color = "black"),       # Optional: Red y-axis labels
     axis.ticks.y = element_line(color = "black"),
     axis.line = element_line(size = 0.8, colour = "black")
   )
 ggsave("recall.eps", width = 8, height = 6, dpi = 400)

 ggplot(aggregated_results, aes(x = factor(alpha), y = mean_type_I, group = simulation)) +
   geom_line(aes(linetype = simulation), color = "black", size = 1) +
   geom_point(aes(shape = simulation), color = "black", size = 3) +
   labs(
     x = TeX("$\\alpha$ Level"),
     y = ""
   ) +
   scale_linetype_manual(values = c(
     "3U,3N" = "solid",
     "3U,4N" = "dashed",
     "3U,5N" = "dotted",
     "3U,6N" = "dotdash"
   )) +
   scale_shape_manual(values = c(
     "3U,3N" = 16,
     "3U,4N" = 17,
     "3U,5N" = 15,
     "3U,6N" = 18
   )) +
   theme_minimal(base_size = 16) +
   theme(
     plot.title = element_text(size = 18, face = "bold"),
     axis.title = element_text(size = 20),
     axis.text = element_text(size = 22),
     legend.title = element_blank(),
     legend.text = element_text(size = 18),
     legend.position = c(0.95, 0.05),
     legend.justification = c("right", "bottom"),
     panel.grid = element_blank(),
     axis.ticks.length.y = unit(-0.25, "cm"),         # Make y ticks point inward
     axis.text.y = element_text(color = "black"),       # Optional: Red y-axis labels
     axis.ticks.y = element_line(color = "black"),
     axis.line = element_line(size = 0.8, colour = "black")
   )
ggsave("typeI.eps", width = 8, height = 6, dpi = 400)
