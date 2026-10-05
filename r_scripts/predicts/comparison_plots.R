library(readr); library(ggplot2); library(dplyr); library(tidyr); library(forcats)

comparison <- read_csv("predicts/results_low_qc/comparison_embeddings_vs_qc.csv")

## --- Palette: 2 sole tinte, per il contrasto che conta --------------------
PAL <- c("QC metrics (4 variables)"   = "#2a78d6",
         "Foundation model embedding" = "#eb6834")

METRICS <- c(f1          = "F1",
             precision   = "Precision",
             sensitivity = "Sensitivity",
             specificity = "Specificity")

SRC_QC  <- "QC metrics (4 variables)"
SRC_EMB <- "Foundation model embedding"

PAL  <- c("#2a78d6", "#eb6834")
names(PAL) <- c(SRC_QC, SRC_EMB)

.src <- function(x) ifelse(x == "QC-only", SRC_QC, SRC_EMB)

## tema condiviso: un solo posto dove cambiare la scala tipografica
theme_fig <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.spacing     = unit(1.2, "lines"),
      strip.text        = element_text(hjust = 0, face = "bold", size = base_size),
      axis.text         = element_text(size = base_size - 1, colour = "#3d3c3a"),
      axis.title        = element_text(size = base_size),
      plot.title        = element_text(face = "bold", size = base_size + 4),
      plot.subtitle     = element_text(colour = "#52514e", size = base_size - 1),
      legend.text       = element_text(size = base_size),
      legend.position   = "top",
      legend.justification = "left",
      plot.background   = element_rect(fill = "white", colour = NA)
    )
}

## =============================== p1 ======================================
dl <- comparison |>
  mutate(label  = ifelse(embedding == "QC-only", model, paste(embedding, model)),
         source = .src(embedding)) |>
  select(label, source, all_of(names(METRICS))) |>
  pivot_longer(all_of(names(METRICS)), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(METRICS[metric], levels = METRICS))

ord <- comparison |>
  mutate(label = ifelse(embedding == "QC-only", model, paste(embedding, model))) |>
  arrange(balanced_accuracy) |> pull(label)
dl$label <- factor(dl$label, levels = ord)

p1 <- ggplot(dl, aes(value, label, colour = source)) +
  geom_segment(aes(x = 0, xend = value, yend = label),
               linewidth = 1, alpha = 0.30, lineend = "round") +
  geom_point(size = 3.6) +
  geom_text(aes(label = sprintf("%.3f", value)),
            hjust = -0.3, size = 4, colour = "#52514e") +
  facet_wrap(~ metric, ncol = 2) +
  scale_x_continuous(limits = c(0, 1.22), breaks = c(0, 0.5, 1),
                     expand = expansion(0)) +
  scale_colour_manual(values = PAL, name = NULL, limits = names(PAL))+
  labs(x = NULL, y = NULL,
       title = "Predicting low_qc cells: embeddings vs QC metrics",
       subtitle = "Stratified 5-fold CV, aggregated out-of-fold predictions (n = 22,767 cells)") +
  theme_fig(14) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#e8e8e4", linewidth = 0.3))

p1
ggsave("predicts/results_low_qc/fig_confronto_metriche.png", p1, width = 10, height = 5.5, dpi = 300)


# fig 2
## =============================== p2 (ROC space) ==========================

p2 <- comparison |>
  mutate(label = ifelse(embedding == "QC-only", model, paste(embedding, model)),
         source = .src(embedding),
         fpr    = 1 - specificity) |>
  ggplot(aes(fpr, sensitivity, colour = source)) +
  geom_abline(slope = 1, intercept = 0, colour = "#c9c9c4", linetype = "22",
              linewidth = 0.5) +
  annotate("text", x = 0.82, y = 0.78, label = "chance", angle = 45,
           size = 4, colour = "#9a9a94") +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = label), size = 4, colour = "#52514e",
                           box.padding = 0.6, seed = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(0.02)) +
  scale_colour_manual(values = PAL, name = NULL) +
  coord_equal() +
  labs(x = "1 − Specificity (false positive rate)",
       y = "Sensitivity (true positive rate)",
       title = "Where the signal is lost",
       subtitle = "Embeddings miss ~35–40% of low quality cells at comparable false positive rate") +
  theme_fig(14)

p2
ggsave("predicts/results_low_qc/fig_sens_spec.png", p2, width = 6, height = 6.2)
