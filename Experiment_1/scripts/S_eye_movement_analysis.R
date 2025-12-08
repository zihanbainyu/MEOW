library(eyesim)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("/Users/bai/Documents/GitHub/MEOW/Experiment_1/scripts")
subj_id <- 501
base_dir <- ".."
results_dir <- file.path(base_dir, "data")
fix_file <- file.path(results_dir, paste0("sub", subj_id), paste0("sub", subj_id, "_fixations.csv"))
data <- read.csv(fix_file)



fixation_groups <- list()
density_maps <- list()
for (tid in unique_trial_ids) {
  trial_data <- oneback_data %>% filter(trial_id == tid)
  if (nrow(trial_data) > 0) {
    fg <- fixation_group(x = trial_data$x, y = trial_data$y, 
                         onset = trial_data$onset, duration = trial_data$duration)
    fixation_groups[[as.character(tid)]] <- fg
    density_maps[[as.character(tid)]] <- eye_density(fg, sigma = 50, xbounds = c(0, 1920), ybounds = c(0, 1080))
  }
}

n_trials <- length(density_maps)
n_cols <- min(3, ceiling(sqrt(n_trials)))
n_rows <- ceiling(n_trials / n_cols)

par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 1))
for (tid in names(density_maps)) {
  plot(density_maps[[tid]], main = paste("Trial", tid))
}
par(mfrow = c(1, 1))

trial_ids <- names(density_maps)
n_trials <- length(trial_ids)
similarity_matrix <- matrix(NA, nrow = n_trials, ncol = n_trials,
                            dimnames = list(trial_ids, trial_ids))

for (i in 1:n_trials) {
  for (j in 1:n_trials) {
    similarity_matrix[i, j] <- similarity(density_maps[[trial_ids[i]]], 
                                          density_maps[[trial_ids[j]]], 
                                          method = 'pearson')
  }
}

sim_df <- melt(similarity_matrix, varnames = c("Trial1", "Trial2"), value.name = "Correlation")

p_heatmap <- ggplot(sim_df, aes(x = Trial1, y = Trial2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = paste("Gaze Pattern Similarity Matrix - Subject", subj_id),
       x = "Trial ID", y = "Trial ID") +
  coord_fixed()

print(p_heatmap)
ggsave(filename = paste0("gaze_similarity_heatmap_sub", subj_id, ".png"),
       plot = p_heatmap, width = 12, height = 10, dpi = 300)

write.csv(similarity_matrix, file = paste0("similarity_matrix_sub", subj_id, ".csv"), row.names = TRUE)


library(eyesim)
library(patchwork)
library(dplyr)
library(eyetrackingR)
library(eyelinker)
library(stringr)
# remotes::install_github("bbuchsbaum/eyesim")
# devtools::install_github("jwdink/eyetrackingR")

setwd("/Users/bai/Documents/GitHub/MEOW/Experiment_1/scripts")
subj_id <- 610
base_dir <- ".."
results_dir <- file.path(base_dir, "data")

fix_file <- file.path(results_dir, paste0("sub", subj_id), paste0("sub", subj_id, "_fixations.csv"))

data <- read.csv(fix_file)

# load 1b and 2b data
# visualize the gaze pattern of each trial
# find which trial belongs to which condition & identity & goal
# compute scanpath similarity between between each trial
# test if episodic memory representation of the original A item 1-back  is reinstated when seeing the lure B item during 2-back, indicated by a higher scanpath similarity/spatial & temporal correlation

```{r, echo=FALSE}
oneback_data <- data %>% filter(task=='1_back')
fg_1b <- fixation_group(x=oneback_data$x, y=oneback_data$y, onset = oneback_data$onset, duration = oneback_data$duration)
plot(fg_1b)

twoback_data <- data %>% filter(task=='2_back')
fg_2b <- fixation_group(x=twoback_data$x, y=twoback_data$y, onset = twoback_data$onset, duration = twoback_data$duration)
plot(fg_2b)

oneback_data %>% group_by(trial_id) %>% summarise(n=n()) %>% 
  {ggplot(., aes(trial_id, n)) + geom_col(fill="steelblue") + 
      geom_hline(yintercept=2, color="red", linetype=2) + theme_classic()}
twoback_data %>% group_by(trial_id) %>% summarise(n=n()) %>% 
  {ggplot(., aes(trial_id, n)) + geom_col(fill="steelblue") + 
      geom_hline(yintercept=2, color="red", linetype=2) + theme_classic()}

A_1_back <- data %>% filter(stim_id=='mst_021_A_l1.png' & task=='1_back')
B_1_back <- data %>% filter(stim_id=='mst_021_B_l1.png' & task=='1_back')
# A_2_back <- data %>% filter(stim_id=='mst_002_A_l1.png' & task=='2_back')
# B_2_back <- data %>% filter(stim_id=='mst_002_B_l1.png' & task=='2_back')

fg_1_back_A <- fixation_group(x=A_1_back$x, y=A_1_back$y, onset = A_1_back$onset, duration = A_1_back$duration)
fg_1_back_B <- fixation_group(x=B_1_back$x, y=B_1_back$y, onset = B_1_back$onset, duration = B_1_back$duration)
plot(fg_1_back_A)
plot(fg_1_back_B)
# fg_2_back_A <- fixation_group(
#   x=A_2_back$x,
#   y=A_2_back$y,
#   onset = A_2_back$onset,
#   duration = A_2_back$duration
# )
# fg_2_back_B <- fixation_group(
#   x=B_2_back$x,
#   y=B_2_back$y,
#   onset = B_2_back$onset,
#   duration = B_2_back$duration
# )
ed_fg_1_back_A <- eye_density(fg_1_back_A, sigma=50, xbounds=c(0,1920), ybounds=c(0,1080), normalize = 1)
ed_fg_1_back_B <- eye_density(fg_1_back_B, sigma=50, xbounds=c(0,1920), ybounds=c(0,1080), normalize = 1)
plot(ed_fg_1_back_A)
plot(ed_fg_1_back_B)
# plot(ed_fg_2_back_A)
# plot(ed_fg_2_back_B)



# ed_fg_2_back_A <- eye_density(fg_2_back_A, sigma=50, xbounds=c(0,1920), ybounds=c(0,1080)) 
# ed_fg_2_back_B <- eye_density(fg_2_back_B, sigma=50, xbounds=c(0,1920), ybounds=c(0,1080)) 
# s1 <- similarity(ed_fg_1_back_A,ed_fg_2_back_A, method = 'pearson')
# s2 <- similarity(ed_fg_1_back_B,ed_fg_2_back_B, method = 'pearson')
s3 <- similarity(ed_fg_1_back_A,ed_fg_1_back_B, method = 'fisherz')
# s4 <- similarity(ed_fg_2_back_A,ed_fg_2_back_B, method = 'pearson')
# s5 <- similarity(ed_fg_1_back_A,ed_fg_2_back_B, method = 'pearson')
# s6 <- similarity(ed_fg_1_back_B,ed_fg_2_back_A, method = 'pearson')

# similarity_results <- tibble::tibble(
#   comparison = c(
#     "A: 1-back vs 2-back",
#     "B: 1-back vs 2-back",
#     "1-back: A vs B",
#     "2-back: A vs B",
#     "Cross: A1 vs B2",
#     "Cross: B1 vs A2"
#   ),
#   pearson = c(s1, s2, s3, s4, s5, s6)
# )
# print(similarity_results)
print(s3)


sp1 <- scanpath(fg_1_back_A)
sp2 <- scanpath(fg_1_back_B)
eyesim:::multi_match(sp1, sp2, c(1920,1080))

example_data <- data %>%
  filter(
    stim_id %in% c("mst_002_A_l1.png", "mst_002_B_l1.png"),
    task %in% c("2_back")
  )
eyetab <- eye_table("x", "y", "duration", "onset", groupvar=c("stim_id"), data=example_data)
eyedens <- density_by(eyetab, groups=c("stim_id"), sigma=100, xbounds=c(0,1920), ybounds=c(0,1080))

p1 <- plot(eyedens$density[[1]])
p2 <- plot(eyedens$density[[2]])


print(p1)
print(p2)


set.seed(1234)
enc_dens <- eyedens %>% filter(stim_id == "mst_002_A_l1.png")
ret_dens <- eyedens %>% filter(stim_id == "mst_002_B_l1.png")

simres1 <- template_similarity(enc_dens, ret_dens, match_on="stim_id", method="fisherz", permutations=50)

t.test(simres1$eye_sim_diff)


print(plot(ed_fg_1_back_A))
print(plot(ed_fg_1_back_B))
print(plot(ed_fg_2_back_A))
print(plot(ed_fg_2_back_B))
