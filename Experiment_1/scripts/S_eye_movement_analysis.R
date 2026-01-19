library(eyesim)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)


setwd('/Users/bai/Documents/GitHub/MEOW/Experiment_1/')

M <- fread("data/eye_movement_data/group_fixations_behavior.csv")
sac <- fread("data/eye_movement_data/group_saccades.csv")
blk <- fread("data/eye_movement_data/group_blinks.csv")

coverage <- M[, .(
  n_fix = .N,
  total_fix_dur = sum(dur),
  mean_fix_dur = mean(dur)
), by=.(subj_id, task, trial_id)]

blk_summary <- blk[, .(n_blk = .N), by=.(subj_id, task, trial_id)]
coverage <- merge(coverage, blk_summary, all.x=TRUE, by=c("subj_id","task","trial_id"))
coverage[is.na(n_blk), n_blk := 0]

coverage[, quality := ifelse(n_fix < 2 | n_blk > 5, "bad", "good")]

cat(sprintf("Good trials: %d (%.1f%%)\n", 
            sum(coverage$quality=="good"), 
            mean(coverage$quality=="good")*100))

M_clean <- M[coverage[quality=="good"], on=.(subj_id, task, trial_id)]
fwrite(M_clean, "data/eye_movement_data/group_fixation_clean.csv")


M <- fread("data/eye_movement_data/group_fixations_clean.csv")

library(data.table)
library(spatstat)
library(ggplot2)

# Load and filter to ONE subject
M <- M[subj_id == 610]  # Test with subject 501

cat(sprintf("Subject 610: %d fixations\n", nrow(M)))
cat(sprintf("1-back trials: %d\n", M[task=="1_back", uniqueN(trial_id)]))
cat(sprintf("2-back trials: %d\n", M[task=="2_back", uniqueN(trial_id)]))

# Setup
screen_w <- 1920; screen_h <- 1080; sigma <- 80
M[, dur_ms := dur / 1000]

# Extract pair info
M[, stim_base := sub("_([AB])_.*", "", stim_id)]
M[, pair_member := sub(".*_([AB])_.*", "\\1", stim_id)]

# Validate extraction
cat("\nStim extraction check:\n")
print(M[1:5, .(stim_id, stim_base, pair_member)])

# Clean coordinates
M <- M[!is.na(x) & !is.na(y) & x >= 0 & x <= screen_w & y >= 0 & y <= screen_h]

# Heatmap function
create_heatmap <- function(x, y, dur, screen_w, screen_h, sigma) {
  if(length(x) < 2) return(matrix(0, 128, 128))
  weights <- dur / sum(dur)
  pp <- ppp(x, y, window=owin(c(0, screen_w), c(0, screen_h)))
  dens <- density(pp, sigma=sigma, weights=weights, dimyx=c(128, 128))
  return(as.matrix(dens))
}

# Compute heatmaps for 1-back (should be fast for 1 subject)
cat("\nComputing 1-back heatmaps...\n")
hm_1b <- M[task=="1_back", {
  cat(".")
  .(heatmap = list(create_heatmap(x, y, dur_ms, screen_w, screen_h, sigma)))
}, by=.(trial_id, stim_base, pair_member, condition, identity)]

cat(sprintf("\n1-back heatmaps: %d\n", nrow(hm_1b)))

# Compute for 2-back
cat("\nComputing 2-back heatmaps...\n")
hm_2b <- M[task=="2_back", {
  cat(".")
  .(heatmap = list(create_heatmap(x, y, dur_ms, screen_w, screen_h, sigma)),
    resp_key = first(resp_key),
    corr_resp = first(corr_resp))
}, by=.(trial_id, stim_base, pair_member, condition, identity, goal)]

hm_2b[, correct := (resp_key == corr_resp)]
cat(sprintf("\n2-back heatmaps: %d\n", nrow(hm_2b)))

# Test matching for ONE 2-back B trial
test_trial <- hm_2b[pair_member=="B"][1]
cat(sprintf("\nTest trial: %s, condition: %s\n", test_trial$stim_base, test_trial$condition))

match_1B <- hm_1b[stim_base==test_trial$stim_base & pair_member=="B"]
cat(sprintf("Found matching 1-back B: %s\n", ifelse(nrow(match_1B)>0, "YES", "NO")))

if(nrow(match_1B) > 0) {
  match_cor <- cor(as.vector(test_trial$heatmap[[1]]), as.vector(match_1B$heatmap[[1]]))
  cat(sprintf("Match correlation: %.3f\n", match_cor))
}
screen_w <- 1920
screen_h <- 1080
sigma <- 80  # as per paper

# Function: duration-weighted heatmap
create_heatmap <- function(x, y, dur, screen_w, screen_h, sigma) {
  if(length(x) < 2) return(matrix(0, 128, 128))
  
  # Normalize durations as weights
  weights <- dur / sum(dur)
  
  # Create point pattern with duration weights
  pp <- ppp(x, y, window=owin(c(0, screen_w), c(0, screen_h)))
  
  # Density with Gaussian smoothing
  dens <- density(pp, sigma=sigma, weights=weights, dimyx=c(128, 128))
  
  return(as.matrix(dens))
}

# Compute heatmaps for 1-back
hm_1b <- M[task=="1-back", .(
  heatmap = list(create_heatmap(x, y, dur, screen_w, screen_h, sigma)),
  stim_id = first(stim_id),
  condition = first(condition),
  identity = first(identity)
), by=.(subj_id, trial_id)]

# Compute heatmaps for 2-back  
hm_2b <- M[task=="2-back", .(
  heatmap = list(create_heatmap(x, y, dur, screen_w, screen_h, sigma)),
  stim_id = first(stim_id),
  condition = first(condition),
  identity = first(identity),
  correct = first(correct)
), by=.(subj_id, trial_id)]



library(data.table)
library(spatstat)
library(ggplot2)


# Extract pair info from stim_id: "mst_195_A_l2.png" → base="mst_195", member="A"
M[, stim_base := sub("_([AB])_.*", "", stim_id)]  # "mst_195"
M[, pair_member := sub(".*_([AB])_.*", "\\1", stim_id)]  # "A" or "B"

# Verify extraction
M[stim_id=="mst_195_A_l2.png", unique(.(stim_base, pair_member))]  # Should be: mst_195, A
M[stim_id=="mst_195_B_l2.png", unique(.(stim_base, pair_member))]  # Should be: mst_195, B

# Screen params
screen_w <- 1920; screen_h <- 1080; sigma <- 80

# Duration-weighted heatmap function
create_heatmap <- function(x, y, dur, screen_w, screen_h, sigma) {
  if(length(x) < 2) return(matrix(0, 128, 128))
  weights <- dur / sum(dur)
  pp <- ppp(x, y, window=owin(c(0, screen_w), c(0, screen_h)))
  dens <- density(pp, sigma=sigma, weights=weights, dimyx=c(128, 128))
  return(as.matrix(dens))
}

# Convert dur to ms
M[, dur_ms := dur / 1000]

# Compute heatmaps
hm_1b <- M[task=="1_back", .(
  heatmap = list(create_heatmap(x, y, dur_ms, screen_w, screen_h, sigma))
), by=.(subj_id, trial_id, stim_base, pair_member, condition, identity)]

hm_2b <- M[task=="2_back", .(
  heatmap = list(create_heatmap(x, y, dur_ms, screen_w, screen_h, sigma)),
  resp_key = first(resp_key),
  corr_resp = first(corr_resp)
), by=.(subj_id, trial_id, stim_base, pair_member, condition, identity, goal)]

# Add correct column
hm_2b[, correct := (resp_key == corr_resp)]

# PRIMARY: 1-back B → 2-back B (same lure)
reinst_primary <- hm_2b[pair_member=="B", {
  match_1B <- hm_1b[subj_id==.BY$subj_id & 
                      stim_base==.BY$stim_base & 
                      pair_member=="B"]
  
  if(nrow(match_1B) > 0) {
    match_cor <- cor(as.vector(heatmap[[1]]), as.vector(match_1B$heatmap[[1]]))
    
    # Mismatch baseline
    pool <- hm_1b[subj_id==.BY$subj_id & stim_base!=.BY$stim_base]
    if(nrow(pool) >= 50) {
      mismatch_cors <- sapply(sample(1:nrow(pool), 50), function(i)
        cor(as.vector(heatmap[[1]]), as.vector(pool$heatmap[[i]])))
      mismatch_cor <- mean(mismatch_cors)
    } else mismatch_cor <- NA
    
    .(type="1B→2B", 
      match_cor, mismatch_cor,
      match_z = atanh(match_cor), 
      mismatch_z = atanh(mismatch_cor),
      reinst = atanh(match_cor) - atanh(mismatch_cor))
  } else .(type="1B→2B", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
}, by=.(subj_id, trial_id, condition, correct, stim_base)]

# SECONDARY: 1-back A → 2-back B (original → lure)
reinst_secondary <- hm_2b[pair_member=="B", {
  match_1A <- hm_1b[subj_id==.BY$subj_id & 
                      stim_base==.BY$stim_base & 
                      pair_member=="A"]
  
  if(nrow(match_1A) > 0) {
    match_cor <- cor(as.vector(heatmap[[1]]), as.vector(match_1A$heatmap[[1]]))
    
    pool <- hm_1b[subj_id==.BY$subj_id & stim_base!=.BY$stim_base]
    if(nrow(pool) >= 50) {
      mismatch_cors <- sapply(sample(1:nrow(pool), 50), function(i)
        cor(as.vector(heatmap[[1]]), as.vector(pool$heatmap[[i]])))
      mismatch_cor <- mean(mismatch_cors)
    } else mismatch_cor <- NA
    
    .(type="1A→2B", 
      match_cor, mismatch_cor,
      match_z = atanh(match_cor), 
      mismatch_z = atanh(mismatch_cor),
      reinst = atanh(match_cor) - atanh(mismatch_cor))
  } else .(type="1A→2B", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
}, by=.(subj_id, trial_id, condition, correct, stim_base)]

# Combine
all_reinst <- rbind(reinst_primary, reinst_secondary)

# Summary
summary_stats <- all_reinst[!is.na(reinst), .(
  n = .N,
  mean_match = mean(match_cor),
  mean_mismatch = mean(mismatch_cor),
  mean_reinst = mean(reinst),
  se_reinst = sd(reinst)/sqrt(.N)
), by=.(type, condition)]

print(summary_stats)

# Visualize
ggplot(all_reinst[!is.na(reinst)], aes(x=condition, y=reinst, fill=condition)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.2) +
  geom_jitter(alpha=0.1, width=0.1) +
  facet_wrap(~type, scales="free_y") +
  labs(title="Gaze Reinstatement", y="Reinstatement Score (Fisher Z)", x="") +
  theme_minimal() +
  theme(legend.position="none")

ggsave("figures/reinstatement_conditions.pdf", width=10, height=5)







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





library(eyesim)
library(data.table)
library(ggplot2)
library(patchwork)

M <- M_ori[subj_id == 612]

# Extract pair info
M[, stim_base := sub("_([AB])_.*", "", stim_id)]
M[, pair_member := sub(".*_([AB])_.*", "\\1", stim_id)]
M[, dur_ms := dur / 1000]

# YOUR EXAMPLE: mst_195
example_base <- "mst_195"

# Create fixation_group objects
fg_1b_A <- with(M[task=="1_back" & stim_base==example_base & pair_member=="A"],
                fixation_group(x=x, y=y, duration=dur_ms, onset=onset))

fg_1b_B <- with(M[task=="1_back" & stim_base==example_base & pair_member=="B"],
                fixation_group(x=x, y=y, duration=dur_ms, onset=onset))

fg_2b_A <- with(M[task=="2_back" & stim_base==example_base & pair_member=="A"],
                fixation_group(x=x, y=y, duration=dur_ms, onset=onset))

fg_2b_B <- with(M[task=="2_back" & stim_base==example_base & pair_member=="B"],
                fixation_group(x=x, y=y, duration=dur_ms, onset=onset))

# Create eye_density objects
ed_1b_A <- eye_density(fg_1b_A, sigma=20, xbounds=c(0,1920), ybounds=c(0,1080), duration_weighted=TRUE)
ed_1b_B <- eye_density(fg_1b_B, sigma=20, xbounds=c(0,1920), ybounds=c(0,1080), duration_weighted=TRUE)
ed_2b_A <- eye_density(fg_2b_A, sigma=20, xbounds=c(0,1920), ybounds=c(0,1080), duration_weighted=TRUE)
ed_2b_B <- eye_density(fg_2b_B, sigma=20, xbounds=c(0,1920), ybounds=c(0,1080), duration_weighted=TRUE)

# Plot in RStudio viewer
plot(ed_1b_A, main="1-back A (trial 631)")
plot(ed_1b_B, main="1-back B (trial 632)")
plot(ed_2b_A, main="2-back A (trial 609)")
plot(ed_2b_B, main="2-back B (trial 611)")

# Compute similarities
cat("\n=== REINSTATEMENT SCORES (Fisher Z) ===\n")
cat(sprintf("1B → 2B (same lure): %.3f\n", similarity(ed_1b_B, ed_2b_B, method="fisherz")))
cat(sprintf("1A → 2B (original→lure): %.3f\n", similarity(ed_1b_A, ed_2b_B, method="fisherz")))
cat(sprintf("1A → 2A (same original): %.3f\n", similarity(ed_1b_A, ed_2b_A, method="fisherz")))
cat(sprintf("1B → 1A (within encoding): %.3f\n", similarity(ed_1b_B, ed_1b_A, method="fisherz")))
