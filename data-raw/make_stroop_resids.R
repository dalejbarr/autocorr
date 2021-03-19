library(tidyverse)
library(lme4)

stroop_ML3 <- read_csv("StroopCleanSet.csv",
		       col_types = "iciicicicciiiic") %>%
  separate(block_name, c("ink_color", "j1", "j2", "j3"), "\\|") %>%
  separate(trial_name, c("word", "k1", "k2", "k3"), "\\|") %>%
  mutate(trial_latency = if_else(trial_error == 0L, NA_integer_, trial_latency),
	 trial_latency = if_else(trial_latency > 10000, NA_integer_, trial_latency),
	 word = tolower(word),
	 cong = (congruent == "Incongruent") - mean(congruent == "Incongruent"),
	 session_id = factor(session_id),
	 ink_color = factor(ink_color),
	 word = factor(word),
	 congruent = factor(congruent),
	 response = factor(trial_response)) %>%
  select(session_id, trial = trial_number,
	 ink_color, word, congruent, cong,
	 response,
	 latency = trial_latency)

mod <- lmer(latency ~ cong + (cong || session_id), stroop_ML3)

allres <- double(nrow(stroop_ML3))
allres[!is.na(stroop_ML3[["latency"]])] <- residuals(mod)

v1 <- round(sqrt(as.numeric(VarCorr(mod)[[1]][1, 1]))) # random intercept
v2 <- round(sqrt(as.numeric(VarCorr(mod)[[2]][1, 1]))) # random slope

stroop_mod <- list(
  fixed = round(fixef(mod)),
  covmx = matrix(c(v1^2, 0,
		   0, v2^2), nrow = 2),
  sigma = round(sigma(mod)),
  resid = split(allres, stroop_ML3[["session_id"]]))

usethis::use_data(stroop_ML3, overwrite = TRUE)
usethis::use_data(stroop_mod, overwrite = TRUE)
