library("RePsychLing")
library("mgcv")
library("tidyverse")

kkl <- as_tibble(KKL)

kkl$Trial <- scale(kkl$trial)
kkl$FirstTrial = kkl$first==1
kkl$Int = interaction(kkl$size, kkl$cardinal)
kkl$Soa = scale(kkl$SOA)

## see Baayen, Vashisth, Kliegl, Bates (2017),
## supplemental materials p. 12
## This is the same model but without by-subject factor smooths
## and without 'AR' correction (AR.start, rho)
## and a by-subject random intercept replacing the factor smooth
mod_nfs2 <- bam(lrt ~ sze * (spt + obj + grv) * orn +
		  s(Trial, by=Int) + 
		  s(subj, bs="re") +
		  s(subj, spt, bs="re") +
		  s(subj, grv, bs="re") +
		  s(subj, obj, bs="re") +
		  s(subj, orn, bs="re") +
		  s(subj, spt_orn, bs="re") +
		  s(Soa),
		data=kkl, method="fREML", discrete=TRUE)

kkl$resid_nfs2 <- residuals(mod_nfs2)

## restore missing trials with NAs
kkl_df <- kkl %>%
  complete(subj, trial) %>%
  select(subj, trial, resid_nfs2) %>%
  rename(resid = resid_nfs2) %>%
  arrange(subj, trial) %>%
  as.data.frame() ## so it's not a tibble

## put it into a matrix for easy resampling
## each row a subject, each column a trial
kkl_mx <- matrix(kkl_df$resid, byrow = TRUE, ncol = 800)

save(kkl_df, file = "data/kkl_df.rdata")
save(kkl_mx, file = "data/kkl_mx.rdata")
