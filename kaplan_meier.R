# Packages
library(xlsx)
library(tidyverse)
library(survival)
library(survminer)

# Import data
base <- read.xlsx("data/big_table_base.xlsx",1)
flup <- read.xlsx("data/big_table_flup.xlsx",1)

# Create dictionary to group diagnosis into ALL, CLL & NHL categories
diag_dict <-
  list("Acute B lymphoblastic leukemia [B-ALL]" = "ALL",
       "Richter's syndrome" = "NHL",
       "Diffuse large B cell lymphoma [DLBCL]" = "NHL",
       "Follicular lymphoma [FL]" = "NHL",
       "Other lymphoma type" = "NHL",
       "Chronic lymphocytic leukemia [CLL]" = "CLL",
       "Other disease" = "Other",
       "FL transformed to DLBCL" = "NHL",
       "Chronic myeloid leukemia [CML]" = "NHL")

# Get progression dates (refractory, relapse or progression) from the flup df
progression_df <-
  flup %>% 
  filter(calculate_response == "Refractory" | calculate_response == "Relapse/Progression") %>% 
  dplyr::mutate(
    progression_date = flup_visit_date
  ) %>% 
  select(reg, progression_date)

# Now we build the OS/PFS df from the previously obtained dfs and dictionary
survival_df <- base %>% 
  # Obtain progression dates from the progression df
  left_join(progression_df,
            by = "reg") %>% 
  # Select baseline and time variables needed
  select(reg, background_legal_frame, background_sex,
         background_diagnosis, infusion_date, progression_date, fu_exitus_date, last_visit) %>% 
  # Create new survival columns
  mutate(
    # create diagnostic groups with previously created dictionary
    diag_group = recode(background_diagnosis, !!!diag_dict),
    # obtain last follow-up date (last visit or death)
    last_date = if_else(!is.na(fu_exitus_date), fu_exitus_date, last_visit),
    # calculate the follow-up time to death or last visit
    surv_time = last_visit - infusion_date,
    # death status at last follow-up day
    died = if_else(!is.na(fu_exitus_date), TRUE, FALSE),
    # PFS
    pfs_last_date = if_else(!is.na(progression_date), progression_date, last_visit),
    pfs_time = pfs_last_date - infusion_date,
    progressed_died = if_else(!is.na(progression_date), TRUE, FALSE)
  ) %>% 
  relocate(diag_group, .before = infusion_date) %>% 
  # Remove patients missing infusion date
  filter(
    # remove patients missing the infusion date
    !is.na(infusion_date),
    # remove diags other than ALL, CLL & NHL
    diag_group != "Other"
    )

# Create Kaplan-Meier curve objects
os <- survfit(Surv(surv_time, died)~1, data=survival_df)
pfs <- survfit(Surv(pfs_time, progressed_died)~1, data=survival_df)
os_diag <- survfit(Surv(surv_time, died) ~ diag_group,
                   data=survival_df)
pfs_diag <- survfit(Surv(surv_time, progressed_died) ~ diag_group,
                    data=survival_df)

# Create Kaplan-Meier plots
plot_os <-
  ggsurvplot(os, data = survival_df,
             risk.table = TRUE,
             censor = TRUE,
             conf.int = F,
             
             break.time.by = 150,  # axis ticks every 150 days
             # xlim = c(0, 100),
             xlab = "Time (days)", ylab = "OS probability",
             
             palette = "#e41a1c", # color for the line
             
             legend = "none"
  )

plot_pfs <-
  ggsurvplot(pfs, data = survival_df,
             risk.table = TRUE,
             censor = TRUE,
             conf.int = F,
             
             break.time.by = 150,  # axis ticks every 150 days
             # xlim = c(0, 100),
             xlab = "Time (days)", ylab = "PFS probability",
             
             palette = "#e41a1c", # color for the line
             
             legend = "none"
  )

plot_os_diag <-
  ggsurvplot(os_diag, data = survival_df,
             risk.table = TRUE,
             censor = TRUE,
             conf.int = F,
             
             break.time.by = 150,  # axis ticks every 150 days
             # xlim = c(0, 100),
             xlab = "Time (days)", ylab = "OS probability",
             
             palette = c("#d8b365", "#5ab4ac","#5e3c99") # colors for the lines
  )

plot_pfs_diag <-
  ggsurvplot(pfs_diag, data = survival_df,
             risk.table = TRUE,
             censor = TRUE,
             conf.int = F,
             
             break.time.by = 150,  # axis ticks every 150 days
             # xlim = c(0, 100),
             xlab = "Time (days)", ylab = "PFS probability",
             
             palette = c("#d8b365", "#5ab4ac","#5e3c99") # colors for the lines
  )

# Save plots

# workaround method to ggsave the list of objects created by ggsurvplot
# (https://github.com/kassambara/survminer/issues/152#issuecomment-938941051)
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave("figures/plot_os.png", plot_os)
ggsave("figures/plot_pfs.png", plot_pfs)
ggsave("figures/plot_os_diag.png", plot_os_diag)
ggsave("figures/plot_pfs_diag.png", plot_pfs_diag)
