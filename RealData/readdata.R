library(Hmisc)
library(haven)   # for read_xpt
library(dplyr)

# Acculturation: relevant in diverse cohorts
ACQ_H <- read_xpt(
  file="RealData/ACQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Alcohol Use: lifestyle confounder
ALQ_H <- read_xpt(
  file="RealData/ALQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Cardiovascular Health: Outcome
CDQ_H <- read_xpt(
  file="RealData/CDQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Diet Behavior & Nutrition: diet quality/energy intake as upstream confounders.
DBQ_H <- read_xpt(
  file="RealData/DBQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Demographics Data: essential baseline confounders
DEMO_H <- read_xpt(
  file="RealData/DEMO_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Diabetes: Mediator
# DIQ_H <- read_xpt(
#   file="RealData/DIQ_H.xpt",
#   col_select = NULL,
#   skip = 0,
#   n_max = Inf,
#   .name_repair = "unique"
# )

# Mental Health - Depression Screener: psychosocial confounder
DPQ_H <- read_xpt(
  file="RealData/DPQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Hospital Utilization & Access to Care: proxy for detection/management
GHB_H <- read_xpt(
  file="RealData/GHB_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Hospital Utilization & Access to Care: proxy for detection/management
HUQ_H <- read_xpt(
  file="RealData/HUQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Health Insurance: proxy for detection/management
HIQ_H <- read_xpt(
  file="RealData/HIQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Income: SES confounder
INQ_H <- read_xpt(
  file="RealData/INQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Medical Conditions: for reference, has questions regarding heart attack, obesity...
MCQ_H <- read_xpt(
  file="RealData/MCQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Physical Activity: lifestyle confounder upstream of obesity/diabetes/CVD.
PAQ_H <- read_xpt(
  file="RealData/PAQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Sleep Disorders: linked to obesity, diabetes, and CVD.
SLQ_H <- read_xpt(
  file="RealData/SLQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Smoking - Cigarette Use: major lifestyle confounder.
SMQ_H <- read_xpt(
  file="RealData/SMQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)

# Cholesterol - LDL & Triglycerides: Mediator
# TRIGLY_H <- read_xpt(
#   file="RealData/TRIGLY_H.xpt",
#   col_select = NULL,
#   skip = 0,
#   n_max = Inf,
#   .name_repair = "unique"
# )

# Weight History: BMI, exposure
WHQ_H <- read_xpt(
  file="RealData/WHQ_H.xpt",
  col_select = NULL,
  skip = 0,
  n_max = Inf,
  .name_repair = "unique"
)


# SEQN   Respondent sequence number
# LBDLDL: LDL-cholesterol (mg/dL)
# WHD010 Current self-reported height (inches)
# WHQ020 Current self-reported weight (pounds)
# MCQ160B Ever told had congestive heart failure
# MCQ160C Ever told you had coronary heart disease
# MCQ160D Ever told you had angina/angina pectoris
# MCQ160E Ever told you had heart attack
# MCQ160F Ever told you had a stroke

merged <- CDQ_H %>%                   # Cardiovascular Health: Outcome
  left_join(WHQ_H, by = "SEQN")  %>%  # Weight History: BMI, exposure
  # left_join(TRIGLY_H, by = "SEQN")  %>% # Cholesterol - LDL: Mediator
  left_join(GHB_H, by = "SEQN")  %>%  #Glycohemoglobin (%): Mediator
  # left_join(DIQ_H, by = "SEQN")  %>%  # Diabetes: Mediator
  left_join(ACQ_H, by = "SEQN")  %>%  # Acculturation: relevant in diverse cohorts
  left_join(ALQ_H, by = "SEQN")  %>%  # Alcohol Use: lifestyle confounder
  left_join(DBQ_H, by = "SEQN")  %>%  # Diet Behavior & Nutrition: diet quality/energy intake as upstream confounders.
  left_join(DEMO_H, by = "SEQN")  %>% # Demographics Data: essential baseline confounders
  left_join(DPQ_H, by = "SEQN")  %>%  # Mental Health - Depression Screener: psychosocial confounder
  left_join(HIQ_H, by = "SEQN")  %>%  # Health Insurance: proxy for detection/management
  left_join(HUQ_H, by = "SEQN")  %>%  # Hospital Utilization & Access to Care: proxy for detection/management
  left_join(INQ_H, by = "SEQN")  %>%  # Income: SES confounder
  left_join(MCQ_H, by = "SEQN")  %>%  # Medical Conditions: for reference, has questions regarding heart attack, obesity...
  left_join(PAQ_H, by = "SEQN")  %>%  # Physical Activity: lifestyle confounder upstream of obesity/diabetes/CVD.
  left_join(SLQ_H, by = "SEQN")  %>%  # Sleep Disorders: linked to obesity, diabetes, and CVD.
  left_join(SMQ_H, by = "SEQN")       # Smoking - Cigarette Use: major lifestyle confounder.

key_vars <- c(
  
  # ID
  "SEQN",
  
  # Demographics
  "RIDAGEYR","RIAGENDR","RIDRETH1","DMDEDUC2",
  
  # Income
  "INDFMIN2",
  #"INDHHIN2",
  #"INDFMPIR",
  
  # Alcohol
  "ALQ101",
  #"ALQ110","ALQ120Q","ALQ120U","ALQ130","ALQ141Q","ALQ141U","ALQ151","ALQ160",
  
  # Smoking
  "SMQ020",
  #"SMQ040","SMD641","SMD650",
  
  # Physical activity
  "PAQ605",
  #"PAQ620","PAQ635","PAQ650","PAQ665",
  #"PAQ610","PAD615",#"PAQ625","PAD630","PAQ640","PAD645","PAQ655","PAD660","PAQ670","PAD675","PAQ680","PAQ706",
  
  # Diet behavior
  "DBQ700",
  #"DBD895","DBD905",
  #"DBD900",
  
  # Sleep
  "SLD010H",
  #"SLQ050","SLQ060",
  
  # Health insurance
  "HIQ011",
  #"HIQ031A","HIQ031B","HIQ031C","HIQ031D","HIQ031E","HIQ031F","HIQ031G","HIQ031H","HIQ031I","HIQ031J",
  
  # Obesity / BMI
  "WHD010","WHD020",
  # "BMI",
  
  # Diabetes
  "LBXGH",
  # "DIQ010","DIQ050","DID040",
  #"diab_prob_type2"
  
  # Cholesterol - LDL: Mediator
  # "LBDLDL", "LBXTR",
  
  # Cardiovascular conditions
  "CDQ001","CDQ010","MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"
  #"CDQ002","CDQ003","CDQ004","CDQ005","CDQ006","CDQ008","CDQ009A","CDQ009B","CDQ009C","CDQ009D","CDQ009E","CDQ009F","CDQ009G","CDQ009H",
  
  # Family history
  #"MCQ300A"
)

analysis_variables <- merged %>%
  dplyr::select(any_of(key_vars))
# colSums(is.na(analysis_variables))

analysis_variables$RIAGENDR <- 
  ifelse(analysis_variables$RIAGENDR %in% c(7, 9),
         NA, 
         analysis_variables$RIAGENDR)

analysis_variables$DMDEDUC2 <- 
  ifelse(analysis_variables$DMDEDUC2 %in% c(7, 9),
         NA, 
         analysis_variables$DMDEDUC2)

analysis_variables$INDFMIN2 <- 
  ifelse(analysis_variables$INDFMIN2 %in% c(77, 99),
         NA, 
         analysis_variables$INDFMIN2)

analysis_variables$ALQ101 <- 
  ifelse(analysis_variables$ALQ101 %in% c(7, 9),
         NA, 
         analysis_variables$ALQ101)

analysis_variables$SMQ020 <- 
  ifelse(analysis_variables$SMQ020 %in% c(7, 9),
         NA, 
         analysis_variables$SMQ020)

analysis_variables$PAQ605 <- 
  ifelse(analysis_variables$PAQ605 %in% c(7, 9),
         NA, 
         analysis_variables$PAQ605)

analysis_variables$DBQ700 <- 
  ifelse(analysis_variables$DBQ700 %in% c(7, 9),
         NA, 
         analysis_variables$DBQ700)

analysis_variables$SLD010H <- 
  ifelse(analysis_variables$SLD010H %in% c(77, 99),
         NA, 
         analysis_variables$SLD010H)

analysis_variables$HIQ011 <- 
  ifelse(analysis_variables$HIQ011 %in% c(7, 9),
         NA, 
         analysis_variables$HIQ011)

analysis_variables$WHD010 <- 
  ifelse(analysis_variables$WHD010 %in% c(7777, 9999),
         NA, 
         analysis_variables$WHD010)

analysis_variables$WHD020 <- 
  ifelse(analysis_variables$WHD020 %in% c(7777, 9999),
         NA, 
         analysis_variables$WHD020)

# analysis_variables$DIQ010 <- 
#   ifelse(analysis_variables$DIQ010 %in% c(7, 9),
#          NA, 
#          analysis_variables$DIQ010)

# analysis_variables$DID040 <- 
#   ifelse(analysis_variables$DID040 %in% c(666,777,999),
#          NA, 
#          analysis_variables$DID040)

# analysis_variables$DIQ050 <- 
#   ifelse(analysis_variables$DIQ050 %in% c(7, 9),
#          NA, 
#          analysis_variables$DIQ050)

analysis_variables$CDQ001 <- 
  ifelse(analysis_variables$CDQ001 %in% c(7, 9),
         NA, 
         analysis_variables$CDQ001)

analysis_variables$CDQ010 <- 
  ifelse(analysis_variables$CDQ010 %in% c(7, 9),
         NA, 
         analysis_variables$CDQ010)

analysis_variables$MCQ160B <- 
  ifelse(analysis_variables$MCQ160B %in% c(7, 9),
         NA, 
         analysis_variables$MCQ160B)

analysis_variables$MCQ160C <- 
  ifelse(analysis_variables$MCQ160C %in% c(7, 9),
         NA, 
         analysis_variables$MCQ160C)

analysis_variables$MCQ160D <- 
  ifelse(analysis_variables$MCQ160D %in% c(7, 9),
         NA, 
         analysis_variables$MCQ160D)

analysis_variables$MCQ160E <- 
  ifelse(analysis_variables$MCQ160E %in% c(7, 9),
         NA, 
         analysis_variables$MCQ160E)

analysis_variables$MCQ160F <- 
  ifelse(analysis_variables$MCQ160F %in% c(7, 9),
         NA, 
         analysis_variables$MCQ160F)


# analysis_variables <- analysis_variables %>%
#   mutate(
#     diab_prob_type2 = case_when(
#       DIQ010 == 1 & (is.na(DID040) | DID040 > 30 | is.na(DIQ050) | DIQ050 != 1) ~ 1,  # Type II
#       DIQ010 == 1 & DID040 <= 30 & DIQ050 == 1 ~ 0,                                   # Type I → not mediator
#       TRUE ~ 0                                                                        # No diabetes
#     )
#   )

analysis_variables <- analysis_variables %>%
  mutate(
    height_m = WHD010 * 0.0254,    # 1 inch = 0.0254 m
    weight_kg = WHD020 * 0.453592, # 1 pound = 0.453592 kg
    BMI = weight_kg / (height_m^2),
    obesity = ifelse(BMI >= 31, 1, ifelse(BMI <= 29, 0, NA))
  )



dataset_vars <- c(
  
  # ID
  #"SEQN",
  
  # Demographics
  "RIDAGEYR","RIAGENDR","RIDRETH1","DMDEDUC2",
  
  # Income
  "INDFMIN2",
  #"INDHHIN2",
  #"INDFMPIR",
  
  # Alcohol
  "ALQ101",
  #"ALQ110","ALQ120Q","ALQ120U","ALQ130","ALQ141Q","ALQ141U","ALQ151","ALQ160",
  
  # Smoking
  "SMQ020",
  #"SMQ040","SMD641","SMD650",
  
  # Physical activity
  "PAQ605",
  #"PAQ620","PAQ635","PAQ650","PAQ665",
  #"PAQ610","PAD615",#"PAQ625","PAD630","PAQ640","PAD645","PAQ655","PAD660","PAQ670","PAD675","PAQ680","PAQ706",
  
  # Diet behavior
  "DBQ700",
  #"DBD895","DBD905",
  #"DBD900",
  
  # Sleep
  "SLD010H",
  #"SLQ050","SLQ060",
  
  # Health insurance
  "HIQ011",
  #"HIQ031A","HIQ031B","HIQ031C","HIQ031D","HIQ031E","HIQ031F","HIQ031G","HIQ031H","HIQ031I","HIQ031J",
  
  # Obesity / BMI
  "obesity",
  #"WHD010","WHD020",
  "BMI",
  
  # Diabetes
  "LBXGH",
  # "diab_prob_type2",
  #"DIQ010","DIQ050","DID040",
  
  # Cholesterol - LDL: Mediator
  # "LBXTR",
  # "LBDLDL",
  
  # Cardiovascular conditions
  "CDQ001","CDQ010","MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"
  #"CDQ002","CDQ003","CDQ004","CDQ005","CDQ006","CDQ008","CDQ009A","CDQ009B","CDQ009C","CDQ009D","CDQ009E","CDQ009F","CDQ009G","CDQ009H",
  
  # Family history
  #"MCQ300A"
)

analysis_dataset <- analysis_variables %>%
  dplyr::select(any_of(dataset_vars))

analysis_dataset_complete <- na.omit(analysis_dataset)
analysis_dataset_complete <- analysis_dataset_complete %>%
  mutate(
    RIAGENDR = ifelse(RIAGENDR==1,1,0),
    ALQ101 = ifelse(ALQ101==1,1,0),
    SMQ020 = ifelse(SMQ020==1,1,0),
    PAQ605 = ifelse(PAQ605==1,1,0),
    HIQ011 = ifelse(HIQ011==1,1,0),
    CDQ001 = ifelse(CDQ001==1,1,0),
    CDQ010 = ifelse(CDQ010==1,1,0),
    MCQ160B = ifelse(MCQ160B==1,1,0),
    MCQ160C = ifelse(MCQ160C==1,1,0),
    MCQ160D = ifelse(MCQ160D==1,1,0),
    MCQ160E = ifelse(MCQ160E==1,1,0),
    MCQ160F = ifelse(MCQ160F==1,1,0)
  )
saveRDS(analysis_dataset_complete, file = "RealData/analysis_dataset_complete.rds")


d_sd <- (mean(analysis_dataset_complete$BMI[analysis_dataset_complete$obesity==1])
         -mean(analysis_dataset_complete$BMI[analysis_dataset_complete$obesity==0]))/sd(analysis_dataset_complete$BMI)
hist(analysis_dataset_complete$BMI,100)

length(analysis_dataset_complete$BMI[analysis_dataset_complete$BMI <= 37 & analysis_dataset_complete$BMI >= 23])
