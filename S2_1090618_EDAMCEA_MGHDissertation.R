################################################################################
##TITLE: CEA for the Electronic Clinical Decision Support for Acute Fever Management
##COURSE: MSc in Modelling for Global Health
##CANDIDATE NUMBER: 1090618
################################################################################
rm(list = ls())
set.seed(100)

##SET WORKING DIRECTORY 
#uncomment the below code:
#setwd(#User modified by setting own working directory, cloned Github repository)
getwd()

##LOAD PACKAGES
install.packages('haven')
library(haven)
install.packages('dplyr')
library(dplyr)
library(ggplot2)
install.packages('data.tree')
library(data.tree)
install.packages('DiagrammeR')
library(DiagrammeR)
install.packages(heemod)
library(heemod)
library(Matrix)
library(lme4)
library(stringr)
library(tidyr)
install.packages("fmsb")
library(fmsb)
install.packages('flextable')
library(flextable)
library(scales)


##UPLOAD NEWEST CLINICAL TRIAL DATASET
#uncomment the below code:
#clinicaldata <- read_dta(#add path to 'S5_EDAM_all_Datasets_21May2025.dta' in Supplementary Information S5)

##CREATE VARIABLE FOR A PATIENT'S DANGER SIGNS
clinicaldata <- clinicaldata %>%
  mutate(
    observed_dangersigns = as.integer(
      (!is.na(ask_unable_drink) & ask_unable_drink == 1) |
        (!is.na(ask_vomit) & ask_vomit == 1) |
        (!is.na(ask_convulsion) & ask_convulsion == 1) |
        (!is.na(look_unable_sit) & look_unable_sit == 1) |
        (!is.na(look_convulsing) & look_convulsing == 1) |
        (!is.na(ask_confused) & ask_confused == 1) |
        (!is.na(mental_state) & mental_state %in% c(2, 3, 4))
    )
  )

clinicaldata <- clinicaldata %>%
  mutate(
    measured_dangersigns = as.integer(
      (!is.na(respiratory_rate) & respiratory_rate >= 25) |
        (!is.na(oxygen_saturation) & oxygen_saturation < 90) |
        (!is.na(heart_rate) & (heart_rate < 40 | heart_rate > 140)) |
        (!is.na(systolic_bp) & systolic_bp < 90) |
        (!is.na(diastolic_bp) & diastolic_bp < 30)
    )
  )

clinicaldata <- clinicaldata %>%
  mutate(
    dangersigns = if_else(
      observed_dangersigns == 1 | measured_dangersigns == 1,
      1L, 0L
    )
  )


################################################################################
##DETERMINE TIME SPENT IN TRIAL (NECESSARY FOR APPT TIME COST INPUTS)
clinicaldata$started_time <- as.POSIXct(clinicaldata$started_time, format = "%Y-%m-%d %H:%M:%S")
clinicaldata$completed_time <- as.POSIXct(clinicaldata$completed_time, format = "%Y-%m-%d %H:%M:%S")

clinicaldata$duration <- clinicaldata$completed_time - clinicaldata$started_time

duration_arm0 <- clinicaldata %>%
  filter(arm == 0, !is.na(duration)) %>%
  select(duration)

duration_arm1 <- clinicaldata %>%
  filter(arm == 1, !is.na(duration)) %>%
  select(duration)

#Convert to minutes
duration_arm0$duration <- as.numeric(duration_arm0$duration, units = "mins")
duration_arm1$duration <- as.numeric(duration_arm1$duration, units = "mins")

duration_arm0 <- duration_arm0 %>%
  mutate(arm = "Arm 0")

duration_arm1 <- duration_arm1 %>%
  mutate(arm = "Arm 1")

duration_combined <- bind_rows(duration_arm0, duration_arm1)

#Visualise
ggplot(duration_combined, aes(x = arm, y = duration, fill = arm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Duration by Arm", y = "Duration (minutes)", x = "Study Arm") +
  scale_fill_manual(values = c("Arm 0" = "#66c2a5", "Arm 1" = "#fc8d62")) +
  coord_cartesian(ylim = c(0, 100)) +  # adjust limits as needed
  theme_minimal()

#Looking at outliers
summary(duration_combined$duration)

t.test(duration_arm0$duration, duration_arm1$duration)

#Investigating time difference between arms
wilcox.test(duration_arm0$duration, duration_arm1$duration)
#Trust Wilcox test due to outliers

#Calculate median percent increase from control to intervention arm
medianduration_arm0 <- median(duration_arm0$duration, na.rm = TRUE)
medianduration_arm1 <- median(duration_arm1$duration, na.rm = TRUE)
durationpercent_increase <- ((medianduration_arm1 - medianduration_arm0) / medianduration_arm0) * 100

#For sensitivity and scenario analyses: Get quartiles for each arm
quartiles_arm0 <- quantile(duration_arm0$duration, probs = c(0.25, 0.75), na.rm = TRUE)
quartiles_arm1 <- quantile(duration_arm1$duration, probs = c(0.25, 0.75), na.rm = TRUE)

#Extract values
q25_arm0 <- quartiles_arm0[1]
q75_arm0 <- quartiles_arm0[2]
q25_arm1 <- quartiles_arm1[1]
q75_arm1 <- quartiles_arm1[2]

#Calculate percent increases
q25_pct_increase <- 100 * (q25_arm1 - q25_arm0) / q25_arm0
q75_pct_increase <- 100 * (q75_arm1 - q75_arm0) / q75_arm0

quantile(duration_arm1$duration, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
################################################################################
##DATA CLEANING & COLLECTING PROBABILITIES FOR DECISION TREE
################################################################################
##CONTROL ARM PROBABILITIES
pControl <- 0.5

##MALARIA RDT RESULT PROBABILITIES: CONTROL
pMalaria_Control <- clinicaldata %>%
  filter(arm == 0) %>%
  summarise(percentage = mean(malaria_rdt == 1, na.rm = TRUE))

pNoMalaria_Control <- clinicaldata %>%
  filter(arm == 0) %>%
  summarise(percentage = mean(malaria_rdt == 2, na.rm = TRUE))

pNoMalariaRDT_Control <- clinicaldata %>% #Those who did not take the RDT
  filter(arm == 0) %>%
  summarise(percentage = mean(malaria_rdt == 3, na.rm = TRUE))

#PRIMARY OUTCOME (ANTIBIOTIC PRESCRIBED) PROBABILITIES: CONTROL
pMalariaNoAB_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 1) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))

pMalariaAB_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 1) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))

pNoMalariaNoAB_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 2) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))

pNoMalariaAB_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 2) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))

pNoMalariaRDTNoAB_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 3) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))

pNoMalariaRDTAB_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 3) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))

##DIAGNOSIS PROBABILITIES: CONTROL (EXPLORATORY PURPOSES ONLY)
Malaria_Count_Control <- nrow(clinicaldata %>%
                            filter(arm == 0, malaria_rdt == 1))

Malaria_Diagnoses_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 1) %>%   # Filter for arm == 0 and malaria_rdt == 1
  group_by(diagnosis_combined) %>%                        # Group by diagnosis
  summarise(count = n(), .groups = 'drop') %>%   # Count occurrences of each diagnosis
  mutate(percentage = (count / Malaria_Count_Control) * 100) %>%   # Calculate percentage relative to the entire group
  arrange(desc(percentage))                      # Sort by percentage in descending order
as.data.frame(Malaria_Diagnoses_Control)

NoMalaria_Count_Control <- nrow(clinicaldata %>%
                                filter(arm == 0, malaria_rdt == 2))

NoMalaria_Diagnoses_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 2) %>%   # Filter for arm == 0 and malaria_rdt == 1
  group_by(diagnosis_combined) %>%                        # Group by diagnosis
  summarise(count = n(), .groups = 'drop') %>%   # Count occurrences of each diagnosis
  mutate(percentage = (count / NoMalaria_Count_Control) * 100) %>%   # Calculate percentage relative to the entire group
  arrange(desc(percentage))                      # Sort by percentage in descending order
as.data.frame(NoMalaria_Diagnoses_Control)

NoMalariaRDT_Count_Control <- nrow(clinicaldata %>%
                                  filter(arm == 0, malaria_rdt == 3))

NoMalariaRDT_Diagnoses_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 3) %>%   # Filter for arm == 0 and malaria_rdt == 1
  group_by(diagnosis_combined) %>%                        # Group by diagnosis
  summarise(count = n(), .groups = 'drop') %>%   # Count occurrences of each diagnosis
  mutate(percentage = (count / NoMalariaRDT_Count_Control) * 100) %>%   # Calculate percentage relative to the entire group
  arrange(desc(percentage))                      # Sort by percentage in descending order
as.data.frame(NoMalariaRDT_Diagnoses_Control)


##ANTIBIOTICS BY TYPE PROBABILITIES: CONTROL
#This is asking: 'of the antibiotics prescribed in each control branch malaria RDT result node, 
#what types of antibiotics were they?' This is needed for weighted cost inputs
AB_WithMalaria_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 1) %>%   
  mutate(which_antibiotics = ifelse(which_antibiotics == "", "None", which_antibiotics)) %>%
  group_by(which_antibiotics) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

as.data.frame(AB_WithMalaria_Control)  

AB_WithoutMalaria_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 2) %>%   
  mutate(which_antibiotics = ifelse(which_antibiotics == "", "None", which_antibiotics)) %>%
  group_by(which_antibiotics) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

as.data.frame(AB_WithoutMalaria_Control)  

AB_WithoutMalariaRDT_Control <- clinicaldata %>%
  filter(arm == 0, malaria_rdt == 3) %>%   
  mutate(which_antibiotics = ifelse(which_antibiotics == "", "None", which_antibiotics)) %>%
  group_by(which_antibiotics) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

as.data.frame(AB_WithoutMalariaRDT_Control)  

##REFERRAL PROBABILITIES: CONTROL
library(haven)
library(dplyr)
clinicaldata$d7_patient_recovered <- as.numeric(clinicaldata$d7_patient_recovered)

pReferredH_Control <- clinicaldata %>%
  filter(arm == 0) %>%
  summarise(percentage = mean(hosp_referral == 1, na.rm = TRUE))


pHospitalRef_MalariaNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 1,
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_MalariaAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 1,
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pHospitalRef_NoMalariaNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 2,
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pHospitalRef_NoMalariaAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 2,
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pHospitalRef_NoMalariaRDTNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 3,
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_NoMalariaRDTAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 3,
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

##RECOVERY AT FOLLOW UP PROBABILITIES: CONTROL

pRecovered_MalariaNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 1,
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_MalariaAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 1,
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoMalariaNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 2,
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoMalariaAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 2,
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)


pRecovered_NoMalariaRDTNoAB_Control <-  clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 3,
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoMalariaRDTAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 3,
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

##FOLLOW UP -- UNSCHEDULED REPRESENTING AT A PHC PROBABILITIES: CONTROL
pRepresent_MalariaNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 1,
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pRepresent_MalariaAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 1,
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_NoMalariaNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 2,
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_NoMalariaAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 2,
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)



pRepresent_NoMalariaRDTNoAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 3,
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_NoMalariaRDTAB_Control <- clinicaldata %>%
  filter(
    arm == 0,
    malaria_rdt == 3,
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


################################################################################
##INTERVENTION ARM PROBABILITIES 
pIntervention <- 0.5

##DANGER SIGN(S) PRESENT PROBABILITY: INTERVENTION 
pDangerSigns_Intervention <- clinicaldata %>%
  filter(arm == 1) %>%
  summarise(percentage = mean(dangersigns == 1, na.rm = TRUE))
pNoDangerSigns_Intervention <- clinicaldata %>%
  filter(arm == 1) %>%
  summarise(percentage = mean(dangersigns == 0, na.rm = TRUE))

##DANGER SIGN(S) PRESENT: (SEVERE) MALARIA RDT RESULT PROBABILITIES: INTERVENTION
pSevereMalaria_Intervention <- clinicaldata %>% #nobody had severe malaria
  filter(arm == 1, dangersigns == 1) %>%
  summarise(percentage = mean(malaria_rdt == 1, na.rm = TRUE))

pNoSevereMalaria_Intervention <- clinicaldata %>% 
  filter(arm == 1, dangersigns == 1) %>%
  summarise(percentage = mean(malaria_rdt == 2, na.rm = TRUE))

pNoSevereMalariaRDT_Intervention <- clinicaldata %>% ##Those who did not take the RDT
  filter(arm == 1, dangersigns == 1) %>%
  summarise(percentage = mean(malaria_rdt == 3, na.rm = TRUE))


##NO DANGER SIGN(S) PRESENT: SYMPTOMS ASSESSMENT & MALARIA RDT RESULT PROBABILITIES: INTERVENTION
#Everyone in the no danger signs branches receives malaria RDT. No clinical trial probabilities exist to account for otherwise
pMalariaRDT_AfterSymptoms_Intervention <- 1 

#Symptoms assessment in intervention arm establishes whether the patient has upper respiratory (UR) symptoms
pURSymptoms_Intervention <- clinicaldata %>% 
  filter(arm == 1, dangersigns == 0) %>%
  summarise(percentage = mean(resp_symptoms != ""))

pNoURSymptoms_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 0) %>%
  summarise(percentage = mean(resp_symptoms == ""))

#PRIMARY OUTCOME (ANTIBIOTIC PRESCRIBED) PROBABILITIES: INTERVENTION
pNoSevereMalariaNoAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 1, malaria_rdt == 2) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))
  
pNoSevereMalariaAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 1, malaria_rdt == 2) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))

pSevereMalariaNoAB_Intervention <- clinicaldata %>% #nobody had severe malaria
  filter(arm == 1, dangersigns == 1, malaria_rdt == 1) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))

pSevereMalariaAB_Intervention <- clinicaldata %>%   #nobody had severe malaria
  filter(arm == 1, dangersigns == 1, malaria_rdt == 1) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))


pNoSevereMalariaRDTNoAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 1, malaria_rdt == 3) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))


pNoSevereMalariaRDTAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 1, malaria_rdt == 3) %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))

  
pNoURSNoAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 0, resp_symptoms == "") %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))

pNoURSAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 0, resp_symptoms == "") %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))

pURSNoAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 0, resp_symptoms != "") %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 0, na.rm = TRUE)))

pURSAB_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 0, resp_symptoms != "") %>%
  summarise(proportion = ifelse(n() == 0, NA, mean(primary_outcome == 1, na.rm = TRUE)))

##ANTIBIOTICS BY TYPE PROBABILITIES: INTERVENTION
#This is asking: 'of the antibiotics prescribed in each intervention branch diagnosis result node, 
#what types of antibiotics were they?' This is needed for weighted cost inputs

#No one had severe malaria so it is not included

AB_WithoutSevereMalaria_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 1, malaria_rdt == 2) %>%   # added malaria_rdt == 2
  mutate(which_antibiotics = ifelse(which_antibiotics == "", "None", which_antibiotics)) %>%
  group_by(which_antibiotics) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

as.data.frame(AB_WithoutSevereMalaria_Intervention)  # View the full data frame in console

AB_WithoutSevereMalariaRDT_Intervention <- clinicaldata %>%
  filter(arm == 1, dangersigns == 1, malaria_rdt == 3) %>%   # added malaria_rdt == 2
  mutate(which_antibiotics = ifelse(which_antibiotics == "", "None", which_antibiotics)) %>%
  group_by(which_antibiotics) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

as.data.frame(AB_WithoutSevereMalariaRDT_Intervention)  # View the full data frame in console

AB_WithURS_Intervention <- clinicaldata %>%
filter(arm == 1, dangersigns == 0, resp_symptoms != "") %>% 
  mutate(which_antibiotics = ifelse(which_antibiotics == "", "None", which_antibiotics)) %>%
  group_by(which_antibiotics) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

as.data.frame(AB_WithURS_Intervention)  # View the full data frame in console

AB_WithoutURS_Intervention <- clinicaldata %>%
filter(arm == 1, dangersigns == 0, resp_symptoms == "") %>% 
  mutate(which_antibiotics = ifelse(which_antibiotics == "", "None", which_antibiotics)) %>%
  group_by(which_antibiotics) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

as.data.frame(AB_WithoutURS_Intervention)  # View the full data frame in console


##REFERRAL PROBABILITIES: INTERVENTION 
pReferredH_Intervention <- clinicaldata %>%
  filter(arm == 1) %>%
  summarise(percentage = mean(hosp_referral == 1, na.rm = TRUE))

pHospitalRef_SevereMalariaNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 1,
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_SevereMalariaAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 1,
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_NoSevereMalariaNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 2,
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_NoSevereMalariaAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 2,
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_NoSevereMalariaRDTNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 3,
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_NoSevereMalariaRDTAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 3,
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_NoURSNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms == "",
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_NoURSAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms == "",
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pHospitalRef_URSNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms != "",
    primary_outcome == 0,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pHospitalRef_URSAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms != "",
    primary_outcome == 1,
    hosp_referral %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(hosp_referral == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


##RECOVERY AT FOLLOW UP PROBABILITIES: INTERVENTION
pRecovered_SevereMalariaNoAB_Intervention <- clinicaldata %>% #nobody had severe malaria
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 1,
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_SevereMalariaAB_Intervention <- clinicaldata %>% #nobody had severe malaria
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 1,
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoSevereMalariaNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 2,
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoSevereMalariaAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 2,
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoSevereMalariaRDTNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 3,
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoSevereMalariaRDTAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 3,
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoURSNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms == "",
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_NoURSAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms == "",
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_URSNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms != "",
    primary_outcome == 0,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)

pRecovered_URSAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms != "",
    primary_outcome == 1,
    d7_patient_recovered %in% c(0, 1)  # Exclude NAs and 2s
  ) %>%
  summarise(proportion_recovered = mean(d7_patient_recovered == 1)) %>%
  pull(proportion_recovered)


##FOLLOW UP -- UNSCHEDULED REPRESENTING AT A PHC PROBABILITIES: INTERVENTION
pRepresent_SevereMalariaNoAB_Intervention <- clinicaldata %>% #nobody had severe malaria
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 1,
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_SevereMalariaAB_Intervention <- clinicaldata %>% #nobody had severe malaria
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 1,
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_NoSevereMalariaNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 2,
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_NoSevereMalariaAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 2,
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pRepresent_NoSevereMalariaRDTNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 3,
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_NoSevereMalariaRDTAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 1,
    malaria_rdt == 3,
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)


pRepresent_NoURSNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms == "",
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pRepresent_NoURSAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms == "",
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pRepresent_URSNoAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms != "",
    primary_outcome == 0,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

pRepresent_URSAB_Intervention <- clinicaldata %>%
  filter(
    arm == 1,
    dangersigns == 0,
    resp_symptoms != "",
    primary_outcome == 1,
    d7_patient_recovered == 0,
    d7_represent_elsewhere %in% c(0, 1) | d7_represent %in% c(0, 1)  # Include only valid binary values
  ) %>%
  summarise(
    proportion = if (n() > 0) {
      mean(d7_represent_elsewhere == 1 | d7_represent == 1, na.rm = TRUE)
    } else {
      0
    }
  ) %>%
  pull(proportion)

##FOLLOW UP -- ADMITTED TO HOSPITAL: BOTH ARMS
pFollowUp_HospAdmit <- clinicaldata %>%
  summarise(percentage = mean(d7_admitted_to_hospital == 1, na.rm = TRUE))

##A list of the base case probabilities for use in the PSA
base_probs <- list(
  pMalaria_Control = pMalaria_Control,
  pNoMalaria_Control = pNoMalaria_Control,
  pNoMalariaRDT_Control = pNoMalariaRDT_Control,
  pMalariaAB_Control = pMalariaAB_Control,
  pMalariaNoAB_Control = pMalariaNoAB_Control,
  pNoMalariaAB_Control = pNoMalariaAB_Control,
  pNoMalariaNoAB_Control = pNoMalariaNoAB_Control,
  pNoMalariaRDTAB_Control = pNoMalariaRDTAB_Control,
  pNoMalariaRDTNoAB_Control = pNoMalariaRDTNoAB_Control,
  pRecovered_MalariaAB_Control = pRecovered_MalariaAB_Control,
  pRecovered_MalariaNoAB_Control = pRecovered_MalariaNoAB_Control,
  pRecovered_NoMalariaAB_Control = pRecovered_NoMalariaAB_Control,
  pRecovered_NoMalariaNoAB_Control = pRecovered_NoMalariaNoAB_Control,
  pRecovered_NoMalariaRDTAB_Control = pRecovered_NoMalariaRDTAB_Control,
  pRecovered_NoMalariaRDTNoAB_Control = pRecovered_NoMalariaRDTNoAB_Control,
  pDangerSigns_Intervention = pDangerSigns_Intervention,
  pSevereMalaria_Intervention = pSevereMalaria_Intervention,
  pNoSevereMalaria_Intervention = pNoSevereMalaria_Intervention,
  pNoSevereMalariaRDT_Intervention = pNoSevereMalariaRDT_Intervention,
  pSevereMalariaAB_Intervention = pSevereMalariaAB_Intervention,
  pSevereMalariaNoAB_Intervention = pSevereMalariaNoAB_Intervention,
  pNoSevereMalariaAB_Intervention = pNoSevereMalariaAB_Intervention,
  pNoSevereMalariaNoAB_Intervention = pNoSevereMalariaNoAB_Intervention,
  pNoSevereMalariaRDTAB_Intervention = pNoSevereMalariaRDTAB_Intervention,
  pNoSevereMalariaRDTNoAB_Intervention = pNoSevereMalariaRDTNoAB_Intervention,
  pNoURSymptoms_Intervention = pNoURSymptoms_Intervention,
  pNoURSAB_Intervention = pNoURSAB_Intervention,
  pNoURSNoAB_Intervention = pNoURSNoAB_Intervention,
  pURSymptoms_Intervention = pURSymptoms_Intervention,
  pURSAB_Intervention = pURSAB_Intervention,
  pURSNoAB_Intervention = pURSNoAB_Intervention,
  pRecovered_SevereMalariaAB_Intervention = pRecovered_SevereMalariaAB_Intervention,
  pRecovered_SevereMalariaNoAB_Intervention = pRecovered_SevereMalariaNoAB_Intervention,
  pRecovered_NoSevereMalariaAB_Intervention = pRecovered_NoSevereMalariaAB_Intervention,
  pRecovered_NoSevereMalariaNoAB_Intervention = pRecovered_NoSevereMalariaNoAB_Intervention,
  pRecovered_NoSevereMalariaRDTAB_Intervention = pRecovered_NoSevereMalariaRDTAB_Intervention,
  pRecovered_NoSevereMalariaRDTNoAB_Intervention = pRecovered_NoSevereMalariaRDTNoAB_Intervention,
  pRecovered_NoURSAB_Intervention = pRecovered_NoURSAB_Intervention,
  pRecovered_NoURSNoAB_Intervention = pRecovered_NoURSNoAB_Intervention,
  pRecovered_URSAB_Intervention = pRecovered_URSAB_Intervention,
  pRecovered_URSNoAB_Intervention = pRecovered_URSNoAB_Intervention
)

base_probs <- lapply(base_probs, function(x) as.numeric(x))
################################################################################
##PARAMETERS
################################################################################
##1: HEALTHCARE SYSTEM PERSPECTIVE COSTS
##1.1: ANTIBIOTIC COSTS

amoxicillinKHR <- 486                    #Cost in Cambodian Riels
amoxicillinUSD <- amoxicillinKHR*0.00025 #Cost in US Dollars
amoxicillin <- amoxicillinUSD * 5        #Cost of course -- multiplied by avg. 
                                         #number of days prescribed for PHCs

penicillinKHR <- 500
penicillinUSD <- penicillinKHR*0.00025
penicillin <- penicillinUSD * 5

ciprofloxacinKHR <- 498.36
ciprofloxacinUSD <- ciprofloxacinKHR*0.00025 
ciprofloxacin <- ciprofloxacinUSD * 7

cotrimoxazoleKHR <- 104.47 
cotrimoxazoleUSD <- cotrimoxazoleKHR*0.00025
cotrimoxazole <- cotrimoxazoleUSD * 5

cloxacillinKHR <- 433.6 
cloxacillinUSD <- cloxacillinKHR*0.00025
cloxacillin <- cloxacillinUSD * 5

metronidazoleUSD <- 0.015 #metronidazole KHR data missing
metronidazole <- metronidazoleUSD * 7

#Hard coding the costs calculated above into a list, to be used below
Antibiotic_Costs <- c(amoxicillin = 0.607, penicillin = 0.625, ciprofloxacin = 0.872,
                      cotrimoxazole = 0.130, cloxacillin = 0.542, metronidazole = 0.105)

##WEIGHTED AB COSTS: CONTROL ARM
#1. AB Quantities: These are hard coded from the 'AB by Type' probability data frames creates earlier
ABQuantities_WithMalaria_Control <- c(amoxicillin = 7, penicillin = 0, ciprofloxacin = 0, 
                                     cotrimoxazole = 1, cloxacillin = 0, metronidazole = 0)

#2. Finding a weighted, 'typical' cost of receiving an antibiotic in the given diagnosis node
ABWeightedCosts_WithMalaria_Control <- sum(Antibiotic_Costs * ABQuantities_WithMalaria_Control) / sum(ABQuantities_WithMalaria_Control)


ABQuantities_WithoutMalaria_Control <- c(amoxicillin = 657, penicillin = 100, ciprofloxacin = 42, 
                                         cotrimoxazole = 113, cloxacillin = 10, metronidazole = 7)

ABWeightedCosts_WithoutMalaria_Control <- sum(Antibiotic_Costs * ABQuantities_WithoutMalaria_Control) / sum(ABQuantities_WithoutMalaria_Control)

ABQuantities_WithoutMalariaRDT_Control <- c(amoxicillin = 397, penicillin = 21, ciprofloxacin = 30, 
                                            cotrimoxazole = 37, cloxacillin = 4, metronidazole = 4)

ABWeightedCosts_WithoutMalariaRDT_Control <- sum(Antibiotic_Costs * ABQuantities_WithoutMalariaRDT_Control) / sum(ABQuantities_WithoutMalariaRDT_Control)

ABTotal_Control <- ABQuantities_WithMalaria_Control + ABQuantities_WithoutMalaria_Control + ABQuantities_WithoutMalariaRDT_Control

print(ABTotal_Control)

#The above values are used to estimate antibiotic costs as accurately as possible. 
#However, to obtain the final number of antibiotic treated patients, 
#we must double check for potential antibiotics prescribed 
#but with no type recorded. True total number below:
sum(clinicaldata$arm == 0 & clinicaldata$primary_outcome == 1, na.rm = TRUE)
ABTotal_Control_Num <-1468


##WEIGHTED AB COSTS: INTERVENTION ARM
##Nobody had severe malaria

ABQuantities_WithoutSevereMalaria_Intervention <- c(amoxicillin = 456, penicillin = 10, ciprofloxacin = 16, 
                                                    cotrimoxazole = 47, cloxacillin = 21, metronidazole = 2)


ABWeightedCosts_WithoutSevereMalaria_Intervention <- sum(Antibiotic_Costs * ABQuantities_WithoutSevereMalaria_Intervention) / sum(ABQuantities_WithoutSevereMalaria_Intervention)


ABQuantities_WithoutSevereMalariaRDT_Intervention <- c(amoxicillin = 173, penicillin = 1, ciprofloxacin = 2, 
                                                      cotrimoxazole = 11, cloxacillin = 5, metronidazole = 1)

ABWeightedCosts_WithoutSevereMalariaRDT_Intervention <- sum(Antibiotic_Costs * ABQuantities_WithoutSevereMalariaRDT_Intervention) / sum(ABQuantities_WithoutSevereMalariaRDT_Intervention)

ABQuantities_NoURS_Intervention <- c(amoxicillin = 16, penicillin = 2, ciprofloxacin = 3, 
                                    cotrimoxazole = 2, cloxacillin = 9, metronidazole = 1)
ABWeightedCosts_NoURS_Intervention <- sum(Antibiotic_Costs * ABQuantities_NoURS_Intervention) / sum(ABQuantities_NoURS_Intervention)

ABQuantities_URS_Intervention <- c(amoxicillin = 449, penicillin = 11, ciprofloxacin = 23, 
                                  cotrimoxazole = 22, cloxacillin = 6, metronidazole = 3)
ABWeightedCosts_URS_Intervention <- sum(Antibiotic_Costs * ABQuantities_URS_Intervention) / sum(ABQuantities_URS_Intervention)
 
ABTotal_Intervention <- ABQuantities_WithoutSevereMalaria_Intervention + ABQuantities_WithoutSevereMalariaRDT_Intervention + ABQuantities_NoURS_Intervention + ABQuantities_URS_Intervention
print(ABTotal_Intervention)

#The above values are used to estimate antibiotic costs as accurately as possible. 
#However, to obtain the final number of antibiotic treated patients, 
#we must double check for potential antibiotics prescribed 
#but with no type recorded. True total number below:
sum(clinicaldata$arm == 1 & clinicaldata$primary_outcome == 1, na.rm = TRUE)
ABTotal_Intervention_Num <-1304



##1.2: COSTS IN FOLLOW UP PERIOD
cambodia_inflation <- 148.3/100 #using CPI ratio from 2010 to 2023 (most recent)
cFollowUp <- 6.19 * cambodia_inflation #inflate from WHO 2010 price estimate

cInpatientPPP <- 22.02   #Average inpatient cost in Cambodia in PPP in 2010
PPPConversion <- 1460.39 #Cambodia specific PPP conversion in 2010

cHospital_Admission <- (cInpatientPPP * PPPConversion * cambodia_inflation) * 0.00025 #2010 PPP --> 2010 KHR --> ~present KHR --> present USD
cRepresentPHC <- 6.19 * cambodia_inflation 


##1.3: CONTROL ARM SPECIFIC COSTS
cRDT <- 1.1183 #Malaria RDT
cAppt_Control <- 6.19 * cambodia_inflation #Assume standard primary care outpatient bed-day cost

##1.4: INTERVENTION ARM SPECIFIC COST
cAppt_Intervention <- cAppt_Control * (durationpercent_increase/100) #Liberal estimate using median appointment duration % increase from time calculations
cEDAM <- ((6000 * 0.031) + 3.5) / (1 * 2 * 365) #Intervention per use cost: 6000 THB tablet --> tablet in USD + $3.5 USD charger --> used by 2 patients per day for 1 year
cCRP <- 0.69 #C-Reactive Protein Test

##1.4.1: VITAL SIGN ASSESSMENT COSTS: CONSERVATIVE ESTIMATES 
cPulseox <- 349.5 / (2 * 2 * 365) #Pulse ox per use cost: 2 year life span, 2 uses per day
cBloodpressure <- 56.75 / (2 * 2 * 365) #Blood pressure machine per use cost: 2 year life span, 2 uses per day

##1.5: WILLINGNESS TO PAY THRESHOLD
#This is the Cambodia GDP per capita times a proportion established by Pichon_Riviere (2023)
Cambodia_WTP <- 2429.75 * 0.52 

##2: SOCIETAL COSTS
##2.1: ANTIMICROBIAL RESISTANCE (AMR) COSTS 
##Convert costs (in Thai Baht) from Shrestha (2018) to present day costs
thai_inflation <- 122.1/112.5 #using CPI ratio from 2018-2023 (most recent)
thai_aminoglycoside <- 12.33 * thai_inflation
thai_BSP <- 10.26 * thai_inflation
thai_carbapenem <- 12.64 * thai_inflation
thai_cephalosporins <- 9.58 * thai_inflation
thai_macrolides <- 0.29 * thai_inflation
thai_NSP <- 2.71 * thai_inflation
thai_quinolones <- 19.15 * thai_inflation

#Now convert to present day costs in Cambodia 
thai_to_cambodia_GDP <- 2429.7/7182  ##GDP adjustment factor, useful when comparing between LMICs
cambodia_aminoglycoside <- thai_aminoglycoside * thai_to_cambodia_GDP
cambodia_BSP <- thai_BSP * thai_to_cambodia_GDP
cambodia_carbapenem <- thai_carbapenem * thai_to_cambodia_GDP
cambodia_cephalosporins <- thai_cephalosporins * thai_to_cambodia_GDP
cambodia_macrolides <- thai_macrolides * thai_to_cambodia_GDP
cambodia_NSP <- thai_NSP * thai_to_cambodia_GDP
cambodia_quinolones <- thai_quinolones * thai_to_cambodia_GDP

#See table 2.2 footnotes in main report for antibiotic class assumptions
AMR_Costs <- c(amoxicillin = cambodia_NSP, penicillin = cambodia_NSP, ciprofloxacin = cambodia_quinolones,
cotrimoxazole = cambodia_quinolones, cloxacillin = cambodia_NSP, 
metronidazole = cambodia_NSP)

print(AMR_Costs)

#1. AMR Quantities: These are hard coded from the 'AB by Type' probability data frames creates earlier
AMRQuantities_WithMalaria_Control <- c(amoxicillin = 7, penicillin = 0, ciprofloxacin = 0, 
                                    cotrimoxazole = 1, cloxacillin = 0, metronidazole = 0)

#2. Finding a weighted, 'typical' cost of AMR by antibiotic consumption in the given diagnosis node
AMRWeightedCosts_WithMalaria_Control <- sum(AMR_Costs * AMRQuantities_WithMalaria_Control) / sum(AMRQuantities_WithMalaria_Control)

AMRQuantities_WithoutMalaria_Control <- c(amoxicillin = 657, penicillin = 100, ciprofloxacin = 42, 
                                          cotrimoxazole = 113, cloxacillin = 10, metronidazole = 7)
AMRWeightedCosts_WithoutMalaria_Control <- sum(AMR_Costs * AMRQuantities_WithoutMalaria_Control) / sum(AMRQuantities_WithoutMalaria_Control)

AMRQuantities_WithoutMalariaRDT_Control <- c(amoxicillin = 397, penicillin = 21, ciprofloxacin = 30, 
                                             cotrimoxazole = 37, cloxacillin = 4, metronidazole = 4)
AMRWeightedCosts_WithoutMalariaRDT_Control <- sum(AMR_Costs * AMRQuantities_WithoutMalariaRDT_Control) / sum(AMRQuantities_WithoutMalariaRDT_Control)

AMRCostsTotal_Control <- mean(AMRWeightedCosts_WithMalaria_Control, AMRWeightedCosts_WithoutMalaria_Control, AMRWeightedCosts_WithoutMalariaRDT_Control)

##WEIGHTED AMR COSTS: INTERVENTION ARM
#Nobody had severe malaria

AMRQuantities_WithoutSevereMalaria_Intervention <- c(amoxicillin = 456, penicillin = 10, ciprofloxacin = 16, 
                                                     cotrimoxazole = 47, cloxacillin = 21, metronidazole = 2)
AMRWeightedCosts_WithoutSevereMalaria_Intervention <- sum(AMR_Costs * AMRQuantities_WithoutSevereMalaria_Intervention) / sum(AMRQuantities_WithoutSevereMalaria_Intervention)

AMRQuantities_WithoutSevereMalariaRDT_Intervention <- c(amoxicillin = 173, penicillin = 1, ciprofloxacin = 2, 
                                                        cotrimoxazole = 11, cloxacillin = 5, metronidazole = 1)
AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention <- sum(AMR_Costs * AMRQuantities_WithoutSevereMalariaRDT_Intervention) / sum(AMRQuantities_WithoutSevereMalariaRDT_Intervention)

AMRQuantities_NoURS_Intervention <- c(amoxicillin = 16, penicillin = 2, ciprofloxacin = 3, 
                                      cotrimoxazole = 2, cloxacillin = 9, metronidazole = 1)
AMRWeightedCosts_NoURS_Intervention <- sum(AMR_Costs * AMRQuantities_NoURS_Intervention) / sum(AMRQuantities_NoURS_Intervention)

AMRQuantities_URS_Intervention <- c(amoxicillin = 449, penicillin = 11, ciprofloxacin = 23, 
                                    cotrimoxazole = 22, cloxacillin = 6, metronidazole = 3)
AMRWeightedCosts_URS_Intervention <- sum(AMR_Costs * AMRQuantities_URS_Intervention) / sum(AMRQuantities_URS_Intervention)

AMRCostsTotal_Intervention <- mean(AMRWeightedCosts_WithoutSevereMalaria_Intervention, AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention, AMRWeightedCosts_NoURS_Intervention, AMRWeightedCosts_URS_Intervention)

##2.2: PRODUCTIVITY LOSS COSTS
dailywage <- 2429.7/365 #GDP per capita/ 365 days in a year
employmentrate <- .7880 #From Consensus

#See main report section 2.4 for days off work assumptions 
cProductivity_failure <- dailywage * 7 * employmentrate #cost of 7 working days lost
cProductivity_success <- dailywage * 3.5 * employmentrate #cost of only 3.5 working days lost


##A list of the most important basecase costs for use in the PSA
base_costs <- list(
  cRDT = cRDT, 
  cCRP = cCRP,
  cAppt_Control = cAppt_Control,
  cAppt_Intervention = cAppt_Intervention,
  cFollowUp = cFollowUp,
  cRepresentPHC = cRepresentPHC, 
  cHospital_Admission = cHospital_Admission,
  cEDAM = cEDAM,
  cPulseox = cPulseox, 
  cBloodpressure = cBloodpressure,
  ABWeightedCosts_WithMalaria_Control = ABWeightedCosts_WithMalaria_Control,
  ABWeightedCosts_WithoutMalaria_Control = ABWeightedCosts_WithoutMalaria_Control,
  ABWeightedCosts_WithoutMalariaRDT_Control = ABWeightedCosts_WithoutMalariaRDT_Control,
  ABWeightedCosts_WithoutSevereMalaria_Intervention = ABWeightedCosts_WithoutSevereMalaria_Intervention,
  ABWeightedCosts_WithoutSevereMalariaRDT_Intervention = ABWeightedCosts_WithoutSevereMalariaRDT_Intervention,
  ABWeightedCosts_URS_Intervention = ABWeightedCosts_URS_Intervention,
  ABWeightedCosts_NoURS_Intervention = ABWeightedCosts_NoURS_Intervention,
  AMRWeightedCosts_WithMalaria_Control = AMRWeightedCosts_WithMalaria_Control, 
  AMRWeightedCosts_WithoutMalaria_Control = AMRWeightedCosts_WithoutMalaria_Control,
  AMRWeightedCosts_WithoutMalariaRDT_Control = AMRWeightedCosts_WithoutMalariaRDT_Control,
  AMRWeightedCosts_WithoutSevereMalaria_Intervention = AMRWeightedCosts_WithoutSevereMalaria_Intervention,
  AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention = AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention,
  AMRWeightedCosts_URS_Intervention = AMRWeightedCosts_URS_Intervention,
  AMRWeightedCosts_NoURS_Intervention = AMRWeightedCosts_NoURS_Intervention
)
################################################################################
##PARAMETERS: HEALTH OUTCOME VALUES
##HEALTH OUTCOMES
moderate_disabilityweight <- 0.051 #From Global Burden of Disease Study

timesick_treatmentsuccess <- 3.5 
  #assume patients who recovered, did so at an average between 1-7 days
timesick_treatmentfailure <- 7 
  #patients not recovered are sick for the full 7 days

##Calculate DALYs
#See Figure 2.2 in main report for DALY calculations
DALY_moderate_success <- moderate_disabilityweight * (timesick_treatmentsuccess/365)
DALY_moderate_failure <- moderate_disabilityweight * (timesick_treatmentfailure/365)


################################################################################
##TERMINAL NODE VALUES
#Summation of costs for each terminal node, prior to building the decision tree
#See Supplemental Information S1 for more explanation on costs included in each node

perspective <- 'societal' #or healthcare

##TERMINAL NODE VALUES: CONTROL ARM
##Nodes: Malaria Positive, No AB
cost1 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control)
health1 <- DALY_moderate_success
cost2 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control) 
health2 <- DALY_moderate_failure
##Nodes: Malaria Positive, AB
cost3 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_success + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * pHospitalRef_MalariaAB_Control)
health3 <- DALY_moderate_success
cost4 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control) 
health4 <- DALY_moderate_failure

##Nodes: Malaria Negative, No AB
cost5 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control)
health5 <- DALY_moderate_success
cost6 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control)
health6 <- DALY_moderate_failure
##Nodes: Malaria Negative, AB
cost7 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_success + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control)
health7 <- DALY_moderate_success
cost8 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control)
health8 <- DALY_moderate_failure

##Nodes: No Malaria RDT, No AB
cost9 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control) + cProductivity_success else cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control)
health9 <- DALY_moderate_success
cost10 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control) + cProductivity_failure else cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control)
health10 <- DALY_moderate_failure
##Nodes: No Malaria RDT, AB
cost11 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_success + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control)
health11 <- DALY_moderate_success
cost12 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control) 
health12 <- DALY_moderate_failure


##TERMINAL NODE VALUES: INTERVENTION ARM, DANGER SIGNS PRESENT 
##Nodes: Severe Malaria, No AB (Nobody had severe malaria)
cost13 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention)
health13 <- DALY_moderate_success
cost14 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention)
health14 <- DALY_moderate_failure
##Nodes: Severe Malaria, AB (Nobody had severe malaria, so not adding AB costs)
cost15 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention)
health15 <- DALY_moderate_success
cost16 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention)
health16 <- DALY_moderate_failure

##Nodes: No Severe Malaria, No AB
cost17 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention)
health17 <- DALY_moderate_success
cost18 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention)
health18 <- DALY_moderate_failure
##Nodes: No Severe Malaria, AB
cost19 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention)
health19 <- DALY_moderate_success
cost20 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention)
health20 <- DALY_moderate_failure

#Nodes: No Severe Malaria RDT, No AB
cost21 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention)
health21 <- DALY_moderate_success
cost22 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention)
health22 <- DALY_moderate_failure
##Nodes: No Severe Malaria RDT, AB
cost23 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention)
health23 <- DALY_moderate_success
cost24 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention)
health24 <- DALY_moderate_failure


##TERMINAL NODE VALUES: INTERVENTION ARM, NO DANGER SIGNS PRESENT
##Nodes: UR Symptoms Present, No AB
cost25 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention)
health25 <- DALY_moderate_success
cost26 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention)
health26 <- DALY_moderate_failure
##Nodes: UR Symptoms Present, AB
cost27 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_success + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * pHospitalRef_URSAB_Intervention)
health27 <- DALY_moderate_success
cost28 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_failure + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention)
health28 <- DALY_moderate_failure

##Nodes: UR Symptoms Not Present, No AB
cost29 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention)
health29 <- DALY_moderate_success
cost30 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention)
health30 <- DALY_moderate_failure
##Nodes: UR Symptoms Not Present, AB
cost31 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_success + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention) 
health31 <- DALY_moderate_success
cost32 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_failure + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention)
health32 <- DALY_moderate_failure
################################################################################
##MODELLING THE DECISION TREE
################################################################################
##ROOT NODE
library(data.tree)

#In case of any confusion in data.tree between meaning of prob versus p, keep these synonymous: 
#Helper function to add child with both prob and p set
AddNodeWithProb <- function(parent, name, prob = NULL) {
  child <- parent$AddChild(name)
  if (!is.null(prob)) {
    child$prob <- prob
    child$p <- prob
  }
  return(child)
}

#Start the tree at the Root Node
tree <- Node$new("Trial Randomisation")

################################################################################
## CONTROL BRANCHES
control <- AddNodeWithProb(tree, "Control", prob = 0.5)

#'prob' will from here on pull a probability value calculated above in the data cleaning stage

## MALARIA RDT RESULT BRANCHES - CONTROL
branchC1 <- AddNodeWithProb(control, "Malaria Positive", prob = pMalaria_Control)
branchC2 <- AddNodeWithProb(control, "Malaria Negative", prob = pNoMalaria_Control)
branchC3 <- AddNodeWithProb(control, "No Malaria RDT", prob = pNoMalariaRDT_Control)

## PRIMARY OUTCOME BRANCHES - CONTROL
branchC1.1 <- AddNodeWithProb(branchC1, "No AB", prob = pMalariaNoAB_Control)
branchC1.2 <- AddNodeWithProb(branchC1, "AB", prob = pMalariaAB_Control)
branchC2.1 <- AddNodeWithProb(branchC2, "No AB", prob = pNoMalariaNoAB_Control)
branchC2.2 <- AddNodeWithProb(branchC2, "AB", prob = pNoMalariaAB_Control)
branchC3.1 <- AddNodeWithProb(branchC3, "No AB", prob = pNoMalariaRDTNoAB_Control)
branchC3.2 <- AddNodeWithProb(branchC3, "AB", prob = pNoMalariaRDTAB_Control)

## TREATMENT SUCCESS BRANCHES & TERMINAL NODES - CONTROL
branchC1.1.1 <- AddNodeWithProb(branchC1.1, "Treatment Success", prob = pRecovered_MalariaNoAB_Control)
branchC1.1.1$cost <- cost1
branchC1.1.1$daly <- health1

branchC1.1.2 <- AddNodeWithProb(branchC1.1, "Treatment Failure", prob = 1 - pRecovered_MalariaNoAB_Control)
branchC1.1.2$cost <- cost2
branchC1.1.2$daly <- health2

branchC1.2.1 <- AddNodeWithProb(branchC1.2, "Treatment Success", prob = pRecovered_MalariaAB_Control)
branchC1.2.1$cost <- cost3
branchC1.2.1$daly <- health3

branchC1.2.2 <- AddNodeWithProb(branchC1.2, "Treatment Failure", prob = 1 - pRecovered_MalariaAB_Control)
branchC1.2.2$cost <- cost4
branchC1.2.2$daly <- health4

branchC2.1.1 <- AddNodeWithProb(branchC2.1, "Treatment Success", prob = pRecovered_NoMalariaNoAB_Control)
branchC2.1.1$cost <- cost5
branchC2.1.1$daly <- health5

branchC2.1.2 <- AddNodeWithProb(branchC2.1, "Treatment Failure", prob = 1 - pRecovered_NoMalariaNoAB_Control)
branchC2.1.2$cost <- cost6
branchC2.1.2$daly <- health6

branchC2.2.1 <- AddNodeWithProb(branchC2.2, "Treatment Success", prob = pRecovered_NoMalariaAB_Control)
branchC2.2.1$cost <- cost7
branchC2.2.1$daly <- health7

branchC2.2.2 <- AddNodeWithProb(branchC2.2, "Treatment Failure", prob = 1 - pRecovered_NoMalariaAB_Control)
branchC2.2.2$cost <- cost8
branchC2.2.2$daly <- health8

branchC3.1.1 <- AddNodeWithProb(branchC3.1, "Treatment Success", prob = pRecovered_NoMalariaRDTNoAB_Control)
branchC3.1.1$cost <- cost9
branchC3.1.1$daly <- health9

branchC3.1.2 <- AddNodeWithProb(branchC3.1, "Treatment Failure", prob = 1 - pRecovered_NoMalariaRDTNoAB_Control)
branchC3.1.2$cost <- cost10
branchC3.1.2$daly <- health10

branchC3.2.1 <- AddNodeWithProb(branchC3.2, "Treatment Success", prob = pRecovered_NoMalariaRDTAB_Control)
branchC3.2.1$cost <- cost11
branchC3.2.1$daly <- health11

branchC3.2.2 <- AddNodeWithProb(branchC3.2, "Treatment Failure", prob = 1 - pRecovered_NoMalariaRDTAB_Control)
branchC3.2.2$cost <- cost12
branchC3.2.2$daly <- health12

################################################################################
## INTERVENTION BRANCHES
intervention <- AddNodeWithProb(tree, "Intervention", prob = 0.5)

## VITAL SIGNS ASSESSMENT RESULTS - INTERVENTION
branchI1 <- AddNodeWithProb(intervention, "Danger Signs Present", prob = pDangerSigns_Intervention)
branchI2 <- AddNodeWithProb(intervention, "No Danger Signs Present", prob = pNoDangerSigns_Intervention)

## DANGER SIGNS PRESENT - MALARIA RDT BRANCHES - INTERVENTION
branchI1.1 <- AddNodeWithProb(branchI1, "Malaria Positive", prob = 0) #nobody had severe malaria
branchI1.2 <- AddNodeWithProb(branchI1, "Malaria Negative", prob = pNoSevereMalaria_Intervention)
branchI1.3 <- AddNodeWithProb(branchI1, "No Malaria RDT", prob = pNoSevereMalariaRDT_Intervention)

## DANGER SIGNS PRESENT - PRIMARY OUTCOME BRANCHES - INTERVENTION
branchI1.1.1 <- AddNodeWithProb(branchI1.1, "No AB", prob = 0)
branchI1.1.2 <- AddNodeWithProb(branchI1.1, "AB", prob = 0)
branchI1.2.1 <- AddNodeWithProb(branchI1.2, "No AB", prob = pNoSevereMalariaNoAB_Intervention)
branchI1.2.2 <- AddNodeWithProb(branchI1.2, "AB", prob = pNoSevereMalariaAB_Intervention)
branchI1.3.1 <- AddNodeWithProb(branchI1.3, "No AB", prob = pNoSevereMalariaRDTNoAB_Intervention)
branchI1.3.2 <- AddNodeWithProb(branchI1.3, "AB", prob = pNoSevereMalariaRDTAB_Intervention)

## DANGER SIGNS PRESENT - TREATMENT SUCCESS BRANCHES & TERMINAL NODES - INTERVENTION
branchI1.1.1.1 <- AddNodeWithProb(branchI1.1.1, "Treatment Success", prob = 0)
branchI1.1.1.1$cost <- cost13
branchI1.1.1.1$daly <- health13

branchI1.1.1.2 <- AddNodeWithProb(branchI1.1.1, "Treatment Failure", prob = 0)
branchI1.1.1.2$cost <- cost14
branchI1.1.1.2$daly <- health14

branchI1.1.2.1 <- AddNodeWithProb(branchI1.1.2, "Treatment Success", prob = 0)
branchI1.1.2.1$cost <- cost15
branchI1.1.2.1$daly <- health15

branchI1.1.2.2 <- AddNodeWithProb(branchI1.1.2, "Treatment Failure", prob = 0)
branchI1.1.2.2$cost <- cost16
branchI1.1.2.2$daly <- health16

branchI1.2.1.1 <- AddNodeWithProb(branchI1.2.1, "Treatment Success", prob = pRecovered_NoSevereMalariaNoAB_Intervention)
branchI1.2.1.1$cost <- cost17
branchI1.2.1.1$daly <- health17

branchI1.2.1.2 <- AddNodeWithProb(branchI1.2.1, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaNoAB_Intervention)
branchI1.2.1.2$cost <- cost18
branchI1.2.1.2$daly <- health18

branchI1.2.2.1 <- AddNodeWithProb(branchI1.2.2, "Treatment Success", prob = pRecovered_NoSevereMalariaAB_Intervention)
branchI1.2.2.1$cost <- cost19
branchI1.2.2.1$daly <- health19

branchI1.2.2.2 <- AddNodeWithProb(branchI1.2.2, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaAB_Intervention)
branchI1.2.2.2$cost <- cost20
branchI1.2.2.2$daly <- health20

branchI1.3.1.1 <- AddNodeWithProb(branchI1.3.1, "Treatment Success", prob = pRecovered_NoSevereMalariaRDTNoAB_Intervention)
branchI1.3.1.1$cost <- cost21
branchI1.3.1.1$daly <- health21

branchI1.3.1.2 <- AddNodeWithProb(branchI1.3.1, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaRDTNoAB_Intervention)
branchI1.3.1.2$cost <- cost22
branchI1.3.1.2$daly <- health22

branchI1.3.2.1 <- AddNodeWithProb(branchI1.3.2, "Treatment Success", prob = pRecovered_NoSevereMalariaRDTAB_Intervention)
branchI1.3.2.1$cost <- cost23
branchI1.3.2.1$daly <- health23

branchI1.3.2.2 <- AddNodeWithProb(branchI1.3.2, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaRDTAB_Intervention)
branchI1.3.2.2$cost <- cost24
branchI1.3.2.2$daly <- health24

## NO DANGER SIGNS PRESENT - SYMPTOMS ASSESSMENT TO RDT - INTERVENTION 
branchI2.1 <- AddNodeWithProb(branchI2, "Symptoms Assessment to Malaria RDT", prob = 1)

## NO DANGER SIGNS PRESENT - MALARIA RDT BRANCHES - INTERVENTION
branchI2.2 <- AddNodeWithProb(branchI2.1, "UR Symptoms", prob = pURSymptoms_Intervention)
branchI2.3 <- AddNodeWithProb(branchI2.1, "No UR Symptoms", prob = pNoURSymptoms_Intervention)

## NO DANGER SIGNS PRESENT - PRIMARY OUTCOME BRANCHES - INTERVENTION
branchI2.2.1 <- AddNodeWithProb(branchI2.2, "No AB", prob = pURSNoAB_Intervention)
branchI2.2.2 <- AddNodeWithProb(branchI2.2, "AB", prob = pURSAB_Intervention)
branchI2.3.1 <- AddNodeWithProb(branchI2.3, "No AB", prob = pNoURSNoAB_Intervention)
branchI2.3.2 <- AddNodeWithProb(branchI2.3, "AB", prob = pNoURSAB_Intervention)

## NO DANGER SIGNS PRESENT - TREATMENT SUCCESS BRANCHES & TERMINAL NODES - INTERVENTION
branchI2.2.1.1 <- AddNodeWithProb(branchI2.2.1, "Treatment Success", prob = pRecovered_URSNoAB_Intervention)
branchI2.2.1.1$cost <- cost25
branchI2.2.1.1$daly <- health25

branchI2.2.1.2 <- AddNodeWithProb(branchI2.2.1, "Treatment Failure", prob = 1 - pRecovered_URSNoAB_Intervention)
branchI2.2.1.2$cost <- cost26
branchI2.2.1.2$daly <- health26

branchI2.2.2.1 <- AddNodeWithProb(branchI2.2.2, "Treatment Success", prob = pRecovered_URSAB_Intervention)
branchI2.2.2.1$cost <- cost27
branchI2.2.2.1$daly <- health27

branchI2.2.2.2 <- AddNodeWithProb(branchI2.2.2, "Treatment Failure", prob = 1 - pRecovered_URSAB_Intervention)
branchI2.2.2.2$cost <- cost28
branchI2.2.2.2$daly <- health28

branchI2.3.1.1 <- AddNodeWithProb(branchI2.3.1, "Treatment Success", prob = pRecovered_NoURSNoAB_Intervention)
branchI2.3.1.1$cost <- cost29
branchI2.3.1.1$daly <- health29

branchI2.3.1.2 <- AddNodeWithProb(branchI2.3.1, "Treatment Failure", prob = 1 - pRecovered_NoURSNoAB_Intervention)
branchI2.3.1.2$cost <- cost30
branchI2.3.1.2$daly <- health30

branchI2.3.2.1 <- AddNodeWithProb(branchI2.3.2, "Treatment Success", prob = pRecovered_NoURSAB_Intervention)
branchI2.3.2.1$cost <- cost31
branchI2.3.2.1$daly <- health31

branchI2.3.2.2 <- AddNodeWithProb(branchI2.3.2, "Treatment Failure", prob = 1 - pRecovered_NoURSAB_Intervention)
branchI2.3.2.2$cost <- cost32
branchI2.3.2.2$daly <- health32
################################################################################
#Function to ensure prob and p are in fact the same value
VerifyProbs <- function(tree) {
  mismatches <- list()
  
  tree$Do(function(node) {
    if (!is.null(node$prob) && !is.null(node$p)) {
      if (!isTRUE(all.equal(node$prob, node$p))) {
        mismatches[[length(mismatches) + 1]] <- paste0("Mismatch at node '", node$name, 
                                                       "': prob = ", node$prob, 
                                                       ", p = ", node$p)
      }
    } else if (!is.null(node$prob) || !is.null(node$p)) {
      mismatches[[length(mismatches) + 1]] <- paste0("Missing prob or p at node '", node$name, "'")
    }
  })
  
  if (length(mismatches) == 0) {
    cat("✅ All nodes have matching `prob` and `p`.\n")
  } else {
    cat("⚠️ Found mismatches:\n")
    for (msg in mismatches) {
      cat(msg, "\n")
    }
  }
}
VerifyProbs(tree)

##FOLD TREE UP
CalculateExpectedValues <- function(node) {
  if (node$isLeaf) {
    # Terminal nodes
    node$expected_cost <- node$cost
    node$expected_daly <- node$daly
  } else {
    # Internal nodes
    total_cost <- 0
    total_daly <- 0
    
    # Loop through 
    for (child in node$children) {
      CalculateExpectedValues(child)
      prob <- child$p
      total_cost <- total_cost + prob * child$expected_cost
      total_daly <- total_daly + prob * child$expected_daly
    }
    
    # Assign expected values to current node
    node$expected_cost <- total_cost
    node$expected_daly <- total_daly
  }
}

##CALCULATE EXPECTED VALUES
CalculateExpectedValues(tree)
CalculateExpectedValues(tree$Control)
CalculateExpectedValues(tree$Intervention)


cost_control <- tree$Control$expected_cost
daly_control <- as.numeric(tree$Control$expected_daly)

cost_intervention <- tree$Intervention$expected_cost
daly_intervention <- as.numeric(tree$Intervention$expected_daly)


##CALCULATE ICER FOR BASECASE
delta_cost <- cost_intervention - cost_control
delta_daly <- daly_control - daly_intervention #DALYs averted!
prescriptions_averted <-ABTotal_Control_Num - ABTotal_Intervention_Num


basecase_ICER <- delta_cost / delta_daly
basecase_ICER_ABs <- delta_cost / prescriptions_averted

#Net Monetary Benefit
incremental_basecase_nmb <- Cambodia_WTP * delta_daly - delta_cost

##ICER Summary Results
is.list(cost_control)
is.list(daly_control)
is.list(ABTotal_Control_Num)

ICER_summary_results <- data.frame(
  Outcome = c("Cost", "DALY", "ICER ($/DALY)", "Antibiotics Averted", "ICER ($/AB)"),
  Control = c(
    round(as.numeric(cost_control), 2),
    round(as.numeric(daly_control), 5),
    NA,
    ABTotal_Control_Num,
    NA
  ),
  Intervention = c(
    round(as.numeric(cost_intervention), 2),
    round(as.numeric(daly_intervention), 5),
    NA,
    ABTotal_Intervention_Num,
    NA
  ),
  Delta = c(
    round(as.numeric(delta_cost), 2),
    round(as.numeric(delta_daly), 5),
    round(as.numeric(basecase_ICER), 2),
    prescriptions_averted,
    round(as.numeric(basecase_ICER_ABs), 2)
  ),
  stringsAsFactors = FALSE
)

options(scipen = 999)  # Prevent scientific notation

################################################################################

##PLOTTING TREE
SetGraphStyle(tree, direction = "LR")

##CUSTOMISE DECISION TREE PLOT
SetNodeStyle(
  tree,
  label = function(node) {
    if (node$isLeaf) {
      paste0(node$name, "\nCost: ", node$cost, "\nDALY: ", node$daly)
    } else {
      node$name 
    }
  },
  shape = "box",
  style = "filled",
  fontsize = function(node) if (node$isLeaf) 12 else 10,
  fillcolor = function(node) {
    if (node$isLeaf) {
      "lightgoldenrod1"  # Terminal nodes
    } else {
      "lightblue"  # Internal nodes
    }
  }
)

graph <- ToDiagrammeRGraph(tree)


render_graph(graph)

##Cost-Effectiveness Plane
library(ggplot2)

# Prepare data
ce_plane_data_deltas <- data.frame(
  strategy = c("Intervention"),
  delta_cost = as.numeric(delta_cost),
  delta_ab = as.numeric(prescriptions_averted),
  delta_daly = as.numeric(delta_daly)
)

max_abs_x <- max(abs(ce_plane_data_deltas$delta_daly)) * 1.1
max_abs_y <- max(abs(ce_plane_data_deltas$delta_cost)) * 1.1

# Clean CE plane
ggplot(ce_plane_data_deltas, aes(x = delta_daly, y = delta_cost)) +
  geom_hline(yintercept = 0, color = "grey80", linetype = "solid", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linetype = "solid", linewidth = 0.3) +
  
  # Main data point
  geom_point(color = "goldenrod", shape = 8, size = 4) +
  geom_text(aes(label = strategy), vjust = -1.5, hjust = 0.5, size = 3.5) +
  
  # WTP line
  geom_abline(slope = 1263.47, intercept = 0, linetype = "dashed", color = "peachpuff4", size = 1) +
  annotate("text", 
           x = 0.6 * max_abs_x, 
           y = 1263.47 * 0.6 * max_abs_x,
           label = "WTP Threshold", 
           color = "peachpuff4", 
           angle = 45,
           hjust = -0.1, 
           size = 3.2) +
  
  # Quadrant labels
  annotate("text", x = 0.75 * max_abs_x, y = 0.85 * max_abs_y,
           label = "More costly\nMore effective", size = 3, hjust = 0.5) +
  annotate("text", x = 0.75 * max_abs_x, y = -0.85 * max_abs_y,
           label = "Less costly\nMore effective", size = 3, hjust = 0.5) +
  annotate("text", x = -0.75 * max_abs_x, y = 0.85 * max_abs_y,
           label = "More costly\nLess effective", size = 3, hjust = 0.5) +
  annotate("text", x = -0.75 * max_abs_x, y = -0.85 * max_abs_y,
           label = "Less costly\nLess effective", size = 3, hjust = 0.5) +
  
  # Labels and axis formatting
  labs(
    x = "DALYs Averted",
    y = "Incremental Cost (USD)"
  ) +
  
  coord_cartesian(xlim = c(-max_abs_x, max_abs_x),
                  ylim = c(-max_abs_y, max_abs_y)) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, hjust = 0.5, face = "italic")
  )


##Zoomed out version of CE Plane
ggplot(ce_plane_data_deltas, aes(x = delta_daly, y = delta_cost)) +
  
  # Axes lines
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  
  # Main point
  geom_point(color = "goldenrod", shape = 8, size = 4) +
  geom_text(aes(label = strategy), vjust = -1.5, hjust = 0.5, size = 3.5) +
  
  # WTP line
  geom_abline(slope = 1263.47, intercept = 0, linetype = "dashed", color = "peachpuff4", size = 1) +
  
  # WTP label (adjusted coordinates)
  annotate("text", 
           x = 0.2, 
           y = 1263.47 * 0.15,  # approx 505
           label = "WTP Threshold", 
           color = "peachpuff4", 
           angle = 45,
           hjust = -0.1, 
           size = 3.2) +
  
  # Labels
  labs(
    x = "DALYs Averted",
    y = "Incremental Cost (USD)"
  ) +
  
  # Updated zoom window
  coord_cartesian(xlim = c(0, 1), ylim = c(-max_abs_y, 550)) +  # Make sure y-axis goes above 500
  
  # Theme settings
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, hjust = 0.5, face = "italic")
  )

##
ce_plane_data <- data.frame(
  strategy = c("Control", "Intervention"),
  cost = c(as.numeric(cost_control), as.numeric(cost_intervention)),
  ab_count = c(as.numeric(ABTotal_Control_Num), as.numeric(ABTotal_Intervention_Num))
)
ce_plane_data$ab_averted <- max(ce_plane_data$ab_count) - ce_plane_data$ab_count

ggplot(ce_plane_data, aes(x = ab_averted, y = cost, label = strategy)) +
  geom_point(size = 4) +
  geom_text(vjust = -1.2, hjust = 0.5) +
  labs(
    title = "Cost vs. Antibiotic Use Reduction",
    x = "Antibiotic Treated Patients Averted",
    y = "Cost (USD)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic")
  ) 

################################################################################
##SENSITIVITY ANALYSES
################################################################################
##ONE WAY SENSITIVITY ANALYSIS 
#Function to recreate costs and tree logic
build_owsa_tree <- function(cost_params, perspective = "societal") {
  with(as.list(cost_params), {
  ##TERMINAL NODE VALUES: CONTROL ARM
  ##Nodes: Malaria Positive, No AB
  cost1 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control)
  health1 <- DALY_moderate_success
  cost2 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control) 
  health2 <- DALY_moderate_failure
  ##Nodes: Malaria Positive, AB
  cost3 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_success + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * pHospitalRef_MalariaAB_Control)
  health3 <- DALY_moderate_success
  cost4 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control) 
  health4 <- DALY_moderate_failure
  
  ##Nodes: Malaria Negative, No AB
  cost5 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control)
  health5 <- DALY_moderate_success
  cost6 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control)
  health6 <- DALY_moderate_failure
  ##Nodes: Malaria Negative, AB
  cost7 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_success + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control)
  health7 <- DALY_moderate_success
  cost8 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control)
  health8 <- DALY_moderate_failure
  
  ##Nodes: No Malaria RDT, No AB
  cost9 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control) + cProductivity_success else cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control)
  health9 <- DALY_moderate_success
  cost10 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control) + cProductivity_failure else cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control)
  health10 <- DALY_moderate_failure
  ##Nodes: No Malaria RDT, AB
  cost11 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_success + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control)
  health11 <- DALY_moderate_success
  cost12 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control) 
  health12 <- DALY_moderate_failure
  
  
  ##TERMINAL NODE VALUES: INTERVENTION ARM, DANGER SIGNS PRESENT 
  ##Nodes: Severe Malaria, No AB (Nobody had severe malaria)
  cost13 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention)
  health13 <- DALY_moderate_success
  cost14 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention)
  health14 <- DALY_moderate_failure
  ##Nodes: Severe Malaria, AB (Nobody had severe malaria, so not adding AB costs)
  cost15 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention)
  health15 <- DALY_moderate_success
  cost16 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention)
  health16 <- DALY_moderate_failure
  
  ##Nodes: No Severe Malaria, No AB
  cost17 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention)
  health17 <- DALY_moderate_success
  cost18 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention)
  health18 <- DALY_moderate_failure
  ##Nodes: No Severe Malaria, AB
  cost19 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention)
  health19 <- DALY_moderate_success
  cost20 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention)
  health20 <- DALY_moderate_failure
  
  #Nodes: No Severe Malaria RDT, No AB
  cost21 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention)
  health21 <- DALY_moderate_success
  cost22 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention)
  health22 <- DALY_moderate_failure
  ##Nodes: No Severe Malaria RDT, AB
  cost23 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention)
  health23 <- DALY_moderate_success
  cost24 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention)
  health24 <- DALY_moderate_failure
  
  
  ##TERMINAL NODE VALUES: INTERVENTION ARM, NO DANGER SIGNS PRESENT
  ##Nodes: UR Symptoms Present, No AB
  cost25 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention)
  health25 <- DALY_moderate_success
  cost26 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention)
  health26 <- DALY_moderate_failure
  ##Nodes: UR Symptoms Present, AB
  cost27 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_success + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * pHospitalRef_URSAB_Intervention)
  health27 <- DALY_moderate_success
  cost28 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_failure + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention)
  health28 <- DALY_moderate_failure
  
  ##Nodes: UR Symptoms Not Present, No AB
  cost29 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention)
  health29 <- DALY_moderate_success
  cost30 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention)
  health30 <- DALY_moderate_failure
  ##Nodes: UR Symptoms Not Present, AB
  cost31 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_success + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention) 
  health31 <- DALY_moderate_success
  cost32 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_failure + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention)
  health32 <- DALY_moderate_failure
  
  tree <- Node$new("Trial Randomisation")
  
  ################################################################################
  ## CONTROL BRANCHES
  control <- AddNodeWithProb(tree, "Control", prob = 0.5)
  
  ## MALARIA RDT RESULT BRANCHES - CONTROL
  branchC1 <- AddNodeWithProb(control, "Malaria Positive", prob = pMalaria_Control)
  branchC2 <- AddNodeWithProb(control, "Malaria Negative", prob = pNoMalaria_Control)
  branchC3 <- AddNodeWithProb(control, "No Malaria RDT", prob = pNoMalariaRDT_Control)
  
  ## PRIMARY OUTCOME BRANCHES - CONTROL
  branchC1.1 <- AddNodeWithProb(branchC1, "No AB", prob = pMalariaNoAB_Control)
  branchC1.2 <- AddNodeWithProb(branchC1, "AB", prob = pMalariaAB_Control)
  branchC2.1 <- AddNodeWithProb(branchC2, "No AB", prob = pNoMalariaNoAB_Control)
  branchC2.2 <- AddNodeWithProb(branchC2, "AB", prob = pNoMalariaAB_Control)
  branchC3.1 <- AddNodeWithProb(branchC3, "No AB", prob = pNoMalariaRDTNoAB_Control)
  branchC3.2 <- AddNodeWithProb(branchC3, "AB", prob = pNoMalariaRDTAB_Control)
  
  ## TREATMENT SUCCESS BRANCHES & TERMINAL NODES - CONTROL
  branchC1.1.1 <- AddNodeWithProb(branchC1.1, "Treatment Success", prob = pRecovered_MalariaNoAB_Control)
  branchC1.1.1$cost <- cost1
  branchC1.1.1$daly <- health1
  
  branchC1.1.2 <- AddNodeWithProb(branchC1.1, "Treatment Failure", prob = 1 - pRecovered_MalariaNoAB_Control)
  branchC1.1.2$cost <- cost2
  branchC1.1.2$daly <- health2
  
  branchC1.2.1 <- AddNodeWithProb(branchC1.2, "Treatment Success", prob = pRecovered_MalariaAB_Control)
  branchC1.2.1$cost <- cost3
  branchC1.2.1$daly <- health3
  
  branchC1.2.2 <- AddNodeWithProb(branchC1.2, "Treatment Failure", prob = 1 - pRecovered_MalariaAB_Control)
  branchC1.2.2$cost <- cost4
  branchC1.2.2$daly <- health4
  
  branchC2.1.1 <- AddNodeWithProb(branchC2.1, "Treatment Success", prob = pRecovered_NoMalariaNoAB_Control)
  branchC2.1.1$cost <- cost5
  branchC2.1.1$daly <- health5
  
  branchC2.1.2 <- AddNodeWithProb(branchC2.1, "Treatment Failure", prob = 1 - pRecovered_NoMalariaNoAB_Control)
  branchC2.1.2$cost <- cost6
  branchC2.1.2$daly <- health6
  
  branchC2.2.1 <- AddNodeWithProb(branchC2.2, "Treatment Success", prob = pRecovered_NoMalariaAB_Control)
  branchC2.2.1$cost <- cost7
  branchC2.2.1$daly <- health7
  
  branchC2.2.2 <- AddNodeWithProb(branchC2.2, "Treatment Failure", prob = 1 - pRecovered_NoMalariaAB_Control)
  branchC2.2.2$cost <- cost8
  branchC2.2.2$daly <- health8
  
  branchC3.1.1 <- AddNodeWithProb(branchC3.1, "Treatment Success", prob = pRecovered_NoMalariaRDTNoAB_Control)
  branchC3.1.1$cost <- cost9
  branchC3.1.1$daly <- health9
  
  branchC3.1.2 <- AddNodeWithProb(branchC3.1, "Treatment Failure", prob = 1 - pRecovered_NoMalariaRDTNoAB_Control)
  branchC3.1.2$cost <- cost10
  branchC3.1.2$daly <- health10
  
  branchC3.2.1 <- AddNodeWithProb(branchC3.2, "Treatment Success", prob = pRecovered_NoMalariaRDTAB_Control)
  branchC3.2.1$cost <- cost11
  branchC3.2.1$daly <- health11
  
  branchC3.2.2 <- AddNodeWithProb(branchC3.2, "Treatment Failure", prob = 1 - pRecovered_NoMalariaRDTAB_Control)
  branchC3.2.2$cost <- cost12
  branchC3.2.2$daly <- health12
  
  ################################################################################
  ## INTERVENTION BRANCHES
  intervention <- AddNodeWithProb(tree, "Intervention", prob = 0.5)
  
  ## VITAL SIGNS ASSESSMENT RESULTS - INTERVENTION
  branchI1 <- AddNodeWithProb(intervention, "Danger Signs Present", prob = pDangerSigns_Intervention)
  branchI2 <- AddNodeWithProb(intervention, "No Danger Signs Present", prob = pNoDangerSigns_Intervention)
  
  ## DANGER SIGNS PRESENT - MALARIA RDT BRANCHES - INTERVENTION
  branchI1.1 <- AddNodeWithProb(branchI1, "Malaria Positive", prob = 0)
  branchI1.2 <- AddNodeWithProb(branchI1, "Malaria Negative", prob = pNoSevereMalaria_Intervention)
  branchI1.3 <- AddNodeWithProb(branchI1, "No Malaria RDT", prob = pNoSevereMalariaRDT_Intervention)
  
  ## DANGER SIGNS PRESENT - PRIMARY OUTCOME BRANCHES - INTERVENTION
  branchI1.1.1 <- AddNodeWithProb(branchI1.1, "No AB", prob = 0)
  branchI1.1.2 <- AddNodeWithProb(branchI1.1, "AB", prob = 0)
  branchI1.2.1 <- AddNodeWithProb(branchI1.2, "No AB", prob = pNoSevereMalariaNoAB_Intervention)
  branchI1.2.2 <- AddNodeWithProb(branchI1.2, "AB", prob = pNoSevereMalariaAB_Intervention)
  branchI1.3.1 <- AddNodeWithProb(branchI1.3, "No AB", prob = pNoSevereMalariaRDTNoAB_Intervention)
  branchI1.3.2 <- AddNodeWithProb(branchI1.3, "AB", prob = pNoSevereMalariaRDTAB_Intervention)
  
  ## DANGER SIGNS PRESENT - TREATMENT SUCCESS BRANCHES & TERMINAL NODES - INTERVENTION
  branchI1.1.1.1 <- AddNodeWithProb(branchI1.1.1, "Treatment Success", prob = 0)
  branchI1.1.1.1$cost <- cost13
  branchI1.1.1.1$daly <- health13
  
  branchI1.1.1.2 <- AddNodeWithProb(branchI1.1.1, "Treatment Failure", prob = 0)
  branchI1.1.1.2$cost <- cost14
  branchI1.1.1.2$daly <- health14
  
  branchI1.1.2.1 <- AddNodeWithProb(branchI1.1.2, "Treatment Success", prob = 0)
  branchI1.1.2.1$cost <- cost15
  branchI1.1.2.1$daly <- health15
  
  branchI1.1.2.2 <- AddNodeWithProb(branchI1.1.2, "Treatment Failure", prob = 0)
  branchI1.1.2.2$cost <- cost16
  branchI1.1.2.2$daly <- health16
  
  branchI1.2.1.1 <- AddNodeWithProb(branchI1.2.1, "Treatment Success", prob = pRecovered_NoSevereMalariaNoAB_Intervention)
  branchI1.2.1.1$cost <- cost17
  branchI1.2.1.1$daly <- health17
  
  branchI1.2.1.2 <- AddNodeWithProb(branchI1.2.1, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaNoAB_Intervention)
  branchI1.2.1.2$cost <- cost18
  branchI1.2.1.2$daly <- health18
  
  branchI1.2.2.1 <- AddNodeWithProb(branchI1.2.2, "Treatment Success", prob = pRecovered_NoSevereMalariaAB_Intervention)
  branchI1.2.2.1$cost <- cost19
  branchI1.2.2.1$daly <- health19
  
  branchI1.2.2.2 <- AddNodeWithProb(branchI1.2.2, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaAB_Intervention)
  branchI1.2.2.2$cost <- cost20
  branchI1.2.2.2$daly <- health20
  
  branchI1.3.1.1 <- AddNodeWithProb(branchI1.3.1, "Treatment Success", prob = pRecovered_NoSevereMalariaRDTNoAB_Intervention)
  branchI1.3.1.1$cost <- cost21
  branchI1.3.1.1$daly <- health21
  
  branchI1.3.1.2 <- AddNodeWithProb(branchI1.3.1, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaRDTNoAB_Intervention)
  branchI1.3.1.2$cost <- cost22
  branchI1.3.1.2$daly <- health22
  
  branchI1.3.2.1 <- AddNodeWithProb(branchI1.3.2, "Treatment Success", prob = pRecovered_NoSevereMalariaRDTAB_Intervention)
  branchI1.3.2.1$cost <- cost23
  branchI1.3.2.1$daly <- health23
  
  branchI1.3.2.2 <- AddNodeWithProb(branchI1.3.2, "Treatment Failure", prob = 1 - pRecovered_NoSevereMalariaRDTAB_Intervention)
  branchI1.3.2.2$cost <- cost24
  branchI1.3.2.2$daly <- health24
  
  ## NO DANGER SIGNS PRESENT - SYMPTOMS ASSESSMENT TO RDT - INTERVENTION 
  branchI2.1 <- AddNodeWithProb(branchI2, "Symptoms Assessment to Malaria RDT", prob = 1)
  
  ## NO DANGER SIGNS PRESENT - MALARIA RDT BRANCHES - INTERVENTION
  branchI2.2 <- AddNodeWithProb(branchI2.1, "UR Symptoms", prob = pURSymptoms_Intervention)
  branchI2.3 <- AddNodeWithProb(branchI2.1, "No UR Symptoms", prob = pNoURSymptoms_Intervention)
  
  ## NO DANGER SIGNS PRESENT - PRIMARY OUTCOME BRANCHES - INTERVENTION
  branchI2.2.1 <- AddNodeWithProb(branchI2.2, "No AB", prob = pURSNoAB_Intervention)
  branchI2.2.2 <- AddNodeWithProb(branchI2.2, "AB", prob = pURSAB_Intervention)
  branchI2.3.1 <- AddNodeWithProb(branchI2.3, "No AB", prob = pNoURSNoAB_Intervention)
  branchI2.3.2 <- AddNodeWithProb(branchI2.3, "AB", prob = pNoURSAB_Intervention)
  
  ## NO DANGER SIGNS PRESENT - TREATMENT SUCCESS BRANCHES & TERMINAL NODES - INTERVENTION
  branchI2.2.1.1 <- AddNodeWithProb(branchI2.2.1, "Treatment Success", prob = pRecovered_URSNoAB_Intervention)
  branchI2.2.1.1$cost <- cost25
  branchI2.2.1.1$daly <- health25
  
  branchI2.2.1.2 <- AddNodeWithProb(branchI2.2.1, "Treatment Failure", prob = 1 - pRecovered_URSNoAB_Intervention)
  branchI2.2.1.2$cost <- cost26
  branchI2.2.1.2$daly <- health26
  
  branchI2.2.2.1 <- AddNodeWithProb(branchI2.2.2, "Treatment Success", prob = pRecovered_URSAB_Intervention)
  branchI2.2.2.1$cost <- cost27
  branchI2.2.2.1$daly <- health27
  
  branchI2.2.2.2 <- AddNodeWithProb(branchI2.2.2, "Treatment Failure", prob = 1 - pRecovered_URSAB_Intervention)
  branchI2.2.2.2$cost <- cost28
  branchI2.2.2.2$daly <- health28
  
  branchI2.3.1.1 <- AddNodeWithProb(branchI2.3.1, "Treatment Success", prob = pRecovered_NoURSNoAB_Intervention)
  branchI2.3.1.1$cost <- cost29
  branchI2.3.1.1$daly <- health29
  
  branchI2.3.1.2 <- AddNodeWithProb(branchI2.3.1, "Treatment Failure", prob = 1 - pRecovered_NoURSNoAB_Intervention)
  branchI2.3.1.2$cost <- cost30
  branchI2.3.1.2$daly <- health30
  
  branchI2.3.2.1 <- AddNodeWithProb(branchI2.3.2, "Treatment Success", prob = pRecovered_NoURSAB_Intervention)
  branchI2.3.2.1$cost <- cost31
  branchI2.3.2.1$daly <- health31
  
  branchI2.3.2.2 <- AddNodeWithProb(branchI2.3.2, "Treatment Failure", prob = 1 - pRecovered_NoURSAB_Intervention)
  branchI2.3.2.2$cost <- cost32
  branchI2.3.2.2$daly <- health32
  
  return(tree)
  })
}

#Function to run the recreated tree with high and low estimates 
perform_owsa <- function(base_costs, change = 0.2) {
  owsa_results <- data.frame(
    Parameter = character(),
    Low_Value = numeric(),
    Low_ICER = numeric(),
    High_Value = numeric(),
    High_ICER = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (param in names(base_costs)) {
    cat("Running OWSA for:", param, "\n")
    
    base_val <- base_costs[[param]]
    
    # Create low and high cost lists by copying base_costs
    costs_low <- base_costs
    costs_high <- base_costs
    
    # Modify only the current parameter
    costs_low[[param]] <- base_val * (1 - change)
    costs_high[[param]] <- base_val * (1 + change)
    
    # Run model with low value
    tree_low <- build_owsa_tree(costs_low)
    CalculateExpectedValues(tree_low)
    ICER_low <- tryCatch({
      delta_cost_low <- tree_low$Intervention$expected_cost - tree_low$Control$expected_cost
      delta_daly_low <- tree_low$Control$expected_daly - tree_low$Intervention$expected_daly 
      delta_cost_low / delta_daly_low
    }, error = function(e) {
      warning(paste("LOW value failed for", param, ":", e$message))
      NA
    })
    
    # Run model with high value
    tree_high <- build_owsa_tree(costs_high)
    CalculateExpectedValues(tree_high)
    ICER_high <- tryCatch({
      delta_cost_high <- tree_high$Intervention$expected_cost - tree_high$Control$expected_cost
      delta_daly_high <- tree_high$Control$expected_daly - tree_high$Intervention$expected_daly
      delta_cost_high / delta_daly_high
    }, error = function(e) {
      warning(paste("HIGH value failed for", param, ":", e$message))
      NA
    })
    
    if (!is.na(ICER_low) && !is.na(ICER_high)) {
      owsa_results <- rbind(owsa_results, data.frame(
        Parameter = param,
        Low_Value = round(costs_low[[param]], 2),
        Low_ICER = round(ICER_low, 2),
        High_Value = round(costs_high[[param]], 2),
        High_ICER = round(ICER_high, 2)
      ))
    } else {
      warning(paste("Skipping", param, "due to NA ICERs"))
    }
  }
  
  return(owsa_results)
}

owsa_results <- perform_owsa(base_costs)

#Order by influence
owsa_results$ICER_range <- abs(owsa_results$percentage - owsa_results$percentage.1)
owsa_results <- owsa_results[order(-owsa_results$ICER_range), ]

##OWSA Visualisation
# Prepare the OWSA results for plotting
owsa_plot_data <- owsa_results
owsa_plot_data$Min_ICER <- pmin(owsa_results$percentage, owsa_results$percentage.1)
owsa_plot_data$Max_ICER <- pmax(owsa_results$percentage, owsa_results$percentage.1)
owsa_results$Parameter <- factor(owsa_results$Parameter,
                                 levels = owsa_results$Parameter[order(owsa_results$ICER_range)])


#Visual-Friendly names to include in Tornado Plot
costs_friendly_names <- c(
  cRDT = "Cost of Malaria RDT",
  cCRP = "Cost of CRP",
  cAppt_Control = "Appointment Cost (Control)",
  cAppt_Intervention = "Appointment Cost (Intervention)",
  cFollowUp = "Follow-Up Visit Cost",
  cRepresentPHC = "Re-Presentation Cost at PHC",
  cHospital_Admission = "Hospital Admission Cost",
  cEDAM = "EDAM Cost",
  cPulseox = "Pulse Oximeter Cost",
  cBloodpressure = "Blood Pressure Monitor Cost",
  
  ABWeightedCosts_WithMalaria_Control = "AB Cost (Malaria, Control)",
  ABWeightedCosts_WithoutMalaria_Control = "AB Cost (No Malaria, Control)",
  ABWeightedCosts_WithoutMalariaRDT_Control = "AB Cost (No Malaria, No RDT, Control)",
  ABWeightedCosts_WithoutSevereMalaria_Intervention = "AB Cost (No Severe Malaria, Intervention)",
  ABWeightedCosts_WithoutSevereMalariaRDT_Intervention = "AB Cost (No Severe Malaria, No RDT, Intervention)",
  ABWeightedCosts_URS_Intervention = "AB Cost (URS, Intervention)",
  ABWeightedCosts_NoURS_Intervention = "AB Cost (No URS, Intervention)",
  
  AMRWeightedCosts_WithMalaria_Control = "AMR Cost (Malaria, Control)",
  AMRWeightedCosts_WithoutMalaria_Control = "AMR Cost (No Malaria, Control)",
  AMRWeightedCosts_WithoutMalariaRDT_Control = "AMR Cost (No Malaria, No RDT, Control)",
  AMRWeightedCosts_WithoutSevereMalaria_Intervention = "AMR Cost (No Severe Malaria, Intervention)",
  AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention = "AMR Cost (No Severe Malaria, No RDT, Intervention)",
  AMRWeightedCosts_URS_Intervention = "AMR Cost (URS, Intervention)",
  AMRWeightedCosts_NoURS_Intervention = "AMR Cost (No URS, Intervention)"
)
owsa_results$Friendly_Name <- costs_friendly_names[owsa_results$Parameter]
owsa_plot_data$Friendly_Name <- costs_friendly_names[owsa_plot_data$Parameter]


#Plot
library(ggplot2)
library(dplyr)
library(scales)

#Keep only top 10 parameters by influence
owsa_top10 <- owsa_plot_data %>%
  mutate(
    ICER_Range = abs(Max_ICER - Min_ICER),
    is_top_param = FALSE  # initialize
  ) %>%
  arrange(desc(ICER_Range)) %>%
  slice(1:10) %>%
  mutate(
    is_top_param = row_number() == 1,
    Friendly_Name = factor(Friendly_Name, levels = rev(Friendly_Name))
  )


# Axis breaks and formatting
min_val <- min(owsa_top10$Min_ICER, na.rm = TRUE)
max_val <- max(owsa_top10$Max_ICER, na.rm = TRUE)
x_breaks <- pretty(c(min_val, max_val), n = 5)  # 5 breaks to reduce clutter

ggplot(owsa_top10, aes(y = Friendly_Name)) +
  
  # Main bars
  geom_segment(aes(x = Min_ICER, xend = Max_ICER, yend = Friendly_Name),
               linewidth = 6, color = "beige") +
  
  # Dots for low/high values
  geom_point(aes(x = percentage, color = "Low Value"), size = 3) +
  geom_point(aes(x = percentage.1, color = "High Value"), size = 3) +
  
  # Min and Max ICER labels for top parameter, inside the dots (on the beige bar)
  geom_text(
    data = owsa_top10 %>% filter(is_top_param),
    aes(x = Min_ICER, label = paste0("$", format(round(Min_ICER, 0), big.mark = ","))),
    hjust = -0.1, vjust = 0.5, size = 2.2, color = "black"
  ) +
  geom_text(
    data = owsa_top10 %>% filter(is_top_param),
    aes(x = Max_ICER, label = paste0("$", format(round(Max_ICER, 0), big.mark = ","))),
    hjust = 1.1, vjust = 0.5, size = 2.2, color = "black"
  ) +
  
  geom_text(
    data = owsa_top10 %>% filter(!is_top_param),
    aes(x = Min_ICER, label = paste0("$", format(round(Min_ICER, 0), big.mark = ","))),
    hjust = 1.2, vjust = 0.5, size = 2, color = "gray30"
  ) +
  geom_text(
    data = owsa_top10 %>% filter(!is_top_param),
    aes(x = Max_ICER, label = paste0("$", format(round(Max_ICER, 0), big.mark = ","))),
    hjust = -0.1, vjust = 0.5, size = 2, color = "gray30"
  ) +
  
  # Reference line
  geom_vline(
    xintercept = mean(c(owsa_results$percentage, owsa_results$percentage.1), na.rm = TRUE),
    linetype = "dashed", color = "gray50"
  ) +
  
  # Scales and formatting
  scale_x_continuous(
    breaks = x_breaks,
    labels = label_number(scale = 1e-6, suffix = " million", accuracy = 0.1)
  ) +
  scale_color_manual(values = c("Low Value" = "peachpuff4", "High Value" = "forestgreen")) +
  
  labs(
    x = "ICER",
    y = "Parameter",
    color = "Input Value"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

################################################################################
##PROBABILISTIC SENSITIVITY ANALYSIS
library(heemod)

N <- 1000 #run for 1000 iterations

set.seed(100)

##PSA PART 1: CREATE FUNCTIONS NECESSARY FOR PSA SAMPLE DISTRIBUTIONS

#One time function to apply CIs to Probability Distributions
get_Probbeta_params <- function(mean, lower, upper) {
  obj <- function(par) {
    a <- par[1]; b <- par[2]
    (a / (a + b) - mean)^2 +
      (qbeta(0.025, a, b) - lower)^2 +
      (qbeta(0.975, a, b) - upper)^2
  }
  fit <- optim(c(1, 1), obj, method = "L-BFGS-B", lower = c(0.001, 0.001))
  fit$par
}

get_beta_from_binomial <- function(prob, n) {
  x <- round(prob * n)
  bt <- binom.test(x, n)
  get_Probbeta_params(bt$estimate, bt$conf.int[1], bt$conf.int[2])
}

#One time function to apply CIs to Dirichlet Probability Distributions (Malaria RDT Outcome ONLY)
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  samples <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    samples[, i] <- rgamma(n, shape = alpha[i], rate = 1)
  }
  samples <- samples / rowSums(samples)
  return(samples)
}

get_dirichlet_params <- function(probs, n_total) {
  alpha <- probs * n_total
  return(alpha)
}
  
#One time function to apply IHME CIs to DALY Distributions
get_DALYbeta_params <- function(mean, lower, upper) {
  obj <- function(par) {
    a <- par[1]; b <- par[2]
    (a / (a + b) - mean)^2 +
      (qbeta(0.025, a, b) - lower)^2 +
      (qbeta(0.975, a, b) - upper)^2
  }
  fit <- optim(c(1, 1), obj, method = "L-BFGS-B", lower = c(0.001, 0.001))
  fit$par
}


##Dirichlet distribution labels needed for malaria RDT outcome
d1 <- rdirichlet(N, get_dirichlet_params(c(base_probs$pMalaria_Control, base_probs$pNoMalaria_Control, base_probs$pNoMalariaRDT_Control), 2428)) #vector of probabilities, sample size
d2 <- rdirichlet(N, get_dirichlet_params(c(base_probs$pSevereMalaria_Intervention, base_probs$pNoSevereMalaria_Intervention, base_probs$pNoSevereMalariaRDT_Intervention), 1302))


##PSA PART 2: CREATE PARAMETERS DISTRIBUTIONS 

##PARAMETER SAMPLES
psa_samples <- data.frame(
  cRDT = rgamma(N, shape = (base_costs$cRDT^2)/(0.2^2), rate = (base_costs$cRDT)/(0.2^2)),
  cCRP = rgamma(N, shape = (base_costs$cCRP^2)/(0.2^2), rate = (base_costs$cCRP)/(0.2^2)), #not final inputs, refer to paper instead
  cAppt_Control = rgamma(N, shape = (base_costs$cAppt_Control^2)/(0.2^2), rate = (base_costs$cAppt_Control)/(0.2^2)),
  cAppt_Intervention = rgamma(N, shape = (base_costs$cAppt_Intervention^2)/(0.2^2), rate = (base_costs$cAppt_Intervention)/(0.2^2)), #not final inputs
  cFollowUp = rgamma(N, shape = (base_costs$cFollowUp^2)/(0.2^2), rate = (base_costs$cFollowUp)/(0.2^2)),
  cRepresentPHC = rgamma(N, shape = (base_costs$cRepresentPHC^2)/(0.2^2), rate = (base_costs$cRepresentPHC)/(0.2^2)),
  cHospital_Admission = rgamma(N, shape = (base_costs$cHospital_Admission^2)/(0.2^2), rate = (base_costs$cHospital_Admission)/(0.2^2)),
  cEDAM = rgamma(N, shape = (base_costs$cEDAM^2)/(0.2^2), rate = (base_costs$cEDAM)/(0.2^2)),
  cPulseox = rgamma(N, shape =(base_costs$cPulseox^2)/(0.2^2), rate = (base_costs$cPulseox)/(0.2^2)),
  cBloodpressure = rgamma(N, shape =(base_costs$cBloodpressure^2)/(0.2^2), rate = (base_costs$cBloodpressure)/(0.2^2)),
  ABWeightedCosts_WithMalaria_Control = rgamma(N, shape =  (base_costs$ABWeightedCosts_WithMalaria_Control^2)/(0.2^2), rate =(base_costs$ABWeightedCosts_WithMalaria_Control)/(0.2^2)),
  ABWeightedCosts_WithoutMalaria_Control = rgamma(N, shape =(base_costs$ABWeightedCosts_WithoutMalaria_Control^2)/(0.2^2), rate = (base_costs$ABWeightedCosts_WithoutMalaria_Control)/(0.2^2)), 
  ABWeightedCosts_WithoutMalariaRDT_Control = rgamma(N, shape = (base_costs$ABWeightedCosts_WithoutMalariaRDT_Control^2)/(0.2^2), rate = (base_costs$ABWeightedCosts_WithoutMalariaRDT_Control)/(0.2^2)), 
  ABWeightedCosts_WithoutSevereMalaria_Intervention = rgamma(N, shape = (base_costs$ABWeightedCosts_WithoutSevereMalaria_Intervention^2)/(0.2^2), rate = (base_costs$ABWeightedCosts_WithoutSevereMalaria_Intervention)/(0.2^2)), 
  ABWeightedCosts_WithoutSevereMalariaRDT_Intervention = rgamma(N, shape = (base_costs$ABWeightedCosts_WithoutSevereMalariaRDT_Intervention^2)/(0.2^2), rate = (base_costs$ABWeightedCosts_WithoutSevereMalariaRDT_Intervention)/(0.2^2)), 
  ABWeightedCosts_URS_Intervention = rgamma(N, shape =  (base_costs$ABWeightedCosts_URS_Intervention^2)/(0.2^2), rate = (base_costs$ABWeightedCosts_URS_Intervention)/(0.2^2)), 
  ABWeightedCosts_NoURS_Intervention = rgamma(N, shape =(base_costs$ABWeightedCosts_NoURS_Intervention^2)/(0.2^2), rate = (base_costs$ABWeightedCosts_NoURS_Intervention)/(0.2^2)),
  AMRWeightedCosts_WithMalaria_Control = rgamma(N, shape =(base_costs$AMRWeightedCosts_WithMalaria_Control^2)/(0.2^2), rate = (base_costs$AMRWeightedCosts_WithMalaria_Control)/(0.2^2)),
  AMRWeightedCosts_WithoutMalaria_Control = rgamma(N, shape = (base_costs$AMRWeightedCosts_WithoutMalaria_Control^2)/(0.2^2), rate = (base_costs$AMRWeightedCosts_WithoutMalaria_Control)/(0.2^2)),
  AMRWeightedCosts_WithoutMalariaRDT_Control = rgamma(N,shape = (base_costs$AMRWeightedCosts_WithoutMalariaRDT_Control^2)/(0.2^2), rate = (base_costs$AMRWeightedCosts_WithoutMalariaRDT_Control)/(0.2^2)),
  AMRWeightedCosts_WithoutSevereMalaria_Intervention = rgamma(N, shape =(base_costs$AMRWeightedCosts_WithoutSevereMalaria_Intervention^2)/(0.2^2), rate = (base_costs$AMRWeightedCosts_WithoutSevereMalaria_Intervention)/(0.2^2)),
  AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention = rgamma(N, shape = (base_costs$AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention^2)/(0.2^2), rate = (base_costs$AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention)/(0.2^2)), 
  AMRWeightedCosts_URS_Intervention = rgamma(N, shape = (base_costs$AMRWeightedCosts_URS_Intervention^2)/(0.2^2), rate = (base_costs$AMRWeightedCosts_URS_Intervention)/(0.2^2)), 
  AMRWeightedCosts_NoURS_Intervention = rgamma(N, shape = (base_costs$AMRWeightedCosts_NoURS_Intervention^2)/(0.2^2), rate = (base_costs$AMRWeightedCosts_NoURS_Intervention)/(0.2^2)),
  #probabilities samples
  #Dirichlet Distribution for 3 branches (ONLY Malaria RDT Outcomes)
  pMalaria_Control = d1[,1], 
  pNoMalaria_Control = d1[,2], 
  pNoMalariaRDT_Control = d1[,3],
 
  pMalariaAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pMalariaAB_Control, n = 9)
    rbeta(N, p[1], p[2])
  },
  pMalariaNoAB_Control =  {
    p <- get_beta_from_binomial(prob = base_probs$pMalariaNoAB_Control, n = 9)
    rbeta(N, p[1], p[2])
  },
  pNoMalariaAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pNoMalariaAB_Control, n = 1487)
    rbeta(N, p[1], p[2])
  },
  pNoMalariaNoAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pNoMalariaNoAB_Control, n = 1487)
    rbeta(N, p[1], p[2])
  },
  pNoMalariaRDTAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pNoMalariaRDTAB_Control, n = 926)
    rbeta(N, p[1], p[2])
  },
  pNoMalariaRDTNoAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pNoMalariaRDTNoAB_Control, n = 926)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_MalariaAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_MalariaAB_Control, n = 8)
    rbeta(N, p[1], p[2])
  },

  pRecovered_MalariaNoAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_MalariaNoAB_Control, n = 1)
    rbeta(N, p[1], p[2])
  },
  pRecovered_NoMalariaAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoMalariaAB_Control, n = 956)
    rbeta(N, p[1], p[2])
  },
  pRecovered_NoMalariaNoAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoMalariaNoAB_Control, n = 531)
    rbeta(N, p[1], p[2])
  },
  pRecovered_NoMalariaRDTAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoMalariaRDTAB_Control, n = 502)
    rbeta(N, p[1], p[2])
  },
  pRecovered_NoMalariaRDTNoAB_Control = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoMalariaRDTNoAB_Control, n = 424)
    rbeta(N, p[1], p[2])
  },
  
  pDangerSigns_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pDangerSigns_Intervention, n = 2324)
    rbeta(N, p[1], p[2])
  },
  #Dirichlet Distribution for 3 branches (ONLY Malaria RDT Outcomes)
  pSevereMalaria_Intervention = d2[,1], 
  pNoSevereMalaria_Intervention = d2[,2], 
  pNoSevereMalariaRDT_Intervention = d2[,3],

  pNoSevereMalariaAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pNoSevereMalariaAB_Intervention, n = 556)
    rbeta(N, p[1], p[2])
  },
  pNoSevereMalariaNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pNoSevereMalariaNoAB_Intervention, n = 374)
    rbeta(N, p[1], p[2])
  },
  
  pNoSevereMalariaRDTAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pNoSevereMalariaRDTAB_Intervention, n = 194)
    rbeta(N, p[1], p[2])
  },
  
  pNoSevereMalariaRDTNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pNoSevereMalariaRDTNoAB_Intervention, n = 178)
    rbeta(N, p[1], p[2])
  },
  
  pNoURSymptoms_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pNoURSymptoms_Intervention, n = 1022)
    rbeta(N, p[1], p[2])
  },
  
  pNoURSAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pNoURSAB_Intervention, n = 162)
    rbeta(N, p[1], p[2])
  },
  
  pNoURSNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pNoURSNoAB_Intervention, n = 162)
    rbeta(N, p[1], p[2])
  },
  
  pURSymptoms_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pURSymptoms_Intervention, n = 1022)
    rbeta(N, p[1], p[2])
  },
  
  pURSAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pURSAB_Intervention, n = 860)
    rbeta(N, p[1], p[2])
  },
  
  pURSNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pURSNoAB_Intervention, n = 860)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_NoSevereMalariaAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoSevereMalariaAB_Intervention, n = 556)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_NoSevereMalariaNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoSevereMalariaNoAB_Intervention, n = 374)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_NoSevereMalariaRDTAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoSevereMalariaRDTAB_Intervention, n = 194)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_NoSevereMalariaRDTNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoSevereMalariaRDTNoAB_Intervention, n = 178)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_NoURSAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoURSAB_Intervention, n = 34)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_NoURSNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_NoURSNoAB_Intervention, n = 128)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_URSAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_URSAB_Intervention, n = 520)
    rbeta(N, p[1], p[2])
  },
  
  pRecovered_URSNoAB_Intervention = {
    p <- get_beta_from_binomial(prob = base_probs$pRecovered_URSNoAB_Intervention, n = 340)
    rbeta(N, p[1], p[2])
  },
  #Add DALYs samples
  moderate_disabilityweight = { p <- get_DALYbeta_params(0.051, 0.032, 0.074); rbeta(N, p[1], p[2]) },
  severe_disabilityweight = { p <- get_DALYbeta_params(0.051, 0.032, 0.074); rbeta(N, p[1], p[2]) }
)

##Empty dataframe to store PSA results in
psa_results <- data.frame(psa_intervention_cost = numeric(0), psa_intervention_daly = numeric(0), 
                          psa_control_cost = numeric(0), psa_control_daly = numeric(0),
                          psa_delta_cost = numeric(0), psa_delta_daly = numeric(0),
                          ab_use_intervention = numeric(0), psa_ab_averted = numeric(0), 
                          psa_ICER = numeric(0),psa_ICER_ab = numeric(0), psa_intervention_nmb = numeric(0), 
                          psa_control_nmb = numeric(0))
psa_nmb <- data.frame(psa_intervention_nmb = numeric(0), psa_control_nmb = numeric(0))

##PSA PART 3: BUILD PSA LOOP
psa_loop <- for(i in 1:N){

  ##PART 3.1: EXTRACT PARAMETERS FOR RUN
  cRDT <- psa_samples$cRDT[i]
  cCRP <- psa_samples$cCRP[i]
  cAppt_Control <- psa_samples$cAppt_Control[i]
  cAppt_Intervention <- psa_samples$cAppt_Intervention[i]
  cFollowUp <- psa_samples$cFollowUp[i]
  cRepresentPHC <- psa_samples$cRepresentPHC[i]
  cHospital_Admission <- psa_samples$cHospital_Admission[i]
  cEDAM <- psa_samples$cEDAM[i]
  cPulseox <- psa_samples$cPulseox[i]
  cBloodpressure <- psa_samples$cBloodpressure[i]
  ABWeightedCosts_WithMalaria_Control <- psa_samples$ABWeightedCosts_WithMalaria_Control[i]
  ABWeightedCosts_WithoutMalaria_Control <- psa_samples$ABWeightedCosts_WithoutMalaria_Control[i]
  ABWeightedCosts_WithoutMalariaRDT_Control <- psa_samples$ABWeightedCosts_WithoutMalariaRDT_Control[i]
  ABWeightedCosts_WithoutSevereMalaria_Intervention <- psa_samples$ABWeightedCosts_WithoutSevereMalaria_Intervention[i]
  ABWeightedCosts_WithoutSevereMalariaRDT_Intervention <- psa_samples$ABWeightedCosts_WithoutSevereMalariaRDT_Intervention[i]
  ABWeightedCosts_URS_Intervention <- psa_samples$ABWeightedCosts_URS_Intervention[i]
  ABWeightedCosts_NoURS_Intervention <- psa_samples$ABWeightedCosts_NoURS_Intervention[i]
  AMRWeightedCosts_WithMalaria_Control <- psa_samples$AMRWeightedCosts_WithMalaria_Control[i]
  AMRWeightedCosts_WithoutMalaria_Control <- psa_samples$AMRWeightedCosts_WithoutMalaria_Control[i]
  AMRWeightedCosts_WithoutMalariaRDT_Control <- psa_samples$AMRWeightedCosts_WithoutMalariaRDT_Control[i]
  AMRWeightedCosts_WithoutSevereMalaria_Intervention <- psa_samples$AMRWeightedCosts_WithoutSevereMalaria_Intervention[i]
  AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention <- psa_samples$AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention[i]
  AMRWeightedCosts_URS_Intervention <- psa_samples$AMRWeightedCosts_URS_Intervention[i]
  AMRWeightedCosts_NoURS_Intervention <- psa_samples$AMRWeightedCosts_NoURS_Intervention[i]
  #Add Probabilities Samples
  pMalaria_Control <- psa_samples$pMalaria_Control[i]
  pNoMalaria_Control <- psa_samples$pNoMalaria_Control[i]
  pNoMalariaRDT_Control <- psa_samples$pNoMalariaRDT_Control[i]
  pMalariaAB_Control <- psa_samples$pMalariaAB_Control[i]
  pMalariaNoAB_Control <- psa_samples$pMalariaNoAB_Control[i]
  pNoMalariaAB_Control <- psa_samples$pNoMalariaAB_Control[i]
  pNoMalariaNoAB_Control <- psa_samples$pNoMalariaNoAB_Control[i]
  pNoMalariaRDTAB_Control <- psa_samples$pNoMalariaRDTAB_Control[i]
  pNoMalariaRDTNoAB_Control <- psa_samples$pNoMalariaRDTNoAB_Control[i]
  pRecovered_MalariaAB_Control <- psa_samples$pRecovered_MalariaAB_Control[i]
  pRecovered_MalariaNoAB_Control <- psa_samples$pRecovered_MalariaNoAB_Control[i]
  pRecovered_NoMalariaAB_Control <- psa_samples$pRecovered_NoMalariaAB_Control[i]
  pRecovered_NoMalariaNoAB_Control <- psa_samples$pRecovered_NoMalariaNoAB_Control[i]
  pRecovered_NoMalariaRDTAB_Control <- psa_samples$pRecovered_NoMalariaRDTAB_Control[i]
  pRecovered_NoMalariaRDTNoAB_Control <- psa_samples$pRecovered_NoMalariaRDTNoAB_Control[i]
  pDangerSigns_Intervention <- psa_samples$pDangerSigns_Intervention[i]
  pSevereMalaria_Intervention <- psa_samples$pSevereMalaria_Intervention[i]
  pNoSevereMalaria_Intervention <- psa_samples$pNoSevereMalaria_Intervention[i]
  pNoSevereMalariaRDT_Intervention <- psa_samples$pNoSevereMalariaRDT_Intervention[i]
  pNoSevereMalariaAB_Intervention <- psa_samples$pNoSevereMalariaAB_Intervention[i]
  pNoSevereMalariaNoAB_Intervention <- psa_samples$pNoSevereMalariaNoAB_Intervention[i]
  pNoSevereMalariaRDTAB_Intervention <- psa_samples$pNoSevereMalariaRDTAB_Intervention[i]
  pNoSevereMalariaRDTNoAB_Intervention <- psa_samples$pNoSevereMalariaRDTNoAB_Intervention[i]
  pNoURSymptoms_Intervention <- psa_samples$pNoURSymptoms_Intervention[i]
  pNoURSAB_Intervention <- psa_samples$pNoURSAB_Intervention[i]
  pNoURSNoAB_Intervention <- psa_samples$pNoURSNoAB_Intervention[i]
  pURSymptoms_Intervention <- psa_samples$pURSymptoms_Intervention[i]
  pURSAB_Intervention <- psa_samples$pURSAB_Intervention[i]
  pURSNoAB_Intervention <- psa_samples$pURSNoAB_Intervention[i]
  pRecovered_NoSevereMalariaAB_Intervention <- psa_samples$pRecovered_NoSevereMalariaAB_Intervention[i]
  pRecovered_NoSevereMalariaNoAB_Intervention <- psa_samples$pRecovered_NoSevereMalariaNoAB_Intervention[i]
  pRecovered_NoSevereMalariaRDTAB_Intervention <- psa_samples$pRecovered_NoSevereMalariaRDTAB_Intervention[i]
  pRecovered_NoSevereMalariaRDTNoAB_Intervention <- psa_samples$pRecovered_NoSevereMalariaRDTNoAB_Intervention[i]
  pRecovered_URSAB_Intervention <- psa_samples$pRecovered_URSAB_Intervention[i]
  pRecovered_URSNoAB_Intervention <- psa_samples$pRecovered_URSNoAB_Intervention[i]
  pRecovered_NoURSAB_Intervention <- psa_samples$pRecovered_NoURSAB_Intervention[i]
  pRecovered_NoURSNoAB_Intervention <- psa_samples$pRecovered_NoURSNoAB_Intervention[i]
  #Add DALYs samples
  moderate_disabilityweight <- psa_samples$moderate_disabilityweight[i]
  severe_disabilityweight <- psa_samples$severe_disabilityweight[i]
  
  # CHECK 1: Are sampled parameters different?
  if(i <= 3) {
    cat("Iteration", i, ": cRDT =", cRDT, ", cCRP =", cCRP, "\n")
  }
  
  
  ##PART 3.2: COPY PARMS INTO THIS PSA FUNCTION
  ##COPY DALYS INTO FUNCTION
  DALY_moderate_success <- moderate_disabilityweight * (3.5/365)
  DALY_moderate_failure <- moderate_disabilityweight * (7/365)
  
  ##COPY COSTS & HEALTH OUTCOMES INTO FUNCTION
  ##Nodes: Malaria Positive, No AB
  cost1 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control)
  health1 <- DALY_moderate_success
  cost2 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control) 
  health2 <- DALY_moderate_failure
  ##Nodes: Malaria Positive, AB
  cost3 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_success + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * pHospitalRef_MalariaAB_Control)
  health3 <- DALY_moderate_success
  cost4 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control) 
  health4 <- DALY_moderate_failure
  
  ##Nodes: Malaria Negative, No AB
  cost5 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control)
  health5 <- DALY_moderate_success
  cost6 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control)
  health6 <- DALY_moderate_failure
  ##Nodes: Malaria Negative, AB
  cost7 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_success + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control)
  health7 <- DALY_moderate_success
  cost8 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control)
  health8 <- DALY_moderate_failure
  
  ##Nodes: No Malaria RDT, No AB
  cost9 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control) + cProductivity_success else cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control)
  health9 <- DALY_moderate_success
  cost10 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control) + cProductivity_failure else cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control)
  health10 <- DALY_moderate_failure
  ##Nodes: No Malaria RDT, AB
  cost11 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_success + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control)
  health11 <- DALY_moderate_success
  cost12 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control) 
  health12 <- DALY_moderate_failure
  
  
  ##TERMINAL NODE VALUES: INTERVENTION ARM, DANGER SIGNS PRESENT 
  ##Nodes: Severe Malaria, No AB (Nobody had severe malaria)
  cost13 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention)
  health13 <- DALY_moderate_success
  cost14 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention)
  health14 <- DALY_moderate_failure
  ##Nodes: Severe Malaria, AB (Nobody had severe malaria, so not adding AB costs)
  cost15 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention)
  health15 <- DALY_moderate_success
  cost16 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention)
  health16 <- DALY_moderate_failure
  
  ##Nodes: No Severe Malaria, No AB
  cost17 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention)
  health17 <- DALY_moderate_success
  cost18 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention)
  health18 <- DALY_moderate_failure
  ##Nodes: No Severe Malaria, AB
  cost19 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention)
  health19 <- DALY_moderate_success
  cost20 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention)
  health20 <- DALY_moderate_failure
  
  #Nodes: No Severe Malaria RDT, No AB
  cost21 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention)
  health21 <- DALY_moderate_success
  cost22 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention)
  health22 <- DALY_moderate_failure
  ##Nodes: No Severe Malaria RDT, AB
  cost23 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention)
  health23 <- DALY_moderate_success
  cost24 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention)
  health24 <- DALY_moderate_failure
  
  
  ##TERMINAL NODE VALUES: INTERVENTION ARM, NO DANGER SIGNS PRESENT
  ##Nodes: UR Symptoms Present, No AB
  cost25 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention)
  health25 <- DALY_moderate_success
  cost26 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention)
  health26 <- DALY_moderate_failure
  ##Nodes: UR Symptoms Present, AB
  cost27 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_success + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * pHospitalRef_URSAB_Intervention)
  health27 <- DALY_moderate_success
  cost28 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_failure + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention)
  health28 <- DALY_moderate_failure
  
  ##Nodes: UR Symptoms Not Present, No AB
  cost29 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention)
  health29 <- DALY_moderate_success
  cost30 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention)
  health30 <- DALY_moderate_failure
  ##Nodes: UR Symptoms Not Present, AB
  cost31 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_success + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention) 
  health31 <- DALY_moderate_success
  cost32 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_failure + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention)
  health32 <- DALY_moderate_failure
  
  
  
  # CHECK 2: Are calculated costs different?
  if(i <= 3) {
    cat("Iteration", i, ": cost1 =", cost1, ", cost3 =", cost3, "\n")
  }
  
  ##CLONE TREE
  cloned_tree <- tree$clone()
  
  cloned_tree$Do(function(node) {
    if (!node$isRoot && !is.null(node$prob)) {
      node$p <- node$prob
    }
  })
  
  ##PART 3.3: REASSIGN COST AND QALY TO TERMINAL NODES (BECAUSE THEY ARE STATIC)
  node1 <- cloned_tree$Control$`Malaria Positive`$`No AB`$`Treatment Success`
  node1$cost <- cost1
  node1$daly <- health1
  
  node2 <- cloned_tree$Control$`Malaria Positive`$`No AB`$`Treatment Failure`
  node2$cost <- cost2 
  node2$daly <- health2
  
  node3 <- cloned_tree$Control$`Malaria Positive`$AB$`Treatment Success`
  node3$cost <- cost3
  node3$daly <- health3
  node3$ab_count <- 8 #hard coded from the AB Type Probabilities data frames
  node3$ab_path <- 0.5 * psa_samples$pMalaria_Control[i] * psa_samples$pMalariaAB_Control[i] *
    psa_samples$pRecovered_MalariaAB_Control[i]
  
  node4 <- cloned_tree$Control$`Malaria Positive`$AB$`Treatment Failure`
  node4$cost <- cost4 
  node4$daly <- health4
  node4$ab_count <- 8
  node4$ab_path <- 0.5 * psa_samples$pMalaria_Control[i] * psa_samples$pMalariaAB_Control[i] *
    (1 - psa_samples$pRecovered_MalariaAB_Control[i])
  
  node5 <- cloned_tree$Control$`Malaria Negative`$`No AB`$`Treatment Success`
  node5$cost <- cost5 
  node5$daly <- health5
  
  node6 <- cloned_tree$Control$`Malaria Negative`$`No AB`$`Treatment Failure`
  node6$cost <- cost6
  node6$daly <- health6
  
  node7 <- cloned_tree$Control$`Malaria Negative`$AB$`Treatment Success`
  node7$cost <- cost7
  node7$daly <- health7
  node7$ab_count <- 956
  node7$ab_path <- 0.5 * psa_samples$pNoMalaria_Control[i] * psa_samples$pNoMalariaAB_Control[i] *
    psa_samples$pRecovered_NoMalariaAB_Control[i]
  
  node8 <- cloned_tree$Control$`Malaria Negative`$AB$`Treatment Failure`
  node8$cost <- cost8
  node8$daly <- health8
  node8$ab_count <- 956
  node8$ab_path <- 0.5 * psa_samples$pNoMalaria_Control[i] * psa_samples$pNoMalariaAB_Control[i] *
    (1 - psa_samples$pRecovered_NoMalariaAB_Control[i])
  
  node9 <- cloned_tree$Control$`No Malaria RDT`$`No AB`$`Treatment Success`
  node9$cost <- cost9
  node9$daly <- health9
  
  node10 <- cloned_tree$Control$`No Malaria RDT`$`No AB`$`Treatment Failure`
  node10$cost <- cost10
  node10$daly <- health10
  
  node11 <- cloned_tree$Control$`No Malaria RDT`$AB$`Treatment Success`
  node11$cost <- cost11
  node11$daly <- health11
  node11$ab_count <- 502
  node11$ab_path <- 0.5 * psa_samples$pNoMalariaRDT_Control[i] * psa_samples$pNoMalariaRDTAB_Control[i] *
    psa_samples$pRecovered_NoMalariaRDTAB_Control[i]
  
  node12 <- cloned_tree$Control$`No Malaria RDT`$AB$`Treatment Failure`
  node12$cost <- cost12
  node12$daly <- health12
  node12$ab_count <- 502
  node12$ab_path <- 0.5 * psa_samples$pNoMalariaRDT_Control[i] * psa_samples$pNoMalariaRDTAB_Control[i] *
    (1 -psa_samples$pRecovered_NoMalariaRDTAB_Control[i])
  
  node13 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Positive`$`No AB`$`Treatment Success`
  node13$cost <- cost13
  node13$daly <- health13
  
  node14 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Positive`$`No AB`$`Treatment Failure`
  node14$cost <- cost14
  node14$daly <- health14
  
  node15 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Positive`$AB$`Treatment Success`
  node15$cost <- cost15
  node15$daly <- health15
  #Not adding AB Count because no Severe Malaria in Intervention Arm
  
  node16 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Positive`$AB$`Treatment Failure`
  node16$cost <- cost16
  node16$daly <- health16
  
  node17 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Negative`$`No AB`$`Treatment Success`
  node17$cost <- cost17
  node17$daly <- health17
  
  node18 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Negative`$`No AB`$`Treatment Failure`
  node18$cost <- cost18 
  node18$daly <- health18
  
  node19 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$`Treatment Success`
  node19$cost <- cost19
  node19$daly <- health19
  node19$ab_count <- 556
  node19$ab_path <- 0.5 * psa_samples$pDangerSigns_Intervention[i] * psa_samples$pNoSevereMalaria_Intervention[i] *
    psa_samples$pNoSevereMalariaAB_Intervention[i] * psa_samples$pRecovered_NoSevereMalariaAB_Intervention[i]
  
  node20 <- cloned_tree$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$`Treatment Failure`
  node20$cost <- cost20
  node20$daly <- health20
  node20$ab_count <- 556
  node20$ab_path <- 0.5 * psa_samples$pDangerSigns_Intervention[i] * psa_samples$pNoSevereMalaria_Intervention[i] *
    psa_samples$pNoSevereMalariaAB_Intervention[i] * (1 - psa_samples$pRecovered_NoSevereMalariaAB_Intervention[i])
  
  node21 <- cloned_tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$`No AB`$`Treatment Success`
  node21$cost <- cost21
  node21$daly <- health21
  
  node22 <- cloned_tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$`No AB`$`Treatment Failure`
  node22$cost <- cost22
  node22$daly <- health22
  
  node23 <- cloned_tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$AB$`Treatment Success`
  node23$cost <- cost23
  node23$daly <- health23
  node23$ab_count <- 194
  node23$ab_path <- 0.5 * psa_samples$pDangerSigns_Intervention[i] * psa_samples$pNoSevereMalariaRDT_Intervention[i] * 
    psa_samples$pNoSevereMalariaRDTAB_Intervention[i] * psa_samples$pRecovered_NoSevereMalariaRDTAB_Intervention[i]
  
  node24 <- cloned_tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$AB$`Treatment Failure`
  node24$cost <- cost24
  node24$daly <- health24
  node24$ab_count <- 194
  node24$ab_path <- 0.5 * psa_samples$pDangerSigns_Intervention[i] * psa_samples$pNoSevereMalariaRDT_Intervention[i] * 
    psa_samples$pNoSevereMalariaRDTAB_Intervention[i] * psa_samples$pRecovered_NoSevereMalariaRDTAB_Intervention[i]
  
  node25 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$`No AB`$`Treatment Success`
  node25$cost <- cost25
  node25$daly <- health25
  
  node26 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$`No AB`$`Treatment Failure`
  node26$cost <- cost26
  node26$daly <- health26
  
  node27 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$AB$`Treatment Success`
  node27$cost <- cost27
  node27$daly <- health27
  node27$ab_count <- 520 
  node27$ab_path <- 0.5 * (1 - psa_samples$pDangerSigns_Intervention[i]) * psa_samples$pURSymptoms_Intervention[i] *
    psa_samples$pNoURSAB_Intervention[i] * psa_samples$pRecovered_URSAB_Intervention[i]
  
  node28 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$AB$`Treatment Failure`
  node28$cost <- cost28
  node28$daly <- health28
  node28$ab_count <- 520
  node28$ab_path <- 0.5 * (1 - psa_samples$pDangerSigns_Intervention[i]) * psa_samples$pURSymptoms_Intervention[i] *
    psa_samples$pNoURSAB_Intervention[i] * (1 - psa_samples$pRecovered_URSAB_Intervention[i])
  
  node29 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$`No AB`$`Treatment Success`
  node29$cost <- cost29
  node29$daly <- health29
  
  node30 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$`No AB`$`Treatment Failure`
  node30$cost <- cost30
  node30$daly <- health30
  
  node31 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$AB$`Treatment Success`
  node31$cost <- cost31
  node31$daly <- health31
  node31$ab_count <- 34 
  node31$ab_path <- 0.5 * (1 - psa_samples$pDangerSigns_Intervention[i]) * psa_samples$pNoURSymptoms_Intervention[i] * 
    psa_samples$pNoURSAB_Intervention[i] * psa_samples$pRecovered_NoURSAB_Intervention[i]
  
  node32 <- cloned_tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$AB$`Treatment Failure`
  node32$cost <- cost32
  node32$daly <- health32
  node32$ab_count <- 34 
  node32$ab_path <- 0.5 * (1 - psa_samples$pDangerSigns_Intervention[i]) * psa_samples$pNoURSymptoms_Intervention[i] * 
    psa_samples$pNoURSAB_Intervention[i] * (1 - psa_samples$pRecovered_NoURSAB_Intervention[i])
  

  # CHECK 3: Are nodes getting updated with different values?
  if(i <= 3) {
    cat("Iteration", i, ": node1$cost =", node1$cost, "\n")
  }
  
  terminal_nodes <- cloned_tree$leaves
  control_nodes <- Filter(function(node) "Control" %in% node$path, terminal_nodes)
  intervention_nodes <- Filter(function(node) "Intervention" %in% node$path, terminal_nodes)
  
  ##PART 3.4: RE-EVALUATE AND STORE
  intervention_branch <- cloned_tree$Get(function(x) x, filterFun = function(x) x$name == "Intervention")[[1]]
  control_branch <- cloned_tree$Get(function(x) x, filterFun = function(x) x$name == "Control")[[1]]
  
  CalculateExpectedValues(intervention_branch)
  CalculateExpectedValues(control_branch)
  
  psa_intervention_cost <- intervention_branch$expected_cost
  psa_intervention_daly <- intervention_branch$expected_daly
  psa_control_cost <- control_branch$expected_cost
  psa_control_daly <- control_branch$expected_daly
  
  ab_use_control <- sum(vapply(control_nodes, function(node) {
    ab <- if (!is.null(node$ab_count) && !is.na(node$ab_count)) node$ab_count else 0
    pr <- if (!is.null(node$ab_path) && !is.na(node$ab_path)) node$ab_path else 0
    ab * pr
  }, numeric(1)))
  
  
  ab_use_intervention <- sum(vapply(intervention_nodes, function(node) {
    ab <- if (!is.null(node$ab_count) && !is.na(node$ab_count)) node$ab_count else 0
    pr <- if (!is.null(node$ab_path) && !is.na(node$ab_path)) node$ab_path else 0
    ab * pr
  }, numeric(1)))
  
  ab_averted <- ab_use_control - ab_use_intervention
  
  ##CALCULATE ICER
  psa_delta_cost <- psa_intervention_cost - psa_control_cost
  psa_delta_daly <- psa_control_daly - psa_intervention_daly
  
  psa_ICER <- psa_delta_cost / psa_delta_daly
  psa_ICER_ab <- psa_delta_cost / ab_averted
  
  ##Net Monetary Benefit for EVPI 
  psa_intervention_nmb <- (Cambodia_WTP * psa_intervention_daly) - psa_intervention_cost
  psa_control_nmb <- (Cambodia_WTP * psa_control_daly) - psa_control_cost
  
   psa_results <- rbind(
    psa_results,
    data.frame(
      psa_intervention_cost = as.numeric(psa_intervention_cost),
      psa_intervention_daly = as.numeric(psa_intervention_daly),
      psa_control_cost = as.numeric(psa_control_cost),
      psa_control_daly = as.numeric(psa_control_daly),
      psa_delta_cost = as.numeric(psa_delta_cost),
      psa_delta_daly = as.numeric(psa_delta_daly),
      ab_use_intervention = as.numeric(ab_use_intervention),
      psa_ab_averted = as.numeric(ab_averted),
      psa_ICER = as.numeric(psa_ICER),
      psa_ICER_ab = as.numeric(psa_ICER_ab),
      psa_intervention_nmb = as.numeric(psa_intervention_nmb),
      psa_control_nmb = as.numeric(psa_control_nmb)
    )
  )
}


################################################################################
##PSA Visualisations
##Histograms and CE Plane for Supplemental Information S1
hist(psa_results$psa_intervention_cost, 
     main = "", 
     xlab = "EDAM Intervention Cost (USD)", 
     col = "peachpuff4", 
     breaks = 30)
title(main = "PSA: Distribution of Intervention Costs", 
      font.main = 3,      # bold
      line = 1,           # distance from top
      cex.main = 1)     # size of the title text


hist(psa_results$psa_ab_averted, 
     main = "", 
     xlab = "Antibiotic Treated Patients Averted", 
     col = "peachpuff4", 
     breaks = 30)
title(main = "PSA: Distribution of Antibiotic Treated Patients Averted with EDAM", 
      font.main = 3,      # bold
      line = 1,           # distance from top
      cex.main = 1)     # size of the title text

hist(psa_results$psa_intervention_nmb, 
     main = "", 
     xlab = " EDAM Net Monetary Benefit", 
     col = "peachpuff4", 
     breaks = 30)
title(main = "PSA: Distribution of EDAM Net Monetary Benefit", 
      font.main = 3,      # bold
      line = 1,           # distance from top
      cex.main = 1)     # size of the title text

##CE Planes
##DALYs
plot(psa_results$psa_delta_daly, psa_results$psa_delta_cost,
     xlab = "Incremental DALYs Averted",
     ylab = "Incremental Cost (USD)",
     main = "PSA: EDAM Cost-Effectiveness Plane (ICER by DALYs Averted)",
     font.main = 3,  # Italic
     pch = 19, col = rgb(0.4, 0.3, 0.1, 0.5) )  # semi-transparent green

abline(h = 0, col = "peachpuff4")  # zero cost line
abline(v = 0, col = "peachpuff4")  # zero effect line

##CE Plane for Main Report Chapter 3
##ABs Averted
# Calculate mean values for highlighting
mean_ab <- mean(psa_results$psa_ab_averted, na.rm = TRUE)
mean_cost <- mean(psa_results$psa_delta_cost, na.rm = TRUE)

# Expand axis ranges slightly for breathing room
x_range <- range(psa_results$psa_ab_averted, na.rm = TRUE)
y_range <- range(psa_results$psa_delta_cost, na.rm = TRUE)

x_pad <- 0.1 * diff(x_range)
y_pad <- 0.1 * diff(y_range)

xlim <- c(x_range[1] - x_pad, x_range[2] + x_pad)
ylim <- c(y_range[1] - y_pad, y_range[2] + y_pad)

# Create plot
plot(
  psa_results$psa_ab_averted, psa_results$psa_delta_cost,
  xlab = "Antibiotic Treated Patient Averted",
  ylab = "Incremental Cost (USD)",
  main = "",
  pch = 19,
  col = rgb(0.4, 0.3, 0.1, 0.4),   # More translucent
  cex.lab = 0.9,                   # Slightly smaller labels
  cex.axis = 0.8,                  # Slightly smaller tick labels
  xlim = xlim,
  ylim = ylim,
  bty = "n"                        # Remove box border
)

# Axes
abline(h = 0, col = "grey80", lwd = 1)
abline(v = 0, col = "grey80", lwd = 1)

# Highlight mean
points(mean_ab, mean_cost, pch = 8, col = "yellow", cex = 1.5)
text(mean_ab, mean_cost, labels = "Mean", pos = 3, cex = 0.8, font = 2, col = "yellow")

##PSA: PSEUDO-TORNADO DIAGRAMS
#Visuals for Supplemental Information S1
psa_regression_df <- cbind(psa_results, psa_samples)


##Looking at Coefficients 
##Repeat this step for whichever parameters are of interest
psa_effect_model <- lm(psa_intervention_nmb ~ AMRWeightedCosts_WithoutSevereMalaria_Intervention
                       + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention + AMRWeightedCosts_NoURS_Intervention + 
                         AMRWeightedCosts_URS_Intervention, data = psa_regression_df)

coefs <- coef(psa_effect_model)
coefs <- coefs[-1]
psatornado_parameter_names <- names(coefs)
psatornado_influences <- abs(coefs)


psatornado_df <- data.frame(
  CostParameter = psatornado_parameter_names,
  CostInfluence = psatornado_influences
)

#Order parameters by influence
psatornado_df <- psatornado_df[order(-psatornado_df$CostInfluence), ]
psatornado_df$CostParameter <- costs_friendly_names[psatornado_df$CostParameter]

ggplot(psatornado_df, aes(x = reorder(CostParameter, CostInfluence), y = CostInfluence)) +
  geom_bar(stat = "identity", fill = "peachpuff4", width = 0.4) +
  coord_flip() +
  labs(
    title = "Pseudo-Tornado Plot: Antimicrobial Resistance Costs Impact on NMB",
    x = "Cost Parameter",
    y = "Regression Coefficient (Absolute Value)"
  ) +
  theme_minimal() +  # Apply minimal theme first
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.background = element_rect(fill = "white"),  # to make background white explicitly
    plot.background = element_rect(fill = "white")    # to be safe for plot bg
  )

coefs_signed <- coef(psa_effect_model)[-1]  # keep signs

################################################################################
##EXPECTED VALUE OF PERFECT INFORMATION
################################################################################
#1: Find the maximum net monetary benefit per iteration 
#This represents the decision with perfect information
max_nmb_per_iter <- pmax(psa_results$psa_control_nmb, psa_results$psa_intervention_nmb)

#2. Find the mean net monetary benefit for both strategies across all iterations
mean_nmb_control <- mean(psa_results$psa_control_nmb)
mean_nmb_intervention <- mean(psa_results$psa_intervention_nmb)

#3. Select the maximum of those two averages: Ie. max mean net monetary benefit
#This represents current best decision under uncertainty
max_mean_nmb <- max(mean_nmb_control, mean_nmb_intervention)

#4. Calculate EVPI:
EVPI <- mean(max_nmb_per_iter) - max_mean_nmb

print(paste("Expected Value of Perfect Information (EVPI):", EVPI))
################################################################################
##SCENARIO ANALYSES
################################################################################
##OPTIMISED EDAM USE SCEANARIO
################################################################################
#Scenario Analysis: Optimised EDAM Use (Varying Compliance & Appointment Length)
#1: List of basecase prescribing probabilities and basecase intervention appointment cost
scenarios_comp <- list(
  basecase_comp = list(
    pNoSevereMalariaAB_Intervention = base_probs$pNoSevereMalariaAB_Intervention,
    pNoSevereMalariaNoAB_Intervention = base_probs$pNoSevereMalariaNoAB_Intervention,
    pNoSevereMalariaRDTAB_Intervention = base_probs$pNoSevereMalariaRDTAB_Intervention,
    pNoSevereMalariaRDTNoAB_Intervention = base_probs$pNoSevereMalariaRDTNoAB_Intervention,
    pNoURSAB_Intervention = base_probs$pNoURSAB_Intervention,
    pNoURSNoAB_Intervention = base_probs$pNoURSNoAB_Intervention,
    pURSAB_Intervention = base_probs$pURSAB_Intervention,
    pURSNoAB_Intervention = base_probs$pURSNoAB_Intervention,
    
    cAppt_Intervention = cAppt_Control * (315/100) ##Liberal estimate from basecase
  ), 

#2: List of optimised prescribing probabilities and reduced intervention appointment cost
#Based on findings from the clinical trial team, 
#Antibiotic prescribing probability should be reduced by 46% (see Methods Chapter) 
#And this 46% should shift to the likelihood of not being prescribed an antibiotic
  scenario_comp1 = list(  
    pNoSevereMalariaAB_Intervention = (.54 * base_probs$pNoSevereMalariaAB_Intervention),
    pNoSevereMalariaNoAB_Intervention = base_probs$pNoSevereMalariaNoAB_Intervention + (.46 * base_probs$pNoSevereMalariaAB_Intervention),
    pNoSevereMalariaRDTAB_Intervention = (.54 * base_probs$pNoSevereMalariaRDTAB_Intervention),
    pNoSevereMalariaRDTNoAB_Intervention = base_probs$pNoSevereMalariaRDTNoAB_Intervention + (.46 * base_probs$pNoSevereMalariaRDTAB_Intervention),
    pNoURSAB_Intervention = (.54 * base_probs$pNoURSAB_Intervention),
    pNoURSNoAB_Intervention = base_probs$pNoURSNoAB_Intervention + (.46 * base_probs$pNoURSAB_Intervention),
    pURSAB_Intervention = (.54 * base_probs$pURSAB_Intervention),
    pURSNoAB_Intervention = base_probs$pURSNoAB_Intervention + (.46 * base_probs$pURSAB_Intervention),
    
    cAppt_Intervention = cAppt_Control * (q25_pct_increase/100) ##Conservative time estimate
  ))

##Function to regenerate costs
generate_scenario_terminal_costs <- function(cAppt_Intervention){
  ##TERMINAL NODE VALUES: CONTROL ARM
  ##Nodes: Malaria Positive, No AB
  list(
  cost1 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaNoAB_Control),
  cost2 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaNoAB_Control),
  ##Nodes: Malaria Positive, AB
  cost3 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_success + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * pHospitalRef_MalariaAB_Control),
  cost4 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control) + ABWeightedCosts_WithMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithMalaria_Control + (cHospital_Admission * (pHospitalRef_MalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_MalariaAB_Control),
  
  ##Nodes: Malaria Negative, No AB
  cost5 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control) + cProductivity_success else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaNoAB_Control),
  cost6 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control) + cProductivity_failure else cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaNoAB_Control),
  ##Nodes: Malaria Negative, AB
  cost7 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_success + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * pHospitalRef_NoMalariaAB_Control),
  cost8 <- if (perspective == 'societal') cRDT + cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control) + ABWeightedCosts_WithoutMalaria_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalaria_Control else cRDT + cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalaria_Control + (cHospital_Admission * (pHospitalRef_NoMalariaAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaAB_Control),
  
  ##Nodes: No Malaria RDT, No AB
  cost9 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control) + cProductivity_success else cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTNoAB_Control),
  cost10 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control) + cProductivity_failure else cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTNoAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTNoAB_Control),
  ##Nodes: No Malaria RDT, AB
  cost11 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_success + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * pHospitalRef_NoMalariaRDTAB_Control),
  cost12 <- if (perspective == 'societal') cAppt_Control + cFollowUp + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control) + ABWeightedCosts_WithoutMalariaRDT_Control + cProductivity_failure + AMRWeightedCosts_WithoutMalariaRDT_Control else cAppt_Control + cFollowUp + ABWeightedCosts_WithoutMalariaRDT_Control + (cHospital_Admission * (pHospitalRef_NoMalariaRDTAB_Control + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoMalariaRDTAB_Control),
  
  
  ##TERMINAL NODE VALUES: INTERVENTION ARM, DANGER SIGNS PRESENT 
  ##Nodes: Severe Malaria, No AB (Nobody had severe malaria)
  cost13 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaNoAB_Intervention),
  cost14 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaNoAB_Intervention),
  ##Nodes: Severe Malaria, AB (Nobody had severe malaria, so not adding AB costs)
  cost15 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_SevereMalariaAB_Intervention),
  cost16 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_SevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_SevereMalariaAB_Intervention),
  
  ##Nodes: No Severe Malaria, No AB
  cost17 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaNoAB_Intervention),
  cost18 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaNoAB_Intervention),
  ##Nodes: No Severe Malaria, AB
  cost19 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaAB_Intervention),
  cost20 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalaria_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalaria_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaAB_Intervention),
  
  #Nodes: No Severe Malaria RDT, No AB
  cost21 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTNoAB_Intervention),
  cost22 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTNoAB_Intervention),
  ##Nodes: No Severe Malaria RDT, AB
  cost23 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_success + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * pHospitalRef_NoSevereMalariaRDTAB_Intervention),
  cost24 <- if (perspective == 'societal') cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + cProductivity_failure + AMRWeightedCosts_WithoutSevereMalariaRDT_Intervention else cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_WithoutSevereMalariaRDT_Intervention + (cHospital_Admission * (pHospitalRef_NoSevereMalariaRDTAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoSevereMalariaRDTAB_Intervention),
  
  
  ##TERMINAL NODE VALUES: INTERVENTION ARM, NO DANGER SIGNS PRESENT
  ##Nodes: UR Symptoms Present, No AB
  cost25 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * pHospitalRef_URSNoAB_Intervention),
  cost26 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + (cHospital_Admission * (pHospitalRef_URSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSNoAB_Intervention),
  ##Nodes: UR Symptoms Present, AB
  cost27 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_success + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * pHospitalRef_URSAB_Intervention),
  cost28 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + cProductivity_failure + AMRWeightedCosts_URS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + cCRP + ABWeightedCosts_URS_Intervention + (cHospital_Admission * (pHospitalRef_URSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_URSAB_Intervention),
  
  ##Nodes: UR Symptoms Not Present, No AB
  cost29 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_success else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * pHospitalRef_NoURSNoAB_Intervention),
  cost30 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention) + cBloodpressure + cPulseox + cEDAM + cProductivity_failure else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + (cHospital_Admission * (pHospitalRef_NoURSNoAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSNoAB_Intervention),
  ##Nodes: UR Symptoms Not Present, AB
  cost31 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_success + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * pHospitalRef_NoURSAB_Intervention),
  cost32 <- if (perspective == 'societal') cRDT + cAppt_Intervention + cFollowUp + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention) + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + cProductivity_failure + AMRWeightedCosts_NoURS_Intervention else cRDT + cAppt_Intervention + cFollowUp + cBloodpressure + cPulseox + cEDAM + ABWeightedCosts_NoURS_Intervention + (cHospital_Admission * (pHospitalRef_NoURSAB_Intervention + pFollowUp_HospAdmit)) + (cRepresentPHC * pRepresent_NoURSAB_Intervention)
 )
}

##Check
print(tree$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$p)

##Function to regenerate Probabilities
update_tree_scenario <- function(tree, params){
  cat("Updating tree with scenario parameters:\n")
  print(params)
  
  tryCatch({
    tree$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$p <- params$pNoSevereMalariaAB_Intervention
    cat("✅ Set pNoSevereMalariaAB_Intervention to", params$pNoSevereMalariaAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pNoSevereMalariaAB_Intervention\n")
  })
  
  tryCatch({
    tree$Intervention$`Danger Signs Present`$`Malaria Negative`$`No AB`$p <- params$pNoSevereMalariaNoAB_Intervention
    cat("✅ Set pNoSevereMalariaNoAB_Intervention to", params$pNoSevereMalariaNoAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pNoSevereMalariaNoAB_Intervention\n")
  })
  
  tryCatch({
    tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$AB$p <- params$pNoSevereMalariaRDTAB_Intervention
    cat("✅ Set pNoSevereMalariaRDTAB_Intervention to", params$pNoSevereMalariaRDTAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pNoSevereMalariaRDTAB_Intervention\n")
  })
  
  tryCatch({
    tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$`No AB`$p <- params$pNoSevereMalariaRDTNoAB_Intervention
    cat("✅ Set pNoSevereMalariaRDTNoAB_Intervention to", params$pNoSevereMalariaRDTNoAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pNoSevereMalariaRDTNoAB_Intervention\n")
  })
  
  tryCatch({
    tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$AB$p <- params$pNoURSAB_Intervention
    cat("✅ Set pNoURSAB_Intervention to", params$pNoURSAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pNoURSAB_Intervention\n")
  })
  
  tryCatch({
    tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$`No AB`$p <- params$pNoURSNoAB_Intervention
    cat("✅ Set pNoURSNoAB_Intervention to", params$pNoURSNoAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pNoURSNoAB_Intervention\n")
  })
  
  tryCatch({
    tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$AB$p <- params$pURSAB_Intervention
    cat("✅ Set pURSAB_Intervention to", params$pURSAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pURSAB_Intervention\n")
  })
  
  tryCatch({
    tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$`No AB`$p <- params$pURSNoAB_Intervention
    cat("✅ Set pURSNoAB_Intervention to", params$pURSNoAB_Intervention, "\n")
  }, error = function(e) {
    cat("❌ Could not update pURSNoAB_Intervention\n")
  })
  
  return(tree)
}

library(data.tree)

#Data frame to store results
scenario_comp_results <- data.frame(
  scenario = character(),
  scenario_cost = numeric(),
  scenario_daly = numeric(),
  stringsAsFactors = FALSE
)


#Run Scenario Analysis
scenario_name <- "scenario_comp1"
params <- scenarios_comp[[scenario_name]] #pull scenario (updated) probabilities and cost 
cAppt_Intervention <- params$cAppt_Intervention

cat("\n=== Running scenario:", scenario_name, "===\n")
print(params)

new_terminal_costs <- generate_scenario_terminal_costs(cAppt_Intervention)
new_terminal_costs_numeric <- data.frame(lapply(new_terminal_costs, function(col) {
  sapply(col, function(x) {
    if (is.data.frame(x) && nrow(x) == 1 && ncol(x) == 1) {
      as.numeric(x[[1,1]])
    } else {
      as.numeric(x)
    }
  })
}))

#Reassign terminal costs generated with the above functions
colnames(new_terminal_costs_numeric) <- paste0("cost", seq_along(new_terminal_costs_numeric))

tree$Intervention$`Danger Signs Present`$`Malaria Positive`$`No AB`$`Treatment Success`$cost <- new_terminal_costs_numeric$cost13
tree$Intervention$`Danger Signs Present`$`Malaria Positive`$`No AB`$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost14
tree$Intervention$`Danger Signs Present`$`Malaria Positive`$AB$`Treatment Success`$cost <- new_terminal_costs_numeric$cost15
tree$Intervention$`Danger Signs Present`$`Malaria Positive`$AB$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost16
tree$Intervention$`Danger Signs Present`$`Malaria Negative`$`No AB`$`Treatment Success`$cost <- new_terminal_costs_numeric$cost17
tree$Intervention$`Danger Signs Present`$`Malaria Negative`$`No AB`$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost18
tree$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$`Treatment Success`$cost <- new_terminal_costs_numeric$cost19
tree$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost20
tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$`No AB`$`Treatment Success`$cost <- new_terminal_costs_numeric$cost21
tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$`No AB`$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost22
tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$AB$`Treatment Success`$cost <- new_terminal_costs_numeric$cost23
tree$Intervention$`Danger Signs Present`$`No Malaria RDT`$AB$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost24
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$`No AB`$`Treatment Success`$cost <- new_terminal_costs_numeric$cost25
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$`No AB`$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost26
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$AB$`Treatment Success`$cost <- new_terminal_costs_numeric$cost27
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`UR Symptoms`$AB$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost28
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$`No AB`$`Treatment Success`$cost <- new_terminal_costs_numeric$cost29
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$`No AB`$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost30
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$AB$`Treatment Success`$cost <- new_terminal_costs_numeric$cost31
tree$Intervention$`No Danger Signs Present`$`Symptoms Assessment to Malaria RDT`$`No UR Symptoms`$AB$`Treatment Failure`$cost <- new_terminal_costs_numeric$cost32


#Clone base case tree
tree_scenario <- tree$clone()

# Reset all node$p to default base probs
tree_scenario$Do(function(node) {
  if (!node$isRoot && !is.null(node$prob)) {
    node$p <- node$prob
  }
})

cat("Before update - base p value for Malaria Negative AB:\n")
print(tree_scenario$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$p)

#Applying scenario updates
tree_scenario <- update_tree_scenario(tree_scenario, params)

# Confirm one update
cat("After update - new p value for Malaria Negative AB:\n")
print(tree_scenario$Intervention$`Danger Signs Present`$`Malaria Negative`$AB$p)

# Calculate expected values
CalculateExpectedValues(tree_scenario)

# Extract result
scenario_cost <- tree_scenario$Intervention$expected_cost
scenario_daly <- tree_scenario$Intervention$expected_daly

# Store results in data frame
scenario_comp_results <- rbind(scenario_comp_results, data.frame(
  scenario = scenario_name,
  scenario_cost = scenario_cost,
  scenario_daly = scenario_daly,
  stringsAsFactors = FALSE
))

# Calculate scenario analysis ICER
scenario_ICER <- ((scenario_cost - cost_control)/ (daly_control - scenario_daly))

scenario_delta_cost <- (scenario_cost - cost_control)
scenario_delta_daly <- (daly_control - scenario_daly)
################################################################################
##Scenario Analysis: Compliance Part 2 -- Antibiotic Prescriptions 
#Hard coded antibiotic counts from AB Type data frames)
scenario_base_ab_counts <- c(
  SevMalariaNeg_Intervention_AB = 556,
  SevMalariaNoRDT_Intervention_AB = 194,
  URS_Intervention_AB = 520,
  NoURS_Intervention_AB = 34
)

#Hard code base case probability of receiving an antibiotic 
#From probabilities in data cleaning stage
scenario_ab_nodes_base <- list(
  SevMalariaNeg_Intervention_AB = 0.5978495,
  SevMalariaNoRDT_Intervention_AB = 0.5215054,
  URS_Intervention_AB =  0.6046512,
  NoUR_Intervention_AB = 0.2098765
)
##Hard code each probability reduced by 46%
scenario_ab_nodes_scenario <- list(
  SevMalariaNeg_Intervention_AB = 0.3228387,
  SevMalariaNoRDT_Intervention_AB = 0.2816129,
  URS_Intervention_AB = 0.3265116,
  NoUR_Intervention_AB = 0.1133333
)

# Calculate percentage change
scenario_AB_percentage_change <- mapply(function(base, scenario) {
  ((scenario - base) / base) * 100
}, base = scenario_ab_nodes_base, scenario = scenario_ab_nodes_scenario)

# Convert to list
percentage_change_list <- as.list(scenario_AB_percentage_change)

# Print result
percentage_change_list

# Compute scenario AB prescriptions by scaling base counts
scenario_ab_counts <- scenario_base_ab_counts * (1 + (scenario_AB_percentage_change / 100))
# Calculate delta from base
delta_ab_counts <- scenario_ab_counts - scenario_base_ab_counts
delta_total_ab <- sum(delta_ab_counts)

# Print results
cat("\nAB prescription change by node:\n")
print(delta_ab_counts)
cat("\nPercentage change in AB prescriptions by node:\n")
print(round(scenario_AB_percentage_change, 2))
cat("Total change in ABs:", sum(delta_ab_counts), "\n")

ab_change_df <- data.frame(
  Node = names(delta_ab_counts),
  AB_Change = as.numeric(delta_ab_counts)
)

scenario_delta_prescriptions <- prescriptions_averted - delta_total_ab
scenario_ICER_ABs <- scenario_delta_cost / scenario_delta_prescriptions
################################################################################
##Scenario Analysis: Visualisations

#Visual friendly names
ab_friendly_names <- c(
  SevMalariaNeg_Intervention_AB = "Danger Signs, Severe Malaria Negative",
  SevMalariaNoRDT_Intervention_AB = "Danger Signs, No Malaria RDT",
  URS_Intervention_AB = "No Danger Signs, UR Symptoms Present",
  NoURS_Intervention_AB = "No Danger Signs, No UR Symptoms",
  NoUR_Intervention_AB = "No Danger Signs, No UR Symptoms"
)


##Visualisation 1: Main Report Chapter 3
ab_waterfall <- ab_change_df %>%
  mutate(
    Type = ifelse(AB_Change >= 0, "Increase", "Reduction"),
    FriendlyNode = recode(Node, !!!ab_friendly_names)
  )

ggplot(ab_waterfall, aes(x = reorder(FriendlyNode, -AB_Change), y = AB_Change, fill = Type)) +
  geom_col(width = 0.3) +  # thinner bars
  coord_flip() +
  scale_fill_manual(values = c("Reduction" = "peachpuff4", "Increase" = "peachpuff")) +
  labs(
       y = "Change in Number of Antibiotic Treated Patients",
       x = "Node") +
  theme_minimal() + 
  theme(plot.title = element_text(face = "italic", size = 14))


##Visualisations for Supplemental Information S1
##Visualisation 2
scenario_percentage_df <- data.frame(
  Node = names(scenario_AB_percentage_change),
  PercentChange = as.numeric(scenario_AB_percentage_change)
) %>%
  mutate(FriendlyNode = recode(Node, !!!ab_friendly_names))
ggplot(scenario_percentage_df, aes(x = reorder(FriendlyNode, PercentChange), y = PercentChange)) +
  geom_col(fill = "peachpuff4") +
  coord_flip() +
  labs(title = "% Change in Antibiotic Treated Patients Prescriptions by Node",
       x = "Node",
       y = "Percent Change (%)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "italic", size = 14))

##Visualisation 3
scenario_icer_data <- tibble::tibble(
  scenario = c("Base case", "Scenario"),
  delta_cost = c(delta_cost$percentage, scenario_delta_cost$percentage),
  delta_daly = c(delta_daly, scenario_delta_daly$percentage),
  prescriptions_averted = c(prescriptions_averted, scenario_delta_prescriptions)
) %>%
  mutate(
    ICER_daly = delta_cost / delta_daly,
    ICER_rx = delta_cost / prescriptions_averted
  )

ggplot(scenario_icer_data, aes(x = ICER_rx, y = ICER_daly, color = scenario)) +
  geom_point(size = 4) +
  labs(
    title = "Scenario Comparison of Incremental Cost-Effectiveness",
    x = "Incremental Cost per Antibiotic Treated Patient Averted (USD)",
    y = "Incremental Cost per DALY (USD)",
    color = "Analysis Scenario"
  ) +
  scale_color_manual(
    values = c("Base case" = "peachpuff2", "Scenario" = "peachpuff4")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )
################################################################################
##STUDY DESIGN AND CLUSTER EFFECT
################################################################################
#Fit the generalized linear mixed model (GLMM)
EDAM_glm_model <- glmer(
  primary_outcome ~ arm + agegroup + sex + (1 | phc_name),
  data = clinicaldata,
  family = binomial(link = "logit")
)

#View the model summary
summary(EDAM_glm_model)

#Exponentiate fixed effects to get odds ratios
odds_ratios <- exp(fixef(EDAM_glm_model))
print(odds_ratios)

m <- mean(table(clinicaldata$phc_name))

#ICC calculation
var_comp <- as.data.frame(VarCorr(EDAM_glm_model))
icc <- var_comp$vcov[1] / (var_comp$vcov[1] + (pi^2 / 3))
print(icc)

#Design effect
deff <- 1 + (m - 1) * icc
print(deff)

##Boostrapping, using variance applied as noise to PSA draws, to calculate ICER 95% CI:
n_boot <- 1000
sd_clustering <- sqrt(var_comp$vcov[1])

#Set seed
set.seed(100)

#Perform one bootstrap replicate
bootstrap_icer <- function() {
  sampled_row <- psa_results[sample(nrow(psa_results), 1, replace = TRUE), ]
  
  # Add clustering noise to the effect
  boot_delta_effect <- sampled_row$psa_ab_averted + rnorm(1, mean = 0, sd = sd_clustering)
  
  boot_delta_cost <- sampled_row$psa_intervention_cost - sampled_row$psa_control_cost
  
  # Return ICER (only if effect is non-zero)
  if (is.finite(boot_delta_effect) && boot_delta_effect != 0) {
    return(boot_delta_cost / boot_delta_effect)
  } else {
    return(NA)
  }
}

#Run the bootstrapping
boot_icers <- replicate(n_boot, bootstrap_icer())

# Remove NA and infinite results
boot_icers_clean <- boot_icers[is.finite(boot_icers)]

# Compute 95% confidence interval
boot_icer_ci <- quantile(boot_icers_clean, probs = c(0.025, 0.975), na.rm = TRUE)
print(boot_icer_ci)

#Histogram for Main Report Chapter 3
hist_data <- hist(boot_icers_clean,
                  breaks = 50,
                  xlab = "ICER",
                  main = "",
                  col = "peachpuff4",
                  cex.lab = 0.8,
                  cex.axis = 0.75,
                  border = "white")

# Add vertical lines for confidence interval
abline(v = boot_icer_ci, col = "forestgreen", lty = 2, lwd = 1)

# Add 95% CI labels at the top of the histogram
y_max <- max(hist_data$counts)

text(x = boot_icer_ci[1],
     y = y_max * 0.95,
     labels = "Lower 95% CI",
     col = "forestgreen",
     pos = 2,    # left of the line
     cex = 0.75,
     font = 2)   # bold

text(x = boot_icer_ci[2],
     y = y_max * 0.95,
     labels = "Upper 95% CI",
     col = "forestgreen",
     pos = 4,    # right of the line
     cex = 0.75,
     font = 2)   # bold

################################################################################
##END OF COST-EFFECTIVENESS ANALYSIS 
##THANK YOU FOR READING




