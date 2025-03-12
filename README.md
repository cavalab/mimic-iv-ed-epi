# mimic-iv-ed
A repository of code used to analyse the MIMIC-IV-ED dataset

These files, run in order, reproduce the analysis we've used to analyse the MIMIC-IV-ED dataset. They assume you have downloaded the dataset to your local machine from Physionet, along with MIMIC-IV.

1. "edstays-extract-time-of-day-and-previous-visits.r": extracts time spent in ER, time of arrival, year group associated with patient, and visit history.
2. "edstays_triage_data_cleaning_mult_visits.r": cleans and categorises vital signs, pain scores, and demographic variables.
3. "mimic-chief-complaint.r": translates the free-text "chief complaint" field to a series of 54 binary indicators associated with each visit.
4. "mimic-chief-complaint.r": obtains patient age at each visit from the MIMIC-IV dataset, and translates ICD codes associated with each stay to a Charlson score. Also obtains insurance, marital, and language info for visits resulting in an admission.
5. "edstays-triage-hists.r": produces histograms showing the distribution of variables across the cohort, allowing at-a-glance identification of reference groups.
6. "mimic_mult_visit_with_vitals_binary_recode.r": reformat data on all visits to ER to a series of binary indicators.
7. "adm_edstays_binary_recode.r": reformat data on all visits to ER *resulting in an admission* to a series of binary indicators, including additional covariates.
8. "edstays_add_vitalsigns_during_stay.r": attach data on vital signs from readings taken during a patient's stay, cleaned, categorised, and reformat as binary indicators.

All other files (starting "ps2") are named according to their exposure and outcome. These fit propensity scores for a given exposure variable, calculate odds ratios using two PS methodologies (covariate adjustment and matching with McNemar's test), save them, and produce OR plots (for categorical exposures).


