# Tutorial on building, calibrating, and applying a discrete-event simulation cancer model

## Simulating ground truth data
To generate oracle data, run `ground_truth/simulate_truth.R`. The only data file needed as an input is `data/background_mortality_data.xlsx`, which is based on the 1980 U.S. male and female life tables from the Human Mortality Database. The model includes healthy (H), preclinical cancer (P), and clinical cancer (C), and death (D) states, with substates for early versus late-stage cancer and death from cancer or from other causes. The outputs include:
- decision model inputs
  - relative survival by stage at cancer diagnosis (`relative_survival_cancer.csv`)
- parameter priors for calibration
  - uniform prior distributions for the unknown parameters (`priors.RData`)
- calibration targets
  - preclinical cancer incidence (`prevalence_asymptomatic_cancer.csv`)
  - incidence of symptomatically-detected cancer without screening (`incidence_symptomatic_cancer.csv`)
  - stage distribution of symptomatically-detected cancer at diagnosis (`stage_distr.csv`).

## Running the model
Currently, the function `load_default_params()` in `R/01_load_model_inputs.R` can take two sequences of states in `v_states` as inputs:
- c('H', 'P', 'C', 'D')
- c('H', 'L', 'P', 'C', 'D')

where 'H' denotes healthy, 'L' precancerous lesion, 'P' preclinical cancer, 'C' clinical cancer, and 'D' death. If the input vector for `v_states` includes 'L', the decision model functions in `R/02_decision_model_functions.R` will model the transition from H to L, then from L to P. Otherwise, simply the transition from H to P will be modeled. For both cases, the next step will be to model the progression of cancer stages in the order of `v_cancer`, which should be inputted in chronological order and include the stages listed in the relative survival data (`relative_survival_cancer.csv`). For each cancer stage before the last stage, the time from preclinical cancer stage i to the next stage will be simulated and compared to the simulated time from preclinical cancer stage i to detection. In the non-screening scenario, the stage of diagnosis will be the first stage at which the time to detection is less than the time to progression to the next stage. If the last preclinical cancer stage is reached, we simply simulate the time to detection. From detection, the time to death from cancer based on the survival curve for the stage at detection. Finally, the time and cause of death is calculated as the earliest of death from cancer and death from other causes.

Add: Updating parameters with `update_param_list()`, `make_param_map()`, and `update_param_from_map()`
