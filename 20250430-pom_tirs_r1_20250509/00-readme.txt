- "./run_08s0_pval.bat":
  - this batch file runs the script in "./08-scripts_oc_pval/*.r"
    to simulate data and analyze them with POM

===
- "./08-scripts_oc_pval":
  - include all simulation scripts
  - "08s0-oc_genda_1_xi_1_a_1_po_*.r" is for Z random with different
    K = 5, 10, and 30
  - "08s0-oc_genda_2_xi_1_a_1_po_*.r" is for Z|U,a with different
    K = 5, 10, and 30

- "./code":
  - include all POM functions
  - see "./code/00-readme.txt" for a brief documentation 

- "./log_pval":
  - the batch file would dump messages to this directory 

- "./results_oc_pval":
  - the batch file would save results to this directory 


===
- For manuscript
  - "./paper_scripts_v4":
    - the script under this directory would generate figures for the POM paper

  - "./paper_plots_v4":
    - all figures would be saved to this directory


===
- For additional configurations (ntrt, noise, and ps)
  - "./paper_supp_scripts_v4"
  - "./paper_supp_plots_v4"

  - with different sample sizes of single-arm study
    - "./run_08s0_pval_ntrt.bat"
    - "./08-scripts_oc_pval_ntrt"
    - "./log_pval_ntrt"
    - "./results_oc_pval_ntrt"

  - with larger outcome noises (error sd = 2) of single-arm study
    - "./run_08s0_pval_noise.bat"
    - "./08-scripts_oc_pval_noise"
    - "./log_pval_noise"
    - "./results_oc_pval_noise"

  - with PS weighting and matching
    - "./run_08s0_pval_ps.bat"
    - "./08-scripts_oc_pval_ps"
    - "./log_pval_ps"
    - "./results_oc_pval_ps"


