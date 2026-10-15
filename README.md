# ancova_mismatch

Code, frozen protocol, and saved results for a simulation study of eight analyses of a categorical factor's effect
with a numeric covariate. Software: R 4.3.3 with car 3.1.2. Run the lines below in R from this folder.

Tables and figures from the saved results:

    Sys.setenv(SIM_OUTPUT_ROOT = "output")   # figures: TIFF at 600 dpi; add FIG_FORMAT = "png" for PNG files
    source("make_tables_ancova.R"); source("make_figures_ancova.R")
    system("python score_extension.py output/ancova_ext/ancova_results.csv")

Analytical calculation and checks:

    Sys.setenv(TC_NDES = "40", TC_NX = "300"); source("theory_checks.R")
    source("theory_limits.R"); source("theory_independent_check.R"); source("validate_ancova.R")
    source("export_for_python.R"); system("python validate_python_ancova.py")   # statsmodels 0.15.0, Python 3.13.16

Rerun the simulations (hours, fresh R session for each):

    Sys.setenv(SIM_OUTPUT_ROOT = "output", SIM_MODE = "full", SIM_SEED = "8419"); source("explore_ancova.R")
    Sys.setenv(SIM_OUTPUT_ROOT = "output", SC_REPS = "500", SC_DESIGNS = "20", SC_SEED = "9300000"); source("scaling_ancova.R")
    Sys.setenv(SIM_OUTPUT_ROOT = "output", SIM_MODE = "ext", SIM_SEED = "7207"); source("explore_ancova_ext.R")

`docs/PROTOCOL_ancova_v1.md` is the frozen protocol, with the criteria and the md5 values of the three simulation
scripts. `output/` holds the saved results.
