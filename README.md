
# Replication files for Inference Based on Time-Varying SVARs Identified with Sign Restrictions
by Jonas E. Arias, Juan F. Rubio-Ramirez, Minchul Shin, and Daniel F. Waggoner

Here is the link to a local copy of the paper, [[Local copy]](https://mcmcs.github.io/papers/TV_SVAR.pdf)


***Important:*** Estimating our model takes a relatively long time. For reference, the estimation was performed using the following hardware specifications: Architecture: x86_64, CPU op-mode(s): 32-bit, 64-bit. Byte
Order: Little Endian. CPU(s): 32. On-line CPU(s) list 0-15. Off-line CPU(s) list:
16-31. Thread(s) per core: 1. Core(s) per socket: 16. Model Name: Intel(R) Xeon (R)
Platinum 8488C. CPU MHz: 2400. BogoMIPS: 4800.00. We use the following software:
Red Hat Enterprise Linux 8.10 and MATLAB 2024a. Under the baseline specification
this code takes about 5 days to run.

To save users time, we also provide pre-estimated parameters (i.e., MCMC draws), which facilitate generating figures and tables directly. To skip the estimation step, please download the following ***four*** files:

[[File 1]](https://www.dropbox.com/scl/fi/0j7tzi5hbglsuom2st26e/temp_results_jointpr_2s_rest97_sam1r997ndraws1000000.mat?rlkey=q7hiclkka5b3k2nzrxi99k1id&dl=0), [[File 2]](https://www.dropbox.com/scl/fi/sqpaw61euvtcu21duz1f1/temp_results_jointpr_rest115r3024ndraws1000000.mat?rlkey=s2zkzmcu3xakxx5q2cow4dnja&dl=0), [[File 3]](https://www.dropbox.com/scl/fi/v8498sipdhcydaz9n16yv/temp_results_jointpr_rest118r997ndraws500000_restALL.mat?rlkey=gyeahkmmombsi2b6235s2lxu5&dl=0), [[File 4]](https://www.dropbox.com/scl/fi/bxohvl7n6xcnegffb4j7d/temp_results_unrest_ndraws_50000_is_restrict_B0.mat?rlkey=5c1l1ylrbxixanpx1o5o1jdwr&dl=0)

and place them in the folder: `Baseline\results\`

***Actual replication steps:***
To replicate the results in the paper proceed as follows:

1. Run the m-file `Baseline\main_jointpr_taylor.m` This file obtains the draws required to produce the results in the paper. 

    Users interested in replicating the figures and the tables in the paper without producing the draws can skip this step. For users' convenience, the draws can be found here: `Baseline\results\temp_results_jointpr_rest115r3024ndraws1000000.mat`

2. Run the m-file `Baseline\replicate_figures_and_table.m` to reproduce the main results of the paper. This file produces Figures 1-6 and Table 1.

3. Run the m-file `Baseline\replicate_figures_in_appendix_part_1.m` to reporduce the figures III.1, III.6, III.7, and III.8. To replicate the figures you must previously run the following two m-files:
    1. `Baseline\main_jointpr_taylor_ood1_R2ALL.m`
    2. `constant_parameters\code\run_mainfile.m`

    Users interested in replicating the figures and the table in the paper without producing the draws can skip these two steps. For users' convenience, the draws can be found here:

    `Baseline\results\temp_results_jointpr_rest118r997ndraws500000_restALL.mat` and `constant_parameters\code\results\results.mat`, respectively.

4. Run the m-file: `Baseline\replicate_figures_in_appendix_part_2.m` to reproduce the figures III.2, III.3, III.4, and III.5.

    To replicate the figures you must previously run the following m-file `Baseline\main_jointpr_2MPshocks_cspi_ood.m`

    Users interested in replicating the figures and the table in the paper without producing the draws can skip this step. For users' convenience, the draws can be found here: `Baseline\results\temp_results_jointpr_2s_rest97_sam1r997ndraws1000000.mat`

5. Run the m-file `Baseline\replicate_figure_V1.m` to reproduce the figure V.1.

    To replicate the figures you must previously run the following m-file: `Baseline\main_jointpr_taylor.m`.

    Users interested in replicating the figures and the table in the paper without rpoducing the draws can skip this step and can use the file described on step 1.

6. Run the m-file `Baseline\online_appendix_V\Baseline_simulated_data\replicate_figure_V2.m` to reproduce the figure V.2.
    To replicate the figures you must previusly run the following m-file: `Baseline\online_appendix_V\Baseline_simulated_data\main_jointpr_taylor.m`

    Users interested in replicating the figures and the table in the paper without producing the draws can skip this step by downloading the following file:

    [[Link to download result file for this step]](https://www.dropbox.com/scl/fi/ya3qqkkdc2rbgsop7m88k/temp_results_jointpr_rest115r3024ndraws200000.mat?rlkey=ce7taeoupvj4pb5g59ji08xjm&dl=0)

    and place it in the following folder: `Baseline\online_appendix_V\Baseline_simulated_data\results\`