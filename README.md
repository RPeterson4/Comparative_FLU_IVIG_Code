# Comparative_FLU_IVIG_Code

Code for second FLU-IVIG paper, comparing the ordinal endpoint and proportional odds model to other models fitted to other endpoints
under simulation.

File descriptions:

analytic_prop_odds: Code used to analytically derive IVIG group distributions on day 7 with odds ratios that approximate the pre-
specified value of 1.77. Defaulted to FLU-IVIG placebo group and its seven treatment effects. Compiled HTML file displays results of
P1-T1.

simulation_prop_odds: Code used to run the simulation of eight models fitted to six endpoints. Using eight cores takes about 2 hours 
per simulation. Compiled HTML file displays results of P1-T1.

data_generate_prop_odds: Code of various functions needed to generate data and calculate power. In order to simulate different placebo
group distributions and treatment effects in simulation_prop_odds, transition matrices and day 0 placebo group distributions must be 
modified in this file.

