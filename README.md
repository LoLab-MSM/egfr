# egfr
PySB models of EGFR receptor driven signaling.  

Four model versions are included:
chen_sorger_2009_egfr: A PySB version of ErbB signaling through phosphorylated Erk and Akt, intended to recreate the model 
originally published in Chen et al (2009).

chen_sorger_2009_egfr_simplified: A somewhat simplified version of the above; most simplifications are in the parameters used 
for internalization of receptors.

ems_egfr_plus_mtor: The simplified version of EGFR signaling expanded to include signaling through mTOR.  Also includes several
different topologies for ErbB2-driven downstream signaling.

ems_egfr_plus_apoptosis: Includes ErbB signaling through Erk and Akt, mTOR signaling, intrinsic apoptosis signaling, and 
crosstalk events between the two.
