from pysb import *
from collections import OrderedDict

# Parameters for EGFR Model - A431 cells
# Obtained from Chen-Sorger Jacobian files unless otherwise specified.
# Commented values prefixed by 'c' are the variables names from Chen Sorger Jacobian files.
parameter_dict = OrderedDict([
    #('erbb1_0',
     #Parameter('erbb1_0', 1.08e6)
     #),
    ('initial_amounts', # Initial values for all starting species (in molecules/cell)
     [Parameter('erbb1_0', 1.08e6), #531
      Parameter('erbb2_0', 4.62e5), #c141
      Parameter('erbb3_0', 6.23e3), #c140
      Parameter('erbb4_0', 7.94e2), #c143
      Parameter('ATP_0',   1.2e9), #c105
      Parameter('DEP_0',   7e4), #c280
      Parameter('CBL_0',   5e3), #c12
      Parameter('GAP_0', 5.35e5), #c14
      Parameter('SHC_0', 1.1e6), #c31
      # Parameter('SHCPase_0', 1000)
      Parameter('GRB2_0', 1264), #c22
      #Parameter('SOS_0', 6.63e4) removed to better represent Chen/Sorger model
      Parameter('RAS_0', 5.81e4), #c26
      Parameter('RAF_0', 7.11e4), #c41
      Parameter('BRAF_0', 7.11e4), #Set to Raf amount
      Parameter('CRAF_0', 7.11e4), #Set to Raf amount
      Parameter('KSR_0', 7.11e4), #Set to Raf amount
      Parameter('MEK_0', 3.02e6), #c47
      Parameter('ERK_0', 6.95e5), #c55
      Parameter('PP1_0', 5e4), #c44
      Parameter('PP2_0', 1.25e5), #c53
      Parameter('PP3_0', 1.69e4), #c60
      Parameter('GRB2_SOS_0', 8.89e7), #This added to better represent Chen Sorger model c30
      Parameter('GAB1_0', 94868.3), #c426
      Parameter('PI3K_0', 3.55656e7), #c287 c455?
      Parameter('SHP2_0', 1e6), #c463
      Parameter('PIP_0',     3.94e5), #c444
      Parameter('PTEN_0',    5.62e4), #c279
      Parameter('SHP_0',     2.21e3), #c461
      Parameter('AKT_0',     9.05e5), #c107
      Parameter('PDK1_0',     3.00416e8), #c109
      Parameter('PP2A_III_0', 4.5e5), #c113
      Parameter('Pase9t_0', 0), #c521
      Parameter('ERKPP_0', 0),
      Parameter('AKTPP_0', 0),
      Parameter('mTOR_0', 2095), #RPPA
      Parameter('TORC1_ptns_0', 2095), # These and mTORC2 ptns have been set to same as mTOR to make mTOR the limiting factor in complex formation (NOT sure if this is actually the case).
      Parameter('TORC2_ptns_0', 2095), # See comment above for mTORC1
      Parameter('AMPK_0', 846), #RPPA
      Parameter('LKB1_0', 1910), #RPPA
      Parameter('TSC_0', 1570), #Based on relative protein amounts in Pezze et al 2012 (AKT amount is 14x the TSC amount -- done in HeLa cells)
      Parameter('S6K_0', 2.3e8), #RPPA - p70S6K
      Parameter('Rheb_0', 1500), #Not sure what to make this -- set to TSC level so that it can be completely inhibited under the right circumstances
      Parameter('FKBP38_0', 6000), #Calculated based on absolute prostate cancer protein values at MOPED ~3x FKBP38 as mTOR (used relative to our RPPA measurement).
      Parameter('rpS6_0', 28742), #RPPA
      Parameter('EIF4EBP1_0', 1835), #RPPA
      Parameter('EIF4E_0', 3600), #Calculated based on 2:1 EIF4E:EIF4EBP1 ratio in HEK293 cells (MOPED) -- Not cancer cell line so not sure how true this holds for other cell lines.
      Parameter('RSK1_0', 1666), #Calculated based on relative amounts of EIF4EBP1 and RSK1 in HEK293 cells -- see caveat above.
      Parameter('ELK1_0', 80000), #Calculated based on relative amounts of ERK in model/HEK293 cells and ELK1 in HEK293 cells -- see caveat above.
      #These parameters taken from the EARM lopez Embedded model (except Bim and Puma)
      Parameter('Bid_0'   , 4.0e4), # Bid
      Parameter('BclxL_0' , 2.0e4), # cytosolic BclxL
      Parameter('Mcl1_0'  , 2.0e4), # Mitochondrial Mcl1
      Parameter('Bcl2_0'  , 2.0e4), # Mitochondrial Bcl2
      Parameter('Bad_0'   , 1.0e3), # Bad
      Parameter('Noxa_0'  , 1.0e3), # Noxa
      Parameter('CytoC_0' , 5.0e5), # cytochrome c
      Parameter('Smac_0'  , 1.0e5), # Smac
      Parameter('Bax_0'   , 0.8e5), # Bax
      Parameter('Bak_0'   , 0.2e5), # Bak
      Parameter('Bim_0', 1.0e3),
      Parameter('Puma_0', 1.0e3),
      Parameter('C8_0'    , 2.0e4), # procaspase-8
      Parameter('FOXO_0', 1e2)
    ]),
    # Parameters ('k' prefixed variables are Chen-Sorger variable names from Jacobian files):
    # Receptor-level rate parameters:
    ('erbb1_0',
     Parameter('erbb1_0_init', 1.08e6)),
    ('EGF_bind_ErbB1',
     [Parameter('EGF_bind_ErbB1kf', 1e7), #k1
      Parameter('EGF_bind_ErbB1kr', .0033) #kd1
      ]),
    ('EGF_bind_ErbB1d',
     [Parameter('EGF_bind_ErbB1dkf', 1e7), #k1
      Parameter('EGF_bind_ErbB1dkr', .0033) #kd1
      ]),
    ('HRG_bind_ErbB3',
     [Parameter('HRG_bind_ErbB3kf', 1e7), #k119
      Parameter('HRG_bind_ErbB3kr', .0103115) #kd119
      ]),
    ('HRG_bind_ErbB4',
     [Parameter('HRG_bind_ErbB4kf', 1e7), #k119
      Parameter('HRG_bind_ErbB4kr', .0103115) #kd119
      ]),
    ('EGF_bind_ErbB2_ErbB3',
     [Parameter('EGF_bind_ErbB2_ErbB3kf', 800), #k1c
      Parameter('EGF_bind_ErbB2_ErbB3kr', 1) #kd1c
      ]),
    ('EGF_bind_ErbB2_ErbB4',
     [Parameter('EGF_bind_ErbB2_ErbB4kf', 518), #k1d
      Parameter('EGF_bind_ErbB2_ErbB4kr', 1e-1), #kd1d
      ]),
    ('EGFE_bind_ErbBE',
     [Parameter('EGFE_bind_ErbBEkf', 5.426e-2), #k10b
      Parameter('EGFE_bind_ErbBEkr', 1.1e-2) #kd10
      ]),
    ('HRGE_bind_ErbBE',
     [Parameter('HRGE_bind_ErbBEkf', 5.426e-2),
      Parameter('HRGE_bind_ErbBEkr', 1.1e-2)
      ]),
    ('kdeg_HRG',
     Parameter('kdeg_HRG', 5.7e-4)
     ),
    ('ErbB1_bind_ErbB1L',
     [Parameter('ErbB1_bind_ErbB1Lkf', 7.44622e-6), #k2
      Parameter('ErbB1_bind_ErbB1Lkr', 1.6e-1) #kd2
      ]),
    ('ErbB1_bind_ErbB1',
     [Parameter('ErbB1_bind_ErbB1kf', 7.44622e-6), #k2
      Parameter('ErbB1_bind_ErbB1kr', 1.6e-1) #kd2
      ]),
    ('ErbB1_bind_ErbB2',
     [Parameter('ErbB1_bind_ErbB2kf', 3.73632e-8), #k2b 
      Parameter('ErbB1_bind_ErbB2kr', 1.6e-2) #kd2b
      ]),
    ('ErbB1_bind_ErbB3',
     [Parameter('ErbB1_bind_ErbB3kf', 3.73632e-8), #k2b
      Parameter('ErbB1_bind_ErbB3kr', 1.6e-2) #kd2b
      ]),
    ('ErbB1_bind_ErbB4',
     [Parameter('ErbB1_bind_ErbB4kf', 3.73632e-8), #k2b
      Parameter('ErbB1_bind_ErbB4kr', 1.6e-2) #kd2b
      ]),
    ('ErbB1L_bind_ErbB1L',
     [Parameter('ErbB1L_bind_ErbB1Lkf', 7.44622e-6), #k2
      Parameter('ErbB1L_bind_ErbB1Lkr', 1.6e-1) #kd2
      ]),
    ('ErbB1L_bind_ErbB2',
     [Parameter('ErbB1L_bind_ErbB2kf', 3.73632e-8), #k2b 
      Parameter('ErbB1L_bind_ErbB2kr', 1.6e-2) #kd2b
      ]),
    ('ErbB2_bind_ErbB2',
     [Parameter('ErbB2_bind_ErbB2kf', 8.36983e-9), #k103
      Parameter('ErbB2_bind_ErbB2kr', 1.6e-2) #kd103
      ]),
    ('ErbB1L_bind_ErbB3',
     [Parameter('ErbB1L_bind_ErbB3kf', 3.73632e-8), #k2b
      Parameter('ErbB1L_bind_ErbB3kr', 1.6e-2) #kd2b
      ]),
    ('ErbB2_bind_ErbB3',
     [Parameter('ErbB2_bind_ErbB3kf', 1.48131e-8), #k120
      Parameter('ErbB2_bind_ErbB3kr', 1e-1) #kd120
      ]),
    ('ErbB1L_bind_ErbB4',
     [Parameter('ErbB1L_bind_ErbB4kf', 3.73632e-8), #k2b
      Parameter('ErbB1L_bind_ErbB4kr', 1.6e-2) #kd2b
      ]),
    ('ErbB2_bind_ErbB4',
     [Parameter('ErbB2_bind_ErbB4kf', 1.48131e-8), #k120
      Parameter('ErbB2_bind_ErbB4kr', 1e-1) #kd120
      ]),
    ('ErbB1_bind_ATP',
     [Parameter('ErbB1_bind_ATPkf', 1.8704e-8), #k122
      Parameter('ErbB1_bind_ATPkr', 1), #kd122
      ]),
    ('ErbB2_bind_ATP',
     [Parameter('ErbB2_bind_ATPkf', 1.8704e-8),
      Parameter('ErbB2_bind_ATPkr', 1) #kd122
      ]),
    ('ErbB23_bind_ATP',
     [Parameter('ErbB3_bind_ATPkf', 1.8704e-8), #k122
      Parameter('ErbB3_bind_ATPkr', 1) #kd122
      ]),
    ('ErbB4_bind_ATP',
     [Parameter('ErbB4_bind_ATPkf', 1.8704e-8), #k122
      Parameter('ErbB4_bind_ATPkr', 1) #kd122
      ]),
      #All unphosphorylated ErbB dimers binding DEP rate constants are set to 0 in Jacobian files (k95). 
    ('ErbBP1_bind_DEP',
     [Parameter('ErbBP1_bind_DEPkf', 5e-5), #k94 or k94b (equal in Jacobian files)
      Parameter('ErbBP1_bind_DEPkr', 1e-2) #kd94
      ]),
    ('ErbBP2_bind_DEP',
     [Parameter('ErbBP2_bind_DEPkf', 5e-5), #k94 or k94b (equal in Jacobian files)
      Parameter('ErbBP2_bind_DEPkr', 1e-2) #kd94
      ]),
    ('ErbBP4_bind_DEP',
     [Parameter('ErbBP4_bind_DEPkf', 5e-5), #k94 or k94b (equal in Jacobian files)
      Parameter('ErbBP4_bind_DEPkr', 1e-2) #kd94
      ]),
      #All phosphorylated ErbB dimers binding ATP rate constants are set to 0 in Jacobian files (k123).
    ('ATP_phos_ErbB', #kd123
     .177828),
    ('DEP_dephos_ErbB', #kd95
     33),
    ('ErbB1P_ErbBXP_bind',
     [Parameter('ErbB1P_ErbBXP_bindkf', 5e-7), #k102
      Parameter('ErbB1P_ErbBXP_bindkr', 5.61009) #kd102
      ]),
    ('ErbB2P_ErbBXP_bind',
     [Parameter('ErbB2P_ErbBXP_bindkf', 8.36983e-9), #k103
      Parameter('ErbB2P_ErbBXP_bindkr', 1.6e-2) #kd103
      ]),
    ('ErbB2_lateralsignal',
      8.36983e-9), #k103
    ('kint_const_1',
     Parameter('kint_const_1', 5e-6)
     ),
    ('kint_const_2',
     Parameter('kint_const_2', 5e-7)
     ),
    ('kint_const_3',
     Parameter('kint_const_3', 1e-6)
     ),
    ('kint_const_4',
     Parameter('kint_const_4', 1e-6)
     ),
    ('ErbB1_Cbl_Y1045_Cyto',
     [Parameter('ErbB1_Cbl_Y1045_Cyto_kf', 6.73e-6),
      Parameter('ErbB1_Cbl_Y1045_Cyto_kr', 1.66e-4)
      ]),
    ('ErbB1_Cbl_Y1045_Endo',
     [Parameter('ErbB1_Cbl_Y1045_Endo_kf', 1e-15),
      Parameter('ErbB1_Cbl_Y1045_Endo_kr', 8.0833e-3)
      ]),
    ('ErbB1_Cbl_Grb2',
     [Parameter('ErbB1_Cbl_Grb2_kf', 6.73e-6),
      Parameter('ErbB1_Cbl_Grb2_kr', 1.66e-4)
      ]),
    ('kint_no_cPP_1', 
     [Parameter('kint_no_cPP_1kf', .013), #k6
      Parameter('kint_no_cPP_1kr', 5e-5) #kd6
      ]),
    ('kint_no_cPP_2',
     [Parameter('kint_no_cPP_2kf', 5e-5), #k7
      Parameter('kint_no_cPP_2kr', 1.38e-4) #kd7
      ]),
    ('Akt_path_intern_1Xdimers',
     [Parameter('Akt_path_intern_1Xdimers_kf', 5e-5),
      Parameter('Akt_path_intern_1Xdimers_kr', 5e-4)
      ]),
    ('Akt_path_intern_2Xdimers',
     [Parameter('Akt_path_intern_2Xdimers_kf', 5e-5),
      Parameter('Akt_path_intern_2Xdimers_kr', 5e-4)
      ]),
    ('Akt_path_intern_23_PI3Kdimers',
     [Parameter('Akt_path_intern_23_PI3Kdimers_kf', 5e-5),
      Parameter('Akt_path_intern_23_PI3Kdimers_kr', 5e-4)
      ]),
      #While Jacobian file contains parameter k4b for CPP binding to non-ErbB1 dimers, this is set to 0.
    ('CPP_bind_ErbB1dimers',
     [Parameter('CPP_bind_ErbB1dimerskf', 6.73e-6), #k4
      Parameter('CPP_bind_ErbB1dimerskr', 1.66e-4) #kd4
      ]),
    ('CPPE_bind_ErbB1dimers',
     [Parameter('CPPE_bind_ErbB1dimerskf', 1e-15), #k5 or k5b - Both are set to 0 across all cell types
      Parameter('CPPE_bind_ErbB1dimerskr', 8.0833e-3) #kd5b
      ]),
    ('CPP_int',
     [Parameter('CPP_intkf', 1.667e-8), #k15 Endo --> Cyto
      Parameter('CPP_intkr', 1e-15) #kd15
      ]),
    ('kint_CM_ErbB1ErbB1',
     Parameter('kint_CM_ErbB1ErbB1', 6.73e-5)
     ),
    ('kint_CM_ErbB1ErbBX',
     Parameter('kint_CM_ErbB1ErbBX', 6.73e-6)
     ),
    ('kint_NCM_ErbB1ErbB1',
     Parameter('kint_NCM_ErbB1ErbB1', 6.73e-8)
     ),
    ('kint_NCM_ErbB1ErbBX',
     Parameter('kint_NCM_ErbB1ErbBX', 6.73e-9)
     ),
    ('kint_NCM_ErbB2ErbBX',
     Parameter('kint_NCM_ErbB2ErbBX', 6.73e-7)
     ),
    ('kdeg_erbb1',
     Parameter('kdeg_erbb1', 2.67e-3)
     ),
    ('kdeg_erbb234',
     Parameter('kdeg_erbb234', 4.16e-4)
     ),
    ('kdeg_1',
     Parameter('kdeg_1', .00266742) #k60
     ),
    ('kdeg_2',
     Parameter('kdeg_2', .0471248) #k60b
     ),
    ('kdeg_3',
     Parameter('kdeg_3', 5.2e-4) #k60c
     ),
    ('kdeg_4',
     Parameter('kdeg_4', 5.7e-4) #k61
     ),
    ('kdeg_5',
     Parameter('kdeg_5', 4.16e-4) #k62b
     ),
    ('kdeg_6',
     Parameter('kdeg_6', 4.16e-4) 
     ),
    ('kdeg_7',
     Parameter('kdeg_7', 4.16e-4)
     ),
    ('kdeg_8',
     Parameter('kdeg_8', 4.16e-4)
     ),
    ('krecyc_ErbB1',
     Parameter('krecyc_ErbB1', 1.66e-4)
     ),
    ('krecyc_ErbB234',
     Parameter('krecyc_ErbB234', 1.66e-3)
     ),
     # MAPK pathway rate parameters:
    ('ErbB_bind_GAP_1',
     [Parameter('ErbB_bind_GAP_1kf', 5.91474e-7), #k8
      Parameter('ErbB_bind_GAP_1kr', 2e-1) #kd8
      ]),
    ('ErbB_bind_GAP_2',
     [Parameter('ErbB_bind_GAP_2kf', 9.34641e-6), #k8b
      Parameter('ErbB_bind_GAP_2kr', 2e-2) #kd8b
      ]),
    ('GAP_bind_SHC',
     [Parameter('GAP_bind_SHCkf', 1.39338e-7), #k22
      Parameter('GAP_bind_SHCkr', 1e-1) #kd22 or kd22b (assigned same value in Jacobian file).
      ]),
    ('SHC_phos',
     Parameter('SHC_phoskf', 6) #k23
      ),
    ('SHC_unbound_dephos',
      Parameter('SHC_unbound_dephoskf', 5e-3) #k36
      ),
    ('GRB2_SOS_bind_SHCP_GAP',
     [Parameter('GRB2_SOS_bind_SHCP_GAPkf', 5e-5), #k41
      Parameter('GRB2_SOS_bind_SHCP_GAPkr', 4.29e-2) #kd41
      ]),
    ('SHCP_bind_GRB2SOS',
     [Parameter('SHCP_bind_GRB2SOSkf', 3.5e-5), #k33
      Parameter('SHCP_bind_GRB2SOSkr', 2e-1) #kd33
      ]),
    ('GAP_bind_SHCP_GRB2_SOS',
     [Parameter('GAP_bind_SHCP_GRB2_SOSkf', 4e-7), #k32
      Parameter('GAP_bind_SHCP_GRB2_SOSkr', 1e-1) #kd32
      ]),
    ('GAP_bind_SHCP',
     [Parameter('GAP_bind_SHCPkf', 1.5e-6), #k37
      Parameter('GAP_bind_SHCPkr', 3e-1) #kd37
      ]),
    ('GRB2_bind_SOS',
     [Parameter('GRB2_bind_SOSkf', 7.5e-6), #k35
      Parameter('GRB2_bind_SOSkr', 1.5e-3) #kd35
      ]),
    ('SOS_bind_GAP_SHCP_GRB2',
     [Parameter('SOS_bind_GAP_SHCP_GRB2kf', 1.67e-5), #k25
      Parameter('SOS_bind_GAP_SHCP_GRB2kr', 2.14e-2) #kd25
      ]),
    ('SOS_bind_SHCP_GRB2',
     [Parameter('SOS_bind_SHCP_GRB2kf', 5e-5), #k40
      Parameter('SOS_bind_SHCP_GRB2kr', 6.4e-2) #kd40
      ]),
    ('SOS_bind_GAP_GRB2',
     [Parameter('SOS_bind_GAP_GRB2kf', 1.67e-5), #k17
      Parameter('SOS_bind_GAP_GRB2kr', .06) #kd17
      ]),
    ('RASGDP_bind_bound_GRB2_SOS',
     [Parameter('RASGDP_bind_bound_GRB2_SOSkf', 2.5e-8), #k18
      Parameter('RASGDP_bind_bound_GRB2_SOSkr', 1.3) #kd18
      ]),
    ('RASGTP_bind_bound_GRB2_SOS',
     [Parameter('RASGTP_bind_bound_GRB2_SOSkf', 1.667e-5), #k19
      Parameter('RASGTP_bind_bound_GRB2_SOSkr', 5e-1) #kd19
      ]),
    ('RASGTPact_bind_bound_GRB2_SOS',
     [Parameter('RASGTPact_bind_bound_GRB2_SOSkf', 1.1068e-5), #k20
      Parameter('RASGTPact_bind_bound_GRB2_SOSkr', 4e-1) #kd20
      ]),
    ('RASGTP_unbind_GRB2_SOS',
     [Parameter('RASGTP_unbind_GRB2_SOSkf', 2.3e-5), #kd21
      Parameter('RASGTP_unbind_GRB2_SOSkr', 3.67e-2) #k21
      ]),
    ('Ras_intrinsic_function',
     [Parameter('Ras_intrinsic_function_GTPase', 2.3e-5),
      Parameter('Ras_intrinsic_function_GDPexchange', 3.67e-10)
      ]),
    ('RASGTP_bind_RAF',
     [Parameter('RASGTP_bind_RAFkf', 5e-6), #k28
      Parameter('RASGTP_bind_RAFkr', 5.3e-3) #kd28
      ]),
    ('RASGTP_RAF_cat',
     Parameter('RASGTP_RAF_catkc', 10),
      ),
    ('RAFP_PP1',
     [Parameter('RAFP_PP1kf', 6e-5), #k42
      Parameter('RAFP_PP1kr', .0141589), #kd42
      Parameter('RAFP_PP1kc', 31.6228) #kd43
      ]),
    ('RAFP_MEK',
     [Parameter('RAFP_MEKkf', 1.07e-5), #k44
      Parameter('RAFP_MEKkr', 3.3e-2), #kd52
      Parameter('RAFP_MEKkc', 1.9) #kd45
      ]),
    ('MEKP_PP2',
     [Parameter('MEKP_PP2kf', 4.74801e-8), #k50
      Parameter('MEKP_PP2kr', .252982), #kd50
      Parameter('MEKP_PP2kc', .112387) #kd49
      ]),
    ('RAFP_MEKP',
     [Parameter('RAFP_MEKPkf', 1.07e-5), #k44
      Parameter('RAFP_MEKPkr', 3.3e-2), #kd52
      Parameter('RAFP_MEKPkc', 8e-1) #kd47
      ]),
    ('MEKPP_PP2',
     [Parameter('MEKPP_PP2kf', 2.37e-5), #k48
      Parameter('MEKPP_PP2kr', .79), #kd48
      Parameter('MEKPP_PP2kc', .112387) #kd49
      ]),
    ('MEKPP_ERK',
     [Parameter('MEKPP_ERKkf', 8.85125e-6), #k52
      Parameter('MEKPP_ERKkr', 1.833e-2), #kd44
      Parameter('MEKPP_ERKkc', .28) #kd53
      ]),
    ('ERKP_PP3',
     [Parameter('ERKP_PP3kf', 8.33e-7), #k58
      Parameter('ERKP_PP3kr', 56.7862), #kd58
      Parameter('ERKP_PP3kc', .0076) #kd57
      ]),
    ('MEKPP_ERKP',
     [Parameter('MEKPP_ERKPkf', 8.85125e-6), #k52
      Parameter('MEKPP_ERKPkr', 1.833e-2), #kd44
      Parameter('MEKPP_ERKPkc', 70.1662) #kd55
      ]),
    ('ERKPP_PP3',
     [Parameter('ERKPP_PP3kf', .000397392), #k56
      Parameter('ERKPP_PP3kr', 5), #kd56
      Parameter('ERKPP_PP3kc', .0076) #kd57
      ]),
    ('PP3_deg',
     Parameter('PP3_degkf', .0150356) #k116
     ),
      # AKT pathway event rates:
    ('GRB2_bind_GAP',
     [Parameter('GRB2_bind_GAPkf', 1.67e-5), #k16
      Parameter('GRB2_bind_GAPkr', 5.5e-1) #kd24
      ]),
    ('GRB2_bind_GAP_2',
     [Parameter('GRB2_bind_GAP_2kf', 1.67e-5), #k16
      Parameter('GRB2_bind_GAP_2kr', 2.75e-1) #kd63
      ]),
    ('GRB2_SOS_bind_GAP',
     [Parameter('GRB2_SOS_bind_GAPkf', 7.5e-6), #k34
      Parameter('GRB2_SOS_bind_GAPkr', 3e-2) #kd34
      ]),
    ('GRB2_bind_GAB1',
     [Parameter('GRB2_bind_GAB1kf', 6.67e-5), #k105
      Parameter('GRB2_bind_GAB1kr', 1e-1) #kd105
      ]),
    ('GAB1_bind_ATP',
     [Parameter('GAB1_bind_ATPkf', 1.8704e-8), #k122
      Parameter('GAB1_bind_ATPkr', 1) #kd122
      ]),
    ('GAB1_phos',
     Parameter('GAB1_phoskc', .177828) #kd123
     ),
    ('SHP2_dephos_GAB1P',
     [Parameter('SHP2_dephos_GAB1Pkf', 3.33e-5), #k107
      Parameter('SHP2_dephos_GAB1Pkr', 1e-1), #kd107
      Parameter('SHP2_dephos_GAB1Pkc', 5) #kd108
      ]),
    ('GAB1_bind_PI3K_1',
     [Parameter('GAB1_bind_PI3Kkf', 1.5e-5), #k66
      Parameter('GAB1_bind_PI3Kkr', 2e-1) #kd66
      ]),
    ('GAB1_bind_PI3K_2',
     [Parameter('GAB1_bind_PI3K_2kf', 5e-5), #k67
      Parameter('GAB1_bind_PI3K_2kr', 2e-2) #kd67
      ]),
    ('PIP2_chain_PI3K',
     [Parameter('PIP2_chain_PI3Kkf', 1.33e-5), #k106
      Parameter('PIP2_chain_PI3Kkr', 1e-1) #kd106
      ]),
    ('PIP2_self_catalysis',
     Parameter('PIP2_self_catalysiskc', 2.05e1) #kd68b
     ),
    ('PIP2_bind_PI3K_1',
     [Parameter('PIP2_bind_PI3K_1kf', 2.63418e-8), #k106b
      Parameter('PIP2_bind_PI3K_1kr', 1e-1) #kd106b
      ]),
    ('PIP2_PI3K_catalysis',
     Parameter('PIP2_PI3K_catalysiskc', 2e-1) #kd68
     ),
    ('ErbB23_bind_PI3K_1',
     [Parameter('ErbB23_bind_PI3K_1kf', 1.33e-8), #k106
      Parameter('ErbB23_bind_PI3K_1kr', 1e-1) #kd106
      ]),
    ('ErbB23_bind_PI3K_2',
     [Parameter('ErbB23_bind_PI3K_2kf', 1.33e-8), #k106
      Parameter('ErbB23_bind_PI3K_2kr', 1e-1) #kd106
      ]),
    ('ErbB23_bind_PI3K_3',
     [Parameter('ErbB23_bind_PI3K_3kf', 1.33e-8), #k106
      Parameter('ErbB23_bind_PI3K_3kr', 1e-1) #kd106
      ]),
    ('ErbB23_bind_PI3K_4',
     [Parameter('ErbB23_bind_PI3K_4kf', 1.33e-8), #k106
      Parameter('ErbB23_bind_PI3K_4kr', 1e-1) #kd106
      ]),
    ('ErbB23_bind_PI3K_5',
     [Parameter('ErbB23_bind_PI3K_5kf', 1.33e-8), #k106
      Parameter('ErbB23_bind_PI3K_5kr', 1e-1) #kd106
      ]),
    ('ErbB23_bind_PI3K_6',
     [Parameter('ErbB23_bind_PI3K_6kf', 1.33e-8), #k106
      Parameter('ErbB23_bind_PI3K_6kr', 1e-1) #kd106
      ]),
    ('ErbB23_PI3K_cat_1',
     Parameter('ErbB23_PI3K_cat_1kc', 2e-5),
      ),
    ('ErbB23_PI3K_cat_2',
     Parameter('ErbB23_PI3K_cat_2kc', 2e-5),
      ),
    ('ErbB23_PI3K_cat_3',
     Parameter('ErbB23_PI3K_cat_3kc', 2e-5),
      ),
    ('ErbB23_PI3K_cat_4',
     Parameter('ErbB23_PI3K_cat_4kc', 2e-5),
      ),
    ('ErbB23_PI3K_cat_5',
     Parameter('ErbB23_PI3K_cat_5kc', 2e-5),
      ),
    ('ErbB23_PI3K_cat_6',
     Parameter('ErbB23_PI3K_cat_6kc', 2e-5),
      ),                        
    ('PIP3_bind_AKT',
     [Parameter('PIP3_bind_AKTkf', 3.33e-5), #k69
      Parameter('PIP3_bind_AKTkr', 1e-1) #kd69
      ]),
    ('PIP3_bind_PDK1',
     [Parameter('PIP3_bind_PDK1kf', 1e-15), #k76
      Parameter('PIP3_bind_PDK1kr', 142.262) #kd76
      ]),
    ('AKT_PIP3_bind_PDK1',
     [Parameter('AKT_PIP3_bind_PDK1kf', 6.67e-7), #k70
      Parameter('AKT_PIP3_bind_PDK1kr', 1e-1) #kd70
      ]),
    ('PDK1_AKT_catalysis',
     Parameter('PDK1_AKT_catalysiskc', 2.52e1) #kd71
     ),
    ('PDK1_AKTP_catalysis',
     Parameter('PDK1_AKTP_catalysiskc', 5.01187) #kd72
     ),
    ('AKTP_dephos',
     [Parameter('AKTP_dephoskf', .00374845), #k73
      Parameter('AKTP_dephoskr', 5e-1), #kd73
      Parameter('AKTP_dephoskc', .00633957) #kd75
      ]),
    ('AKTPP_dephos',
     [Parameter('AKTPP_dephoskf', 6.36184e-7), #k74
      Parameter('AKTPP_dephoskr', .355656), #kd74
      Parameter('AKTPP_dephoskc', .00633957) #kd75
      ]),
    ('PIP3_dephos',
     [Parameter('PIP3_dephoskf', 5e-6), #k109
      Parameter('PIP3_dephoskr', 1e-1), #kd109
      Parameter('PIP3_dephoskc', 2e-1) #kd104
      ]),
      # Crosstalk event rates:
    ('ERKPP_phos_GAB1P',
     [Parameter('ERKPP_phos_GAB1Pkf', 3.33e-4), #k110
      Parameter('ERKPP_phos_GAB1Pkr', 1e-1), #kd110
      Parameter('ERKPP_phos_GAB1Pkc', 6.57) #kd111
      ]),
    ('Pase9t_dephos_GAB1PP',
     [Parameter('Pase9t_dephos_GAB1PPkf', 8.33e-8), #k117
      Parameter('Pase9t_dephos_GAB1PPkr', 1e-1), #kd117
      Parameter('Pase9t_dephos_GAB1PPkc', 3e-2) #kd118
      ]),
    ('ERKPP_phos_SOS',
     [Parameter('ERKPP_phos_SOSkf', 1.67e-5), #k64
      Parameter('ERKPP_phos_SOSkr', 3e-1), #kd64
      Parameter('ERKPP_phos_SOSkc', 2e-1) #kd65
      ]),
    ('SOSP_bind_GRB2',
     [Parameter('SOSP_bind_GRB2kf', 8.33e-7), #k101
      Parameter('SOSP_bind_GRB2kr', .03) #kd101
      ]),
    ('AKTPP_phos_RAFP',
     [Parameter('AKTPP_phos_RAFPkf', 4.98816e-6), #k114
      Parameter('AKTPP_phos_RAFPkr', 1e-1), #kd114
      Parameter('AKTPP_phos_RAFPkc', 1) #kd115
      ]),
    ('RASGDP_bind_PI3K',
     [Parameter('RASGDP_bind_PI3Kkf', .0047067), #k112
      Parameter('RASGDP_bind_PI3Kkr', 1e-1), #kd112
      Parameter('RASGDP_bind_PI3Kkc', 177.828) #kd113
      ]),
    ('RASGTP_bind_PI3K',
     [Parameter('RASGTP_bind_PI3Kkf', 1e-8),
      Parameter('RASGTP_bind_PI3Kkr', 1e-1)
      ]),
    ('RAS_PI3K_cat_PIP',
     [Parameter('RAS_PI3K_cat_PIPkf', 2.6e-8),
      Parameter('RAS_PI3K_cat_PIPkr', 1e-1),
      Parameter('RAS_PI3K_cat_PIPkc', 2e-1)
      ]),
    # All parameters below are for mTOR pathway -- used generic Chen Sorger rate constants unless otherwise specified.
    ('mTOR_bind_TORC1ptns',
     [Parameter('mTOR_bind_TORC1ptnskf', 1e-5),
      Parameter('mTOR_bind_TORC1ptnskr', 1e-1)
      ]),
    ('mTOR_bind_TORC2ptns',
     [Parameter('mTOR_bind_TORC2ptnskf', 1e-5),
      Parameter('mTOR_bind_TORC2ptnskr', 1e-1)
      ]),
    ('mTORC2_bind_AKTP',
     [Parameter('mTORC2_bind_AKTPkf', 1e-5),
      Parameter('mTORC2_bind_AKTPkr', 1e-1)
      ]),
    ('mTORC2_cat_AKTPP',
     Parameter('mTORC2_cat_AKTPPkc', 1e-1)
     ),
    ('AKTPP_bind_TSC2',
     [Parameter('AKTPP_bind_TSC2kf', 1e-5),
      Parameter('AKTPP_bind_TSC2kr', 1e-1),
      ]),
    ('AKTPP_phosS939_TSC2',
      Parameter('AKTPP_phosS939_TSC2kc', 1e-1)
      ),
    ('AKTPP_phosS981_TSC2',
      Parameter('AKTPP_phosS981_TSC2kc', 1e-1)
      ),
    ('AKTPP_phosT1462_TSC2',
      Parameter('AKTPP_phosT1462_TSC2kc', 1e-1)
      ),
    ('ERKPP_phos_TSC2',
     [Parameter('ERKPP_phos_TSC2kf', 1e-5),
      Parameter('ERKPP_phos_TSC2kr', 1e-1),
      Parameter('ERKPP_phos_TSC2kc', 1e-1)
      ]),
    ('ERKPP_phos_RSK1',
     [Parameter('ERKPP_phos_RSK1kf', 1e-5),
      Parameter('ERKPP_phos_RSK1kr', 1e-1),
      Parameter('ERKPP_phos_RSK1kc', 1e-1)
      ]),
    ('RSK1_autocat',
     Parameter('RSK1_autocatkc', 1e-1)
     ),
    ('PDK1_phos_RSK1',
     [Parameter('PDK1_phos_RSK1kf', 1e-5),
     Parameter('PDK1_phos_RSK1kr', 1e-1),
     Parameter('PDK1_phos_RSK1kc', 1e-1)
     ]),
    ('RSK1_phos_TSC2',
     [Parameter('RSK1_phos_TSC2kf', 1e-5),
      Parameter('RSK1_phos_TSC2kr', 1e-1),
      Parameter('RSK1_phos_TSC2kc', 1e-1)
      ]),
    ('RSK1_phos_LKB1',
     [Parameter('RSK1_phos_LKB1kf', 1e-5),
      Parameter('RSK1_phos_LKB1kr', 1e-1),
      Parameter('RSK1_phos_LKB1kc', 1e-1)
      ]),
    ('LKB1_phos_AMPK',
     [Parameter('LKB1_phos_AMPKkf', 1e-5),
      Parameter('LKB1_phos_AMPKkr', 1e-1),
      Parameter('LKB1_phos_AMPKkc', 1e-1)
      ]),
    ('AMPK_phos_TSC2',
     [Parameter('AMPK_phos_TSC2kf', 1e-5),
      Parameter('AMPK_phos_TSC2kr', 1e-1),
      Parameter('AMPK_phos_TSC2kc', 1e-1)
      ]),
    ('Rheb_intrinsic_GTPase',
     Parameter('Rheb_intrinsic_GTPasekc', 6.7e-4) #From Marshall et al 2009
     ),
    ('Rheb_GDP_GTP',
     Parameter('Rheb_GDP_GTPkf', 1e10) #Chosen to not be rate limiting
     ),
    ('TSC2_bind_Rheb',
     [Parameter('TSC2_bind_Rhebkf', 1e-5),
      Parameter('TSC2_bind_Rhebkr', 1e-1)
      ]),
    ('TSC2_Rheb_GTPase',
     Parameter('TSC2_Rheb_GTPasekc', .0335) #From Marshall et al 2009
     ),
    ('TSC2p_Rheb_GTPase',
     Parameter('TSC2p_Rheb_GTPasekc', .00335) #Made a factor of 10 lower.
     ),
    ('TSC2pS1387_Rheb_GTPase',
     Parameter('TSC2pS1387_Rheb_GTPasekc', .335) #Made a factor of 10 higher.
     ),
    ('Rheb_bind_mTORC1',
     [Parameter('Rheb_bind_mTORC1kf', 1e-5),
      Parameter('Rheb_bind_mTORC1kr', 1e-1)
      ]),
    ('RhebGTP_phos_mTOR',
     Parameter('RhebGTP_phos_mTORkc', 1e-1)
     ),
    ('FKBP38_bind_mTOR',
     [Parameter('FKBP38_bind_mTORkf', 1e-5),
      Parameter('FKBP38_bind_mTORkr', 1e-1)
      ]),
    ('Rheb_bind_FKBP38',
     [Parameter('Rheb_bind_FKBP38kf', 1e-5),
      Parameter('Rheb_bind_FKBP38kr', 1e-1)
      ]),
    ('mTORC1_phos_S6K',
     [Parameter('mTORC1_phos_S6Kkf', 1e-5),
      Parameter('mTORC1_phos_S6Kkr', 1e-1),
      Parameter('mTORC1_phos_S6Kkc', 1e-1)
      ]),
    ('PDK1_phos_S6K',
     [Parameter('PDK1_phos_S6Kkf', 1e-5),
      Parameter('PDK1_phos_S6Kkr', 1e-1),
      Parameter('PDK1_phos_S6Kkc', 1e-1)
      ]),
    ('S6K_phos_rpS6',
     [Parameter('S6K_phos_rpS6kf', 1e-5),
      Parameter('S6K_phos_rpS6kr', 1e-1),
      Parameter('S6K_phos_rpS6kc', 1e-1)
      ]),
    ('EIF4EBP1_bind_EIF4E',
     [Parameter('EIF4EBP1_bind_EIF4Ekf', 1e-5),
      Parameter('EIF4EBP1_bind_EIF4Ekr', 1e-1)
      ]),
    ('mTORC1_phos_EIF4EBP1',
     [Parameter('mTORC1_bind_EIF4EBP1kf', 1e-5),
      Parameter('mTORC1_bind_EIF4EBP1kr', 1e-1),
      Parameter('mTORC1_bind_EIF4EBP1kc', 1e-1)
      ]),
    ('ERKPP_to_nucleus',
     [Parameter('ERKPP_to_nucleuskf', 1.3e-3),
      Parameter('ERKPP_to_nucleuskr', 5e-5)
      ]),
    ('ERKPP_phos_ELK1',
     [Parameter('ERKPP_phos_ELK1kf', 1e-5),
      Parameter('ERKPP_phos_ELK1kr', 1e-1),
      Parameter('ERKPP_phos_ELK1kc', 1e-1)
    ]),
    ('Bad_cyto_to_mito',
     [Parameter('Bad_cyto_to_mitokf', 1e-1),
      Parameter('Bad_cyto_to_mitokr', 1e-3)
      ]),
    ('Noxa_cyto_to_mito',
     [Parameter('Noxa_cyto_to_mitokf', 1e-1),
      Parameter('Noxa_cyto_to_mitokr', 1e-3)
      ]),
    ('Bim_cyto_to_mito',
     [Parameter('Bim_cyto_to_mitokf', 1e-1),
      Parameter('Bim_cyto_to_mitokr', 1e-3)
      ]),
    ('Puma_cyto_to_mito',
     [Parameter('Puma_cyto_to_mitokf', 1e-1),
      Parameter('Puma_cyto_to_mitokr', 1e-3)
      ]),
    ('Bim_bind_Bcl2',
     [Parameter('Bim_bind_Bcl2kf', 1e-6),
      Parameter('Bim_bind_Bcl2kr', 1e4)
      ]),
    ('Bim_bind_BclXL',
     [Parameter('Bim_bind_BclXLkf', 1e-6),
      Parameter('Bim_bind_BclXLkr', 1e4)
      ]),
    ('Bim_bind_Mcl1',
     [Parameter('Bim_bind_Mcl1kf', 1e-6),
      Parameter('Bim_bind_Mcl1kr', 1e4)
      ]),
    ('Puma_bind_Bcl2',
     [Parameter('Puma_bind_Bcl2kf', 1e-6),
      Parameter('Puma_bind_Bcl2kr', 1e4)
      ]),
    ('Puma_bind_BclXL',
     [Parameter('Puma_bind_BclXLkf', 1e-6),
      Parameter('Puma_bind_BclXLkr', 1e4)
      ]),
    ('Puma_bind_Mcl1',
     [Parameter('Puma_bind_Mcl1kf', 1e-6),
      Parameter('Puma_bind_Mcl1kr', 1e4)
      ]),
    ('Bim_activate_Bax',
     [Parameter('Bim_activate_Baxkf', 1e-7),
      Parameter('Bim_activate_Baxkr', 1e-3),
      Parameter('Bim_activate_Baxkc', 1)
      ]),
    ('FOXO_cyto_to_nucleus',
     [Parameter('FOXO_cyto_to_nucleuskf', 1.3e-3),
      Parameter('FOXO_cyto_to_nucleuskr', 5e-5)
      ]),
    ('Akt_prevent_FOXO_transport',
     [Parameter('Akt_prevent_FOXO_transportkf', 1e-7),
      Parameter('Akt_prevent_FOXO_transportkr', 1e-3),
      Parameter('Akt_prevent_FOXO_transportkc', 1),
      ]),
    ('FOXO_tf_Puma',
     Parameter('FOXO_tf_Puma', 1e-6)
     ),
    ('Erk_phos_Bim',
     [Parameter('Erk_phos_Bimkf', 1e-7),
      Parameter('Erk_phos_Bimkr', 1e-3),
      Parameter('Erk_phos_Bimkc', 1)
      ]),
    ('kdeg_Bim',
     Parameter('kdeg_Bim', 1e-4)
     ),
    ('S6K_phos_Bad',
     [Parameter('S6K_phos_Badkf', 1e-7),
      Parameter('S6K_phos_Badkr', 1e-3),
      Parameter('S6K_phos_Badkc', 1)
      ]),
    ('Akt_phos_Bad',
     [Parameter('Akt_phos_Badkf', 1e-7),
      Parameter('Akt_phos_Badkr', 1e-3),
      Parameter('Akt_phos_Badkc', 1)
      ]),
    #Extra Raf parameters for work with JA.
    ('Raf_Raf_dimerization',
     [Parameter('Raf_Raf_dimerizationkf', 3.3e-6), #Original value = 1e6 /M*s --> Converted assuming cell volume = 2 pL
      Parameter('Raf_Raf_dimerizationkr', .05)]), #Set to BRAF/CRAF/KSR dimerization (kBBf/kBBr) rate from JA modeling work
    ('MEK_phos_CRAF',
     [Parameter('MEK_phos_CRAFkf', 1e-7),
      Parameter('MEK_phos_CRAFkr', 1e-3),
      Parameter('MEK_phos_CRAFkc', 1)
      ])
    ])
