from pysb import *
from collections import OrderedDict

# Parameters for EGFR Model - BT474 cells
# Obtained from Chen-Sorger Jacobian files unless otherwise specified.
# Commented values prefixed by 'c' are the variables names from Chen Sorger Jacobian files.
parameter_dict = OrderedDict([
    ('initial_amounts', # Initial values for all starting species (in molecules/cell)
     [Parameter('erbb1_0', 3274), #RPPA
      Parameter('erbb2_0', 3.81e5), #RPPA
      Parameter('erbb3_0', 5e3), #Average of other cell types
      Parameter('erbb4_0', 4e2), #Average of other cell types
      Parameter('ATP_0',   1.2e9), #Average of other cell types
      Parameter('DEP_0',   7e6), #Average of other cell types
      Parameter('CPP_0',   8.1e4), #Schoeberl model
      Parameter('GAP_0', 1.2e4), #Schoeberl model
      Parameter('SHC_0', 1.01e6), #Schoeberl model
      # Parameter('SHCPase_0', 1000)
      Parameter('GRB2_0', 4700), #Average of other cell types
      #Parameter('SOS_0', 6.63e4) removed to better represent Chen/Sorger model
      Parameter('RAS_0', 1.14e7), #Schoeberl model
      Parameter('RAF_0', 4e4), #Schoeberl model
      Parameter('MEK_0', 2495), #RPPA
      Parameter('ERK_0', 2.10e7), #Schoeberl model
      Parameter('PP1_0', 4e4), #Schoeberl model
      Parameter('PP2_0', 4e4), #Schoeberl model
      Parameter('PP3_0', 1e7), #Schoeberl model
      Parameter('GRB2_SOS_0', 4.4e7), #This added to better represent Chen Sorger model c30
      Parameter('GAB1_0', 40000.3), #Average of other cell types
      Parameter('PI3K_0', 1132), #RPPA
      Parameter('SHP2_0', 1e6), #Average of other cell types
      Parameter('PIP_0',     7e5), #Average of other cell types
      Parameter('PTEN_0',    21648), #RPPA
      Parameter('SHP_0',     3200), #Average of other cell types
      Parameter('AKT_0',     22224), #RPPA
      Parameter('PDK1_0',     7224), #RPPA
      Parameter('PP2A_III_0', 4.5e6), #Average of other cell types
      Parameter('Pase9t_0', 0), #c521
      Parameter('ERKPP_0', 1152), #RPPA
      Parameter('AKTPP_0', 1536), #RPPA
      Parameter('AKTP_0', 1847), #RPPA
      Parameter('ErbB1P_0', 1130), #RPPA
      Parameter('ErbB2P_0', 2838), #RPPA
    ]),
    # Parameters ('k' prefixed variables are Chen-Sorger variable names from Jacobian files):
    # Receptor-level rate parameters:
    ('LAP_transport_MC',
     [Parameter('LAP_transport_MCkf', 1.3e-3), #from generic internalization rate constant in Chen Sorger paper
      Parameter('LAP_transport_MCkr', 5e-5)
      ]),
    ('LAP_bind_ErbB1',
     [Parameter('LAP_bind_ErbB1kf', 1e7),
      Parameter('LAP_bind_ErbB1kr', 3e-2) #Calculated from Ki estimate ~3 nM/L given in Konecny et al Cancer Research 2006
      ]),
    ('LAP_bind_ErbB2',
     [Parameter('LAP_bind_ErbB2kf', 1e7),
      Parameter('LAP_bind_ErbB2kr', .13) #Calculated from Ki estimate ~13 nM/L given in Konecny et al Cancer Research 2006
      ]),
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
    ('ErbB21_bind_ATP',
     [Parameter('ErbB21_bind_ATPkf', 1.8704e-8),
      Parameter('ErbB21_bind_ATPkr', 1),
      ]),
    ('ErbB22_bind_ATP',
     [Parameter('ErbB2_bind_ATPkf', 1.8704e-8),
      Parameter('ErbB2_bind_ATPkr', 1) #kd122
      ]),
    ('ErbB23_bind_ATP',
     [Parameter('ErbB3_bind_ATPkf', 1.8704e-8), #k122
      Parameter('ErbB3_bind_ATPkr', 1) #kd122
      ]),
    ('ErbB24_bind_ATP',
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
    ('kint_no_cPP_1', 
     [Parameter('kint_no_cPP_1kf', .013), #k6
      Parameter('kint_no_cPP_1kr', 5e-5) #kd6
      ]),
    ('kint_no_cPP_2',
     [Parameter('kint_no_cPP_2kf', 5e-5), #k7
      Parameter('kint_no_cPP_2kr', 1.38e-4) #kd7
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
    ('LAP_intern',
     [Parameter('LAP_internkf', 1.3e-3), #From generic internalization constant in Chen Sorger 2009
      Parameter('LAP_internkr', 5e-5)
      ]),
    ('LAP_ErbB1_intern',
     [Parameter('LAP_ErbB1_internkf', 1.3e-3),
      Parameter('LAP_ErbB1_internkr', 5e-5)
      ]),
    ('LAP_ErbB1d_intern',
     [Parameter('LAP_ErbB1d_internkf', 1.3e-3),
      Parameter('LAP_ErbB1d_internkr', 5e-5)
      ]),
    ('LAP_ErbB2_intern',
     [Parameter('LAP_ErbB2_internkf', 1.3e-3),
      Parameter('LAP_ErbB2_internkr', 5e-5)
      ]),
    ('LAP_ErbB2d_intern',
     [Parameter('LAP_ErbB2d_internkf', 1.3e-3),
      Parameter('LAP_ErbB2d_internkr', 5e-5)
      ]),
    ('LAP2_ErbB_intern',
     [Parameter('LAP2_ErbB_internkf', 1.3e-3),
      Parameter('LAP2_ErbB_internkr', 5e-5)
      ]),
    ('LAP_deg',
     Parameter('LAP_deg', 5.7e-4) #Don't know what this should be -- starting from EGF degradation rate
      ),
    ('LAP_ErbB1_deg',
     Parameter('LAP_ErbB1_deg', .00266742)
     ),
    ('LAP_ErbB1d_deg',
     Parameter('LAP_ErbB1d_deg', 4.16e-4)
     ),
    ('LAP_ErbB2_deg',
     Parameter('LAP_ErbB2_deg', .0471248)
     ),
    ('LAP_ErbB2d_deg',
     Parameter('LAP_ErbB2d_deg', 4.16e-4)
     ),
    ('LAP2_ErbB_deg',
     Parameter('LAP2_ErbB_deg', 4.16e-4)
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
     [Parameter('SHC_phoskf', 6), #k23
      Parameter('SHC_phoskr', 6e-2) #kd23
      ]),
    ('SHC_unbound_phos',
     [Parameter('SHC_unbound_phoskf', 1e-15), #kd36 - Chen-Sorger model written P --> U
      Parameter('SHC_unbound_phoskr', 5e-3) #k36
      ]),
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
     [Parameter('RASGDP_bind_bound_GRB2_SOSkf', 2.5e-5), #k18
      Parameter('RASGDP_bind_bound_GRB2_SOSkr', 1.3) #kd18
      ]),
    ('RASGTP_bind_bound_GRB2_SOS',
     [Parameter('RASGTP_bind_bound_GRB2_SOSkf', 1.667e-7), #k19
      Parameter('RASGTP_bind_bound_GRB2_SOSkr', 5e-1) #kd19
      ]),
    ('RASGTPact_bind_bound_GRB2_SOS',
     [Parameter('RASGTPact_bind_bound_GRB2_SOSkf', 1.1068e-5), #k20
      Parameter('RASGTPact_bind_bound_GRB2_SOSkr', 4e-1) #kd20
      ]),
    ('RASGTP_unbind_GRB2_SOS',
     [Parameter('RASGTP_unbind_GRB2_SOSkf', 2.3e-1), #kd21
      Parameter('RASGTP_unbind_GRB2_SOSkr', 3.67e-7) #k21
      ]),
    ('RASGTP_bind_RAF',
     [Parameter('RASGTP_bind_RAFkf', 5e-6), #k28
      Parameter('RASGTP_bind_RAFkr', 5.3e-3) #kd28
      ]),
    ('RASGTP_RAF_cat',
     [Parameter('RASGTP_RAF_catkf', 1.17e-6), #k29
      Parameter('RASGTP_RAF_catkr', 3.1) #kd29
      ]),
    ('RAFP_PP1',
     [Parameter('RAFP_PP1kf', 6e-5), #k42
      Parameter('RAFP_PP1kr', .0141589), #kd42
      Parameter('RAFP_PP1kc', 31.6228) #kd43
      ]),
    ('MEK_syn',
     Parameter('MEK_syn', .14)
     ),
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
    ('PI3K_syn',
     Parameter('PI3K_syn', .03)
     ),
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
    ('ErbB23_bind_PI3K',
     [Parameter('ErbB23_bind_PI3Kkf', 5e-5), #k67
      Parameter('ErbB23_bind_PI3Kkr', 2e-2) #kd67
      ]),
    ('ErbB23_PI3K_bind_PIP2',
     [Parameter('ErbB23_PI3K_bind_PIP2kf', 1.33e-5), #k106
      Parameter('ErbB23_PI3K_bind_PIP2kr', 1e-1) #kd106
      ]),
    ('ErbB23_PI3K_cat_PIP3',
     Parameter('ErbB23_PI3K_cat_PIP3kc', 123)
     ),
    ('AKT_syn',
     Parameter('AKT_syn', 6.6)
     ),
    ('PDK1_syn',
     Parameter('PDK1_syn', 2)
     ),
    ('PTEN_syn',
     Parameter('PTEN_syn', .31)
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
      ])
    ])
