# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 17:06:50 2014

@author: Erin
"""

from pysb import Model, Observable

Model()
import chen_modules
from earm import lopez_modules
from earm import albeck_modules

# Receptor-level events with EGF ligand
chen_modules.rec_monomers()
chen_modules.rec_monomers_lig_EGF()
chen_modules.rec_monomers_scaffold_proteins()
chen_modules.rec_events_lig_EGF()
chen_modules.rec_events_scaffold_protein_binding_grb2()
chen_modules.rec_events_scaffold_protein_binding_shc()
chen_modules.rec_events_scaffold_protein_binding_gab1()
chen_modules.receptor_dimerization()
chen_modules.receptor_phosphorylation()
chen_modules.receptor_dephosphorylation()
chen_modules.receptor_erbb2_lateral_signaling()
chen_modules.receptor_cbl_interactions_erbb1()
chen_modules.receptor_internalization_constitutive()
chen_modules.receptor_internalization_erbb1_clathrin_med()
chen_modules.receptor_internalization_erbb1_clathrin_indepen()
chen_modules.receptor_internalization_erbb234()
chen_modules.receptor_recycling_erbb1()
chen_modules.receptor_recycling_erbb234()
chen_modules.rec_initial()
chen_modules.rec_initial_lig_hEGF()
chen_modules.rec_initial_scaffold_proteins()

# MAPK pathway
chen_modules.mapk_monomers()
chen_modules.mapk_initial()
chen_modules.mapk_events()

# Erk nuclear localization
chen_modules.erk_nuclear_monomers()
chen_modules.erk_nuclear_initial()
chen_modules.erk_nuclear_events()

# AKT pathway
chen_modules.akt_monomers()
chen_modules.akt_initial()
chen_modules.akt_events()

# mTOR signaling
chen_modules.mtor_monomers()
chen_modules.mtor_initial()
chen_modules.mtor_complex_formation()
chen_modules.mtorc1_signaling_monomers()
chen_modules.mtorc1_signaling_initial()
chen_modules.mtorc1_signaling()
chen_modules.mtorc2_signaling()
chen_modules.tsc2_monomers()
chen_modules.tsc2_initial()
chen_modules.tsc2_inhibition_by_akt()
chen_modules.tsc2_inhibition_by_erk()
chen_modules.tsc2_activation_by_erk()
chen_modules.tsc2_gap_function()

# Apoptotic signaling
albeck_modules.ligand_to_c8_monomers()
chen_modules.apoptosis_monomers()
chen_modules.apoptosis_initial()
albeck_modules.apaf1_to_parp_monomers()
lopez_modules.translocate_tBid_Bax_BclxL()
lopez_modules.tBid_activates_Bax_and_Bak()
lopez_modules.effector_auto_activation()
lopez_modules.tBid_binds_all_anti_apoptotics()
lopez_modules.effectors_bind_anti_apoptotics()
lopez_modules.sensitizers_bind_anti_apoptotics()
lopez_modules.lopez_pore_formation()
chen_modules.apoptosis_sensitizer_translocation()
chen_modules.apoptosis_bim_and_puma_bind_anti_apoptotics()
chen_modules.apoptosis_bim_activate_bax()
albeck_modules.pore_to_parp()

# Crosstalk between MAPK and AKT pathways
chen_modules.crosstalk_mapk_akt_monomers()
chen_modules.crosstalk_mapk_akt_initial()
chen_modules.crosstalk_mapk_akt_events()

# Crosstalk between ErbB signaling and apoptotic signaling
chen_modules.crosstalk_erbb_apoptosis_monomers()
chen_modules.crosstalk_erbb_apoptosis_initial()
chen_modules.crosstalk_erbb_apoptosis_events()

# Observables
Observable('obsAKTPP', AKT(bpip3=None, bpdk1=None, S='PP'))
Observable('obsErbB1_P_CE', erbb(ty='1', st='P'))
Observable('obsERKPP', ERK(st='PP'))
Observable('active_mTORC1', mTOR(S2448='P'))
Observable('S6K_PP', S6K(T252='P', T412='P'))
Observable('mBid',  Bid(state='M'))
Observable('aSmac', Smac(state='A'))
Observable('cPARP', PARP(state='C'))
Observable('nuclear_FOXO', FOXO(loc='N'))
Observable('mito_Puma', Puma(state='M'))
Observable('mito_Bad', Bad(state='M'))
Observable('mito_Bim', Bim(state='M'))
