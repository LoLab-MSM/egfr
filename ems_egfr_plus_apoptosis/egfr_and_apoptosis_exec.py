# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 17:06:50 2014

@author: Erin
"""

from pysb import Model, Observable

Model()
from . import receptor_modules
from . import mapk_modules
from . import akt_modules
from . import mtor_modules
from . import nuclear_events_modules
from . import apoptosis_modules
from . import crosstalk_modules
from earm import lopez_modules
from earm import albeck_modules

# Receptor-level events with EGF ligand
receptor_modules.rec_monomers()
receptor_modules.rec_monomers_lig_EGF()
receptor_modules.rec_monomers_inh_ERL()
receptor_modules.rec_monomers_scaffold_proteins()
receptor_modules.rec_events_lig_EGF()
receptor_modules.rec_events_inh_ERL()
receptor_modules.rec_events_scaffold_protein_binding_grb2()
receptor_modules.rec_events_scaffold_protein_binding_shc()
#receptor_modules.rec_events_scaffold_protein_binding_gab1()
receptor_modules.receptor_dimerization()
receptor_modules.receptor_phosphorylation()
receptor_modules.receptor_dephosphorylation()
receptor_modules.receptor_erbb2_lateral_signaling()
#receptor_modules.receptor_cbl_interactions_erbb1()
#receptor_modules.receptor_internalization_constitutive()
#receptor_modules.receptor_internalization_erbb1_clathrin_med()
#receptor_modules.receptor_internalization_erbb1_clathrin_indepen()
#receptor_modules.receptor_internalization_erbb234()
#receptor_modules.receptor_recycling_erbb1()
#receptor_modules.receptor_recycling_erbb234()
receptor_modules.rec_initial()
receptor_modules.rec_initial_lig_hEGF()
receptor_modules.rec_initial_inh_ERL()
receptor_modules.rec_initial_scaffold_proteins()

# MAPK pathway
mapk_modules.mapk_monomers(simplified_raf=True, raf_dimers=False)
mapk_modules.mapk_initial(simplified_raf=True, raf_dimers=False)
mapk_modules.mapk_events(simplified_raf=True, raf_dimers=False)

# Erk nuclear localization
nuclear_events_modules.erk_nuclear_monomers()
nuclear_events_modules.erk_nuclear_initial()
nuclear_events_modules.erk_nuclear_events()

# AKT pathway
akt_modules.akt_monomers()
akt_modules.akt_initial()
akt_modules.akt_events()
#
## mTOR signaling
mtor_modules.mtor_monomers()
mtor_modules.mtor_initial()
mtor_modules.mtor_complex_formation()
mtor_modules.mtorc1_signaling_monomers()
mtor_modules.mtorc1_signaling_initial()
mtor_modules.mtorc1_signaling()
mtor_modules.mtorc2_signaling()
mtor_modules.tsc2_monomers()
mtor_modules.tsc2_initial()
mtor_modules.tsc2_inhibition_by_akt()
mtor_modules.tsc2_inhibition_by_erk()
mtor_modules.tsc2_activation_by_erk()
mtor_modules.tsc2_gap_function()
#
## Apoptotic signaling
albeck_modules.ligand_to_c8_monomers()
apoptosis_modules.apoptosis_monomers()
apoptosis_modules.apoptosis_initial()
albeck_modules.apaf1_to_parp_monomers()
lopez_modules.translocate_tBid_Bax_BclxL()
lopez_modules.tBid_activates_Bax_and_Bak()
lopez_modules.effector_auto_activation()
lopez_modules.tBid_binds_all_anti_apoptotics()
lopez_modules.effectors_bind_anti_apoptotics()
lopez_modules.sensitizers_bind_anti_apoptotics()
lopez_modules.lopez_pore_formation()
apoptosis_modules.apoptosis_sensitizer_translocation()
apoptosis_modules.apoptosis_bim_and_puma_bind_anti_apoptotics()
apoptosis_modules.apoptosis_bim_activate_bax()
albeck_modules.pore_to_parp()

## Crosstalk between MAPK and AKT pathways
crosstalk_modules.crosstalk_mapk_akt_monomers()
crosstalk_modules.crosstalk_mapk_akt_initial()
crosstalk_modules.crosstalk_mapk_akt_events()
#
## Crosstalk between ErbB signaling and apoptotic signaling
crosstalk_modules.crosstalk_erbb_apoptosis_monomers()
crosstalk_modules.crosstalk_erbb_apoptosis_initial()
crosstalk_modules.crosstalk_erbb_apoptosis_events()
#
## Observables
#Observable('obsAKTPP', AKT(bpip3=None, bpdk1=None, S='PP'))
#Observable('obsErbB1_P_CE', erbb(ty='1', st='P'))
#Observable('obsERKPP', ERK(st='PP'))
#Observable('active_mTORC1', mTOR(S2448='P'))
#Observable('S6K_PP', S6K(T252='P', T412='P'))
Observable('mBid',  Bid(state='M'))
Observable('aSmac', Smac(state='A'))
Observable('cPARP', PARP(state='C'))
Observable('obsPARP', PARP(state='U'))
#Observable('nuclear_FOXO', FOXO(loc='N'))
#Observable('mito_Puma', Puma(state='M'))
#Observable('mito_Bad', Bad(state='M'))
#Observable('mito_Bim', Bim(state='M'))
