# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:29:46 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from parameter_dict_A431 import parameter_dict as par

def crosstalk_mapk_akt_monomers():
    Monomer('Pase9t', ['bgab1'])
    
def crosstalk_mapk_akt_initial():
    alias_model_components()
    Initial(Pase9t(bgab1=None), Pase9t_0)

def crosstalk_mapk_akt_events():
    """Defines crosstalk events between MAPK and AKT pathways."""
    
    alias_model_components()
    #ERK:P:P phosphorylates ErbBdimer-GRB2-GAB1:P (making it unable to bind PI3K)
    catalyze_state(ERK(st='PP', loc='C'), 'b', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None), 'bERKPP', 'S', 'P', 'PP', (par['ERKPP_phos_GAB1P']))

    #ErbBdimer-GRB2-GAB1:P:P is dephosphorylated by Pase9t
    catalyze_state(Pase9t(), 'bgab1', GAB1(bgrb2=ANY), 'bPase9t', 'S', 'PP', 'P', (par['Pase9t_dephos_GAB1PP']))

    #AKT:PP phosphorylates RAF:P at Ser295, preventing MEK phosphorylation.
    catalyze_state(AKT(S='PP', bpip3=None), 'bpdk1', RAF(st='P'), 'b', 'ser259', 'U', 'P', (par['AKTPP_phos_RAFP']))

    #RAS-GTP binds PI3K and activates PI3K catalytic function.
    bind(RAS(bsos=None, braf=None, st='GTP'), 'bpi3k', PI3K(bgab1=None, bpip=None, berb=None), 'bras', par['RASGTP_bind_PI3K'])
    
    catalyze_state(PI3K(bras=ANY, bgab1=None, berb=None), 'bpip', PIP(bpdk1=None, bakt=None), 'bpi3k', 'S', 'PIP2', 'PIP3', par['RAS_PI3K_cat_PIP'])
    
def crosstalk_erbb_apoptosis_monomers():
    Monomer('FOXO', ['active', 'loc', 'b'], {'active':['Y', 'N'], 'loc':['C', 'N']}) #FOXO 1/3 transcription factors

def crosstalk_erbb_apoptosis_initial():
    alias_model_components()
    Initial(FOXO(active='Y', loc='C', b=None), FOXO_0)

def crosstalk_erbb_apoptosis_events(akt_puma=True, erk_bim=True, s6k_bad=True, akt_bad=True):
    """Crosstalk interactions between pathways downstream of ErbB signaling and apoptotic signaling."""
    
    alias_model_components()
    if akt_puma:
        Rule('FOXO_cyto_to_nucleus',
             FOXO(active='Y', loc='C', b=None) <>
             FOXO(active='Y', loc='N', b=None),
             *par['FOXO_cyto_to_nucleus'])
        
        catalyze_state(AKT(S='PP', bpdk1=None), 'bpip3', FOXO(loc='C'), 'b', 'active', 'Y', 'N', par['Akt_prevent_FOXO_transport'])
        
        #This is almost certainly not the best way to handle this but it's somewhere to start
        Rule('FOXO_tf',
             FOXO(loc='N') >>
             FOXO(loc='N') + Puma(bf=None, state='C'),
             par['FOXO_tf_Puma'])
        
    if erk_bim:
        catalyze_state(ERK(st='PP', loc='C'), 'b', Bim(state='C'), 'bf', 'S69', 'U', 'P', par['Erk_phos_Bim'])
        
        degrade(Bim(state='C', bf=None, S69='P'), par['kdeg_Bim'])
        
    if s6k_bad:
        catalyze_state(S6K(T252='P', T412='P'), 'b', Bad(state='C'), 'bf', 'S75', 'U', 'P', par['S6K_phos_Bad'])
    
    if akt_bad:
        catalyze_state(AKT(S='PP', bpdk1=None), 'bpip3', Bad(state='C'), 'bf', 'S99', 'U', 'P', par['Akt_phos_Bad'])