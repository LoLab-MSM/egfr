# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:27:29 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from .parameter_dict_A431 import parameter_dict as par

def apoptosis_monomers():
    # **Activators**.
    # Bid, states: Untruncated, Truncated, truncated and Mitochondrial
    Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})

    # **Effectors**
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']})

    # Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']})

    # **Anti-Apoptotics**
    Monomer('Bcl2', ['bf', 'state'], {'state': ['C', 'M']})
    Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('Mcl1', ['bf', 'state'], {'state':['M', 'C']})

    # **Sensitizers**
    Monomer('Bad', ['bf', 'state', 'S75', 'S99'], {'state':['C', 'M'], 'S75': ['U', 'P'], 'S99':['U', 'P']}) #S75: phosphorylation site for S6K; S99: major site of AKT phosphorylation
    Monomer('Noxa', ['bf', 'state'], {'state': ['C', 'M']})
    Monomer('Bim', ['bf', 'state', 'S69'], {'state': ['C', 'M'], 'S69':['U', 'P']}) #S69: phosphorylation site for Erk1/Erk2
    Monomer('Puma', ['bf', 'state'], {'state': ['C', 'M']})

    # **Cytochrome C and Smac**
    Monomer('CytoC', ['bf', 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', ['bf', 'state'], {'state':['M', 'C', 'A']})

def apoptosis_initial():
    alias_model_components()
    Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bad(bf=None, state='C', S75='U', S99='U'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None, state='C'), Bcl2_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None, state='M'), Mcl1_0)
    Initial(Noxa(bf=None, state='C'), Noxa_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)
    Initial(Bim(bf=None, state='C', S69='U'), Bim_0)
    Initial(Puma(bf=None, state='C'), Puma_0)
    Initial(C8(bf=None, state='pro'), C8_0)
    
    
def apoptosis_sensitizer_translocation():
    """Translocation of apoptotic sensitizers from cytosol to mitochondria."""
    
    equilibrate(Bad(state='C', S75='U', S99='U', bf=None), Bad(state='M', S75='U', S99='U', bf=None), par['Bad_cyto_to_mito'])
    equilibrate(Noxa(state='C', bf=None), Noxa(state='M', bf=None), par['Noxa_cyto_to_mito'])
    equilibrate(Bim(state='C', S69='U', bf=None), Bim(state='M', S69='U', bf=None), par['Bim_cyto_to_mito'])
    equilibrate(Puma(state='C', bf=None), Puma(state='M', bf=None), par['Puma_cyto_to_mito'])
    
def apoptosis_bim_and_puma_bind_anti_apoptotics():
    """Binding of Bim and Puma to anti-apoptotic proteins Bcl2, Bcl-XL, and Mcl1."""
    
    bind_table([[                        Bcl2,                      BclxL(state='M'),           Mcl1(state='M')],
                [Bim(state='M'),        (par['Bim_bind_Bcl2']),     (par['Bim_bind_BclXL']),    (par['Bim_bind_Mcl1'])],
                [Puma(state='M'),       (par['Puma_bind_Bcl2']),    (par['Puma_bind_BclXL']),   (par['Puma_bind_Mcl1'])]],
               'bf', 'bf')

def apoptosis_bim_activate_bax():
    """Isoforms BimS and Bim-alpha3 can interact with Bax and activate it."""
    
    catalyze_state(Bim(state='M'), 'bf', Bax(s1=None, s2=None), 'bf', 'state', 'M', 'A', (par['Bim_activate_Bax']))