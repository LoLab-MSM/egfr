# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:25:38 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from parameter_dict_A431 import parameter_dict as par

def erk_nuclear_monomers():
    Monomer('ELK1', ['b', 'S383'], {'S383':['U','P']})
         
def erk_nuclear_initial():
    alias_model_components()    
    Initial(ELK1(b=None, S383='U'), ELK1_0)
        
def erk_nuclear_events():
    """Activity of ERK in the nucleus."""
    
    alias_model_components()
    # ERK:PP can translocate to the nucleus and activate ELK-1 by phosphorylating S383 and S389.
    equilibrate(ERK(st='PP', loc='C', b=None), ERK(st='PP', loc='N', b=None), par['ERKPP_to_nucleus'])
    
    catalyze_state(ERK(st='PP', loc='N'), 'b', ELK1(), 'b', 'S383', 'U', 'P', (par['ERKPP_phos_ELK1']))