# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:15:05 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from parameter_dict_A431 import parameter_dict as par

#Calculated molecules for a given extracellular concentration
#These assume a cellular volume = 2 pL
#Cell radius then = .78 pm
#Volume of cube with same radius = (2r)^3 = 3.8 pL
#Which gives an approximate extracellular volume of 3.8 - 2 = 1.8 pL

#Based on above assumptions, molecules/cell for various extracellular ligand/inhibitor concentrations = CAV where:
# conc = C
# Avogadro's number = A
# Volume of extracellular from above = V

# .010 nM = 11 molecules
# 1.00 nM = 1080 molecules
# 5.00 nM = 5420 molecules
# 1 microM = 1.08 * 10^6 molecules
# 3 microM = 3.25 * 10^6 molecules
        
# Receptor level interactions.  
# This includes ligand binding, binding of drugs that target receptors, receptor dimerization events, and receptor internalization and recycling.

def rec_monomers():

    """ Declares the ErbB receptor interaction monomers (except ligands).
    """
    Monomer('erbb', ['bl', 'bd', 'b', 'Y1045', 'ty', 'st', 'loc', 'pi3k1', 'pi3k2', 'pi3k3', 'pi3k4', 'pi3k5', 'pi3k6'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E']}) 
    # ErbB receptors, types 1-4    
    # Sites: bl: lig, bd: dimer, 
    # b: binding site for ATP during phosphorylation and for scaffolding proteins, 
    # Y1045: binding site for the Cbl, a protein which catalyzes ErbB1 ubiquitinylation.  Binding of Cbl to this site is required for ErbB1 degradation (though not for internalization).
    # ty: rec type, st: (U)n(P)hosphorylated, 
    # loc: (C)yto 'brane or (E)ndosome 'brane, 
    # pi3k 1-6: PI3K binding sites (these 6 sites occur on ErbB3)      
    Monomer('DEP', ['b']) # A generic phosphorylase for dephosphorylation of receptors
    Monomer('ATP', ['b'])
    Monomer('ADP') # ATP and ADP
    Monomer('CBL', ['b']) # An E3 ubiquitin ligase responsible for ErbB1 ubiquitination (Roepstorff 2008)

def rec_monomers_lig_EGF():
    """Declares the monomer for the ligand EGF."""
    Monomer('EGF', ['b', 'st'], {'st':['M', 'E']}) # Epidermal Growth Factor ligand

def rec_monomers_lig_HRG():
    """Declares the monomer for the ligand heregulin."""
    Monomer('HRG', ['b', 'st'], {'st': ['M', 'E']}) # Heregulin ligand

def rec_monomers_inh_ERL():
    """Declares the monomer for the inhibitor erlotinib."""
    Monomer('ERL', ['b'])
    
def rec_monomers_scaffold_proteins():
    """Declares the monomers for scaffolding proteins that may be shared between pathways."""
    Monomer('SHC', ['bgap', 'bgrb', 'batp', 'st'], {'st':['U','P']})
    Monomer('GRB2', ['b', 'bsos', 'bgap', 'bgab1', 'bcbl'])
    Monomer('GAB1', ['bgrb2', 'bshp2', 'bpi3k', 'batp','bERKPP','bPase9t','S'],{'S':['U','P','PP']})

def rec_initial_lig_hEGF():
    """Declares the initial conditions for high EGF (5 nM)."""
    Parameter('EGF_0',     5420) # c1 5 nM EGF = 5420 molec/cell; .01 nM EGF = 11 molec/cell
    alias_model_components()
    Initial(EGF(b=None, st='M'), EGF_0)
    
def rec_initial_lig_lEGF():
    """Declares the initial conditions for low EGF (.01 nM)."""
    Parameter('EGF_0',      11) # c1 5 nM EGF = 5420 molec/cell; .01 nM EGF = 11 molec/cell
    alias_model_components()
    Initial(EGF(b=None, st='M'), EGF_0)

def rec_initial_lig_hHRG():
    """Declares the initial conditions for high heregulin (5 nM)."""
    Parameter('HRG_0',         5420) # c514 5 nM HRG = 5420 molec/cell; .01 nM EGF = 11 molec/cell
    alias_model_components()
    Initial(HRG(b=None, st='M'), HRG_0)

def rec_initial_lig_lHRG():
    """Declares the initial conditions for low heregulin (.01 nM)."""
    Parameter('HRG_0',         11) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    Initial(HRG(b=None, st='M'), HRG_0)

def rec_initial_inh_ERL():
    """Declares the initial conditions for erlotinib inhibition (3 microM)."""
    Parameter('ERL_0', 3.25e6)
    alias_model_components()
    Initial(ERL(b=None), ERL_0)
    
def rec_initial():
    """Declares the initial conditions for all monomers in receptor interactions except ligands."""
    # # Initial concentrations (except DEP1) for all cell types taken from Chen et al 2009 -- see Jacobian files
    alias_model_components()
    Initial(erbb(bl=None, bd=None, b=None, ty='1', st='U', loc='C', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), par['erbb1_0'])
    Initial(erbb(bl=None, bd=None, b=None, ty='2', st='U', loc='C', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), erbb2_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='3', st='U', loc='C', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), erbb3_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='4', st='U', loc='C', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), erbb4_0)
    Initial(ATP(b=None), ATP_0)
    Initial(DEP(b=None), DEP_0)
    Initial(CBL(b=None), CBL_0)
    
def rec_initial_scaffold_proteins():
    """Declares initial conditions for scaffolding proteins that may be shared between pathways."""
    alias_model_components()
    Initial(SHC(bgap=None, bgrb=None, batp=None, st='U'), SHC_0)
    Initial(GRB2(b=None, bsos=None, bgap=None, bgab1=None, bcbl=None), GRB2_0)
    Initial(GAB1(bgrb2=None, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'), GAB1_0)

def rec_events_lig_EGF():
    """ Receptor events involving the ligand EGF."""
    
    alias_model_components()
    
    # Binding of EGF to undimerized ErbB1
    bind_table([[                                                                                                                           EGF(st='M')],
                [erbb(ty='1', bl=None, bd=None, b=None, st='U', loc='C'),                                                                   (par['EGF_bind_ErbB1'])],
                [erbb(ty='3', b=None, bd=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),    None],
                [erbb(ty='4', b=None, bd=None, st='U', loc='C'),                                                                            None]],
                'bl', 'b')
                
    # Binding of EGF to dimerized ErbB1
    
    bind_complex(erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C') % erbb(bl=None, bd=1, b=None, st='U', loc='C'), 'bl', EGF(st='M', b=None), 'b', par['EGF_bind_ErbB1d'], m1=erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C'))
    
    bind_complex(erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C') % erbb(bl=ANY, bd=1, b=None, st='U', loc='C'), 'bl', EGF(st='M', b=None), 'b', par['EGF_bind_ErbB1d'], m1=erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C'))
    
    # EGF binding/unbinding from endosomal receptors
    
    Rule('EGFE_bind_ErbBE',
         erbb(ty='1', bl=None, loc='E') + EGF(st='E', b=None) <>
         erbb(ty='1', bl=1, loc='E') % EGF(st='M', b=1),
         *par['EGFE_bind_ErbBE'])
    
    # degradation of endosomal EGF
    degrade(EGF(b=None, st='E'), par['kdeg_4'])
    
def rec_events_lig_HRG():
    """ Receptor events involving the ligand heregulin."""
    alias_model_components()
    # Binding of HRG to undimerized ErbB3 and ErbB4
    bind_table([[                                                                                                                           HRG(st='M')],
                [erbb(ty='1', bl=None, bd=None, b=None, st='U', loc='C'),                                                                   None],
                [erbb(ty='3', b=None, bd=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),    (par['HRG_bind_ErbB3'])],
                [erbb(ty='4', b=None, bd=None, st='U', loc='C'),                                                                            (par['HRG_bind_ErbB4'])]],
                'bl', 'b')
                
    # Binding of HRG to dimerized ErbB3 and ErbB4
    for i in ['3', '4']:
        bind_complex(erbb(ty=i, bl=None, bd=1, b=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(bl=None, bd=1, b=None, st='U', loc='C'), 'bl', HRG(st='M', b=None), 'b', par['HRG_bind_ErbB'+i], m1=erbb(ty=i, bl=None, bd=1, b=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    #HRG binding/unbinding from endosomal receptors

    for i in ['3', '4']:
        Rule('HRGE_bind_ErbBE'+i,
             erbb(ty=i, bl=None, loc='E') + HRG(st='E', b=None) <>
             erbb(ty=i, bl=1, loc='E') % HRG(st='M', b=1),
             *par['HRGE_bind_ErbBE'])
    
    # degradation of endosomal HRG     
    degrade(HRG(b=None, st='E'), par['kdeg_HRG'])

def rec_events_inh_ERL():
    """Receptor events involving the EGFR kinase inhibitor erlotinib.  Binds in the ATP binding pocket."""
    alias_model_components()
    
    #Binding of erlotinib to EGFR
    #Assumption here: erlotinib only binds to dimers and its ligand binding status doesn't matter (need to check this)
    bind(erbb(ty='1', bd=ANY, st='U', loc='C'), 'b', ERL(), 'b', par['EGFR_bind_ERL'])    

def rec_events_scaffold_protein_binding_shc():
    """Binding of scaffold proteins (which are often used in multiple 'pathways') to phosphorylated dimers.  This module covers events with the scaffold protein Shc that aren't pathway specific."""
    
    # SHC binds to ErbB dimers
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),'b', SHC(batp=None, st='U', bgrb=None), 'bgap', par['GAP_bind_SHC'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))
        
    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),'b', SHC(batp=None, st='U', bgrb=None), 'bgap', par['GAP_bind_SHC'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    # Bound and unbound SHC phosphorylation:
    Rule('SHC_phos',
         erbb(bd=ANY, st='P', b=1) % SHC(bgap=1, bgrb=None, batp=None, st='U') >>
         erbb(bd=ANY, st='P', b=1) % SHC(bgap=1, bgrb=None, batp=None, st='P'),
         par['SHC_phos'])

    #SHC:P binds/unbinds ErbB dimers
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=None), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=None), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    #SHC:P-GRB2 binds ErbB dimers
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcbl=None, b=2), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcbl=None, b=2), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))

    # Unbound SHC dephosphorylation  
    Rule('SHC_unbound_dephos',
         SHC(bgap=None, bgrb=None, batp=None, st='P') >>
         SHC(bgap=None, bgrb=None, batp=None, st='U'),
         par['SHC_unbound_dephos'])
    
def rec_events_scaffold_protein_binding_grb2():    
    """Binding of scaffold proteins (which are often used in multiple 'pathways') to phosphorylated dimers.  This module covers events with the scaffold protein Grb2 that aren't pathway specific."""
    
    # GRB2 binds to ErbBdimer-SHC:P without SOS:
    bind(SHC(batp=None, st='P', bgap=ANY), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=None, bcbl=None), 'b', par['GRB2_bind_GAP'])
    
    #GRB2 binds ErbBdimer-complex (without requiring SHC bound to complex):
    #Bind GRB2 without SOS already bound:
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bgap=None), 'bgap', par['GRB2_bind_GAP_2'], m1=erbb(ty='1', bd=1, st='P', b=None))
    
    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcbl=None, bgap=None), 'bgap', par['GRB2_bind_GAP_2'], m1=erbb(ty='2', bd=1, st='P', b=None))

def rec_events_scaffold_protein_binding_gab1():
    """Binding of scaffold proteins (which are often used in multiple 'pathways') to phosphorylated dimers.  This module covers events with the scaffold protein Gab1 that aren't pathway specific."""     
     
    #GAB1 binds ErbB:ErbB-GRB2. 
    bind_complex(erbb(bd=1) % erbb(bd=1) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None, bcbl=None), 'bgab1', GAB1(bgrb2=None, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'), 'bgrb2', par['GRB2_bind_GAB1'])
       
     
def receptor_dimerization():
    """ ErbB receptor dimerization. 
    """
    
    # ErbB dimerization
    # ErbB1 is not required to contain a ligand in order to dimerize (3 and 4 are)
    # Rates for ErbB1 dimerization with and without ligand are different
    # Assumptions: 
    # ErbB3/ErbB3 and ErbB3/ErbB4 dimers neglected (assumed to be at very low concentration since both components likely are)
    # ErbB3 and ErbB4 have the same dimerization rate independent of ligand
    # If we ever add multiple ligand types per receptor, the rates for dimerization as defined currently will not depend on the ligand type
    # These rules only cover binding/unbinding of unphosphorylated receptors.  Binding/unbinding of phosphorylated receptors is currently in ErbB2 lateral signaling section 
    
    alias_model_components()    
    
    erbb1 = erbb(ty='1', bl=None, b=None, st='U', loc='C')
    erbb1Lig = erbb(ty='1', bl=ANY, b=None, st='U', loc='C')
    erbb2Lig = erbb(ty='2', b=None, st='U', loc='C')
    erbb3Lig = erbb(ty='3', b=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None)
    erbb4Lig = erbb(ty='4', b=None, st='U', loc='C')
    bind_table([[                          erbb1,                      erbb1Lig,                     erbb2Lig,                    erbb3Lig, erbb4Lig],
                [erbb1,                    (par['ErbB1_bind_ErbB1']),  None,                         None,                        None,     None],
                [erbb1Lig,                 (par['ErbB1_bind_ErbB1L']), (par['ErbB1L_bind_ErbB1L']),  None,                        None,     None],
                [erbb2Lig,                 (par['ErbB1_bind_ErbB2']),  (par['ErbB1L_bind_ErbB2']),   (par['ErbB2_bind_ErbB2']),   None,     None],
                [erbb3Lig,                 (par['ErbB1_bind_ErbB3']),  (par['ErbB1L_bind_ErbB3']),   (par['ErbB2_bind_ErbB3']),   None,     None],
                [erbb4Lig,                 (par['ErbB1_bind_ErbB4']),  (par['ErbB1L_bind_ErbB4']),   (par['ErbB2_bind_ErbB4']),   None,     None]],
                'bd', 'bd')

    #alias_model_components()

def receptor_phosphorylation():
    """ ErbB receptor phosphorylation.
    """
    
    # ATP binding: 
    # Assumption: ATP only binds to dimers
    # Assumption: ATP only binds to dimers that contain at least one ligand.
    # ErbB1, ErbB2, and ErbB4 have kinase domains and can bind ATP.
    
    bind_complex(erbb(ty='1', st='U', loc='C', bl=None, b=None, bd=1) % erbb(st='U', loc='C', bl=ANY, b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', ATP(b=None), 'b', par['ErbB1_bind_ATP'], m1=erbb(ty='1', st='U', loc='C', bl=None, b=None, bd=1))
    
    bind_complex(erbb(ty='2', st='U', loc='C', b=None, bd=1) % erbb(st='U', loc='C', bl=ANY, b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', ATP(b=None), 'b', par['ErbB2_bind_ATP'], m1=erbb(ty='2', st='U', loc='C', b=None, bd=1))
    
    bind_complex(erbb(ty='4', st='U', loc='C', bl=ANY, b=None, bd=1) % erbb(st='U', loc='C', b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', ATP(b=None), 'b', par['ErbB4_bind_ATP'], m1=erbb(ty='4', st='U', loc='C', bl=ANY, b=None, bd=1))   

    bind_complex(erbb(ty='4', st='U', loc='C', bl=None, b=None, bd=1) % erbb(st='U', loc='C', bl=ANY, b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', ATP(b=None), 'b', par['ErbB4_bind_ATP'], m1=erbb(ty='4', st='U', loc='C', bl=None, b=None, bd=1))

    # Cross phosphorylation: only ErbB1, 2, and 4 have ATP, and they can cross-phosphorylate any other receptor
    # ErbB2:ErbB2 pairs only happen by dissociation of phosphorylated monomers

    # Assumption: Both dimers become phosphorylated/dephosphorylated concurrently (unrealistic)

    for i in ['1','2','4']:
        for j in ['1','2','3','4']:
            Rule("cross_phospho_"+i+"_"+j,
                 ATP(b=1) % erbb(ty=i, b=1,    bd=2, st='U') % erbb(ty=j, bd=2, b=None, st='U', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) >>
                 ADP()    + erbb(ty=i, b=None, bd=2, st='P') % erbb(ty=j, bd=2, b=None, st='P', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
                 Parameter("kcp"+i+j, par['ATP_phos_ErbB']))

def receptor_dephosphorylation():
    """ErbB receptor dephosphorylation.
    """
    # Receptor Dephosphorylation
    # DEPHOSPHORYLATION: 
    #  * Density enhanced phosphatase1 (DEP1) dephosphorylates ERB1 (at the cell-membrane)
    #  * Protein Tyrosine Phosphatase1b (PTP1b) dephosphorylates all RTKs (at the endo-membrane) 
    #  Berset, TA, Hoier, EF, Hajnal, A: Genes Dev. 19:1328-1340 (2005)
    #  Haj, FG, Verver, PJ, Squire, A, Neel, BG, Bastiaens, PI: Science 295:1708-1711 (2002)
    
    # Assumption: Both dimers are dephosphorylated concurrently.
    # Assumption: Bound ligand is not necessary for phosphatase to bind.
    # Assumption: Phosphatase binds at kinase domain (so can only bind to ErbB1, ErbB2, and ErbB4).

    for i in ['1', '2', '4']:
        bind_complex(erbb(ty=i, st='P', loc='C', b=None, Y1045=None, bd=1) % erbb(st='P', loc='C', b=None, Y1045=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', DEP(b=None), 'b', par['ErbBP'+i+'_bind_DEP'], m1=erbb(ty=i, st='P', loc='C', Y1045=None, b=None, bd=1))

    for i in ['1','2','4']:
        for j in ['1','2','3','4']:
            Rule("cross_DEphospho_"+i+"_"+j,
                 DEP(b=1)   %  erbb(ty=i, b=1,    bd=2, st='P', Y1045=None) % erbb(ty=j, bd=2, b=None, st='P', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) >>
                 DEP(b=None) + erbb(ty=i, b=None, bd=2, st='U', Y1045=None) % erbb(ty=j, bd=2, b=None, st='U', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
                 Parameter("kcd"+i+j, par['DEP_dephos_ErbB']))

def receptor_erbb2_lateral_signaling():
    """Unbinding of phosphorylated ErbB dimers.  This allows phosphorylated ErbB2 to unbind from ligand-bound, phosphorylated ErbB1 aftr an EGF signal.  
    The monomeric phosphorylated ErbB2 can then phosphorylate ErbB1, ErbB3 or ErbB4, allowing 'lateral' transmission of the EGF signal through ErbB2/ErbB3 and ErbB2/ErbB4 dimers.
    """
    
    #ErbB2 lateral signaling - ErbB2P-ErbB2P dimers can only form by the dissociation of ligand-containing, phosphorylated dimers containing ErbB2.  
    # The monomeric activated ErbB2 can then bind and activate other monomers (ErbB1, 3, or 4) 
    # This allows an EGF signal to be transmitted by ErbB2/ErbB3 and ErbB2/ErbB4 complexes, even though 3 and 4 can't bind EGF. 
    # This also allows an HRG signal to be transmitted by ErbB1/ErbB2 dimers, even though ErbB1 can't bind HRG.
    # It also allows the formation of active ErbB2/ErbB2 dimers (that still require an EGF or HRG signal to initialize signaling).
    
    # Rules for binding/unbinding of phosphorylated receptors.

    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C', Y1045=None), 'bd', erbb(ty='1', bd=None, st='P', b=None, loc='C', Y1045=None), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C', Y1045=None), 'bd', erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C', Y1045=None), 'bd', erbb(ty='3', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C', Y1045=None), 'bd', erbb(ty='4', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])

    bind(erbb(bl=None, ty='2', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'bd', erbb(bl=None, ty='2', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'bd', par['ErbB2P_ErbBXP_bind'])
    bind(erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='3', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'bd', par['ErbB2P_ErbBXP_bind'])
    bind(erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='4', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB2P_ErbBXP_bind'])
    
    # Phosphorylation of ErbB1, ErbB3 or ErbB4 by a phosphorylated ErbB2 receptor to yield phosphorylated ErbB1/ErbB2, ErbB2/ErbB3 or ErbB2/ErbB4 receptor dimers.  
    
    for i in ['1', '3', '4']:
        Rule('ErbB2_lateralsignal_'+i,
             erbb(ty='2', bd=None, st='P', b=None, loc='C') + erbb(ty=i, bd=None, st='U', b=None, loc='C', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) >>
             erbb(ty='2', bd=1, st='P', b=None, loc='C') % erbb(ty=i, bd=1, st='P', b=None, loc='C', Y1045=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
             Parameter('ErbB2_lateralsignal_k'+i, par['ErbB2_lateralsignal']))

def receptor_cbl_interactions_erbb1():
    """Interactions of ErbB1 with the protein Cbl.  Cbl can either bind directly to Y1045P of ErbB1, or indirectly to the adaptor protein Grb2.
    The first interaction is required for ErbB1 degradation and the second for ErbB1 internalization (Roepstorff 2008). """

    # Direct interaction via Y1045P on ErbB1
    bind_complex(erbb(ty='1', st='P', loc='C') % erbb(st='P', loc='C', Y1045=None), 'Y1045', CBL(), 'b', par['ErbB1_Cbl_Y1045_Cyto'], m1=erbb(ty='1', st='P', loc='C'))
    bind_complex(erbb(ty='1', st='P', loc='E') % erbb(st='P', loc='E', Y1045=None), 'Y1045', CBL(), 'b', par['ErbB1_Cbl_Y1045_Endo'], m1=erbb(ty='1', st='P', loc='E'))
    
    # Indirect interaction via Grb2
    bind_complex(erbb(loc='C', ty='1') % erbb(loc='C', ty='1') % GRB2(bgap=ANY), 'bcbl', CBL(), 'b', par['ErbB1_Cbl_Grb2'])
    bind_complex(erbb(loc='E', ty='1') % erbb(loc='E', ty='1') % GRB2(bgap=ANY), 'bcbl', CBL(), 'b', par['ErbB1_Cbl_Grb2'])
    

def receptor_internalization_constitutive():
    """Constitutive internalization of unphosphorylated, inactive receptors."""

    Rule('rec_intern_constit_1',
         erbb(bd=None, b=None, loc='C', st='U', ty='1') >>
         erbb(bd=None, b=None, loc='E', st='U', ty='1'),
         par['kint_const_1'])
    
    Rule('recd_intern_constit_1',
         erbb(bd=1, b=None, loc='C', st='U', ty='1') % erbb(bd=1, b=None, loc='C', st='U') >>
         erbb(bd=1, b=None, loc='E', st='U', ty='1') % erbb(bd=1, b=None, loc='E', st='U'),
         par['kint_const_1'])
         
    Rule('rec_intern_constit_2',
         erbb(bd=None, b=None, loc='C', st='U', ty='2') >>
         erbb(bd=None, b=None, loc='E', st='U', ty='2'),
         par['kint_const_2'])
    
    Rule('recd_intern_constit_2',
         erbb(bd=1, b=None, loc='C', st='U', ty='2') % erbb(bd=1, b=None, loc='C', st='U') >>
         erbb(bd=1, b=None, loc='E', st='U', ty='2') % erbb(bd=1, b=None, loc='E', st='U'),
         par['kint_const_2'])
         
    Rule('rec_intern_constit_3',
         erbb(bd=None, b=None, loc='C', st='U', ty='3') >>
         erbb(bd=None, b=None, loc='E', st='U', ty='3'),
         par['kint_const_3'])
    
    Rule('rec_intern_constit_4',
         erbb(bd=None, b=None, loc='C', st='U', ty='4') >>
         erbb(bd=None, b=None, loc='E', st='U', ty='4'),
         par['kint_const_4'])

def receptor_internalization_erbb1_clathrin_med():
    """Internalization of activated EGFR that is clathrin-mediated and requires a GRB2:Cbl interaction (Sorkin and Goh 2009)."""
    
    # Two rates implemented: one for ErbB1-ErbB1 dimers and one for ErbB1 heterodimers.    
    
    Rule('rec_intern_CM_ErbB1ErbB1',
         erbb(loc='C', ty='1') % erbb(loc='C', ty='1') % GRB2(bgap=2, bcbl=ANY) % CBL(b=ANY) >> 
         erbb(loc='E', ty='1') % erbb(loc='E', ty='1') % GRB2(bgap=2, bcbl=ANY) % CBL(b=ANY),
         par['kint_CM_ErbB1ErbB1'])
         
    for i in ['2', '3', '4']:
        Rule('rec_intern_CM_ErbB1ErbB'+i,
             erbb(loc='C', ty='1') % erbb(loc='C', ty=i) % GRB2(bgap=2, bcbl=ANY) % CBL(b=ANY) >> 
             erbb(loc='E', ty='1') % erbb(loc='E', ty=i) % GRB2(bgap=2, bcbl=ANY) % CBL(b=ANY),
             par['kint_CM_ErbB1ErbBX'])

def receptor_internalization_erbb1_clathrin_indepen():
    """Internalization of activated EGFR that is clathrin independent.  
       Clathrin-independent endocytosis of EGFR is typically experimentally observed when high (non-physiological) concentrations of EGF are used.
       The dominant pathway in-vivo is therefore thought to be clathrin-mediated (Sorkin and Goh 2009).
       However, these normally non-physiologically high concentrations of EGF can occur in tumors (Roepstorff 2008)."""
       
    #Two rates implemented: one for ErbB1-ErbB1 dimers and one for ErbB1 heterodimers.
    Rule('rec_intern_NCM_ErbB1ErbB1_Grb2_1',
             erbb(loc='C', ty='1') % erbb(loc='C', ty='1') % GRB2(bgap=2, bcbl=None) >> 
             erbb(loc='E', ty='1') % erbb(loc='E', ty='1') % GRB2(bgap=2, bcbl=None),
             par['kint_NCM_ErbB1ErbB1'])   
    
    Rule("rec_intern_NCM_ErbB1ErbB1_Shc_1",
              erbb(loc='C', ty='1') % erbb(loc='C', ty='1') % SHC(bgap=2) >>
              erbb(loc='E', ty='1') % erbb(loc='E', ty='1') % SHC(bgap=2),
             par['kint_NCM_ErbB1ErbB1'])
     
    Rule("rec_intern_NCM_ErbB1ErbB1_1",
              erbb(loc='C', ty='1', st='P', b=None) % erbb(loc='C', ty='1', st='P', b=None) >>
              erbb(loc='E', ty='1', st='P', b=None) % erbb(loc='E', ty='1', st='P', b=None),
             par['kint_NCM_ErbB1ErbB1']) 
    
       
    for i in ['2', '3', '4']:
        
        Rule('rec_intern_NCM_ErbB1ErbB1_Grb2_'+i,
             erbb(loc='C', ty='1') % erbb(loc='C', ty=i) % GRB2(bgap=2, bcbl=None) >> 
             erbb(loc='E', ty='1') % erbb(loc='E', ty=i) % GRB2(bgap=2, bcbl=None),
             par['kint_NCM_ErbB1ErbBX'])
         
        Rule("rec_intern_NCM_ErbB1ErbB1_Shc_"+i,
              erbb(loc='C', ty='1') % erbb(loc='C', ty=i) % SHC(bgap=2) >>
              erbb(loc='E', ty='1') % erbb(loc='E', ty=i) % SHC(bgap=2),
             par['kint_NCM_ErbB1ErbBX'])
         
        Rule("rec_intern_NCM_ErbB1ErbB1_"+i,
              erbb(loc='C', ty='1', st='P', b=None) % erbb(loc='C', ty=i, st='P', b=None) >>
              erbb(loc='E', ty='1', st='P', b=None) % erbb(loc='E', ty=i, st='P', b=None),
             par['kint_NCM_ErbB1ErbBX'])

def receptor_internalization_erbb234():
    """Internalization of Erb2, ErbB3, and ErbB4 complexes that don't contain ErbB1."""

    for i in ['2', '3', '4']:
        Rule('rec_intern_ErbB2ErbB_Grb2_'+i,
             erbb(loc='C', ty='2', st='P') % erbb(loc='C', ty=i, st='P') % GRB2(bgap=2, bcbl=None) >> 
             erbb(loc='E', ty='2', st='P') % erbb(loc='E', ty=i, st='P') % GRB2(bgap=2, bcbl=None),
             par['kint_NCM_ErbB2ErbBX'])
        
        Rule('rec_intern_ErbB2ErbB_Shc_'+i,
             erbb(loc='C', ty='2', st='P') % erbb(loc='C', ty=i, st='P') % SHC(bgap=2) >> 
             erbb(loc='E', ty='2', st='P') % erbb(loc='E', ty=i, st='P') % SHC(bgap=2),
             par['kint_NCM_ErbB2ErbBX'])
        
        Rule('rec_intern_ErbB2ErbB_'+i,
             erbb(loc='C', ty='2', st='P', b=None) % erbb(loc='C', ty=i, st='P', b=None) >> 
             erbb(loc='E', ty='2', st='P', b=None) % erbb(loc='E', ty=i, st='P', b=None),
             par['kint_NCM_ErbB2ErbBX'])
             
def receptor_degradation_erbb1():
    """Degradation of ErbB1 after endocytosis.  Degradation requires that ligand be bound to ErbB1 and that Cbl be bound to Y1045 of ErbB1 (Roepstorff 2008)."""
    
    degrade(erbb(loc='E', ty='1', bl=ANY, Y1045=ANY), par['kdeg_erbb1'])
    
def receptor_degradation_erbb234():
    """Degradation of ErbB2, ErbB3, and ErbB4 containing complexes."""
    
    for i in ['2', '3', '4']:
        degrade(erbb(loc='E', ty=i, bd=None), par['kdeg_erbb234'])
        degrade(erbb(loc='E', ty='2') % erbb(loc='E', ty=i), par['kdeg_erbb234'])
        
def receptor_recycling_erbb1():
    """Recycling of ErbB1-containing complexes to plasma membrane from endosomes."""
    
    Rule('rec_recycling_ErbB1',
         erbb(loc='E', ty='1', Y1045=None, bl=None, bd=None) >>
         erbb(loc='C', ty='1', Y1045=None, bl=None, bd=None),
         par['krecyc_ErbB1'])    
    
    Rule('rec_recycling_ErbB1d',
         erbb(loc='E', ty='1', Y1045=None, bl=None) % erbb(loc='E', Y1045=None, bl=None) >>
         erbb(loc='C', ty='1', Y1045=None, bl=None) % erbb(loc='C', Y1045=None, bl=None),
         par['krecyc_ErbB1'])

def receptor_recycling_erbb234():
    """Recycling of ErbB2, ErbB3, and ErbB4 complexes to plasma membrane from endosomes."""

    for i in ['2', '3', '4']:
        Rule('rec_recycling_ErbB'+i+'_monomer',
             erbb(loc='E', ty=i, bd=None) >>
             erbb(loc='C', ty=i, bd=None),
             par['krecyc_ErbB234'])
         
        Rule('rec_recycling_ErbB'+i+'_dimer',
             erbb(loc='E', ty='2', bd=1) % erbb(loc='E', ty=i, bd=1) >>
             erbb(loc='C', ty='2', bd=1) % erbb(loc='C', ty=i, bd=1),
             par['krecyc_ErbB234'])