"""
Overview:
========

PySB implementation of the ErbB related MAPK and AKT signaling
pathways originally published in [Chen2009]_.

This file contains functions that implement the ErbB execution
pathway in four modules:

- Receptor layer events, taking into account all ErbB1-4 interactions with ligand.
- AKT pathway
- MAPK pathway
- crosstalk between AKT/MAPK pathways.
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

# Receptor Layer, 0

# Rates obtained from Chen et al. Table 1 pg. 5
KF = 1e-5
KR = 1e-1
KCP = 1e-1
KCD = 1e-2
KDIMF = 1.6e-6
KDIMR = 1.6e-1
KINTF = 1.3e-3
KINTR = 5.0e-5
KDEG = .1

from parameter_dict_A431 import parameter_dict as par
        
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
    
def rec_monomers_scaffold_proteins():
    """Declares the monomers for scaffolding proteins that may be shared between pathways."""
    Monomer('SHC', ['bgap', 'bgrb', 'batp', 'st'], {'st':['U','P']})
    Monomer('GRB2', ['b', 'bsos', 'bgap', 'bgab1', 'bcbl'])
    Monomer('GAB1', ['bgrb2', 'bshp2', 'bpi3k', 'batp','bERKPP','bPase9t','S'],{'S':['U','P','PP']})

def rec_initial_lig_hEGF():
    """Declares the initial conditions for high EGF (5 nM)."""
    Parameter('EGF_0',     3.01e12) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    Initial(EGF(b=None, st='M'), EGF_0)
    
def rec_initial_lig_lEGF():
    """Declares the initial conditions for low EGF (.01 nM)."""
    Parameter('EGF_0',      6.02e9) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    Initial(EGF(b=None, st='M'), EGF_0)

def rec_initial_lig_hHRG():
    """Declares the initial conditions for high heregulin (5 nM)."""
    Parameter('HRG_0',         3.01e12) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    Initial(HRG(b=None, st='M'), HRG_0)

def rec_initial_lig_lHRG():
    """Declares the initial conditions for low heregulin (.01 nM)."""
    Parameter('HRG_0',         6.02e9) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    Initial(HRG(b=None, st='M'), HRG_0)
    
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
    
    print bind_complex(erbb(ty='1', st='U', loc='C', bl=ANY, b=None, bd=1) % erbb(st='U', loc='C', b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', ATP(b=None), 'b', par['ErbB1_bind_ATP'], m1=erbb(ty='1', st='U', loc='C', bl=ANY, b=None, bd=1))  
    
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

def mapk_monomers():
    Monomer('SOS', ['bgrb', 'bras', 'bERKPP', 'st'], {'st':['U', 'P']})
    Monomer('RAS', ['bsos', 'braf', 'bpi3k', 'st'], {'st':['GDP', 'GTP']})
    Monomer('RAF', ['b', 'st', 'ser295'], {'st':['U', 'P'], 'ser295':['U', 'P']})
    Monomer('PP1', ['b'])
    Monomer('PP2', ['b'])
    Monomer('PP3', ['b'])
    Monomer('MEK', ['b', 'st'], {'st':['U', 'P', 'PP']})
    Monomer('ERK', ['b', 'st', 'loc'], {'st':['U', 'P', 'PP'], 'loc':['C', 'N']})

def mapk_initial():
    # Initial values declared in parameter dictionary for given cell type.
    alias_model_components()

    Initial(RAS(bsos=None, braf=None, bpi3k=None, st='GDP'), RAS_0)
    Initial(RAF(b=None, st='U', ser295='U'), RAF_0)
    Initial(MEK(b=None, st='U'), MEK_0)
    Initial(ERK(b=None, st='U', loc='C'), ERK_0)
    Initial(PP1(b=None), PP1_0)
    Initial(PP2(b=None), PP2_0)
    Initial(PP3(b=None), PP3_0)
    Initial(GRB2(b=None, bsos=1, bgap=None, bgab1=None, bcbl=None) % SOS(bgrb=1, bras=None, bERKPP=None, st='U'), GRB2_SOS_0)
    Initial(ERK(b=None, st='PP', loc='C'), ERKPP_0)

def mapk_events():

    # =====================
    # Alias model components for names in present namespace
    alias_model_components()

    #Bind GRB2 to ErbB dimers with SOS already bound:
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=2, bgab1=None, bgap=None) % SOS(bras=None, bERKPP=None, bgrb=2), 'bgap', par['GRB2_SOS_bind_GAP'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=2, bgab1=None, bgap=None, bcbl=None) % SOS(bras=None, bERKPP=None, bgrb=2), 'bgap', par['GRB2_SOS_bind_GAP'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    # Shc:P bound to ErbB dimers can bind Grb2-Sos
    bind_complex(SHC(batp=None, st='P', bgrb=None, bgap=ANY), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=1, b=None, bcbl=None) % SOS(bras=None, bERKPP=None, bgrb=1), 'b', par['GRB2_SOS_bind_SHCP_GAP'])

    # SHC:P can bind GRB2-SOS without being attached to an ErbB dimer: 
    bind_complex(SHC(batp=None, st='P', bgrb=None, bgap=None), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=1, b=None, bcbl=None) % SOS(bras=None, bERKPP=None, bgrb=1), 'b', par['SHCP_bind_GRB2SOS'])

    # ErbB dimers can bind the free SHC:P-GRB2-SOS complex:
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=3, b=2, bcbl=None) % SOS(bras=None, bERKPP=None, bgrb=3), 'bgap', par['GAP_bind_SHCP_GRB2_SOS'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=3, b=2, bcbl=None) % SOS(bras=None, bERKPP=None, bgrb=3), 'bgap', par['GAP_bind_SHCP_GRB2_SOS'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))

    # GRB2 and SOS bind/disassociate:
    bind(GRB2(bgap=None, bgab1=None, b=None), 'bsos', SOS(bras=None, bERKPP=None), 'bgrb', par['GRB2_bind_SOS'])

    #Although no free SOS is present initially in Chen Sorger model, GRB2-SOS can disassociate (see above), so these are necessary.
    # SOS binds to ErbBdimer-SHC:P-GRB2  
    bind_complex(SHC(batp=None, st='P', bgrb=1, bgap=ANY) % GRB2(b=1, bgap=None, bgab1=None, bcbl=None), 'bsos', SOS(bras=None, bERKPP=None), 'bgrb', par['SOS_bind_GAP_SHCP_GRB2'])

    #SOS binds SHC:P-GRB2 without complex
    bind_complex(SOS(bras=None, bERKPP=None, bgrb=None), 'bgrb', SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, b=1, bcbl=None), 'bsos', par['SOS_bind_SHCP_GRB2'])

    # SOS also binds ErbBdimer-GRB2
    bind_complex(GRB2(bgap=ANY, bgab1=None, b=None, bsos=None), 'bsos', SOS(bras=None, bgrb=None, bERKPP=None), 'bgrb', par['SOS_bind_GAP_GRB2'])

    # ErbBdimer-GRB2-SOS and ErbBdimer-SHC:P-GRB2-SOS can bind either Ras-GDP or Ras-GTP: k1, k1r
    bind_complex(GRB2(bgap=ANY, bgab1=None, b=None, bsos=1) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GDP', bpi3k=None), 'bsos', par['RASGDP_bind_bound_GRB2_SOS'])

    bind_complex(SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcbl=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GDP', bpi3k=None), 'bsos', par['RASGDP_bind_bound_GRB2_SOS'])

    bind_complex(GRB2(bgap=ANY, bgab1=None, b=None, bsos=1) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GTP', bpi3k=None), 'bsos', par['RASGTP_bind_bound_GRB2_SOS'])

    bind_complex(SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcbl=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GTP', bpi3k=None), 'bsos', par['RASGTP_bind_bound_GRB2_SOS'])

    #Ras-GDP --> Ras-GTP transition (GDP exchange) occurs at a different (faster) rate when bound to Sos (a guanine exchange factor - GEF) than when unbound
    #Ras-GTP --> Ras-GDP transition (GTP hydrolysis) is also covered by these rules.  A GAP (GTPase activating protein) would theoretically affect this rate.

    Rule('RASGTP_to_GDP_SOS_GRB2_SHCP_complex',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcbl=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', bpi3k=None) <>
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcbl=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', bpi3k=None),
         *par['RASGTP_unbind_GRB2_SOS'])

    Rule('RASGTP_to_GDP_SOS_GRB2_GAP_complex',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', bpi3k=None) <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', bpi3k=None),
         *par['RASGTP_unbind_GRB2_SOS'])
    
    #Ras has its own intrinsic (slower) GTPase and GDP exchange rates when it is unbound.
    equilibrate(RAS(braf=None, bsos=None, st='GTP', bpi3k=None), RAS(braf=None, bsos=None, st='GDP', bpi3k=None), par['Ras_intrinsic_function'])

    # Activation of RAF -> RAF:P by RAS-GTP
    catalyze_state(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', RAF(ser295='U'), 'b', 'st', 'U', 'P', par['RASGTP_bind_RAF']+par['RASGTP_RAF_cat'])
    
    # Deactivation of RAF:P -> RAF by PP1
    catalyze(PP1(), 'b', RAF(st='P', ser295='U'), 'b', RAF(st='U', ser295='U'),
             (par['RAFP_PP1']))

    # Activation of MEK -> MEK:P by activated RAF
    catalyze(RAF(st='P', ser295='U'), 'b', MEK(st='U'), 'b', MEK(st='P'),
             (par['RAFP_MEK']))

    # Deactivation of MEK:P -> MEK by PP2
    catalyze(PP2(), 'b', MEK(st='P'), 'b', MEK(st='U'),
             (par['MEKP_PP2']))
    
    # Activation of MEK:P -> MEK:P:P by activated RAF
    catalyze(RAF(st='P', ser295='U'), 'b', MEK(st='P'), 'b', MEK(st='PP'),
             (par['RAFP_MEKP']))

    # Deactivation of MEK:P:P -> MEK:P by PP2
    catalyze(PP2(), 'b', MEK(st='PP'), 'b', MEK(st='P'),
             (par['MEKPP_PP2']))
    
    # Activation of ERK -> ERK:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='U', loc='C'), 'b', ERK(st='P', loc='C'),
             (par['MEKPP_ERK']))

    # Deactivation of ERK:P -> ERK by PP3
    catalyze(PP3(), 'b', ERK(st='P', loc='C'), 'b', ERK(st='U', loc='C'),
             (par['ERKP_PP3']))

    # Activation of ERK:P -> ERK:P:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='P', loc='C'), 'b', ERK(st='PP', loc='C'),
             (par['MEKPP_ERKP']))

    # Deactivation of ERK:P:P -> ERK:P by PP3
    catalyze(PP3(), 'b', ERK(st='PP', loc='C'), 'b', ERK(st='P', loc='C'),
             (par['ERKPP_PP3']))

    # Degradation of PP3
    degrade(PP3(b=None), par['PP3_deg'])

def mapk_negative_feedback():
    """Negative feedback events within the MAPK pathway."""
    
    #ERK:P:P phosphorylates GRB2-SOS, preventing RAS-GDP->RAS-GTP conversion
    catalyze_state(ERK(st='PP', loc='C'), 'b', SOS(bras=None), 'bERKPP', 'st', 'U', 'P', (par['ERKPP_phos_SOS']))

def akt_monomers():
    """ This is the akt part of the pathway from the Chen et al. 2009 paper.  Initial rules for all binding reactions were generated and then coded again using macros and higher order macros.  Initial parameters and conditions were taken from Chen et al. 2009 paper and supplementary, but were later modified in order to get the model working correctly.  This pathway follows AKT from its initial state to a phosphorylated and then double phosphorylated state before returning to unphosphorylated AKT.  The model works correctly, but parameters and rates may need to be modified in order to get best fit.  Parameters and rates included are from trial and error of what best fit the model.  
"""
    #This pathway originally coded by Tim O'Brien.
    Monomer('PI3K',['bgab1','bpip', 'bras', 'berb'])
    Monomer('SHP2',['bgab1'])
    Monomer('PIP', ['bakt', 'bpdk1', 'S', 'bpi3k'], {'S':['PIP2', 'PIP3']})
    Monomer('PTEN', ['bpip3'])
    Monomer('SHP', ['bpip3'])
    Monomer('AKT', ['bpip3', 'bpdk1', 'S'], {'S':['U', 'P', 'PP']})
    Monomer('PDK1', ['bakt', 'bpip3'])
    Monomer('PP2A_III', ['bakt'])

def akt_initial():
    # See parameter dictionary files for given cell type for initial values.
    
    alias_model_components()
    
    # Initial conditions 
    Initial(PI3K(bgab1=None, bpip=None, bras=None, berb=None), PI3K_0)
    Initial(SHP2(bgab1=None), SHP2_0)
    Initial(PIP(bakt=None, bpdk1=None, S='PIP2', bpi3k=None), PIP_0)
    Initial(PTEN(bpip3=None), PTEN_0)
    Initial(SHP(bpip3=None), SHP_0)
    Initial(AKT(bpip3=None, bpdk1=None, S='U'), AKT_0)
    Initial(PDK1(bakt=None, bpip3=None), PDK1_0)
    Initial(PP2A_III(bakt=None), PP2A_III_0)
    Initial(AKT(bpip3=None, bpdk1=None, S='PP'), AKTPP_0)
    
def akt_events():

    #ErbBdimer-GRB2-GAB1 phosphorylation - Rates from Table p. 5 Chen et al 2009
    bind_complex(erbb(bd=1) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U'), 'batp', ATP(b=None), 'b', par['GAB1_bind_ATP'])

    Rule('GAB1_phos',
         erbb(bd=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1) >>
         erbb(bd=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + ADP(),
         par['GAB1_phos'])

    #SHP2 can desphosphorylate GAB1-P
    catalyze_state(SHP2(), 'bgab1', GAB1(bgrb2=ANY, bpi3k=None, batp=None, bERKPP=None, bPase9t=None), 'bshp2', 'S', 'P', 'U', (par['SHP2_dephos_GAB1P']))
   
    #After GAB1 phosphorylation, all receptor dimer combinations can bind a single PI3K
    #Chen/Sorger model gives two rate constant sets for different receptor dimers:
    #Rate 1: ErbB1/ErbB1, ErbB1/ErbB2, ErbB1/ErbB4, and ErbB2/ErbB4 dimers:
    
    for i in ['1', '2', '4']:
        bind_complex(erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P'), 'bpi3k', PI3K(bpip=None, bgab1=None, bras=None, berb=None), 'bgab1', par['GAB1_bind_PI3K_1'])
                     
    bind_complex(erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='4') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P'), 'bpi3k', PI3K(bpip=None, bgab1=None, bras=None, berb=None), 'bgab1', par['GAB1_bind_PI3K_1'])

    #Rate 2: ErbB1/ErbB3, ErbB2/ErbB2, and ErbB2/ErbB3 dimers:
    for i in ['1', '2']:
        bind_complex(erbb(bd=1, ty=i) % erbb(bd=1, ty='3') % GRB2(b=None, bsos=None, bgap=None, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P'), 'bpi3k', PI3K(bpip=None, bgab1=None, bras=None, berb=None), 'bgab1', par['GAB1_bind_PI3K_2'])

    bind_complex(erbb(bd=1, ty='2') % erbb(bd=1, ty='2') % GRB2(b=None, bsos=None, bgap=None, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P'), 'bpi3k', PI3K(bpip=None, bgab1=None, bras=None, berb=None), 'bgab1', par['GAB1_bind_PI3K_2']) 
    
    #ErbB2-ErbB3 dimers contain 6 binding domains for PI3K (don't need to be bound to adaptor complex).  This set of reactions assumes sequential binding to 6 sites, which is almost certainly not biologically true.  However, this drastically lowers the number of possible species.  Individual parameters are assigned to each sequential binding so that effective rate parameters can be fit.
    """This is the improved ErbB2/ErB3-PI3K sequence of events -- different from Chen Sorger 2009"""

    bind_complex(erbb(bd=1, ty='3', b=None, st='P', pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(bd=1, ty='2', b=None, st='P'), 'pi3k1', PI3K(bgab1=None, bpip=None, bras=None), 'berb', par['ErbB23_bind_PI3K_1'], m1=erbb(bd=1, ty='3', b=None, st='P', pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))
    bind_complex(erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(bd=1, ty='2', b=None, st='P'), 'pi3k2', PI3K(bgab1=None, bpip=None, bras=None), 'berb', par['ErbB23_bind_PI3K_2'], m1=erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))
    bind_complex(erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(bd=1, ty='2', b=None, st='P'), 'pi3k3', PI3K(bgab1=None, bpip=None, bras=None), 'berb', par['ErbB23_bind_PI3K_3'], m1=erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k4=None, pi3k5=None, pi3k6=None))
    bind_complex(erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k5=None, pi3k6=None) % erbb(bd=1, ty='2', b=None, st='P'), 'pi3k4', PI3K(bgab1=None, bpip=None, bras=None), 'berb', par['ErbB23_bind_PI3K_4'], m1=erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k5=None, pi3k6=None))
    bind_complex(erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k6=None) % erbb(bd=1, ty='2', b=None, st='P'), 'pi3k5', PI3K(bgab1=None, bpip=None, bras=None), 'berb', par['ErbB23_bind_PI3K_5'], m1=erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k6=None))
    bind_complex(erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY) % erbb(bd=1, ty='2', b=None, st='P'), 'pi3k6', PI3K(bgab1=None, bpip=None, bras=None), 'berb', par['ErbB23_bind_PI3K_6'], m1=erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY))

    #PI3K bound directly to ErbB2-ErbB3 dimers (at any of 6 sites) can catalyze PIP2->PIP3:
    
    Rule('ErbB23_PI3K_cat_1',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_1'])
    
    Rule('ErbB23_PI3K_cat_2',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_2'])
    
    Rule('ErbB23_PI3K_cat_3',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_3'])
    
    Rule('ErbB23_PI3K_cat_4',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=None, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_4'])
    
    Rule('ErbB23_PI3K_cat_5',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=None) + PIP(bakt=None, bpdk1=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_5'])
    
    Rule('ErbB23_PI3K_cat_6',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=ANY) + PIP(bakt=None, bpdk1=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=ANY) + PIP(bakt=None, bpdk1=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_6'])
    
    #PI3K bound to complex catalyzes PIP2 -> PIP3
    #Two rate sets for initial binding in Chen/Sorger model:
    #Rate 1: ErbB1/ErbBX dimers:
    
    bind_complex(erbb(bd=ANY, ty='1') % erbb(bd=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None), 'bpip', PIP(S='PIP2', bpdk1=None, bakt=None, bpi3k=None), 'bpi3k', par['PIP2_bind_PI3K_1'])

    #Rate 2: ErbB2/ErbBX dimers, X=2, 3, 4:
    for i in ['2', '3', '4']:
        bind_complex(erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None), 'bpip', PIP(S='PIP2', bpdk1=None, bakt=None, bpi3k=None), 'bpi3k', par['PIP2_chain_PI3K'])
    
    #Two catalysis rates in Chen/Sorger model:
    #Rate 1: ErbB2/ErbB3 dimers:
    Rule('PIP2_PI3K_catalysis_1',
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='3') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', bpdk1=None, bakt=None, bpi3k=1) >>
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='3') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', bpdk1=None, bakt=None, bpi3k=None),
         par['PIP2_self_catalysis'])

    #Rate 2: All other dimers:
    for i in ['1', '2', '3', '4']:
        Rule('PIP2_PI3K_catalysis_2_'+i,
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', bpdk1=None, bakt=None, bpi3k=1) >>
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', bpdk1=None, bakt=None, bpi3k=None),
             par['PIP2_PI3K_catalysis'])

    for i in ['2', '4']:
        Rule('PIP2_PI3K_catalysis_3_'+i,
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', bpdk1=None, bakt=None, bpi3k=1) >>
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', bpdk1=None, bakt=None, bpi3k=None),
             par['PIP2_PI3K_catalysis'])
             
    # Setting up the binding reactions necessary for AKT to be phosphorylated and move through the pathway
    bind_table([[                                                 AKT(S='U', bpdk1=None),       AKT(S='P', bpdk1=None)],
                [PIP(S='PIP3', bpdk1=None, bpi3k=None),       (par['PIP3_bind_AKT']),     (par['PIP3_bind_AKT'])]],
                'bakt', 'bpip3')
    
    # AKT-PIP3 is phosphorylated by PDK1 to AKTP; PDK1-PIP3 and AKTP are released
    bind(PDK1(bpip3=None), 'bakt', AKT(bpip3=ANY, S='U'), 'bpdk1', par['AKT_PIP3_bind_PDK1'])
    
    Rule('PDK1_AKT_catalysis',
         PDK1(bpip3=None, bakt=1) % AKT(bpip3=2, S='U', bpdk1=1) % PIP(S='PIP3', bpdk1=None, bpi3k=None, bakt=2) >>
         PDK1(bpip3=3, bakt=None) % PIP(S='PIP3', bpdk1=3, bpi3k=None, bakt=None) + AKT(bpip3=None, S='P', bpdk1=None),
         par['PDK1_AKT_catalysis'])

    # PIP3 unbinds PDK1
    bind(PIP(S='PIP3', bakt=None, bpi3k=None), 'bpdk1', PDK1(bakt=None), 'bpip3', par['PIP3_bind_PDK1'])

    # AKTP-PIP3 is phosphorylated by PDK1 to AKTPP
    bind(PDK1(bpip3=None), 'bakt', AKT(bpip3=ANY, S='P'), 'bpdk1', par['AKT_PIP3_bind_PDK1'])
    
    Rule('PDK1_AKTP_catalysis',
         PDK1(bpip3=None, bakt=1) % AKT(bpip3=2, S='P', bpdk1=1) % PIP(S='PIP3', bpdk1=None, bpi3k=None, bakt=2) >>
         PDK1(bpip3=3, bakt=None) % PIP(S='PIP3', bpdk1=3, bpi3k=None, bakt=None) + AKT(bpip3=None, S='PP', bpdk1=None),
         par['PDK1_AKTP_catalysis'])

    # AKTP is dephosphorylated by PP2A-III back to AKT
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'bpdk1', 'S', 'P', 'U',(par['AKTP_dephos']))
   
    # AKTPP is dephosphorylated by PP2A-III back to AKTP
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'bpdk1', 'S', 'PP', 'P',(par['AKTPP_dephos']))

    # PIP3 is dephosphorylated by PTEN to PIP2
    catalyze_state(PTEN, 'bpip3', PIP(bakt=None, bpi3k=None), 'bpdk1', 'S', 'PIP3', 'PIP2', (par['PIP3_dephos']))

    # PIP3 is dephosphorylated by SHP to PIP2
    catalyze_state(SHP, 'bpip3', PIP(bakt=None, bpi3k=None), 'bpdk1', 'S', 'PIP3', 'PIP2', (par['PIP3_dephos']))

def mtor_monomers():
    Monomer('mTOR', ['bcomplex', 'bcat', 'S2448', 'bFKBP38'], {'S2448':['U','P']})
    Monomer('TORC1_ptns', ['bmTOR']) #A complex composed of Raptor, mLST8, PRAS40, and DEPTOR.  With mTOR becomes mTORC1 (mTOR complex 1)
    Monomer('TORC2_ptns', ['bmTOR']) #A complex composed of Rictor, GbetaL, and mSIN1, mLST8, and DEPTOR. With mTOR becomes mTORC2 (mTOR complex 2)    
      
def mtor_initial():
    alias_model_components()
    Initial(mTOR(bcomplex=None, bcat=None, bFKBP38=None, S2448='U'), mTOR_0)
    Initial(TORC1_ptns(bmTOR=None), TORC1_ptns_0)
    Initial(TORC2_ptns(bmTOR=None), TORC2_ptns_0)

def mtor_complex_formation():
    """Formation of mTOR complexes 1 and 2. """    
    
    alias_model_components()
    # mTOR can bind either the proteins in mTORC1 or mTORC2:
    bind(mTOR(bcat=None, S2448='U', bFKBP38=None), 'bcomplex', TORC1_ptns(), 'bmTOR', (par['mTOR_bind_TORC1ptns']))
    
    bind(mTOR(bcat=None, S2448='U', bFKBP38=None), 'bcomplex', TORC2_ptns(), 'bmTOR', (par['mTOR_bind_TORC2ptns']))
         
def tsc2_monomers():
    Monomer('AMPK', ['T172', 'b'], {'T172':['U', 'P']}) #Phosphorylated at T172 by LKB1
    Monomer('LKB1', ['b', 'S431'], {'S431':['U', 'P']})
    Monomer('TSC', ['b', 'S664', 'S1798', 'S1387', 'S939', 'S981', 'T1462'], {'T1462': ['U', 'P'], 'S664':['U', 'P'], 'S1798':['U','P'], 'S1387':['U','P'], 'S939':['U', 'P'], 'S981':['U', 'P']}) #Composed of both TSC1 and TSC2
    Monomer('RSK1', ['b', 'T573', 'S380', 'S221'], {'T573':['U','P'], 'S380':['U','P'], 'S221':['U','P']})

def tsc2_initial():
    alias_model_components()
    Initial(AMPK(T172='U', b=None), AMPK_0)
    Initial(LKB1(b=None, S431='U'), LKB1_0)
    Initial(TSC(b=None, S664='U', S1798='U', S1387='U', S939='U', S981='U', T1462='U'), TSC_0)
    Initial(RSK1(b=None, T573='U', S380='U', S221='U'), RSK1_0)

def tsc2_inhibition_by_akt():
    """Inhibition of TSC2 by AKT:PP."""
    
    alias_model_components()
    # AKT:PP phosphorylates TSC2 at S939, S981, and T1462, which causes it to translocate from the membrane to the cytosol and stops TSC2's GAP activity on Rheb (Cai 2006, Journal of Cell Biology).
    bind(AKT(S='PP', bpip3=None), 'bpdk1', TSC(S939='U', S981='U', T1462='U'), 'b', par['AKTPP_bind_TSC2'])
    
    for site in ['S939', 'S981', 'T1462']:
        Rule('AKTPP_phos_TSC2_'+site,
        AKT(S='PP', bpip3=None, bpdk1=1) % TSC({'b':1, site:'U'}) >>
        AKT(S='PP', bpip3=None, bpdk1=None) + TSC({'b':None, site:'P'}),
        par['AKTPP_phos'+site+'_TSC2'])

def tsc2_inhibition_by_erk():
    """Inhibition of TSC2 by ERK:PP.  This includes both direct phosphorylation of TSC2 by ERK and also indirect inhibition through ERK phosphorylation of RSK1."""
    
    alias_model_components()
    # ERK:PP can also phosphorylate TSC2 at S664, again inhibiting its GAP activity.
    catalyze_state(ERK(st='PP', loc='C'), 'b', TSC(S939='U', S981='U', T1462='U'), 'b', 'S664', 'U', 'P', (par['ERKPP_phos_TSC2']))

    # ERK:PP phosphorylates RSK1 at T573. RSK1 then autocatalyzes phosphorylation at S380, which allows binding of PDK1, which phosphorylates S221, giving fully active RSK1.
    catalyze_state(ERK(st='PP', loc='C'), 'b', RSK1(S380='U', S221='U'), 'b', 'T573', 'U', 'P', (par['ERKPP_phos_RSK1']))

    Rule('RSK1_autocatalysis',
         RSK1(S380='U', S221='U', T573='P') >>
         RSK1(S380='P', S221='U', T573='P'),
         par['RSK1_autocat'])

    catalyze_state(PDK1(bakt=None), 'bpip3', RSK1(S380='P', T573='P'), 'b', 'S221', 'U', 'P', (par['PDK1_phos_RSK1']))

    # Active RSK1 can phosphorylate TSC2 at S1798, inhibiting its GAP activity.
    catalyze_state(RSK1(S380='P', T573='P', S221='P'), 'b', TSC(S939='U', S981='U', T1462='U'), 'b', 'S1798', 'U', 'P', (par['RSK1_phos_TSC2']))
    
def tsc2_activation_by_erk():
    """Indirect activation of TSC2 GTPase activating protein function by ERK:PP through ERK:PP->RSK1->LKB1->AMPK kinase cascade."""
    
    alias_model_components()
    # Active RSK1 phosphorylates LKB1 at S431, activating LKB1.
    catalyze_state(RSK1(S380='P', T573='P', S221='P'), 'b', LKB1(), 'b', 'S431', 'U', 'P', (par['RSK1_phos_LKB1']))

    # Active LKB1 phosphorylates AMPK at T172, activating it.
    catalyze_state(LKB1(S431='P'), 'b', AMPK(), 'b', 'T172', 'U', 'P', (par['LKB1_phos_AMPK']))

    # Active AMPK phosphorylates S1387 on TSC2, activating its GAP activity.
    catalyze_state(AMPK(T172='P'), 'b', TSC(S939='U', S981='U', T1462='U'), 'b', 'S1387', 'U', 'P', (par['AMPK_phos_TSC2'])) 

def tsc2_gap_function():
    """GTPase Activating Protein activity of TSC2 on Rheb GTPase."""
    
    alias_model_components()
    # TSC2 can bind Rheb, inhibiting its GTPase activity if TSC2 is phosphorylated on S664 or S1798 (if phosphorylated on S939, S981, or T1462, TSC2 translocates from the membrane to the cytosol and can't bind Rheb at all).
    bind(TSC(S939='U', S981='U', T1462='U'), 'b', Rheb(S='GTP', bmTOR=None), 'bTSC', par['TSC2_bind_Rheb'])

    Rule('TSC2_Rheb',
         TSC(b=1, S1387='U', S664='U', S1798='U') % Rheb(S='GTP', bTSC=1, bmTOR=None) >>
         TSC(b=1, S1387='U', S664='U', S1798='U') % Rheb(S='GDP', bTSC=1, bmTOR=None),
         par['TSC2_Rheb_GTPase'])
    
    for site in ['S664', 'S1798']:
        Rule('TSC2_'+site+'_Rheb',
             TSC({'b':1, site:'P'}) % Rheb(S='GTP', bTSC=1, bmTOR=None) >>
             TSC({'b':1, site:'P'}) % Rheb(S='GDP', bTSC=1, bmTOR=None),
             par['TSC2p_Rheb_GTPase'])

    # If TSC2 is phosphorylated on S1387, its GAP activity is increased.
    Rule('TSC2_S1387_Rheb',
         TSC(b=1, S1387='P') % Rheb(S='GTP', bTSC=1, bmTOR=None) >>
         TSC(b=1, S1387='P') % Rheb(S='GDP', bTSC=1, bmTOR=None),
         par['TSC2pS1387_Rheb_GTPase'])

def mtorc1_signaling_monomers():
    Monomer('Rheb', ['bTSC', 'bFKBP38', 'S', 'bmTOR'], {'S':['GDP', 'GTP']})
    Monomer('FKBP38', ['b'])
    Monomer('S6K', ['T252', 'T412', 'b'], {'T252':['U', 'P'], 'T412':['U','P']})
    Monomer('rpS6', ['b', 'S'], {'S':['U', 'P']})
    Monomer('EIF4EBP1', ['bEIF4E', 'bmTOR', 'S'], {'S':['U','P']})
    Monomer('EIF4E', ['b'])

def mtorc1_signaling_initial():
    alias_model_components()    
    Initial(Rheb(bTSC=None, bFKBP38=None, S='GTP', bmTOR=None), Rheb_0)
    Initial(FKBP38(b=None), FKBP38_0)
    Initial(S6K(b=None, T252='U', T412='U'), S6K_0)
    Initial(rpS6(b=None, S='U'), rpS6_0)
    Initial(EIF4EBP1(bEIF4E=None, bmTOR=None, S='U'), EIF4EBP1_0)
    Initial(EIF4E(b=None), EIF4E_0)

def mtorc1_signaling():
    """Activation of mTORC1 by Rheb GTPase and downstream signaling (mTORC1->S6K->rpS6 and mTORC1->EIF4EBP1)."""
    
    alias_model_components()
    # Rheb possesess its own intrinsic GTPase activity.
    Rule('Rheb_GTPase',
         Rheb(bTSC=None, S='GTP', bFKBP38=None, bmTOR=None) >>
         Rheb(bTSC=None, S='GDP', bFKBP38=None, bmTOR=None),
         par['Rheb_intrinsic_GTPase'])

    # Replacement of GDP with GTP in Rheb site (to give cycle) -- this rxn should not be rate limiting:
    Rule('Rheb_GDP_GTP',
         Rheb(S='GDP') >>
         Rheb(S='GTP'),
         par['Rheb_GDP_GTP'])

        # Rheb-GTP phosphorylates mTOR in mTORC1 on S2448, activating it.
    bind_complex(TORC1_ptns(bmTOR=1) % mTOR(bcomplex=1, S2448='U', bFKBP38=None), 'bcat', Rheb(S='GTP', bTSC=None), 'bmTOR', par['Rheb_bind_mTORC1'])

    Rule('RhebGTP_cat_mTOR',
         mTOR(bcomplex=ANY, bcat=2, S2448='U', bFKBP38=None) % Rheb(S='GTP', bTSC=None, bmTOR=2) >>
         mTOR(bcomplex=ANY, bcat=None, S2448='P', bFKBP38=None) + Rheb(S='GTP', bTSC=None, bmTOR=None),
         par['RhebGTP_phos_mTOR'])

    # FKBP38 can bind mTOR in mTORC1, preventing its activity.
    bind_complex(TORC1_ptns(bmTOR=1) % mTOR(bcomplex=1, bcat=None), 'bFKBP38', FKBP38(), 'b', par['FKBP38_bind_mTOR'])

    # Rheb-GTP can also bind FKBP38, keeping it from inhibiting mTOR.
    bind(Rheb(S='GTP', bmTOR=None), 'bFKBP38', FKBP38(), 'b', par['Rheb_bind_FKBP38'])
    
    # Active mTORC1 phosphorylates S6K at T412, allowing PDK1 to phosphorylate S6K at T252.
    catalyze_state(mTOR(bcomplex=ANY, bFKBP38=None, S2448='P'), 'bcat', S6K(T252='U'), 'b', 'T412', 'U', 'P', (par['mTORC1_phos_S6K']))

    catalyze_state(PDK1(bakt=None), 'bpip3', S6K(T412='P'), 'b', 'T252', 'U', 'P', (par['PDK1_phos_S6K']))

    # Active S6K phosphorylates rpS6.
    catalyze_state(S6K(T252='P', T412='P'), 'b', rpS6(), 'b', 'S', 'U', 'P', (par['S6K_phos_rpS6']))

    # Unphosphorylated EIF4EBP1 binds EIF4E (EIF4E is necessary for mRNA translation, which this binding interaction prevents).
    bind(EIF4EBP1(bmTOR=None, S='U'), 'bEIF4E', EIF4E(), 'b', par['EIF4EBP1_bind_EIF4E'])

    # Active mTORC1 phosphorylates EIF4EBP1, preventing its interaction with EIF4E and activating mRNA translation.
    catalyze_state(mTOR(bcomplex=ANY, bFKBP38=None, S2448='P'), 'bcat', EIF4EBP1(bEIF4E=None), 'bmTOR', 'S', 'U', 'P', (par['mTORC1_phos_EIF4EBP1']))

def mtorc2_signaling():
    """Phosphorylation of Akt:P->Akt:PP by mTORC2."""
    
    alias_model_components()
    # mTORC2 phosphorylates AKT:P -> AKT:PP
    bind_complex(TORC2_ptns(bmTOR=1) % mTOR(bcomplex=1), 'bcat', AKT(S='P', bpip3=None), 'bpdk1', par['mTORC2_bind_AKTP'])

    Rule('mTOR_AKT_cat',
         TORC2_ptns(bmTOR=1) % mTOR(bcomplex=1, bcat=2) % AKT(S='P', bpdk1=2, bpip3=None) >>
         TORC2_ptns(bmTOR=1) % mTOR(bcomplex=1, bcat=None) + AKT(S='PP', bpdk1=None, bpip3=None),
         par['mTORC2_cat_AKTPP'])

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
    Monomer('Bcl2', ['bf'])
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
    Initial(Bcl2(bf=None), Bcl2_0)
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
    catalyze_state(AKT(S='PP', bpip3=None), 'bpdk1', RAF(st='P'), 'b', 'ser295', 'U', 'P', (par['AKTPP_phos_RAFP']))

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
        
        
