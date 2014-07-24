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
        
# Monomer declarations
# ====================

def rec_monomers():

    """ Declares the ErbB receptor interaction monomers (except ligands).
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    Monomer('erbb', ['bl', 'bd', 'b', 'ty', 'st', 'loc', 'pi3k1', 'pi3k2', 'pi3k3', 'pi3k4', 'pi3k5', 'pi3k6', 'cpp'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E'], 'cpp':['Y', 'N']}) # bl: lig, bd: dimer, b: binding, ty: rec type, st: (U)n(P)hosphorylated, loc: (C)yto 'brane or (E)ndosome 'brane, cpp: No real biophysical meaning; useful model marker for presence of CPP bound downstream.
    Monomer('DEP', ['b'])
    Monomer('ATP', ['b'])
    Monomer('ADP')
    Monomer('CPP', ['b', 'loc'], {'loc':['C', 'E']})

def rec_monomers_lig_EGF():
    """Declares the monomer for the ligand EGF."""
    Monomer('EGF', ['b', 'st'], {'st':['M', 'E']}) # Epidermal Growth Factor ligand

def rec_monomers_lig_HRG():
    """Declares the monomer for the ligand heregulin."""
    Monomer('HRG', ['b']) # Heregulin ligand

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
    Initial(HRG(b=None), HRG_0)

def rec_initial_lig_lHRG():
    """Declares the initial conditions for low heregulin (.01 nM)."""
    Parameter('HRG_0',         6.02e9) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    Initial(HRG(b=None), HRG_0)
    
def rec_initial():
    """Declares the initial conditions for all monomers in receptor interactions except ligands."""
    # # Initial concentrations (except DEP1) for all cell types taken from Chen et al 2009 -- see Jacobian files
    alias_model_components()
    Initial(erbb(bl=None, bd=None, b=None, ty='1', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb1_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='2', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb2_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='3', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb3_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='4', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb4_0)
    Initial(ATP(b=None), ATP_0)
    Initial(DEP(b=None), DEP_0)
    Initial(CPP(b=None, loc='C'), CPP_0)

def rec_events_lig_EGF():
    """ Receptor events involving the ligand EGF."""
    bind_table([[                                                                                                                           EGF(st='M')],
                [erbb(ty='1', bl=None, bd=None, b=None, st='U', loc='C'),                                                                   (par['EGF_bind_ErbB1'])],
                [erbb(ty='3', b=None, bd=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),    None],
                [erbb(ty='4', b=None, bd=None, st='U', loc='C'),                                                                            None]],
                'bl', 'b')
    bind_complex(erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C') % erbb(bl=None, ty='1', bd=1, b=None, st='U', loc='C'), 'bl', EGF(st='M', b=None), 'b', par['EGF_bind_ErbB1d'], m1=erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C'))
    
    bind_complex(erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C') % erbb(bl=ANY, ty='1', bd=1, b=None, st='U', loc='C'), 'bl', EGF(st='M', b=None), 'b', par['EGF_bind_ErbB1d'], m1=erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C'))
    
    # EGF binding/unbinding from endosomal receptors (consistent with Chen/Sorger model, only uncomplexed ErbB1 can bind/release EGF:
    Rule('EGFE_bind_ErbBE',
         erbb(ty='1', bd=None, b=None, st='U', loc='E') + EGF(st='E', b=None) <>
         erbb(ty='1', bd=None, b=None, st='U', loc='E') % EGF(st='M', b=None),
         *par['EGFE_bind_ErbBE'])
    
    # Rate 4: degradation of EGF
    degrade(EGF(b=None, st='E'), par['kdeg_4'])
    
    alias_model_components()
    
def rec_events_lig_HRG():
    """ Receptor events involving the ligand heregulin."""
    bind_table([[                                                                                                                           HRG],
                [erbb(ty='1', bl=None, bd=None, b=None, st='U', loc='C'),                                                                   None],
                [erbb(ty='3', b=None, bd=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),    (par['HRG_bind_ErbB3'])],
                [erbb(ty='4', b=None, bd=None, st='U', loc='C'),                                                                            (par['HRG_bind_ErbB4'])]],
                'bl', 'b')
      
def rec_events():
    """ Describe receptor-level events (except ligand interactions) here. 
    """
    
    # Parameter definitions
    # =====================
    # Alias model components for names in present namespace
    alias_model_components()
    # EGF / HRG binding to receptors -- see specific ligand condition above
    # EGF / HRG receptor binding rates obtained from Chen et al Jacobian files
    
    # ErbB dimerization
    # Dimerization rates obtained from Chen et al Jacobian files
    # erbb1 is not required to contain a ligand in order to dimerize (3 and 4 are)
    erbb1 = erbb(ty='1', bl=None, b=None, st='U', loc='C')
    erbb1Lig = erbb(ty='1', bl=ANY, b=None, st='U', loc='C')
    erbb2Lig = erbb(ty='2', b=None, st='U', loc='C')
    erbb3Lig = erbb(ty='3', bl=ANY, b=None, st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None)
    erbb4Lig = erbb(ty='4', bl=ANY, b=None, st='U', loc='C')
    bind_table([[                          erbb1,                      erbb1Lig,                     erbb2Lig,                    erbb3Lig, erbb4Lig],
                [erbb1,                    (par['ErbB1_bind_ErbB1']),  None,                         None,                        None,     None],
                [erbb1Lig,                 (par['ErbB1_bind_ErbB1L']), (par['ErbB1L_bind_ErbB1L']),  None,                        None,     None],
                [erbb2Lig,                 (par['ErbB1_bind_ErbB2']),  (par['ErbB1L_bind_ErbB2']),   (par['ErbB2_bind_ErbB2']),   None,     None],
                [erbb3Lig,                 (par['ErbB1_bind_ErbB3']),  (par['ErbB1L_bind_ErbB3']),   (par['ErbB2_bind_ErbB3']),   None,     None],
                [erbb4Lig,                 (par['ErbB1_bind_ErbB4']),  (par['ErbB1L_bind_ErbB4']),   (par['ErbB2_bind_ErbB4']),   None,     None]],
                'bd', 'bd')

    # ATP binding: ATP only binds to dimers
    # ATP binding rates obtained from Chen et al (Supplementary)
    # include DEP binding here since they both bind to the same site
    for i in ['3', '4']:
        bind_complex(erbb(ty='2', st='U', loc='C', b=None, bd=1) % erbb(ty=i, st='U', loc='C', bl=ANY, b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', ATP(b=None), 'b', par['ErbB2'+i+'_bind_ATP'], m1=erbb(ty='2', st='U', loc='C', b=None, bd=1))
    
    bind_complex(erbb(ty='1', st='U', loc='C', bl=ANY, b=None, bd=1) % erbb(st='U', loc='C', b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', ATP(b=None), 'b', par['ErbB1_bind_ATP'], m1=erbb(ty='1', st='U', loc='C', bl=ANY, b=None, bd=1))

    for i in ['1', '2', '4']:
        bind_complex(erbb(ty=i, st='P', loc='C', b=None, bd=1) % erbb(st='P', loc='C', b=None, bd=1, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', DEP(b=None), 'b', par['ErbBP'+i+'_bind_DEP'], m1=erbb(ty=i, st='P', loc='C', b=None, bd=1))

    # Cross phosphorylation: only erbb1, 2, and 4 have ATP, and they can cross-phosphorylate any other receptor
    # erbb2:erbb2 pairs only happen by dissociation of phosphorylated monomers
    # kcat phosphorylation obtained from Chen et al Table I pg. 5

    # Both dimers become phosphorylated/dephosphorylated to agree better w Chen/Schoeberl model
    # Receptor Dephosphorylation
    # DEPHOSPHORYLATION: 
    #  * Density enhanced phosphatase1 (DEP1) dephosphorylates ERB1 (at the cell-membrane)
    #  * Protein Tyrosine Phosphatase1b (PTP1b) dephosphorylates all RTKs (at the endo-membrane) 
    #  Berset, TA, Hoier, EF, Hajnal, A: Genes Dev. 19:1328-1340 (2005)
    #  Haj, FG, Verver, PJ, Squire, A, Neel, BG, Bastiaens, PI: Science 295:1708-1711 (2002)

    for i in ['1','2','4']:
        for j in ['1','2','3','4']:
            Rule("cross_phospho_"+i+"_"+j,
                 ATP(b=1) % erbb(ty=i, b=1,    bd=2, st='U') % erbb(ty=j, bd=2, b=None, st='U', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) >>
                 ADP()    + erbb(ty=i, b=None, bd=2, st='P') % erbb(ty=j, bd=2, b=None, st='P', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
                 Parameter("kcp"+i+j, par['ATP_phos_ErbB']))
            Rule("cross_DEphospho_"+i+"_"+j,
                 DEP(b=1)   %  erbb(ty=i, b=1,    bd=2, st='P') % erbb(ty=j, bd=2, b=None, st='P', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) >>
                 DEP(b=None) + erbb(ty=i, b=None, bd=2, st='U') % erbb(ty=j, bd=2, b=None, st='U', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
                 Parameter("kcd"+i+j, par['DEP_dephos_ErbB']))


    #ErbB2 lateral signaling - ErbB2P-ErbB2P dimers can only form by the dissociation of ligand-containing, phosphorylated dimers containing ErbB2.  The monomeric activated ErbB2 can then bind and activate other monomers (ErbB1, 3, or 4 -- allows EGF signal to be transmitted by ErbB2/ErbB3 and ErbB2/ErbB4 complexes, even though 3 and 4 can't bind EGF) or bind another phosphorylated ErbB2 to form an active complex (that still requires an EGF signal to get started)

    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='3', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='4', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])

    bind(erbb(bl=None, ty='2', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), 'bd', erbb(bl=None, ty='2', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), 'bd', par['ErbB2P_ErbBXP_bind'])
    bind(erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='3', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'bd', par['ErbB2P_ErbBXP_bind'])
    bind(erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='4', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB2P_ErbBXP_bind'])
    
    for i in ['3', '4']:
        Rule('ErbB2_lateralsignal_'+i,
             erbb(ty='2', bd=None, st='P', b=None, loc='C') + erbb(ty=i, bd=None, st='U', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) >>
             erbb(ty='2', bd=1, st='P', b=None, loc='C') % erbb(ty=i, bd=1, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
             Parameter('ErbB2_lateralsignal_k'+i, par['ErbB2_lateralsignal']))

    # Receptor internalization
    # This internalizes receptors (with/without complexes) after binding to CPP (coated pit protein) as well as without CPP (first 2 sets of rules).  Only ErbB1/ErbB1 dimers and ErbB1/ErbBX:GRB2:SOS:RAS-GTP can bind and be internalized by CPP.
    # The Chen/Sorger model implements different internalization rates for different receptor combinations/complexes:
    # Rate 1: The first four rules are to internalize all ErbB1/ErbB1 complexes in the MAPK pathway (i.e. ErbB1/ErbB1:GRB2:SOS:RAS-GDP and ErbB1/ErbB1:SHC-P:GRB2:SOS:RAS-GDP and all intermediates in their formation.  The last one internalizes undimerized ErbB1.)
    Rule("rec_intern_1",
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % GRB2(bgap=2, bgab1=None, b=None, bcpp=None) <> erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % GRB2(bgap=2, bgab1=None, b=None, bcpp=None),
         *par['kint_no_cPP_1'])

    Rule("rec_intern_2",
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % SHC(bgap=2, batp=None, bgrb=None) <>
         erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % SHC(bgap=2, batp=None, bgrb=None),
         *par['kint_no_cPP_1'])

    Rule('rec_intern_4a',
         erbb(bl=None, bd=1, loc='C', cpp='N', ty='1', st='P', b=None) % erbb(bl=None, bd=1, loc='C', cpp='N', ty='1', st='P', b=None) <>
         erbb(bl=None, bd=1, loc='E', cpp='N', ty='1', st='P', b=None) % erbb(bl=None, bd=1, loc='E', cpp='N', ty='1', st='P', b=None),
         *par['kint_no_cPP_1'])
    
    Rule('rec_intern_4b',
         erbb(bl=ANY, bd=1, loc='C', cpp='N', ty='1', st='P', b=None) % erbb(bl=ANY, bd=1, loc='C', cpp='N', ty='1', st='P', b=None) <>
         erbb(bl=ANY, bd=1, loc='E', cpp='N', ty='1', st='P', b=None) % erbb(bl=ANY, bd=1, loc='E', cpp='N', ty='1', st='P', b=None),
         *par['kint_no_cPP_1'])
    
    Rule('rec_intern_4c',
         erbb(bl=ANY, bd=1, loc='C', cpp='N', ty='1', st='P', b=None) % erbb(bl=None, bd=1, loc='C', cpp='N', ty='1', st='P', b=None) <>
         erbb(bl=ANY, bd=1, loc='E', cpp='N', ty='1', st='P', b=None) % erbb(bl=None, bd=1, loc='E', cpp='N', ty='1', st='P', b=None),
         *par['kint_no_cPP_1'])

    Rule('rec_intern_5',
         erbb(bd=None, loc='C', cpp='N', ty='1', b=None) <>
         erbb(bd=None, loc='E', cpp='N', ty='1', b=None),
         *par['kint_no_cPP_1'])

    # Rate 2: Set to 0 in Chen/Sorger files and not implemented.  Would have internalized single ErbB2, 3, and 4, as well as ErbB2/3,4:SHC complexes (phos/unphos).
    # Rate 3: These rules internalize ErbB1/ErbBX dimers, ErbB2/ErbB2:SHC complexes and intermediates, and ErbB2/ErbB3 and ErbB2/ErbB4 dimers.
    for i in ['2', '3', '4']:
        Rule('rec_intern_6_'+i,
             erbb(bd=1, loc='C', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='C', cpp='N', b=None, ty=i, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) <>
             erbb(bd=1, loc='E', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='E', cpp='N', b=None, ty=i, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
            *par['kint_no_cPP_2'])

    Rule('rec_intern_9',
         erbb(bd=1, loc='C', cpp='N', ty='2') % erbb(bd=1, loc='C', cpp='N', ty='2') % SHC(bgap=2, batp=None, bgrb=None) <>
         erbb(bd=1, loc='E', cpp='N', ty='2') % erbb(bd=1, loc='E', cpp='N', ty='2') % SHC(bgap=2, batp=None, bgrb=None),
         *par['kint_no_cPP_2'])

    for i in ['2', '3', '4']:
        Rule('rec_intern_10_'+i,
             erbb(bl=None, bd=1, loc='C', cpp='N', ty='2', st='P', b=None) % erbb(bl=None, bd=1, loc='C', cpp='N', b=None, ty=i, st='P', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) <>
             erbb(bl=None, bd=1, loc='E', cpp='N', ty='2', st='P', b=None) % erbb(bl=None, bd=1, loc='E', cpp='N', b=None, ty=i, st='P', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),
            *par['kint_no_cPP_2'])
    
    #Modification to original Chen 2007 model -- Internalization of dimers in the PI3K-Akt pathway.
    #Three rates: one for ErbB1/ErbBX dimers and one for ErbB2/ErbBX dimers with Gab1 and any other downstream proteins bound, and one for ErbB2/ErbB3 dimers with PI3K directly bound.  ErbB1/ErbB2 dimers are under the first rate.
    for i in ['1', '2', '3', '4']:
        Rule('Akt_path_intern_1Xdimers'+i,
            erbb(bd=1, loc='C', cpp='N', ty='1', st='P', b=2) % erbb(bd=1, loc='C', cpp='N', ty=i, st='P', b=None) % GRB2(bgap=2, bgab1=3) % GAB1(bgrb2=3) <>
            erbb(bd=1, loc='E', cpp='N', ty='1', st='P', b=2) % erbb(bd=1, loc='E', cpp='N', ty=i, st='P', b=None) % GRB2(bgap=2, bgab1=3) % GAB1(bgrb2=3),
            *par['Akt_path_intern_1Xdimers'])
    
    for i in ['2', '3', '4']:
        Rule('Akt_path_intern_2Xdimers'+i,
            erbb(bd=1, loc='C', cpp='N', ty='2', st='P', b=2) % erbb(bd=1, loc='C', cpp='N', ty=i, st='P', b=None) % GRB2(bgap=2, bgab1=3) % GAB1(bgrb2=3) <>
            erbb(bd=1, loc='E', cpp='N', ty='2', st='P', b=2) % erbb(bd=1, loc='E', cpp='N', ty=i, st='P', b=None) % GRB2(bgap=2, bgab1=3) % GAB1(bgrb2=3),
            *par['Akt_path_intern_2Xdimers'])
     
    Rule('Akt_path_intern_23_PI3Kdimers',
          erbb(bd=1, loc='C', cpp='N', ty='2', st='P', b=2) % erbb(bd=1, loc='C', cpp='N', ty=3, st='P', b=None) % PI3K(berb=ANY) <>
          erbb(bd=1, loc='E', cpp='N', ty='2', st='P', b=2) % erbb(bd=1, loc='E', cpp='N', ty=3, st='P', b=None) % PI3K(berb=ANY),
            *par['Akt_path_intern_23_PI3Kdimers'])
        
    # CPP bound to receptors can catalyze their internalization (when they are bound to any complex containing GRB2, except GAB1 complex):
    # Binding to CPP and internalization rates are conflated in order to better match Chen-Sorger model.
    Rule('CPP_bind_GAP_GRB2',
         CPP(loc='C', b=None) + erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(ty='1', bd=1, loc='C', cpp='N') % GRB2(bgap=2, bgab1=None, b=None, bcpp=None) <>
         erbb(ty='1', bd=1, loc='E', cpp='Y') % erbb(ty='1', bd=1, loc='E', cpp='Y') % GRB2(bgap=2, bcpp=3, bgab1=None, b=None) % CPP(loc='E', b=3),
         *par['CPP_bind_ErbB1dimers'])

    Rule('CPP_bind_SHC_GRB2',
         erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(ty='1', bd=1, loc='C', cpp='N') % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None, bgap=None) + CPP(loc='C', b=None) <>
         erbb(ty='1', bd=1, loc='E', cpp='Y') % erbb(ty='1', bd=1, loc='E', cpp='Y') % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None, bgap=None) % CPP(loc='E', b=4),
         *par['CPP_bind_ErbB1dimers'])

    for i in ['1', '2', '3', '4']:
        Rule('CPP_bind_ErbB1_RASGTP_complex_'+i,
        erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(bd=1, ty=i, loc='C', cpp='N') % GRB2(bgap=2, bgab1=None, b=None, bcpp=None, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') + CPP(loc='C', b=None) <>
        erbb(ty='1', bd=1, loc='E', cpp='N') % erbb(bd=1, ty=i, loc='E', cpp='N') % GRB2(bgap=2, bgab1=None, b=None, bcpp=5, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') % CPP(loc='E', b=5),
        *par['CPP_bind_ErbB1dimers'])

    Rule('CPPE_bind_GAP_GRB2',
         erbb(bd=1, ty='1', loc='E', cpp='N') % erbb(bd=1, ty='1', loc='E', cpp='N') % GRB2(bgap=2, bcpp=None, bgab1=None, b=None) + CPP(loc='E', b=None) <>
         erbb(bd=1, ty='1', loc='E', cpp='Y') % erbb(bd=1, ty='1', loc='E', cpp='Y') % GRB2(bgap=2, bcpp=3, bgab1=None, b=None) % CPP(loc='E', b=3),
         *par['CPPE_bind_ErbB1dimers'])

    Rule('CPPE_bind_SHC_GRB2',
         erbb(bd=1, ty='1', loc='E', cpp='N') % erbb(bd=1, ty='1', loc='E', cpp='N') % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None, bgap=None) + CPP(loc='E', b=None) <>
         erbb(bd=1, ty='1', loc='E', cpp='Y') % erbb(bd=1, ty='1', loc='E', cpp='Y') % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None, bgap=None) % CPP(loc='E', b=4),
         *par['CPPE_bind_ErbB1dimers'])
    
    Rule("CPP_intern",
         CPP(loc='E', b=None) <> CPP(loc='C', b=None),
         *par['CPP_int'])
         
    # Receptor degradation
    # This degrades all receptor combos within an endosome
    # The Chen/Sorger model implements different degradation rates for different species:
    # Rate 1: These rules degrade all 2EGF:ErbB1/ErbB1 complexes in the MAPK pathway (i.e. ErbB1/ErbB1:GRB2:SOS:RAS-GDP and ErbB1/ErbB1:SHC-P:GRB2:SOS:RAS-GDP and all intermediates in their formation.), as well as single ErbB1.  
    degrade(erbb(bd=1, loc='E', ty='1', bl=ANY) % erbb(bd=1, loc='E', ty='1', bl=ANY) % GRB2(bgap=2, bcpp=None, bgab1=None, b=None), par['kdeg_1'])

    degrade(erbb(bd=1, loc='E', ty='1', bl=ANY) % erbb(bd=1, loc='E', ty='1', bl=ANY) % SHC(bgap=2, batp=None), par['kdeg_1'])

    degrade(erbb(bd=1, loc='E', ty='1', bl=ANY, b=None) % erbb(bd=1, loc='E', ty='1', bl=ANY, b=None), par['kdeg_1'])

    degrade(erbb(bd=None, loc='E', ty='1'), par['kdeg_1'])

    # Rate 2: These rules degrade all ErbB1/ErbBX species and all ErbB2/ErbB2 species in MAPK pathway.  Chen/Sorger model also included degradation of single ErbB2, 3, and 4 under this constant, but as these are never internalized by Chen/Sorger rule set, these degradation rxns were ignored.

    # Rate 3: These rules degrade all ErbB2/ErbB3 and all ErbB2/ErbB4 complexes in MAPK pathway.

    # Rate 4 -- EGF degradation -- see ligand condition above

    # Rate 5: Degrades ErbB1/ErbBX, ErbB2/ErbB3, and ErbB2/ErbB4 dimers (when no complex attached).
    for i in ['2', '3', '4']:
        degrade(erbb(bd=1, loc='E', ty='1', b=None) % erbb(bd=1, loc='E', ty=i, b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), par['kdeg_5'])

    for i in ['3', '4']:
        degrade(erbb(bd=1, loc='E', ty='2', b=None) % erbb(bd=1, loc='E', ty=i, b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), par['kdeg_5'])

    #Modification to Chen 2007 model: Added degradation of complexes in the AKT/PI3K pathway.
    degrade(erbb(bd=1, loc='E', ty='1', b=2) % erbb(bd=1, loc='E', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % GRB2(bgap=2, bgab1=3) % GAB1(bgrb2=3), par['kdeg_6'])
    degrade(erbb(bd=1, loc='E', ty='2', b=2) % erbb(bd=1, loc='E', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % GRB2(bgap=2, bgab1=3) % GAB1(bgrb2=3), par['kdeg_7'])
    degrade(erbb(bd=1, loc='E', ty='2', b=None) % erbb(bd=1, loc='E', ty='3') % PI3K(berb=ANY), par['kdeg_8'])

def mapk_monomers():
    Monomer('SHC', ['bgap', 'bgrb', 'batp', 'st'], {'st':['U','P']})
    Monomer('GRB2', ['b', 'bsos', 'bgap', 'bgab1', 'bcpp'])
    Monomer('SOS', ['bgrb', 'bras', 'bERKPP', 'st'], {'st':['U', 'P']})
    Monomer('RAS', ['bsos', 'braf', 'bpi3k', 'st'], {'st':['GDP', 'GTP']})
    Monomer('RAF', ['b', 'st', 'ser295'], {'st':['U', 'P'], 'ser295':['U', 'P']})
    Monomer('PP1', ['b'])
    Monomer('PP2', ['b'])
    Monomer('PP3', ['b'])
    Monomer('MEK', ['b', 'st'], {'st':['U', 'P', 'PP']})
    Monomer('ERK', ['b', 'st'], {'st':['U', 'P', 'PP']})

def mapk_initial():
    # Initial values declared in parameter dictionary for given cell type.
    alias_model_components()

    Initial(SHC(bgap=None, bgrb=None, batp=None, st='U'), SHC_0)
    Initial(GRB2(b=None, bsos=None, bgap=None, bgab1=None, bcpp=None), GRB2_0)
    Initial(RAS(bsos=None, braf=None, bpi3k=None, st='GDP'), RAS_0)
    Initial(RAF(b=None, st='U', ser295='U'), RAF_0)
    Initial(MEK(b=None, st='U'), MEK_0)
    Initial(ERK(b=None, st='U'), ERK_0)
    Initial(PP1(b=None), PP1_0)
    Initial(PP2(b=None), PP2_0)
    Initial(PP3(b=None), PP3_0)
    Initial(GRB2(b=None, bsos=1, bgap=None, bgab1=None, bcpp=None) % SOS(bgrb=1, bras=None, bERKPP=None, st='U'), GRB2_SOS_0)
    Initial(ERK(b=None, st='PP'), ERKPP_0)

    
def mapk_events():

    # =====================
    # Alias model components for names in present namespace
    alias_model_components()

    # SHC binds to ErbB dimers
    # Chen-Sorger model assigns 2 sets of rate constants to different dimer combinations.  The kf is the same variable; two different kr variables are used but are assigned the same values in the Jacobian files.  These have been combined into one set in this model.
    
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),'b', SHC(batp=None, st='U', bgrb=None), 'bgap', par['GAP_bind_SHC'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))
        
    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None),'b', SHC(batp=None, st='U', bgrb=None), 'bgap', par['GAP_bind_SHC'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    #SHC:P binds ErbB dimers
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=None), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=None), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    #SHC:P-GRB2 binds ErbB dimers
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=2), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=2), 'bgap', par['GAP_bind_SHCP'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))

    # Bound and unbound SHC phosphorylation - These are represented by two kf, kr pairs in the Chen-Sorger model:
    
    Rule('SHC_phos',
         erbb(bd=ANY, st='P', b=1) % SHC(bgap=1, bgrb=None, batp=None, st='U') <>
         erbb(bd=ANY, st='P', b=1) % SHC(bgap=1, bgrb=None, batp=None, st='P'),
         *par['SHC_phos'])

    # Forward rate is set to 0 (only dephosphorylation is occurring).  
    Rule('SHC_unbound_phos',
         SHC(bgap=None, bgrb=None, batp=None, st='U') <>
         SHC(bgap=None, bgrb=None, batp=None, st='P'),
         *par['SHC_unbound_phos'])
    
    # GRB2 binds to ErbBdimer-SHC:P with or without SOS:
    bind_complex(SHC(batp=None, st='P', bgrb=None, bgap=ANY), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1), 'b', par['GRB2_SOS_bind_SHCP_GAP'])

    bind(SHC(batp=None, st='P', bgap=ANY), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None), 'b', par['GRB2_bind_GAP'])

    # SHC:P can bind GRB2-SOS without being attached to an ErbB dimer:
    
    bind_complex(SHC(batp=None, st='P', bgrb=None, bgap=None), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1), 'b', par['SHCP_bind_GRB2SOS'])

    # ErbB dimers can bind the free SHC:P-GRB2-SOS complex:
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=3, bcpp=None, b=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=3), 'bgap', par['GAP_bind_SHCP_GRB2_SOS'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=3, bcpp=None, b=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=3), 'bgap', par['GAP_bind_SHCP_GRB2_SOS'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), m2=SHC(batp=None, st='P', bgrb=2, bgap=None))


    # GRB2 and SOS bind/disassociate:
    bind(GRB2(bgap=None, bgab1=None, b=None, bcpp=None), 'bsos', SOS(bras=None, bERKPP=None, st='U'), 'bgrb', par['GRB2_bind_SOS'])

    #Although no free SOS is present initially in Chen Sorger model, GRB2-SOS can disassociate (see above), so these are necessary.
    # SOS binds to ErbBdimer-SHC:P-GRB2  
    bind_complex(SHC(batp=None, st='P', bgrb=1, bgap=ANY) % GRB2(b=1, bgap=None, bgab1=None, bcpp=None), 'bsos', SOS(bras=None, st='U', bERKPP=None), 'bgrb', par['SOS_bind_GAP_SHCP_GRB2'])

    #SOS binds SHC:P-GRB2 without complex
    bind_complex(SOS(bras=None, st='U', bERKPP=None, bgrb=None), 'bgrb', SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=1), 'bsos', par['SOS_bind_SHCP_GRB2'])

    # SOS also binds ErbBdimer-GRB2
    bind_complex(GRB2(bgap=ANY, bgab1=None, b=None, bsos=None, bcpp=None), 'bsos', SOS(bras=None, bgrb=None, bERKPP=None, st='U'), 'bgrb', par['SOS_bind_GAP_GRB2'])

    # ErbBdimer-GRB2-SOS and ErbBdimer-SHC:P-GRB2-SOS can bind either Ras-GDP or Ras-GTP: k1, k1r
    bind_complex(GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GDP', bpi3k=None), 'bsos', par['RASGDP_bind_bound_GRB2_SOS'])

    bind_complex(SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GDP', bpi3k=None), 'bsos', par['RASGDP_bind_bound_GRB2_SOS'])

    bind_complex(GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GTP', bpi3k=None), 'bsos', par['RASGTP_bind_bound_GRB2_SOS'])

    bind_complex(SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U'), 'bras', RAS(braf=None, bsos=None, st='GTP', bpi3k=None), 'bsos', par['RASGTP_bind_bound_GRB2_SOS'])

    #Ras-GDP --> Ras-GTP transition (GDP exchange) occurs at a different (faster) rate when bound to Sos (a guanine exchange factor - GEF) than when unbound
    #Ras-GTP --> Ras-GDP transition (GTP hydrolysis) is also covered by these rules.  A GAP (GTPase activating protein) would theoretically affect this rate.

    Rule('RASGTP_to_GDP_SOS_GRB2_SHCP_complex',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', bpi3k=None) <>
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', bpi3k=None),
         *par['RASGTP_unbind_GRB2_SOS'])

    Rule('RASGTP_to_GDP_SOS_GRB2_GAP_complex',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', bpi3k=None) <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', bpi3k=None),
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
    catalyze(MEK(st='PP'), 'b', ERK(st='U'), 'b', ERK(st='P'),
             (par['MEKPP_ERK']))

    # Deactivation of ERK:P -> ERK by PP3
    catalyze(PP3(), 'b', ERK(st='P'), 'b', ERK(st='U'),
             (par['ERKP_PP3']))

    # Activation of ERK:P -> ERK:P:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='P'), 'b', ERK(st='PP'),
             (par['MEKPP_ERKP']))

    # Deactivation of ERK:P:P -> ERK:P by PP3
    catalyze(PP3(), 'b', ERK(st='PP'), 'b', ERK(st='P'),
             (par['ERKPP_PP3']))

    # Degradation of PP3
    degrade(PP3(b=None), par['PP3_deg'])

def akt_monomers():
    """ This is the akt part of the pathway from the Chen et al. 2009 paper.  Initial rules for all binding reactions were generated and then coded again using macros and higher order macros.  Initial parameters and conditions were taken from Chen et al. 2009 paper and supplementary, but were later modified in order to get the model working correctly.  This pathway follows AKT from its initial state to a phosphorylated and then double phosphorylated state before returning to unphosphorylated AKT.  The model works correctly, but parameters and rates may need to be modified in order to get best fit.  Parameters and rates included are from trial and error of what best fit the model.  
"""
    #This pathway coded by Tim O'Brien.
    Monomer('GAB1', ['bgrb2', 'bshp2', 'bpi3k', 'batp','bERKPP','bPase9t','S'],{'S':['U','P','PP']})
    Monomer('PI3K',['bgab1','bpip', 'bras', 'berb'])
    Monomer('SHP2',['bgab1'])
    Monomer('PIP', ['bakt', 'both', 'S', 'bpi3k_self','bself2'], {'S':['PIP2', 'PIP3']})
    Monomer('PTEN', ['bpip3', 'both'])
    Monomer('SHP', ['bpip3', 'both'])
    Monomer('AKT', ['bpip3', 'both', 'S'], {'S':['U', 'P', 'PP']})
    Monomer('PDK1', ['bakt', 'both'])
    Monomer('PP2A_III', ['bakt', 'both'])

def akt_initial():
    # See parameter dictionary files for given cell type for initial values.
    
    alias_model_components()
    
    # Initial conditions 
    Initial(GAB1(bgrb2=None, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'), GAB1_0)
    Initial(PI3K(bgab1=None, bpip=None, bras=None, berb=None), PI3K_0)
    Initial(SHP2(bgab1=None), SHP2_0)
    Initial(PIP(bakt=None, both=None, S='PIP2', bpi3k_self=None, bself2=None), PIP_0)
    Initial(PTEN(bpip3=None, both=None), PTEN_0)
    Initial(SHP(bpip3=None, both=None), SHP_0)
    Initial(AKT(bpip3=None, both=None, S='U'), AKT_0)
    Initial(PDK1(bakt=None, both=None), PDK1_0)
    Initial(PP2A_III(bakt=None, both=None), PP2A_III_0)
    Initial(AKT(bpip3=None, both=None, S='PP'), AKTPP_0)
    
def akt_events():
    #GRB2 binds ErbBdimer-complex (without requiring SHC bound to complex):
    #Bind GRB2 without SOS already bound (two Chen-Sorger rate constants for different receptor dimers):
    #Rate 1: ErbB1/ErbB1 dimers (endosomal and plasma membrane), ErbB2/ErbB2 dimers (endosomal and plasma membrane), ErbB2/ErbB3 dimers (plasma membrane), and ErbB2/ErbB4 dimers (endosomal membrane):
    
    bind_complex(erbb(bd=1, st='P', ty='1', b=None) % erbb(bd=1, st='P', ty='1', b=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None), 'bgap', par['GRB2_bind_GAP_2'], m1=erbb(bd=1, st='P', ty='1', b=None))
    
    bind_complex(erbb(bd=1, st='P', ty='2', b=None) % erbb(bd=1, st='P', ty='2', b=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None), 'bgap', par['GRB2_bind_GAP_2'], m1=erbb(bd=1, st='P', ty='2', b=None))
    
    bind_complex(erbb(bd=1, st='P', ty='2', loc='C', b=None) % erbb(bd=1, st='P', ty='3', loc='C', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None), 'bgap', par['GRB2_bind_GAP_2'], m1=erbb(bd=1, st='P', ty='2', loc='C', b=None))
    
    bind_complex(erbb(bd=1, st='P', ty='2', loc='E', b=None) % erbb(bd=1, st='P', ty='4', loc='E', b=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None), 'bgap', par['GRB2_bind_GAP_2'], m1=erbb(bd=1, st='P', ty='2', loc='E', b=None))

    #Rate 2: ErbB1/ErbBX, X=2, 3, 4 (endosomal and plasma membrane), ErbB2/ErbB3 dimers (endosomal membrane), and ErbB2/ErbB4 dimers (plasma membrane):
    for i in ['2', '3', '4']:
        bind_complex(erbb(bd=1, st='P', ty='1', b=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None), 'bgap', par['GRB2_bind_GAP'], m1=erbb(bd=1, st='P', ty='1', b=None))

    bind_complex(erbb(bd=1, st='P', ty='2', loc='E', b=None) % erbb(bd=1, st='P', ty='3', loc='E', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None), 'bgap', par['GRB2_bind_GAP'], m1=erbb(bd=1, st='P', ty='2', loc='E', b=None))

    bind_complex(erbb(bd=1, st='P', ty='2', loc='C', b=None) % erbb(bd=1, st='P', ty='4', loc='C', b=None), 'b', GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None), 'bgap', par['GRB2_bind_GAP'], m1=erbb(bd=1, st='P', ty='2', loc='C', b=None))

    #Bind GRB2 to ErbB dimers with SOS already bound (one rate constant set for all dimer combinations):
    for i in ['1', '2', '3', '4']:
        bind_complex(erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=2, bgab1=None, bcpp=None, bgap=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=2), 'bgap', par['GRB2_SOS_bind_GAP'], m1=erbb(ty='1', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))

    for i in ['2', '3', '4']:
        bind_complex(erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) % erbb(ty=i, bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None), 'b', GRB2(b=None, bsos=2, bgab1=None, bcpp=None, bgap=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=2), 'bgap', par['GRB2_SOS_bind_GAP'], m1=erbb(ty='2', bd=1, st='P', b=None, pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None))


    #GAB1 binds ErbB:ErbB-GRB2. 
    
    bind_complex(erbb(bd=1) % erbb(bd=1) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None, bcpp=None), 'bgab1', GAB1(bgrb2=None, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'), 'bgrb2', par['GRB2_bind_GAB1'])
    
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
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_1'])
    
    Rule('ErbB23_PI3K_cat_2',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_2'])
    
    Rule('ErbB23_PI3K_cat_3',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=None, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_3'])
    
    Rule('ErbB23_PI3K_cat_4',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=None, pi3k6=None) + PIP(bakt=None, both=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_4'])
    
    Rule('ErbB23_PI3K_cat_5',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=None) + PIP(bakt=None, both=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=None) + PIP(bakt=None, both=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_5'])
    
    Rule('ErbB23_PI3K_cat_6',
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=ANY) + PIP(bakt=None, both=None, S='PIP2') >>
         erbb(bd=1, ty='2', b=None, st='P') % erbb(bd=1, ty='3', b=None, st='P', pi3k1=ANY, pi3k2=ANY, pi3k3=ANY, pi3k4=ANY, pi3k5=ANY, pi3k6=ANY) + PIP(bakt=None, both=None, S='PIP3'), 
         par['ErbB23_PI3K_cat_6'])
    
    #PI3K bound to complex catalyzes PIP2 -> PIP3
    #Two rate sets for initial binding in Chen/Sorger model:
    #Rate 1: ErbB1/ErbBX dimers:
    
    bind_complex(erbb(bd=ANY, ty='1') % erbb(bd=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None), 'bpip', PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=None), 'bpi3k_self', par['PIP2_bind_PI3K_1'])

    #Rate 2: ErbB2/ErbBX dimers, X=2, 3, 4:
    for i in ['2', '3', '4']:
        bind_complex(erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None), 'bpip', PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=None), 'bpi3k_self', par['PIP2_chain_PI3K'])
    
    #Two catalysis rates in Chen/Sorger model:
    #Rate 1: ErbB2/ErbB3 dimers:
    Rule('PIP2_PI3K_catalysis_1',
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='3') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1) >>
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='3') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', both=None, bakt=None, bself2=None, bpi3k_self=None),
         par['PIP2_self_catalysis'])

    #Rate 2: All other dimers:
    for i in ['1', '2', '3', '4']:
        Rule('PIP2_PI3K_catalysis_2_'+i,
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1) >>
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', both=None, bakt=None, bself2=None, bpi3k_self=None),
             par['PIP2_PI3K_catalysis'])

    for i in ['2', '4']:
        Rule('PIP2_PI3K_catalysis_3_'+i,
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1) >>
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', both=None, bakt=None, bself2=None, bpi3k_self=None),
             par['PIP2_PI3K_catalysis'])
             
    # Setting up the binding reactions necessary for AKT to be phosphorylated and move through the pathway
    bind_table([[                                                 AKT(S='U', both=None),       AKT(S='P', both=None)],
                [PIP(S='PIP3', both=None, bpi3k_self=None),       (par['PIP3_bind_AKT']),     (par['PIP3_bind_AKT'])]],
                'bakt', 'bpip3')
    
    # AKT-PIP3 is phosphorylated by PDK1 to AKTP; PDK1-PIP3 and AKTP are released
    bind(PDK1(both=None), 'bakt', AKT(bpip3=ANY, S='U'), 'both', par['AKT_PIP3_bind_PDK1'])
    
    Rule('PDK1_AKT_catalysis',
         PDK1(both=None, bakt=1) % AKT(bpip3=2, S='U', both=1) % PIP(S='PIP3', both=None, bpi3k_self=None, bakt=2) >>
         PDK1(both=3, bakt=None) % PIP(S='PIP3', both=3, bpi3k_self=None, bakt=None) + AKT(bpip3=None, S='P', both=None),
         par['PDK1_AKT_catalysis'])

    # PIP3 unbinds PDK1
    bind(PIP(S='PIP3', bakt=None, bpi3k_self=None), 'both', PDK1(bakt=None), 'both', par['PIP3_bind_PDK1'])

    # AKTP-PIP3 is phosphorylated by PDK1 to AKTPP
    bind(PDK1(both=None), 'bakt', AKT(bpip3=ANY, S='P'), 'both', par['AKT_PIP3_bind_PDK1'])
    
    Rule('PDK1_AKTP_catalysis',
         PDK1(both=None, bakt=1) % AKT(bpip3=2, S='P', both=1) % PIP(S='PIP3', both=None, bpi3k_self=None, bakt=2) >>
         PDK1(both=3, bakt=None) % PIP(S='PIP3', both=3, bpi3k_self=None, bakt=None) + AKT(bpip3=None, S='PP', both=None),
         par['PDK1_AKTP_catalysis'])

    # AKTP is dephosphorylated by PP2A-III back to AKT
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'both', 'S', 'P', 'U',(par['AKTP_dephos']))
   
    # AKTPP is dephosphorylated by PP2A-III back to AKTP
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'both', 'S', 'PP', 'P',(par['AKTPP_dephos']))

    # PIP3 is dephosphorylated by PTEN to PIP2
    catalyze_state(PTEN, 'bpip3', PIP(bakt=None, bpi3k_self=None), 'both', 'S', 'PIP3', 'PIP2', (par['PIP3_dephos']))

    # PIP3 is dephosphorylated by SHP to PIP2
    catalyze_state(SHP, 'bpip3', PIP(bakt=None, bpi3k_self=None), 'both', 'S', 'PIP3', 'PIP2', (par['PIP3_dephos']))

def crosstalk_monomers():
    Monomer('Pase9t', ['bgab1'])
    alias_model_components()
    
def crosstalk_initial():
    Initial(Pase9t(bgab1=None), Pase9t_0)

def crosstalk_events():
    #ERK:P:P phosphorylates ErbBdimer-GRB2-GAB1:P (making it unable to bind PI3K)
    catalyze_state(ERK(st='PP'), 'b', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None), 'bERKPP', 'S', 'P', 'PP', (par['ERKPP_phos_GAB1P']))

    #ErbBdimer-GRB2-GAB1:P:P is dephosphorylated by Pase9t
    catalyze_state(Pase9t(), 'bgab1', GAB1(bgrb2=ANY), 'bPase9t', 'S', 'PP', 'P', (par['Pase9t_dephos_GAB1PP']))

    #ERK:P:P phosphorylates GRB2-SOS, preventing RAS-GDP->RAS-GTP conversion
    #To conform with Chen/Sorger model, this only effects ErbB1/ErbB1 dimers containing SOS and free SOS, and phosphorylated SOS can only bind ErbB1/ErbB1 complexes, not free GRB2:
    catalyze_state(ERK(st='PP'), 'b', SOS(bgrb=None, bras=None), 'bERKPP', 'st', 'U', 'P', (par['ERKPP_phos_SOS']))

    bind_complex(erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=None, st='U'), 'bERKPP', ERK(st='PP', b=None), 'b', par['ERKPP_phos_SOS'][0:2])

    bind_complex(erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='U', bgrb=ANY), 'bERKPP', ERK(st='PP', b=None), 'b', par['ERKPP_phos_SOS'][0:2])

    Rule('ERKPP_phos_SOS_1',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=1, st='U') % ERK(st='PP', b=1) >>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=None, st='P') + ERK(st='PP', b=None),
         par['ERKPP_phos_SOS'][2])

    Rule('ERKPP_phos_SOS_2',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=1, st='U', bgrb=ANY) % ERK(st='PP', b=1) >>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='P', bgrb=ANY) + ERK(st='PP', b=None),
         par['ERKPP_phos_SOS'][2])

    bind_complex(erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None), 'bsos', SOS(bras=None, bgrb=None, bERKPP=None, st='P'), 'bgrb', par['SOSP_bind_GRB2'])

    bind_complex(erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=ANY), 'bsos', SOS(bras=None, bERKPP=None, st='P', bgrb=None), 'bgrb', par['SOSP_bind_GRB2'])

    #AKT:PP phosphorylates RAF:P at Ser295, preventing MEK phosphorylation.
    catalyze_state(AKT(S='PP', bpip3=None), 'both', RAF(st='P'), 'b', 'ser295', 'U', 'P', (par['AKTPP_phos_RAFP']))

    #RAS-GTP binds PI3K and activates PI3K catalytic function.
    bind(RAS(bsos=None, braf=None, st='GTP'), 'bpi3k', PI3K(bgab1=None, bpip=None, berb=None), 'bras', par['RASGTP_bind_PI3K'])
    
    catalyze_state(PI3K(bras=ANY, bgab1=None, berb=None), 'bpip', PIP(both=None, bakt=None, bself2=None), 'bpi3k_self', 'S', 'PIP2', 'PIP3', par['RAS_PI3K_cat_PIP'])
