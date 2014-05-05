"""
Overview:
========

PySB implementation of the ErbB related MAPK and AKT signaling
pathways originally published in [Chen2009]_.

This file contains functions that implement the ErbB execution
pathway in three modules:

- Receptor layer events, taking into account all ErbB1-4 interactions with ligand.
- AKT pathway
- MAPK pathway

"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components
# from egfr.shared import * # modified model aliases

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
#FIXME: What is Inh in reaction list?
        
# Monomer declarations
# ====================

def rec_monomers():

    """ Declares the ErbB receptor interactions.
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    Monomer('EGF', ['b', 'st'], {'st':['M', 'E']}) # Epidermal Growth Factor ligand
    Monomer('HRG', ['b']) # Heregulin ligand
    Monomer('erbb', ['bl', 'bd', 'b', 'ty', 'st', 'loc', 'pi3k1', 'pi3k2', 'pi3k3', 'pi3k4', 'pi3k5', 'pi3k6', 'cpp'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E'], 'cpp':['Y', 'N']}) # bl: lig, bd: dimer, b: binding, ty: rec type, st: (U)n(P)hosphorylated, loc: (C)yto 'brane or (E)ndosome 'brane, cpp: No real biophysical meaning; useful model marker for presence of CPP bound downstream.
    Monomer('DEP', ['b'])
    Monomer('ATP', ['b'])
    Monomer('ADP')
    Monomer('CPP', ['b', 'loc'], {'loc':['C', 'E']})

def rec_initial_lig_hEGF():
    Parameter('EGF_0',     3.01e12) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',         0) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    
def rec_initial_lig_lEGF():
    Parameter('EGF_0',      6.02e9) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',         0) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()

def rec_initial_lig_hHRG():
    Parameter('EGF_0',      0) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',       3.01e12) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()

def rec_initial_lig_lHRG():
    Parameter('EGF_0',      0) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',         6.02e9) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    
def rec_initial():
    # # Initial concentrations (except DEP1) for all cell types taken from Chen et al 2009 -- see Jacobian files
    alias_model_components()

    Initial(EGF(b=None, st='M'), EGF_0)
    Initial(HRG(b=None), HRG_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='1', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb1_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='2', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb2_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='3', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb3_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='4', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb4_0)
    Initial(ATP(b=None), ATP_0)
    Initial(DEP(b=None), DEP_0)
    Initial(CPP(b=None, loc='C'), CPP_0)
            
def rec_events():
    """ Describe receptor-level events here. 
    """
    
    # Parameter definitions
    # =====================
    # Alias model components for names in present namespace
    alias_model_components()
    # EGF / HRG binding to receptors
    # EGF / HRG receptor binding rates obtained from Chen et al Jacobian files
    bind_table([[                                                                                    EGF(st='M'),                                   HRG],
                [erbb(ty='1', bl=None, bd=None, b=None, st='U', loc='C'),                            (par['EGF_bind_ErbB1']),                       None],
                [erbb(ty='3', b=None, bd=None, st='U', loc='C'),                                     None,                                          (par['HRG_bind_ErbB3'])],
                [erbb(ty='4', b=None, bd=None, st='U', loc='C'),                                     None,                                          (par['HRG_bind_ErbB4'])]],
                'bl', 'b')

    Rule('EGF_bind_ErbB1dimers_1',
         erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C') % erbb(bl=None, ty='1', bd=1, b=None, st='U', loc='C') + EGF(st='M', b=None) <>
         erbb(ty='1', bl=2, bd=1, b=None, st='U', loc='C') % erbb(bl=None, ty='1', bd=1, b=None, st='U', loc='C') % EGF(st='M', b=2),
         *par['EGF_bind_ErbB1d'])
    
    Rule('EGF_bind_ErbB1dimers_2',
         erbb(ty='1', bl=None, bd=1, b=None, st='U', loc='C') % erbb(bl=ANY, ty='1', bd=1, b=None, st='U', loc='C') + EGF(st='M', b=None) <>
         erbb(ty='1', bl=2, bd=1, b=None, st='U', loc='C') % erbb(bl=ANY, ty='1', bd=1, b=None, st='U', loc='C') % EGF(st='M', b=2),
         *par['EGF_bind_ErbB1d'])
        
    # EGF binding/unbinding from endosomal receptors (consistent with Chen/Sorger model, only uncomplexed ErbB1 can bind/release EGF:
    Rule('EGFE_bind_ErbBE',
         erbb(ty='1', bd=None, b=None, st='U', loc='E') + EGF(st='E', b=None) <>
         erbb(ty='1', bd=None, b=None, st='U', loc='E') % EGF(st='M', b=None),
         *par['EGFE_bind_ErbBE'])
    
    # ErbB dimerization
    # Dimerization rates obtained from Chen et al Jacobian files
    # erbb1 is not required to contain a ligand in order to dimerize (3 and 4 are)
    erbb1 = erbb(ty='1', bl=None, b=None, st='U', loc='C')
    erbb1Lig = erbb(ty='1', bl=ANY, b=None, st='U', loc='C')
    erbb2Lig = erbb(ty='2', b=None, st='U', loc='C')
    erbb3Lig = erbb(ty='3', bl=ANY, b=None, st='U', loc='C')
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
        Rule('ATP_bind_ErbB2'+i,
             erbb(ty='2', st='U', loc='C', b=None, bd=1) % erbb(ty=i, st='U', loc='C', b=None, bd=1) + ATP(b=None) <>
             erbb(ty='2', st='U', loc='C', b=2, bd=1) % erbb(ty=i, st='U', loc='C', b=None, bd=1) % ATP(b=2),
             *par['ErbB2'+i+'_bind_ATP'])

    Rule('ATP_bind_ErbB1',
             erbb(ty='1', st='U', loc='C', bl=ANY, b=None, bd=1) % erbb(st='U', loc='C', b=None, bd=1) + ATP(b=None) <>
             erbb(ty='1', st='U', loc='C', bl=ANY, b=2, bd=1) % erbb(st='U', loc='C', b=None, bd=1) % ATP(b=2),
             *par['ErbB1_bind_ATP'])

    for i in ['1', '2', '4']:
        Rule('DEP_bind_ErbB'+i,
             erbb(ty=i, st='P', loc='C', b=None, bd=1) % erbb(st='P', loc='C', b=None, bd=1) + DEP(b=None) <>
             erbb(ty=i, st='P', loc='C', b=2, bd=1) % erbb(st='P', loc='C', b=None, bd=1) % DEP(b=2),
             *par['ErbBP'+i+'_bind_DEP'])

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
                 ATP(b=1) % erbb(ty=i, b=1,    bd=2, st='U') % erbb(ty=j, bd=2, b=None, st='U') >>
                 ADP()    + erbb(ty=i, b=None, bd=2, st='P') % erbb(ty=j, bd=2, b=None, st='P'),
                 Parameter("kcp"+i+j, par['ATP_phos_ErbB']))
            Rule("cross_DEphospho_"+i+"_"+j,
                 DEP(b=1)   %  erbb(ty=i, b=1,    bd=2, st='P') % erbb(ty=j, bd=2, b=None, st='P') >>
                 DEP(b=None) + erbb(ty=i, b=None, bd=2, st='U') % erbb(ty=j, bd=2, b=None, st='U'),
                 Parameter("kcd"+i+j, par['DEP_dephos_ErbB']))


    #ErbB2 lateral signaling - ErbB2P-ErbB2P dimers can only form by the dissociation of ligand-containing, phosphorylated dimers containing ErbB2.  The monomeric activated ErbB2 can then bind and activate other monomers (ErbB1, 3, or 4 -- allows EGF signal to be transmitted by ErbB2/ErbB3 and ErbB2/ErbB4 complexes, even though 3 and 4 can't bind EGF) or bind another phosphorylated ErbB2 to form an active complex (that still requires an EGF signal to get started)

    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='3', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])
    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='4', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])

    bind(erbb(bl=None, ty='2', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), 'bd', erbb(bl=None, ty='2', bd=None, st='P', b=None, loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), 'bd', par['ErbB2P_ErbBXP_bind'])
    bind(erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='3', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB2P_ErbBXP_bind'])
    bind(erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', erbb(ty='4', bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB2P_ErbBXP_bind'])
    
    for i in ['3', '4']:
        Rule('ErbB2_lateralsignal_'+i,
             erbb(ty='2', bd=None, st='P', b=None, loc='C') + erbb(ty=i, bd=None, st='U', b=None, loc='C') >>
             erbb(ty='2', bd=1, st='P', b=None, loc='C') % erbb(ty=i, bd=1, st='P', b=None, loc='C'),
             Parameter('ErbB2_lateralsignal_k'+i, par['ErbB2_lateralsignal']))


    alias_model_components()
    # Receptor internalization
    # This internalizes receptors (with/without complexes) after binding to CPP (coated pit protein) as well as without CPP (first 2 sets of rules).  Only ErbB1/ErbB1 dimers and ErbB1/ErbBX:GAP:GRB2:SOS:RAS-GTP can bind and be internalized by CPP.
    # The Chen/Sorger model implements different internalization rates for different receptor combinations/complexes:
    # Rate 1: The first four rules are to internalize all ErbB1/ErbB1 complexes in the MAPK pathway (i.e. ErbB1/ErbB1:GAP:GRB2:SOS:RAS-GDP and ErbB1/ErbB1:GAP:SHC-P:GRB2:SOS:RAS-GDP and all intermediates in their formation.  The last one internalizes undimerized ErbB1.)
    Rule("rec_intern_1",
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None) <> erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None),
         *par['kint_no_cPP_1'])

    Rule("rec_intern_2",
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, bgrb=None) <>
         erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, bgrb=None),
         *par['kint_no_cPP_1'])

    Rule('rec_intern_3',
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % GAP(bd=ANY, b=None, bgrb2=None) <>
         erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % GAP(bd=ANY, b=None, bgrb2=None),
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

# Rate 2: Set to 0 in Chen/Sorger files and not implemented.  Would have internalized single ErbB2, 3, and 4, as well as ErbB2/3,4:GAP:SHC complexes (phos/unphos).
    # Rate 3: These rules internalize ErbB1/ErbBX dimers, ErbB2/ErbB2:GAP:SHC complexes and intermediates, and ErbB2/ErbB3 and ErbB2/ErbB4 dimers.
    for i in ['2', '3', '4']:
        Rule('rec_intern_6_'+i,
             erbb(bd=1, loc='C', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='C', cpp='N', b=None, ty=i) <>
             erbb(bd=1, loc='E', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='E', cpp='N', b=None, ty=i),
            *par['kint_no_cPP_2'])

    Rule('rec_intern_8',
         erbb(bd=1, loc='C', cpp='N', ty='2') % erbb(bd=1, loc='C', cpp='N', ty='2') % GAP(bd=ANY, b=None, bgrb2=None) <>
         erbb(bd=1, loc='E', cpp='N', ty='2') % erbb(bd=1, loc='E', cpp='N', ty='2') % GAP(bd=ANY, b=None, bgrb2=None),
         *par['kint_no_cPP_2'])

    Rule('rec_intern_9',
         erbb(bd=1, loc='C', cpp='N', ty='2') % erbb(bd=1, loc='C', cpp='N', ty='2') % GAP(bd=ANY, b=2, bgrb2=None) % SHC(bgap=2, batp=None, bgrb=None) <>
         erbb(bd=1, loc='E', cpp='N', ty='2') % erbb(bd=1, loc='E', cpp='N', ty='2') % GAP(bd=ANY, b=2, bgrb2=None) % SHC(bgap=2, batp=None, bgrb=None),
         *par['kint_no_cPP_2'])

    for i in ['2', '3', '4']:
        Rule('rec_intern_10_'+i,
             erbb(bl=None, bd=1, loc='C', cpp='N', ty='2', st='P', b=None) % erbb(bl=None, bd=1, loc='C', cpp='N', b=None, ty=i, st='P') <>
             erbb(bl=None, bd=1, loc='E', cpp='N', ty='2', st='P', b=None) % erbb(bl=None, bd=1, loc='E', cpp='N', b=None, ty=i, st='P'),
            *par['kint_no_cPP_2'])
        
    # CPP bound to receptors can catalyze their internalization (when they are bound to any complex containing GRB2, except GAB1 complex):
    # Binding to CPP and internalization rates are conflated in order to better match Chen-Sorger model.
    Rule('CPP_bind_GAP_GRB2',
         CPP(loc='C', b=None) + erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(ty='1', bd=1, loc='C', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None) <>
         erbb(ty='1', bd=1, loc='E', cpp='Y') % erbb(ty='1', bd=1, loc='E', cpp='Y') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3, bgab1=None, b=None) % CPP(loc='E', b=3),
         *par['CPP_bind_ErbB1dimers'])

    Rule('CPP_bind_SHC_GRB2',
         erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(ty='1', bd=1, loc='C', cpp='N') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None, bgap=None) + CPP(loc='C', b=None) <>
         erbb(ty='1', bd=1, loc='E', cpp='Y') % erbb(ty='1', bd=1, loc='E', cpp='Y') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None, bgap=None) % CPP(loc='E', b=4),
         *par['CPP_bind_ErbB1dimers'])

    for i in ['1', '2', '3', '4']:
        Rule('CPP_bind_ErbB1_RASGTP_complex_'+i,
        erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(bd=1, ty=i, loc='C', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') + CPP(loc='C', b=None) <>
        erbb(ty='1', bd=1, loc='E', cpp='N') % erbb(bd=1, ty=i, loc='E', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=5, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') % CPP(loc='E', b=5),
        *par['CPP_bind_ErbB1dimers'])

    Rule('CPPE_bind_GAP_GRB2',
         erbb(bd=1, ty='1', loc='E', cpp='N') % erbb(bd=1, ty='1', loc='E', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=None, bgab1=None, b=None) + CPP(loc='E', b=None) <>
         erbb(bd=1, ty='1', loc='E', cpp='Y') % erbb(bd=1, ty='1', loc='E', cpp='Y') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3, bgab1=None, b=None) % CPP(loc='E', b=3),
         *par['CPPE_bind_ErbB1dimers'])

    Rule('CPPE_bind_SHC_GRB2',
         erbb(bd=1, ty='1', loc='E', cpp='N') % erbb(bd=1, ty='1', loc='E', cpp='N') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None, bgap=None) + CPP(loc='E', b=None) <>
         erbb(bd=1, ty='1', loc='E', cpp='Y') % erbb(bd=1, ty='1', loc='E', cpp='Y') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None, bgap=None) % CPP(loc='E', b=4),
         *par['CPPE_bind_ErbB1dimers'])
    
    Rule("CPP_intern",
         CPP(loc='E', b=None) <> CPP(loc='C', b=None),
         *par['CPP_int'])
         
    # Receptor degradation
    # This degrades all receptor combos within an endosome
    # The Chen/Sorger model implements different degradation rates for different species:
    # Rate 1: These rules degrade all 2EGF:ErbB1/ErbB1 complexes in the MAPK pathway (i.e. ErbB1/ErbB1:GAP:GRB2:SOS:RAS-GDP and ErbB1/ErbB1:GAP:SHC-P:GRB2:SOS:RAS-GDP and all intermediates in their formation.), as well as single ErbB1.  
    degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3) % erbb(bd=1, loc='E', ty='1', bl=4) % GAP(bd=ANY, bgrb2=2, b=None) % GRB2(bgap=2, bcpp=None, bgab1=None, b=None), par['kdeg_1'])

    degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3) % erbb(bd=1, loc='E', ty='1', bl=4) % GAP(bd=ANY, b=2, bgrb2=None) % SHC(bgap=2, batp=None), par['kdeg_1'])

    degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3) % erbb(bd=1, loc='E', ty='1', bl=4) % GAP(bd=ANY, b=None, bgrb2=None), par['kdeg_1'])

    degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3, b=None) % erbb(bd=1, loc='E', ty='1', bl=4, b=None), par['kdeg_1'])

    degrade(erbb(bd=None, loc='E', ty='1'), par['kdeg_1'])

    # Rate 2: These rules degrade all ErbB1/ErbBX species and all ErbB2/ErbB2 species in MAPK pathway.  Chen/Sorger model also included degradation of single ErbB2, 3, and 4 under this constant, but as these are never internalized by Chen/Sorger rule set, these degradation rxns were ignored.
    for i in ['2', '3', '4']:
        degrade(erbb(bd=1, loc='E', ty='1') % erbb(bd=1, loc='E', ty=i) % GAP(bd=ANY), par['kdeg_2'])

    degrade(erbb(bd=1, loc='E', ty='2') % erbb(bd=1, loc='E', ty='2') % GAP(bd=ANY), par['kdeg_2'])

    # Rate 3: These rules degrade all ErbB2/ErbB3 and all ErbB2/ErbB4 complexes in MAPK pathway.
    for i in ['3', '4']:
        degrade(erbb(bd=1, loc='E', ty='2') % erbb(bd=1, loc='E', ty=i) % GAP(bd=ANY), par['kdeg_3'])

    # Rate 4: degradation of EGF
    degrade(EGF(b=None, st='E'), par['kdeg_4'])

    # Rate 5: Degrades ErbB1/ErbBX, ErbB2/ErbB3, and ErbB2/ErbB4 dimers (when no complex attached).
    for i in ['2', '3', '4']:
        degrade(erbb(bd=1, loc='E', ty='1', b=None) % erbb(bd=1, loc='E', ty=i, b=None), par['kdeg_5'])

    for i in ['3', '4']:
        degrade(erbb(bd=1, loc='E', ty='2', b=None) % erbb(bd=1, loc='E', ty=i, b=None), par['kdeg_5'])

def mapk_monomers():
    Monomer('GAP', ['bd', 'b', 'bgrb2'])
    Monomer('SHC', ['bgap', 'bgrb', 'batp', 'st'], {'st':['U','P']})
    # Monomer('SHCPase', ['b'])
    Monomer('GRB2', ['b', 'bsos', 'bgap', 'bgab1', 'bcpp'])
    Monomer('SOS', ['bgrb', 'bras', 'bERKPP', 'st'], {'st':['U', 'P']})
    Monomer('RAS', ['bsos', 'braf', 'bpi3k', 'st', 'act'], {'st':['GDP', 'GTP'], 'act':['N', 'Y']})
    Monomer('RAF', ['b', 'st', 'ser295'], {'st':['U', 'P'], 'ser295':['U', 'P']})
    Monomer('PP1', ['b'])
    Monomer('PP2', ['b'])
    Monomer('PP3', ['b'])
    Monomer('MEK', ['b', 'st'], {'st':['U', 'P', 'PP']})
    Monomer('ERK', ['b', 'st'], {'st':['U', 'P', 'PP']})

def mapk_initial():
    # Initial values declared in parameter dictionary for given cell type.
    alias_model_components()

    Initial(GAP(bd=None, b=None, bgrb2=None), GAP_0)
    Initial(SHC(bgap=None, bgrb=None, batp=None, st='U'), SHC_0)
    Initial(GRB2(b=None, bsos=None, bgap=None, bgab1=None, bcpp=None), GRB2_0)
    #Initial(SOS(bgrb=None, bras=None, bERKPP=None, st='U'), SOS_0)
    Initial(RAS(bsos=None, braf=None, bpi3k=None, st='GDP', act='N'), RAS_0)
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

    # GAP binds to phosphorylated dimers
    # in the present we use MatchOnce to insure correct representation of the binding
    # similar to Chen et al
    # Chen/Sorger model implements 2 rate constants for this rxn, one for ErbB1/ErbB1, ErbB2/ErbB2, ErbB2/ErbB3, and ErbB2/ErbB4 dimers, and one for ErbB1/ErbBX dimers.
    # Rate 1: ErbB1/ErbB1, ErbB2/ErbB2, ErbB2/ErbB3, and ErbB2/ErbB4 dimers.
    Rule("GAP_binding_1",
         MatchOnce(erbb(bd=1, b=None, st='P', ty='1') % erbb(bd=1, b=None, st='P', ty='1')) + GAP(bd=None, b=None, bgrb2=None) <>
         MatchOnce(erbb(bd=1, b=2,    st='P', ty='1') % erbb(bd=1, b=None, st='P', ty='1') % GAP(bd=2, b=None, bgrb2=None)),
         *par['ErbB_bind_GAP_1'])

    for i in ['2', '3', '4']:
        Rule('GAP_binding_2'+i,
             MatchOnce(erbb(bd=1, b=None, st='P', ty='2') % erbb(bd=1, b=None, st='P', ty=i)) + GAP(bd=None, b=None, bgrb2=None) <>
             MatchOnce(erbb(bd=1, b=2, st='P', ty='2') % erbb(bd=1, b=None, st='P', ty=i) % GAP(bd=2, b=None, bgrb2=None)),
             *par['ErbB_bind_GAP_1'])

    # Rate 2: ErbB1/ErbBX, X=2, 3, 4  Note: In Chen/Sorger rxn list, plasma membrane ErbB1/ErbB2 dimers are assigned Rate 1 (above); however the other 5 ErbB1/ErbBX combinations (plasma and endosomal membranes) are assigned Rate 2.  ErbB1/ErbB2 was assigned the latter in this model under the assumption that this was accidental.
    for i in ['2', '3', '4']:
        Rule('GAP_binding_3'+i,
             MatchOnce(erbb(bd=1, b=None, st='P', ty='1') % erbb(bd=1, b=None, st='P', ty=i)) + GAP(bd=None, b=None, bgrb2=None) <>
             MatchOnce(erbb(bd=1, b=2, st='P', ty='1') % erbb(bd=1, b=None, st='P', ty=i) % GAP(bd=2, b=None, bgrb2=None)),
             *par['ErbB_bind_GAP_2'])

        Rule('GAP_binding_3_2'+i,
             MatchOnce(erbb(bd=1, b=None, st='P', ty=i) % erbb(bd=1, b=None, st='P', ty='1')) + GAP(bd=None, b=None, bgrb2=None) <>
             MatchOnce(erbb(bd=1, b=2, st='P', ty=i) % erbb(bd=1, b=None, st='P', ty='1') % GAP(bd=2, b=None, bgrb2=None)),
             *par['ErbB_bind_GAP_2'])
    
    # SHC binds to GAP-complex
    # Chen-Sorger model assigns 2 sets of rate constants to different dimer combinations.  The kf is the same variable; two different kr variables are used but are assigned the same values in the Jacobian files.  These have been combined into one set in this model.
    
    bind(GAP(bd=ANY, bgrb2=None), 'b', SHC(batp=None, st='U', bgrb=None), 'bgap', par['GAP_bind_SHC'])

    #SHC:P binds GAP
    bind(GAP(bd=ANY, bgrb2=None), 'b', SHC(batp=None, st='P', bgrb=None), 'bgap', par['GAP_bind_SHCP'])

    #SHC:P-GRB2 binds GAP
    Rule('GAP_bind_SHCP_GRB2',
         GAP(bd=ANY, b=None, bgrb2=None) + SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=1) <>
         GAP(bd=ANY, b=2, bgrb2=None) % SHC(batp=None, st='P', bgrb=1, bgap=2) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=1),
         *par['GAP_bind_SHCP'])

    # Bound and unbound SHC phosphorylation - These are represented by two kf, kr pairs in the Chen-Sorger model:
    Rule('SHC_phos',
         GAP(bd=ANY, b=1, bgrb2=None) % SHC(bgap=1, bgrb=None, batp=None, st='U') <>
         GAP(bd=ANY, b=1, bgrb2=None) % SHC(bgap=1, bgrb=None, batp=None, st='P'),
         *par['SHC_phos'])

    # Forward rate is set to 0 (only dephosphorylation is occurring).  
    Rule('SHC_unbound_phos',
         SHC(bgap=None, bgrb=None, batp=None, st='U') <>
         SHC(bgap=None, bgrb=None, batp=None, st='P'),
         *par['SHC_unbound_phos'])
    
    # The two rules below can be used if a binding kf,kr and a kc are desired:
    # Rule("Shc_bind_ATP",
    #      GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=None, st='U') + ATP(b=None) <>
    #      GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=2, st='U') % ATP(b=2),
    #      Parameter("ShcATPf",KF), Parameter("ShcATPr",KR))
    
    # Rule("Shc_phos",
    #      GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=2, st='U') % ATP(b=2) >>
    #      GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=None, st='P') + ADP(),
    #      Parameter("ShcPhosc", KCP))
    
    # GRB2 binds to GAP-SHC:P with or without SOS:
    Rule('GRB2_bind_GAP_SHCP_1',
         SHC(batp=None, st='P', bgrb=None, bgap=ANY) + GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) <>
         SHC(batp=None, st='P', bgrb=2, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         *par['GRB2_SOS_bind_SHCP_GAP'])

    bind(SHC(batp=None, st='P', bgap=ANY), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None), 'b', par['GRB2_bind_GAP'])
    
    # Can use this simpler version if GRB2-SOS complex isn't present alone:
    #    bind(SHC(batp=None, st='P'), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None), 'b', [KF, KR])

    # SHC:P can bind GRB2-SOS without being attached to GAP:
    Rule('SHCP_bind_GRB2SOS',
         SHC(batp=None, st='P', bgrb=None, bgap=None) + GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) <>
         SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         *par['SHCP_bind_GRB2SOS'])

    # GAP can bind the free SHC:P-GRB2-SOS complex:
    Rule('GAP_bind_SHCP_GRB2_SOS',
         GAP(bd=ANY, b=None, bgrb2=None) + SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=2, bcpp=None, b=1) % SOS(bras=None, bERKPP=None, st='U', bgrb=2) <>
         GAP(bd=ANY, b=3, bgrb2=None) % SHC(batp=None, st='P', bgrb=1, bgap=3) % GRB2(bgap=None, bgab1=None, bsos=2, bcpp=None, b=1) % SOS(bras=None, bERKPP=None, st='U', bgrb=2),
         *par['GAP_bind_SHCP_GRB2_SOS'])

    # GRB2 and SOS bind/disassociate:
    bind(GRB2(bgap=None, bgab1=None, b=None, bcpp=None), 'bsos', SOS(bras=None, bERKPP=None, st='U'), 'bgrb', par['GRB2_bind_SOS'])

    #Although no free SOS is present initially in Chen Sorger model, GRB2-SOS can disassociate (see above), so these are necessary.
    # SOS binds to GAP-SHC:P-GRB2  
    bind(GRB2(b=ANY, bgap=None, bgab1=None, bcpp=None), 'bsos', SOS(bras=None, st='U', bERKPP=None), 'bgrb', par['SOS_bind_GAP_SHCP_GRB2'])

    #SOS binds SHC:P-GRB2 without complex
    Rule('SOS_bind_SHCP_GRB2',
         SOS(bras=None, st='U', bERKPP=None, bgrb=None) + SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=1) <>
         SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=2, bcpp=None, b=1) % SOS(bras=None, st='U', bERKPP=None, bgrb=2),
         *par['SOS_bind_SHCP_GRB2'])

    # SOS also binds GAP-GRB2
    Rule("GAP_GRB2_bind_SOS",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=None, bcpp=None) + SOS(bras=None, bgrb=None, bERKPP=None, st='U') <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, st='U', bERKPP=None),
         *par['SOS_bind_GAP_GRB2'])

    # GAP-GRB2-SOS and GAP-SHC:P-GRB2-SOS catalyze RAS-GDP->RAS-GTP:
    Rule("GAP_GRB2_SOS_bind_RASGDP",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None) <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGDP_bind_bound_GRB2_SOS'])

    Rule("GAP_SHCP_GRB2_SOS_bind_RASGDP",
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None) <>
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGDP_bind_bound_GRB2_SOS'])

    # Instead of a one-way catalytic process, the Chen-Sorger model implements this as a bidirectional process, as below:
    Rule('GAP_GRB2_SOS_bind_RASGTP',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='N', bpi3k=None) <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_bind_bound_GRB2_SOS'])

    Rule('GAP_SHCP_GRB2_SOS_bind_RASGTP',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='N', bpi3k=None) <>
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_bind_bound_GRB2_SOS'])

    # If a catalytic process is desired instead, use these rules:
    # Rule("GAP_GRB2_SOS_catRAS",
    #      GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP') >>
    #      GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP'),
    #      Parameter("GAP_GRB2_SOS_catRASc", KCP))

    # Rule("GAP_SHCP_GRB2_SOS_catRAS",
    #      GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP') >>
    #      GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP'),
    #      Parameter("GAP_SHCP_GRB2_SOS_catRASc", KCP))

    #FIXME: Need to add re-binding of RAS-GTP to complexes, deal with RAS active/unactive? and separate pools of internalized/un-internalized RAS (and RAF, MEK, and ERK).
         
    #can use this simpler implementation of above if GRB2-SOS isn't present on its own as in Chen Sorger model:
    #catalyze_state(SOS(bgrb=ANY, st='U', bERKPP=None), 'bras', RAS(braf=None), 'bsos', 'st', 'GDP', 'GTP', (KF, KR, KCD))

    # Recycling of activated RAS-GTP --> RAS-GDP.  In Chen/Sorger model, activated RAS-GTP is produced upon Raf phosphorylation.
    Rule('RASGTPact_bind_SOS_SHCP_complex',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='Y', bpi3k=None) <>
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=None, st='GTP', act='N', bpi3k=None),
         *par['RASGTPact_bind_bound_GRB2_SOS'])

    Rule('RASGTPact_bind_SOS_GRB2_GAP_complex',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='Y', bpi3k=None) <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=None, st='GTP', act='N', bpi3k=None),
         *par['RASGTPact_bind_bound_GRB2_SOS'])

    Rule('RASGTP_unbind_SOS_GRB2_SHCP_complex',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', act='N', bpi3k=None) <>
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_unbind_GRB2_SOS'])

    Rule('RASGTP_unbind_SOS_GRB2_GAP_complex',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', act='N', bpi3k=None) <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_unbind_GRB2_SOS'])

    # Activation of RAF -> RAF:P by RAS-GTP
    Rule('RASGTP_bind_RAF',
         RAS(bsos=None, bpi3k=None, st='GTP', act='N', braf=None) + RAF(st='U', ser295='U', b=None) <>
         RAS(bsos=None, bpi3k=None, st='GTP', act='N', braf=1) % RAF(st='U', ser295='U', b=1),
         *par['RASGTP_bind_RAF'])

    Rule('RASGTP_RAF_cat',
         RAS(bsos=None, bpi3k=None, st='GTP', act='Y', braf=None) + RAF(st='P', ser295='U', b=None) <>
         RAS(bsos=None, bpi3k=None, st='GTP', act='N', braf=1) % RAF(st='U', ser295='U', b=1),
         *par['RASGTP_RAF_cat'])
    
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
    """ This is the akt part of the pathway from the Chen et al. 2009 paper.  Initial rules for all binding reactions were generated and then coded again using macros and higher order macros.  Initial parameters and conditions were taken from Chen et al. 2009 paper and supplementary, but were later modified in order to get the model working correctly.  This pathway follows AKT from its initial state to a phosphorylated and then double phosphorylated state before returning to unphosphorylated AKT.  The model works correctly, but parameters and rates may need to be modified in order to get best fit.  Parameters and rates included are from trial and error of what best fit the model.  The last unbinding reactions may not be needed because of the catalyze_state macros, but were left in just in case these are needed later.  
"""
    #This pathway coded by Tim O'Brien.
    Monomer('GAB1', ['bgrb2', 'bshp2', 'bpi3k', 'bpi3k2', 'bpi3k3', 'bpi3k4', 'bpi3k5', 'bpi3k6','batp','bERKPP','bPase9t','S'],{'S':['U','P','PP']})
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
    Initial(GAB1(bgrb2=None, bshp2=None, bpi3k=None, bpi3k2=None, bpi3k3=None, bpi3k4=None, bpi3k5=None, bpi3k6=None, batp=None, bERKPP=None, bPase9t=None, S='U'), GAB1_0)
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
    #GRB2 binds GAP-complex (without requiring SHC bound to complex):
    #Bind GRB2 without SOS already bound (two Chen-Sorger rate constants for different receptor dimers):
    #Rate 1: ErbB1/ErbB1 dimers (endosomal and plasma membrane), ErbB2/ErbB2 dimers (endosomal and plasma membrane), ErbB2/ErbB3 dimers (plasma membrane), and ErbB2/ErbB4 dimers (endosomal membrane):
    Rule('GRB2_bind_GAP_2',
         erbb(bd=1, ty='1') % erbb(bd=1, ty='1') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) <>
         erbb(bd=1, ty='1') % erbb(bd=1, ty='1') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP_2'])

    Rule('GRB2_bind_GAP_3',
         erbb(bd=1, ty='2') % erbb(bd=1, ty='2') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) <>
         erbb(bd=1, ty='2') % erbb(bd=1, ty='2') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP_2'])
    
    Rule('GRB2_bind_GAP_4',
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='3', loc='C') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) <>
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='3', loc='C') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP_2'])

    Rule('GRB2_bind_GAP_5',
          erbb(bd=1, ty='4', loc='E') % erbb(bd=1, ty='2', loc='E') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) <>
          erbb(bd=1, ty='4', loc='E') % erbb(bd=1, ty='2', loc='E') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
          *par['GRB2_bind_GAP_2'])

    #Rate 2: ErbB1/ErbBX, X=2, 3, 4 (endosomal and plasma membrane), ErbB2/ErbB3 dimers (endosomal membrane), and ErbB2/ErbB4 dimers (plasma membrane):
    for i in ['2', '3', '4']:
        Rule('GRB2_bind_GAP_6_'+i,
        erbb(bd=1, ty='1') % erbb(bd=1, ty=i) % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) <>
        erbb(bd=1, ty='1') % erbb(bd=1, ty=i) % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
        *par['GRB2_bind_GAP'])

    Rule('GRB2_bind_GAP_7',
         erbb(bd=1, ty='2', loc='E') % erbb(bd=1, ty='3', loc='E') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) <>
         erbb(bd=1, ty='2', loc='E') % erbb(bd=1, ty='3', loc='E') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP'])

    Rule('GRB2_bind_GAP_8',
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='4', loc='C') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) <>
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='4', loc='C') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP'])

    #Bind GRB2 to GAP with SOS already bound (one rate constant set for all dimer combinations):
    Rule('GRB2_bind_GAP_1',
         GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=1, bgab1=None, bcpp=None, bgap=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) <>
         GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=1, bgab1=None, bcpp=None, bgap=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         *par['GRB2_SOS_bind_GAP'])

    #GAB1 binds GAP-GRB2. Specify plasma membrane complexes in order to prevent complex building on endosomal receptors, so that degradation rxns (above in receptor events) can be simplified -- GAB1 complexes are not degraded as per Chen/Sorger model 

    Rule('GRB2_bind_GAB1',
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None, bcpp=None) + GAB1(bgrb2=None, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U') <>
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=2, bcpp=None) % GAB1(bgrb2=2, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'),
         *par['GRB2_bind_GAB1'])
    
    #GAP-GRB2-GAB1 phosphorylation - Rates from Table p. 5 Chen et al 2009
    Rule('GAB1_bind_ATP',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') + ATP(b=None) <>
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1),
         *par['GAB1_bind_ATP'])

    Rule('GAB1_phos',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1) >>
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + ADP(),
         par['GAB1_phos'])

    #SHP2 can desphosphorylate GAB1-P
    catalyze_state(SHP2(), 'bgab1', GAB1(bgrb2=ANY, bpi3k=None, batp=None, bERKPP=None, bPase9t=None), 'bshp2', 'S', 'P', 'U', (par['SHP2_dephos_GAB1P']))
   
    #After GAB1 phosphorylation, all receptor dimer combinations can bind a single PI3K
    #Chen/Sorger model gives two rate constant sets for different receptor dimers:
    #Rate 1: ErbB1/ErbB1, ErbB1/ErbB2, ErbB1/ErbB4, and ErbB2/ErbB4 dimers:
    for i in ['1', '2', '4']:
        Rule('GAB1_bind_PI3K_1_'+i,
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) <>
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_1'])

    Rule('GAB1_bind_PI3K_2',
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='4') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) <>
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='4') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_1'])

    #Rate 2: ErbB1/ErbB3, ErbB2/ErbB2, and ErbB2/ErbB3 dimers:
    for i in ['1', '2']:
        Rule('GAB1_bind_PI3K_3_'+i,
             erbb(bd=ANY, ty=i) % erbb(bd=ANY, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) <>
             erbb(bd=ANY, ty=i) % erbb(bd=ANY, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_2'])

    Rule('GAB1_bind_PI3K_4',
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='2') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) <>
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='2') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_2'])

    #GAB1-PI3K bound to complex containing ErbB2/ErbB3 binds 1-6 PIP2 (creates chains; doesn't necessarily represent biology but accurately reproduces Chen Sorger 2009 model).
    #First bind a single PIP2 to PI3K complex - this rule created by catalyze_state below

    #Then create chains of up to 6 PIP2 molecules attached to a single PI3K:
    assemble_chain_sequential_base(PI3K(berb=None, bras=None, bgab1=ANY, bpip=None), 'bpip', PIP(bakt=None, both=None, S='PIP2'), 'bpi3k_self', 'bself2', 6, [par['PIP2_chain_PI3K']]*5, \
                                   erbb(bd=ANY, ty='2', b=ANY, loc='C') % erbb(bd=ANY, ty='3', b=None, loc='C') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P'))
    
    #To accurately reproduce Chen Sorger model, allow PIP2 to catalyze PIP2->PIP3 conversion of final chain unit.
    Rule('PIP2_self_catalysis_1',
         PIP(bakt=None, both=None, S='PIP2', bpi3k_self=ANY, bself2=1) % PIP(bakt=None, both=None, S='PIP2', bpi3k_self=1, bself2=None) >>
         PIP(bakt=None, both=None, S='PIP2', bpi3k_self=ANY, bself2=None) + PIP(bakt=None, both=None, S='PIP3', bpi3k_self=None, bself2=None),
          par['PIP2_self_catalysis'])
    
    #ErbB2-ErbB3 dimers contain 6 binding domains for PI3K (don't need to be bound to adaptor complex).
    """This is the improved ErbB2/ErB3-PI3K sequence of events -- different from Chen Sorger 2009"""
    #for i in range(1,6):
    #   rfix=Rule('ErbB2_3bindPI3K_'+str(i),
    #            erbb(bd=1, ty='2') % erbb(bd=1, ty='3') + PI3K(bpip=None, berb=None) <>  erbb(bd=1, ty='2') % erbb(bd=1, ty='3') + PI3K(bpip=None, berb=1),
    #            Parameter('GAB1PI3Kf'+str(i), 1e-5),
    #              Parameter('GAB1PI3Kr'+str(i), 1e-1))
    #   m=rfix.reactant_pattern.complex_patterns[0].monomer_patterns[1]
    #  m.site_conditions = {'bd' : 1, 'ty':'3', 'pi3k'+str(i):None}
    #  m=rfix.product_pattern.complex_patterns[0].monomer_patterns[1]
    #  m.site_conditions = {'bd' : 1, 'ty':'3', 'pi3k'+str(i):1}

    #PI3K bound directly to ErbB2-ErbB3 dimers (at any of 6 sites) can catalyze PIP2->PIP3:
    #catalyze_state(PI3K(berb=ANY), 'bpip', PIP(bakt=None, both=None), 'bpi3k', 'S', 'PIP2', 'PIP3', (1e-5, 1e-1, 1e-1))
    

    #PI3K bound to complex catalyzes PIP2 -> PIP3
    #Two rate sets for initial binding in Chen/Sorger model:
    #Rate 1: ErbB1/ErbBX dimers:
    for i in ['1', '2', '3', '4']:
        Rule('PIP2_bind_PI3K_1_'+i,
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=None) <>
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1),
             *par['PIP2_bind_PI3K_1'])

    #Rate 2: ErbB2/ErbBX dimers, X=2, 3, 4:
    #FIXME: What is up with v701 in reaction list?
    for i in ['2', '3', '4']:
        Rule('PIP2_bind_PI3K_2_'+i,
        erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=None) <>
        erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1),
             *par['PIP2_chain_PI3K'])
    
    #Two catalysis rates in Chen/Sorger model:
    #Rate 1: ErbB2/ErbB3 dimers:
    Rule('PIP2_PI3K_catalysis_1',
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1) >>
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', both=None, bakt=None, bself2=None, bpi3k_self=None),
         par['PIP2_self_catalysis'])

    #Rate 2: All other dimers:
    for i in ['1', '2', '3', '4']:
        Rule('PIP2_PI3K_catalysis_2_'+i,
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1) >>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', both=None, bakt=None, bself2=None, bpi3k_self=None),
         par['PIP2_PI3K_catalysis'])

    for i in ['2', '4']:
        Rule('PIP2_PI3K_catalysis_3_'+i,
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1) >>
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP3', both=None, bakt=None, bself2=None, bpi3k_self=None),
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
    #ERK:P:P phosphorylates GAP-GRB2-GAB1:P (making it unable to bind PI3K)
    catalyze_state(ERK(st='PP'), 'b', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None, bpi3k2=None, bpi3k3=None, bpi3k4=None, bpi3k5=None, bpi3k6=None), 'bERKPP', 'S', 'P', 'PP', (par['ERKPP_phos_GAB1P']))

    #GAP-GRB2-GAB1:P:P is dephosphorylated by Pase9t
    catalyze_state(Pase9t(), 'bgab1', GAB1(bgrb2=ANY), 'bPase9t', 'S', 'PP', 'P', (par['Pase9t_dephos_GAB1PP']))

    #ERK:P:P phosphorylates GRB2-SOS, preventing RAS-GDP->RAS-GTP conversion
    #To conform with Chen/Sorger model, this only effects ErbB1/ErbB1 dimers containing SOS and free SOS, and phosphorylated SOS can only bind ErbB1/ErbB1 complexes, not free GRB2:
    catalyze_state(ERK(st='PP'), 'b', SOS(bgrb=None, bras=None), 'bERKPP', 'st', 'U', 'P', (par['ERKPP_phos_SOS']))

    Rule('ERKPP_bind_SOS_1',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=None, st='U') + ERK(st='PP', b=None) <>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=1, st='U') % ERK(st='PP', b=1),
        *par['ERKPP_phos_SOS'][0:2])

    Rule('ERKPP_bind_SOS_2',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='U', bgrb=ANY) + ERK(st='PP', b=None) <>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=1, st='U', bgrb=ANY) % ERK(st='PP', b=1),
        *par['ERKPP_phos_SOS'][0:2])

    Rule('ERKPP_phos_SOS_1',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=1, st='U') % ERK(st='PP', b=1) >>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=None, st='P') + ERK(st='PP', b=None),
         par['ERKPP_phos_SOS'][2])

    Rule('ERKPP_phos_SOS_2',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=1, st='U', bgrb=ANY) % ERK(st='PP', b=1) >>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='P', bgrb=ANY) + ERK(st='PP', b=None),
         par['ERKPP_phos_SOS'][2])

    Rule('SOSP_bind_GRB2_1',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) + SOS(bras=None, bgrb=None, bERKPP=None, st='P') <>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=1, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='P'),
         *par['SOSP_bind_GRB2'])

    Rule('SOSP_bind_GRB2_2',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=ANY) + SOS(bras=None, bERKPP=None, st='P', bgrb=None) <>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='P', bgrb=1),
         *par['SOSP_bind_GRB2'])

    #AKT:P:P phosphorylates RAF:P at Ser295, preventing MEK phosphorylation.
    catalyze_state(AKT(S='PP', bpip3=None), 'both', RAF(st='P'), 'b', 'ser295', 'U', 'P', (par['AKTPP_phos_RAFP']))

    #RAS-GDP binds PI3K and PI3K catalyzes GDP --> GTP transformation. (Note: Matches model, but is this actually what's happening biologically?)
    catalyze_state(PI3K(bgab1=ANY, bpip=None), 'bras', RAS(bsos=None, braf=None), 'bpi3k', 'st', 'GDP', 'GTP', (par['RASGDP_bind_PI3K']))

