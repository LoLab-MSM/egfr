"""
Overview:
========

PySB implementation of the ErbB related MAPK and AKT signaling
pathways originally published in [Chen2009]_.

This file containst functions that implement the ErbB execution
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

from .parameter_dict_BT474 import parameter_dict as par
#FIXME: What is Inh in reaction list?
        
# Monomer declarations
# ====================

def rec_monomers():

    """ Declares the ErbB receptor interactions.
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    # Monomer('EGF', ['b', 'st'], {'st':['M', 'E']}) # Epidermal Growth Factor ligand
    # Monomer('HRG', ['b']) # Heregulin ligand
    Monomer('erbb', ['bl', 'bd', 'b', 'ty', 'st', 'loc', 'pi3k1', 'pi3k2', 'pi3k3', 'pi3k4', 'pi3k5', 'pi3k6', 'cpp'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E'], 'cpp':['Y', 'N']}) # bl: lig, bd: dimer, b: binding, ty: rec type, st: (U)n(P)hosphorylated, loc: (C)yto 'brane or (E)ndosome 'brane, cpp: No real biophysical meaning; useful model marker for presence of CPP bound downstream.

    Monomer('DEP', ['b'])
    Monomer('ATP', ['b'])
    Monomer('ADP')
    Monomer('CPP', ['b', 'loc'], {'loc':['C', 'E']})
    Monomer('LAP', ['b', 'loc'], {'loc':['M', 'C', 'E']})
    Monomer('BKM120', ['b', 'loc'], {'loc':['M', 'C']})

def rec_initial_lig_hEGF():
    Parameter('EGF_0',      3.01e12) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',         0) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()
    
def rec_initial_lig_lEGF():
    Parameter('EGF_0',      6.02e9) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',         0) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()

def rec_initial_lig_pEGF():
    Parameter('EGF_0', 6.02e8) # 1 pm EGF = 6.02e8 molec/cell

def rec_initial_lig_hHRG():
    Parameter('EGF_0',      0) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',       3.01e12) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()

def rec_initial_lig_lHRG():
    Parameter('EGF_0',      0) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    Parameter('HRG_0',         6.02e9) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell
    alias_model_components()

def rec_initial_inhib_LAP():
    Parameter('LAP_0', 6.02e14) # 1 microM lapatinib = 6.02e14 molec/cell

def rec_initial_inhib_PI3K():
    Parameter('BKM120_0', 0) # 1 microM BKM120 = 6.02e14 molec/cell

def rec_initial():
    # # Initial concentrations (except DEP1) for all cell types taken from Chen et al 2009 -- see Jacobian files
    alias_model_components()

    # Initial(EGF(b=None, st='M'), EGF_0)
    # Initial(HRG(b=None), HRG_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='1', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb1_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='2', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb2_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='3', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb3_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='4', st='U', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), erbb4_0)
    Initial(ATP(b=None), ATP_0)
    Initial(DEP(b=None), DEP_0)
    Initial(CPP(b=None, loc='C'), CPP_0)
    Initial(LAP(b=None, loc='M'), LAP_0)
    #Initial(BKM120(b=None, loc='M'), BKM120_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='1', st='P', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), ErbB1P_0)
    Initial(erbb(bl=None, bd=None, b=None, ty='2', st='P', loc='C', pi3k1=None, pi3k2=None, pi3k3=None, pi3k4=None, pi3k5=None, pi3k6=None, cpp='N'), ErbB2P_0)
            
def rec_events():
    """ Describe receptor-level events here. 
    """
    
    # Parameter definitions
    # =====================
    # Alias model components for names in present namespace
    alias_model_components()
    # EGF / HRG binding to receptors
    # EGF / HRG receptor binding rates obtained from Chen et al Jacobian files
    # bind_table([[                                                          EGF(st='M'),                                   HRG],
    #             [erbb(ty='1', bd=None, b=None, st='U', loc='C'),           (par['EGF_bind_ErbB1']),               None],
    #             [erbb(ty='3', bd=None, b=None, st='U', loc='C'),            None,                                 (par['HRG_bind_ErbB3'])],
    #             [erbb(ty='4', bd=None, b=None, st='U', loc='C'),            None,                                 (par['HRG_bind_ErbB4'])]],
    #             'bl', 'b')
    
    # ErbB dimerization
    # Dimerization rates obtained from Chen et al Jacobian files

    # MODIFICATION - Added lapatinib transport into the cell.
    # FURTHER MOD - Added BKM120 (PI3K inhibitor) transport into cell.

    equilibrate(LAP(loc='M', b=None), LAP(loc='C', b=None), par['LAP_transport_MC'])

    #equilibrate(BKM120(loc='M', b=None), BKM120(loc='C', b=None), par['BKM120_transport_MC'])
    
    # MODIFICATION for Rexer model -- Monomers are not required to contain a ligand to form dimers but won't form active phosphorylated dimers without ligand (except ErbB2 containing dimers).
    erbb1Lig = erbb(ty='1', b=None, st='U', loc='C')
    erbb2Lig = erbb(ty='2', bl=None, b=None, st='U', loc='C')
    erbb3Lig = erbb(ty='3', b=None, st='U', loc='C')
    erbb4Lig = erbb(ty='4', b=None, st='U', loc='C')
    bind_table([[                          erbb1Lig,                    erbb2Lig,                     erbb3Lig, erbb4Lig],
                [erbb1Lig,                 (par['ErbB1_bind_ErbB1']),   None,                         None,     None],
                [erbb2Lig,                 (par['ErbB1_bind_ErbB2']),   (par['ErbB2_bind_ErbB2']),    None,     None],
                [erbb3Lig,                 (par['ErbB1_bind_ErbB3']),   (par['ErbB2_bind_ErbB3']),    None,     None],
                [erbb4Lig,                 (par['ErbB1_bind_ErbB4']),   (par['ErbB2_bind_ErbB4']),    None,     None]],
        'bd', 'bd')

    # MODIFICATION for Rexer model -- Added lapatinib binding to ErbB1 and ErbB2

    bind_table([[             LAP(loc='C')],
                [erbb1Lig,    (par['LAP_bind_ErbB1'])],
                [erbb2Lig,    (par['LAP_bind_ErbB2'])]],
                'b', 'b')

    # MODIFICATION for Rexer data: All ErbB2 containing dimers can be directly phosphorylated -- allows signal without ligand
    # Also: All other dimers are only phosphorylated by interaction with dissassociated ErbB2P -- no phosphorylation of others w/o ligand present.
    # Also: Phosphorylation will not occur in lapatinib is bound.

    for i in ['1', '2', '3', '4']:
        Rule('ATP_bind_ErbB2'+i,
             erbb(ty='2', st='U', loc='C', b=None, bd=1) % erbb(ty=i, st='U', loc='C', b=None, bd=1) + ATP(b=None) |
             erbb(ty='2', st='U', loc='C', b=2, bd=1) % erbb(ty=i, st='U', loc='C', b=None, bd=1) % ATP(b=2),
             *par['ErbB2'+i+'_bind_ATP'])
    
    for i in ['1', '2', '4']:
        Rule('DEP_bind_ErbB'+i,
             erbb(ty=i, st='P', loc='C', b=None, bd=1) % erbb(st='P', loc='C', b=None, bd=1) + DEP(b=None) |
             erbb(ty=i, st='P', loc='C', b=2, bd=1) % erbb(st='P', loc='C', b=None, bd=1) % DEP(b=2),
             *par['ErbBP'+i+'_bind_DEP'])

    # Cross phosphorylation: only erbb1, 2, and 4 have ATP, and they can cross-phosphorylate any other receptor (once activated by ligand or in Rexer model, 2 doesn't have to be activated)
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


    #ErbB2 lateral signaling - ErbB2P-ErbB2P dimers can only form by the dissociation of ligand-containing, phosphorylated dimers containing ErbB2 (in Rexer model, these can form independently of any other ErbB receptor; in normal model, they need another ligand-activated receptor).  The monomeric activated ErbB2 can then bind and activate other monomers (ErbB1, 3, or 4 -- allows EGF signal to be transmitted by ErbB2/ErbB3 and ErbB2/ErbB4 complexes, even though 3 and 4 can't bind EGF) or bind another phosphorylated ErbB2 to form an active complex (that still requires an EGF signal to get started)

    bind(erbb(ty='1', bd=None, st='P', b=None, loc='C'), 'bd', erbb(bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB1P_ErbBXP_bind'])

    bind(erbb(ty='2', bd=None, st='P', b=None, loc='C'), 'bd', erbb(bd=None, st='P', b=None, loc='C'), 'bd', par['ErbB2P_ErbBXP_bind'])

    for i in ['1', '3', '4']:
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
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None) | erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None),
         *par['kint_no_cPP_1'])

    Rule("rec_intern_2",
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None) |
         erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None),
         *par['kint_no_cPP_1'])

    Rule('rec_intern_3',
         erbb(bd=1, loc='C', cpp='N', ty='1') % erbb(bd=1, loc='C', cpp='N', ty='1') % GAP(bd=ANY, b=None, bgrb2=None) |
         erbb(bd=1, loc='E', cpp='N', ty='1') % erbb(bd=1, loc='E', cpp='N', ty='1') % GAP(bd=ANY, b=None, bgrb2=None),
         *par['kint_no_cPP_1'])

    Rule('rec_intern_4',
         erbb(bd=1, loc='C', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='C', cpp='N', ty='1', b=None) |
         erbb(bd=1, loc='E', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='E', cpp='N', ty='1', b=None),
         *par['kint_no_cPP_1'])

    Rule('rec_intern_5',
         erbb(bd=None, loc='C', cpp='N', ty='1', b=None) |
         erbb(bd=None, loc='E', cpp='N', ty='1', b=None),
         *par['kint_no_cPP_1'])

    # Rate 2: Set to 0 in Chen/Sorger files and not implemented.  Would have internalized single ErbB2, 3, and 4, as well as ErbB2/3,4:GAP:SHC complexes (phos/unphos).
    # Rate 3: These rules internalize ErbB1/ErbBX dimers, ErbB2/ErbB2:GAP:SHC complexes and intermediates, and ErbB2/ErbB3 and ErbB2/ErbB4 dimers.
    for i in ['2', '3', '4']:
        Rule('rec_intern_6_'+i,
             erbb(bd=1, loc='C', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='C', cpp='N', b=None, ty=i) |
             erbb(bd=1, loc='E', cpp='N', ty='1', st='P', b=None) % erbb(bd=1, loc='E', cpp='N', b=None, ty=i),
            *par['kint_no_cPP_2'])

        Rule('rec_intern_7_'+i,
             erbb(bd=1, loc='C', cpp='N', ty=i, b=None) % erbb(bd=1, loc='C', cpp='N', b=None, st='P', ty='1') |
             erbb(bd=1, loc='E', cpp='N', ty=i, b=None) % erbb(bd=1, loc='E', cpp='N', b=None, st='P', ty='1'),
            *par['kint_no_cPP_2'])

    Rule('rec_intern_8',
         erbb(bd=1, loc='C', cpp='N', ty='2') % erbb(bd=1, loc='C', cpp='N', ty='2') % GAP(bd=ANY, b=None, bgrb2=None) |
         erbb(bd=1, loc='E', cpp='N', ty='2') % erbb(bd=1, loc='E', cpp='N', ty='2') % GAP(bd=ANY, b=None, bgrb2=None),
         *par['kint_no_cPP_2'])

    Rule('rec_intern_9',
         erbb(bd=1, loc='C', cpp='N', ty='2') % erbb(bd=1, loc='C', cpp='N', ty='2') % GAP(bd=ANY, b=2, bgrb2=None) % SHC(bgap=2, batp=None, bgrb=None) |
         erbb(bd=1, loc='E', cpp='N', ty='2') % erbb(bd=1, loc='E', cpp='N', ty='2') % GAP(bd=ANY, b=2, bgrb2=None) % SHC(bgap=2, batp=None, bgrb=None),
         *par['kint_no_cPP_2'])

    for i in ['2', '3', '4']:
        Rule('rec_intern_10_'+i,
             erbb(bd=1, loc='C', cpp='N', ty='2', st='P', b=None) % erbb(bd=1, loc='C', cpp='N', b=None, ty=i) |
             erbb(bd=1, loc='E', cpp='N', ty='2', st='P', b=None) % erbb(bd=1, loc='E', cpp='N', b=None, ty=i),
            *par['kint_no_cPP_2'])

        Rule('rec_intern_11_'+i,
             erbb(bd=1, loc='C', cpp='N', ty=i, b=None) % erbb(bd=1, loc='C', cpp='N', b=None, st='P', ty='2') |
             erbb(bd=1, loc='E', cpp='N', ty=i, b=None) % erbb(bd=1, loc='E', cpp='N', b=None, st='P', ty='2'),
            *par['kint_no_cPP_2'])
        
    # CPP bound to receptors can catalyze their internalization (when they are bound to any complex containing GRB2, except GAB1 complex):
    # Binding to CPP and internalization rates are conflated in order to better match Chen-Sorger model.
    Rule('CPP_bind_GAP_GRB2',
         CPP(loc='C', b=None) + erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(ty='1', bd=1, loc='C', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None) |
         erbb(ty='1', bd=1, loc='E', cpp='Y') % erbb(ty='1', bd=1, loc='E', cpp='Y') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3, bgab1=None, b=None) % CPP(loc='E', b=3),
         *par['CPP_bind_ErbB1dimers'])

    Rule('CPP_bind_SHC_GRB2',
         erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(ty='1', bd=1, loc='C', cpp='N') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None, bgap=None) + CPP(loc='C', b=None) |
         erbb(ty='1', bd=1, loc='E', cpp='Y') % erbb(ty='1', bd=1, loc='E', cpp='Y') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None, bgap=None) % CPP(loc='E', b=4),
         *par['CPP_bind_ErbB1dimers'])

    for i in ['1', '2', '3', '4']:
        Rule('CPP_bind_ErbB1_RASGTP_complex_'+i,
        erbb(ty='1', bd=1, loc='C', cpp='N') % erbb(bd=1, ty=i, loc='C', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') + CPP(loc='C', b=None) |
        erbb(ty='1', bd=1, loc='E', cpp='N') % erbb(bd=1, ty=i, loc='E', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=5, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') % CPP(loc='E', b=5),
        *par['CPP_bind_ErbB1dimers'])

        Rule('CPP_bind_ErbB1_RASGTP_complex2_'+i,
         erbb(ty=i, bd=1, loc='C', cpp='N') % erbb(bd=1, ty='1', loc='C', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') + CPP(loc='C', b=None) |
         erbb(ty=i, bd=1, loc='E', cpp='N') % erbb(bd=1, ty='1', loc='E', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=5, bsos=3) % SOS(bgrb=3, bERKPP=None, bras=4) % RAS(bsos=4, braf=None, bpi3k=None, st='GTP') % CPP(loc='E', b=5),
         *par['CPP_bind_ErbB1dimers'])
    
    Rule('CPPE_bind_GAP_GRB2',
         erbb(bd=1, ty='1', loc='E', cpp='N') % erbb(bd=1, ty='1', loc='E', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=None, bgab1=None, b=None) + CPP(loc='E', b=None) |
         erbb(bd=1, ty='1', loc='E', cpp='Y') % erbb(bd=1, ty='1', loc='E', cpp='Y') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3, bgab1=None, b=None) % CPP(loc='E', b=3),
         *par['CPPE_bind_ErbB1dimers'])

    Rule('CPPE_bind_SHC_GRB2',
         erbb(bd=1, ty='1', loc='E', cpp='N') % erbb(bd=1, ty='1', loc='E', cpp='N') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None, bgap=None) + CPP(loc='E', b=None) |
         erbb(bd=1, ty='1', loc='E', cpp='Y') % erbb(bd=1, ty='1', loc='E', cpp='Y') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None, bgap=None) % CPP(loc='E', b=4),
         *par['CPPE_bind_ErbB1dimers'])
    
    Rule("CPP_intern",
         CPP(loc='E', b=None) | CPP(loc='C', b=None),
         *par['CPP_int'])

    # MODIFICATION for Rexer model -- Added lapatinib internalization and degradation, and the internalization and degradation of lapatinib-containing complexes.

    equilibrate(LAP(b=None, loc='C'), LAP(b=None, loc='E'), par['LAP_intern'])

    Rule('LAP_ErbB1_intern',
         LAP(b=1, loc='C') % erbb(bd=None, b=1, ty='1', loc='C', cpp='N') |
         LAP(b=1, loc='E') % erbb(bd=None, b=1, ty='1', loc='E', cpp='N'),
         *par['LAP_ErbB1_intern'])

    Rule('LAP_ErbB1d_intern',
         LAP(b=1, loc='C') % erbb(bd=2, b=1, ty='1', loc='C', cpp='N') % erbb(bd=2, b=None, loc='C', cpp='N') |
         LAP(b=1, loc='E') % erbb(bd=2, b=1, ty='1', loc='E', cpp='N') % erbb(bd=2, b=None, loc='E', cpp='N'),
         *par['LAP_ErbB1d_intern'])

    Rule('LAP_ErbB2_intern',
         LAP(b=1, loc='C') % erbb(bd=None, b=1, ty='2', loc='C', cpp='N') |
         LAP(b=1, loc='E') % erbb(bd=None, b=1, ty='2', loc='E', cpp='N'),
         *par['LAP_ErbB2_intern'])

    Rule('LAP_ErbB2d_intern',
         LAP(b=1, loc='C') % erbb(bd=2, b=1, ty='2', loc='C', cpp='N') % erbb(bd=2, b=None, loc='C', cpp='N') |
         LAP(b=1, loc='E') % erbb(bd=2, b=1, ty='2', loc='E', cpp='N') % erbb(bd=2, b=None, loc='E', cpp='N'),
         *par['LAP_ErbB2d_intern'])

    Rule('LAP2_ErbB_intern',
         LAP(b=1, loc='C') % LAP(b=2, loc='C') % erbb(bd=3, b=1, loc='C', cpp='N') % erbb(bd=3, b=2, loc='C', cpp='N') |
         LAP(b=1, loc='E') % LAP(b=2, loc='E') % erbb(bd=3, b=1, loc='E', cpp='N') % erbb(bd=3, b=2, loc='E', cpp='N'),
         *par['LAP2_ErbB_intern'])

    degrade(LAP(b=None, loc='E'), par['LAP_deg'])

    degrade(LAP(b=1, loc='E') % erbb(bd=None, b=1, ty='1', loc='E', cpp='N'), par['LAP_ErbB1_deg'])

    degrade(LAP(b=1, loc='E') % erbb(bd=2, b=1, ty='1', loc='E', cpp='N') % erbb(bd=2, b=None, loc='E', cpp='N'), par['LAP_ErbB1d_deg'])

    degrade(LAP(b=1, loc='E') % erbb(bd=None, b=1, ty='2', loc='E', cpp='N'), par['LAP_ErbB2_deg'])

    degrade(LAP(b=1, loc='E') % erbb(bd=2, b=1, ty='2', loc='E', cpp='N') % erbb(bd=2, b=None, loc='E', cpp='N'), par['LAP_ErbB2d_deg'])

    degrade(LAP(b=1, loc='E') % LAP(b=2, loc='E') % erbb(bd=3, b=1, loc='E', cpp='N') % erbb(bd=3, b=2, loc='E', cpp='N'), par['LAP2_ErbB_deg'])
         
    # Receptor degradation
    # This degrades all receptor combos within an endosome
    # The Chen/Sorger model implements different degradation rates for different species:
    # Rate 1: These rules degrade all 2EGF:ErbB1/ErbB1 complexes in the MAPK pathway (i.e. ErbB1/ErbB1:GAP:GRB2:SOS:RAS-GDP and ErbB1/ErbB1:GAP:SHC-P:GRB2:SOS:RAS-GDP and all intermediates in their formation.), as well as single ErbB1.  
    # degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3) % erbb(bd=1, loc='E', ty='1', bl=4) % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=None, bgab1=None, b=None), par['kdeg_1'])

    # degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3) % erbb(bd=1, loc='E', ty='1', bl=4) % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None), par['kdeg_1'])

    # degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3) % erbb(bd=1, loc='E', ty='1', bl=4) % GAP(bd=ANY, b=None, bgrb2=None), par['kdeg_1'])

    # degrade(EGF(b=3) % EGF(b=4) % erbb(bd=1, loc='E', ty='1', bl=3) % erbb(bd=1, loc='E', ty='1', bl=4), par['kdeg_1'])

    degrade(erbb(bd=None, loc='E', ty='1'), par['kdeg_1'])

    # Rate 2: These rules degrade all ErbB1/ErbBX species and all ErbB2/ErbB2 species in MAPK pathway.  Chen/Sorger model also included degradation of single ErbB2, 3, and 4 under this constant, but as these are never internalized by Chen/Sorger rule set, these degradation rxns were ignored.
    for i in ['2', '3', '4']:
        degrade(erbb(bd=1, loc='E', ty='1') % erbb(bd=1, loc='E', ty=i) % GAP(bd=ANY), par['kdeg_2'])
        degrade(erbb(bd=1, loc='E', ty=i) % erbb(bd=1, loc='E', ty='1') % GAP(bd=ANY), par['kdeg_2'])

    degrade(erbb(bd=1, loc='E', ty='2') % erbb(bd=1, loc='E', ty='2') % GAP(bd=ANY), par['kdeg_2'])

    # Rate 3: These rules degrade all ErbB2/ErbB3 and all ErbB2/ErbB4 complexes in MAPK pathway.
    for i in ['3', '4']:
        degrade(erbb(bd=1, loc='E', ty='2') % erbb(bd=1, loc='E', ty=i) % GAP(bd=ANY), par['kdeg_3'])
        degrade(erbb(bd=1, loc='E', ty=i) % erbb(bd=1, loc='E', ty='2') % GAP(bd=ANY), par['kdeg_3'])

    # Rate 4: degradation of EGF
    # degrade(EGF(b=None, st='E'), par['kdeg_4'])

    # Rate 5: Degrades ErbB1/ErbBX, ErbB2/ErbB3, and ErbB2/ErbB4 dimers (when no complex attached).
    for i in ['2', '3', '4']:
        degrade(erbb(bd=1, loc='E', ty='1', b=None) % erbb(bd=1, loc='E', ty=i, b=None), par['kdeg_5'])
        degrade(erbb(bd=1, loc='E', ty=i, b=None) % erbb(bd=1, loc='E', ty='1', b=None), par['kdeg_5'])

    for i in ['3', '4']:
        degrade(erbb(bd=1, loc='E', ty='2', b=None) % erbb(bd=1, loc='E', ty=i, b=None), par['kdeg_5'])
        degrade(erbb(bd=1, loc='E', ty=i, b=None) % erbb(bd=1, loc='E', ty='2', b=None), par['kdeg_5'])

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
    Monomer('ERK', ['b', 'st', 'loc'], {'st':['U', 'P', 'PP'], 'loc':['C', 'N']})

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
    Initial(ERK(b=None, st='U', loc='C'), ERK_0)
    Initial(PP1(b=None), PP1_0)
    Initial(PP2(b=None), PP2_0)
    Initial(PP3(b=None), PP3_0)
    Initial(GRB2(b=None, bsos=1, bgap=None, bgab1=None, bcpp=None) % SOS(bgrb=1, bras=None, bERKPP=None, st='U'), GRB2_SOS_0)
    Initial(ERK(b=None, st='PP', loc='C'), ERKPP_0)

    
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
         MatchOnce(erbb(bd=1, b=None, st='P', ty='1') % erbb(bd=1, b=None, st='P', ty='1')) + GAP(bd=None, b=None, bgrb2=None) |
         MatchOnce(erbb(bd=1, b=2,    st='P', ty='1') % erbb(bd=1, b=None, st='P', ty='1') % GAP(bd=2, b=None, bgrb2=None)),
         *par['ErbB_bind_GAP_1'])

    for i in ['2', '3', '4']:
        Rule('GAP_binding_2'+i,
             MatchOnce(erbb(bd=1, b=None, st='P', ty='2') % erbb(bd=1, b=None, st='P', ty=i)) + GAP(bd=None, b=None, bgrb2=None) |
             MatchOnce(erbb(bd=1, b=2, st='P', ty='2') % erbb(bd=1, b=None, st='P', ty=i) % GAP(bd=2, b=None, bgrb2=None)),
             *par['ErbB_bind_GAP_1'])

        Rule('GAP_binding_2_2'+i,
             MatchOnce(erbb(bd=1, b=None, st='P', ty=i) % erbb(bd=1, b=None, st='P', ty='2')) + GAP(bd=None, b=None, bgrb2=None) |
             MatchOnce(erbb(bd=1, b=2, st='P', ty=i) % erbb(bd=1, b=None, st='P', ty='2') % GAP(bd=2, b=None, bgrb2=None)),
             *par['ErbB_bind_GAP_1'])

    # Rate 2: ErbB1/ErbBX, X=2, 3, 4  Note: In Chen/Sorger rxn list, plasma membrane ErbB1/ErbB2 dimers are assigned Rate 1 (above); however the other 5 ErbB1/ErbBX combinations (plasma and endosomal membranes) are assigned Rate 2.  ErbB1/ErbB2 was assigned the latter in this model under the assumption that this was accidental.
    for i in ['2', '3', '4']:
        Rule('GAP_binding_3'+i,
             MatchOnce(erbb(bd=1, b=None, st='P', ty='1') % erbb(bd=1, b=None, st='P', ty=i)) + GAP(bd=None, b=None, bgrb2=None) |
             MatchOnce(erbb(bd=1, b=2, st='P', ty='1') % erbb(bd=1, b=None, st='P', ty=i) % GAP(bd=2, b=None, bgrb2=None)),
             *par['ErbB_bind_GAP_2'])

        Rule('GAP_binding_3_2'+i,
             MatchOnce(erbb(bd=1, b=None, st='P', ty=i) % erbb(bd=1, b=None, st='P', ty='1')) + GAP(bd=None, b=None, bgrb2=None) |
             MatchOnce(erbb(bd=1, b=2, st='P', ty=i) % erbb(bd=1, b=None, st='P', ty='1') % GAP(bd=2, b=None, bgrb2=None)),
             *par['ErbB_bind_GAP_2'])
    
    # SHC binds to GAP-complex
    # Chen-Sorger model assigns 2 sets of rate constants to different dimer combinations.  The kf is the same variable; two different kr variables are used but are assigned the same values in the Jacobian files.  These have been combined into one set in this model.
    
    bind(GAP(bd=ANY, bgrb2=None), 'b', SHC(batp=None, st='U', bgrb=None), 'bgap', par['GAP_bind_SHC'])

    #SHC:P binds GAP
    bind(GAP(bd=ANY, bgrb2=None), 'b', SHC(batp=None, st='P', bgrb=None), 'bgap', par['GAP_bind_SHCP'])

    #SHC:P-GRB2 binds GAP
    Rule('GAP_bind_SHCP_GRB2',
         GAP(bd=ANY, b=None, bgrb2=None) + SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=1) |
         GAP(bd=ANY, b=2, bgrb2=None) % SHC(batp=None, st='P', bgrb=1, bgap=2) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=1),
         *par['GAP_bind_SHCP'])

    # Bound and unbound SHC phosphorylation - These are represented by two kf, kr pairs in the Chen-Sorger model:
    Rule('SHC_phos',
         GAP(bd=ANY, b=1, bgrb2=None) % SHC(bgap=1, bgrb=None, batp=None, st='U') |
         GAP(bd=ANY, b=1, bgrb2=None) % SHC(bgap=1, bgrb=None, batp=None, st='P'),
         *par['SHC_phos'])

    # Forward rate is set to 0 (only dephosphorylation is occurring).  
    Rule('SHC_unbound_phos',
         SHC(bgap=None, bgrb=None, batp=None, st='U') |
         SHC(bgap=None, bgrb=None, batp=None, st='P'),
         *par['SHC_unbound_phos'])
    
    # GRB2 binds to GAP-SHC:P with or without SOS:
    Rule('GRB2_bind_GAP_SHCP_1',
         SHC(batp=None, st='P', bgrb=None, bgap=ANY) + GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) |
         SHC(batp=None, st='P', bgrb=2, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         *par['GRB2_SOS_bind_SHCP_GAP'])

    bind(SHC(batp=None, st='P', bgap=ANY), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None), 'b', par['GRB2_bind_GAP'])

    # SHC:P can bind GRB2-SOS without being attached to GAP:
    Rule('SHCP_bind_GRB2SOS',
         SHC(batp=None, st='P', bgrb=None, bgap=None) + GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) |
         SHC(batp=None, st='P', bgrb=2, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         *par['SHCP_bind_GRB2SOS'])

    # GAP can bind the free SHC:P-GRB2-SOS complex:
    Rule('GAP_bind_SHCP_GRB2_SOS',
         GAP(bd=ANY, b=None, bgrb2=None) + SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=2, bcpp=None, b=1) % SOS(bras=None, bERKPP=None, st='U', bgrb=2) |
         GAP(bd=ANY, b=3, bgrb2=None) % SHC(batp=None, st='P', bgrb=1, bgap=3) % GRB2(bgap=None, bgab1=None, bsos=2, bcpp=None, b=1) % SOS(bras=None, bERKPP=None, st='U', bgrb=2),
         *par['GAP_bind_SHCP_GRB2_SOS'])

    # GRB2 and SOS bind/disassociate:
    bind(GRB2(bgap=None, bgab1=None, b=None, bcpp=None), 'bsos', SOS(bras=None, bERKPP=None, st='U'), 'bgrb', par['GRB2_bind_SOS'])

    #Although no free SOS is present initially in Chen Sorger model, GRB2-SOS can disassociate (see above), so these are necessary.
    # SOS binds to GAP-SHC:P-GRB2  
    bind(GRB2(b=ANY, bgap=None, bgab1=None, bcpp=None), 'bsos', SOS(bras=None, st='U', bERKPP=None), 'bgrb', par['SOS_bind_GAP_SHCP_GRB2'])

    #SOS binds SHC:P-GRB2 without complex
    Rule('SOS_bind_SHCP_GRB2',
         SOS(bras=None, st='U', bERKPP=None, bgrb=None) + SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=1) |
         SHC(batp=None, st='P', bgrb=1, bgap=None) % GRB2(bgap=None, bgab1=None, bsos=2, bcpp=None, b=1) % SOS(bras=None, st='U', bERKPP=None, bgrb=2),
         *par['SOS_bind_SHCP_GRB2'])

    # SOS also binds GAP-GRB2
    Rule("GAP_GRB2_bind_SOS",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=None, bcpp=None) + SOS(bras=None, bgrb=None, bERKPP=None, st='U') |
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, st='U', bERKPP=None),
         *par['SOS_bind_GAP_GRB2'])

    # GAP-GRB2-SOS and GAP-SHC:P-GRB2-SOS catalyze RAS-GDP->RAS-GTP:
    Rule("GAP_GRB2_SOS_bind_RASGDP",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None) |
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGDP_bind_bound_GRB2_SOS'])

    Rule("GAP_SHCP_GRB2_SOS_bind_RASGDP",
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None) |
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGDP_bind_bound_GRB2_SOS'])

    # Instead of a one-way catalytic process, the Chen-Sorger model implements this as a bidirectional process, as below:
    Rule('GAP_GRB2_SOS_bind_RASGTP',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='N', bpi3k=None) |
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_bind_bound_GRB2_SOS'])

    Rule('GAP_SHCP_GRB2_SOS_bind_RASGTP',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='N', bpi3k=None) |
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_bind_bound_GRB2_SOS'])

    # Recycling of activated RAS-GTP --> RAS-GDP.  In Chen/Sorger model, activated RAS-GTP is produced upon Raf phosphorylation.
    Rule('RASGTPact_bind_SOS_SHCP_complex',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='Y', bpi3k=None) |
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', act='N', bpi3k=None),
         *par['RASGTPact_bind_bound_GRB2_SOS'])

    Rule('RASGTPact_bind_SOS_GRB2_GAP_complex',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP', act='Y', bpi3k=None) |
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', act='N', bpi3k=None),
         *par['RASGTPact_bind_bound_GRB2_SOS'])

    Rule('RASGTP_unbind_SOS_GRB2_SHCP_complex',
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', act='N', bpi3k=None) |
         SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_unbind_GRB2_SOS'])

    Rule('RASGTP_unbind_SOS_GRB2_GAP_complex',
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GTP', act='N', bpi3k=None) |
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP', act='N', bpi3k=None),
         *par['RASGTP_unbind_GRB2_SOS'])

    # Activation of RAF -> RAF:P by RAS-GTP
    Rule('RASGTP_bind_RAF',
         RAS(bsos=None, bpi3k=None, st='GTP', act='N', braf=None) + RAF(st='U', ser295='U', b=None) |
         RAS(bsos=None, bpi3k=None, st='GTP', act='N', braf=1) % RAF(st='U', ser295='U', b=1),
         *par['RASGTP_bind_RAF'])

    Rule('RASGTP_RAF_cat',
         RAS(bsos=None, bpi3k=None, st='GTP', act='Y', braf=None) + RAF(st='P', ser295='U', b=None) |
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
    Initial(AKT(bpip3=None, both=None, S='P'), AKTP_0)
    
def akt_events():
    #GRB2 binds GAP-complex (without requiring SHC bound to complex):
    #Bind GRB2 without SOS already bound (two Chen-Sorger rate constants for different receptor dimers):
    #Rate 1: ErbB1/ErbB1 dimers (endosomal and plasma membrane), ErbB2/ErbB2 dimers (endosomal and plasma membrane), ErbB2/ErbB3 dimers (plasma membrane), and ErbB2/ErbB4 dimers (endosomal membrane):
    Rule('GRB2_bind_GAP_2',
         erbb(bd=1, ty='1') % erbb(bd=1, ty='1') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) |
         erbb(bd=1, ty='1') % erbb(bd=1, ty='1') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP_2'])

    Rule('GRB2_bind_GAP_3',
         erbb(bd=1, ty='2') % erbb(bd=1, ty='2') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) |
         erbb(bd=1, ty='2') % erbb(bd=1, ty='2') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP_2'])
    
    Rule('GRB2_bind_GAP_4',
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='3', loc='C') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) |
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='3', loc='C') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP_2'])

    Rule('GRB2_bind_GAP_5',
          erbb(bd=1, ty='4', loc='E') % erbb(bd=1, ty='2', loc='E') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) |
          erbb(bd=1, ty='4', loc='E') % erbb(bd=1, ty='2', loc='E') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
          *par['GRB2_bind_GAP_2'])

    #Rate 2: ErbB1/ErbBX, X=2, 3, 4 (endosomal and plasma membrane), ErbB2/ErbB3 dimers (endosomal membrane), and ErbB2/ErbB4 dimers (plasma membrane):
    for i in ['2', '3', '4']:
        Rule('GRB2_bind_GAP_6_'+i,
        erbb(bd=1, ty='1') % erbb(bd=1, ty=i) % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) |
        erbb(bd=1, ty='1') % erbb(bd=1, ty=i) % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
        *par['GRB2_bind_GAP'])

    Rule('GRB2_bind_GAP_7',
         erbb(bd=1, ty='2', loc='E') % erbb(bd=1, ty='3', loc='E') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) |
         erbb(bd=1, ty='2', loc='E') % erbb(bd=1, ty='3', loc='E') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP'])

    Rule('GRB2_bind_GAP_8',
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='4', loc='C') % GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=None) |
         erbb(bd=1, ty='2', loc='C') % erbb(bd=1, ty='4', loc='C') % GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=None, bgab1=None, bcpp=None, bgap=2),
         *par['GRB2_bind_GAP'])

    #Bind GRB2 to GAP with SOS already bound (one rate constant set for all dimer combinations):
    Rule('GRB2_bind_GAP_1',
         GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=1, bgab1=None, bcpp=None, bgap=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) |
         GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=1, bgab1=None, bcpp=None, bgap=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         *par['GRB2_SOS_bind_GAP'])

    #GAB1 binds GAP-GRB2. Specify plasma membrane complexes in order to prevent complex building on endosomal receptors, so that degradation rxns (above in receptor events) can be simplified -- GAB1 complexes are not degraded as per Chen/Sorger model 

    Rule('GRB2_bind_GAB1',
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None, bcpp=None) + GAB1(bgrb2=None, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U') |
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=2, bcpp=None) % GAB1(bgrb2=2, bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'),
         *par['GRB2_bind_GAB1'])
    
    #GAP-GRB2-GAB1 phosphorylation - Rates from Table p. 5 Chen et al 2009
    Rule('GAB1_bind_ATP',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') + ATP(b=None) |
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1),
         *par['GAB1_bind_ATP'])

    Rule('GAB1_phos',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1) >>
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + ADP(),
         par['GAB1_phos'])

    #SHP2 can desphosphorylate GAB1-P
    catalyze_state(SHP2(), 'bgab1', GAB1(bgrb2=ANY, bpi3k=None, batp=None, bERKPP=None, bPase9t=None), 'bshp2', 'S', 'P', 'U', (par['SHP2_dephos_GAB1P']))
    
    # MODFICIATION for Rexer model: Added binding of inhibitor BKM120 to PI3K
    # FURTHER MODIFICATION: Added 'mutated' PI3K (PI3K that is catalytically active regardless of what it's bound to)

    #catalyze_state(PI3K, 'bpip', PIP(both=None, bakt=None, bself2=None), 'bpi3k_self', 'S', 'PIP2', 'PIP3', par['PIP2_PI3Kmut_catalysis'])

    #bind(PI3K(), 'bpip', BKM120(loc='C'), 'b', par['BKM120_bind_PI3K'])
    
    #After GAB1 phosphorylation, all receptor dimer combinations can bind a single PI3K
    #Chen/Sorger model gives two rate constant sets for different receptor dimers:
    #Rate 1: ErbB1/ErbB1, ErbB1/ErbB2, ErbB1/ErbB4, and ErbB2/ErbB4 dimers:
    for i in ['1', '2', '4']:
        Rule('GAB1_bind_PI3K_1_'+i,
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) |
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_1'])

    Rule('GAB1_bind_PI3K_2',
         erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='4') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) |
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='4') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_1'])

    #Rate 2: ErbB1/ErbB3, ErbB2/ErbB2, and ErbB2/ErbB3 dimers:
    for i in ['1', '2']:
        Rule('GAB1_bind_PI3K_3_'+i,
             erbb(bd=ANY, ty=i) % erbb(bd=ANY, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) |
             erbb(bd=ANY, ty=i) % erbb(bd=ANY, ty='3') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_2'])

    Rule('GAB1_bind_PI3K_4',
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='2') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) |
             erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty='2') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
             *par['GAB1_bind_PI3K_2'])
    
    #PI3K bound to complex catalyzes PIP2 -> PIP3
    #Two rate sets for initial binding in Chen/Sorger model:
    #Rate 1: ErbB1/ErbBX dimers:
    for i in ['1', '2', '3', '4']:
        Rule('PIP2_bind_PI3K_1_'+i,
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=None) |
             erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=1, bgab1=ANY, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=1),
             *par['PIP2_bind_PI3K_1'])

    #Rate 2: ErbB2/ErbBX dimers, X=2, 3, 4:
    #FIXME: What is up with v701 in reaction list?
    for i in ['2', '3', '4']:
        Rule('PIP2_bind_PI3K_2_'+i,
        erbb(bd=ANY, ty='2') % erbb(bd=ANY, ty=i) % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=ANY, bras=None) + PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=None) |
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


    #ErbB2/ErbB3 dimers contain 6 separate PI3K binding sites.  Instead of representing this as 6 separate sites on ErbB3 (which would create a lot of extra species), it is represented as 1 virtual site with a catalytic rate 6 times faster than the normal rate.

    Rule('ErbB23_bind_PI3K',
         erbb(bd=1, ty='2', loc='C', st='P', b=None) % erbb(bd=1, ty='3', loc='C', st='P', b=None) + PI3K(bpip=None, bgab1=None, bras=None) |
         erbb(bd=1, ty='2', loc='C', st='P', b=None) % erbb(bd=1, ty='3', loc='C', st='P', b=2) % PI3K(bpip=None, bgab1=2, bras=None),
         *par['ErbB23_bind_PI3K'])

    Rule('ErbB23_PI3K_bind_PIP2',
         erbb(bd=1, ty='2', loc='C', st='P', b=None) % erbb(bd=1, ty='3', loc='C', st='P', b=2) % PI3K(bpip=None, bgab1=2, bras=None) + PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=None) |
         erbb(bd=1, ty='2', loc='C', st='P', b=None) % erbb(bd=1, ty='3', loc='C', st='P', b=2) % PI3K(bpip=3, bgab1=2, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=3),
         *par['ErbB23_PI3K_bind_PIP2'])

    Rule('ErbB23_PI3K_cat_PIP3',
         erbb(bd=1, ty='2', loc='C', st='P', b=None) % erbb(bd=1, ty='3', loc='C', st='P', b=2) % PI3K(bpip=3, bgab1=2, bras=None) % PIP(S='PIP2', both=None, bakt=None, bself2=None, bpi3k_self=3) >>
         erbb(bd=1, ty='2', loc='C', st='P', b=None) % erbb(bd=1, ty='3', loc='C', st='P', b=2) % PI3K(bpip=None, bgab1=2, bras=None) + PIP(S='PIP3', both=None, bakt=None, bself2=None, bpi3k_self=None),
         par['ErbB23_PI3K_cat_PIP3'])
             
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

    # MODIFICATION for Rexer model: AKT:P -> AKT:PP phosphorylation is done by mTORC2, not by PDK1 (see downstream_signaling_events module below)

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
    catalyze_state(ERK(st='PP', loc='C'), 'b', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None, bpi3k2=None, bpi3k3=None, bpi3k4=None, bpi3k5=None, bpi3k6=None), 'bERKPP', 'S', 'P', 'PP', (par['ERKPP_phos_GAB1P']))

    #GAP-GRB2-GAB1:P:P is dephosphorylated by Pase9t
    catalyze_state(Pase9t(), 'bgab1', GAB1(bgrb2=ANY), 'bPase9t', 'S', 'PP', 'P', (par['Pase9t_dephos_GAB1PP']))

    #ERK:P:P phosphorylates GRB2-SOS, preventing RAS-GDP->RAS-GTP conversion
    #To conform with Chen/Sorger model, this only effects ErbB1/ErbB1 dimers containing SOS and free SOS, and phosphorylated SOS can only bind ErbB1/ErbB1 complexes, not free GRB2:
    catalyze_state(ERK(st='PP', loc='C'), 'b', SOS(bgrb=None, bras=None), 'bERKPP', 'st', 'U', 'P', (par['ERKPP_phos_SOS']))

    Rule('ERKPP_bind_SOS_1',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=None, st='U') + ERK(st='PP', b=None, loc='C') |
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=1, st='U') % ERK(st='PP', b=1, loc='C'),
        *par['ERKPP_phos_SOS'][0:2])

    Rule('ERKPP_bind_SOS_2',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='U', bgrb=ANY) + ERK(st='PP', b=None, loc='C') |
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=1, st='U', bgrb=ANY) % ERK(st='PP', b=1, loc='C'),
        *par['ERKPP_phos_SOS'][0:2])

    Rule('ERKPP_phos_SOS_1',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=1, st='U') % ERK(st='PP', b=1, loc='C') >>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=ANY, bERKPP=None, st='P') + ERK(st='PP', b=None, loc='C'),
         par['ERKPP_phos_SOS'][2])

    Rule('ERKPP_phos_SOS_2',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=1, st='U', bgrb=ANY) % ERK(st='PP', b=1, loc='C') >>
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='P', bgrb=ANY) + ERK(st='PP', b=None, loc='C'),
         par['ERKPP_phos_SOS'][2])

    Rule('SOSP_bind_GRB2_1',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=None) + SOS(bras=None, bgrb=None, bERKPP=None, st='P') |
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=1, bgap=ANY, bgab1=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='P'),
         *par['SOSP_bind_GRB2'])

    Rule('SOSP_bind_GRB2_2',
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None, b=ANY) + SOS(bras=None, bERKPP=None, st='P', bgrb=None) |
         erbb(bd=ANY, ty='1') % erbb(bd=ANY, ty='1') % GAP(bd=ANY, b=ANY, bgrb2=None) % SHC(batp=None, st='P', bgrb=ANY, bgap=ANY) % GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=ANY) % SOS(bras=None, bERKPP=None, st='P', bgrb=1),
         *par['SOSP_bind_GRB2'])

    #AKT:P:P phosphorylates RAF:P at Ser295, preventing MEK phosphorylation.
    catalyze_state(AKT(S='PP', bpip3=None), 'both', RAF(st='P'), 'b', 'ser295', 'U', 'P', (par['AKTPP_phos_RAFP']))

    #RAS-GDP binds PI3K and PI3K catalyzes GDP --> GTP transformation. (Note: Matches model, but is this actually what's happening biologically?)
    #MODIFICATION for Rexer model: Reversed direction of this interaction (RAS now activates PI3K) - Castellano Genes and Cancer 2011
    #catalyze_state(PI3K(bgab1=ANY, bpip=None), 'bras', RAS(bsos=None, braf=None), 'bpi3k', 'st', 'GDP', 'GTP', (par['RASGDP_bind_PI3K']))
    bind(RAS(bsos=ANY, braf=None, st='GTP', act='N'), 'bpi3k', PI3K(bgab1=None, bpip=None), 'bras', par['RAS_bind_PI3K'])
    
    catalyze_state(PI3K(bgab1=None, bras=ANY), 'bpip', PIP(bakt=None, both=None), 'bpi3k_self', 'S', 'PIP2', 'PIP3', (par['RAS_cat_PI3K']))

#These are AKT and ERK downstream signaling events suggested as additions by Brent Rexer.

def downstream_signaling_monomers():
    Monomer('mTOR', ['bcomplex', 'bcat', 'S2448', 'bFKBP38'], {'S2448':['U','P']})
    Monomer('TORC1_ptns', ['bmTOR']) #A complex composed of Raptor, mLST8, PRAS40, and DEPTOR.  With mTOR becomes mTORC1 (mTOR complex 1)
    Monomer('TORC2_ptns', ['bmTOR']) #A complex composed of Rictor, GbetaL, and mSIN1, mLST8, and DEPTOR. With mTOR becomes mTORC2 (mTOR complex 2)
    Monomer('AMPK', ['T172', 'b'], {'T172':['U', 'P']}) #Phosphorylated at T172 by LKB1
    Monomer('LKB1', ['b', 'S431'], {'S431':['U', 'P']})
    Monomer('TSC', ['b', 'S664', 'S1798', 'S1387', 'S939', 'S981', 'T1462'], {'T1462': ['U', 'P'], 'S664':['U', 'P'], 'S1798':['U','P'], 'S1387':['U','P'], 'S939':['U', 'P'], 'S981':['U', 'P']}) #Composed of both TSC1 and TSC2
    Monomer('S6K', ['T252', 'T412', 'b'], {'T252':['U', 'P'], 'T412':['U','P']})
    Monomer('Rheb', ['bTSC', 'bFKBP38', 'S', 'bmTOR'], {'S':['GDP', 'GTP']})
    Monomer('FKBP38', ['b'])
    Monomer('rpS6', ['b', 'S'], {'S':['U', 'P']})
    Monomer('EIF4EBP1', ['bEIF4E', 'bmTOR', 'S'], {'S':['U','P']})
    Monomer('EIF4E', ['b'])
    Monomer('RSK1', ['b', 'T573', 'S380', 'S221'], {'T573':['U','P'], 'S380':['U','P'], 'S221':['U','P']})
    Monomer('ELK1', ['b', 'S383'], {'S383':['U','P']})
    alias_model_components()
    
def downstream_signaling_initial():
    Initial(mTOR(bcomplex=None, bcat=None, bFKBP38=None, S2448='U'), mTOR_0)
    Initial(TORC1_ptns(bmTOR=None), TORC1_ptns_0)
    Initial(TORC2_ptns(bmTOR=None), TORC2_ptns_0)
    Initial(AMPK(T172='U', b=None), AMPK_0)
    Initial(LKB1(b=None, S431='U'), LKB1_0)
    Initial(TSC(b=None, S664='U', S1798='U', S1387='U', S939='U', S981='U', T1462='U'), TSC_0)
    Initial(S6K(b=None, T252='U', T412='U'), S6K_0)
    Initial(Rheb(bTSC=None, bFKBP38=None, S='GTP', bmTOR=None), Rheb_0)
    Initial(FKBP38(b=None), FKBP38_0)
    Initial(rpS6(b=None, S='U'), rpS6_0)
    Initial(EIF4EBP1(bEIF4E=None, bmTOR=None, S='U'), EIF4EBP1_0)
    Initial(EIF4E(b=None), EIF4E_0)
    Initial(RSK1(b=None, T573='U', S380='U', S221='U'), RSK1_0)
    Initial(ELK1(b=None, S383='U'), ELK1_0)
    
def downstream_signaling_events():
    
    # mTOR can bind either the proteins in mTORC1 or mTORC2:
    bind(mTOR(bcat=None, S2448='U', bFKBP38=None), 'bcomplex', TORC1_ptns(), 'bmTOR', (par['mTOR_bind_TORC1ptns']))
    
    bind(mTOR(bcat=None, S2448='U', bFKBP38=None), 'bcomplex', TORC2_ptns(), 'bmTOR', (par['mTOR_bind_TORC2ptns']))

    # mTORC2 phosphorylates AKT:P -> AKT:PP
    bind_complex(TORC2_ptns(bmTOR=1) % mTOR(bcomplex=1), 'bcat', AKT(S='P', bpip3=None), 'both', par['mTORC2_bind_AKTP'])

    Rule('mTOR_AKT_cat',
         TORC2_ptns(bmTOR=1) % mTOR(bcomplex=1, bcat=2) % AKT(S='P', both=2, bpip3=None) >>
         TORC2_ptns(bmTOR=1) % mTOR(bcomplex=1, bcat=None) + AKT(S='PP', both=None, bpip3=None),
         par['mTORC2_cat_AKTPP'])

    # AKT:PP phosphorylates TSC2 at S939, S981, and T1462, which causes it to translocate from the membrane to the cytosol and stops TSC2's GAP activity on Rheb (Cai 2006, Journal of Cell Biology).
    # ERK:PP can also phosphorylate TSC2 at S664, again inhibiting its GAP activity.

    bind(AKT(S='PP', bpip3=None), 'both', TSC(S939='U', S981='U', T1462='U'), 'b', par['AKTPP_bind_TSC2'])
    
    for site in ['S939', 'S981', 'T1462']:
        Rule('AKTPP_phos_TSC2_'+site,
        AKT(S='PP', bpip3=None, both=1) % TSC({'b':1, site:'U'}) >>
        AKT(S='PP', bpip3=None, both=None) + TSC({'b':None, site:'P'}),
        par['AKTPP_phos'+site+'_TSC2'])

    # Rule('AKTPP_phos_TSC2_S939',
    #     AKT(S='PP', bpip3=None, both=1) % TSC(b=1, S939='U', S981='U', T1462='U') >>
    #     AKT(S='PP', bpip3=None, both=None) + TSC(b=None, S939='P', S981='U', T1462='U'),
    #     par['AKTPP_phosS939_TSC2'])

    # Rule('AKTPP_phos_TSC2_S981',
    #     AKT(S='PP', bpip3=None, both=1) % TSC(b=1, S939='U', S981='U', T1462='U') >>
    #     AKT(S='PP', bpip3=None, both=None) + TSC(b=None, S939='U', S981='P', T1462='U'),
    #     par['AKTPP_phosS981_TSC2'])

    # Rule('AKTPP_phos_TSC2_T1462',
    #     AKT(S='PP', bpip3=None, both=1) % TSC(b=1, S939='U', S981='U', T1462='U') >>
    #     AKT(S='PP', bpip3=None, both=None) + TSC(b=None, S939='U', S981='U', T1462='P'),
    #     par['AKTPP_phosT1462_TSC2'])

    catalyze_state(ERK(st='PP', loc='C'), 'b', TSC(S939='U', S981='U', T1462='U'), 'b', 'S664', 'U', 'P', (par['ERKPP_phos_TSC2']))

    # ERK:PP phosphorylates RSK1 at T573. RSK1 then autocatalyzes phosphorylation at S380, which allows binding of PDK1, which phosphorylates S221, giving fully active RSK1.
    catalyze_state(ERK(st='PP', loc='C'), 'b', RSK1(S380='U', S221='U'), 'b', 'T573', 'U', 'P', (par['ERKPP_phos_RSK1']))

    Rule('RSK1_autocatalysis',
         RSK1(S380='U', S221='U', T573='P') >>
         RSK1(S380='P', S221='U', T573='P'),
         par['RSK1_autocat'])

    catalyze_state(PDK1(bakt=None), 'both', RSK1(S380='P', T573='P'), 'b', 'S221', 'U', 'P', (par['PDK1_phos_RSK1']))

    # Active RSK1 can phosphorylate TSC2 at S1798, inhibiting its GAP activity.
    catalyze_state(RSK1(S380='P', T573='P', S221='P'), 'b', TSC(S939='U', S981='U', T1462='U'), 'b', 'S1798', 'U', 'P', (par['RSK1_phos_TSC2']))

    # Active RSK1 phosphorylates LKB1 at S431, activating LKB1.
    catalyze_state(RSK1(S380='P', T573='P', S221='P'), 'b', LKB1(), 'b', 'S431', 'U', 'P', (par['RSK1_phos_LKB1']))

    # Active LKB1 phosphorylates AMPK at T172, activating it.
    catalyze_state(LKB1(S431='P'), 'b', AMPK(), 'b', 'T172', 'U', 'P', (par['LKB1_phos_AMPK']))

    # Active AMPK phosphorylates S1387 on TSC2, activating its GAP activity.
    catalyze_state(AMPK(T172='P'), 'b', TSC(S939='U', S981='U', T1462='U'), 'b', 'S1387', 'U', 'P', (par['AMPK_phos_TSC2']))

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

    catalyze_state(PDK1(bakt=None), 'both', S6K(T412='P'), 'b', 'T252', 'U', 'P', (par['PDK1_phos_S6K']))

    # Active S6K phosphorylates rpS6.
    catalyze_state(S6K(T252='P', T412='P'), 'b', rpS6(), 'b', 'S', 'U', 'P', (par['S6K_phos_rpS6']))

    # Unphosphorylated EIF4EBP1 binds EIF4E (EIF4E is necessary for mRNA translation, which this binding interaction prevents).
    bind(EIF4EBP1(bmTOR=None, S='U'), 'bEIF4E', EIF4E(), 'b', par['EIF4EBP1_bind_EIF4E'])

    # Active mTORC1 phosphorylates EIF4EBP1, preventing its interaction with EIF4E and activating mRNA translation.
    catalyze_state(mTOR(bcomplex=ANY, bFKBP38=None, S2448='P'), 'bcat', EIF4EBP1(bEIF4E=None), 'bmTOR', 'S', 'U', 'P', (par['mTORC1_phos_EIF4EBP1']))

    # ERK:PP can translocate to the nucleus and activate ELK-1 by phosphorylating S383 and S389.
    equilibrate(ERK(st='PP', loc='C', b=None), ERK(st='PP', loc='N', b=None), par['ERKPP_to_nucleus'])
    
    catalyze_state(ERK(st='PP', loc='N'), 'b', ELK1(), 'b', 'S383', 'U', 'P', (par['ERKPP_phos_ELK1']))
