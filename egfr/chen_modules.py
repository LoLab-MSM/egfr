"""
Overview
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

# Monomer declarations
# ====================

def rec_monomers():

    """ Declares the ErbB receptor interactions.
    'bf' is the default site to be used for all binding/catalysis reactions.
    """
    Monomer('EGF', ['b']) # Epidermal Growth Factor ligand
    Monomer('HRG', ['b']) # Heregulin ligand
    Monomer('erbb', ['bl', 'bd', 'b', 'ty', 'st', 'loc', 'pi3k1', 'pi3k2', 'pi3k3', 'pi3k4', 'pi3k5', 'pi3k6', 'cpp'], {'ty':['1','2','3','4'], 'st':['U','P'], 'loc':['C','E'], 'cpp':['Y', 'N']}) # bl: lig, bd: dimer, b: binding, ty: rec type, st: (U)n(P)hosphorylated, loc: (C)yto 'brane or (E)ndosome 'brane, cpp: No real biophysical meaning; useful model marker for presence of CPP bound downstream.

    Monomer('DEP', ['b'])
    Monomer('ATP', ['b'])
    Monomer('ADP')
    Monomer('CPP', ['b', 'loc'], {'loc':['C', 'E']})

def rec_initial():
    # # Initial concentrations (except DEP1) for all cell types taken from Chen et al 2009 -- see Jacobian files
    Parameter('EGF_0',      3.01e12) # c1 5 nM EGF = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell 
    Parameter('HRG_0',         0) # c514 5 nM HRG = 3.01e12 molec/cell; .01 nM EGF = 6.02e9 molec/cell

    # # These values are for A431 cells
    Parameter('erbb1_0',  1.08e6) #c531
    Parameter('erbb2_0',  4.62e5) #c141
    Parameter('erbb3_0',  6.23e3) #c140
    Parameter('erbb4_0',  7.94e2) #c143
    Parameter('ATP_0',     1.2e9) #c105
    Parameter('DEP_0',       7e4) #c280
    Parameter('CPP_0', 5e3) #c12

    # These values are for H1666 cells
    # Parameter('erbb1_0',  1.60e5) #c531
    # Parameter('erbb2_0',  6.83e3) #c141
    # Parameter('erbb3_0',  6.05e3) #c140
    # Parameter('erbb4_0',  2.59e1) #c143
    # Parameter('ATP_0', 1.2e9) #c105
    # Parameter('DEP_0', 6.23876e6) #c280
    # Parameter('CPP_0', 1.59621e6) #c12

    # These values are for H3255 cells
    # Parameter('erbb1_0',  1.29e6) #c531
    # Parameter('erbb2_0',  3.16e4) #c141
    # Parameter('erbb3_0',  4.48e4) #c140
    # Parameter('erbb4_0',  2.58e1) #c143
    # Parameter('ATP_0',     1.2e9) #c105
    # Parameter('DEP_0',       1.2448e9) #c280
    # Parameter('CPP_0', 4.49873e6) #c12

    alias_model_components()

    Initial(EGF(b=None), EGF_0)
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
    # EGF / HRG receptor binding rates obtained from Chen et al (Supplementary)
    bind_table([[                                                          EGF,             HRG],
                [erbb(ty='1', bd=None, b=None, st='U', loc='C'),   (1e7, 3.3e-2),          None],
                [erbb(ty='3', bd=None, b=None, st='U', loc='C'),            None,   (1e7, 7e-2)],
                [erbb(ty='4', bd=None, b=None, st='U', loc='C'),            None,   (1e7, 7e-2)]],
                'bl', 'b')
    
    # ErbB dimerization
    # Dimerization rates obtained from Chen et al (Supplementary)
    # erbb's is required to containe a ligand (except for erbb2 which cannot bind a ligand)
    erbb1Lig = erbb(ty='1', bl=ANY, b=None, st='U', loc='C')
    erbb2Lig = erbb(ty='2', bl=None, b=None, st='U', loc='C')
    erbb3Lig = erbb(ty='3',bl=ANY, b=None, st='U', loc='C')
    erbb4Lig = erbb(ty='4',bl=ANY, b=None, st='U', loc='C')
    bind_table([[                          erbb1Lig,            erbb2Lig, erbb3Lig, erbb4Lig],
                [erbb1Lig,        (7.45e-6, 1.6e-1),                None,     None,     None],
                [erbb2Lig,        (3.74e-8, 1.6e-2),  (1.67e-10, 1.6e-2),     None,     None],
                [erbb3Lig,        (3.74e-8, 1.6e-2),  (1.67e-10, 1.6e-2),     None,     None],
                [erbb4Lig,        (3.74e-8, 1.6e-2),  (1.67e-10, 1.6e-2),     None,     None]],
        'bd', 'bd')

    # ATP binding: ATP only binds to dimers
    # ATP binding rates obtained from Chen et al (Supplementary)
    # include DEP binding here since they both bind to the same site
    bind_table([[                                                ATP,  DEP],
                [erbb(ty='1', bd=ANY, st='U', loc='C'), (1.87e-8, 1), (5e-5, 1e-2)],
                [erbb(ty='2', bd=ANY, st='U', loc='C'), (1.87e-8, 1), (5e-5, 1e-2)],
                [erbb(ty='4', bd=ANY, st='U', loc='C'), (1.87e-8, 1), (5e-5, 1e-2)]],
        'b', 'b')

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
                 ATP(b=1) % erbb(ty=i, b=1,    bd=2, st='U') % erbb(ty=j, bd=2, st='U') >>
                 ADP()    + erbb(ty=i, b=None, bd=2, st='P') % erbb(ty=j, bd=2, st='P'),
                 Parameter("kcp"+i+j, 1e-1))
            Rule("cross_DEphospho_"+i+"_"+j,
                 DEP(b=1)   %  erbb(ty=i, b=1,    bd=2, st='P') % erbb(ty=j, bd=2, st='P') >>
                 DEP(b=None) + erbb(ty=i, b=None, bd=2, st='U') % erbb(ty=j, bd=2, st='U'),
                 Parameter("kcd"+i+j, 1e-1))

    # Receptor internalization
    # This internalizes all receptor combos after binding to CPP (coated pit protein).
    # FIXME: this should just be a state transformation (i.e. the rule looks noisy)
    # Internalization rates taken from Chen et al Table I pg. 5
    Rule("rec_intern",
         erbb(bd=1, loc='C', cpp='N') % erbb(bd=1, loc='C', cpp='N') <> erbb(bd=1, loc='E', cpp='N') % erbb(bd=1, loc='E', cpp='N'),
         Parameter("kintf", 1.3e-3), Parameter("kintr", 5e-5))

    # CPP bound to receptors can catalyze their internalization (when they are bound to any complex containing GRB2, except GAB1 complex):
    Rule('CPP_bind_GAP_GRB2',
         CPP(loc='C', b=None) + erbb(bd=1, loc='C', cpp='N') % erbb(bd=1, loc='C', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bgab1=None, b=None, bcpp=None) <>
         erbb(bd=1, loc='C', cpp='Y') % erbb(bd=1, loc='C', cpp='Y') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3, bgab1=None, b=None) % CPP(loc='C', b=3),
         Parameter('cpp_gap_grb2_bind_f', KF),
         Parameter('cpp_gap_grb2_bind_r', KR))

    Rule('CPP_bind_SHC_GRB2',
         erbb(bd=1, loc='C', cpp='N') % erbb(bd=1, loc='C', cpp='N') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None, bgap=None) + CPP(loc='C', b=None) <>
         erbb(bd=1, loc='C', cpp='Y') % erbb(bd=1, loc='C', cpp='Y') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None, bgap=None) % CPP(loc='C', b=4),
         Parameter('cpp_shc_grb2_bind_f', KF),
         Parameter('cpp_shc_grb2_bind_r', KR))

    Rule("CPP_rec_int_1",
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3) % CPP(loc='C', b=3) <>
         erbb(bd=1, loc='E') % erbb(bd=1, loc='E') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3) % CPP(loc='E', b=3),
         Parameter('kcppintf_1', 1.3e-3),
         Parameter('kcppintr_1', 5e-5))

    Rule("CPP_rec_int_2",
         erbb(bd=1, loc='C') % erbb(bd=1, loc='C') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4) % CPP(loc='C', b=4) <>
         erbb(bd=1, loc='E') % erbb(bd=1, loc='E') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4) % CPP(loc='E', b=4),
         Parameter('kcppintf_2', 1.3e-3),
         Parameter('kcppintr_2', 5e-5))

    Rule('CPPE_bind_GAP_GRB2',
         erbb(bd=1, loc='E', cpp='N') % erbb(bd=1, loc='E', cpp='N') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=None, bgab1=None) + CPP(loc='E', b=None) <>
         erbb(bd=1, loc='E', cpp='Y') % erbb(bd=1, loc='E', cpp='Y') % GAP(bd=ANY, bgrb2=2) % GRB2(bgap=2, bcpp=3, bgab1=None) % CPP(loc='E', b=3),
         Parameter('cppe_gap_grb2_bind_f', KF),
         Parameter('cppe_gap_grb2_bind_r', KR))

    Rule('CPPE_bind_SHC_GRB2',
         erbb(bd=1, loc='E', cpp='N') % erbb(bd=1, loc='E', cpp='N') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=None, bgab1=None) + CPP(loc='E', b=None) <>
         erbb(bd=1, loc='E', cpp='Y') % erbb(bd=1, loc='E', cpp='Y') % GAP(bd=ANY, b=2) % SHC(bgap=2, batp=None, st='P', bgrb=3) % GRB2(b=3, bcpp=4, bgab1=None) % CPP(loc='E', b=4),
         Parameter('cppe_shc_grb2_bind_f', KF),
         Parameter('cppe_shc_grb2_bind_r', KR))
    
    Rule("CPP_intern",
         CPP(loc='E', b=None) <> CPP(loc='C', b=None),
         Parameter('cppintf', 1.3e-3),
         Parameter('cppintr', 5e-5))
         
    # Receptor degradation
    # This degrades all receptor combos within an endosome
    # Should this have a forward and reverse rate - k62b
    degrade(erbb(bd=1, loc='E') % erbb(bd=1, loc='E'), Parameter("kdeg", 4.16e-4))

def mapk_monomers():
    Monomer('GAP', ['bd', 'b', 'bgrb2'])
    Monomer('SHC', ['bgap', 'bgrb', 'batp', 'st'], {'st':['U','P']})
    # Monomer('SHCPase', ['b'])
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

    # These values are for A431 cells (obtained from Chen et al 2009 Jacobian files) -- values prefixed by 'c' are Chen/Sorger variable names:
    Parameter('GAP_0', 5.35e5) #c14
    Parameter('SHC_0', 1.1e6) #c31
    # Parameter('SHCPase_0', 1000)
    Parameter('GRB2_0', 1264) #c22
    #Parameter('SOS_0', 6.63e4) removed to better represent Chen/Sorger model
    Parameter('RAS_0', 5.81e4) #c26
    Parameter('RAF_0', 7.11e4) #c41
    Parameter('MEK_0', 3.02e6) #c47
    Parameter('ERK_0', 6.95e5) #c55
    Parameter('PP1_0', 5e4) #c44
    Parameter('PP2_0', 1.25e5) #c53
    Parameter('PP3_0', 1.69e4) #c60
    Parameter('GRB2_SOS_0', 8.89e7) #This added to better represent Chen Sorger model c30

    # These values are for H1666 cells:
    # Parameter('GAP_0', 1.06697e9) #c14
    # Parameter('SHC_0', 1.1e6) #c31
    # # Parameter('SHCPase_0', 1000)
    # Parameter('GRB2_0', 12649.1) #c22
    # #Parameter('SOS_0', 6.63e4) removed to better represent Chen/Sorger model
    # Parameter('RAS_0', 183713) #c26
    # Parameter('RAF_0', 7113.12) #c41
    # Parameter('MEK_0', 3.02e6) #c47
    # Parameter('ERK_0', 6.95e5) #c55
    # Parameter('PP1_0', 28117.1) #c44
    # Parameter('PP2_0', 39363.9) #c53
    # Parameter('PP3_0', 168702) #c60
    # Parameter('GRB2_SOS_0', 2.81171e6) #This added to better represent Chen Sorger model c30

    #These values are for H3255 cells:
    # Parameter('GAP_0', 3.37405e6) #c14
    # Parameter('SHC_0', 1.1e6) #c31
    # # Parameter('SHCPase_0', 1000)
    # Parameter('GRB2_0', 400) #c22
    # #Parameter('SOS_0', 6.63e4) removed to better represent Chen/Sorger model
    # Parameter('RAS_0', 58095.2) #c26
    # Parameter('RAF_0', 1.26491e6) #c41
    # Parameter('MEK_0', 3.02e6) #c47
    # Parameter('ERK_0', 6.95e5) #c55
    # Parameter('PP1_0', 28117.1) #c44
    # Parameter('PP2_0', 39363.9) #c53
    # Parameter('PP3_0', 5.33484e6) #c60
    # Parameter('GRB2_SOS_0', 5e7) #This added to better represent Chen Sorger model c30
    
    alias_model_components()

    Initial(GAP(bd=None, b=None, bgrb2=None), GAP_0)
    Initial(SHC(bgap=None, bgrb=None, batp=None, st='U'), SHC_0)
    Initial(GRB2(b=None, bsos=None, bgap=None, bgab1=None, bcpp=None), GRB2_0)
    #Initial(SOS(bgrb=None, bras=None, bERKPP=None, st='U'), SOS_0)
    Initial(RAS(bsos=None, braf=None, bpi3k=None, st='GDP'), RAS_0)
    Initial(RAF(b=None, st='U', ser295='U'), RAF_0)
    Initial(MEK(b=None, st='U'), MEK_0)
    Initial(ERK(b=None, st='U'), ERK_0)
    Initial(PP1(b=None), PP1_0)
    Initial(PP2(b=None), PP2_0)
    Initial(PP3(b=None), PP3_0)
    Initial(GRB2(b=None, bsos=1, bgap=None, bgab1=None, bcpp=None) % SOS(bgrb=1, bras=None, bERKPP=None, st='U'), GRB2_SOS_0)

    
def mapk_events():

    # =====================
    # Alias model components for names in present namespace
    alias_model_components()

    # GAP binds to phosphorylated dimers
    # in the present we use MatchOnce to insure correct representation of the binding
    # similar to Chen et al
    Rule("GAP_binding",
         MatchOnce(erbb(bd=1, b=None, st='P') % erbb(bd=1, b=None, st='P')) + GAP(bd=None, b=None, bgrb2=None) <>
         MatchOnce(erbb(bd=1, b=2,    st='P') % erbb(bd=1, b=None, st='P')  % GAP(bd=2, b=None, bgrb2=None)),
         Parameter("kerbb_dim_GAPf", KF), Parameter("kerbb_dim_GAPr", KR))
    
    # SHC binds to GAP-complex
    bind(GAP(bd=ANY, bgrb2=None), 'b', SHC(batp=None, st='U'), 'bgap', [KF, KR])

    # SHC phosphorylation
    Rule("Shc_bind_ATP",
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=None, st='U') + ATP(b=None) <>
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=2, st='U') % ATP(b=2),
         Parameter("ShcATPf",KF), Parameter("ShcATPr",KR))
    
    Rule("Shc_phos",
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=2, st='U') % ATP(b=2) >>
         GAP(bd=ANY, b=1) % SHC(bgap=1, bgrb=None, batp=None, st='P') + ADP(),
         Parameter("ShcPhosc", KCP))
    
    # GRB2 binds to GAP-SHC:P with or without SOS:
    Rule('GRB2_bind_GAP_SHCP_1',
         SHC(batp=None, st='P', bgrb=None) + GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) <>
         SHC(batp=None, st='P', bgrb=2) % GRB2(bgap=None, bgab1=None, bsos=1, bcpp=None, b=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         Parameter('GRB2_bind_GAP_SHCP_f', KF),
         Parameter('GRB2_bind_GAP_SHCP_r', KR))

    bind(SHC(batp=None, st='P'), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=None, bcpp=None), 'b', [KF,KR])
    
    # Can use this simpler version if GRB2-SOS complex isn't present alone:
    #    bind(SHC(batp=None, st='P'), 'bgrb', GRB2(bgap=None, bgab1=None, bsos=ANY, bcpp=None), 'b', [KF, KR])

    # GRB2 and SOS bind/disassociate:
    bind(GRB2(bgap=None, bgab1=None, b=None, bcpp=None), 'bsos', SOS(bras=None, bERKPP=None, st='U'), 'bgrb', [7.5e-6, 1.5])

    #Although no free SOS is present initially in Chen Sorger model, GRB2-SOS can disassociate (see above), so these are necessary.
    # SOS binds to GAP-SHC:P-GRB2 - rate obtained from Chen et al. pg 5. 
    bind(GRB2(b=ANY, bgap=None, bgab1=None, bcpp=None), 'bsos', SOS(bras=None, st='U', bERKPP=None), 'bgrb', [7.5e-6, 1.5])

    # SOS also binds GAP-GRB2
    Rule("GAP_GRB2_bind_SOS",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=None, bcpp=None) + SOS(bras=None, bgrb=None, bERKPP=None, st='U') <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, st='U', bERKPP=None),
         Parameter("GAP_GRB2_bind_SOSf", 7.5e-6),
         Parameter("GAP_GRB2_bind_SOSr", 1.5))

    # GAP-GRB2-SOS and GAP-SHC:P-GRB2-SOS catalyze RAS-GDP->RAS-GTP:
    Rule("GAP_GRB2_SOS_bind_RASGDP",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP') <>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP'),
         Parameter("GAP_GRB2_SOS_bind_RASf", KF),
         Parameter("GAP_GRB2_SOS_bind_RASr", KR))

    Rule("GAP_SHCP_GRB2_SOS_bind_RASGDP",
         GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GDP') <>
         GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP'),
         Parameter("GAP_SHCP_GRB2_SOS_bind_RASGDPf", KF),
         Parameter("GAP_SHCP_GRB2_SOS_bind_RASGDPr", KR))

    Rule("GAP_GRB2_SOS_catRAS",
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP') >>
         GRB2(bgap=ANY, bgab1=None, b=None, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP'),
         Parameter("GAP_GRB2_SOS_catRASc", KCP))

    Rule("GAP_SHCP_GRB2_SOS_catRAS",
         GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=2, bgrb=1, bERKPP=None, st='U') % RAS(braf=None, bsos=2, st='GDP') >>
         GRB2(bgap=None, bgab1=None, b=ANY, bsos=1, bcpp=None) % SOS(bras=None, bgrb=1, bERKPP=None, st='U') + RAS(braf=None, bsos=None, st='GTP'),
         Parameter("GAP_SHCP_GRB2_SOS_catRASc", KCP))
         
    #can use this simpler implementation of above if GRB2-SOS isn't present on its own as in Chen Sorger model:
    #catalyze_state(SOS(bgrb=ANY, st='U', bERKPP=None), 'bras', RAS(braf=None), 'bsos', 'st', 'GDP', 'GTP', (KF, KR, KCD))

    # Activation of RAF -> RAF:P by RAS-GTP 
    catalyze(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', RAF(st='U', ser295='U'), 'b', RAF(st='P', ser295='U'),
             (KF,KR,KCP))

    # Deactivation of RAF:P -> RAF by PP1
    catalyze(PP1(), 'b', RAF(st='P', ser295='U'), 'b', RAF(st='U', ser295='U'),
             (KF,KR,KCD))

    # Activation of MEK -> MEK:P by activated RAF
    catalyze(RAF(st='P', ser295='U'), 'b', MEK(st='U'), 'b', MEK(st='P'),
             (KF,KR,KCP))

    # Deactivation of MEK:P -> MEK by PP2
    catalyze(PP2(), 'b', MEK(st='P'), 'b', MEK(st='U'),
             (KF,KR,KCD))
    
    # Activation of MEK:P -> MEK:P:P by activated RAF
    catalyze(RAF(st='P', ser295='U'), 'b', MEK(st='P'), 'b', MEK(st='PP'),
             (KF,KR,KCP))

    # Deactivation of MEK:P:P -> MEK:P by PP2
    catalyze(PP2(), 'b', MEK(st='PP'), 'b', MEK(st='P'),
             (KF,KR,KCD))
    
    # Activation of ERK -> ERK:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='U'), 'b', ERK(st='P'),
             (KF,KR,KCP))

    # Deactivation of ERK:P -> ERK by PP3
    catalyze(PP3(), 'b', ERK(st='P'), 'b', ERK(st='U'),
             (KF,KR,KCD))

    # Activation of ERK:P -> ERK:P:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='P'), 'b', ERK(st='PP'),
             (KF,KR,KCP))

    # Deactivation of ERK:P:P -> ERK:P by PP3
    catalyze(PP3(), 'b', ERK(st='PP'), 'b', ERK(st='P'),
             (KF,KR,KCD))

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
    # Initial concentrations from Chen Sorger Jacobian files - 'c' prefixed terms are their variable names
    # These values are for A431 cells:
    Parameter('GAB1_0', 94868.3) #c426
    Parameter('PI3K_0', 3.55656e7) #c287 c455?
    Parameter('SHP2_0', 1e6) #c463
    Parameter('PIP_0',     3.94e5) #c444
    Parameter('PTEN_0',    5.62e4) #c279
    Parameter('SHP_0',     2.21e3) #c461
    Parameter('AKT_0',     9.05e5) #c107
    Parameter('PDK1_0',     3.00416e8) #c109
    Parameter('PP2A_III_0', 4.5e5) #c113

    # These values are for H1666 cells:
    # Parameter('GAB1_0', 5.33484e7) #c426
    # Parameter('PI3K_0', 3.55656e5) #c287 c455?
    # Parameter('SHP2_0', 1.12202e5) #c463
    # Parameter('PIP_0',     1.2448e6) #c444
    # Parameter('PTEN_0',    5000) #c279
    # Parameter('SHP_0',     7000) #c461
    # Parameter('AKT_0',     9.05e5) #c107
    # Parameter('PDK1_0',     1.8955e6) #c109
    # Parameter('PP2A_III_0', 401063) #c113

    # These values are for H3255 cells:
    # Parameter('GAB1_0', 30000) #c426
    # Parameter('PI3K_0', 2e9) #c287 c455?
    # Parameter('SHP2_0', 3.16228e6) #c463
    # Parameter('PIP_0',     700000) #c444
    # Parameter('PTEN_0',    158114) #c279
    # Parameter('SHP_0',     700) #c461
    # Parameter('AKT_0',     9.05e5) #c107
    # Parameter('PDK1_0',     3.00416e8) #c109
    # Parameter('PP2A_III_0', 2.53054e7) #c113
    
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
def akt_events():
    #GRB2 binds GAP-complex (without requiring SHC bound to complex and with or without SOS already bound) k16, kd24
    Rule('GRB2_bind_GAP_1',
         GAP(bd=ANY, b=None, bgrb2=None) + GRB2(b=None, bsos=1, bgab1=None, bcpp=None, bgap=None) % SOS(bras=None, bERKPP=None, st='U', bgrb=1) <>
         GAP(bd=ANY, b=None, bgrb2=2) % GRB2(b=None, bsos=1, bgab1=None, bcpp=None, bgap=2) % SOS(bras=None, bERKPP=None, st='U', bgrb=1),
         Parameter('GRB2_bind_GAP_1_f', KF),
         Parameter('GRB2_bind_GAP_1_r', KR))

    bind(GRB2(b=None, bsos=None, bgab1=None, bcpp=None), 'bgap', GAP(bd=ANY, b=None), 'bgrb2', [KF, KR])

    #Can use the simpler version below if GRB2-SOS isn't present by itself:
    #bind(GRB2(b=None, bsos=WILD, bgab1=None, bcpp=None), 'bgap', GAP(bd=ANY, b=None), 'bgrb2', [KF, KR])

    #GAB1 binds GAP-GRB2 k105, kd105
    bind(GRB2(b=None, bsos=None, bgap=ANY, bcpp=None), 'bgab1', GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, S='U'), 'bgrb2', [KF, KR])
    #GAP-GRB2-GAB1 phosphorylation - Rates from Table p. 5 Chen et al 2009
    Rule('GAB1_bind_ATP',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') + ATP(b=None) <>
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1),
         Parameter('GAB1ATPf', KF),
         Parameter('GAB1ATPr', KR))

    Rule('GAB1_phos',
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=1, bERKPP=None, bPase9t=None, bgrb2=ANY, S='U') % ATP(b=1) >>
         GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + ADP(),
         Parameter('GAB1Phosc', KCP))

    #SHP2 can desphosphorylate GAB1-P
    catalyze_state(SHP2(), 'bgab1', GAB1(bgrb2=ANY, bpi3k=None, batp=None, bERKPP=None, bPase9t=None), 'bshp2', 'S', 'P', 'U', (KF, KR, KCD))
   
    #After GAB1 phosphorylation, all receptor dimer combinations can bind a single PI3K
    Rule('GAB1_bind_PI3K_1',
         GAB1(bshp2=None, bpi3k=None, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') + PI3K(bpip=None, bgab1=None, bras=None) <>
         GAB1(bshp2=None, bpi3k=1, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P') % PI3K(bpip=None, bgab1=1, bras=None),
         Parameter('GAB1PI3Kf', KF),
         Parameter('GAB1PI3Kr', KR))

    #GAB1-PI3K bound to complex containing ErbB2/ErbB3 binds 1-6 PIP2 (creates chains; doesn't necessarily represent biology but accurately reproduces Chen Sorger 2009 model).
    #First bind a single PIP2 to PI3K complex - this rule created by catalyze_state below

    #Then create chains of up to 6 PIP2 molecules attached to a single PI3K:
    assemble_chain_sequential_base(PI3K(berb=None, bras=None, bgab1=ANY, bpip=None), 'bpip', PIP(bakt=None, both=None, S='PIP2'), 'bpi3k_self', 'bself2', 6, [[1e-5, 1e-1]]*5, \
                                   erbb(bd=ANY, ty='3', b=ANY, loc='C') % erbb(bd=ANY, ty='2', b=None, loc='C') % GAP(bd=ANY, b=None, bgrb2=ANY) % GRB2(b=None, bsos=None, bgap=ANY, bgab1=ANY) % GAB1(bshp2=None, bpi3k=ANY, batp=None, bERKPP=None, bPase9t=None, bgrb2=ANY, S='P'))
    
    #To accurately reproduce Chen Sorger model, allow PIP2 to catalyze PIP2->PIP3 conversion of final chain unit.
    Rule('PIP2_self_catalysis_1',
         PIP(bakt=None, both=None, S='PIP2', bpi3k_self=ANY, bself2=1) % PIP(bakt=None, both=None, S='PIP2', bpi3k_self=1, bself2=None) >>
         PIP(bakt=None, both=None, S='PIP2', bpi3k_self=ANY, bself2=None) + PIP(bakt=None, both=None, S='PIP3', bpi3k_self=None, bself2=None),
          Parameter('PIP2_self_catalysis_1_kf', KCP))
    
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
    catalyze_state(PI3K(bgab1=ANY), 'bpip', PIP(bakt=None, both=None, bself2=None), 'bpi3k_self', 'S', 'PIP2', 'PIP3', (KF, KR, KCP))
             
     # Setting up the binding reactions necessary for AKT to be phosphorylated and move through the pathway
    bind_table([[                       AKT(S='U', both=None),       AKT(S='P', both=None),      AKT(S='PP', both=None)],
                [PIP(S='PIP3', both=None, bpi3k_self=None),       (KF, KR),     (KF, KR),    (KF, KR)]],
                'bakt', 'bpip3')

    # AKT-PIP3 is phosphorylated by PDK1 to AKTP
    catalyze(PDK1, 'bakt', AKT(bpip3=ANY, S='U'), 'both', AKT(bpip3=ANY, S='P'), (KF, KR, KCP))

    # AKTP-PIP3 is phosphorylated by PDK1 to AKTPP
    catalyze(PDK1, 'bakt', AKT(bpip3=ANY, S='P'), 'both', AKT(bpip3=ANY, S='PP'), (KF, KR, KCP))

    # AKTP is dephosphorylated by PP2A-III back to AKT
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'both', 'S', 'P', 'U',(KF, KR, KCD))
   
    # AKTPP is dephosphorylated by PP2A-III back to AKTP
    catalyze_state(PP2A_III, 'bakt', AKT(bpip3=None), 'both', 'S', 'PP', 'P',(KF, KR, KCD))

    # PIP3 is dephosphorylated by PTEN to PIP2
    catalyze_state(PTEN, 'bpip3', PIP(bakt=None, bpi3k_self=None), 'both', 'S', 'PIP3', 'PIP2', (KF, KR, KCD))

    # PIP3 is dephosphorylated by SHP to PIP2
    catalyze_state(SHP, 'bpip3', PIP(bakt=None, bpi3k_self=None), 'both', 'S', 'PIP3', 'PIP2', (KF, KR, KCD))

def crosstalk_monomers():
    Monomer('Pase9t', ['bgab1'])
    alias_model_components()
def crosstalk_initial():
    Parameter('Pase9t_0', 0) #c521

def crosstalk_events():
    #ERK:P:P phosphorylates GAP-GRB2-GAB1:P (making it unable to bind PI3K)
    catalyze_state(ERK(st='PP'), 'b', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None, bpi3k2=None, bpi3k3=None, bpi3k4=None, bpi3k5=None, bpi3k6=None), 'bERKPP', 'S', 'P', 'PP', (KF, KR, KCP))

    #GAP-GRB2-GAB1:P:P is dephosphorylated by Pase9t
    catalyze_state(Pase9t(), 'bgab1', GAB1(bgrb2=ANY), 'bPase9t', 'S', 'PP', 'P', (KF, KR, KCD))

    #ERK:P:P phosphorylates GRB2-SOS, preventing RAS-GDP->RAS-GTP conversion
    catalyze_state(ERK(st='PP'), 'b', SOS(bgrb=ANY, bras=None), 'bERKPP', 'st', 'U', 'P', (KF, KR, KCP))

    #AKT:P:P phosphorylates RAF:P at Ser295, preventing MEK phosphorylation.
    catalyze_state(AKT(S='PP', bpip3=None), 'both', RAF(st='P'), 'b', 'ser295', 'U', 'P', (KF, KR, KCP))

    #RAS-GTP binds PI3K
    bind(RAS(bsos=None, braf=None, st='GTP'), 'bpi3k', PI3K(bgab1=None, bpip=None), 'bras', [KF, KR])

    #RAS-GTP-PI3K binds PIP2
    bind(PI3K(bgab1=None, bpip=None, bras=ANY), 'bpip', PIP(bakt=None, both=None, bself2=None, S='PIP2'), 'bpi3k_self',
         [KF, KR])

    #RAS-GTP-PI3K-PIP2 disassociates to give PIP3
    Rule('RASGTPPI3KcatPIP',
         RAS(bsos=None, braf=None, st='GTP', bpi3k=1) % PI3K(bgab1=None, bpip=2, bras=1) % PIP(bakt=None, both=None, bpi3k_self=2, S='PIP2') >>
         RAS(bsos=None, braf=None, st='GTP', bpi3k=1) % PI3K(bgab1=None, bpip=None, bras=1) + PIP(bakt=None, both=None, bpi3k_self=None, S='PIP3'),
         Parameter('RASGTPPI3K_kc', KCP))
