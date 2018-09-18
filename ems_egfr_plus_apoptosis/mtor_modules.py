# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:23:07 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from .parameter_dict_A431 import parameter_dict as par

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