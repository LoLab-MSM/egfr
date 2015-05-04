# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:20:24 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from parameter_dict_A431 import parameter_dict as par

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