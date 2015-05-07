# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:18:08 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from parameter_dict_A431 import parameter_dict as par

def mapk_monomers(simplified_raf=True, raf_dimers=True, mek_dimers=True, braf=False, craf=False, ksr=False):
    Monomer('SOS', ['bgrb', 'bras', 'bERKPP', 'st'], {'st':['U', 'P']})
    Monomer('RAS', ['bsos', 'braf', 'bpi3k', 'st'], {'st':['GDP', 'GTP']})
    if simplified_raf == True:
        if raf_dimers == True:
            Monomer('RAF', ['b', 'dim', 'st', 'ser259'], {'st':['U', 'P'], 'ser259':['U', 'P']})
        else:
            Monomer('RAF', ['b', 'st', 'ser259'], {'st':['U', 'P'], 'ser259':['U', 'P']})
    else:
        if raf_dimers == True:
            if braf == True:
                Monomer('BRAF', ['b', 'dim', 'st'], {'st':['U', 'P']})
            if craf == True:
                Monomer('CRAF', ['b', 'dim', 'st', 'ser259', 'ser338'], {'st':['U', 'P'], 'ser259':['U', 'P'], 'ser338':['U', 'P']})
            if ksr == True:
                Monomer('KSR', ['b', 'dim'])
        else:
            if braf == True:
                Monomer('BRAF', ['b', 'st'], {'st':['U', 'P']})
            if craf == True:
                Monomer('CRAF', ['b', 'st', 'ser259', 'ser338'], {'st':['U', 'P'], 'ser259':['U', 'P'], 'ser338':['U', 'P']})
            if ksr == True:
                Monomer('KSR', ['b'])
    Monomer('PP1', ['b'])
    Monomer('PP2', ['b'])
    Monomer('PP3', ['b'])
    if mek_dimers == True:
        Monomer('MEK', ['b', 'dim', 'st'], {'st':['U', 'P', 'PP']})
    else:
        Monomer('MEK', ['b', 'st'], {'st':['U', 'P', 'PP']})    
    Monomer('ERK', ['b', 'st', 'loc'], {'st':['U', 'P', 'PP'], 'loc':['C', 'N']})

def mapk_initial(simplified_raf=True, raf_dimers=True, mek_dimers=True, braf=False, craf=False, ksr=False):
    # Initial values declared in parameter dictionary for given cell type.
    alias_model_components()

    Initial(RAS(bsos=None, braf=None, bpi3k=None, st='GDP'), RAS_0)
    if simplified_raf == True:
        if raf_dimers == True:
            Initial(RAF(b=None, dim=None, st='U', ser259='U'), RAF_0)
        else:
            Initial(RAF(b=None, st='U', ser259='U'), RAF_0)
    else:
        if raf_dimers == True:
            if braf == True:
                Initial(BRAF(b=None, dim=None, st='U'), BRAF_0)
            if craf == True:
                Initial(CRAF(b=None, dim=None, st='U', ser259='U', ser338='U'), CRAF_0)
            if ksr == True:
                Initial(KSR(b=None, dim=None), KSR_0)
        
        else:
            if braf == True:
                Initial(BRAF(b=None, st='U'), BRAF_0)
            if craf == True:
                Initial(CRAF(b=None, st='U', ser259='U', ser338='U'), CRAF_0)
            if ksr == True:
                Initial(KSR(b=None), KSR_0)
    
    if mek_dimers == True:    
        Initial(MEK(b=None, dim=None, st='U'), MEK_0)
    
    else:
        Initial(MEK(b=None, st='U'), MEK_0)        
        
    Initial(ERK(b=None, st='U', loc='C'), ERK_0)
    Initial(PP1(b=None), PP1_0)
    Initial(PP2(b=None), PP2_0)
    Initial(PP3(b=None), PP3_0)
    Initial(GRB2(b=None, bsos=1, bgap=None, bgab1=None, bcbl=None) % SOS(bgrb=1, bras=None, bERKPP=None, st='U'), GRB2_SOS_0)

def mapk_events(simplified_raf=True, raf_dimers=True, mek_dimers=True, braf=False, craf=False, ksr=False):

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

    if simplified_raf == True:
    
        # Deactivation of RAF:P -> RAF by PP1
        catalyze(PP1(), 'b', RAF(st='P', ser259='U'), 'b', RAF(st='U', ser259='U'),
                 (par['RAFP_PP1']))
                 
        if raf_dimers == True:
            #Raf can dimerize.  
            #As per Freeman Mol Cell 2013, 1) Raf dimerization depends on Ras interaction, and 2) the catalytic function of WT Raf depends on its dimerization.
            # To match these results the following assumptions have been made: 
            # Only dimeric Raf is catalytically active.
            # Only Raf that has been bound  by Ras can dimerize.

            # RAS-GTP binds RAF
            bind(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', RAF(dim=None), 'b', par['RASGTP_bind_RAF'])          
        
            bind(RAF(b=ANY), 'dim', RAF(b=None), 'dim', par['Raf_Raf_dimerization'])
            
            #RAF is activated by the other RAF kinase in a hetero/homo dimer (Hu 2013 Cell).  In the simplified model we have only homodimers and no positive feedback from MEK.
            Rule('RAF_RAF_transphosphorylation',
                 RAF(b=ANY, dim=ANY) % RAF(b=None, dim=ANY, st='U') >> RAF(b=ANY, dim=ANY) % RAF(b=None, dim=ANY, st='P'),
                 par['RASGTP_RAF_cat'])    
            
            # Activation of MEK -> MEK:P by activated RAF
            catalyze(RAF(st='P', ser259='U', dim=ANY), 'b', MEK(st='U'), 'b', MEK(st='P'), (par['RAFP_MEK']))
            
            # Activation of MEK:P -> MEK:P:P by activated RAF
            catalyze(RAF(st='P', ser259='U', dim=ANY), 'b', MEK(st='P'), 'b', MEK(st='PP'), (par['RAFP_MEKP']))
        
        else:
            # Activation of RAF -> RAF:P by RAS-GTP
            catalyze_state(RAS(bpi3k=None, st='GTP'), 'braf', RAF(), 'b', 'st', 'U', 'P', par['RASGTP_bind_RAF']+[par['RASGTP_RAF_cat']])
            
            # Activation of MEK -> MEK:P by activated RAF
            catalyze(RAF(st='P', ser259='U'), 'b', MEK(st='U'), 'b', MEK(st='P'), (par['RAFP_MEK']))
            
            # Activation of MEK:P -> MEK:P:P by activated RAF
            catalyze(RAF(st='P', ser259='U'), 'b', MEK(st='P'), 'b', MEK(st='PP'), (par['RAFP_MEKP']))
    
    else:
        if braf == True:
            
            #Deactivation of BRAF:P -> BRAF by PP1
            catalyze(PP1(), 'b', BRAF(st='P'), 'b', BRAF(st='U'), (par['RAFP_PP1']))
        
            if raf_dimers == True:
                # RAS-GTP binds BRAF (undimerized or dimerized -- separate to prevent formation of complexes with two RAS molecules)
                bind(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', BRAF(dim=None), 'b', par['RASGTP_bind_RAF'])

                bind_complex(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', BRAF(dim=1, b=None) % BRAF(dim=1, b=None), 'b', par['RASGTP_bind_RAF'], m2=BRAF(dim=1, b=None))                
                
                #BRAF exhibits significant constitutive homodimerization prior to mitogen stimulation (Freeman Mol Cell 2013)
                bind(BRAF(), 'dim', BRAF(b=None), 'dim', par['Raf_Raf_dimerization'])
                
                #BRAF is activated by the other RAF kinase in a hetero/homo dimer (Hu 2013 Cell).  This activation depends on phosphorylation in the N-terminus, which is constitutive for BRAF.
                Rule('BRAF_BRAF_transphosphorylation',
                    BRAF(b=ANY, dim=ANY) % BRAF(b=None, dim=ANY, st='U') >> BRAF(b=ANY, dim=ANY) % BRAF(b=None, dim=ANY, st='P'),
                    par['RASGTP_RAF_cat'])                
                
                # Activation of MEK -> MEK:P by activated BRAF
                catalyze(BRAF(st='P', dim=ANY), 'b', MEK(st='U'), 'b', MEK(st='P'), (par['RAFP_MEK']))
                
                # Activation of MEK:P -> MEK:P:P by activated RAF
                catalyze(BRAF(st='P', dim=ANY), 'b', MEK(st='P'), 'b', MEK(st='PP'), (par['RAFP_MEKP']))
                
                if craf == True:
                    bind(BRAF(b=ANY), 'dim', CRAF(b=None), 'dim', par['Raf_Raf_dimerization'])
                    
                    Rule('BRAF_CRAF_transphosphorylation',
                    BRAF(b=ANY, dim=ANY) % CRAF(b=None, dim=ANY, st='U') >> BRAF(b=ANY, dim=ANY) % CRAF(b=None, dim=ANY, st='P'),
                    par['RASGTP_RAF_cat'])    
                    
                    bind(BRAF(b=None), 'dim', CRAF(b=ANY), 'dim', par['Raf_Raf_dimerization'])
                    
                    Rule('CRAF_BRAF_transphosphorylation',
                    CRAF(b=ANY, dim=ANY, ser338='P') % BRAF(b=None, dim=ANY, st='U') >> CRAF(b=ANY, dim=ANY, ser338='P') % BRAF(b=None, dim=ANY, st='P'),
                    par['RASGTP_RAF_cat']) 
                    
                if ksr == True:
                    bind(BRAF(b=ANY), 'dim', KSR(), 'dim', par['Raf_Raf_dimerization'])
            
            else:
                
                #Instead of BRAF phosphorylation of MEK depending on dimerization and transactivation, its activation is approximated by phosphorylation by RAS.
                catalyze_state(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', BRAF(), 'b', 'st', 'U', 'P', par['RASGTP_bind_RAF']+par['RASGTP_RAF_cat'])                
                
                # Activation of MEK -> MEK:P by activated BRAF
                catalyze(BRAF(st='P'), 'b', MEK(st='U'), 'b', MEK(st='P'), (par['RAFP_MEK']))
                
                # Activation of MEK:P -> MEK:P:P by activated BRAF
                catalyze(BRAF(st='P'), 'b', MEK(st='P'), 'b', MEK(st='PP'), (par['RAFP_MEKP']))
        
        if craf == True:
            
            #Deactivation of CRAF:P -> CRAF by PP1
            catalyze(PP1(), 'b', CRAF(st='P', ser259='U'), 'b', CRAF(st='U', ser259='U'), (par['RAFP_PP1']))
            
            if raf_dimers == True:
                
                # RAS-GTP binds CRAF
                bind(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', CRAF(dim=None), 'b', par['RASGTP_bind_RAF'])      
            
                bind(CRAF(b=ANY), 'dim', CRAF(b=None), 'dim', par['Raf_Raf_dimerization'])
                
                #CRAF is activated by the other RAF kinase in a hetero/homo dimer (Hu 2013 Cell).  This activation depends on phosphorylation in the N-terminus, which is mediated by MEK at Ser338 for CRAF.
                Rule('CRAF_CRAF_transphosphorylation',
                    CRAF(b=ANY, dim=ANY, ser338='P') % CRAF(b=None, dim=ANY, st='U') >> CRAF(b=ANY, dim=ANY, ser338='P') % CRAF(b=None, dim=ANY, st='P'),
                    par['RASGTP_RAF_cat'])
                
                # Activation of MEK -> MEK:P by activated CRAF
                catalyze(CRAF(st='P', dim=ANY, ser259='U'), 'b', MEK(st='U'), 'b', MEK(st='P'), (par['RAFP_MEK']))
                
                # Activation of MEK:P -> MEK:P:P by activated CRAF
                catalyze(CRAF(st='P', dim=ANY, ser259='U'), 'b', MEK(st='P'), 'b', MEK(st='PP'), (par['RAFP_MEKP']))
                
                #Phosphorylation of CRAF at S338 by MEK
                
                catalyze_state(MEK(st='PP'), 'b', CRAF(ser259='U'), 'b', 'ser338', 'U', 'P', (par['MEK_phos_CRAF']))
                
                if ksr == True:
                   bind(CRAF(st='P'), 'dim', KSR(), 'dim', par['Raf_Raf_dimerization']) 
            
            else:
                #Instead of CRAF phosphorylation of MEK depending on dimerization and transactivation, its activation is approximated by phosphorylation by RAS.                 
                catalyze_state(RAS(bsos=None, bpi3k=None, st='GTP'), 'braf', CRAF(ser259='U'), 'b', 'st', 'U', 'P', par['RASGTP_bind_RAF']+par['RASGTP_RAF_cat'])
                           
                # Activation of MEK -> MEK:P by activated CRAF
                catalyze(CRAF(st='P', ser259='U'), 'b', MEK(st='U'), 'b', MEK(st='P'), (par['RAFP_MEK']))
                
                # Activation of MEK:P -> MEK:P:P by activated CRAF
                catalyze(CRAF(st='P', ser259='U'), 'b', MEK(st='P'), 'b', MEK(st='PP'), (par['RAFP_MEKP']))
        
    # Deactivation of MEK:P -> MEK by PP2
    catalyze(PP2(), 'b', MEK(st='P'), 'b', MEK(st='U'), (par['MEKP_PP2']))

    # Deactivation of MEK:P:P -> MEK:P by PP2
    catalyze(PP2(), 'b', MEK(st='PP'), 'b', MEK(st='P'), (par['MEKPP_PP2']))
    
    # Activation of ERK -> ERK:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='U', loc='C'), 'b', ERK(st='P', loc='C'), (par['MEKPP_ERK']))

    # Deactivation of ERK:P -> ERK by PP3
    catalyze(PP3(), 'b', ERK(st='P', loc='C'), 'b', ERK(st='U', loc='C'), (par['ERKP_PP3']))

    # Activation of ERK:P -> ERK:P:P by activated MEK:P:P
    catalyze(MEK(st='PP'), 'b', ERK(st='P', loc='C'), 'b', ERK(st='PP', loc='C'), (par['MEKPP_ERKP']))

    # Deactivation of ERK:P:P -> ERK:P by PP3
    catalyze(PP3(), 'b', ERK(st='PP', loc='C'), 'b', ERK(st='P', loc='C'), (par['ERKPP_PP3']))

    # Degradation of PP3
    degrade(PP3(b=None), par['PP3_deg'])

def mapk_negative_feedback():
    """Negative feedback events within the MAPK pathway."""
    
    #ERK:P:P phosphorylates GRB2-SOS, preventing RAS-GDP->RAS-GTP conversion
    catalyze_state(ERK(st='PP', loc='C'), 'b', SOS(bras=None), 'bERKPP', 'st', 'U', 'P', (par['ERKPP_phos_SOS']))