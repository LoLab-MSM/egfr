# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 16:45:18 2015

@author: Erin
"""

import pymc as pm
import numpy as np
import pysb.integrate
from mpi4py import MPI
import uuid
import os
import scipy
import theano
from theano import tensor as t
import shutil
import pickle
from pymc.backends import text
from ems_egfr_plus_apoptosis.egfr_and_apoptosis_exec import model as egfr

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

u = None
if rank == 0:
    u = str(uuid.uuid4())
    for rankn in range(1,size):
        comm.send(u, dest=rankn, tag=1)

else:
    u = comm.recv(source=0, tag=1)

basetmp = '/tmp/shockle'
catalog_dir = os.path.join(basetmp, 'pythoncompiled',  u+'-'+str(rank))
intermediate_dir = os.path.join(basetmp, 'pythonintermediate',  u+'-'+str(rank))

os.makedirs(catalog_dir, mode=0o700)
os.makedirs(intermediate_dir, mode=0o700)

#monkeypatching the catalog and intermediate_dir
scipy.weave.inline_tools.function_catalog = scipy.weave.catalog.catalog(catalog_dir)
scipy.weave.catalog.intermediate_dir = lambda: intermediate_dir

model = pm.Model()

tspan = np.linspace(0, 10000, num=10000)

solver = pysb.integrate.Solver(egfr, tspan, integrator='lsoda', nsteps=10000)

parp = egfr.parameters['PARP_0'].value

@theano.compile.ops.as_op(itypes=[t.dvector],otypes=[t.dscalar, t.dscalar]) #to use gpu use type t.fvector for all inputs/outputs
def likelihood(param_vector):
    # Sub in parameter values for current location in parameter space and simulate
    for i in range(len(param_vector)):
        egfr.parameters_rules()[name_dict[i]].value = 10**param_vector[i]

    #First run is with erlotinib present: apoptosis should happen
    solver.run()
    
    #We want to maximize this fraction (i.e. apoptosis should occur under this condition)
    erl_cparp_frac = np.max(solver.yobs['cPARP'])/parp
    
    erl_error = np.log10(erl_cparp_frac)    
    
    if np.isnan(erl_error):
        erl_error = -np.inf
    
    #Second run is with no erlotinib: apoptosis should not happen
    egfr.parameters['ERL_0'].value = 0
    
    solver.run()
    
    #This time we want to maximize the fraction of pARP that remains uncleaved
    no_erl_parp_frac = np.min(solver.yobs['obsPARP'])/parp
    
    no_erl_error = np.log10(no_erl_parp_frac)
    
    if np.isnan(no_erl_error):
        no_erl_error = -np.inf    
    
    print('Erl error = ',erl_error,'No erl error = ',no_erl_error)
    
    return np.array(erl_error), np.array(no_erl_error) #to use gpu add .astype('float32') to end of first two arrays

#Setting up PyMC model
with model:
    #Create dictionary of parameter locations in vector and names for use in likelihood function    
    name_dict = {i: param.name for i, param in enumerate([param for param in egfr.parameters_rules()])}  
    
    params = pm.Normal('params', mu=[np.log10(param.value) for param in egfr.parameters_rules()], sd=np.array([2.0]*len(egfr.parameters_rules())), shape=(len(egfr.parameters_rules())))
    
    erl, no_erl = likelihood(model.params)   
    
    erl_like = pm.ArbLikelihood('erl_output', erl)
    no_erl_like = pm.ArbLikelihood('no_erl_output', no_erl)    
    
    erl = pm.Deterministic('erl', erl)
    no_erl = pm.Deterministic('no_erl', no_erl)

    nseedchains = 10*len(egfr.parameters_rules())

    step = pm.Dream_mpi(variables=[model.params], nseedchains=nseedchains, blocked=True, multitry=5, start_random=True, save_history=True, parallel=False, adapt_crossover=True)
    
    #old_trace = text.load('2015_04_29_earm_direct_mtdreamzs_normal_prior')
    trace = pm.sample(100, step, njobs=5, use_mpi=True) #pass njobs=None to start multiple chains on different cpus
    
    if rank == 0:
        text.dump('2015_06_07_egfr_qualitative_calibration', trace)    
    
        dictionary_to_pickle = {}

        for dictionary in trace:
            for var in dictionary:
                dictionary_to_pickle[var] = trace[var] 
    
        pickle.dump(dictionary_to_pickle, open('2015_06_07_egfr_qualitative_calibration.p', 'wb'))
    
        from helper_fxns import convert_param_vec_dict_to_param_dict
        from helper_fxns import merge_traces
        from helper_fxns import print_convergence_summary

        #old_traces = pickle.load(open('2015_04_29_earm_direct_mtdreamzs_normal_prior_merged_traces_80000.p'))
        #trace_list = [old_traces, dictionary_to_pickle]
        #merged_traces = merge_traces(trace_list)
    
        #pickle.dump(merged_traces, open('2015_04_30_earm_direct_mtdreamzs_normal_prior_merged_traces_95000.p', 'wb'))
    
        trace_just_params = dictionary_to_pickle
        #trace_just_params = merged_traces
        del trace_just_params['erl_output']
        del trace_just_params['no_erl_output']
        del trace_just_params['erl']
        del trace_just_params['no_erl']
        del trace_just_params['momp']
        param_vec_dict = convert_param_vec_dict_to_param_dict(trace_just_params, egfr.parameters_rules())
        print_convergence_summary(param_vec_dict)

shutil.rmtree(catalog_dir)
shutil.rmtree(intermediate_dir)  