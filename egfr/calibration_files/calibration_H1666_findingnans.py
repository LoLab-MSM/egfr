# Fits EGFR model against Chen/Sorger 2009 experimental data.

import bayessb
import pysb.integrate
import numpy
import matplotlib.pyplot as plt
import os
import itertools
import pickle
import sys

from erbb_exec import model

with open('calibration_H1666_highEGF_fittedparams_unnorm_splined_fromhighEGFval_priorvar1.5', 'rb') as handle:
    fittedparams = pickle.loads(handle.read())

for i in range(len(model.parameters)):
    if model.parameters[i].name in fittedparams:
        model.parameters[i].value = fittedparams[model.parameters[i].name]

def normalize(trajectories):
    """Rescale a matrix of model trajectories to 0-1"""
    ymin = trajectories.min(0)
    ymax = trajectories.max(0)
    return ((trajectories - ymin) / (ymax - ymin))

def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return numpy.vstack([recarray[name] for name in names]).T

def likelihood(mcmc, position):
    """Distance between model trajectories and experimental data"""
    ysim = mcmc.simulate(position, observables=True)
    ysim_array_he = extract_records(ysim, obs_names)
    ysim_slice_he = ysim_array_he[exptimepts]
   # ysim_slice_he = ysim_array_he[splined_exp_pts]
    #FIXME: exp_var is really a sdev (and so is prior_var)
    # Column 0 = AKTPP, column 1 = ErbB1P, column 2 = ERKPP
    try:
        objAKT[mcmc.iter] = numpy.sum((yslice_exp_he[:,0] - ysim_slice_he[:,0]) ** 2 / (2 * yslice_exp_var_he[:,0] ** 2))
        objErb[mcmc.iter] = numpy.sum((yslice_exp_he[:,1] - ysim_slice_he[:,1]) ** 2 / (2 * yslice_exp_var_he[:,1] ** 2))
        objERK[mcmc.iter] = numpy.sum((yslice_exp_he[:,2] - ysim_slice_he[:,2]) ** 2 / (2 * yslice_exp_var_he[:,2] ** 2))
        # mcmc.positions[mcmc.iter,:] = mcmc.test_position
        # if numpy.isnan(objAKT[mcmc.iter]) or numpy.isnan(objErb[mcmc.iter]) or numpy.isnan(objERK[mcmc.iter]) is True:
        #     #print some information about the maximum-likelihood estimate parameter set
        #     print
        #     print '%-10s %-12s %-12s %s' % ('parameter', 'actual', 'fitted', 'log10(fit/actual)')
        #     fitted_values = mcmc.cur_params()[mcmc.estimate_idx]
        #     param_dict = {}
        #     for param, new_value in zip(opts.estimate_params, fitted_values):
        #         change = numpy.log10(new_value / param.value)
        #         values = (param.name, param.value, new_value, change)
        #         print '%-10s %-12.2g %-12.2g %-+6.2f' % values
        #         param_dict[param.name] = new_value

        #         with open('calibration_H1666_highEGF_findingnans7', 'wb') as handle:
        #             pickle.dump(param_dict, handle)

        #             name = [p.name for p in opts.estimate_params]
        #             oldvalues = [p.value for p in opts.estimate_params]
        #             newvalues = mcmc.cur_params()[mcmc.estimate_idx]
        #             oldvalues_array = numpy.array(oldvalues)
        #             newvalues_array = numpy.array(newvalues)
        #             change = numpy.log10(newvalues_array / oldvalues_array)
        #             combined = numpy.column_stack((name, oldvalues, newvalues, change))
        #             numpy.savetxt('calibration_H1666_highEGF_paramchange_findingnans7.txt', combined, delimiter=' ', fmt='%s')
        #             numpy.save('calibration_alltestedparams_H1666_highEGF_findingnans7.npy', mcmc.positions)
        #             numpy.save('calibration_allpositions_H1666_highEGF_findingnans7.npy', mcmc.get_mixed_accepts(burn=opts.nsteps/10))
        #             numpy.save('calibration_fittedparams_H1666_highEGF_findingnans7.npy', zip(opts.estimate_params, mcmc.cur_params()[mcmc.estimate_idx]))
        #             numpy.save('calibration_likelihoods_H1666_highEGF_findingnans7.npy', mcmc.likelihoods)
        #             numpy.save('calibration_AKTlikelihoods_H1666_highEGF_findingnans7.npy', objAKT)
        #             numpy.save('calibration_Erblikelihoods_H1666_highEGF_findingnans7.npy', objErb)
        #             numpy.save('calibration_ERKlikelihoods_H1666_highEGF_findingnans7.npy', objERK)
        #             sys.exit()
        return objAKT[mcmc.iter] + objErb[mcmc.iter] + objERK[mcmc.iter]
    except AttributeError:
         return  numpy.sum((yslice_exp_he - ysim_slice_he) ** 2 / (2 * yslice_exp_var_he ** 2))

def prior(mcmc, position):
    """Distance to original parameter values"""
    return numpy.sum((position - prior_mean) ** 2 / ( 2 * prior_var))

def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 1 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

# data is already scaled to 0-1
splined_exp_pts = numpy.array([0,20,41,62,82,124,249,374,499,999])
data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_H1666_highEGF_unnorm_splined_nonzeroinit.npy')
ydata_exp_he = numpy.load(data_filename)
yslice_exp_he = ydata_exp_he[splined_exp_pts]
var_data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_var_H1666_highEGF_unnorm_splined_unweighted_nonzeroinit.npy')
exp_var_he = numpy.load(var_data_filename) #Standard deviation was calculated from the mean by assuming a coefficient of variation of .25; sdev's equal to 0 were set to 1 to avoid division by 0 errors
yslice_exp_var_he = exp_var_he[splined_exp_pts]

tspan = numpy.linspace(0,7200, num=144000)
exptimepts = numpy.array([0, 2999, 5999, 8999, 11999, 17999, 35999, 53999, 71999, 143999])

obs_names = ['obsAKTPP', 'obsErbB1_ErbB_P_CE', 'obsERKPP']

opts = bayessb.MCMCOpts()
opts.model = model
opts.tspan = tspan
opts.integrator = 'vode'
opts.nsteps = 50000
opts.order = 5
opts.method = 'bdf'
# opts.ixpr = True
opts.anneal_length = 12500
opts.norm_step_size = .4

scenario = 1

#Initialize arrays for recording objective function values for each variable (AKT, ErbB1, and ERK)
objAKT = numpy.zeros(opts.nsteps)
objErb = numpy.zeros(opts.nsteps)
objERK = numpy.zeros(opts.nsteps)

# A few estimation scenarios:
if scenario == 1:
    # estimate rates only (not initial conditions)
    opts.estimate_params = model.parameters_rules()
elif scenario == 2:
    # use hessian
    opts.estimate_params = model.parameters_rules()
    # Warning: hessian-guidance is expensive when fitting many parameters -- the
    # time to calculate the hessian increases with the square of the number of
    # parameters to fit!
    opts.use_hessian = True
    opts.hessian_period = opts.nsteps / 6
else:
    raise RuntimeError("unknown scenario number")

# values for prior calculation
prior_mean = [numpy.log10(p.value) for p in opts.estimate_params]
# prior_var is set to 3.0 so that (since calc is in log space) parameters can vary within 3 orders of magnitude and not be penalized.
prior_var =  4


opts.likelihood_fn = likelihood
opts.prior_fn = prior
opts.step_fn = step
opts.seed = 1
#opts.atol=numpy.load('calibration_atol_values.npy')
opts.atol=1e-3
opts.rtol=1e-3
opts.intsteps = 5000
mcmc = bayessb.MCMC(opts)

mcmc.run()

#print some information about the maximum-likelihood estimate parameter set
print
print '%-10s %-12s %-12s %s' % ('parameter', 'actual', 'fitted', 'log10(fit/actual)')
fitted_values = mcmc.cur_params()[mcmc.estimate_idx]
param_dict = {}
for param, new_value in zip(opts.estimate_params, fitted_values):
    change = numpy.log10(new_value / param.value)
    values = (param.name, param.value, new_value, change)
    print '%-10s %-12.2g %-12.2g %-+6.2f' % values
    param_dict[param.name] = new_value

with open('calibration_H1666_highEGF_fittedparams_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit', 'wb') as handle:
    pickle.dump(param_dict, handle)

name = [p.name for p in opts.estimate_params]
oldvalues = [p.value for p in opts.estimate_params]
newvalues = mcmc.cur_params()[mcmc.estimate_idx]
oldvalues_array = numpy.array(oldvalues)
newvalues_array = numpy.array(newvalues)
change = numpy.log10(newvalues_array / oldvalues_array)
combined = numpy.column_stack((name, oldvalues, newvalues, change))
numpy.savetxt('calibration_H1666_highEGF_paramchange_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.txt', combined, delimiter=' ', fmt='%s')
numpy.save('calibration_alltestedparams_H1666_highEGF_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.npy', mcmc.positions)
numpy.save('calibration_allpositions_H1666_highEGF_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.npy', mcmc.get_mixed_accepts(burn=opts.nsteps/10))
numpy.save('calibration_fittedparams_H1666_highEGF_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.npy', zip(opts.estimate_params, mcmc.cur_params()[mcmc.estimate_idx]))
numpy.save('calibration_likelihoods_H1666_highEGF_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.npy', mcmc.likelihoods)
numpy.save('calibration_AKTlikelihoods_H1666_highEGF_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.npy', objAKT)
numpy.save('calibration_Erblikelihoods_H1666_highEGF_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.npy', objErb)
numpy.save('calibration_ERKlikelihoods_H1666_highEGF_unnorm_splined_fromhighEGFval_priorvar4_somenans_nonzeroinit.npy', objERK)
