# Fits EGFR model against Chen/Sorger 2009 experimental data.

import bayessb
import pysb.integrate
import numpy
import matplotlib.pyplot as plt
import os
import itertools

from erbb_exec import model as hem 

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
    ysim_he = mcmc.simulate(position, observables=True)
    ysim_he_array = extract_records(ysim_he, obs_names)
    ysim_he_norm = normalize(ysim_he_array)
    mcmc.solver.run(mcmc.cur_params(position), y0_lEGF)
    ysim_le = mcmc.solver.yobs
    ysim_le_array = extract_records(ysim_le, obs_names)
    ysim_le_norm = normalize(ysim_le_array)
    mcmc.solver.run(mcmc.cur_params(position), y0_hHRG)
    ysim_hh = mcmc.solver.yobs
    ysim_hh_array = extract_records(ysim_hh, obs_names)
    ysim_hh_norm = normalize(ysim_hh_array)
    mcmc.solver.run(mcmc.cur_params(position), y0_lHRG)
    ysim_lh = mcmc.solver.yobs
    ysim_lh_array = extract_records(ysim_lh, obs_names)
    ysim_lh_norm = normalize(ysim_lh_array)
    #FIXME: exp_var is really a sdev (and so is prior_var)
    return numpy.sum((ydatahEGF_norm - ysim_he_norm) ** 2 / (2 * exp_var_hEGF ** 2)) + numpy.sum((ydatalEGF_norm - ysim_le_norm) ** 2 / (2 * exp_var_lEGF ** 2)) \
        + numpy.sum((ydatahHRG_norm - ysim_hh_norm) ** 2 / (2 * exp_var_hHRG ** 2)) + numpy.sum((ydatalHRG_norm - ysim_lh_norm) ** 2 / (2 * exp_var_lHRG ** 2))

def prior(mcmc, position):
    """Distance to original parameter values"""
    return numpy.sum((position - prior_mean) ** 2 / ( 2 * prior_var))

def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

# data is already scaled to 0-1
data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_A431_highEGF.npy')
ydatahEGF_norm = numpy.load(data_filename)
data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_A431_lowEGF.npy')
ydatalEGF_norm = numpy.load(data_filename)
data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_A431_highHRG.npy')
ydatahHRG_norm = numpy.load(data_filename)
data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_A431_lowHRG.npy')
ydatalHRG_norm = numpy.load(data_filename)
var_data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_var_A431_highEGF.npy')
exp_var_hEGF = numpy.load(var_data_filename) #Standard deviation was calculated from the mean by assuming a coefficient of variation of .25; sdev's equal to 0 were set to 1 to avoid division by 0 errors.
var_data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_var_A431_lowEGF.npy')
exp_var_lEGF = numpy.load(var_data_filename) #Standard deviation was calculated from the mean by assuming a coefficient of variation of .25; sdev's equal to 0 were set to 1 to avoid division by 0 errors.
var_data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_var_A431_highHRG.npy')
exp_var_hHRG = numpy.load(var_data_filename) #Standard deviation was calculated from the mean by assuming a coefficient of variation of .25; sdev's equal to 0 were set to 1 to avoid division by 0 errors.
var_data_filename = os.path.join(os.path.dirname(__file__), 'experimental_data_var_A431_lowHRG.npy')
exp_var_lHRG = numpy.load(var_data_filename) #Standard deviation was calculated from the mean by assuming a coefficient of variation of .25; sdev's equal to 0 were set to 1 to avoid division by 0 errors.
# Convergance criteria: posterior should be less than 8.36 (the sum of all experimental variances).

# Initial values for different conditions:
data_filename = os.path.join(os.path.dirname(__file__), 'y0_lEGF.npy')
y0_lEGF = numpy.load(data_filename)
data_filename = os.path.join(os.path.dirname(__file__), 'y0_hHRG.npy')
y0_hHRG = numpy.load(data_filename)
os.path.join(os.path.dirname(__file__), 'y0_lHRG.npy')
y0_lHRG = numpy.load(data_filename)

tspan = numpy.array([0., 150., 300., 450., 600., 900., 1800., 2700., 3600., 7200.]) #10 unevenly spaced time points

obs_names = ['obsAKTPP', 'obsErbB1_ErbB_P_CE', 'obsERKPP']

opts = bayessb.MCMCOpts()
opts.model = hem
opts.tspan = tspan
opts.integrator = 'vode'
opts.nsteps = 50000

scenario = 1

# A few estimation scenarios:
if scenario == 1:
    # estimate rates only (not initial conditions)
    opts.estimate_params = hem.parameters_rules()
elif scenario == 2:
    # use hessian
    opts.estimate_params = hem.parameters_rules()
    # Warning: hessian-guidance is expensive when fitting many parameters -- the
    # time to calculate the hessian increases with the square of the number of
    # parameters to fit!
    opts.use_hessian = True
    opts.hessian_period = opts.nsteps / 6
else:
    raise RuntimeError("unknown scenario number")

# values for prior calculation
prior_mean = [numpy.log10(p.value) for p in opts.estimate_params]
# prior_var is set to 6.0 so that (since calc is in log space) parameters can vary within 6 orders of magnitude and not be penalized.
prior_var =  6.0


opts.likelihood_fn = likelihood
opts.prior_fn = prior
opts.step_fn = step
opts.seed = 1
opts.atol=1e-6
opts.rtol=1e-3
opts.intsteps = 5000
opts.with_jacobian = True
mcmc = bayessb.MCMC(opts)

mcmc.run()

numpy.save('calibration_allpositions_A431.npy', mcmc.get_mixed_accepts(burn=opts.nsteps/10))
numpy.save('calibration_fittedparams_A431.npy', zip(opts.estimate_params, mcmc.cur_params()[mcmc.estimate_idx]))
#print some information about the maximum-likelihood estimate parameter set
print
print '%-10s %-12s %-12s %s' % ('parameter', 'actual', 'fitted', 'log10(fit/actual)')
fitted_values = mcmc.cur_params()[mcmc.estimate_idx]
for param, new_value in zip(opts.estimate_params, fitted_values):
    change = numpy.log10(new_value / param.value)
    values = (param.name, param.value, new_value, change)
    print '%-10s %-12.2g %-12.2g %-+6.2f' % values

# plot data and simulated trajectories before and after the fit
colors = ('r', 'g', 'b')
patterns = ('k--', 'k', 'k:x')
# generate a legend with the deconvolved colors and styles
for style in colors + patterns:
    plt.plot(numpy.nan, numpy.nan, style)
plt.legend(('pAKT', 'pErbB1', 'pERK', 'initial', 'final', 'data'), loc='lower right')
# simulate initial and final trajectories and plot those along with the data
yinitial = pysb.integrate.odesolve(model, tspan)
yinitial_array = extract_records(yinitial, obs_names)
yfinal_array = extract_records(mcmc.simulate(observables=True), obs_names)
initial_lines = plt.plot(tspan, normalize(yinitial_array))
data_lines = plt.plot(tspan, ydata_norm)
final_lines = plt.plot(tspan, normalize(yfinal_array))
for il, dl, fl, c in zip(initial_lines, data_lines, final_lines, colors):
    il.set_color(c)
    dl.set_color(c)
    fl.set_color(c)
    il.set_linestyle('--')
    dl.set_linestyle(':')
    dl.set_marker('x')
plt.show()

