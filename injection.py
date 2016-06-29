from pylab import *
import nest
import nest.voltage_trace
import nest.raster_plot
from NeuroTools import signals

hh_defaults = {
	'g_L': 20.0,
	'E_L':-56.0,
	'g_Na': 7000.0,
	'E_Na':55.0,
	'g_K':1800.0 
}
## Basic Current Injection

# # Env Setup
# nest.ResetKernel()

# ## Set Model Defaults
# nest.SetDefaults('hh_cond_exp_traub',hh_defaults)
# nest.SetKernelStatus({ 'resolution': 0.01 })
# nest.SetKernelStatus({'overwrite_files':True})

# hh_neuron = nest.Create('hh_cond_exp_traub', n=1)
# nest.SetStatus(hh_neuron, {'I_e': 2000.})

# voltmeter = nest.Create('voltmeter')
# nest.SetStatus(voltmeter, {'to_file': True, 'interval': 0.1 })

# nest.Connect(voltmeter,hh_neuron)

# nest.Simulate(100)

# nest.voltage_trace.from_device(voltmeter)
# nest.voltage_trace.show()

## Current vs FR Curve

# Env Setup
nest.ResetKernel()
nest.SetDefaults('hh_cond_exp_traub',hh_defaults)
nest.SetKernelStatus({ 'resolution': 0.01 })
nest.SetKernelStatus({'overwrite_files':True})

hh_neurons = nest.Create('hh_cond_exp_traub', n=20)

for k in range(20):
	nest.SetStatus([hh_neurons[k]], {'I_e': k * 200.})

sd = nest.Create('spike_detector')
nest.SetStatus(sd,{'to_file': True})
nest.ConvergentConnect(hh_neurons,sd)

nest.Simulate(1000)

figure()
nest.raster_plot.from_device(sd, hist=True)


