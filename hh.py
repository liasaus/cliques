#! /usr/bin/env python

import nest
import nest.voltage_trace
from pylab import *

nest.ResetKernel()
nest.SetKernelStatus({
		"overwrite_files": True,
		"resolution": 0.001
	})

neuron = nest.Create("hh_cond_exp_traub", 1, [{
		"E_L": -65.#,
		#"I_e": 100.
	}])
print nest.GetStatus(neuron)
#noise = nest.Create("poisson_generator", 2)

voltmeter = nest.Create("voltmeter")
# Spike Detectors
spikes = nest.Create('spike_detector', 1)

nest.ConvergentConnect(neuron, spikes)

# nest.SetStatus(noise,[{"rate": 80000.0}, {"rate": 15000.0}])

nest.SetStatus(voltmeter,[{
		"to_file": True, 
		"withtime": True,
		"interval": 0.1
	}])

#nest.ConvergentConnect(noise, neuron, [1.2, -1.], [1.0, 1.0])
nest.Connect(voltmeter, neuron)

nest.Simulate(100.0)

nest.voltage_trace.from_device(voltmeter)
plt.gca().set_ylim([-90., 0.])

nest.raster_plot.from_device(spikes)

show()
