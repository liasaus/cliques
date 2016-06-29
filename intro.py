import nest
import nest.voltage_trace
import pylab
neuron = nest.Create('iaf_neuron')
sine = nest.Create('ac_generator', 1, {
					'amplitude': 100.0,
					'frequency': 2.0})
noise = nest.Create('poisson_generator', 2, 
					[{ 'rate': 70000.0},
					{'rate':20000.0}])
voltmeter = nest.Create('voltmeter',1,{'withgid':True})
nest.Connect(sine,neuron)
nest.Connect(voltmeter,neuron)
nest.ConvergentConnect(noise,neuron,[1.0, -1.0],1.0)
nest.Simulate(1000.0)
nest.voltage_trace.from_device(voltmeter)
