from pylab import *
import nest
import nest.raster_plot
import nest.voltage_trace
import pylab

def get_spike_data(detector):
	data = nest.GetStatus(detector,'events')
	senders = np.array([])
	times = np.array([])
	for subdevice in data:
		senders = np.append(senders,subdevice['senders'])
		times = np.append(times,subdevice['times'])
	return np.asarray(np.transpose(np.asmatrix([senders,times])))

def extract_events(data, time=None, sel=None):
		"""
		rewrite of nest.raster_plot.extract_events()
		"""

		if time:
			t_max = time[-1]
			if len(time) > 1:
					t_min = time[0]
			else:
					t_min = 0

			t_matches = np.logical_and((data[:,1] >= t_min), (data[:,1] < t_max))
			data = data[t_matches]

		if sel is not None:
			gid_match_acc = zeros(np.shape(data[:,0]),dtype=bool)

			for gid in sel:
				gid_matches = (data[:,0] == gid)
				gid_match_acc = logical_or(gid_match_acc, gid_matches)

			data = data[gid_match_acc, :]
		return data

def clique_activation(data, clique, times):
	clique_data = extract_events(data, time=times, sel=clique)
	if np.size(clique_data) is 0:
		return 0
	unique_gids = np.unique(clique_data[:,0])
	activation_ratio = np.size(unique_gids) / (np.size(clique) + 0.0)
	return activation_ratio

# Nest Setup
nest.ResetKernel()
nest.SetKernelStatus({
	'print_time': True,
	'overwrite_files':True
})

# Simulation Paramters
ORDER = 2000
N_RECORDING = 100
DURATION = 100

# Poisson Parameters
POISSON_START = 60.
POISSON_DURATION = 5.
POISSON_RATE_ON = 100.
POISSON_RATE_OFF = 10.
PROBABILITY_ON = 0.3

# Network Parameters
N_E = 4 * ORDER
N_I = ORDER
N_neurons = N_E + N_I
C_E = N_E / 10
C_I = N_I / 10

# Clique Parameters
CLIQUE_SIZE = 30
CLIQUE_DENSITY = 10 # when 1, each exc neuron is in ~1 clique
N_CLIQUES = int(CLIQUE_DENSITY * (N_E / CLIQUE_SIZE))
CLIQUE_SYNAPSE_WEIGHT = 0.
CLIQUE_DURATION = 25
TRANSIENT_DURATION = 3

# Synapse Weights
g = 0
J_E = 0.5
synapses = {
	'inh': 					-g * J_E,
	'exc': 					J_E,
	'exc_E_to_E':		CLIQUE_SYNAPSE_WEIGHT * J_E,
	'exc_IN_to_I':	0.15,
	'exc_IN_to_E':	J_E,
	'exc_E_to_I':		J_E
}

# Compute Poisson Firing Rates

# Legacy Code:
# nu_ex = POISSON_DYNAMIC_SCALAR * ((V_th - V_resting) / (J_E * C_E * tau_m))
# poisson_rate = 1000.0 * nu_ex * C_E

input_off_rate = 1000.0 * POISSON_RATE_OFF
input_on_rate = 1000.0 * POISSON_RATE_ON

## Create Network Nodes & Synapses
nest.SetDefaults("izhikevich", {
		"d": 8. # big d causes cell to reach regular spiking faster
	})
nodes = nest.Create("izhikevich", N_neurons)
nodes_E = nodes[:N_E]
nodes_I = nodes[N_E:]

N_E_ON = int(np.around(np.size(nodes_E) * PROBABILITY_ON))
nodes_shuffled = list(np.random.permutation(nodes_E))
nodes_E_ON = nodes_shuffled[0:N_E_ON]
nodes_E_OFF = nodes_shuffled[N_E_ON:]
# @TODO turn the above lines into a simple function

# Create Custom Synapses
for key in synapses:
	nest.CopyModel('static_synapse_hom_wd', key, {
		'weight': synapses[key],
		'delay': delay
	})

# Poisson Inputs (off- and on- states)
IN_off = nest.Create('poisson_generator', 1, {
	'rate': input_off_rate 
})
IN_on_pre = nest.Create('poisson_generator', 1, {
	'rate': input_off_rate, 
	'start': 0., 
	'stop': POISSON_START 
})
IN_on = nest.Create('poisson_generator', 1, {
	'rate': input_on_rate, 
	'start': POISSON_START, 
	'stop': POISSON_START + POISSON_DURATION  
})
IN_on_post = nest.Create('poisson_generator', 1, {
	'rate': input_off_rate, 
	'start': POISSON_START + POISSON_DURATION, 
	'stop': DURATION + 0.
})

## Wire It Up!

nest.RandomConvergentConnect(nodes_E, nodes_I, C_E, model='exc_E_to_I')

nest.RandomConvergentConnect(nodes_I, nodes_E, C_I, model='inh')

nest.DivergentConnect(IN_off, nodes_E_OFF, model='exc_IN_to_E')

nest.DivergentConnect(IN_on_pre, nodes_I, model='exc_IN_to_I')
nest.DivergentConnect(IN_on, nodes_I, model='exc_IN_to_I')
nest.DivergentConnect(IN_on_post, nodes_I, model='exc_IN_to_I')

nest.DivergentConnect(IN_on_pre, nodes_E_ON, model='exc_IN_to_E')
nest.DivergentConnect(IN_on, nodes_E_ON, model='exc_IN_to_E')
nest.DivergentConnect(IN_on_post, nodes_E_ON, model='exc_IN_to_E')

# Connect Cliques
cliques = []
for i in range(N_CLIQUES):
	cliques.append(np.random.permutation(nodes_E)[0:CLIQUE_SIZE].tolist())
	nest.ConvergentConnect(cliques[i],cliques[i],model='exc_E_to_E')

# a la Brunel
# nest.RandomConvergentConnect(nodes_E, nodes, C_E, model='excitatory')
# nest.RandomConvergentConnect(nodes_I, nodes, C_I, model='inhibitory')
# nest.DivergentConnect(noise, nodes, model='excitatory')

## Recorders

# Spike Detectors
spikes = nest.Create('spike_detector', 6, 
										[{'label': 'brunel-py-ex', 'to_file': True},
										 {'label': 'brunel-py-in', 'to_file': True},
										 {'label': 'clique-one', 'to_file':True},
										 {'label': 'clique-two', 'to_file':True},
										 {'label': 'clique-three', 'to_file':True},
										 {'label': 'clique-four', 'to_file':True}])
spikes_E = spikes[:1]
spikes_I = spikes[1:2]
spikes_CLIQUE1 = spikes[2:3]
spikes_CLIQUE2 = spikes[3:4]
spikes_CLIQUE3 = spikes[4:5]
spikes_CLIQUE4 = spikes[5:6]
nest.ConvergentConnect(nodes_E[:N_RECORDING], spikes_E)
nest.ConvergentConnect(nodes_I[:N_RECORDING], spikes_I)
nest.ConvergentConnect(cliques[1], spikes_CLIQUE1)
nest.ConvergentConnect(cliques[2], spikes_CLIQUE2)
nest.ConvergentConnect(cliques[3], spikes_CLIQUE3)
nest.ConvergentConnect(cliques[4], spikes_CLIQUE4)

spikes_complete = nest.Create('spike_detector', 1, [{
		'label': 'all-neurons',
		'to_file': True,
		'start': POISSON_START,
		'stop':POISSON_START + CLIQUE_DURATION
	}])
nest.ConvergentConnect(nodes, spikes_complete)

# Voltmeters
N_VOLTMETERS_ON = 3
N_VOLTMETERS_OFF = 3
N_VOLTMETERS = N_VOLTMETERS_ON + N_VOLTMETERS_OFF

voltmeters = nest.Create('voltmeter', N_VOLTMETERS, {
		'withgid': True,
		'start': 55.,
		'stop': 85.,
		'interval': 0.2
	})

clique1_ON = np.intersect1d(cliques[1], nodes_E_ON)
clique1_OFF = np.intersect1d(cliques[1], nodes_E_OFF)
voltmeter_gids_ON = np.random.choice(clique1_ON, N_VOLTMETERS_ON)
voltmeter_gids_OFF = np.random.choice(clique1_OFF, N_VOLTMETERS_OFF)
voltmeter_gids = np.concatenate((voltmeter_gids_ON, voltmeter_gids_OFF))
nest.Connect(voltmeters, voltmeter_gids)

## Run Simulation
nest.Simulate(DURATION)

## Post-Processing
events = nest.GetStatus(spikes, 'n_events')
rate_ex = (events[0] / (N_RECORDING + 0.0)) / (DURATION / 1000.0)
rate_in = (events[1] / (N_RECORDING + 0.0)) / (DURATION / 1000.0)

# Measure Clique Activation
spike_data = get_spike_data(spikes_complete)

activations_transient = []
activations_recurrent = []

transient_times = [POISSON_START, POISSON_START + TRANSIENT_DURATION]
recurrent_times = [POISSON_START, POISSON_START + CLIQUE_DURATION]
for clique in cliques:
	activations_transient.append(clique_activation(spike_data, clique, transient_times))
	activations_recurrent.append(clique_activation(spike_data, clique, recurrent_times))
np.savetxt('activations_transient.dat',activations_transient)
np.savetxt('activations_recurrent.dat',activations_recurrent)


# Print Output
print 'N_excitatory cells:   %d' % N_E
print 'N_inhibitory cells:   %d' % N_I
print 'Clique count          %d' % N_CLIQUES
print 'Clique size           %d' % CLIQUE_SIZE
print 'Poisson OFF rate:     %3.1f s/s' % (input_off_rate / 1000.0)
print 'Poisson ON rate:      %3.1f s/s' % (input_on_rate / 1000.0)
print 'Excitatory rate:      %3.1f s/s' % rate_ex
print 'Inhibitory rate:      %3.1f s/s' % rate_in



nest.raster_plot.from_device(spikes_E, title='Excitatory Cells')

#nest.raster_plot.from_device(spikes_I, title='Inhibitory Cells')

nest.raster_plot.from_device(spikes_CLIQUE1, title='Clique #1')

nest.raster_plot.from_device(spikes_CLIQUE2, title='Clique #2')

nest.raster_plot.from_device(spikes_CLIQUE3, title='Clique #3')

nest.raster_plot.from_device(spikes_CLIQUE4, title='Clique #4')


figure(figsize=(6,12))
for i in range(N_VOLTMETERS):
	plt.subplot(N_VOLTMETERS, 1, i+1)
	nest.voltage_trace.from_device([voltmeters[i]],
			title='Neuron in Clique #1, gid=%d' % voltmeter_gids[i])
	plt.gca().label_outer()
	plt.gca().set_ylabel("mv")
	plt.gca().set_xlabel("")

figure()
plt.scatter(activations_transient, activations_recurrent)
plt.xlim(0,1)
plt.ylim(0,1)
plt.gca().set_ylabel("p recurrent spike")
plt.gca().set_xlabel("p transient spike")

show()
