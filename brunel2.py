import nest
import nest.raster_plot
import nest.voltage_trace
import pylab
nest.ResetKernel()

order = 2000

g = 5.0
eta = 2.0
delay = 1.0
tau_m = 20.0
V_resting = -70
V_reset = -60
V_th = -50
N_E = 4 * order
N_I = order
N_neurons = N_E + N_I
C_E = N_E / 10
C_I = N_I / 10
J_E = 0.1
J_I = -g * J_E
nu_ex = eta * ((V_th - V_resting) / (J_E * C_E * tau_m))
p_rate = 1000.0 * nu_ex * C_E
nest.SetKernelStatus({'print_time': True})

## Create Network Nodes & Synapses

nest.SetDefaults('iaf_psc_delta', {
    'C_m': 1.0,
    'tau_m': tau_m,
    't_ref': 2.0,
    'E_L': V_resting + 0.0,
    'V_th': V_th + 0.0,
    'V_reset': V_reset + 0.0
    })
nodes = nest.Create('iaf_psc_delta', N_neurons)
nodes_E = nodes[:N_E]
nodes_I = nodes[N_E:]
nest.CopyModel('static_synapse_hom_wd', 'excitatory', 
               {'weight': J_E, 'delay': delay})
nest.CopyModel('static_synapse_hom_wd', 'inhibitory', 
              {'weight': J_I, 'delay': delay})
noise = nest.Create('poisson_generator', 1, 
                    {'rate': p_rate, 'origin': 0., 'start': 0., 'stop': 200. })

## Wire It Up!

nest.RandomConvergentConnect(nodes_E, nodes, C_E, model='excitatory')
nest.RandomConvergentConnect(nodes_I, nodes, C_I, model='inhibitory')
nest.DivergentConnect(noise, nodes, model='excitatory')

## Recorders

# Spike Detectors
spikes = nest.Create('spike_detector', 2, 
                    [{'label': 'brunel-py-ex'},
                     {'label': 'brunel-py-in'}])
spikes_E = spikes[:1]
spikes_I = spikes[1:]
N_rec = 50
nest.ConvergentConnect(nodes_E[:N_rec], spikes_E)
nest.ConvergentConnect(nodes_I[:N_rec], spikes_I)

# Voltmeter
voltmeter = nest.Create('voltmeter', 1, {'withgid': True})
nest.Connect(voltmeter, nodes_E[:1])

## Run Simulation
simtime = 200
nest.Simulate(simtime)

## Post-Processing

print 'Poisson Rate:        %.2f s/s' % (p_rate / 1000.0)

events = nest.GetStatus(spikes, 'n_events')
print events
rate_ex = (events[0] / (N_rec + 0.0)) / (simtime / 1000.0)
print 'Excitatory rate:     %.2f s/s' % rate_ex
rate_in = (events[1] / (N_rec + 0.0)) / (simtime / 1000.0)
print 'Inhibitory rate:     %.2f s/s' % rate_in

pylab.figure()
nest.raster_plot.from_device(spikes_E, hist=True)

pylab.figure()
nest.raster_plot.from_device(spikes_I, hist=True)

pylab.figure()
nest.voltage_trace.from_device(voltmeter)
pylab.show()
