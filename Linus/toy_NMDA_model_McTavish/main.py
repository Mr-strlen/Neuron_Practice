from neuron import h, gui

# Load mod file containing McTavish et al. NMDA mechanism
h.nrn_load_dll("models/mctavish_syncbylocation/src/x86_64/libnrnmech.dylib") # for Mac

# Create a section
soma = h.Section()
soma.L = soma.diam = 20
soma.insert('pas')
soma.g_pas = 1e-5
soma.e_pas = -70

soma2 = h.Section()
soma2.L = soma2.diam = 20
soma2.insert('pas')
soma2.g_pas = 1e-5
soma2.e_pas = -70

# Insert McTavish et al. NMDA receptor mechanism
syn_nmda = h.AmpaNmda(0.5, sec=soma)
#syn_nmda.gnmda = 0.1
#syn_nmda.mg = 1
# Set NMDA receptor properties
#soma.gnmda = 0.01
#soma.mg = 1.2

# Set stimulus
stim = h.NetStim()
stim.number = 1
stim.start = 10
stim.interval = 0
stim.noise = 0

# Add an NMDA synapse to the soma
syn = h.Exp2Syn(0.5, sec=soma)
syn.tau1 = 0.1
syn.tau2 = 1
syn.e = 0

syn2 = h.Exp2Syn(0.5, sec=soma2)
syn2.tau1 = 0.1
syn2.tau2 = 1
syn2.e = 0
#syn.gmax = 0.1

stim_to_syn = h.NetCon(stim, syn_nmda)
stim_to_syn2 = h.NetCon(stim, syn2)
#stim_to_syn = h.NetCon(stim, syn_nmda)
stim_to_syn.weight[0] = 4000
stim_to_syn2.weight[0] = 0.01

# Record voltage
v_vec = h.Vector().record(soma(0.5)._ref_v)
v_vec2 = h.Vector().record(soma2(0.5)._ref_v)

# Run simulation
#h.load_file('stdrun.hoc')
h.init()
h.tstop = 300
h.dt = 0.1
h.run()

# Plot results
import matplotlib.pyplot as plt
plt.plot(v_vec)
plt.plot(v_vec2+1)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.show()
