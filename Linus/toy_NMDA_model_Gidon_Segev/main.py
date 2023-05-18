# Here's an example code that includes both AMPA and NMDA synapses using the exp2syn mechanism for AMPA and the exp2synnmda mechanism from Gidon Segev 2012 for NMDA

import os
os.chdir("model")

import time

time_start = time.time()  # 记录开始时间

import numpy as np
import matplotlib.pyplot as plt
from neuron import h, gui

# Load mod files
#h.nrn_load_dll("nrnmech.dll")
#h.nrn_load_dll('C:\\path\\to\\mod\\files\\mod_func.dll')
h.load_file("win_gidon_segev/set_t_table.hoc")

# Create the cell
soma = h.Section(name='soma')
soma.L = soma.diam = 20
soma.insert('hh')

# Create dendrite
dend = h.Section(name='dend')
dend.L = 1000
dend.diam = 1.0
dend.insert('pas')

# Connect the sections
dend.connect(soma)

# Add synapses
syn1 = h.ExpSyn(dend(0.5))
syn2 = h.NMDA(dend(0.5))
syn1.e = 0  # AMPA reversal potential
syn2.e = 0  # NMDA reversal potential
syn1.tau = 0.5
syn2.tau1 = 0.5
syn2.tau2 = 50

# Add stimulus
stim = h.NetStim()
stim.interval = 50
stim.number = 1
stim.start = 100

# Connect stimulus to synapse
#nc1 = h.NetCon(stim, syn1)
#nc1.delay = 1
#nc1.weight[0] = 0.01

nc2 = h.NetCon(stim, syn2)
nc2.delay = 1
nc2.weight[0] = 0.01

# Record variables
soma_v = h.Vector().record(soma(0.5)._ref_v)
dend_v = h.Vector().record(dend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)

# Run the simulation
h.tstop = 1000
h.run()

time_end = time.time()

print("Time: " + str(time_end-time_start))
# Plot the results
plt.plot(t, soma_v, label='soma')
plt.plot(t, dend_v, label='dend')
plt.legend()
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.show()


