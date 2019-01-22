import numpy as np
from scipy import constants as cst
import gf
from sample import Sample
from source import Source
import update_eqs as ue
from materials import Material
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt


class Simulation:

	def run_simulation(self, vacuum_height):
		#make sample, make sure we're unit independent
		dim = np.array([10.,10.])
		rel_eps = 2
		dx = 1
		cond = 1
		mcond = 1
		substrate = Sample(dim=dim, vacuum_height=vacuum_height, rel_eps=rel_eps, dx=dx, cond=cond, mcond=mcond)
		bdim = np.array([5.,5.])
		center = np.array([0.,0.])
		print('center (nm)', center)
		substrate.add_block(block_dim=bdim, center=center, rel_eps=rel_eps, cond=cond, mcond=mcond)

		# print(substrate.eps_space)
		# print(substrate.cond_space)
		# print(substrate.mcond_space)

		# plt.figure(dpi=100)
		# plt.imshow(sample.eps_space)
		# plt.axis('off')
		# plt.show()

		#make source
		#run simulation

		# if block is not none:
		# 	space = block['existing space']
		# 	relative_epsilon = substrate['substrate eps_r']
		# 	delta_x = substrate['delta_x']
		# 	conductivity = substrate['conductivity']
		# 	conductivity_m = substrate['conductivity_m']
		# 	sample.add_block(space, dim, center, block_height, rel_eps, conduct, conduct_m, delta_x)

if __name__=='__main__':
	sim = Simulation()
	sim.run_simulation(vacuum_height=5)
