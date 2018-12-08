import numpy as np
from scipy import constants as cst
from GeneralFunctions import general_functions as gf
from BuildSample import build_sample
from Source import source
from UpdateEquations import update_function
from dispersive_materials import Materials


class Simulation:

	def run_simulation(self, height, width, substrate_mat, nano_mat, nano_shape, source_energy, incident_angle):
		dim = 2 #just keep like this for now
		vacuum_height = height
		relative_epsilon = substrate['substrate eps_r']
		delta_x = substrate['delta_x']
		conductivity = substrate['conductivity']
		conductivity_m = substrate['conductivity_m']

		eps_space, cond_space, condm_space = build_sample.create_substrate(dim, vacuum_height, relative_epsilon, delta_x, conductivity, conductivity_m)
		# if block is not none:
		# 	space = block['existing space']
		# 	relative_epsilon = substrate['substrate eps_r']
		# 	delta_x = substrate['delta_x']
		# 	conductivity = substrate['conductivity']
		# 	conductivity_m = substrate['conductivity_m']
		# 	sample.add_block(space, dim, center, block_height, rel_eps, conduct, conduct_m, delta_x)
