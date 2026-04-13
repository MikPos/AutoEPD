import ap
import itertools
import mod

#############################
# Different Exercises Input #
#############################
# Exercises from the Levy book on Arrow Pushing.
# Within a script, instanciate the classes to run, and run them with the ".run()" method.
#
# Exercises I, J and K are skipped so far due to lack of Radical support.

class ExA(ap.Instance):
	def __init__(self):
		exercise_a_start_1 = mod.Graph.fromSMILES("N#[C-]", "a start 1")
		exercise_a_start_2 = mod.Graph.fromSMILES("CI", "a start 2")
		exercise_a_end_1 = mod.Graph.fromSMILES("N#CC", "a end 1")
		exercise_a_end_2 = mod.Graph.fromSMILES("[I-]",  "a end 2")
		a_start_1 = ap.termFromGraph(exercise_a_start_1)
		a_start_2 = ap.termFromGraph(exercise_a_start_2)
		a_end_1 = ap.termFromGraph(exercise_a_end_1)
		a_end_2 = ap.termFromGraph(exercise_a_end_2)

		self.nameSet = set([a_start_1, a_start_2, a_end_1, a_end_2])
		self.sources = [a_start_1, a_start_2]
		self.targets = [a_end_1, a_end_2]
		self.name = "Exercise A"

class ExB(ap.Instance):
	def __init__(self):
		exercise_b_start_1 = mod.Graph.fromSMILES("C[NH-]", "b start 1")
		exercise_b_start_2 = mod.Graph.fromSMILES("CC(OC)=O", "b start 2")
		exercise_b_end = mod.Graph.fromSMILES("CNC(OC)(C)[O-]", "b end")
		b_start_1 = ap.termFromGraph(exercise_b_start_1)
		b_start_2 = ap.termFromGraph(exercise_b_start_2)
		b_end = ap.termFromGraph(exercise_b_end)

		self.nameSet = set([b_start_1, b_start_2, b_end])
		self.sources = [b_start_1, b_start_2]
		self.targets = [b_end]
		self.name = "Exercise B"

class ExC(ap.Instance):
	def __init__(self):
		exercise_c_start = mod.Graph.fromSMILES("CNC(OC)(C)[O-]", "c start")
		exercise_c_end_1 = mod.Graph.fromSMILES("CC(=O)NC", "c end 1")
		exercise_c_end_2 = mod.Graph.fromSMILES("[O-]C", "c end 2")
		c_start = ap.termFromGraph(exercise_c_start)
		c_end_1 = ap.termFromGraph(exercise_c_end_1)
		c_end_2 = ap.termFromGraph(exercise_c_end_2)

		self.nameSet = set([c_start, c_end_1, c_end_2])
		self.sources = [c_start]
		self.targets = [c_end_1, c_end_2]
		self.name = "Exercise C"

class ExD(ap.Instance):
	def __init__(self):
		exercise_d_start_1 = mod.Graph.fromSMILES("CN", "d start 1")
		exercise_d_start_2 = mod.Graph.fromSMILES("CCCCl", "d start 2")
		exercise_d_end_1 = mod.Graph.fromSMILES("CCC[NH2+]C", "d end 1")
		exercise_d_end_2 = mod.Graph.fromSMILES("[Cl-]")
		d_start_1 = ap.termFromGraph(exercise_d_start_1)
		d_start_2 = ap.termFromGraph(exercise_d_start_2)
		d_end_1 = ap.termFromGraph(exercise_d_end_1)
		d_end_2 = ap.termFromGraph(exercise_d_end_2)

		self.nameSet = set([d_start_1, d_start_2, d_end_1, d_end_2])
		self.sources = [d_start_1, d_start_2]
		self.targets = [d_end_1, d_end_2]
		self.name = "Exercise D"

class ExE(ap.Instance):
	def __init__(self):
		exercise_e_start_1 = mod.Graph.fromSMILES("CC(O)[CH2-]", "e start 1")
		exercise_e_start_2 = mod.Graph.fromSMILES("O=CC", "e start 2")
		exercise_e_end = mod.Graph.fromSMILES("CC(O)CC([O-])C", "e end")
		e_start_1 = ap.termFromGraph(exercise_e_start_1)
		e_start_2 = ap.termFromGraph(exercise_e_start_2)
		e_end = ap.termFromGraph(exercise_e_end)

		self.nameSet = set([e_start_1, e_start_2, e_end])
		self.sources = [e_start_1, e_start_2]
		self.targets = [e_end]
		self.name = "Exercise E"

class ExF(ap.Instance):
	def __init__(self):
		exercise_f_start_1 = mod.Graph.fromSMILES("CC(C)(C)O", "f start 1")
		exercise_f_start_2 = mod.Graph.fromSMILES("[H+]", "f start 2")
		exercise_f_end = mod.Graph.fromSMILES("CC(C)(C)[OH2+]", "f end")
		f_start_1 = ap.termFromGraph(exercise_f_start_1)
		f_start_2 = ap.termFromGraph(exercise_f_start_2)
		f_end = ap.termFromGraph(exercise_f_end)

		self.nameSet = set([f_start_1, f_start_2, f_end])
		self.sources = [f_start_1, f_start_2]
		self.targets = [f_end]
		self.name = "Exercise F"

class ExG(ap.Instance):
	def __init__(self):
		exercise_g_start = mod.Graph.fromSMILES("CC(C)(C)[OH2+]", "g start")
		exercise_g_end_1 = mod.Graph.fromSMILES("C[C+](C)C", "g end 1")
		exercise_g_end_2 = mod.Graph.fromSMILES("O", "g end 2")
		g_start = ap.termFromGraph(exercise_g_start)
		g_end_1 = ap.termFromGraph(exercise_g_end_1)
		g_end_2 = ap.termFromGraph(exercise_g_end_2)

		self.nameSet = set([g_start, g_end_1, g_end_2])
		self.sources = [g_start]
		self.targets = [g_end_1, g_end_2]
		self.name = "Exercise G"

class ExH(ap.Instance):
	def __init__(self):
		exercise_h_start_1 = mod.Graph.fromSMILES("COC(=O)C(C(=O)OC)CCl", "h start 1")
		exercise_h_start_2 = mod.Graph.fromSMILES("CC(C)(C)[O-]", "h start 2")
		exercise_h_end_1 = mod.Graph.fromSMILES("COC(=O)[C-](C(=O)OC)CCl", "h end 2")
		exercise_h_end_2 = mod.Graph.fromSMILES("CC(C)(C)O", "h end 2")
		h_start_1 = ap.termFromGraph(exercise_h_start_1)
		h_start_2 = ap.termFromGraph(exercise_h_start_2)
		h_end_1 = ap.termFromGraph(exercise_h_end_1)
		h_end_2 = ap.termFromGraph(exercise_h_end_2)

		self.nameSet = set([h_start_1, h_start_2, h_end_1, h_end_2])
		self.sources = [h_start_1, h_start_2]
		self.targets = [h_end_1, h_end_2]
		self.name = "Exercise H"

class ExM(ap.Instance):
	def __init__(self):
		exercise_m_start = mod.Graph.fromSMILES("COC(=C)OCC=C", "m start")
		exercise_m_end = mod.Graph.fromSMILES("COC(=O)CCC=C", "n end")
		m_start = ap.termFromGraph(exercise_m_start)
		m_end = ap.termFromGraph(exercise_m_end)

		self.nameSet = set([m_start, m_end])
		self.sources = [m_start]
		self.targets = [m_end]
		self.name = "Exercise M"
		# self.size_limit = 1
		# self.iteration_limit = 3

class ExN(ap.Instance):
	def __init__(self):
		exercise_n_start = mod.Graph.fromSMILES("C1C[CH+]C(CCC=CCCC2=CCCCC2)CC1", "n start")
		exercise_n_end = mod.Graph.fromSMILES("C1CC2C(CC1)CCC3C2CC[C+]4C3CCCC4", "n end")
		n_start = ap.termFromGraph(exercise_n_start)
		n_end = ap.termFromGraph(exercise_n_end)

		self.nameSet = set([n_start, n_end])
		self.sources = [n_start]
		self.targets = [n_end]
		self.name = "Exercise N"

exAll = [ExA, ExB, ExC, ExD, ExE, ExF, ExG, ExH, ExM] #, ExN]
