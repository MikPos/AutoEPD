import ap_new as ap
import itertools
import mod

#############################
# Different Exercises Input #
#############################
# Only one input can be used at a time.
# Comment in or out the exercise you would like to test before running the code.
# Make sure only one is active.
#
# Exercises I, J and K are skipped so far due to lack of Radical support.

class ExA(ap.Instance):
	def __init__(self):
		exercise_1_start_1 = mod.Graph.fromSMILES("N#[C-]", "a start 1")
		exercise_1_start_2 = mod.Graph.fromSMILES("CI", "a start 2")
		exercise_1_end_1 = mod.Graph.fromSMILES("N#CC", "a end 1")
		exercise_1_end_2 = mod.Graph.fromSMILES("[I-]",  "a end 2")
		a_start_1 = ap.termFromGraph(exercise_1_start_1)
		a_start_2 = ap.termFromGraph(exercise_1_start_2)
		a_end_1 = ap.termFromGraph(exercise_1_end_1)
		a_end_2 = ap.termFromGraph(exercise_1_end_2)

		self.nameSet = set([a_start_1, a_start_2, a_end_1, a_end_2])
		self.sources = [a_start_1, a_start_2]
		self.targets = [a_end_1, a_end_2]
		self.name = "Exercise A"

class ExB(ap.Instance):
	def __init__(self):
		exercise_2_start_1 = mod.Graph.fromSMILES("C[NH-]")
		exercise_2_start_2 = mod.Graph.fromSMILES("CC(OC)=O")
		exercise_2_end = mod.Graph.fromSMILES("CNC(OC)(C)[O-]")
		b_start_1 = ap.termFromGraph(exercise_2_start_1)
		b_start_2 = ap.termFromGraph(exercise_2_start_2)
		b_end = ap.termFromGraph(exercise_2_end)

		self.nameSet = set([b_start_1, b_start_2, b_end])
		self.sources = [b_start_1, b_start_2]
		self.targets = [b_end]
		self.name = "Exercise B"

class ExC(ap.Instance):
	def __init__(self):
		exercise_3_start = mod.Graph.fromSMILES("CNC(OC)(C)[O-]")
		exercise_3_end_1 = mod.Graph.fromSMILES("CC(=O)NC")
		exercise_3_end_2 = mod.Graph.fromSMILES("[O-]C")
		c_start = ap.termFromGraph(exercise_3_start)
		c_end_1 = ap.termFromGraph(exercise_3_end_1)
		c_end_2 = ap.termFromGraph(exercise_3_end_2)

		self.nameSet = set([c_start, c_end_1, c_end_2])
		self.sources = [c_start]
		self.targets = [c_end_1, c_end_2]
		self.name = "Exercise C"

class ExD(ap.Instance):
	def __init__(self):
		exercise_4_start_1 = mod.Graph.fromSMILES("CN", "ds1")
		exercise_4_start_2 = mod.Graph.fromSMILES("CCCCl", "ds2")
		exercise_4_end_1 = mod.Graph.fromSMILES("CCC[NH2+]C", "de1")
		exercise_4_end_2 = mod.Graph.fromSMILES("[Cl-]")
		d_start_1 = ap.termFromGraph(exercise_4_start_1)
		d_start_2 = ap.termFromGraph(exercise_4_start_2)
		d_end_1 = ap.termFromGraph(exercise_4_end_1)
		d_end_2 = ap.termFromGraph(exercise_4_end_2)

		self.nameSet = set([d_start_1, d_start_2, d_end_1, d_end_2])
		self.sources = [d_start_1, d_start_2]
		self.targets = [d_end_1, d_end_2]
		self.name = "Exercise D"

class ExE(ap.Instance):
	def __init__(self):
		exercise_5_start_1 = mod.Graph.fromSMILES("CC(O)[CH2-]")
		exercise_5_start_2 = mod.Graph.fromSMILES("O=CC")
		exercise_5_end = mod.Graph.fromSMILES("CC(O)CC([O-])C")
		e_start_1 = ap.termFromGraph(exercise_5_start_1)
		e_start_2 = ap.termFromGraph(exercise_5_start_2)
		e_end = ap.termFromGraph(exercise_5_end)

		self.nameSet = set([e_start_1, e_start_2, e_end])
		self.sources = [e_start_1, e_start_2]
		self.targets = [e_end]
		self.name = "Exercise E"

class ExF(ap.Instance):
	def __init__(self):
		exercise_6_start_1 = mod.Graph.fromSMILES("CC(C)(C)O", "f start 1")
		exercise_6_start_2 = mod.Graph.fromSMILES("[H+]", "f start 2")
		exercise_6_end = mod.Graph.fromSMILES("CC(C)(C)[OH2+]", "f end")
		f_start_1 = ap.termFromGraph(exercise_6_start_1)
		f_start_2 = ap.termFromGraph(exercise_6_start_2)
		f_end = ap.termFromGraph(exercise_6_end)

		self.nameSet = set([f_start_1, f_start_2, f_end])
		self.sources = [f_start_1, f_start_2]
		self.targets = [f_end]
		self.name = "Exercise F"

class ExG(ap.Instance):
	def __init__(self):
		exercise_7_start = mod.Graph.fromSMILES("CC(C)(C)[OH2+]")
		exercise_7_end_1 = mod.Graph.fromSMILES("C[C+](C)C")
		exercise_7_end_2 = mod.Graph.fromSMILES("O")
		g_start = ap.termFromGraph(exercise_7_start)
		g_end_1 = ap.termFromGraph(exercise_7_end_1)
		g_end_2 = ap.termFromGraph(exercise_7_end_2)

		self.nameSet = set([g_start, g_end_1, g_end_2])
		self.sources = [g_start]
		self.targets = [g_end_1, g_end_2]
		self.name = "Exercise G"

class ExH(ap.Instance):
	def __init__(self):
		exercise_8_start_1 = mod.Graph.fromSMILES("COC(=O)C(C(=O)OC)CCl")
		exercise_8_start_2 = mod.Graph.fromSMILES("CC(C)(C)[O-]")
		exercise_8_end_1 = mod.Graph.fromSMILES("COC(=O)[C-](C(=O)OC)CCl")
		exercise_8_end_2 = mod.Graph.fromSMILES("CC(C)(C)O")
		h_start_1 = ap.termFromGraph(exercise_8_start_1)
		h_start_2 = ap.termFromGraph(exercise_8_start_2)
		h_end_1 = ap.termFromGraph(exercise_8_end_1)
		h_end_2 = ap.termFromGraph(exercise_8_end_2)

		self.nameSet = set([h_start_1, h_start_2, h_end_1, h_end_2])
		self.sources = [h_start_1, h_start_2]
		self.targets = [h_end_1, h_end_2]
		self.name = "Exercise H"

class ExM(ap.Instance):
	def __init__(self):
		exercise_13_start = mod.Graph.fromSMILES("COC(=C)OCC=C")
		exercise_13_end = mod.Graph.fromSMILES("COC(=O)CCC=C", "end")
		m_start = ap.termFromGraph(exercise_13_start)
		m_end = ap.termFromGraph(exercise_13_end)

		self.nameSet = set([m_start, m_end])
		self.sources = [m_start]
		self.targets = [m_end]
		self.name = "Exercise M"

class ExN(ap.Instance):
	def __init__(self):
		exercise_14_start = mod.Graph.fromSMILES("C1C[CH+]C(CCC=CCCC2=CCCCC2)CC1")
		exercise_14_end = mod.Graph.fromSMILES("C1CC2C(CC1)CCC3C2CC[C+]4C3CCCC4")
		n_start = ap.termFromGraph(exercise_14_start)
		n_end = ap.termFromGraph(exercise_14_end)

		self.nameSet = set([n_start, n_end])
		self.sources = [n_start]
		self.targets = [n_end]
		self.name = "Exercise N"
