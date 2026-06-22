import ap
import time
import itertools
import mod

class ExDaniel(ap.Instance):
	def __init__(self):
		daniel_start = mod.Graph.fromSMILES("BrC1C=C(C2C=CC=CC=2)N=CN=1")
		# daniel_start2 = mod.Graph.fromSMILES("[NH2-]")
		daniel_start3 = mod.Graph.fromSMILES("N")
		daniel_end = mod.Graph.fromSMILES("[NH2-]C1=CC(C2C=CC=CC=2)NC=N1")
		daniel_end2 = mod.Graph.fromSMILES("Br")
		ds = ap.termFromGraph(daniel_start)
		# ds2 = ap.termFromGraph(daniel_start2)
		ds3 = ap.termFromGraph(daniel_start3)
		de = ap.termFromGraph(daniel_end)
		de2 = ap.termFromGraph(daniel_end2)

		self.nameSet = set([ds, ds3, de, de2])
		self.sources = [ds, ds3]
		self.targets = [de, de2]
		self.name = "DanielTest"
		self.iteration_limit = 3

exAll = [ExDaniel]
