import ap_new as ap
import time
import itertools
import mod

class ExLongKlaus1(ap.Instance):
	def __init__(self):
		# Long and Klaus Example 1 (Fig 7B)
		start1 = mod.Graph.fromSMILES("COC1CCC(CBr)CC1")
		start2 = mod.Graph.fromSMILES("I")
		end1 = mod.Graph.fromSMILES("OC1CCC(CBr)CC1")
		end2 = mod.Graph.fromSMILES("CI")
		s1 = ap.termFromGraph(start1)
		s2 = ap.termFromGraph(start2)
		e1 = ap.termFromGraph(end1)
		e2 = ap.termFromGraph(end2)

		self.nameSet = set([s1,s2,e1,e2])
		self.sources = [s1,s2]
		self.targets = [e1,e2]
		self.name = "LongKlaus1"
		self.iteration_limit = 3

class ExLongKlaus2(ap.Instance):
	def __init__(self):
		# Long and Klaus Example 2 (Fig 7C)
		start1 = mod.Graph.fromSMILES("C=C(C)C(=O)OC")
		start2 = mod.Graph.fromSMILES("O")
		end1 = mod.Graph.fromSMILES("C=C(C)C(=O)O")
		end2 = mod.Graph.fromSMILES("CO")
		s1 = ap.termFromGraph(start1)
		s2 = ap.termFromGraph(start2)
		e1 = ap.termFromGraph(end1)
		e2 = ap.termFromGraph(end2)

		self.nameSet = set([s1,s2,e1,e2])
		self.sources = [s1,s2]
		self.targets = [e1,e2]
		self.name = "LongKlaus2"
		self.iteration_limit = 3

exAll = [ExLongKlaus1, ExLongKlaus2]
