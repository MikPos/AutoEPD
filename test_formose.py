import ap_new as ap
import time
import itertools
import mod

class ExFormoseKetoEnol(ap.Instance):
	def __init__(self):
		formaldehyde = mod.Graph.fromSMILES("C=O", name="Formaldehyde")
		glycolaldehyde = mod.Graph.fromSMILES( "OCC=O", name="Glycolaldehyde")
		post_keto = mod.Graph.fromSMILES("OC=CO", "PostKeto")
		post_aldol = mod.Graph.fromSMILES("OC(C=O)CO", "PostAldol")
		f = ap.termFromGraph(formaldehyde)
		g = ap.termFromGraph(glycolaldehyde)
		pk = ap.termFromGraph(post_keto)
		pa = ap.termFromGraph(post_aldol)

		self.nameSet = set([g,pk])
		self.sources = [g]
		self.targets = [pk]
		self.name = "Keto-Enol"
		self.iteration_limit = 4

class ExFormoseAldol(ap.Instance):
	def __init__(self):
		formaldehyde = mod.Graph.fromSMILES("C=O", name="Formaldehyde")
		glycolaldehyde = mod.Graph.fromSMILES( "OCC=O", name="Glycolaldehyde")
		post_keto = mod.Graph.fromSMILES("OC=CO", "PostKeto")
		post_aldol = mod.Graph.fromSMILES("OC(C=O)CO", "PostAldol")
		f = ap.termFromGraph(formaldehyde)
		g = ap.termFromGraph(glycolaldehyde)
		pk = ap.termFromGraph(post_keto)
		pa = ap.termFromGraph(post_aldol)

		self.nameSet = set([f,pk,pa])
		self.sources = [f,pk]
		self.targets = [pa]
		self.name = "Aldol"
		self.iteration_limit = 4

exAll = [ExFormoseKetoEnol, ExFormoseAldol]
