import mod
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import ap

class ExM_incorrect(ap.Instance):
	def __init__(self):
		exercise_13_start = mod.Graph.fromSMILES("COC(=C)OCC=C")
		exercise_13_end = mod.Graph.fromSMILES("COC(=O)CCC=C", "end")
		m_start = ap.termFromGraph(exercise_13_start)
		m_end = ap.termFromGraph(exercise_13_end)
		
		self.nameSet = set([m_start, m_end])
		self.sources = [m_start]
		self.targets = [m_end]
		self.name = "Exercise M Incorrect"
		
ExM_incorrect.run(prettyPrint=True)