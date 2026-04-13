import itertools
import subprocess
import random
import mod
from mod import LabelSettings, LabelType, LabelRelation, BondType, causality

from typing import Optional

ls = LabelSettings(LabelType.Term, LabelRelation.Specialisation)
lsString = LabelSettings(LabelType.String, LabelRelation.Specialisation)

termBondFromBondType = {
	BondType.Invalid: "__error1",
	BondType.Single: "p(0)",
	BondType.Double: "p(p(0))",
	BondType.Triple: "p(p(p(0)))",
	BondType.Aromatic: "__error2"
}

cycle_of_3 = mod.Graph.fromDFS("[*]1{*}[*]{*}[*]{*}1")
cycle_of_4 = mod.Graph.fromDFS("[*]1{*}[*]{*}[*]{*}[*]{*}1")
double_double_bond = mod.Graph.fromDFS("[*]{e(p(p(0)))}[*]{e(p(p(0)))}[*]")
non_chemical_patterns = [cycle_of_3, cycle_of_4, double_double_bond]
forbidden_sub_graphs = set()


def termFromGraph(g):
	"""
	Takes a mod graph in direct formatting as an argument.
	Returns the same graph with labels converted to term formatting.
	"""
	s = "graph [\n"
	for v in g.vertices:
		s += 'node [ id %d label "a(%s, %d)" ]' % (v.id, v.atomId.symbol, v.charge)
	for e in g.edges:
		s += 'edge [ source %d target %d label "e(%s)" ]' % (e.source.id, e.target.id, termBondFromBondType[e.bondType])
	s +="]\n"
	return mod.Graph.fromGMLString(s, name=g.name + ", term", add=False)

def decodeVertexLabel(l):
	"""
	Take a string representing an vertex in term formatting as an argument.
	Returns the string decoded from term to direct formatting.
	"""
	assert l.startswith("a(")
	assert l.endswith(")")
	l = l[2:-1].split(", ")
	lab = l[0]
	c = int(l[1])
	if c > 0:
		lab += str(c) + "+"
	elif c < 0:
		lab += str(abs(c)) + "-"
	return lab

def decodeEdgeLabel(l):
	"""
	Take a string representing an edge in term formatting as an argument.
	Returns the string decoded from term to direct formatting.
	"""
	assert l.startswith("e(")
	assert l.endswith(")")
	l = l[2:-1]
	if l == "0":
		assert False
	elif l == "p(0)":
		bt = "-"
	elif l == "p(p(0))":
		bt = "="
	elif l == "p(p(p(0)))":
		bt = "#"
	else:
		bt = None
	return bt

def graphFromTerm(g):
	"""
	Takes a mod graph in term formatting as an argument.
	Returns the same graph with labels converted to direct formatting.
	"""
	s = "graph [\n"
	for v in g.vertices:
		s += 'node [ id %d label "%s" ]\n' % (v.id, decodeVertexLabel(v.stringLabel))
	for e in g.edges:
		s += 'edge [ source %d target %d label "%s" ]\n' % (e.source.id, e.target.id, decodeEdgeLabel(e.stringLabel))
	s += "]\n"
	return mod.Graph.fromGMLString(s, name=g.name + ", string", add=False)

def makeRuleFromVertexMap(vMap):
	"""
	Help Text
	"""
	left = ""
	right = ""
	# Vertices
	for v in vMap.rule.vertices:
		# Left side of rule
		if not v.left:
			continue
		vG = vMap.match[v.left]
		if not vG:
			continue
		vl = vG.vertex
		left  += '\t\tnode [ id %d label "%s" ]\n' % (v.id, decodeVertexLabel(vl.stringLabel))
		# Rigth side of rule
		if not v.right:
			continue
		vH = vMap.comatch[v.right]
		if not vH:
			continue
		vr = vH.vertex
		if not vr.isNull():
			right += '\t\tnode [ id %d label "%s" ]\n' % (v.id, decodeVertexLabel(vr.stringLabel))
	# Edges
	for e in vMap.rule.edges:
		# Left side of edge
		if not e.left:
			continue
		esG = vMap.match[e.left.source]
		etG = vMap.match[e.left.target]
		if not esG or not etG:
			continue
		el = None
		for eCand in esG.incidentEdges:
			if eCand.target == etG:
				el = eCand
				break
		assert el is not None
		left += '\t\tedge [ source %d target %d label "%s" ]\n' % (e.source.id, e.target.id, decodeEdgeLabel(el.stringLabel))
		# Right side of edge
		if not e.right:
			continue
		esH = vMap.comatch[e.right.source]
		etH = vMap.comatch[e.right.target]
		if not esH or not etH:
			continue
		er = None
		for eCand in esH.incidentEdges:
			if eCand.target == etH:
				er = eCand
				break
		assert er is not None
		right += '\t\tedge [ source %d target %d label "%s" ]\n' % (e.source.id, e.target.id, decodeEdgeLabel(er.stringLabel))
	s = "rule [\n\tleft [\n%s\t]\n\tright [\n%s\t]\n]\n" % (left, right)
	return s

def dgProject(dg, gMap, withEdges=True):
	"""
	Takes a DiGraph and a Graph Map, both containing vertices and edges in term formatting as arguments.
	Optional Argument: withEdges -> Boolean. Default to True. Can be disabled to not create an Edge Map.

	Returns a new DiGraph along with a new Graph Map and Edge Map all with direct label formatting.

	The function rebuilds the DiGraph and Graph Maps with direct labels rather than the origianl term labels.
	It also constructs an Edge Map that uses the decoded labels.

	"""
	graphs = {}
	for v in dg.vertices:
		graph = gMap(v.graph)
		if graph:
			graphs[v.graph] = graph
	dgNew = mod.DG(graphDatabase=[], labelSettings=lsString)
	edgeMap = {}
	rules = {}
	with dgNew.build() as b:
		for hyperedge in dg.edges:
			if not withEdges:
				d = mod.Derivation()
				d.left = [graphs[v.graph] for v in hyperedge.sources]
				d.right = [graphs[v.graph] for v in hyperedge.targets]
			else:
				rule_strings = set()
				vms = mod.DGVertexMapper(hyperedge, upToIsomorphismGDH = True, rightLimit = 1)
				for vm in vms:
					rule_strings.add(makeRuleFromVertexMap(vm))
				rs = []
				for s in rule_strings:
					if s not in rules:
						rules[s] = mod.Rule.fromGMLString(s, add=False)
					rs.append(rules[s])
				for r in rs:
					d = mod.Derivation()
					d.left = [graphs[v.graph] for v in hyperedge.sources]
					d.right = [graphs[v.graph] for v in hyperedge.targets]
					d.rule = r
				new_edge = b.addDerivation(d)
			edgeMap[new_edge] = hyperedge
	graphMap = {}
	for old, new in graphs.items():
		graphMap[new] = old
	return dgNew, graphMap, edgeMap


def allPartitions(graph_set):
	"""
	Takes a state within a graph expansion state space as an argument.
	Returns a list of all possible pairs of subsets of the state, where each pair of subsets contains exactly all graphs within the state.
	"""
	def split(graph_set, i, subset_1, subset_2):
		if i == len(graph_set):
			yield subset_1, subset_2
			return
		subset_1.append(graph_set[i])
		yield from split(graph_set, i + 1, subset_1, subset_2)
		subset_1.pop()
		subset_2.append(graph_set[i])
		yield from split(graph_set, i + 1, subset_1, subset_2)
		subset_2.pop()
	yield from split(graph_set, 0, [], [])



def computeNextIteration(prev, rules, build, dg, sizeLimit, seen, direction):
	"""
	Help Text
	"""
	result = []
	def addSeenGraphs(graph_set):
		t = tuple(sorted(graph_set, key=lambda g: g.id))
		seen.add(t)
	def have_graphs_been_seen(graph_set):
		t = tuple(sorted(graph_set, key=lambda g: g.id))
		return t in seen
	for graph_set in prev:
		if have_graphs_been_seen(graph_set):
			continue
		if sizeLimit is not None and len(graph_set) > sizeLimit:
			continue
		for reactant_graph_subset, rest_of_graph_set in allPartitions(graph_set):
			if len(reactant_graph_subset) == 0:
				continue
			if have_graphs_been_seen(reactant_graph_subset):
				for hyperedge in getOutEdges(dg, reactant_graph_subset):
					new_graph_set = [v.graph for v in hyperedge.targets] + rest_of_graph_set
					if not have_graphs_been_seen(new_graph_set):
						result.append(new_graph_set)
			else:
				for rule in rules:
					edge_set = build.apply(reactant_graph_subset, rule)
					for hyperedge in edge_set:
						contains_non_chemical_pattern = False
						parts_of_product = [x for x in hyperedge.targets]
						if len(parts_of_product) > sizeLimit:
							contains_non_chemical_pattern = True						
						for v in hyperedge.targets:
							if contains_non_chemical_pattern:
								break
							for pattern in non_chemical_patterns:
								if pattern.monomorphism(v.graph, labelSettings=ls) > 0:
									contains_non_chemical_pattern = True
									forbidden_sub_graphs.add(v.graph)
									break
								
						if contains_non_chemical_pattern:
							continue
						product_graph_subset = [v.graph for v in hyperedge.targets]
						new_graph_set = [v.graph for v in hyperedge.targets] + rest_of_graph_set
						if direction == "reverse":
							reversed_rule = reverse_rule_map[rule.id]
							derivation = mod.Derivation()
							derivation.left = product_graph_subset
							derivation.rule = reversed_rule
							derivation.right = reactant_graph_subset
							build.addDerivation(derivation)
						if not have_graphs_been_seen(new_graph_set):
							result.append(new_graph_set)
				addSeenGraphs(reactant_graph_subset)
		addSeenGraphs(graph_set)
	return result



def makeDG(*, rules, sources, graphDatabase, sizeLimit=None, iterationLimit=None, withEdgeProjection=True, targets=None, msgName: Optional[str]):
	"""
	Help Text
	"""
	msgPrefix = None if msgName is None else f"   makeDG({msgName}):"
	seen = set()
	dg = mod.DG(graphDatabase=graphDatabase, labelSettings=ls)
	with dg.build() as b:
		previous_sources = [sources]
		previous_targets = [targets]
		iteration = 0
		while len(previous_sources) > 0 or len(previous_targets) > 0:
			iteration += 1
			if iterationLimit is not None and iteration > iterationLimit:
				break
			if msgPrefix is not None:
				print(msgPrefix, "Iteration %d on %d source sets and %d target sets" % (iteration, len(previous_sources), len(previous_targets)))
			next_sources = computeNextIteration(previous_sources, rules, b, dg, sizeLimit, seen, "forward")
			previous_sources = next_sources
			next_targets = computeNextIteration(previous_targets, rules, b, dg, sizeLimit, seen, "reverse")
			previous_targets = next_targets

	if msgPrefix is not None:
		print(msgPrefix, "Projecting term DG to string DG (%d vertices, %d edges)." % (dg.numVertices, dg.numEdges))
	dgString, graphMap, edgeMap = dgProject(dg, graphFromTerm, withEdgeProjection)
	if msgPrefix is not None:
		print(msgPrefix, "Projection done")
	class DGData:
		pass
	res = DGData()
	res.dg = dg
	res.dgString = dgString
	res.graphMap = graphMap
	res.edgeMap = edgeMap
	
	return res

#############################
# Loads the partial charges #
#############################

def loadPartialCharges(dgData):
	"""
	Help Text
	"""
	dgString = dgData.dgString
	graphMap = dgData.graphMap
	partial_charge_graph_map_string = {}

	for vertex in dgString.vertices:
		all_partial_charges = computeAllCharges(vertex.graph)
		partial_charge_graph_map_string[vertex.graph] = all_partial_charges

	partial_charge_graph_map = {}
	for graph_string, partial_charge_string in partial_charge_graph_map_string.items():
		g = graphMap[graph_string]
		partial_charges_vertex_map = {}
		for v in g.vertices:
			vString = graph_string.getVertexFromExternalId(v.id)
			partial_charges_vertex_map[v] = partial_charge_string[str(vString.id)]
		partial_charge_graph_map[g] = partial_charges_vertex_map
	class AtomValues:
		pass
	res = AtomValues()
	res.atomVal = partial_charge_graph_map
	res.atomValString = partial_charge_graph_map_string
	res.lowerBound = 200
	res.scale = 100
	return res


###
# Not yet edited:
###
def calcPathways(*, ruleData, dgData, sources, targets, useComplexObjFunction, maxNumSplits=None):
	"""
	Help Text
	"""
	dg = dgData.dg
	valLowerBound = ruleData.atomVals.lowerBound
	valScale = ruleData.atomVals.scale
	
	valMap = {}
	for e in dg.edges:
		vms = mod.DGVertexMapper(e, upToIsomorphismGDH = True, rightLimit = 1)
		vals = []
		for vm in vms:
			vals.append(ruleData.eval(vm.rule ,vm))
		valMap[e] = min(vals)

	flow = mod.hyperflow.Model(dg, ilpSolver="CPLEX")
	for g in forbidden_sub_graphs:
		if dg.findVertex(g):
			flow.addConstraint(mod.vertex[g] == 0)
	if maxNumSplits is not None:
		expr = mod.FlowLinExp()
		for e in dg.edges:
			if e.numTargets > 1:
				expr += mod.isEdgeUsed[e]
		flow.addConstraint(expr <= maxNumSplits)

	s = sum(a.numVertices for a in sources)
	t = sum(a.numVertices for a in targets)
	print(f"s: {s}, t: {t}")
	assert s == t
	srcsUnique = set(sources)
	tarsUnique = set(targets)
	for a in srcsUnique:
		flow.addSource(a)
		flow.addConstraint(mod.inFlow[a] >= sources.count(a))
	for a in tarsUnique:
		flow.addSink(a)
	obj = mod.hyperflow.LinExp()
	for e in dg.edges:
		if not e.inverse.isNull():
			flow.addConstraint(mod.isBothReverseUsed[e] == 0)
		m = valMap[e]
		if useComplexObjFunction:
			m += valLowerBound
			m *= valScale
			assert m >= 0
			m = int(m)

		obj += mod.isEdgeUsed[e] * m
	flow.objectiveFunction = obj
	class FlowData(object):
		pass
	res = FlowData()
	res.flow = flow
	res.valMap = valMap
	return res


def printSolutions(*, ruleData, dgData, flowData, prettyPrint):
	"""
	Help Text
	"""
	atomValString = ruleData.atomVals.atomValString
	dgString = dgData.dgString
	graphMap = dgData.graphMap
	edgeMap = dgData.edgeMap
	flow = flowData.flow
	valMap = flowData.valMap
	res = []
	for s in flow.solutions:
		p = mod.DGPrinter()
		p.graphvizPrefix = """
			layout = "dot"
		"""
		pGraph = p.graphPrinter
		pGraph.setMolDefault()
		pGraph.collapseHydrogens = False if prettyPrint == False else True
		pGraph.simpleCarbons = False if prettyPrint == False else True
		p.withInlineGraphs = True
		vVis = lambda v: s.eval(mod.vertex[graphMap[v.graph]]) != 0
		p.pushVertexVisible(vVis)
		p.pushEdgeVisible(lambda e: s.eval(mod.edgeFlow[edgeMap[e]]) != 0)
		if prettyPrint == False:
			p.pushEdgeLabel(lambda e: ", ".join(r.name for r in edgeMap[e].rules))
			p.pushEdgeLabel(lambda e: "%.2f" % valMap[edgeMap[e]])
		fDG, fCoords = dgString.print(p)

		fName = "out/%d_%d_aux.tex" % (flow.id, s.id)
		with open(fName, "w") as f:
			f.write("\\resizebox{\\paperwidth}{!}{%%\n")
			f.write("\\input{%s}\n" % fDG[:-4])
			f.write("\\begin{tikzpicture}[remember picture, overlay]\n")
			for v in dgString.vertices:
				g = v.graph
				if not vVis(v): continue
				for vGraph in g.vertices:
					label = "\\node[inner sep=1, at=(v-%d-0-v-%d.-20), anchor=160] {\\tiny $%.2f$};\n"
					f.write(label % (v.id, vGraph.id, atomValString[g][str(vGraph.id)]))
			f.write("\\end{tikzpicture}\n")
			f.write("}%%\n")
		mod.post.command("compileTikz \"%s\" \"%s\" 3" % (fName[:-4], fCoords[:-4]))
		res.append(fName[:-3] + "pdf")
	return res


##########################
# Electron Pushing Rules #
##########################

breakSingleBond = mod.Rule.fromGMLString("""rule [
	ruleID "breakSingleBond"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(_S1, 0)" ]
		edge [ source 0 target 1 label "e(p(0))" ]
	]
	right [
		node [ id 0 label "a(_S0, 1)" ]
		node [ id 1 label "a(_S1, -1)" ]
	]
]""", add=False)
formSingleBond = breakSingleBond.makeInverse()
breakDoubleBond = mod.Rule.fromGMLString("""rule [
	ruleID "breakDoubleBond"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(_S1, 0)" ]
		edge [ source 0 target 1 label "e(p(p(_E)))" ]
	]
	right [
		node [ id 0 label "a(_S0, 1)" ]
		node [ id 1 label "a(_S1, -1)" ]
		edge [ source 0 target 1 label "e(p(_E))" ]
	]
]""", add=False)
formDoubleBond = breakDoubleBond.makeInverse()

reverse_rule_map = {
	formSingleBond.id: breakSingleBond,
	breakSingleBond.id: formSingleBond,
	formDoubleBond.id: breakDoubleBond,
	breakDoubleBond.id: formDoubleBond
}

def chargeSeparation():
	"""
	Help Text?
	"""
	class RuleData(object):
		def __init__(self, rules):
			self.rules = rules
			self.atomVals = None
		def eval(self, r, vMap):
			assert self.atomVals is not None
			atomVal = self.atomVals.atomVal
			def check(a, b):
				assert not a.isNull()
				assert not b.isNull()
				assert a.stringLabel[:4] == b.stringLabel[:4]
			v0 = r.getVertexFromExternalId(0)
			v1 = r.getVertexFromExternalId(1)
			v0l = vMap.match[v0.left].vertex
			v0r = vMap.comatch[v0.right].vertex
			v1l = vMap.match[v1.left].vertex
			v1r = vMap.comatch[v1.right].vertex
			check(v0l, v0r)
			check(v1l, v1r)
			val0 = atomVal[v0l.graph][v0l]
			val1 = atomVal[v1l.graph][v1l]
			if r in [breakSingleBond, breakDoubleBond]:
				return val1 - val0
			elif r in [formSingleBond, formDoubleBond]:
				return -(val0 - val1)
			else:
				assert False
	return RuleData([breakSingleBond, formSingleBond, breakDoubleBond, formDoubleBond])

def getOutEdges(dg, molecules):
  """
	Help Text
	"""
  assert dg
  molecules = list(sorted(molecules))
  result = []
  for edge in dg.edges:
    graph_set_candidate = list(sorted(vertex.graph for vertex in edge.sources))
    if graph_set_candidate == molecules:
      result.append(edge)
  return result


###############################
# Partial Charge Calculation: #
###############################
def computeAllCharges(molecule):
	"""
	Takes a mod graph as an argument.
	Returns a dictionary containing the partial charges associated with each vertex in the graph.
	Dictionary format:
	- vertex id -> partial charge
	"""
	charges = {}
	for vertex in molecule.vertices:
		charges[str(vertex.id)] = round(calcGasteigerCharge(vertex), 2)
	return charges

def calcGasteigerCharge(vertex):
	"""
	Help Text
	"""
	total = 0

	for iteration in range(1,7):
		neighbors_of_higher_electronegativity, neighbors_of_lower_electronegativity = seperateNeighbors(vertex.incidentEdges)

		sum = 0
		for neighbor in neighbors_of_higher_electronegativity:
			ion_potential = getIonPotential(vertex.stringLabel) / 0.75
			difference_in_electronegativity = getElectronegativity(neighbor.stringLabel) - getElectronegativity(vertex.stringLabel)
			sum += (1/ion_potential) * difference_in_electronegativity

		for neighbor in neighbors_of_lower_electronegativity:
			ion_potential = getIonPotential(neighbor.stringLabel) / 0.75
			difference_in_electronegativity = getElectronegativity(neighbor.stringLabel) - getElectronegativity(vertex.stringLabel)
			sum += (1/ion_potential) * difference_in_electronegativity

		sum *= 0.5**(iteration-1)
		total += sum
	charge = 0
	if "+" in vertex.stringLabel:
		charge = 1
	if "-" in vertex.stringLabel:
		charge = -1
	total += charge
	return total

def getElectronegativity(label):
	"""
	Takes an atomic label as an argument and returns the associated electronnegativity.
	During the process it removes any charge indication from the label string eg.: "+" or "-"
	"""
	sum = electronegativity[label.replace("1-", "").replace("1+", "")]
	return sum

def getIonPotential(label):
	"""
	Takes an atomic label as an argument and returns the associated ion potential.
	During the process it removes any charge indication from the label string eg.: "+" or "-"
	"""
	return ionPotential[label.replace("1-", "").replace("1+", "")]

def seperateNeighbors(edges):
	"""
	Takes a range of incident edges from a mod vertex as an argument.
	Returns two lists of vertex id's, one with id's of vertices with higher electronegativity, and and one with id's of vertices with lower electronegativity.
	"""
	lower_electronegativity_neighbors = []
	higher_electronegativity_neighbors = []
	for edge in edges:
		if getElectronegativity(edge.target.stringLabel) < getElectronegativity(edge.source.stringLabel):
			if edge.stringLabel == "=":
				lower_electronegativity_neighbors.append(edge.target)
			lower_electronegativity_neighbors.append(edge.target)
		elif getElectronegativity(edge.target.stringLabel) > getElectronegativity(edge.source.stringLabel):
			if edge.stringLabel == "=":
				higher_electronegativity_neighbors.append(edge.target)
			higher_electronegativity_neighbors.append(edge.target)
		else:
			continue
	return higher_electronegativity_neighbors, lower_electronegativity_neighbors

electronegativity = {
	"H": 2.20,
	"C": 2.55,
	"O": 3.44,
	"N": 3.04,
	"P": 2.19,
	"Cl": 3.16,
	"Br": 2.96,
	"F": 3.98,
	"Si": 1.90,
	"I": 2.66,
}

ionPotential = {
	"H": 13.59,
	"C": 11.26,
	"O": 13.61,
	"N": 14.53,
	"P": 10.48,
	"Cl": 13.01,
	"Br": 11.84,
	"F": 17.42,
	"Si": 8.15,
	"I": 10.45,
}

def checkRealisable(s: mod.hyperflow.Solution) -> bool:
	"""
	Help Text
	"""
	return causality.RealisabilityQuery(s.model.dg).findDAG(s) is not None


###################
# Instance Class: #
###################

class Instance:
	def run(self, prettyPrint=False, useComplexObjectiveFunction=False) -> None:
		"""
		Runs an instance of the problem, based on SMILES strings defined in the class being run.
		
		Optional Arguments:
		- prettyPrint -> Boolean. False by default. Removes additional print statements for testing and performance checking.
		- useComplexObjectiveFunction -> Boolean. False by default. Changes the objective function used to find a flow solution within the hypergraph.
		"""
		msgPrefix = f"Instance.run({self.name}):"
		SIZE_LIMIT = getattr(self, "size_limit", 3)
		ITERATION_LIMIT = getattr(self, "iteration_limit", 3)
		print(msgPrefix, f"size_limit={SIZE_LIMIT}, iteration_limit={ITERATION_LIMIT}")

		import time

		ruleData = chargeSeparation()

		print("=" * 80)
		print(msgPrefix, "making DG")
		print("=" * 80)
		timeStart = time.perf_counter()
		dgData = makeDG(
			rules=ruleData.rules,
			sources=self.sources,
			targets=self.targets,
			graphDatabase=self.sources + self.targets,
			sizeLimit=SIZE_LIMIT,
			iterationLimit=ITERATION_LIMIT,
			msgName=self.name
		)
		# rename graphs
		for v in dgData.dgString.vertices:
			g = v.graph
			for a in self.nameSet:
				if a.isomorphism(g, labelSettings=lsString) > 0:
					g.name = a.name
					break
		timeDG = time.perf_counter()
		print("-" * 80)
		print(msgPrefix, f"time {timeDG - timeStart} DG")
		print(msgPrefix, f"|V| = {dgData.dg.numVertices}")
		print(msgPrefix, f"|E| = {dgData.dg.numEdges}")
		print("=" * 80)
		print(msgPrefix, "getting partial charges")
		ruleData.atomVals = loadPartialCharges(dgData)
		timePartChg = time.perf_counter()
		print("-" * 80)
		print(msgPrefix, f"time {timePartChg - timeDG} partial charges")
		print("=" * 80)
		print(msgPrefix, "finding flow solutions")
		flowData = calcPathways(
			ruleData=ruleData, dgData=dgData,
			sources=self.sources, targets=self.targets,
			useComplexObjFunction=useComplexObjectiveFunction
		)
		flow = flowData.flow
		flow.addEnumerationVar(mod.isEdgeUsed)

		flow.findSolutions(verbosity=1, maxNumSolutions=1)
		flow.solutions.list()

		for solution in flow.solutions:
			checkRealisable(solution)
		
		timeFlow = time.perf_counter()
		print("-" * 80)
		print(msgPrefix, f"time {timeFlow - timePartChg} flow")
		print("=" * 80)
		print(msgPrefix, "printing solutions")
		files = printSolutions(ruleData=ruleData, dgData=dgData, flowData=flowData, prettyPrint=prettyPrint)
		for f in files:
			print(msgPrefix, "files:", f)
		timePrint = time.perf_counter()
		print("-" * 80)
		print(msgPrefix, f"time {timePrint - timeFlow} printing")
		print("=" * 80)


		print(msgPrefix, f"time {timePrint - timeStart} total")
