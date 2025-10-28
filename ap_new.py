import itertools
import subprocess
import random
import mod
from mod import LabelSettings, LabelType, LabelRelation, BondType, ruleGMLString, graphGMLString

ls = LabelSettings(LabelType.Term, LabelRelation.Unification)
lsString = LabelSettings(LabelType.String, LabelRelation.Unification)

termBondFromBondType = {
	BondType.Invalid: "__error1",
	BondType.Single: "p(0)",
	BondType.Double: "p(p(0))",
	BondType.Triple: "p(p(p(0)))",
	BondType.Aromatic: "__error2"
}

cycle3 = mod.graphDFS("[*]1{*}[*]{*}[*]{*}1")
cycle4 = mod.graphDFS("[*]1{*}[*]{*}[*]{*}[*]{*}1")
badList = [cycle3, cycle4]
for i in [-1, 1]:
	break
	g = mod.graphDFS("[a(_a1, %d)]{*}[a(_a2, %d)]" % (i, i))
	badList.append(g)
doubleDouble = mod.graphDFS("[*]{e(p(p(0)))}[*]{e(p(p(0)))}[*]")
badList.append(doubleDouble)
badMols = set()


def termFromGraph(g):
	s = "graph [\n"
	for v in g.vertices:
		s += 'node [ id %d label "a(%s, %d)" ]' % (v.id, v.atomId.symbol, v.charge)
	for e in g.edges:
		s += 'edge [ source %d target %d label "e(%s)" ]' % (e.source.id, e.target.id, termBondFromBondType[e.bondType])
	s +="]\n"
	return graphGMLString(s, name=g.name + ", term", add=False)

def decodeVertexLabel(l):
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
		# raise ValueError("Can not convert edge term '%s' to a molecule." % l)
	return bt

def graphFromTerm(g):
	s = "graph [\n"
	for v in g.vertices:
		s += 'node [ id %d label "%s" ]\n' % (v.id, decodeVertexLabel(v.stringLabel))
	for e in g.edges:
		s += 'edge [ source %d target %d label "%s" ]\n' % (e.source.id, e.target.id, decodeEdgeLabel(e.stringLabel))
	s += "]\n"
	return graphGMLString(s, name=g.name + ", string", add=False)

def _haxAtomMapEvaluateToMakeRules(vMap):
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

rules = {}  # due to a bug, the rules must be kept alive outside the DG
# def dgProject(dg, gMap, path, withEdges=True): # Why are we even doing this? Atom Maps?
def dgProject(dg, gMap, withEdges=True): # Why are we even doing this? Atom Maps?
	graphs = {}
	rule_strings = set()
	for v in dg.vertices:
		graph = gMap(v.graph)
		if graph:
			graphs[v.graph] = graph
	dgNew = mod.DG(graphDatabase=[], labelSettings=lsString)
	eMap = {}
	with dgNew.build() as b:
		for e in dg.edges:
			if not withEdges:
				d = mod.Derivation()
				d.left = [graphs[v.graph] for v in e.sources]
				d.right = [graphs[v.graph] for v in e.targets]
			else:
				rule_strings.clear()
				vms = mod.DGVertexMapper(e, upToIsomorphismGDH = True, rightLimit = 1)
				for vm in vms:
					rule_strings.add(_haxAtomMapEvaluateToMakeRules(vm))
				rs = []
				for s in rule_strings:
					if s not in rules:
						rules[s] = ruleGMLString(s, add=False)
					rs.append(rules[s])
				for r in rs: #WHY??????
					d = mod.Derivation()
					d.left = [graphs[v.graph] for v in e.sources]
					d.right = [graphs[v.graph] for v in e.targets]
					d.rule = r
				eNew = b.addDerivation(d) ## <- This is always just the last rule found
			eMap[eNew] = e
	graphMap = {}
	for old, new in graphs.items():
		graphMap[new] = old
	return dgNew, graphMap, eMap


def allPartitions(molecule):
	def split(molecule, i, component1, component2):
		if i == len(molecule):
			yield component1, component2
			return
		component1.append(molecule[i])
		yield from split(molecule, i + 1, component1, component2)
		component1.pop()
		component2.append(molecule[i])
		yield from split(molecule, i + 1, component1, component2)
		component2.pop()
	yield from split(molecule, 0, [], [])



# def computeNextIteration(prev, rules, build, dg, sizeLimit, seen, direction, digraph):
def computeNextIteration(prev, rules, build, dg, sizeLimit, seen, direction):
	result = []
	def addSeenGraphs(gs):
		t = tuple(sorted(gs, key=lambda g: g.id))
		seen.add(t)
	def hasBeenSeen(gs):
		t = tuple(sorted(gs, key=lambda g: g.id)) # This meeses with the comparison of the states.... They are not the same order when they are bidirectional.
		return t in seen
	for gs in prev:
		if hasBeenSeen(gs):
			continue
		if sizeLimit is not None and len(gs) > sizeLimit:
			continue
		for gsSub, gsOther in allPartitions(gs):
			if len(gsSub) == 0:
				continue
			if hasBeenSeen(gsSub):
				for e in getOutEdges(dg, gsSub):
					gsNew = [v.graph for v in e.targets] + gsOther
					# Should we even add it here?
					# if direction == "forward":
					# 	digraph.add_edge(tuple(sorted(gs, key=lambda g: g.id)), tuple(sorted(gsNew, key=lambda g: g.id)))
					# else:
					# 	digraph.add_edge(tuple(sorted(gsNew, key=lambda g: g.id)), tuple(sorted(gs, key=lambda g: g.id)))
					if not hasBeenSeen(gsNew):
						result.append(gsNew)
			else:
				for r in rules:
					es = build.apply(gsSub, r)
					for e in es:
						isBad = False
						targets_list = [x for x in e.targets]
						if len(targets_list) > sizeLimit:
							isBad = True						
						for v in e.targets:
							if isBad:
								break
							for pattern in badList:
								if pattern.monomorphism(v.graph, labelSettings=ls) > 0:
									isBad = True
									badMols.add(v.graph)
									break
								
						if isBad:
							continue
						result_gsSub = [v.graph for v in e.targets]
						gsNew = [v.graph for v in e.targets] + gsOther
						# if direction == "forward":
						# 	digraph.add_edge(tuple(sorted(gsSub, key=lambda g: g.id)), tuple(sorted(result_gsSub, key=lambda g: g.id)))
						# else:
						if direction == "reverse":
							inverseRule = ruleDB[r.id]
							derivation = mod.Derivation()
							derivation.left = result_gsSub
							derivation.rule = inverseRule
							derivation.right = gsSub
							build.addDerivation(derivation)
							# build.apply(gsNew, inverseRule) # Add Inverse Rule
							# digraph.add_edge(tuple(sorted(result_gsSub, key=lambda g: g.id)), tuple(sorted(gsSub, key=lambda g: g.id)))
						if not hasBeenSeen(gsNew):
							result.append(gsNew)
				addSeenGraphs(gsSub)
		addSeenGraphs(gs)
	return result



def makeDG(*, rules, sources, graphDatabase, sizeLimit=None, iterationLimit=None, withEdgeProjection=True, targets=None):
	# nGraph = networkx.DiGraph() # This is fine.
	seen = set()
	print(rules)
	dg = mod.DG(graphDatabase=graphDatabase, labelSettings=ls)
	with dg.build() as b:
		prev_sources = [sources]
		prev_targets = [targets]
		# nGraph.add_node(tuple(sorted(sources, key=lambda g: g.id))) # Not ideal way to add them.
		# nGraph.add_node(tuple(sorted(targets, key=lambda g: g.id)))
		iteration = 0
		while len(prev_sources) > 0 or len(prev_targets) > 0:
			iteration += 1
			if iterationLimit is not None and iteration > iterationLimit: # Should be changed so it exits at an overlap in seen_sources and seen_targets.
				break
			print("Iteration %d on %d source sets and %d target sets" % (iteration, len(prev_sources), len(prev_targets)))
			# next_sources = computeNextIteration(prev_sources, rules, b, dg, sizeLimit, seen, "forward", nGraph)
			next_sources = computeNextIteration(prev_sources, rules, b, dg, sizeLimit, seen, "forward")
			prev_sources = next_sources
			# next_targets = computeNextIteration(prev_targets, rules, b, dg, sizeLimit, seen, "reverse", nGraph)
			next_targets = computeNextIteration(prev_targets, rules, b, dg, sizeLimit, seen, "reverse")
			prev_targets = next_targets

	# path = networkx.dijkstra_path(nGraph, source=tuple(sorted(sources, key=lambda g: g.id)), target=tuple(sorted(targets, key=lambda g: g.id)))
	# print(path)

	print("Projecting term DG to string DG (%d vertices, %d edges)." % (dg.numVertices, dg.numEdges))
	# dgString, graphMap, edgeMap = dgProject(dg, graphFromTerm, path, withEdgeProjection)
	dgString, graphMap, edgeMap = dgProject(dg, graphFromTerm, withEdgeProjection)
	print("Projection done")
	class DGData(object):
		pass
	res = DGData()
	res.dg = dg
	res.dgString = dgString
	res.graphMap = graphMap
	res.edgeMap = edgeMap
	
	# return res, path
	return res

#############################
# Loads the partial charges #
#############################

def loadPartialCharges(dgData):
	dgString = dgData.dgString
	graphMap = dgData.graphMap
	atomValString = {}

	for vertex in dgString.vertices:
		partialCharges = computeAllCharges(vertex.graph)
		atomValString[vertex.graph] = partialCharges

	atomVal = {}
	for gString, pcString in atomValString.items():
		g = graphMap[gString]
		pc = {}
		for v in g.vertices:
			vString = gString.getVertexFromExternalId(v.id)
			pc[v] = pcString[str(vString.id)]
		atomVal[g] = pc
	class AtomValues(object):
		pass
	res = AtomValues()
	res.atomVal = atomVal
	res.atomValString = atomValString
	res.lowerBound = 200
	res.scale = 100
	return res


# def calcPathways(*, ruleData, dgData, sources, targets, maxNumSplits=None, pathGraphs):
def calcPathways(*, ruleData, dgData, sources, targets, maxNumSplits=None):
	dg = dgData.dg
	valLowerBound = ruleData.atomVals.lowerBound
	valScale = ruleData.atomVals.scale
	
	valMap = {}
	for e in dg.edges:
		vms = mod.DGVertexMapper(e, upToIsomorphismGDH = True, rightLimit = 1)
		vals = []
		# print(e)
		# for s in e.sources:
		# 	print(s.graph.getGMLString())
		for vm in vms:
			vals.append(ruleData.eval(vm.rule ,vm))
		valMap[e] = min(vals)

	flow = mod.Flow(dg, ilpSolver="CPLEX")
	for g in badMols:
		flow.addConstraint(mod.vertex(g) == 0)
	if maxNumSplits is not None:
		expr = mod.FlowLinExp()
		for e in dg.edges:
			if e.numTargets > 1:
				expr += mod.isEdgeUsed[e]
		flow.addConstraint(expr <= maxNumSplits)
	
	# for vertex in dg.vertices:
	# 	if vertex.graph not in pathGraphs:
	# 		flow.exclude(vertex)

	s = sum(a.numVertices for a in sources)
	t = sum(a.numVertices for a in targets)
	print(f"s: {s}, t: {t}")
	assert s == t
	srcsUnique = set(sources)
	tarsUnique = set(targets)
	for a in srcsUnique:
		flow.addSource(a)
		flow.addConstraint(mod.inFlow(a) >= sources.count(a))
	for a in tarsUnique:
		flow.addSink(a)
	obj = mod.FlowLinExp()
	for e in dg.edges:
		if not e.inverse.isNull():
			flow.addConstraint(mod.isBothReverseUsed(e) == 0)
		m = valMap[e] + valLowerBound
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


def printSolutions(*, ruleData, dgData, flowData):
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
		pGraph.collapseHydrogens = False
		pGraph.simpleCarbons = False
		p.withInlineGraphs = True
		vVis = lambda g, dg: s.eval(mod.vertex(graphMap[g])) != 0
		p.pushEdgeLabel(lambda e: ", ".join(r.name for r in edgeMap[e].rules))
		p.pushVertexVisible(vVis)
		p.pushEdgeVisible(lambda e: s.eval(mod.edge(edgeMap[e])) != 0)
		p.pushEdgeLabel(lambda e: "%.2f" % valMap[edgeMap[e]])
		fDG, fCoords = dgString.print(p)

		fName = "out/%d_%d_aux.tex" % (flow.id, s.id)
		with open(fName, "w") as f:
			f.write("\\resizebox{\\paperwidth}{!}{%%\n")
			f.write("\\input{%s}\n" % fDG[:-4])
			f.write("\\begin{tikzpicture}[remember picture, overlay]\n")
			for v in dgString.vertices:
				g = v.graph
				if not vVis(g, dgString): continue
				for vGraph in g.vertices:
					# print(vGraph)
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

breakSingleBond = ruleGMLString("""rule [
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
breakDoubleBond = ruleGMLString("""rule [
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

########################
# Potential New Rules: #
########################

# Moving Positve Charge:
formPositiveCharge = ruleGMLString("""rule [
	ruleID "pushPositiveCharge"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(_S1, 1)" ]
		edge [ source 0 target 1 label "e(p(0))" ]
	]
	right [
		node [ id 0 label "a(_S0, 1)" ]
		node [ id 1 label "a(_S1, 0)" ]

	]
]""", add=False)
breakPositiveCharge = formPositiveCharge.makeInverse()
formDoublePositiveCharge = ruleGMLString("""rule [
	ruleID "pushDoublePositiveCharge"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(_S1, 1)" ]
		edge [ source 0 target 1 label "e(p(p(_E)))" ]
	]
	right [
		node [ id 0 label "a(_S0, 1)" ]
		node [ id 1 label "a(_S1, 0)" ]
		edge [ source 0 target 1 label "e(p(_E))" ]
	]
]""", add=False)
breakDoublePositiveCharge = formDoublePositiveCharge.makeInverse()

# Moving Negative Charge:
formNegativeCharge = ruleGMLString("""rule [
	ruleID "formNegativeCharge"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(_S1, -1)" ]
	]
	right [
		node [ id 0 label "a(_S0, -1)" ]
		node [ id 1 label "a(_S1, 0)" ]
		edge [ source 0 target 1 label "e(p(0))" ]
	]
]""", add=False)
breakNegativeCharge = formNegativeCharge.makeInverse()
formDoubleNegativeCharge = ruleGMLString("""rule [
	ruleID "pushDoubleNegativeCharge"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(_S1, -1)" ]
		edge [ source 0 target 1 label "e(p(p(_E)))" ]
	]
	right [
		node [ id 0 label "a(_S0, -1)" ]
		node [ id 1 label "a(_S1, 0)" ]
		edge [ source 0 target 1 label "e(p(_E))" ]
	]
]""", add=False)
breakDoubleNegativeCharge = formDoubleNegativeCharge.makeInverse()

# Lone Electron Pairs used as nucleophile:
oxygenLoneElectronPairs = ruleGMLString("""rule [
	ruleID "oxygenLoneElectronPairs"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 1)" ]
		node [ id 1 label "a(O, 0)" ]
	]
	right [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(O, 1)" ]
		edge [ source 0 target 1 label "e(p(0))" ]

	]
]""", add=False)
# reverseOxygenElectronPairs = oxygenLoneElectronPairs.makeInverse()
nitrogenLoneElectronPairs = ruleGMLString("""rule [
	ruleID "nitrogenLoneElectronPairs"
	labelType "term"
	left [
		node [ id 0 label "a(_S0, 1)" ]
		node [ id 1 label "a(N, 0)" ]
	]
	right [
		node [ id 0 label "a(_S0, 0)" ]
		node [ id 1 label "a(N, 1)" ]
		edge [ source 0 target 1 label "e(p(0))" ]

	]
]""", add=False)
# reverseNitrogenElectronPairs = nitrogenLoneElectronPairs.makeInverse()

ruleDB = {
	formNegativeCharge.id: breakNegativeCharge,
	breakNegativeCharge.id: formNegativeCharge,
	formDoubleNegativeCharge.id: breakDoubleNegativeCharge,
	breakDoubleNegativeCharge.id: formDoubleNegativeCharge,
	formPositiveCharge.id: breakPositiveCharge,
	breakPositiveCharge.id: formPositiveCharge,
	formDoublePositiveCharge.id: breakDoublePositiveCharge,
	breakDoublePositiveCharge.id: formDoublePositiveCharge,
	formSingleBond.id: breakSingleBond,
	breakSingleBond.id: formSingleBond,
	formDoubleBond.id: breakDoubleBond,
	breakDoubleBond.id: formDoubleBond,
	nitrogenLoneElectronPairs.id: breakSingleBond,
	oxygenLoneElectronPairs.id: breakSingleBond
}

def chargeSeparation():
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
			# Make similar to _haxAtomMapEvaluate...
			v0l = vMap.match[v0.left].vertex
			v0r = vMap.comatch[v0.right].vertex
			v1l = vMap.match[v1.left].vertex
			v1r = vMap.comatch[v1.right].vertex
			# v0l, v0r = vMap(v0)
			# v1l, v1r = vMap(v1)
			check(v0l, v0r)
			check(v1l, v1r)
			val0 = atomVal[v0l.graph][v0l]
			val1 = atomVal[v1l.graph][v1l]
			if r in [breakSingleBond, breakDoubleBond]:
				return val1 - val0
			elif r in [formSingleBond, formDoubleBond]:
				return -(val0 - val1)
			# elif r in [formNegativeCharge, formDoubleNegativeCharge, breakNegativeCharge, breakDoubleNegativeCharge,
			# 	  	   formPositiveCharge, formDoublePositiveCharge, breakPositiveCharge, breakDoublePositiveCharge,
			# 		   oxygenLoneElectronPairs, nitrogenLoneElectronPairs]:
			# 	return 0.0
			else:
				assert False
	return RuleData([breakSingleBond, formSingleBond, breakDoubleBond, formDoubleBond,
				#   pushPositiveCharge, pushNegativeCharge, pushDoublePositiveCharge, pushDoubleNegativeCharge,
				#   formNegativeCharge, formDoubleNegativeCharge, breakNegativeCharge, breakDoubleNegativeCharge,
				#   formPositiveCharge, formDoublePositiveCharge, breakPositiveCharge, breakDoublePositiveCharge,
				#   oxygenLoneElectronPairs, nitrogenLoneElectronPairs
	])


# Haxed Functions:
def getOutEdges(dg, molecules):
  assert dg
  molecules = list(sorted(molecules))
  result = []
  for edge in dg.edges:
    gsCandidate = list(sorted(vertex.graph for vertex in edge.sources))
    if gsCandidate == molecules:
      result.append(edge)
  return result


#########################
# Haxed Partial Charge: #
#########################
def computePartialCharge():
	return random.random()

def computeAllCharges(molecule):
	charges = {}
	for vertex in molecule.vertices:
		charges[str(vertex.id)] = round(calcGasteigerCharge(vertex), 2)
	return charges

def calcGasteigerCharge(vertex):
	total = 0

	for alpha in range(1,7):
		listHigher, listLower = seperateNeighbors(vertex.incidentEdges)

		sum = 0
		for j in listHigher: # List of neighboring vertices with higher electronegativity.
			maxDiff = getIonPotential(vertex.stringLabel) / 0.75#/ ionPotential[j.stringLabel]
			electronDiff = getElectronegativity(j.stringLabel) - getElectronegativity(vertex.stringLabel)
			sum += (1/maxDiff) * electronDiff

		for k in listLower: # List of neighboring vertices with lower electronegativity.
			maxDiff = getIonPotential(k.stringLabel) / 0.75#/ ionPotential[vertex.stringLabel]
			electronDiff = getElectronegativity(k.stringLabel) - getElectronegativity(vertex.stringLabel)
			sum += (1/maxDiff) * electronDiff

		sum *= 0.5**(alpha-1)
		total += sum
	charge = 0
	if "+" in vertex.stringLabel:
		charge = 1
	if "-" in vertex.stringLabel:
		charge = -1
	total += charge
	return total

def getElectronegativity(label):
	sum = electronegativity[label.replace("1-", "").replace("1+", "")]
	return sum

def getIonPotential(label):
	return ionPotential[label.replace("1-", "").replace("1+", "")]



def seperateNeighbors(edges):
	lower = []
	higher = []
	for edge in edges:
		if getElectronegativity(edge.target.stringLabel) < getElectronegativity(edge.source.stringLabel):
			if edge.stringLabel == "=":
				lower.append(edge.target)
			lower.append(edge.target)
		elif getElectronegativity(edge.target.stringLabel) > getElectronegativity(edge.source.stringLabel):
			if edge.stringLabel == "=":
				higher.append(edge.target)
			higher.append(edge.target)
		else:
			continue
	return higher, lower

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
