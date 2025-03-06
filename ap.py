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
		raise ValueError("Can not convert edge term '%s' to a molecule." % l)
	return bt

def graphFromTerm(g):
	s = "graph [\n"
	for v in g.vertices:
		s += 'node [ id %d label "%s" ]\n' % (v.id, decodeVertexLabel(v.stringLabel))
	for e in g.edges:
		s += 'edge [ source %d target %d label "%s" ]\n' % (e.source.id, e.target.id, decodeEdgeLabel(e.stringLabel))
	s += "]\n"
	return graphGMLString(s, name=g.name + ", string", add=False)


_haxRuleStrings = set()
def _haxAtomMapEvaluateToMakeRules(r, vMap):
	left = ""
	right = ""
	for v in r.vertices:
		vl, vr = vMap(v)
		if not vl.isNull():
			left  += '\t\tnode [ id %d label "%s" ]\n' % (v.id, decodeVertexLabel(vl.stringLabel))
		if not vr.isNull():
			right += '\t\tnode [ id %d label "%s" ]\n' % (v.id, decodeVertexLabel(vr.stringLabel))
	for e in r.edges:
		vsl, vsr = vMap(e.source)
		vtl, vtr = vMap(e.target)
		if not e.left.isNull():
			el = None
			for eCand in vsl.incidentEdges:
				if eCand.target == vtl:
					el = eCand
					break
			assert el is not None
			left += '\t\tedge [ source %d target %d label "%s" ]\n' % (e.source.id, e.target.id, decodeEdgeLabel(el.stringLabel))
		if not e.right.isNull():
			er = None
			for eCand in vsr.incidentEdges:
				if eCand.target == vtr:
					er = eCand
					break
			assert er is not None
			right += '\t\tedge [ source %d target %d label "%s" ]\n' % (e.source.id, e.target.id, decodeEdgeLabel(er.stringLabel))
	s = "rule [\n\tleft [\n%s\t]\n\tright [\n%s\t]\n]\n" % (left, right)
	_haxRuleStrings.add(s)
	return 0

rules = {}  # due to a bug, the rules must be kept alive outside the DG
def dgProject(dg, gMap, withEdges=True):
	graphs = {}
	for v in dg.vertices:
		graphs[v.graph] = gMap(v.graph)
	dgNew = mod.DG(graphDatabase=[], labelSettings=lsString)
	eMap = {}
	print(graphs)
	with dgNew.build() as b:
		for e in dg.edges:
			if not withEdges:
				d = mod.Derivation()
				d.left = [graphs[v.graph] for v in e.sources]
				d.right = [graphs[v.graph] for v in e.targets]
			else:
				_haxRuleStrings.clear()
				print("e is:")
				print(e)
				mod.atomMapEvaluate(e, _haxAtomMapEvaluateToMakeRules) # Missing the AtomMapEvaluate function.
				print("done with edge evaluation")
				rs = []
				for s in _haxRuleStrings:
					if s not in rules:
						rules[s] = ruleGMLString(s, add=False)
					rs.append(rules[s])
				for r in rs:
					d = mod.Derivation()
					d.left = [graphs[v.graph] for v in e.sources]
					d.right = [graphs[v.graph] for v in e.targets]
					d.rule = r
			eNew = b.addDerivation(d)
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


def makeDG(*, rules, sources, graphDatabase, sizeLimit=None, iterationLimit=None, withEdgeProjection=True, targets=None):
	seen = set()
	print("graphDatabase: " + str(graphDatabase))
	def addSeenGraphs(gs):
		t = tuple(sorted(gs, key=lambda g: g.id))
		seen.add(t)
	def hasBeenSeen(gs): 
		t = tuple(sorted(gs, key=lambda g: g.id))
		return t in seen
	dg = mod.DG(graphDatabase=graphDatabase, labelSettings=ls)
	with dg.build() as b:
		prev = [sources]
		if targets:
			prev.append(targets)
		iteration = 0
		while len(prev) > 0:
			iteration += 1
			if iterationLimit is not None and iteration > iterationLimit:
				break
			print("Iteration %d on %d sets" % (iteration, len(prev)))
			next_ = []
			count = 0
			for gs in prev:
				def makePretty(a):
					s = a.graphDFS
					s = s.replace("{e(p(0))}", "")
					s = s.replace("{e(p(p(0)))}", "=")
					s = s.replace("[a(C, 0)]", "C")
					s = s.replace("[a(C, 1)]", "[C+]")
					s = s.replace("[a(C, -1)]", "[C-]")
					s = s.replace("[a(O, 0)]", "O")
					s = s.replace("[a(O, 1)]", "[O+]")
					s = s.replace("[a(O, -1)]", "[O-]")
					return s
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
							if not hasBeenSeen(gsNew):
								next_.append(gsNew)
					else:
						for r in rules:
							if count == 0:
								print(r)
								for v in r.right.vertices:
									print(v.stringLabel)
									for e in v.incidentEdges:
										print(e.stringLabel)
							es = b.apply(gsSub, r)
							for e in es:
								isBad = False
								for v in e.targets:
									for pattern in badList:
										if pattern.monomorphism(v.graph, labelSettings=ls) > 0:
											isBad = True
											badMols.add(v.graph)
											break
									if isBad:
										break
								if isBad:
									continue

								gsNew = [v.graph for v in e.targets] + gsOther
								print("has not been seen")
								print(gsNew)
								if count < 1:
									for g in gsNew:
										for v in g.vertices:
											print(v.stringLabel)
											for e in v.incidentEdges:
												print(e.stringLabel)
									count += 1
								if not hasBeenSeen(gsNew):
									next_.append(gsNew)
						addSeenGraphs(gsSub)
				addSeenGraphs(gs)
			prev = next_

	print("Projecting term DG to string DG (%d vertices, %d edges)." % (dg.numVertices, dg.numEdges))
	dgString, graphMap, edgeMap = dgProject(dg, graphFromTerm, withEdgeProjection)
	print("Projection done")
	class DGData(object):
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
	


def calcPathways(*, ruleData, dgData, sources, targets, maxNumSplits=None):
	dg = dgData.dg
	valLowerBound = ruleData.atomVals.lowerBound
	valScale = ruleData.atomVals.scale

	valMap = {}
	for e in dg.edges:
		vals = mod.atomMapEvaluate(e, ruleData.eval)
		valMap[e] = min(vals)

	flow = mod.Flow(dg)
	for g in badMols:
		flow.addConstraint(mod.vertex(g) == 0)
	if maxNumSplits is not None:
		expr = mod.FlowLinExp()
		for e in dg.edges:
			if e.numTargets > 1:
				expr += mod.isEdgeUsed[e]
		flow.addConstraint(expr <= maxNumSplits)

	s = sum(a.numVertices for a in sources)
	t = sum(a.numVertices for a in targets)
	assert s == t
	srcsUnique = set(sources)
	tarsUnique = set(targets)
	for a in srcsUnique:
		flow.addSource(a)
		flow.addConstraint(mod.inFlow(a) >= sources.count(a))
	for a in tarsUnique:
		print(a)
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
					print(vGraph)
					label = "\\node[inner sep=1, at=(v-%d-0-v-%d.-20), anchor=160] {\\tiny $%.2f$};\n"
					f.write(label % (v.id, vGraph.id, atomValString[g][str(vGraph.id)]))
			f.write("\\end{tikzpicture}\n")
			f.write("}%%\n")
		mod.post("compileTikz \"%s\" \"%s\" 3" % (fName[:-4], fCoords[:-4]))
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

# pushPositiveCharge = ruleGMLString("""rule [
# 	ruleID "pushPositiveCharge"
# 	labelType "term"
# 	left [
# 		node [ id 0 label "a(_S0, 0)" ]
# 		node [ id 1 label "a(_S1, 1)" ]
# 		edge [ source 0 target 1 label "e(p(0))" ]
# 	]
# 	right [
# 		node [ id 0 label "a(_S0, 1)" ]
# 		node [ id 1 label "a(_S1, 0)" ]
		
# 	]
# ]""", add=False)

# pushNegativeCharge = ruleGMLString("""rule [
# 	ruleID "pushNegativeCharge"
# 	labelType "term"
# 	left [
# 		node [ id 0 label "a(_S0, 0)" ]
# 		node [ id 1 label "a(_S1, -1)" ]
# 	]
# 	right [
# 		node [ id 0 label "a(_S0, -1)" ]
# 		node [ id 1 label "a(_S1, 0)" ]
# 		edge [ source 0 target 1 label "e(p(0))" ]
# 	]
# ]""", add=False)

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
			v0l, v0r = vMap(v0)
			v1l, v1r = vMap(v1)
			check(v0l, v0r)
			check(v1l, v1r)
			val0 = atomVal[v0l.graph][v0l]
			val1 = atomVal[v1l.graph][v1l]
			if r in [breakSingleBond, breakDoubleBond]:
				return val1 - val0
			elif r in [formSingleBond, formDoubleBond]:
				return -(val0 - val1)
			elif r in [pushPositiveCharge, pushNegativeCharge]:
				return 0.0
			else:
				assert False
	return RuleData([breakSingleBond, formSingleBond, breakDoubleBond, formDoubleBond, 
				#   pushPositiveCharge, pushNegativeCharge])
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