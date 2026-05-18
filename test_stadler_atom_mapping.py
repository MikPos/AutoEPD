import math

import ap
import mod


TOLERANCE = 1e-12


class AtomValues:
	def __init__(self, atomVal):
		self.atomVal = atomVal

def find_expected(vertex_map, atom_values):
    expected = 0.0
    id_maps = {}
    for domain_vertex in vertex_map.map.domain.vertices:
     codomain_vertex = vertex_map.map[domain_vertex]
     id_maps[str(domain_vertex.vertex.id) + "" + domain_vertex.vertex.stringLabel] = str(codomain_vertex.vertex.id) + "" + codomain_vertex.vertex.stringLabel
     expected += abs(
			atom_values[domain_vertex.vertex.graph][domain_vertex.vertex]
			- atom_values[codomain_vertex.vertex.graph][codomain_vertex.vertex]
		)
     
    return expected,id_maps

def get_atom_values(graph_charges):
    atom_values = {}
    for graph, charges in graph_charges.items():
     atom_values[graph] = {
			vertex: charges[str(vertex.id)]
			for vertex in graph.vertices
		}
     
    return atom_values


def test_stadler_break_single_bond_on_larger_molecule():
	reactant = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C" ]
		node [ id 1 label "O" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]
		node [ id 8 label "H" ]
		node [ id 9 label "H" ]
		node [ id 10 label "H" ]
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 0 target 4 label "-" ]
		edge [ source 0 target 5 label "-" ]
		edge [ source 0 target 6 label "-" ]
		edge [ source 2 target 7 label "-" ]
		edge [ source 2 target 8 label "-" ]
		edge [ source 3 target 9 label "-" ]
		edge [ source 3 target 10 label "-" ]
	]""", name="break single reactant", add=False)
	product_left = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C+" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
		node [ id 6 label "H" ]
		edge [ source 0 target 4 label "-" ]
		edge [ source 0 target 5 label "-" ]
		edge [ source 0 target 6 label "-" ]
	]""", name="break single product left", add=False)
	product_right = mod.Graph.fromGMLString("""graph [
		node [ id 1 label "O-" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		node [ id 7 label "H" ]
		node [ id 8 label "H" ]
		node [ id 9 label "H" ]
		node [ id 10 label "H" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 2 target 7 label "-" ]
		edge [ source 2 target 8 label "-" ]
		edge [ source 3 target 9 label "-" ]
		edge [ source 3 target 10 label "-" ]
	]""", name="break single product right", add=False)
	rule = mod.Rule.fromGMLString("""rule [
		ruleID "string breakSingleBond"
		left [
			node [ id 0 label "C" ]
			node [ id 1 label "O" ]
			edge [ source 0 target 1 label "-" ]
		]
		right [
			node [ id 0 label "C+" ]
			node [ id 1 label "O-" ]
		]
	]""", add=False)

	dg = mod.DG(graphDatabase=[reactant, product_left, product_right])
	with dg.build() as build:
		derivation = mod.Derivation()
		derivation.left = [reactant]
		derivation.right = [product_left, product_right]
		derivation.rule = rule
		build.addDerivation(derivation)

	edge = next(iter(dg.edges))
	vertex_maps = list(mod.DGVertexMapper(edge, upToIsomorphismGDH=True, rightLimit=1))
	assert len(vertex_maps) == 1
	vertex_map = vertex_maps[0]

	graph_charges = {
		reactant: ap.computeAllCharges(reactant),
		product_left: ap.computeAllCharges(product_left),
		product_right: ap.computeAllCharges(product_right),
	}
	
	atom_values = get_atom_values(graph_charges)
	expected, id_maps = find_expected(vertex_map, atom_values)
	print(id_maps)
	print(graph_charges)
	rule_data = ap.chargeSeparation()
	rule_data.atomVals = AtomValues(atom_values)
	actual = rule_data.eval(vertex_map.rule, vertex_map, "stadler")

	assert math.isclose(actual, expected, abs_tol=TOLERANCE), (actual, expected)


def test_stadler_form_single_bond_between_larger_fragments():
	reactant_left = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C+" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
		node [ id 6 label "H" ]
		edge [ source 0 target 4 label "-" ]
		edge [ source 0 target 5 label "-" ]
		edge [ source 0 target 6 label "-" ]
	]""", name="form single reactant left", add=False)
	reactant_right = mod.Graph.fromGMLString("""graph [
		node [ id 1 label "O-" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		node [ id 7 label "H" ]
		node [ id 8 label "H" ]
		node [ id 9 label "H" ]
		node [ id 10 label "H" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 2 target 7 label "-" ]
		edge [ source 2 target 8 label "-" ]
		edge [ source 3 target 9 label "-" ]
		edge [ source 3 target 10 label "-" ]
	]""", name="form single reactant right", add=False)
	product = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C" ]
		node [ id 1 label "O" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]
		node [ id 8 label "H" ]
		node [ id 9 label "H" ]
		node [ id 10 label "H" ]
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 0 target 4 label "-" ]
		edge [ source 0 target 5 label "-" ]
		edge [ source 0 target 6 label "-" ]
		edge [ source 2 target 7 label "-" ]
		edge [ source 2 target 8 label "-" ]
		edge [ source 3 target 9 label "-" ]
		edge [ source 3 target 10 label "-" ]
	]""", name="form single product", add=False)
	rule = mod.Rule.fromGMLString("""rule [
		ruleID "string formSingleBond"
		left [
			node [ id 0 label "C+" ]
			node [ id 1 label "O-" ]
		]
		right [
			node [ id 0 label "C" ]
			node [ id 1 label "O" ]
			edge [ source 0 target 1 label "-" ]
		]
	]""", add=False)

	dg = mod.DG(graphDatabase=[reactant_left, reactant_right, product])
	with dg.build() as build:
		derivation = mod.Derivation()
		derivation.left = [reactant_left, reactant_right]
		derivation.right = [product]
		derivation.rule = rule
		build.addDerivation(derivation)

	edge = next(iter(dg.edges))
	vertex_maps = list(mod.DGVertexMapper(edge, upToIsomorphismGDH=True, rightLimit=1))
	assert len(vertex_maps) == 1
	vertex_map = vertex_maps[0]

	graph_charges = {
		reactant_left: ap.computeAllCharges(reactant_left),
		reactant_right: ap.computeAllCharges(reactant_right),
		product: ap.computeAllCharges(product),
	}
	atom_values = get_atom_values(graph_charges)
	expected, id_maps = find_expected(vertex_map, atom_values)
	print(id_maps)
	print(graph_charges)
	rule_data = ap.chargeSeparation()
	rule_data.atomVals = AtomValues(atom_values)
	actual = rule_data.eval(vertex_map.rule, vertex_map, "stadler")

	assert math.isclose(actual, expected, abs_tol=TOLERANCE), (actual, expected)


def test_stadler_break_double_bond_in_exm_like_chain():
	reactant = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C" ]
		node [ id 1 label "O" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "C" ]
		node [ id 5 label "C" ]
		node [ id 6 label "C" ]
		node [ id 7 label "O" ]
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "=" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
		edge [ source 6 target 0 label "-" ]
		edge [ source 6 target 7 label "=" ]
		node [ id 10 label "H" ]
		node [ id 11 label "H" ]
		node [ id 12 label "H" ]
		node [ id 13 label "H" ]
		node [ id 14 label "H" ]
		node [ id 15 label "H" ]
		node [ id 16 label "H" ]
		node [ id 17 label "H" ]
		edge [ source 0 target 10 label "-" ]
		edge [ source 0 target 11 label "-" ]
		edge [ source 2 target 12 label "-" ]
		edge [ source 3 target 13 label "-" ]
		edge [ source 4 target 14 label "-" ]
		edge [ source 4 target 15 label "-" ]
		edge [ source 5 target 16 label "-" ]
		edge [ source 5 target 17 label "-" ]
	]""", name="ExM-like break double reactant", add=False)
	product = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C" ]
		node [ id 1 label "O" ]
		node [ id 2 label "C+" ]
		node [ id 3 label "C-" ]
		node [ id 4 label "C" ]
		node [ id 5 label "C" ]
		node [ id 6 label "C" ]
		node [ id 7 label "O" ]
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
		edge [ source 6 target 0 label "-" ]
		edge [ source 6 target 7 label "=" ]
		node [ id 10 label "H" ]
		node [ id 11 label "H" ]
		node [ id 12 label "H" ]
		node [ id 13 label "H" ]
		node [ id 14 label "H" ]
		node [ id 15 label "H" ]
		node [ id 16 label "H" ]
		node [ id 17 label "H" ]
		edge [ source 0 target 10 label "-" ]
		edge [ source 0 target 11 label "-" ]
		edge [ source 2 target 12 label "-" ]
		edge [ source 3 target 13 label "-" ]
		edge [ source 4 target 14 label "-" ]
		edge [ source 4 target 15 label "-" ]
		edge [ source 5 target 16 label "-" ]
		edge [ source 5 target 17 label "-" ]
	]""", name="ExM-like break double product", add=False)
	rule = mod.Rule.fromGMLString("""rule [
		ruleID "string breakDoubleBond"
		left [
			node [ id 2 label "C" ]
			node [ id 3 label "C" ]
			edge [ source 2 target 3 label "=" ]
		]
		right [
			node [ id 2 label "C+" ]
			node [ id 3 label "C-" ]
			edge [ source 2 target 3 label "-" ]
		]
	]""", add=False)

	dg = mod.DG(graphDatabase=[reactant, product])
	with dg.build() as build:
		derivation = mod.Derivation()
		derivation.left = [reactant]
		derivation.right = [product]
		derivation.rule = rule
		build.addDerivation(derivation)

	edge = next(iter(dg.edges))
	vertex_maps = list(mod.DGVertexMapper(edge, upToIsomorphismGDH=True, rightLimit=1))
	assert len(vertex_maps) == 1
	vertex_map = vertex_maps[0]

	graph_charges = {
		reactant: ap.computeAllCharges(reactant),
		product: ap.computeAllCharges(product),
	}
	atom_values = get_atom_values(graph_charges)
	expected, id_maps = find_expected(vertex_map, atom_values)
	print(id_maps)
	print(graph_charges)
	rule_data = ap.chargeSeparation()
	rule_data.atomVals = AtomValues(atom_values)
	actual = rule_data.eval(vertex_map.rule, vertex_map, "stadler")

	assert math.isclose(actual, expected, abs_tol=TOLERANCE), (actual, expected)




def test_stadler_form_double_bond_in_exm_like_chain():
	reactant = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C" ]
		node [ id 1 label "O" ]
		node [ id 2 label "C+" ]
		node [ id 3 label "C-" ]
		node [ id 4 label "C" ]
		node [ id 5 label "C" ]
		node [ id 6 label "C" ]
		node [ id 7 label "O" ]
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
		edge [ source 6 target 0 label "-" ]
		edge [ source 6 target 7 label "=" ]
		node [ id 10 label "H" ]
		node [ id 11 label "H" ]
		node [ id 12 label "H" ]
		node [ id 13 label "H" ]
		node [ id 14 label "H" ]
		node [ id 15 label "H" ]
		node [ id 16 label "H" ]
		node [ id 17 label "H" ]
		edge [ source 0 target 10 label "-" ]
		edge [ source 0 target 11 label "-" ]
		edge [ source 2 target 12 label "-" ]
		edge [ source 3 target 13 label "-" ]
		edge [ source 4 target 14 label "-" ]
		edge [ source 4 target 15 label "-" ]
		edge [ source 5 target 16 label "-" ]
		edge [ source 5 target 17 label "-" ]
	]""", name="ExM-like form double reactant", add=False)
	product = mod.Graph.fromGMLString("""graph [
		node [ id 0 label "C" ]
		node [ id 1 label "O" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "C" ]
		node [ id 5 label "C" ]
		node [ id 6 label "C" ]
		node [ id 7 label "O" ]
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "=" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
		edge [ source 6 target 0 label "-" ]
		edge [ source 6 target 7 label "=" ]
		node [ id 10 label "H" ]
		node [ id 11 label "H" ]
		node [ id 12 label "H" ]
		node [ id 13 label "H" ]
		node [ id 14 label "H" ]
		node [ id 15 label "H" ]
		node [ id 16 label "H" ]
		node [ id 17 label "H" ]
		edge [ source 0 target 10 label "-" ]
		edge [ source 0 target 11 label "-" ]
		edge [ source 2 target 12 label "-" ]
		edge [ source 3 target 13 label "-" ]
		edge [ source 4 target 14 label "-" ]
		edge [ source 4 target 15 label "-" ]
		edge [ source 5 target 16 label "-" ]
		edge [ source 5 target 17 label "-" ]
	]""", name="ExM-like form double product", add=False)
	rule = mod.Rule.fromGMLString("""rule [
		ruleID "string formDoubleBond"
		left [
			node [ id 2 label "C+" ]
			node [ id 3 label "C-" ]
			edge [ source 2 target 3 label "-" ]
		]
		right [
			node [ id 2 label "C" ]
			node [ id 3 label "C" ]
			edge [ source 2 target 3 label "=" ]
		]
	]""", add=False)

	dg = mod.DG(graphDatabase=[reactant, product])
	with dg.build() as build:
		derivation = mod.Derivation()
		derivation.left = [reactant]
		derivation.right = [product]
		derivation.rule = rule
		build.addDerivation(derivation)

	edge = next(iter(dg.edges))
	vertex_maps = list(mod.DGVertexMapper(edge, upToIsomorphismGDH=True, rightLimit=1))
	assert len(vertex_maps) == 1
	vertex_map = vertex_maps[0]

	graph_charges = {
		reactant: ap.computeAllCharges(reactant),
		product: ap.computeAllCharges(product),
	}
	atom_values = get_atom_values(graph_charges)
	expected, id_maps = find_expected(vertex_map, atom_values)
	print(id_maps)
	print(graph_charges)
	rule_data = ap.chargeSeparation()
	rule_data.atomVals = AtomValues(atom_values)
	actual = rule_data.eval(vertex_map.rule, vertex_map, "stadler")

	assert math.isclose(actual, expected, abs_tol=TOLERANCE), (actual, expected)



def test_custom_reaction_with_multiple_mappings_possible():
	reactant_one = mod.Graph.fromSMILES("[CH2+]O[CH2+]")
	reactant_two = mod.Graph.fromSMILES("[CH2-]CC(=O)[CH2-]")
	product = mod.Graph.fromSMILES("[CH2+]OCCCC(=O)[CH2-]")
	rule = mod.Rule.fromGMLString("""rule [
		ruleID "string formDoubleBond"
		left [
			node [ id 2 label "C+" ]
			node [ id 3 label "C-" ]
		]
		right [
			node [ id 2 label "C" ]
			node [ id 3 label "C" ]
			edge [ source 2 target 3 label "-" ]
		]
	]""", add=False)

	dg = mod.DG(graphDatabase=[reactant_one, reactant_two, product])
	with dg.build() as build:
		derivation = mod.Derivation()
		derivation.left = [reactant_one, reactant_two]
		derivation.right = [product]
		derivation.rule = rule
		build.addDerivation(derivation)

	edge = next(iter(dg.edges))
	vertex_maps = list(mod.DGVertexMapper(edge, upToIsomorphismGDH=True, rightLimit=1))
	assert len(vertex_maps) == 1
	vertex_map = vertex_maps[0]

	graph_charges = {
		reactant_one: ap.computeAllCharges(reactant_one),
		reactant_two: ap.computeAllCharges(reactant_two),
		product: ap.computeAllCharges(product),
	}
	atom_values = get_atom_values(graph_charges)
	expected, id_maps = find_expected(vertex_map, atom_values)
	print(id_maps)
	print(graph_charges)
	rule_data = ap.chargeSeparation()
	rule_data.atomVals = AtomValues(atom_values)
	actual = rule_data.eval(vertex_map.rule, vertex_map, "stadler")

	assert math.isclose(actual, expected, abs_tol=TOLERANCE), (actual, expected)




##################
# RUNNING TESTS: #
##################

test_stadler_break_single_bond_on_larger_molecule()
test_stadler_form_single_bond_between_larger_fragments()
test_stadler_break_double_bond_in_exm_like_chain()
test_stadler_form_double_bond_in_exm_like_chain()
test_custom_reaction_with_multiple_mappings_possible()
print("All tests complete without fail!")