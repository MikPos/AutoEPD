import math

import ap
import mod
import test_rdkit_conversion as rd


TOLERANCE = 1e-12


def graph_from_gml(name, nodes, edges):
	node_lines = "\n".join(
		'node [ id %d label "%s" ]' % (node_id, label)
		for node_id, label in nodes
	)
	edge_lines = "\n".join(
		'edge [ source %d target %d label "%s" ]' % (source, target, label)
		for source, target, label in edges
	)
	return mod.Graph.fromGMLString("""graph [
%s
%s
]""" % (node_lines, edge_lines), name=name, add=False)


def assert_charge_dict_close(actual, expected):
	assert set(actual) == set(expected)
	for atom_id, expected_charge in expected.items():
		assert math.isclose(actual[atom_id], expected_charge, abs_tol=TOLERANCE), (
			atom_id,
			actual[atom_id],
			expected_charge,
		)


def assert_total_charge(actual, expected_total):
	assert math.isclose(sum(actual.values()), expected_total, abs_tol=TOLERANCE)


def methane():
	return graph_from_gml(
		"methane",
		[(1, "C"), (2, "H"), (3, "H"), (4, "H"), (5, "H")],
		[(1, 2, "-"), (1, 3, "-"), (1, 4, "-"), (1, 5, "-")],
	)


def water():
	return graph_from_gml(
		"water",
		[(1, "O"), (2, "H"), (3, "H")],
		[(1, 2, "-"), (1, 3, "-")],
	)


def formaldehyde():
	return graph_from_gml(
		"formaldehyde",
		[(1, "C"), (2, "O"), (3, "H"), (4, "H")],
		[(1, 2, "="), (1, 3, "-"), (1, 4, "-")],
	)


def test_methane_partial_charges_match_expected_values():
	charges = ap.computeAllCharges(methane())
	assert_charge_dict_close(charges, {
		"0": -0.079,
		"1": 0.020,
		"2": 0.020,
		"3": 0.020,
		"4": 0.020,
	})
	assert_total_charge(charges, 0.001)


def test_water_partial_charges_match_expected_values():
	charges = ap.computeAllCharges(water())
	assert_charge_dict_close(charges, {
		"0": -0.339,
		"1": 0.169,
		"2": 0.169,
	})
	assert_total_charge(charges, -0.001)


def test_formaldehyde_partial_charges_match_expected_values():
	charges = ap.computeAllCharges(formaldehyde())
	assert_charge_dict_close(charges, {
		"0": 0.050,
		"1": -0.221,
		"2": 0.086,
		"3": 0.086,
	})
	assert_total_charge(charges, 0.001)

##################
# RUNNING TESTS: #
##################

test_methane_partial_charges_match_expected_values()
test_water_partial_charges_match_expected_values()
test_formaldehyde_partial_charges_match_expected_values()
print("All tests complete without fail!")