import mod

breakSingleBond = mod.Rule.fromGMLString("""rule [
	ruleID "breakSingleBond"
	left [
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		edge [ source 0 target 1 label "-" ]
	]
	right [
		node [ id 0 label "C-" ]
		node [ id 1 label "C+" ]
	]
]""", add=False)

breakDoubleBond = mod.Rule.fromGMLString("""rule [
	ruleID "breakDoubleBond"
	left [
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		edge [ source 0 target 1 label "=" ]
	]
	right [
		node [ id 0 label "C-" ]
		node [ id 1 label "C+" ]
        edge [ source 0 target 1 label "-" ]
	]
]""", add=False)

breakSingleBond.print()
breakDoubleBond.print()