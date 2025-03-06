import ap
import mod

#############################
# Different Exercises Input #
#############################
# Only one input can be used at a time.
# Comment in or out the exercise you would like to test before running the code.
# Make sure only one is active.
#
# Exercises I, J and K are skipped so far due to lack of Radical support.


# Exercise A
exercise_1_start_1 = mod.smiles("N#[C-]", "a_start_1")
exercise_1_start_2 = mod.smiles("CI", "a_start_2")
exercise_1_end_1 = mod.smiles("N#CC", "a_end_1")
exercise_1_end_2 = mod.smiles("[I-]",  "a_end_2")
a_start_1 = ap.termFromGraph(exercise_1_start_1)
a_start_2 = ap.termFromGraph(exercise_1_start_2)
a_end_1 = ap.termFromGraph(exercise_1_end_1)
a_end_2 = ap.termFromGraph(exercise_1_end_2)
nameSet = set([a_start_1, a_start_2, a_end_1, a_end_2])
sources = [a_start_1, a_start_2]
targets = [a_end_1, a_end_2]

# Exercise B
# exercise_2_start_1 = mod.smiles("C[NH-]")
# exercise_2_start_2 = mod.smiles("CC(OC)=O")
# exercise_2_end = mod.smiles("CNC(OC)(C)[O-]", "Din Mor")
# b_start_1 = ap.termFromGraph(exercise_2_start_1)
# b_start_2 = ap.termFromGraph(exercise_2_start_2)
# b_end = ap.termFromGraph(exercise_2_end)
# nameSet = set([b_start_1, b_start_2, b_end])
# sources = [b_start_1, b_start_2]
# targets = [b_end]

# Exercise C
# exercise_3_start = mod.smiles("CNC(OC)(C)[O-]")
# exercise_3_end_1 = mod.smiles("CC(=O)NC")
# exercise_3_end_2 = mod.smiles("[O-]C")
# c_start = ap.termFromGraph(exercise_3_start)
# c_end_1 = ap.termFromGraph(exercise_3_end_1)
# c_end_2 = ap.termFromGraph(exercise_3_end_2)
# nameSet = set([c_start, c_end_1, c_end_2])
# sources = [c_start]
# targets = [c_end_1, c_end_2]

# Exercise D
# # exercise_4_start_1 = mod.smiles("CN", "ds1")
# exercise_4_start_1 = mod.smiles("C[NH-]") # This is the wrong SMILES string, but the original one fails.
# exercise_4_start_2 = mod.smiles("CCCCl", "ds2")
# # exercise_4_end_1 = mod.smiles("CCC[NH2+]C", "de1")
# exercise_4_end_1 = mod.smiles("CCCNC") # Same thing here. The reason it fails is that it cannot add a bond between two atoms wothout both of them having a charge, and also can't remove unless they don't.
# exercise_4_end_2 = mod.smiles("[Cl-]")
# d_start_1 = ap.termFromGraph(exercise_4_start_1)
# d_start_2 = ap.termFromGraph(exercise_4_start_2)
# d_end_1 = ap.termFromGraph(exercise_4_end_1)
# d_end_2 = ap.termFromGraph(exercise_4_end_2)
# nameSet = set([d_start_1, d_start_2, d_end_1, d_end_2])
# sources = [d_start_1, d_start_2]
# targets = [d_end_1, d_end_2]

# Exercise E
# exercise_5_start_1 = mod.smiles("CC(O)[CH2-]")
# exercise_5_start_2 = mod.smiles("O=CC")
# exercise_5_end = mod.smiles("CC(O)CC([O-])C")
# e_start_1 = ap.termFromGraph(exercise_5_start_1)
# e_start_2 = ap.termFromGraph(exercise_5_start_2)
# e_end = ap.termFromGraph(exercise_5_end)
# nameSet = set([e_start_1, e_start_2, e_end])
# sources = [e_start_1, e_start_2]
# targets = [e_end]

# Exercise F
# # exercise_6_start_1 = mod.smiles("CC(C)(C)[O-]", "ex6_start1") # This can't be done, since we only handle forming bonds between two atoms that are both charged.
# exercise_6_start_1 = mod.smiles("CC(C)(C)O", "ex6_start1")
# exercise_6_start_2 = mod.smiles("[H+]", "ex6_start2")
# # exercise_6_end = mod.smiles("CC(C)(C)O", "ex6_end")
# exercise_6_end = mod.smiles("CC(C)(C)[OH2+]", "ex6_end")
# f_start_1 = ap.termFromGraph(exercise_6_start_1)
# f_start_2 = ap.termFromGraph(exercise_6_start_2)
# f_end = ap.termFromGraph(exercise_6_end)
# nameSet = set([f_start_1, f_start_2, f_end])
# sources = [f_start_1, f_start_2]
# targets = [f_end]


# Exercise G
# exercise_7_start = mod.smiles("CC(C)(C)[OH2+]")
# exercise_7_end_1 = mod.smiles("C[C+](C)C")
# exercise_7_end_2 = mod.smiles("O")
# g_start = ap.termFromGraph(exercise_7_start)
# g_end_1 = ap.termFromGraph(exercise_7_end_1)
# g_end_2 = ap.termFromGraph(exercise_7_end_2)
# nameSet = set([g_start, g_end_1, g_end_2])
# sources = [g_start]
# targets = [g_end_1, g_end_2]

# Exercise H
# exercise_8_start_1 = mod.smiles("COC(=O)C(C(=O)OC)CCl")
# exercise_8_start_2 = mod.smiles("CC(C)(C)[O-]")
# exercise_8_end_1 = mod.smiles("COC(=O)[C-](C(=O)OC)CCl")
# exercise_8_end_2 = mod.smiles("CC(C)(C)O")
# h_start_1 = ap.termFromGraph(exercise_8_start_1)
# h_start_2 = ap.termFromGraph(exercise_8_start_2)
# h_end_1 = ap.termFromGraph(exercise_8_end_1)
# h_end_2 = ap.termFromGraph(exercise_8_end_2)
# nameSet = set([h_start_1, h_start_2, h_end_1, h_end_2])
# sources = [h_start_1, h_start_2]
# targets = [h_end_1, h_end_2]

# Exercise M
# exercise_13_start = mod.smiles("COC(=C)OCC=C")
# exercise_13_end = mod.smiles("COC(=O)CCC=C", "end")
# m_start = ap.termFromGraph(exercise_13_start)
# m_end = ap.termFromGraph(exercise_13_end)
# nameSet = set([m_start, m_end])
# sources = [m_start]
# targets = [m_end]

# Exercise N
# exercise_14_start = mod.smiles("C1C[C+]C(CCC=CCCC2=CCCCC2)CC1")
# exercise_14_end = mod.smiles("C1CC2C(CC1)CCC3C(2)CC[C+]4C(3)CCCC4")
# n_start = ap.termFromGraph(exercise_14_start)
# n_end = ap.termFromGraph(exercise_14_end)
# nameSet = set([n_start, n_end])
# sources = [n_start]
# targets = [n_end]


####################
# Code for Testing #
####################

config.flow.computeOldSolutions = False

print(mod.getAvailableILPSolvers())


ruleData = ap.chargeSeparation()
dgData = ap.makeDG(
    rules=ruleData.rules,
    sources=sources,
    graphDatabase=mod.inputGraphs + sources + targets,
    sizeLimit=3,
    iterationLimit=2
    )

for v in dgData.dgString.vertices:
    g = v.graph
    for a in nameSet:
        if a.isomorphism(g, labelSettings=ap.lsString) > 0:
            g.name = a.name
            break
print("|V|:", dgData.dg.numVertices)
print("|E|:", dgData.dg.numEdges)
for v in dgData.dgString.vertices:
    if sum(1 for a in v.graph.vertices if a.atomId == mod.AtomIds.C) != 4: continue
    if sum(1 for a in v.graph.vertices if a.atomId == mod.AtomIds.H) != 6: continue
    if v.graph.numVertices == v.graph.numEdges: continue


ruleData.atomVals = ap.loadPartialCharges(dgData)

flowData = ap.calcPathways(
    ruleData=ruleData, dgData=dgData,
    sources=sources, targets=targets
)
flow = flowData.flow
flow.addEnumerationVar(mod.isEdgeUsed)

flow.findSolutions(verbosity=1, maxNumSolutions=1)
flow.solutions.list()

files = ap.printSolutions(ruleData=ruleData, dgData=dgData, flowData=flowData)
for f in files:
    print("DG:", f)