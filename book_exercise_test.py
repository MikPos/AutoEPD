import ap_new as ap
import time
import itertools
import mod

SIZE_LIMIT = 3
ITERATION_LIMIT = 2
FILE_NAME = "networkx.txt"

timing_dictionary = {}
with open(FILE_NAME, "r") as f:
    for line in f:
        (key, value) = line.split(" : ")
        timing_dictionary[key] = value

overall_time_start = time.perf_counter()

#############################
# Different Exercises Input #
#############################
# Only one input can be used at a time.
# Comment in or out the exercise you would like to test before running the code.
# Make sure only one is active.
#
# Exercises I, J and K are skipped so far due to lack of Radical support.


# Exercise A
# exercise_1_start_1 = mod.smiles("N#[C-]", "a start 1")
# exercise_1_start_2 = mod.smiles("CI", "a start 2")
# exercise_1_end_1 = mod.smiles("N#CC", "a end 1")
# exercise_1_end_2 = mod.smiles("[I-]",  "a end 2")
# a_start_1 = ap.termFromGraph(exercise_1_start_1)
# a_start_2 = ap.termFromGraph(exercise_1_start_2)
# a_end_1 = ap.termFromGraph(exercise_1_end_1)
# a_end_2 = ap.termFromGraph(exercise_1_end_2)
# nameSet = set([a_start_1, a_start_2, a_end_1, a_end_2])
# sources = [a_start_1, a_start_2]
# targets = [a_end_1, a_end_2]
# timing_keyword = "Exercise A"


# Exercise B
# exercise_2_start_1 = mod.smiles("C[NH-]")
# exercise_2_start_2 = mod.smiles("CC(OC)=O")
# exercise_2_end = mod.smiles("CNC(OC)(C)[O-]")
# b_start_1 = ap.termFromGraph(exercise_2_start_1)
# b_start_2 = ap.termFromGraph(exercise_2_start_2)
# b_end = ap.termFromGraph(exercise_2_end)
# nameSet = set([b_start_1, b_start_2, b_end])
# sources = [b_start_1, b_start_2]
# targets = [b_end]
# timing_keyword = "Exercise B"

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
# timing_keyword = "Exercise C"

# Exercise D
# exercise_4_start_1 = mod.smiles("CN", "ds1")
# exercise_4_start_2 = mod.smiles("CCCCl", "ds2")
# exercise_4_end_1 = mod.smiles("CCC[NH2+]C", "de1")
# exercise_4_end_2 = mod.smiles("[Cl-]")
# d_start_1 = ap.termFromGraph(exercise_4_start_1)
# d_start_2 = ap.termFromGraph(exercise_4_start_2)
# d_end_1 = ap.termFromGraph(exercise_4_end_1)
# d_end_2 = ap.termFromGraph(exercise_4_end_2)
# nameSet = set([d_start_1, d_start_2, d_end_1, d_end_2])
# sources = [d_start_1, d_start_2]
# targets = [d_end_1, d_end_2]
# timing_keyword = "Exercise D"

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
# timing_keyword = "Exercise E"

# Exercise F
# exercise_6_start_1 = mod.smiles("CC(C)(C)O", "f start 1")
# exercise_6_start_2 = mod.smiles("[H+]", "f start 2")
# exercise_6_end = mod.smiles("CC(C)(C)[OH2+]", "f end")
# f_start_1 = ap.termFromGraph(exercise_6_start_1)
# f_start_2 = ap.termFromGraph(exercise_6_start_2)
# f_end = ap.termFromGraph(exercise_6_end)
# nameSet = set([f_start_1, f_start_2, f_end])
# sources = [f_start_1, f_start_2]
# targets = [f_end]
# timing_keyword = "Exercise F"


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
# timing_keyword = "Exercise G"

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
# timing_keyword = "Exercise H"

# Exercise M
# exercise_13_start = mod.smiles("COC(=C)OCC=C")
# exercise_13_end = mod.smiles("COC(=O)CCC=C", "end")
# m_start = ap.termFromGraph(exercise_13_start)
# m_end = ap.termFromGraph(exercise_13_end)
# nameSet = set([m_start, m_end])
# sources = [m_start]
# targets = [m_end]
# timing_keyword = "Exercise M"

# Exercise N
exercise_14_start = mod.smiles("C1C[CH+]C(CCC=CCCC2=CCCCC2)CC1")
exercise_14_end = mod.smiles("C1CC2C(CC1)CCC3C2CC[C+]4C3CCCC4")
n_start = ap.termFromGraph(exercise_14_start)
n_end = ap.termFromGraph(exercise_14_end)
nameSet = set([n_start, n_end])
sources = [n_start]
targets = [n_end]
timing_keyword = "Exercise N"


####################
# Code for Testing #
####################


# Get size of molecules
molecule_vertices = 0
molecules_edges = 0
for graph in sources:

    for v in graph.vertices:
        molecule_vertices += 1
        for e in v.incidentEdges:
            if e.target.id < e.source.id:
                molecules_edges += 1

# config.flow.computeOldSolutions = False

print(mod.getAvailableILPSolvers())


ruleData = ap.chargeSeparation()

# Testing Runtime of Constructing DG:
dg_start_time = time.perf_counter()

# dgData, nxPath = ap.makeDG(
dgData = ap.makeDG(
    rules=ruleData.rules,
    sources=sources,
    targets=targets,
    graphDatabase=mod.inputGraphs + sources + targets,
    sizeLimit=SIZE_LIMIT,
    iterationLimit=ITERATION_LIMIT
    )

dg_stop_time = time.perf_counter()

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


# allowed_graphs = list(itertools.chain(*nxPath))
# print(allowed_graphs)

# Test runtime of partial charges:
ruleData.atomVals = ap.loadPartialCharges(dgData)

# Test Runtime of Flow Pathway consruction:
flow_start_time = time.perf_counter()

flowData = ap.calcPathways(
    ruleData=ruleData, dgData=dgData,
    sources=sources, targets=targets,
    # pathGraphs=allowed_graphs
)
flow = flowData.flow
flow.addEnumerationVar(mod.isEdgeUsed)

# Why does it keep being shit?
flow.findSolutions(verbosity=1, maxNumSolutions=1)
flow.solutions.list()

flow_stop_time = time.perf_counter()

overall_time_stop = time.perf_counter()

files = ap.printSolutions(ruleData=ruleData, dgData=dgData, flowData=flowData)
for f in files:
    print("DG:", f)

# for sol in flow.solutions:
#     sol.print()
# p = mod.GraphPrinter()
# usedGraphs = []
# for solution in flow.solutions:
#     for e in dgData.dg.edges:
#         if solution.eval(edge[e]) != 0:
#             for rule in e.rules:
#                 rule.print(p)
            # vms = mod.DGVertexMapper(e, upToIsomorphismGDH = True, rightLimit = 1)
            # for vm in vms:
# for g in usedGraphs:
#     print(g)


# print(ap.rules)
# Compute the runtime:
overall_time = overall_time_stop - overall_time_start
dg_time = dg_stop_time - dg_start_time
flow_time = flow_stop_time - flow_start_time

# Update data with new times:
# print(timing_dictionary)
# exercise_timing_list = timing_dictionary[timing_keyword].strip().split(",")
exercise_timing_list = []

exercise_timing_list.append(f"{molecule_vertices},{molecules_edges},{dgData.dg.numVertices},{dgData.dg.numEdges},{overall_time},{dg_time},{flow_time}")
# exercise_timing_list.append(f"{overall_time},{dg_time},{fp_time},{fs_time}")
exercise_timing_string = "".join(exercise_timing_list)
timing_dictionary[timing_keyword] = exercise_timing_string + "\n"

file = open(FILE_NAME, "w")

for key, value in timing_dictionary.items():
    # print(f"{key}, {value}")
    file.write(f"{key} : {value}")

file.close()