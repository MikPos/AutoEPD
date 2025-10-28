import ap_new as ap
import time
import itertools
import mod

SIZE_LIMIT = 3
ITERATION_LIMIT = 4
FILE_NAME = "non_daniel.txt"

timing_dictionary = {}
with open(FILE_NAME, "r") as f:
    for line in f:
        (key, value) = line.split(" : ")
        timing_dictionary[key] = value

overall_time_start = time.perf_counter()

# Test Example
formaldehyde = mod.smiles("C=O", name="Formaldehyde")
glycolaldehyde = mod.smiles( "OCC=O", name="Glycolaldehyde")

post_keto = mod.smiles("OC=CO")

post_aldol = mod.smiles("OC(C=O)CO")

f = ap.termFromGraph(formaldehyde)
g = ap.termFromGraph(glycolaldehyde)
pk = ap.termFromGraph(post_keto)
pa = ap.termFromGraph(post_aldol)

# Keto
# nameSet = set([g,pk])
# sources = [g]
# targets = [pk]
# timing_keyword = "Keto-Enol"

#Aldol
nameSet = set([f,pk,pa])
sources = [f,pk]
targets = [pa]
timing_keyword = "Aldol"

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
# flow.addEnumerationVar(mod.isEdgeUsed) # Ask Jakob why?

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