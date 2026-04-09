import mod

molecule_1 = mod.Graph.fromSMILES("N#[C-]")
molecule_2 = mod.Graph.fromSMILES("CI")

molecule_1.print()
molecule_2.print()

