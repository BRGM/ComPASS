import ComPASS

simulation = ComPASS.load_eos("water2ph")
print("Phases:", simulation.phases())
print("Components:", simulation.components())
print("Physical contexts:", simulation.contexts())
