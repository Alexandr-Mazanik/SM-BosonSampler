from strawberryfields.ops import *
import strawberryfields as sf

hom_effect = sf.Program(2)

with hom_effect.context as q:
    Fock(1) | q[0]
    Fock(1) | q[1]

    BSgate() | (q[0], q[1])

    MeasureFock() | q

results = []
j = 0
batch_size = 1000
eng = sf.Engine(backend="fock", backend_options={"cutoff_dim": 3})
for i in range(batch_size):
    results.append(eng.run(hom_effect))
    if results[i].samples[0][1] == 2:
        j += 1
print("prob = ", j / batch_size)

# eng.print_applied()
# probs = results.state.all_fock_probs()
# print("\nprobs:\n", results.state.all_fock_probs())
# print("\nmeasurement:\n", results.samples)

# print("\nResults:")
# print("2, 0 -- ", probs[2, 0])
# print("0, 2 -- ", probs[0, 2])
# print("1, 1 -- ", probs[1, 1])
