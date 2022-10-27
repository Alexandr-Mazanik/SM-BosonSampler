import numpy as np
import strawberryfields as sf
from strawberryfields.ops import *

np.random.seed(42)

hom_effect = sf.Program(2)

with hom_effect.context as q:
    Fock(1) | q[0]
    Fock(1) | q[1]

    BSgate() | (q[0], q[1])

    # MeasureFock() | q
    # Interferometer(U) | q

eng = sf.Engine(backend="fock", backend_options={"cutoff_dim": 3})
results = eng.run(hom_effect)

eng.print_applied()
probs = results.state.all_fock_probs()
print("\nprobs:\n", results.state.all_fock_probs())
print("\nmeasurement:\n", results.samples)

print("\nResults:")
print("2, 0 -- ", probs[2, 0])
print("0, 2 -- ", probs[0, 2])
print("1, 1 -- ", probs[1, 1])
