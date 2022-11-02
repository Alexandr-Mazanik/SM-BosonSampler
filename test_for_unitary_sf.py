import main
from strawberryfields.ops import *
import strawberryfields as sf
import numpy as np


def test():
    num_of_modes = 3

    # calculate unitary using strawberry fields
    test_unitary = sf.Program(num_of_modes)
    with test_unitary.context as q:
        BSgate(theta=1.12333, phi=2.31234) | (q[0], q[1])
        BSgate(theta=np.pi/10, phi=0.12334) | (q[0], q[2])

    test_unitary_compiled = test_unitary.compile(compiler="gaussian_unitary")

    s = test_unitary_compiled.circuit[0].op.p[0]
    # print(s)
    u = s[:num_of_modes, :num_of_modes] + 1j * s[num_of_modes:, :num_of_modes]
    print("\nstrawberry fields U:\n", np.round(u, 3))

    # and our unitary
    sampler_unitary_test = main.BosonSampler(num_of_modes)

    sampler_unitary_test.add_BS_gate((1, 2), theta=1.12333, phi_rho=2.31234)
    sampler_unitary_test.add_BS_gate((1, 3), theta=np.pi/10, phi_rho=0.12334)

    sampler_unitary_test.calc_system_unitary()
    sampler_u = sampler_unitary_test.unitary
    print("\nsampler U:\n", np.round(sampler_u, 3))

    # compare them
    if np.array_equal(np.round(sampler_u, 4), np.round(u, 4)):
        print("\nTrue\n")
    else:
        print("\nFalse\n")


if __name__ == '__main__':
    test()
