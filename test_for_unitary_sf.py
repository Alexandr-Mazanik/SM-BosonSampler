import main
from strawberryfields.ops import *
import strawberryfields as sf
import numpy as np


def test():
    num_of_modes = 6

    # calculate unitary using strawberry fields
    test_unitary = sf.Program(num_of_modes)
    with test_unitary.context as q:
        BSgate() | (q[0], q[1])
        BSgate() | (q[1], q[2])
        BSgate() | (q[0], q[2])
        BSgate() | (q[2], q[3])
        BSgate() | (q[1], q[3])
        BSgate() | (q[3], q[4])
        BSgate() | (q[2], q[4])
        BSgate() | (q[4], q[5])
        BSgate() | (q[3], q[5])
        BSgate() | (q[0], q[5])

    test_unitary_compiled = test_unitary.compile(compiler="gaussian_unitary")

    s = test_unitary_compiled.circuit[0].op.p[0]
    # print(s)
    u = s[:6, :6] + 1j * s[6:, :6]
    print("\nstrawberry fields U:\n", np.round(u, 3))

    # and our unitary
    sampler_unitary_test = main.BosonSampler(num_of_modes)

    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (1, 2))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (2, 3))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (1, 3))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (3, 4))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (2, 4))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (4, 5))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (3, 5))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (5, 6))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (4, 6))
    sampler_unitary_test.add_beam_splitter(np.pi / 4, 0, 0, (1, 6))

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
