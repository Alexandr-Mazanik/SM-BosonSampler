import main
from strawberryfields.ops import *
import strawberryfields as sf
import numpy as np
import time


def to_sf_BSgate(q, modes, theta=np.pi/4, phi_rho=0., phi_tau=0.):
    Rgate(phi_tau) | q[modes[0] - 1]
    Rgate(-phi_tau) | q[modes[1] - 1]
    BSgate(theta, phi_rho - phi_tau) | (q[modes[0] - 1], q[modes[1] - 1])


def upload_sf_system(q, f_scheme):
    for f_beam_splitter in f_scheme:
        mode1, mode2 = list(map(int, f_beam_splitter.split()[0:2]))
        theta, phi_rho, phi_tau = list(map(float, f_beam_splitter.split()[2:]))
        to_sf_BSgate(q, (mode1, mode2), theta, phi_rho, phi_tau)


def test():
    # calculate unitary using strawberry fields
    time_sf_start = time.time()

    f_scheme = open('scheme/curr_scheme.txt', 'r')
    num_of_modes = int(f_scheme.readline())

    test_unitary = sf.Program(num_of_modes)
    with test_unitary.context as q:
        upload_sf_system(q, f_scheme)
        f_scheme.close()

    test_unitary_compiled = test_unitary.compile(compiler="gaussian_unitary")

    s = test_unitary_compiled.circuit[0].op.p[0]
    sf_unitary = s[:num_of_modes, :num_of_modes] + 1j * s[num_of_modes:, :num_of_modes]

    time_sf_end = time.time()

    # calculate unitary using boson sampler
    time_bs_start = time.time()

    sampler_unitary = main.BosonSampler(num_of_modes)
    sampler_unitary.upload_scheme_from_file()
    sampler_unitary.calc_system_unitary()
    sampler_unitary = sampler_unitary.unitary

    time_bs_end = time.time()

    # compare them
    print("\nThe time of execution of SF is :", (time_sf_end - time_sf_start) * 10 ** 3, "ms")
    print("\nThe time of execution of Boson_Sampler is :", (time_bs_end - time_bs_start) * 10 ** 3, "ms")

    if np.array_equal(np.round(sampler_unitary, 5), np.round(sf_unitary, 5)):
        print("\nTrue\n")
    else:
        print("\nFalse\n")


if __name__ == '__main__':
    test()
