import numpy as np
import time
from numpy.linalg import multi_dot


class BeamSplitter:
    def __init__(self, theta, phi_rho, phi_tau):
        self._t = np.exp(phi_tau * 1j) * np.cos(theta)
        self._r = np.exp(phi_rho * 1j) * np.sin(theta)

        self._mode1 = self._mode2 = None
        self.unitary = None

    def calc_unitary(self, number_of_modes, modes):
        self._mode1 = modes[0] - 1
        self._mode2 = modes[1] - 1
        self.unitary = np.identity(number_of_modes, dtype=complex)

        self.unitary[self._mode1, self._mode1] = self._t
        self.unitary[self._mode1, self._mode2] = -np.conj(self._r)
        self.unitary[self._mode2, self._mode1] = self._r
        self.unitary[self._mode2, self._mode2] = np.conj(self._t)


class BosonSampler:
    def __init__(self, modes_num):
        self.number_of_modes = modes_num
        self._beam_splitters = []
        self.unitary = np.identity(self.number_of_modes)

    def add_BS_gate(self, modes, theta=np.pi/4, phi_rho=0., phi_tau=0.):
        bs = BeamSplitter(theta, phi_rho, phi_tau)
        bs.calc_unitary(self.number_of_modes, modes)
        self._beam_splitters.append(bs)

    def calc_system_unitary(self):
        # TODO remove the excess identity matrix
        self.unitary = multi_dot([np.identity(self.number_of_modes)] +
                                 [bs.unitary for bs in np.flip(self._beam_splitters)])

    def upload_scheme_from_file(self):
        with open('scheme/curr_scheme.txt', 'r') as f_scheme:
            self.number_of_modes = int(f_scheme.readline())
            for f_beam_splitter in f_scheme:
                mode1, mode2 = list(map(int, f_beam_splitter.split()[0:2]))
                theta, phi_rho, phi_tau = list(map(float, f_beam_splitter.split()[2:]))
                self.add_BS_gate((mode1, mode2), theta, phi_rho, phi_tau)

    def print_system_unitary(self):
        print("U:\n", np.round(self.unitary, 3))


def is_unitary(matrix, dim):
    matrix_dagger = np.conj(matrix.transpose())
    # print(np.round(np.dot(matrix, matrix_dagger)))
    if (np.round(np.dot(matrix, matrix_dagger), 10) == np.identity(dim)).all():
        print("\nis_unitary: True")
    else:
        print("\nis_unitary: False")


def main():
    time_start = time.time()

    sampler = BosonSampler(1)

    sampler.upload_scheme_from_file()

    #sampler.add_BS_gate((1, 2), theta=0.233, phi_rho=23.231, phi_tau=0.232)
    #sampler.add_BS_gate((2, 3), theta=1.234, phi_rho=1.231, phi_tau=0.873)

    sampler.calc_system_unitary()

    sampler.print_system_unitary()
    is_unitary(sampler.unitary, sampler.number_of_modes)

    time_end = time.time()
    print("\nThe time of execution is :",
          (time_end - time_start) * 10 ** 3, "ms")


if __name__ == '__main__':
    main()
