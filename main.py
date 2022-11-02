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
        time_unit_start = time.time()
        self.unitary = multi_dot([np.identity(self.number_of_modes)] +
                                 [bs.unitary for bs in np.flip(self._beam_splitters)])
        time_unit_end = time.time()
        print("--> The time for dot product is :", (time_unit_end - time_unit_start) * 10 ** 3, "ms")

    def upload_scheme_from_file(self):
        with open('scheme/curr_scheme.txt', 'r') as f_scheme:
            self.number_of_modes = int(f_scheme.readline())
            for f_beam_splitter in f_scheme:
                mode1, mode2 = list(map(int, f_beam_splitter.split('\t')[0:2]))
                theta, phi_rho, phi_tau = list(map(float, f_beam_splitter.split('\t')[2:]))
                self.add_BS_gate((mode1, mode2), theta, phi_rho, phi_tau)
        print("--> Scheme was successfully uploaded")

    def export_system_unitary(self):
        with open('scheme/scheme_unitary.txt', 'w') as f_out:
            for line in np.round(self.unitary, 6):
                for element in line:
                    f_out.write(str(element.real) + '\t' + str(element.imag) + '\t')
                f_out.write('\n')
        print("--> Unitary was successfully exported")

    def print_system_unitary(self):
        print("--> U:\n", np.round(self.unitary, 3))


def is_unitary(matrix, dim):
    matrix_dagger = np.conj(matrix.transpose())
    # print(np.round(np.dot(matrix, matrix_dagger)))
    if (np.round(np.dot(matrix, matrix_dagger), 10) == np.identity(dim)).all():
        print("--> is_unitary: True")
    else:
        print("--> is_unitary: False")


def main():
    time_start = time.time()

    sampler = BosonSampler(1)
    sampler.upload_scheme_from_file()
    sampler.calc_system_unitary()
    sampler.export_system_unitary()

    is_unitary(sampler.unitary, sampler.number_of_modes)

    time_end = time.time()
    print("\n--> The time of execution is :", (time_end - time_start) * 10 ** 3, "ms")


if __name__ == '__main__':
    main()
