import numpy as np
import time
from scipy import sparse


class BeamSplitter:
    def __init__(self, theta, phi_rho, phi_tau):
        self._t = np.exp(phi_tau * 1j) * np.cos(theta)
        self._r = np.exp(phi_rho * 1j) * np.sin(theta)

        self._mode1 = self._mode2 = None
        self.unitary = None

    def calc_unitary(self, number_of_modes, modes):
        list_modes = np.array([modes[0] - 1, modes[1] - 1])
        list_data = np.array([self._t - 1, np.conj(self._t) - 1, -np.conj(self._r), self._r])

        row_ind = np.concatenate([np.array(range(number_of_modes)), list_modes, list_modes])
        col_ing = np.concatenate([np.array(range(number_of_modes)), list_modes, np.flip(list_modes)])
        data = np.concatenate([np.ones((number_of_modes, ), dtype=complex), list_data])

        self.unitary = sparse.csr_matrix((data, (row_ind, col_ing)), shape=(number_of_modes, number_of_modes))


class BosonSampler:
    def __init__(self, modes_num):
        self.number_of_modes = modes_num
        self._beam_splitters = []
        self.unitary = sparse.csr_matrix(np.identity(self.number_of_modes))

    def add_BS_gate(self, modes, theta=np.pi/4, phi_rho=0., phi_tau=0.):
        bs = BeamSplitter(theta, phi_rho, phi_tau)
        bs.calc_unitary(self.number_of_modes, modes)
        self._beam_splitters.append(bs)

    def calc_system_unitary(self):
        time_unit_start = time.time()
        self.unitary = sparse.csr_matrix(np.identity(self.number_of_modes))
        for bs in np.flip(self._beam_splitters):
            self.unitary = sparse.csr_matrix.dot(self.unitary, bs.unitary)
        self.unitary = self.unitary.toarray()
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
