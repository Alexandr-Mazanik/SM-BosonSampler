import numpy as np
import time
import itertools
import os
from scipy import sparse
from thewalrus import perm


class BeamSplitter:
    def __init__(self, theta, phi_rho, phi_tau):
        self._t = np.exp(phi_tau * 1j) * np.cos(theta)
        self._r = np.exp(phi_rho * 1j) * np.sin(theta)

        self._mode1 = self._mode2 = None
        self.bs_matrix = None

    def calc_bs_matrix(self, number_of_modes, modes):
        list_modes = np.array([modes[0] - 1, modes[1] - 1])
        list_data = np.array([self._t - 1, np.conj(self._t) - 1, -np.conj(self._r), self._r])

        row_ind = np.concatenate([np.array(range(number_of_modes)), list_modes, list_modes])
        col_ing = np.concatenate([np.array(range(number_of_modes)), list_modes, np.flip(list_modes)])
        data = np.concatenate([np.ones((number_of_modes, ), dtype=complex), list_data])

        self.bs_matrix = sparse.csr_matrix((data, (row_ind, col_ing)), shape=(number_of_modes, number_of_modes))


class Block:
    def __init__(self, modes_num):
        self.number_of_modes = modes_num
        self._beam_splitters = []
        self.block_matrix = sparse.csr_matrix(np.identity(self.number_of_modes))

    def add_BS_gate(self, modes, theta=np.pi/4, phi_rho=0., phi_tau=0.):
        bs = BeamSplitter(theta, phi_rho, phi_tau)
        bs.calc_bs_matrix(self.number_of_modes, modes)
        self._beam_splitters.append(bs)

    def calc_block_matrix(self):
        time_unit_start = time.time()
        self.block_matrix = sparse.csr_matrix(np.identity(self.number_of_modes))
        for bs in np.flip(self._beam_splitters):
            self.block_matrix = sparse.csr_matrix.dot(self.block_matrix, bs.bs_matrix)
        self.block_matrix = self.block_matrix.toarray()
        time_unit_end = time.time()
        print("--> The time for dot product is :", (time_unit_end - time_unit_start) * 10 ** 3, "ms")


class Scheme:
    def __init__(self):
        self._blocks = []

    def upload_scheme_from_file(self, file_name='curr_scheme.txt'):
        ref_mode = 1
        with open(os.path.join('scheme', file_name), 'r') as f_scheme:
            for f_line in f_scheme:
                if f_line.strip() == "BLOCK":
                    modes_num, repetitions_num = list(map(int, f_scheme.readline().split()))
                    block = Block(modes_num)
                    for f_beam_splitter in f_scheme:
                        if f_beam_splitter.strip() == "END":
                            block.calc_block_matrix()
                            for _ in range(repetitions_num):
                                self._blocks.append({'block': block,
                                                     'ref_mode': ref_mode})
                                ref_mode += 1
                            break
                        mode1, mode2 = list(map(int, f_beam_splitter.split('\t')[0:2]))
                        theta, phi_rho, phi_tau = list(map(float, f_beam_splitter.split('\t')[2:]))
                        block.add_BS_gate((mode1, mode2), theta, phi_rho, phi_tau)

        print("--> Scheme was successfully uploaded")


def main():
    scheme = Scheme()
    scheme.upload_scheme_from_file()


if __name__ == '__main__':
    main()
