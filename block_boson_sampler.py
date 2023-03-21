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
        data = np.concatenate([np.ones((number_of_modes,), dtype=complex), list_data])

        self.bs_matrix = sparse.csr_matrix((data, (row_ind, col_ing)), shape=(number_of_modes, number_of_modes))


class Block:
    def __init__(self, modes_num):
        self.number_of_modes = modes_num
        self._beam_splitters = []
        self.block_matrix = sparse.csr_matrix(np.identity(self.number_of_modes))

    def add_BS_gate(self, modes, theta=np.pi / 4, phi_rho=0., phi_tau=0.):
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
        self.blocks = []

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
                                self.blocks.append({'block': block,
                                                    'ref_mode': ref_mode})
                                ref_mode += 1
                            break
                        mode1, mode2 = list(map(int, f_beam_splitter.split('\t')[0:2]))
                        theta, phi_rho, phi_tau = list(map(float, f_beam_splitter.split('\t')[2:]))
                        block.add_BS_gate((mode1, mode2), theta, phi_rho, phi_tau)

        print("--> Scheme was successfully uploaded")


class BosonSampler:
    def __init__(self, scheme, init_config):
        self._scheme = scheme
        self._basis = None
        self._state = None
        self.init_config = init_config

    @staticmethod
    def create_fock_basis(ph_num, modes_num):
        basis = []
        slots_num = ph_num + modes_num
        all_comb_bars = list(itertools.combinations(range(1, slots_num), modes_num - 1))
        for bars in all_comb_bars:
            bars = list(bars)
            bars.append(slots_num)
            bars.insert(0, 0)
            basis_vec = []
            for i in range(modes_num):
                basis_vec.append(bars[i + 1] - bars[i] - 1)
            basis.append(basis_vec)

        return basis

    @staticmethod
    def block_evolve(block, input_state, basis):
        def find_submatrix():
            def get_indices(vec):
                indices = []
                for i in range(block['block'].number_of_modes):
                    if vec[i] > 0:
                        for j in range(vec[i]):
                            indices.append(i)
                return indices

            column_indices = get_indices(basis[vec_in])
            row_indices = get_indices(basis[vec_out])

            return block['block'].block_matrix[:, column_indices][row_indices]

        dim = len(input_state)

        input_state = np.array(input_state)
        transform_matrix = np.empty([dim, dim], dtype=complex)

        for vec_in in range(dim):
            for vec_out in range(dim):
                norm = 1
                for num in basis[vec_in]:
                    norm *= np.math.sqrt(np.math.factorial(num))
                for num in basis[vec_out]:
                    norm *= np.math.sqrt(np.math.factorial(num))

                transform_matrix[vec_in][vec_out] = \
                    perm(find_submatrix(), method="ryser") / norm

        output_state = np.matmul(input_state, transform_matrix)

        return output_state

    @staticmethod
    def collapse_state(state, basis, observation):
        norm_sq = 0
        new_basis = []
        new_state = []
        for i, vec in enumerate(basis):
            if vec[0] == observation:
                new_basis.append(vec)
                new_state.append(state[i])
                norm_sq += abs(state[i]) ** 2
        norm = np.sqrt(norm_sq)
        for i in range(len(new_state)):
            new_state[i] = new_state[i] / norm

        return new_state, new_basis

    def calculate_one_sample(self):
        sample = []
        ph_num = 0
        prev_block = None
        self._basis = []
        self._state = []
        for block_num, block in enumerate(self._scheme.blocks):
            if block_num == 0:
                ph_num = sum(self.init_config[:block['block'].number_of_modes])
            else:
                prev_block = self._scheme.blocks[block_num - 1]
                ph_num = (ph_num -
                          sample[block['ref_mode'] - 2] +
                          sum(self.init_config[prev_block['block'].number_of_modes + prev_block['ref_mode'] - 1:
                                               block['block'].number_of_modes + block['ref_mode'] - 1]))

            basis = self.create_fock_basis(ph_num, block['block'].number_of_modes)

            if block_num == 0:
                state = [1 if i == basis.index(list(self.init_config[:block['block'].number_of_modes])) else 0
                         for i in range(len(basis))]
            else:
                state = [0 for _ in range(len(basis))]
                for i, vec in enumerate(self._basis):
                    new_vec = vec[1:] + self.init_config[
                                        prev_block['block'].number_of_modes + prev_block['ref_mode'] - 1:
                                        block['block'].number_of_modes + block['ref_mode'] - 1]

                    j = basis.index(new_vec)
                    state[j] = self._state[i]

            state = self.block_evolve(block, state, basis)

            prob = [abs(ampl) ** 2 for ampl in state]
            if block_num == len(self._scheme.blocks)-1:
                out_index = np.random.choice(len(basis), p=prob)
                sample += basis[out_index]
            else:
                n_ph_observation = np.zeros([ph_num + 1], dtype=float)
                for i, vec in enumerate(basis):
                    n_ph_observation[vec[0]] += prob[i]

                observation = np.random.choice(np.arange(ph_num + 1), p=n_ph_observation)
                self._state, self._basis = self.collapse_state(state, basis, observation)
                sample.append(observation)

        print("SAMPLE: ", sample)
        return sample

    def sample(self, batch_size=100, file_name='block_sample.txt'):
        with open(os.path.join('samples', file_name), 'w') as f_out:
            for _ in range(batch_size):
                sample = self.calculate_one_sample()
                f_out.write(str(sample) + '\n')
        print("--> Samples successfully exported")


def main():
    scheme = Scheme()
    scheme.upload_scheme_from_file()

    modes_num = len(scheme.blocks) + scheme.blocks[-1]['block'].number_of_modes - 1
    ph_num = 7

    init_config = [1 if i < ph_num else 0 for i in range(modes_num)]

    sampler = BosonSampler(scheme, init_config)
    sampler.sample(batch_size=100000)


if __name__ == '__main__':
    main()
