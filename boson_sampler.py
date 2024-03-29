import numpy as np
import strawberryfields as sf
import time
import itertools
import os
from tqdm import tqdm
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


class Scheme:
    def __init__(self, modes_num=1):
        self.number_of_modes = modes_num
        self._beam_splitters = []
        self.scheme_matrix = sparse.csr_matrix(np.identity(self.number_of_modes))

    def add_bs_gate(self, modes, theta=np.pi / 4, phi_rho=0., phi_tau=0.):
        bs = BeamSplitter(theta, phi_rho, phi_tau)
        bs.calc_bs_matrix(self.number_of_modes, modes)
        self._beam_splitters.append(bs)

    def calc_scheme_matrix(self):
        time_unit_start = time.time()
        self.scheme_matrix = sparse.csr_matrix(np.identity(self.number_of_modes))
        for bs in np.flip(self._beam_splitters):
            self.scheme_matrix = sparse.csr_matrix.dot(self.scheme_matrix, bs.bs_matrix)
        self.scheme_matrix = self.scheme_matrix.toarray()
        time_unit_end = time.time()
        print("--> The time for dot product is :", (time_unit_end - time_unit_start) * 10 ** 3, "ms")

    def haar_random_matrix(self):
        self.scheme_matrix = sf.utils.random_interferometer(self.number_of_modes)

    def upload_scheme_from_file(self, file_name):
        with open(os.path.join('scheme', file_name), 'r') as f_scheme:
            self.number_of_modes = int(f_scheme.readline())
            for f_beam_splitter in f_scheme:
                mode1, mode2 = list(map(int, f_beam_splitter.split('\t')[0:2]))
                theta, phi_rho, phi_tau = list(map(float, f_beam_splitter.split('\t')[2:]))
                self.add_bs_gate((mode1, mode2), theta, phi_rho, phi_tau)
        print("--> Scheme was successfully uploaded")
        print("--> Number of modes: ", self.number_of_modes)

    def upload_unitary(self, file_name):
        self.scheme_matrix = np.empty([self.number_of_modes, self.number_of_modes], dtype=complex)
        with open(os.path.join('scheme', file_name), 'r') as f_scheme:
            for i, f_line in enumerate(f_scheme):
                f_line = f_line.strip()
                values = list(map(float, f_line.split('\t')))
                m = 0
                j = 0
                for _ in range(self.number_of_modes):
                    self.scheme_matrix[i][m] = complex(values[j], values[j+1])
                    j += 2
                    m += 1

    def export_scheme_matrix(self, file_name):
        with open(os.path.join('scheme', file_name), 'w') as f_out:
            for line in np.round(self.scheme_matrix, 6):
                for element in line:
                    f_out.write(str(element.real) + '\t' + str(element.imag) + '\t')
                f_out.write('\n')
        print("--> Unitary was successfully exported")

    def print_scheme_matrix(self):
        print("--> U:\n", np.round(self.scheme_matrix, 3))


class BosonSampler:
    def __init__(self, scheme, init_config):
        self._scheme = scheme
        self._photons_number = sum(init_config)
        self._basis = self.create_fock_basis()
        self._dim = len(self._basis)
        self._init_config = self._basis.index(list(init_config))

        time_unit_start = time.time()
        self._transform_matrix = self.calc_transform_matrix()
        time_unit_end = time.time()
        self.time_one_sample = (time_unit_end - time_unit_start) * 10 ** 3

    def create_fock_basis(self):
        basis = []
        slots_num = self._photons_number + self._scheme.number_of_modes
        all_comb_bars = list(itertools.combinations(range(1, slots_num), self._scheme.number_of_modes - 1))
        for bars in all_comb_bars:
            bars = list(bars)
            bars.append(slots_num)
            bars.insert(0, 0)
            basis_vec = []
            for i in range(self._scheme.number_of_modes):
                basis_vec.append(bars[i+1] - bars[i] - 1)
            basis.append(basis_vec)

        return basis

    def calc_transform_matrix(self):
        def find_submatrix():
            def get_indices(vec):
                indices = []
                for i in range(self._scheme.number_of_modes):
                    if vec[i] > 0:
                        for j in range(vec[i]):
                            indices.append(i)
                return indices

            column_indices = get_indices(self._basis[vec_in])
            row_indices = get_indices(self._basis[vec_out])

            return self._scheme.scheme_matrix[:, column_indices][row_indices]

        transform_matrix = np.empty(self._dim, dtype=float)
        vec_in = self._init_config
        for vec_out in tqdm(range(self._dim), desc="Computing...", disable=False):
            norm = 1
            for num in self._basis[vec_in]:
                norm *= np.math.factorial(num)
            for num in self._basis[vec_out]:
                norm *= np.math.factorial(num)

            transform_matrix[vec_out] = \
                abs(perm(find_submatrix(), method="ryser")) ** 2 / norm

        return transform_matrix

    def sample(self, batch_size, file_name):
        prob_distribution = self._transform_matrix
        choices = [np.random.choice(range(self._dim),
                                    p=prob_distribution) for _ in tqdm(range(batch_size), desc="Sampling...")]

        with open(os.path.join('samples', file_name), 'w') as f_out:
            for result in tqdm(choices, desc="Writing to file..."):
                f_out.write(str(self._basis[result]) + '\n')

        print("--> Samples was successfully exported")

    def export_ground_truth(self, output_file_name='ground_truth'):
        output_file_name = output_file_name + '.txt'
        with open(os.path.join('samples', output_file_name), 'w') as f_out:
            for i, vec in enumerate(self._basis):
                f_out.write(str(vec) + '\t' +
                            str(np.round(self._transform_matrix[i], 6)) + '\n')

    def print_transform_matrix(self):
        print("--> H:\n", self._transform_matrix)


class DistinguishableSampler:
    def __init__(self, scheme, init_config):
        self._scheme = scheme
        self._photons_number = sum(init_config)

    def get_prob_list(self, j):
        prob_list = []
        for i in range(self._scheme.number_of_modes):
            prob_list.append(abs(self._scheme.scheme_matrix[i][j]) ** 2)

        return prob_list

    def calculate_one_sample(self):
        sample = []
        h = []
        for j in range(self._photons_number):
            h.append(np.random.choice(np.arange(self._scheme.number_of_modes), p=self.get_prob_list(j)))

        for i in range(self._scheme.number_of_modes):
            ph_num = 0
            for j, hj in enumerate(h):
                if hj == i:
                    ph_num += 1
            sample.append(ph_num)

        return sample

    def sample(self, batch_size=100, file_name='distinguishable_sample', samples_num=1):
        for i in range(samples_num):
            file_name_full = file_name + str(i + 1) + '.txt'
            with open(os.path.join('samples', file_name_full), 'w') as f_out:
                for _ in tqdm(range(batch_size), desc="Sampling..."):
                    sample = self.calculate_one_sample()
                    f_out.write(str(sample) + '\n')
                print("--> Distinguishable Samples successfully exported")


class UniformSampler:
    def __init__(self, scheme, init_config):
        self._scheme = scheme
        self._photons_number = sum(init_config)
        self._basis = self.create_fock_basis()
        self._dim = len(self._basis)

    def create_fock_basis(self):
        basis = []
        slots_num = self._photons_number + self._scheme.number_of_modes
        all_comb_bars = list(itertools.combinations(range(1, slots_num), self._scheme.number_of_modes - 1))
        for bars in all_comb_bars:
            bars = list(bars)
            bars.append(slots_num)
            bars.insert(0, 0)
            basis_vec = []
            for i in range(self._scheme.number_of_modes):
                basis_vec.append(bars[i+1] - bars[i] - 1)
            basis.append(basis_vec)

        return basis

    def sample(self, batch_size=100, file_name='uniform_sample', samples_num=1):
        for i in range(samples_num):
            choices = [np.random.choice(range(self._dim)) for _ in tqdm(range(batch_size), desc="Sampling...")]
            file_name_full = file_name + str(i + 1) + '.txt'
            with open(os.path.join('samples', file_name_full), 'w') as f_out:
                for result in tqdm(choices, desc="Sampling..."):
                    f_out.write(str(self._basis[result]) + '\n')
                print("--> Uniform Samples successfully exported")


def is_unitary(matrix, dim):
    matrix_dagger = np.conj(matrix.transpose())
    if (np.round(np.dot(matrix, matrix_dagger), 10) == np.identity(dim)).all():
        print("--> is_unitary: True")
    else:
        print("--> is_unitary: False")


def main():
    time_start = time.time()

    ph_num = 4
    modes_num = 16
    init_config = [1 if _ < ph_num else 0 for _ in range(modes_num)]

    scheme = Scheme(modes_num)
    scheme.haar_random_matrix()
    # scheme.upload_unitary('matrix.txt')
    # is_unitary(scheme.scheme_matrix, 25)
    # scheme.calc_scheme_matrix()
    scheme.export_scheme_matrix('scheme_unitary_4_20.txt')

    # sampler = BosonSampler(scheme, init_config)
    # sampler.export_ground_truth()
    # sampler.sample(batch_size=1, file_name='sample.txt')

    # un_sampler = UniformSampler(scheme, init_config)
    # un_sampler.sample(batch_size=10000, file_name='uniform_sample.txt')

    '''
    distinguishable_sampler = DistinguishableSampler(scheme, (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0))
    distinguishable_sampler.sample(batch_size=10000)
    '''
    time_end = time.time()

    print("\n--> The time of execution is :", (time_end - time_start) * 10 ** 3, "ms")


if __name__ == '__main__':
    main()
