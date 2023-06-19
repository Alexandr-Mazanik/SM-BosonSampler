from strawberryfields.ops import *
from thewalrus import perm
import strawberryfields as sf
import numpy as np
import math
import os


def calc_prob_uniform(modes_num, injected_photons_num):
    return 1 / math.comb(injected_photons_num + modes_num - 1, modes_num - 1)


def read_u(modes_num):
    with open('scheme/scheme_unitary.txt', 'r') as f_unitary:
        u = np.empty((modes_num, modes_num), dtype=complex)
        for line_num, f_line in enumerate(f_unitary):
            f_line = f_line.strip()
            string = list(map(float, f_line.split('\t')))
            for i in range(0, 2 * modes_num, 2):
                item = complex(string[i], string[i + 1])
                u[line_num][int(i/2)] = item
    print("--> Scheme was successfully uploaded")

    return u


def export_to_file(results, prob_perm, file_name):
    with open(os.path.join('samples', file_name), 'w') as f_out:
        for i, result in enumerate(results):
            f_out.write(str(result.samples[0]) + '\t' + str(np.round(prob_perm[i], 4)) + '\n')
    print("--> Samples was successfully exported")


def find_submatrix(sample, u, ph_num):
    column_indexes = [i for i in range(ph_num)]
    row_indexes = []
    for i in range(len(sample)):
        if sample[i] > 0:
            for j in range(sample[i]):
                row_indexes.append(i)

    return u[:, column_indexes][row_indexes]


def prob_using_perm(ph_number, sample, u):
    norm = 1
    for num in sample:
        norm *= np.math.factorial(num)

    return abs(perm(find_submatrix(sample, u, ph_number), method="ryser"))**2 / norm


def boson_sampling(modes_num, injected_photons_num, batch_size, file_name):
    boson_sampler = sf.Program(modes_num)
    fock_states_prob = sf.Program(modes_num)

    u = read_u(modes_num)

    with boson_sampler.context as q:
        for i in range(injected_photons_num):
            Fock(1) | q[i]

        Interferometer(u) | q

        MeasureFock() | q

    with fock_states_prob.context as q:
        for i in range(injected_photons_num):
            Fock(1) | q[i]

        Interferometer(u) | q

    eng = sf.Engine(backend="fock", backend_options={"cutoff_dim": injected_photons_num + 1})

    print("--> Probability for uniform distribution: ", calc_prob_uniform(modes_num, injected_photons_num))

    results = []
    prob_perm = []
    for i in range(batch_size):
        result = eng.run(boson_sampler)
        sample = result.samples[0]
        results.append(result)

        prob_perm.append(prob_using_perm(injected_photons_num, sample, u))

        if i % 10 == 0:
            print("--> complete: ", i, "/", batch_size)
    export_to_file(results, prob_perm, file_name=file_name)


def main():
    boson_sampling(modes_num=4, injected_photons_num=3, batch_size=700, file_name='sf_sample.txt')


if __name__ == '__main__':
    main()
