from strawberryfields.ops import *
import strawberryfields as sf
import numpy as np
import math


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


def export_to_file(results, prob):
    with open('samples/sf_sample.txt', 'w') as f_out:
        for i, result in enumerate(results):
            f_out.write(str(result.samples[0]) + '\t' + str(np.round(prob[i], 4)) + '\n')
    print("--> Samples was successfully exported")


def boson_sampling(modes_num, injected_photons_num, batch_size):
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
    prob = []
    probs = eng.run(fock_states_prob)
    for i in range(batch_size):
        results.append(eng.run(boson_sampler))
        prob.append(probs.state.fock_prob(results[i].samples[0]))
        if i % 10 == 0:
            print("--> complete: ", i, "/", batch_size)
    export_to_file(results, prob)


def main():
    boson_sampling(modes_num=4, injected_photons_num=3, batch_size=100)


if __name__ == '__main__':
    main()
