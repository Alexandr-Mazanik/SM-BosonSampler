from strawberryfields.ops import *
import strawberryfields as sf
import numpy as np


def read_u(modes_num):
    with open('scheme/scheme_unitary.txt', 'r') as f_unitary:
        u = np.empty((modes_num, modes_num), dtype=complex)
        for line_num, f_line in enumerate(f_unitary):
            f_line = f_line.strip()
            string = list(map(float, f_line.split('\t')))
            for i in range(0, 2 * modes_num, 2):
                item = complex(string[i], string[i + 1])
                u[line_num][int(i/2)] = item
        print(u)
    print("--> Scheme was successfully uploaded")

    return u


def boson_sampling(modes_num, injected_photons_num, batch_size):
    boson_sampler = sf.Program(modes_num)

    with boson_sampler.context as q:
        for i in range(injected_photons_num):
            Fock(1) | q[i]

        u = read_u(modes_num)
        Interferometer(u) | q

        MeasureFock() | q

    eng = sf.Engine(backend="fock", backend_options={"cutoff_dim": 3})

    results = []
    j = 0
    for i in range(batch_size):
        results.append(eng.run(boson_sampler))
        if results[i].samples[0][1] == 2:
            j += 1
    print("prob = ", j / batch_size)


def main():
    boson_sampling(modes_num=2, injected_photons_num=2, batch_size=1000)


if __name__ == '__main__':
    main()
