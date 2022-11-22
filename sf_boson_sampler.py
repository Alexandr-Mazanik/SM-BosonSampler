import numpy as np
from strawberryfields.ops import *
import strawberryfields as sf


def read_u(modes_num):
    

def boson_sampling(modes_num, injected_photons_num, batch_size):
    boson_sampler = sf.Program(modes_num)

    with boson_sampler.context as q:
        for i in range(injected_photons_num):
            Fock(1) | q[i]

        u = read_u(modes_num)
        Interferometer(u) | q

        MeasureFock() | q

    results = []
    j = 0
    eng = sf.Engine(backend="fock", backend_options={"cutoff_dim": 3})
    for i in range(batch_size):
        results.append(eng.run(boson_sampler))
        if results[i].samples[0][1] == 2:
            j += 1
    print("prob = ", j / batch_size)


def main():
    boson_sampling(modes_num=2, injected_photons_num=2, batch_size=10)


if __name__ == '__main__':
    main()
