import numpy as np
from numpy.linalg import multi_dot


class BeamSplitter:
    def __init__(self, theta, phi, alpha, modes, num_of_modes):
        self.modes = modes
        self._number_of_modes = num_of_modes

        self._t = np.cos(theta)
        self._r = np.exp(phi * 1j) * np.sin(theta)

        self.unitary = self.calc_unitary()

    def calc_unitary(self):
        unitary = np.identity(self._number_of_modes, dtype=complex)

        unitary[self.modes[0]-1, self.modes[0]-1] = self._t
        unitary[self.modes[0]-1, self.modes[1]-1] = -np.conj(self._r)
        unitary[self.modes[1]-1, self.modes[0]-1] = self._r
        unitary[self.modes[1]-1, self.modes[1]-1] = self._t

        return unitary


class BosonSampler:
    def __init__(self, modes_num):
        self.number_of_modes = modes_num
        self._beam_splitters = []
        self.unitary = np.identity(self.number_of_modes)

    def add_beam_splitter(self, theta, phi, alpha, modes):
        self._beam_splitters.append(BeamSplitter(theta, phi, alpha, modes, self.number_of_modes))

    def calc_unitary(self):
        self._beam_splitters = np.flip(self._beam_splitters)
        self.unitary = multi_dot([np.identity(self.number_of_modes)] +
                                 [bs.unitary for bs in self._beam_splitters])

    def get_some_info(self):
        print(self.unitary)


def main():
    sampler = BosonSampler(7)

    for i in range(6):
        sampler.add_beam_splitter(np.pi / 4, 0, 0, (i, i+1))

    sampler.calc_unitary()

    sampler.get_some_info()


if __name__ == '__main__':
    main()
