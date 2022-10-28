import numpy as np
from numpy.linalg import multi_dot


class BeamSplitter:
    def __init__(self, theta, phi, alpha):
        #TODO add phase shifting
        self._t = np.cos(theta)
        self._r = np.exp(phi * 1j) * np.sin(theta)

        self._mode1 = self._mode2 = None
        self.unitary = None

    def calc_unitary(self, number_of_modes, modes):
        self._mode1 = modes[0] - 1
        self._mode2 = modes[1] - 1
        self.unitary = np.identity(number_of_modes, dtype=complex)

        #TODO make it beauty
        self.unitary[self._mode1, self._mode1] = self._t
        self.unitary[self._mode1, self._mode2] = -np.conj(self._r)
        self.unitary[self._mode2, self._mode1] = self._r
        self.unitary[self._mode2, self._mode2] = self._t


class BosonSampler:
    def __init__(self, modes_num):
        self.number_of_modes = modes_num
        self._beam_splitters = []
        self.unitary = np.identity(self.number_of_modes)

    def add_beam_splitter(self, theta, phi, alpha, modes):
        bs = BeamSplitter(theta, phi, alpha)
        bs.calc_unitary(self.number_of_modes, modes)
        self._beam_splitters.append(bs)

    def calc_system_unitary(self):
        #TODO remove the excess identity matrix
        self.unitary = multi_dot([np.identity(self.number_of_modes)] +
                                 [bs.unitary for bs in np.flip(self._beam_splitters)])

    def get_some_info(self):
        #TODO change reading method
        print(self.unitary)


def main():
    #TODO make a test for unitary
    #TODO compare with SF
    #TODO faster
    sampler = BosonSampler(3)

    sampler.add_beam_splitter(np.pi / 4, 0, 0, (1, 2))
    sampler.add_beam_splitter(np.pi / 4, 0, 0, (2, 3))

    sampler.calc_system_unitary()

    sampler.get_some_info()


if __name__ == '__main__':
    main()
