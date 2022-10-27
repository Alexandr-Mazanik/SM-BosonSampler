import numpy as np


class BeamSplitter:
    def __init__(self, theta, phi, alpha, mods):
        self.mods = mods
        self._t = np.cos(theta)
        self._r = np.exp(phi * 1j) * np.sin(theta)

        self.unitary = self.calc_unitary()

    def calc_unitary(self):
        return np.array([[self._t, -np.conj(self._r)], [self._r, self._t]])


class BosonSampler:
    def __init__(self, m):
        self.number_of_modes = m
        self._beam_splitters = []

    def add_beam_splitter(self, theta, phi, alpha, mods):
        self._beam_splitters.append(BeamSplitter(theta, phi, alpha, mods))

    def get_some_info(self):
        print(self._beam_splitters)


def main():
    sampler = BosonSampler(6)

    sampler.add_beam_splitter(np.pi / 4, 0, 0, (1, 2))

    sampler.get_some_info()


if __name__ == '__main__':
    main()
