import random
import numpy as np
import os


def generate_random_scheme(number_of_modes, bs_gate_number, file_name='curr_scheme.txt'):
    with open(os.path.join('scheme', file_name), 'w') as f_scheme:
        f_scheme.write(str(number_of_modes) + '\n')
        for _ in range(bs_gate_number):
            set_of_modes = set()
            while len(set_of_modes) < 2:
                set_of_modes.add(random.randint(1, number_of_modes))
            set_of_modes = list(set_of_modes)

            theta, phi_rho, phi_tau = [2 * np.pi * random.random() for _ in range(3)]

            f_scheme.write(str(set_of_modes[0]) + '\t' + str(set_of_modes[1]) + '\t' + str(theta)
                           + '\t' + str(phi_rho) + '\t' + str(phi_tau) + '\n')


def main():
    number_of_modes = 7
    bs_gate_number = 49

    generate_random_scheme(number_of_modes, bs_gate_number)


if __name__ == '__main__':
    main()
