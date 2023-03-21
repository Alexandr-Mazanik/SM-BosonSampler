import random
import numpy as np
import os


def generate_random_block_scheme(blocks_list, file_name='curr_scheme.txt'):
    with open(os.path.join('scheme', file_name), 'w') as f_scheme:
        for block in blocks_list:
            f_scheme.write("BLOCK" + '\n')
            f_scheme.write(str(block[0]) + ' ' + str(block[2]) + '\n')
            for _ in range(block[1]):
                set_of_modes = set()
                while len(set_of_modes) < 2:
                    set_of_modes.add(random.randint(1, block[0]))
                set_of_modes = list(set_of_modes)

                theta, phi_rho, phi_tau = [2 * np.pi * random.random() for _ in range(3)]

                f_scheme.write(str(set_of_modes[0]) + '\t' + str(set_of_modes[1]) + '\t' + str(theta)
                               + '\t' + str(phi_rho) + '\t' + str(phi_tau) + '\n')
            f_scheme.write("END" + '\n')


def to_simple_scheme(block_scheme_file_name='curr_scheme.txt', simple_scheme_file_name='curr_scheme_simple.txt'):
    block_number = 0
    bs_gates = []
    with open(os.path.join('scheme', block_scheme_file_name), 'r') as f_scheme:
        for f_line in f_scheme:
            if f_line.strip() == "BLOCK":
                block = []
                modes_num, repetitions_num = list(map(int, f_scheme.readline().split()))
                bs_line = f_scheme.readline()
                while bs_line.strip() != "END":
                    mode1, mode2 = list(map(int, bs_line.split('\t')[0:2]))
                    theta, phi_rho, phi_tau = list(map(float, bs_line.split('\t')[2:]))

                    block.append([mode1 + block_number, mode2 + block_number, theta, phi_rho, phi_tau])
                    bs_line = f_scheme.readline()

                for _ in range(repetitions_num):
                    for bs_gate in block:
                        bs_gates.append(bs_gate[:])
                        bs_gate[0] += 1
                        bs_gate[1] += 1
                    block_number += 1

    number_of_modes = block_number + modes_num - 1

    with open(os.path.join('scheme', simple_scheme_file_name), 'w') as f_scheme:
        f_scheme.write(str(number_of_modes) + '\n')
        for bs in bs_gates:
            f_scheme.write(str(bs[0]) + '\t' + str(bs[1]) + '\t' +
                           str(bs[2]) + '\t' + str(bs[3]) + '\t' + str(bs[4]) + '\n')

    return number_of_modes


def main():
    blocks_list = [[3, 9, 4]]
    generate_random_block_scheme(blocks_list)

    to_simple_scheme()


if __name__ == '__main__':
    main()
