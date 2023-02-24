def main():
    data_file_name = 'sf_sample.txt'

    all_samples = []
    occurrences = []
    theory_prob = []

    with open(data_file_name, 'r') as f_data:
        for f_sample in f_data:
            f_sample = f_sample.strip()
            sample, probability = f_sample.split('\t')
            if sample in all_samples:
                occurrences[all_samples.index(sample)] += 1
            else:
                all_samples.append(sample)
                theory_prob.append(probability)
                occurrences.append(1)

    batch_size = sum(occurrences)
    experimental_prob = [i / batch_size for i in occurrences]

    with open('analysis/sample_analysis.txt', 'w') as f_sampler_analysis:
        for i, sample in enumerate(all_samples):
            f_sampler_analysis.write(sample + '\t' + theory_prob[i] + '\t' + str(experimental_prob[i]) + '\n')


if __name__ == '__main__':
    main()
