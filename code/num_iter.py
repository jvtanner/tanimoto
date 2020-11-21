import random
import subprocess
import statistics as stat
import tanimoto


def main():

    # Select 100 random seeds
    seeds = []
    for _ in range(100):
        seeds.append(random.randint(0, 999))

    iterations = ['100', '500', '1000']

    results = []
    for iteration in iterations:
        out = open(f'{iteration}_iteration_output.txt', 'w')
        for seed in seeds:
            completed_process = subprocess.run(["python", "code/pvalue.py", "-n " + iteration,
                                      "-r " + str(seed), 'data/drugs.csv', 'data/targets.csv',
                                      'P54577', 'Q7RTX0'], capture_output=True, text=True)
            p_value = completed_process.stdout
            out.write(p_value)
            results.append(float(p_value.strip()))

        stdev = stat.stdev(results)
        mean = stat.mean(results)
        out.write('Standard Deviation: {:.6f}\n'.format(stdev))
        out.write('Mean: {:.6f}'.format(mean))
        tanimoto.histo(results, iteration)
        print(f'Max p-value for {iteration}: {max(results)}\n')
        out.close()


if __name__ == '__main__':
    main()
