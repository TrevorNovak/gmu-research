import random
import sys
import argparse

###############################################################################
# Simple script used to setup simulation scenarios either within a text file  #
# or through the command line.                                                #
#                                                                             #
# Author: Trevor Novak                                                        #
###############################################################################

def parse_command():
    if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Generate new branch file.')
        parser.add_argument('-b', type=int, nargs='?', default=31,
                            help='set number of branches for this simulation. [Default is 31.]')
        parser.add_argument('-l', type=float, nargs='?', default=.05,
                            help='set the probabilty for random selection of hardened lines for this simulation. [Default is randint(0, 1)]')
        parser.add_argument('-w', type=float, nargs='?', default=.05,
                            help='set the probability for random selection of line located in weather zone for this simulation. [Default is randint(0, 1)]')
        parser.add_argument('-f', type=float, nargs='?', default=.05,
                            help='set the probability for random selection of failure rates for each branch for this simulation. [Default is uniform(0, 1)]')
        parser.add_argument('-fr', type=float, nargs='+', default=[0, 1],
                            help='set the range for random number generation of failure rates. [Default is range(0, 1)')
        parser.add_argument('-o', nargs='?', default="data.txt",
                            help='set the output file. [Default is data.txt]')
        parser.add_argument('-i', nargs='?',
                            help='set the input file. [Default is input.txt]')
        # parser.add_argument('--gen', nargs='?',
        #                     help='generate a new simulation file')

        args = parser.parse_args(sys.argv[1:])
        if(args.i):
            read_file(args)
        else:
            write_file(args)

def read_file(args):
    indata = args.i
    outdata = args.o
    zones = []
    max_harden = 0
    failure_rate_prob = 0
    num_branches = args.b

    with open(indata, "r") as f:
        zones = list(map(int, f.readline().split(',')))
        max_harden = int(f.readline())
        failure_rate_prob = float(f.readline())

    harden_count = 0
    zone_count = 0
    zone_pointer = 0
    with open(outdata, "w") as f:
        for i in range(num_branches):
            if harden_count < max_harden:
                line_harden_val = random.randint(0, 1)
                if(line_harden_val == 1):
                    harden_count += 1
            else:
                line_harden_val = 0

            if zone_pointer < len(zones) and zone_count == zones[zone_pointer]:
                weather_zone_val = 1
                zone_pointer += 1
            else:
                weather_zone_val = 0

            zone_count += 1
            f.write(str(weather_zone_val) + ' ' + str(line_harden_val) + ' ' + "{0:.2f} ".format(failure_rate_prob) + '\n')
            print(str(weather_zone_val) + ' ' + str(line_harden_val) + ' ' + "{0:.2f} ".format(failure_rate_prob))
        print("\nFile " + outdata + " has been generated successfully!\n")
        print("Failure Rate: " + str(failure_rate_prob))
        print("Number of Branches: " + str(num_branches))
        print("Weather Zones: " + str(zones))
        print("Max Hardened Lines: " + str(max_harden))


def write_file(args):
    fn = args.o
    num_branches = args.b
    line_harden_prob = args.l
    weather_zone_prob = args.w
    failure_rate_prob = 0
    fmin = args.fr[0]
    fmax = args.fr[1]

    with open(fn, "w+") as f:
        for i in range(num_branches):
            line_harden_val = random.randint(0, 1)
            weather_zone_val = random.randint(0, 1)
            failure_rate_prob = random.uniform(fmin, fmax)
            f.write(str(weather_zone_val) + ' ' + str(line_harden_val) + ' ' + "{0:.2f} ".format(failure_rate_prob) + '\n')
            print(str(weather_zone_val) + ' ' + str(line_harden_val) + ' ' + "{0:.2f} ".format(failure_rate_prob))

    print("\nFile " + fn + " has been generated successfully!\n")
    print("Number of Branches: " + str(num_branches) + "\n")

def calculate():
    answer = 0
    for i in range(100):
        answer += i
        print(answer)

parse_command()
