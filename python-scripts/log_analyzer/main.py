from analyzer import *
import re
import sys
import argparse
import glob


def parse_command():
    """ NOT CURRENTLY IMPLEMENTED """

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

def main():
    filename = "results.csv"
    myanalyzer = Analyzer()
    text = ""

    logs = glob.glob("*.txt")               # Grab list of .txt files in directory
    logs.sort()
    #print(logs)
    for log in logs:                        # Check for output logs
        if re.match("[a-zA-Z _]+_\d.txt", log):
            text += myanalyzer.read(log)    # Build text string for tokenizer

    myanalyzer.run(text)
    print("The operation was successful. Please check the output file " + filename + " for analysis.")

main()
