from analyzer import *
from util import *
import re
import sys
import argparse
import glob


def parse_command():
    """
    Parses the command line arguments passed in with main.py which are used to
    drive behavior and/or values of the program.
    """
    if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Generate new branch file.')
        parser.add_argument('-o', type=str, nargs='?', default='output.csv',
                            help='Provide an output file. Default is results.csv')
        parser.add_argument('-r', type=int, nargs=2,
                            help='Provide a range of files to process. E.g. "1 10" Default searches current directory for all "filename_#.txt" files')
        parser.add_argument('-l', type=str, nargs='?',
                            help='Provide a list of input files to analyze. Default searches current directory for all "filename_#.txt" files')
        parser.add_argument('-v', type=int, nargs=1, default=0,
                            help='Turns verbose output on. Prints all generated tokens.')

        args = parser.parse_args(sys.argv[1:])
        main(args)

def print_summary(logs, outfiles, flag):
    """
    Prints an output summary based on the analysis that was performed.
    """
    width = 38
    print("\n\n********** ANALYSIS SUMMARY **********")

    print("\nResult of Analysis: " + flag)

    count = 1
    if logs:
        print("\nLogs that were analyzed: ")
        for log in logs:
            print(str(count) + ":  " + str(log))
            count += 1
    else:
        print("\nNo logs were found and analyzed based on input parameters.")

    count = 1
    print("\nOutput File(s): ")
    for outfile in outfiles:
        print(str(count) + ":  "  + str(outfile))
        count += 1

    print("\n", end="")
    for i in range(width):
        print("*", end="")
    print("\n")

def main(args):
    """
    Main function used to coordinate component pieces. Drives the program.
    """
    outfiles = [args.o, "matrix.txt"]
    infile = "input_files/linemap.csv"
    myanalyzer = Analyzer()
    text = ""
    logs = []
    config_values = read_config_file("input_files/config.txt")

    temp_logs = [log for log in glob.glob("*.txt") if re.match('[a-zA-Z _]+_[0-9]+.txt', str(log))]
    temp_logs.sort()

    if args.r:
        start = args.r[0]
        end = args.r[1]+1
        i = start
        for log in temp_logs:
            for i in range(start, end):
                if re.match('[a-zA-Z _-]+_'+str(i)+'.txt', str(log)):
                    logs.append(log)
    elif args.l:
        numbers = args.l.split()
        for log in temp_logs:
            for number in numbers:
                if re.match('[a-zA-Z _-]+_'+number+'.txt', str(log)):
                    logs.append(log)
    elif config_values[RANGE] != "":
        ranges = config_values[RANGE].split(" ")
        start = int(ranges[0])
        end = int(ranges[1])
        i = start
        for log in temp_logs:
            for i in range(start, end):
                if re.match('[a-zA-Z _-]+_'+str(i)+'.txt', str(log)):
                    logs.append(log)
    elif config_values[LIST] != "":
        numbers = config_values[LIST].split()
        for log in temp_logs:
            for number in numbers:
                if re.match('[a-zA-Z _-]+_'+number+'.txt', str(log)):
                    logs.append(log)
    else:
        logs = temp_logs

    if logs:
        for log in logs:
            text += myanalyzer.read(log)    # Build text string for tokenizer

    if logs and len(outfiles) > 0:
        myanalyzer.run(text, outfiles, infile, config_values, args.v)
        print_summary(logs, outfiles, 'SUCCESS')
    else:
        print_summary(logs, outfiles, 'FAILURE')

parse_command()
