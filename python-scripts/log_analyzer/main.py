from analyzer import *
import re
import sys
import argparse
import glob


def parse_command():

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
    width = 38
    print("\n\n********** ANALYSIS SUMMARY **********")

    print("\nResult of Analysis: " + flag)

    count = 1
    if logs:
        print("\nLogs that were analyed: ")
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
    outfiles = [args.o]
    myanalyzer = Analyzer()
    text = ""
    logs = []

    temp_logs = [log for log in glob.glob("*.txt") if re.match('[a-zA-Z _]+_[0-9]+.txt', str(log))]
    temp_logs.sort()

    if args.r:
        start = args.r[0]
        end = args.r[1]+1
        i = start
        for log in temp_logs:
            for i in range(start, end):
                if re.match('[a-zA-Z _]+_'+str(i)+'.txt', str(log)):
                    logs.append(log)
    elif args.l:
        numbers = args.l.split()
        for log in temp_logs:
            for number in numbers:
                if re.match('[a-zA-Z _]+_'+number+'.txt', str(log)):
                    logs.append(log)
    else:
        logs = temp_logs

    if logs:
        for log in logs:
            text += myanalyzer.read(log)    # Build text string for tokenizer

    if logs and len(outfiles) == 1:
        myanalyzer.run(text, outfiles[0], args.v)
        print_summary(logs, outfiles, 'SUCCESS')
    else:
        print_summary(logs, outfiles, 'FAILURE')

parse_command()

#Elapsed time is 154.527286 seconds.
