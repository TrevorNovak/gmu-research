from analyzer import *

def main():
    myanalyzer = Analyzer()
    text = myanalyzer.read("cas_results_5.txt")
    myanalyzer.run(text)

main()
