# gmu-research
A code repository for Dr. Jie Xu's Power Systems Optimization Research.

## Log Analyzer

  The log analyzer is a tool used to tokenize, process, and analyze simulation log files. The results are printed to a .csv file. 
  
### Command Line Arguments:

1. **Verbose**: 

Verbose mode is turned on using the `-v` argument followed by an integer which determines the verbosity level. Currently, 0 turns verbose output off, while 1 simply prints each token. There are plans to add additional verbosity levels, primarily used to check processing against output and for debugging.

####Example:
`python3 main.py -v 1`
