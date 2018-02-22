# gmu-research
A code repository for Dr. Jie Xu's Power Systems Optimization Research.

## Log Analyzer

  The log analyzer is a tool used to tokenize, process, and analyze simulation log files. The results are printed to a .csv file. 
  
### Command Line Arguments:

1. **Verbose**: 

Verbose mode is turned on using the `-v` argument followed by an integer which determines the verbosity level. Currently, 0 turns verbose output off, while 1 simply prints each token. There are plans to add additional verbosity levels, primarily used to check processing against output and for debugging.

##### Example:

`python3 main.py -v 1`

2. **Output File**: 

The default output file is *output.csv*, however, you can specify an output file by using the `-o` argument followed by the file name.

#### Example:

`python3 main.py -o newoutputfile.csv`

3. **Input File**: 

The default input method is to scan the directory of the program and read all files which match the format *filename_#.txt (e.g. cas_results_#.txt)*. 

The regular expression is '[a-zA-Z_ -]+_[0-9]+.txt' which matches any combination of letters, numbers, spaces hyphens, and underscores so long as the last two characters are an underscored followed by a number and the extension .txt). However, there are two additional ways to specify file inputs, either as a range or a list.

   - **Range**

In order to specify a range simply use the `-r` argument followed by two integers specifying the range, inclusively. 

#### Example:

`python3 main.py -r 1 10`

Or for a single file

`python3 main.py -r 10 10`

   - **List**

In order to specify a list simply use the `-l` argument followed by the list of file numbers, separated by a space and wrapped in single or double quotation marks.

#### Example:

`python3 main.py -l "1 3 5 129"`
