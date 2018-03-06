# gmu-research
A code repository for Dr. Jie Xu's Power Systems Optimization Research.

## Log Analyzer

  The log analyzer is a tool used to tokenize, process, and analyze simulation log files. The results are printed to a .csv file. 
  
To run the analyzer:

`python3 main.py`
  
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

The regular expression is '[a-zA-Z_ -]+_[0-9]+.txt' which matches any combination of letters, numbers, spaces, hyphens, and underscores so long as the end of the filename consists of an underscore, followed by a number, and the extension .txt). 

However, there are two additional ways to specify file inputs, either as a range or a list.

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


## Config File

In addition to passing in command line arguments, it is also possible to specify the same information inside of a config file. The file is located in `/input_files/config.txt`. `numlines` can never be left blank. It is always required. All other fields can be left blank and will use their default values.

numlines=

range=

list=

delimiter=

verbose=

###Numlines:

Specify the number of lines (e.g. 233 lines): `numlines=233` 

###Range:

Specify the range of files you wish to process (e.g. 1 through 10): `range=1 10`

###List:

Specify a lost of files you wish to process (e.g. 1, 2, 4, 9): `list=1,2,4,9`

###Delimiter:

Specifies which type of matrix delimiter output you prefer (comma separated value or white space). 

For comma separated value (complete with csv headers): `delimiter=,`

For white space (no headers): `delimiter=sp`

If delimiter is left blank, it will output both matrix types.



