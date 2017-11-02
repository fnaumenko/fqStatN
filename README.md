# fqStatN 
Calculates the statistics of occurrence of ambiguous character 'N', 
and patterns of reads including 'N' in the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file.<br>
These statistics helps to better evaluate the quality of the sequencer output.

The program runs on the command line under Linux and Windows.

## Installation
### Executable file

**Linux**<br>
Go to the desire directory and type commands:<br>
```wget -O fqStatN.gz https://github.com/fnaumenko/fqStatN/releases/download/1.0/fqStatN-Linux-x64.gz```<br>
```gzip -d fqStatN.gz```<br>
```chmod +x fqStatN```

**Windows**<br>
Download archive from [here](https://github.com/fnaumenko/fqStatN/releases/download/1.0/fqStatN-Windows-x64.zip) 
and unzip by any archiver, for instance [WinRar](https://www.win-rar.com/download.html?&L=0).

### Compiling in Linux
Required libraries:<br>
g++<br>
zlib (optionally)

Go to the desired directory and type commands:<br>
```wget -O fqStatN.zip https://github.com/fnaumenko/fqStatN/archive/1.0.zip```<br>
```unzip fqStatN.zip```<br>
```cd fqStatN-1.0```<br>
```make```

If **zlib** is not installed on your system, a linker message will be displayed.<br>
In that case you can compile the program without the ability to work with .gz files: 
open *makefile* in any text editor, uncomment last macro in the second line, comment third line, save *makefile*, and try ```make``` again.<br>
To be sure about **zlib** on your system, type ```whereis zlib```.

## Usage
```
fqStatN [options] file.fq
```

### Help
```
Options:
  -o|--out      duplicate standard output to fqStatN_out.txt file
  -t|--time     print run time
  -h|--help     print usage information and exit
```

## Details

### Input
Primary DNA sequence in [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format.<br>
Compressed file in gzip format (.gz) is acceptable.

### Output
The program displays the frequency of occurrence of 'N' in the position in the read, 
as well as the frequency of the pattern of the reads containing 'N'.<br>
Here is a sample of output:
```
'N' POSITION STATISTICS
pos     count   % of total 'N'
------------------------------
 0      29588   71.1%
 3      4314    10.4%
20      69      0.166%
21      49      0.118%
. . .
47      256     0.615%
48      33      0.0793%
49      1904    4.58%

READ PATTERN STATISTICS
position  10        20        30        40        50    count
01234567890123456789012345678901234567890123456789
----------------------------------------------------------------------
N.................................................    29566     0.172%
...N..............................................     4313     0.0251%
..............................NN..................       83     <0.001%
.............................NNNN.NNNN.........N.N       10     <0.001%
. . .
...............................N..NN..............        3     <0.001%
...............................N..N...............        2     <0.001%
.........................N.....N..NN.N...........N        1     <0.001%

'N' relative to the total number of nucleotides: 0.00483%
Reads that include 'N' relative to the total number of reads: 0.211%
```

##
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
