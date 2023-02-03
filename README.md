# NARDS
NARDS stands for **N**onrandom **A**rrangements of **R**esidues in **D**isorderd **S**huffles. This project is pretty simple, it is a Julia implementation of the NARDINI algorithm. It gives users access to more features 
and runst efficiently using as many threads as there are available.  

## Usage 
It is possible to run the program on the command line. You can obtain a pretty good help page by just running the program by using the ```-h``` flag when you call the program. 
```
% julia NARDS.jl -h
usage: NARDINI.jl [-g WINDOW_SIZE] [-o OUTOUT_PATH]
                  [-n NUM_NULL_MODELS] [--title TITLE] [-h] input

positional arguments:
  input                 The sequence (or fasta) to be analyzed

optional arguments:
  -g, --window_size WINDOW_SIZE
                        The size of the sliding window to use (type:
                        Int64, default: 5)
  -o, --outout_path OUTOUT_PATH
                        The path to the file where the output will be
                        stored (default: ".")
  -n, --num_null_models NUM_NULL_MODELS
                        The number of sequence scrambles to use in the
                        null model (type: Int64, default: 100000)
  --title TITLE         The title of the output plot (default: "")
  -h, --help            show this help message and exit
```
Never the less, here is a more detailed description of each parameter: 
| Parameter                             | Description |
|---------------------------------------|-------------|
| ```-g, --window_size WINDOW_SIZE```         | The size of the sliding window to use for the asymmetry calculation. |
| ```-o, --outout_path OUTOUT_PATH```         | The path to the directory where the output of the algorithm (CSVs and PNGs) will be stored. |
| ```-n, --num_null_models NUM_NULL_MODELS``` | The number of sequences to use in the null model for the parameters. |
| ```--title TITLE```                         | The title of the sequence. Not used if a FASTA is used as input. |
| ```-h, --help```                            | Prints the help message. |

The ```input``` of the  program can be either a sequence or a path to a FASTA file. If it is just a sequence, then the program will run the NARDINI algorithm on that sequence and produce a heatmap result as well as a CSV containing all of the z-scores from the analysis. If the ```input``` is a FASTA file, then the program will run each sequence through the NARDINI algorithm, and the output will be saved in the ```output``` directory. Then, it will produce a figure that gives an overall heatmap comparing all of the sequences present in the FASTA file. 

Note that for large FASTA files, using more threads will significantly speed up the computation time. The number of threads can be changed when launching the program using the Julia ```-t``` flag. For example, to use 8 threads, you would use:
```
% julia -t 8 NARDS.jl ...
```

## Refrences 
Original NARDINI paper: 
```
@article{COHAN2022167373,
    title = {Uncovering Non-random Binary Patterns Within Sequences of Intrinsically Disordered Proteins},
    journal = {Journal of Molecular Biology},
    volume = {434},
    number = {2},
    pages = {167373},
    year = {2022},
    issn = {0022-2836},
    doi = {https://doi.org/10.1016/j.jmb.2021.167373},
    url = {https://www.sciencedirect.com/science/article/pii/S0022283621006100},
    author = {Megan C. Cohan and Min Kyung Shinn and Jared M. Lalmansingh and Rohit V. Pappu},
    keywords = {intrinsically disordered proteins/regions, binary patterns, NARDINI, CIDER}
}
```