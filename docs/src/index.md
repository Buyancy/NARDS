# NARDS.jl

Below is the documentation for if you wanted to use the NARDS (NARDINI) algorithm in a program you wrote in Julia. However, you can also run the script from the command line to analyze sequences or FASTA files containing multiple sequences. 

You can access the help message for the program by running it with the `-h` flag: 
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

NARDS uses as many threads as possible to process largeer FASTA files. In general, doubling the number of threads will approximately half the time it takes to process the file. 


```@meta
CurrentModule = NARDS
```

```@docs
nards
```

```@docs
DEFAULT_CATEGORIES
```

```@docs
DEFAULT_CATEGORY_ORDER
```

```@docs
DEFAULT_CATEGORY_LABELS
```