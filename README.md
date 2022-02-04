# Improved (Faster) Recursive Graph Bisection

This repo contains some extensions of the code corresponding to the SIGIR 2021 short paper 
*Faster Index Reordering with Bipartite Graph Partitioning* by Joel Mackenzie,
Matthias Petri, and Alistair Moffat. This repo is based directly off of the
[faster-graph-bisection](http://github.com/mpetri/faster-graph-bisection) toolkit
which was proposed alongside that SIGIR paper.

## Extensions
The extensions will be outlined in more detail soon. If you are interested, have a look
at the commits and the code in `src/rgb.rs`. We have added a purely iterative implementation
which adds a synchronization step at each level of the "recursion" tree, reducing contention
and resulting in faster processing.

## Citation Information
If you use this code in your own work or research, please consider citing
our prior work:
```
@inproceedings{mpm21-sigir,
 title = {Faster Index Reordering with Bipartite Graph Partitioning},
 author = {J. Mackenzie and M. Petri and A. Moffat},
 booktitle = {Proc. SIGIR},
 pages = {1910--1914},
 year = {2021},
}
```
The paper can be found at the following DOI: https://doi.org/10.1145/3404835.3462991

## Acknowledgements
This work was built on previous work from Dhulipala et. al:
[Compressing Graphs and Indexes with Recursive Graph Bisection](http://www.kdd.org/kdd2016/papers/files/rpp0883-dhulipalaAemb.pdf), 
[ACM Proceedings](https://dl.acm.org/citation.cfm?id=2939862).

We also used the reproducibility study from Mackenzie et. al:
[Compressing Inverted Indexes with Recursive Graph Bisection: A Reproducibility Study](http://engineering.nyu.edu/~suel/papers/bp-ecir19.pdf),
[Springer Proceedings](https://link.springer.com/chapter/10.1007/978-3-030-15712-8_22).

Our codebase is based on the implementation found in the [PISA](https://github.com/pisa-engine/pisa/) search engine, which
corresponds to the reproducibility study discussed above. The codebase works with the
[Common Index File Format](https://github.com/osirrc/ciff), an open-source index exchange format for information
retrieval experimentation.


## Building the code
You can build the code using Cargo:
```
cargo build --release
```

However, if you follow the command above, running the code will give an error:
```
./target/release/create-rgb
03:09:14 [INFO] Error: A gain function needs to be passed at compile time via the environment variable `GAIN` -- Please recompile...
```

The explanation is that, since we experimented with three different gain functions, the desired gain function must be passed in
at compile time via an environment variable. The valid options are `default`, `approx_1`, or `approx_2`. So, recompile as such:
```
GAIN=approx_1 cargo build --release
```

## Running the code
[Usage is the same as in the previous tool](https://github.com/mpetri/faster-graph-bisection#running-the-code).

## Settings and Configuration

A full suite of settings can be found using the `--help` flag, and are listed as follows:
```
create-rgb 0.1.0
Reorders documents using recursive graph bisection and ciff files.

USAGE:
    create-rgb [FLAGS] [OPTIONS] --input <input>

FLAGS:
    -h, --help         Prints help information
    -l, --loggap       Show loggap cost
        --sort-leaf    Sort leaf by identifier
    -V, --version      Prints version information

OPTIONS:
    -c, --cutoff-frequency <cutoff-frequency>
            Maximum length to consider in percentage of documents in the index [default: 0.1]

    -i, --input <input>                          Input file ciff file
        --input-fidx <input-fidx>                Read forward index
        --max-depth <max-depth>                  Maximum depth [default: 100]
    -m, --min-len <min-len>                      Minimum number of occurrences to consider [default: 4096]
    -o, --output-ciff <output-ciff>              Output ciff file
        --output-fidx <output-fidx>              Output forward index
        --output-mapping <output-mapping>        Dump the document map
    -p, --parallel-switch <parallel-switch>
            Depth where we switch from parallel processing to sequential processing [default: 10]

    -r, --recursion-stop <recursion-stop>        Min partition size [default: 16]
    -s, --swap-iterations <swap-iterations>      Swap iterations [default: 20]
```

For example, you can save a forward index using the `--output-fidx` command, and can read a saved forward index
with the `--input-fidx` flag. If you only wish to dump the reordered document map, use the `--output-mapping`
flag. 

## Support
Feel free to raise issues here, we'll do our best to assist.
