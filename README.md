# Metagenome Graph Project

MetaGraph is a tool for scalable construction of annotated genome graphs and sequence-to-graph alignment.

The default index representations in MetaGraph are extremely scalable and support building graphs with trillions of nodes and millions of annotation labels.
At the same time, the provided workflows and their careful implementation, combined with low-level optimizations of the core data structures, enable exceptional query and alignment performance.

#### Main features:
* Large-scale indexing of sequences
* Python API for querying in the server mode
* Support for representing k-mer counts or expression values
* Sequence alignment against very large annotated graphs
* Scalable cleaning of very large de Bruijn graphs (to remove sequencing errors)
* Support for custom alphabets (e.g., {A,C,G,T,N} or amino acids)
* Algorithms for differential assembly

#### Design choices in MetaGraph:
* Use of succinct data structures and efficient representation schemes for extremely high scalability
* Algorithmic choices that work efficiently with succinct data structures (e.g., always prefer batched operations)
* Modular support of different graph and annotation representations
* Use of generic and extensible interfaces to support adding custom index representations / algorithms with little code overhead.

## Install

See [docs online](https://metagraph.ethz.ch/static/docs/index.html).


## Typical workflow
1. Build de Bruijn graph from Fasta files, FastQ files, or [KMC k-mer counters](https://github.com/refresh-bio/KMC/):\
`./metagraph build`
2. Annotate graph using the column compressed annotation:\
`./metagraph annotate`
3. Transform the built annotation to a different annotation scheme:\
`./metagraph transform_anno`
4. Query annotated graph\
`./metagraph query`

### Example
```
DATA="../tests/data/transcripts_1000.fa"

./metagraph build -k 12 -o transcripts_1000 $DATA

./metagraph annotate -i transcripts_1000.dbg --anno-filename -o transcripts_1000 $DATA

./metagraph query -i transcripts_1000.dbg -a transcripts_1000.column.annodbg $DATA

./metagraph stats -a transcripts_1000.column.annodbg transcripts_1000.dbg
```

### Print usage
`./metagraph`

### Build graph

* #### Simple build
```bash
./metagraph build -v --parallel 30 -k 20 --mem-cap-gb 10 \
                        -o <GRAPH_DIR>/graph <DATA_DIR>/*.fasta.gz \
2>&1 | tee <LOG_DIR>/log.txt
```

* #### Build with disk swap (use to limit the RAM usage)
```bash
./metagraph build -v --parallel 30 -k 20 --mem-cap-gb 10 --disk-swap <GRAPH_DIR> \
                        -o <GRAPH_DIR>/graph <DATA_DIR>/*.fasta.gz \
2>&1 | tee <LOG_DIR>/log.txt
```

#### Build from k-mers filtered with KMC
```bash
K=20
./KMC/kmc -ci5 -t4 -k$K -m5 -fm <FILE>.fasta.gz <FILE>.cutoff_5 ./KMC
./metagraph build -v -p 4 -k $K --mem-cap-gb 10 -o graph <FILE>.cutoff_5.kmc_pre
```

### Annotate graph
```bash
./metagraph annotate -v --anno-type row --fasta-anno \
                           -i primates.dbg \
                           -o primates \
                           ~/fasta_zurich/refs_chimpanzee_primates.fa
```

### Convert annotation to Multi-BRWT
1) Cluster columns
```bash
./metagraph transform_anno -v --linkage --greedy \
                           -o linkage.txt \
                           --subsample R \
                           -p NCORES \
                           primates.column.annodbg
```
Requires `N*R/8 + 6*N^2` bytes of RAM, where `N` is the number of columns and `R` is the number of rows subsampled.

2) Construct Multi-BRWT
```bash
./metagraph transform_anno -v -p NCORES --anno-type brwt \
                           --linkage-file linkage.txt \
                           -o primates \
                           --parallel-nodes V \
                           -p NCORES \
                           primates.column.annodbg
```
Requires `M*V/8 + Size(BRWT)` bytes of RAM, where `M` is the number of rows in the annotation and `V` is the number of nodes merged concurrently.

### Query graph
```bash
./metagraph query -v -i <GRAPH_DIR>/graph.dbg \
                        -a <GRAPH_DIR>/annotation.column.annodbg \
                        --discovery-fraction 0.8 --labels-delimiter ", " \
                        query_seq.fa
```

### Align to graph
```bash
./metagraph align -v -i <GRAPH_DIR>/graph.dbg query_seq.fa
```

### Assemble sequences
```bash
./metagraph assemble -v <GRAPH_DIR>/graph.dbg \
                        -o assembled.fa \
                        --unitigs
```

### Assemble differential sequences
```bash
./metagraph assemble -v <GRAPH_DIR>/graph.dbg \
                        --unitigs \
                        -a <GRAPH_DIR>/annotation.column.annodbg \
                        --label-mask-in LABEL_1 \
                        --label-mask-in LABEL_2 \
                        --label-mask-out LABEL_3 \
                        -o diff_assembled.fa
```

### Get stats
Stats for graph
```bash
./metagraph stats graph.dbg
```
Stats for annotation
```bash
./metagraph stats -a annotation.column.annodbg
```
Stats for both
```bash
./metagraph stats -a annotation.column.annodbg graph.dbg
```

## Developer Notes

### Makefile

The `Makefile` in the top level source directory can be used to build and test `metagraph` more conveniently. The following
arguments are supported:
* `env`: environment in which to compile/run (`""`: on the host, `docker`: in a docker container)
* `alphabet`: compile metagraph for a certain alphabet (e.g. `DNA` or `Protein`, default `DNA`)
* `additional_cmake_args`: additional arguments to pass to cmake.

Examples:

```
# compiles metagraph in a docker container for the `DNA` alphabet
make build-metagraph env=docker alphabet=DNA
```


## License
Metagraph is distributed under the GPLv3 License (see LICENSE).
Please find further information in the AUTHORS and COPYRIGHTS files.
