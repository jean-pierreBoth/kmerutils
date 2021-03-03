# Some Kmer counting utilities

This package (currently in development) provides the following tools :

* kmer counting tools.

* Sketching of sequences with up to date sensitive hashing related to the **Probability Jaccard index** see the *jaccardweight* module.  

* A quality server.
  The binary executable *qualityloader* loads qualities from a Fastq file and runs as a server, answering
  basic requests as returning a quality sequence of given seuqnce number , or a block of a sequence or simply the base at a given pos  in a sequence given its rank. This enables storing qualities on a different machine.

* some basic statistics dumps  such as base distributions, read length distributions.

The package has a Julia companion providing interactive access to dumped statistics or interactive inspection of sequences
of bases and qualities.

## Kmer Compression and Counting

The bases are presently encoded on 2 bits.  
Kmer can be stored 32-bit or 64-bit words thus providing compressed representation up to 32 bases with the 2-bit alphabet.  
Kmer and compressed Kmer are represented respectively by trait *KmerT* and *CompressedKmerT*.
A kmer is identified with its reverse complement in the counting methods.  

Kmer counting is multi-threaded and filters unique kmer in a cuckoo filter to spare memory.
Unique kmers are dumped in a separate file with the coordinates (sequence and position in sequence).
Multiple kmers, stored in a Bloom filter, are dumped in another file with their multiplicity. See module *kmercount*

## Hashing and Sketching of data

Similarity between sequences can be estimated by counting common Kmers (i.e by estimating a Jaccard index) between sequences.

* Estimators of the standard Jaccard index (without taking into account Kmer multiplicity) is provided by the Superminhash algorithm.  

* A probability Jaccard index taking into account Kmer multiplicity is also provided with the Probminhash family algorithm.
Probminhash and superminhash are provided by the crate probminhash and are interfaced with kmer generation in module **jaccardweight.rs**.
(Cf [probminhash](https://github.com/jean-pierreBoth/probminhash)).

* The probminhash algorithm is used to provide a complete sketching of a datafile where each sequence has its signature
dumped in a file. This file can be reprocessed to examine neighborhood of a read in term of the Probability Jaccard index. see module *jaccarweight* or *seqblocksketch*.  
For example it takes 141s on a 4-i7 (2.7Ghz) core laptop, to read , generate 8 base kmers and sketch 746333 long reads from a 4.38 Gbases ONT fastq file (Cf [FAB49164_rel3](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-genome/rel_3_4.md)), asking for 200 sketches by read.

* The signatures obtained can be sent in an Ann to study read proximity according to the Jaccard Probability metric.
  See executable *datasketcher*in this crate and the crate *hnsw_rs*

Some others standard tools such :

* Nthash : This is a recursive hashing described in: **"ntHash: recursive nucleotide hashing"**  
     Mohamadi Chu Birol BioInformatics 2016.
It is implemented on all our compressed kmer types.

* Minhash : A generic Minhash implementation based on BinaryHeap and HashMap

## Some basic statistics on sequences

1. Read length distributions.  
    A file giving the number of reads in function of length.  

2. Base distributions.  
    a matrix (100, 4) giving for row i and column j in (1,2,3,4) the number of reads
    where a base (a,c,g,t) corresponding to column j in this order occurs at percentage i.

This file can be reloaded by Julia package Genomics (cf BaseDistribution.jl)

## Quality

Qualities are re-mapped to values between in [0..7] so that they need only 3 bits of storage and are
stored in a wavelet matrix.
The mapping is non uniform and maps the range  [0x25,0x37] to  [1,6].  
The quality part of data are stored in a process serving quality requests described below:

### Quality Server

The server is launched on the server machine by the command:  
 **qualityloader -f filename [ -p portnum] [ --wavelet]**.

The server listens by default to port 4766, the option "--wavelet" asks for wavelet compression.
