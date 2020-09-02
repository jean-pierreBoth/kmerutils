# Some Kmer counting utilities

This package provides basic tools :

* kmer counting tools
  
* some statistics dumps :
  
  1. base distributions
  2. read length distributions
  3. number of bases between to occurences of A/T (or C/G) bases (return times).

* A quality server.
  This binary executable loads qualities from a Fastq file and runs as a server, and answer
  basic requests as returning a sequence of given number , or a block of a sequence or simply the base at a given pos  in a sequence given its rank. This enables storing qualities on a different machine.
  
The package has a Julia companion providing interactive access to dumped statistics or interactive inspection of sequences
of bases and qualities.

## Kmer Compression

The bases can be encoded on 2-bit or 4-bit alphabet. The Kmer can be compressed on these two alphabet (at present time only 2-bit compression is implemented).  
Kmer can be stored in 16-bit, 32-bit or 64-bit words thus providing compressed sequence up to 32 bases.

## Hashing

* Superminhash

An implementation of Superminhash :  
**A new minwise Hashing Algorithm for Jaccard Similarity Estimation**
Otmar Ertl 2017-2018 Cf [superminhash](https://arxiv.org/abs/1706.05698)

The hash values are computed by the sketch method or can be computed before entering SuperMinHash methods.
In this case (pre-hashed values) so the structure just computes permutation according to the paper.

It runs in one pass on data so it can be used in streaming.

* Minhash

A generic Minhash implementation based on BinaryHeap and HashMap

* Inversible Hash

The  module invhash.rs provides a rust implementation for inversible hash in 32bit and 64 bits version
See Thomas Wang's invertible integer hash functions.
at [Wang](https://gist.github.com/lh3/59882d6b96166dfc3d8d) for a snapshot.

See also:  
[Thomas Wang's 32 Bit Mix Function](http://www.cris.com/~Ttwang/tech/inthash.htm)
or [blink](https://chromium.googlesource.com/chromium/blink/+/master/Source/wtf/HashFunctions.h)
and [chris foster blog](http://c42f.github.io/2015/09/21/inverting-32-bit-wang-hash.html)

* Nthash

this is a recursive hashing described in:
    **"ntHash: recursive nucleotide hashing"**  
     Mohamadi Chu Birol BioInformatics 2016.
It is implemented on all our compressed kmer types.

## Counting Kmer

Kmer counting is multi-threaded and filters unique kmer in a cuckoo filter to spare memory.
Unique kmers are dumped in a separate file with the coordinates (sequence and position in sequence).
Multiple kmer are dumped in another file with their multiplicity.

* command to use:

## Some statistics on sequences

### Read length distributions
  
### Base distribution

We dump in a file the following information:

1. number of non acgt (first line of file)
2. a matrix (100, 4) giving for row i and column j in (1,2,3,4)the number of reads
where a base (a,c,g,t) corresponding to column j in this order occurs at percentage i.

This file can be reloaded by Julia package Genomics (cf BaseDistribution.jl)

### Return times

We estimate return times between occurrences of A/T bases or C/G bases. These statistics are reverse complement invariant.
The statistic is computed by default by dividing sequence into small (less than 255 bases length, and by default of length 190) blocks. So we cannot estimate return times between A/T (or C/G) greater than 255 but should be sufficient.

command to run : **./target/release/parsefastq -b 2 -f fastqfile ret -b'A'**

Output file name is $fastqfile.ret_times.bin.

File format: Format is documented in statutils.rs see **function dump_binary** related to struct ReturnTimesLR.

This file can be reloaded from Julia package Genomics (Cf ReturnTimesLRFile).

## Quality

Qualities are re-mapped to values between in [0..7] so that they need only 3 bits of storage and are
stored in a wavelet matrix.
The mapping is non uniform and maps the range  [0x25,0x37] to  [1,6].

## Quality Server

The server is launched on the server machine by the command:  
 **qualityloader -f filename [ -p portnum] [ --wavelet]**.

The server listens by default to port 4766, the option "--wavelet" asks for wavelet compression.

## Julia interface
