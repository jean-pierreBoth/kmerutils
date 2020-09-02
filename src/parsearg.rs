/// this module constains structures
/// describing arguments to different functions : parsing , compression , return times sampling

// All these structures implement clone as it will only rarely be used.


#[derive(Clone)]
pub struct ReturnTimesArgs {
    // true if return times must be computed
    pub to_do: bool,
    /// searched base,
    pub searched_base: u8,
    // size consecutive base regions to analyze
    pub window_size: u8,
    /// if non null gives the name of file to reload ret times from
    file:Option<String>,
}


impl Default for ReturnTimesArgs {
    fn default() -> ReturnTimesArgs {
        ReturnTimesArgs{to_do: false,
                        searched_base: b'A',
                        window_size: 190,
                        file: None}
    } // end of function default
}

//============================================================================


/// This enum stores what is asked for Kmer processing
#[derive(Clone)]
pub enum KmerProcessing {
    /// Nothing is asked
    None,
    /// We must count kmers (heavy)
    Counting,
    /// We just keep track of unique kmers
    Unicity,
}


/// This structure stores all parameters related to Kmers option
#[derive(Clone)]
pub struct KmerArgs {
    /// what is to be done
    pub kmer_task:KmerProcessing,
    /// filename to dump result in
    file:Option<String>,
    /// number of threads to use
    pub nb_threads: usize,
    /// kmer size
    pub kmer_size: usize,
    /// counter size (in bits) 8 or 16 for very high coverage (ONT)
    pub counter_size: usize,   
}



impl Default for KmerArgs {
    fn default() -> KmerArgs {
        KmerArgs {
            kmer_task: KmerProcessing::None,
            file: None,
            nb_threads: 2,   // every body have 2 cores , no?
            kmer_size: 16,
            counter_size: 8,
        }
    } // end of function default
}


//============================================================================

/// gathers all command line options.
#[derive(Default, Clone)]
pub struct ParseFastqArgs {
    /// compression
    pub nb_bits_by_base: u8,
    /// if we have returnTimes computation
    pub ret_times_args_opt: Option<ReturnTimesArgs>,
    /// args for lmer processing
    pub kmer_args: KmerArgs,
    /// data file
    pub filename:String,
}
