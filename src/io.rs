//! this file gathers io , serialization of data to sent to julia ...

// for needletail


extern crate time;
use self::time::*;

use std::path::Path;

use crate::sequence::*;
use crate::parsearg::*;


pub fn parse_with_needletail(parsed_args: ParseFastqArgs) ->  std::result::Result<Vec<Sequence>, &'static str > {
    //
    println!("parsing with needletail file : {} ", parsed_args.filename);
    let path = Path::new(&parsed_args.filename);
    let f_info_res = path.metadata();
    let filesize:u64;
    match f_info_res {
        Ok(meta) => {
            filesize = meta.len();           
        },
        Err(_e) => {
            println!("file does not exist: {:?}", parsed_args.filename);
            return Err("file does not exist");
        },
    }
    //
    let nb_bits:u8 = parsed_args.nb_bits_by_base;
    println!("using nb bits by base : {} ", nb_bits);
    // default size is set to half file size (beccause of quality) divided by 1000
    let default_len = filesize/(2*1_000);
    let mut seq_array : Vec<Sequence > = Vec::with_capacity(default_len as usize);
    let mut n_bases = 0;
    let mut nb_bad_bases = 0;
    let mut nb_bad_read = 0;

    
    
    let start_t = Instant::now();
    let mut reader = needletail::parse_fastx_file(&path).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        n_bases += seqrec.num_bases();
        let nb_bad = count_non_acgt(&seqrec.seq());
        nb_bad_bases = nb_bad_bases + nb_bad;
        if nb_bad == 0 {
            let newseq = Sequence::new(&seqrec.seq(), nb_bits);
            seq_array.push(newseq);
        }
        else {
            nb_bad_read = nb_bad_read+1;
        }
        if seq_array.capacity() <= seq_array.len() + 100  {
            let old_len = seq_array.len() as f64;
            seq_array.reserve((old_len * 1.5) as usize);
        }
        if seq_array.len() % 2000000 == 0 {
            println!(" nb rec read = {} ", seq_array.len());
        }
        //
    }
    //
    println!(" shrinking");
    seq_array.shrink_to_fit();
    println!(" shrinked ");
    //
    let elapsed_t = start_t.elapsed().whole_seconds();
    println!(" elapsed time (s) in parse_with_needletail {} ", elapsed_t);
    //
    println!(" nb rec loaded = {} ", seq_array.len());
    println!("nb_bases {:?}", n_bases);
    println!("nb_bad_bases {:?}", nb_bad_bases);
    println!("nb_bad_read {:?}", nb_bad_read);
    //
    return Ok(seq_array);
}  // end of parse_with_needletail

