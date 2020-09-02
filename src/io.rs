//! this file gathers io , serialization of data to sent to julia ...

// for needletail


extern crate time;
use self::time::*;

use std::path::Path;

use crate::statutils::*;
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

    let mut return_times_opt: Option<ReturnTimesWriter> = None;
    if parsed_args.ret_times_args_opt.is_some() {
        // now we check for return times options
        let mut rt_file_name = parsed_args.filename.clone();
        rt_file_name.push_str(".ret_times.bin");
        if let Some(pos) = rt_file_name.rfind('/') {
            rt_file_name = rt_file_name.split_off(pos+1);
        }
        println!("return times in file : {} ", rt_file_name);
        let ret_times_arg:ReturnTimesArgs = parsed_args.ret_times_args_opt.expect("no return times option").clone();
        let w_size: u8 = ret_times_arg.window_size;
        if w_size == 0 {
            println!(" window size for return times stats must be between 1 and 255");
            panic!();
        }
        let searched_base = ret_times_arg.searched_base;
        return_times_opt = Some(ReturnTimesWriter::new(rt_file_name, searched_base, w_size).unwrap());
    }
    
    let start_t = Instant::now();
    let mut reader = needletail::parse_fastx_file(&path).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        n_bases += seqrec.num_bases();
        let nb_bad = count_non_acgt(&seqrec.seq());
        nb_bad_bases = nb_bad_bases + nb_bad;
        if nb_bad == 0 {
            match return_times_opt {
                // we must capture a mutable reference in the match!!
                Some(ref mut ret) => {
                    if let Err(_) = ret.analyze(seq_array.len(), &seqrec.seq()) {
                        println!("error analyzing sequence");
                    }
                },
                None => {},
            }
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

