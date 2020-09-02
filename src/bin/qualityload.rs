//! loads quality compressed in a wavelet matrix
//! provides a server to quality data.
//! statistics provided for free by Wavelet matrix




#[macro_use]
extern crate lazy_static;

// for logging (debug mostly, switched at compile time in cargo.toml)
#[macro_use]
extern crate log;
extern crate simple_logger;


use ::clap::{App, Arg};

#[doc(no_inline)]
use ::std::process;


use ::std::io::{Cursor,Write};
use ::xdr_codec::{Pack};
use ::zmq::{Message};

use ::std::io;
use ::std::cmp;
use ::std::mem;


// our modules
extern crate kmerutils;
use kmerutils::quality::*;
use kmerutils::qserverclient::*;

const DEBUG : u64 = 0;

lazy_static! {
    #[allow(dead_code)]
    static ref LOG: u64 = {
        let res = init_log();
        res
    };
}

// install a logger facility
fn init_log() -> u64 {
    simple_logger::init().unwrap();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}




///
///        Server statistics
///

#[allow(dead_code)]
struct ServerStatistics {
    bytes_sent: u64,
    ///
    nb_request: u64,
}

//
// Dialog:
// Client sends:  request_handle , request_code. then arguments dependant on request type.
// Server responds : request_handle , status code : then a quality vector if status code is OK. 


struct QualityServer<'a> {
    _myport : u64,
    ///
    // in fact a REP socket. Note Socket keeps an Arc to Context which it comes from.
    socket : &'a zmq::Socket,
    ///
    qualities: &'a Vec<QSequenceWM>, 
}


impl<'a> QualityServer<'a> {
    fn new(port:u64, s: &'a zmq::Socket, qseqvec: &'a Vec<QSequenceWM>) -> QualityServer<'a> {
        QualityServer{_myport: port, socket: s, qualities: qseqvec}
    } // end of new

    
    fn get_request_type(&self, cursor_in : &mut Cursor<Vec<u8>> ) -> Result<RequestId, StatusCode> {
        info!("server entering get_request_type");
        // get id
        let handle : Qhandle;
        let rescode: xdr_codec::Result<(u64, usize)>  = xdr_codec::Unpack::unpack(cursor_in);
        if let Ok(decode) = rescode {
            handle = decode.0;
            println!("server get_quality request handle  {}", handle);
        }
        else {
            // real bad case, could not even get handle
            return Err(StatusCode::ErrXdr);
        }
        // get code for request
        let opcode_u64;
        // get value consumed and length in bytes, recall xdr_codec::Result<T> = result::Result<T, Error>;
        let rescode: xdr_codec::Result<(u64, usize)>  = xdr_codec::Unpack::unpack(cursor_in);
        match rescode {
            Ok(decode) => {opcode_u64 = decode.0},
            Err(_) => { return Err(StatusCode::ErrXdr); }
        }        
        // we have request type
        match opcode_u64 {
            0 => return Ok(RequestId(handle, RequestCode::Zero)),
            1 => return Ok(RequestId(handle, RequestCode::GetQRead)),
            2 => return Ok(RequestId(handle, RequestCode::GetQBlock)),
            9 => return Ok(RequestId(handle, RequestCode::Exit)),           
            _ => {
                println!("server get_opcode, got unknown opcode : {}", opcode_u64);
                return Err(StatusCode::ErrArg);
            },          
        }
    } // end of get_request_type


    fn reply_error(&self, err: ResponseStatus) {
        // make a cursor for response. It must contains request id and quals.len() and then quals
        let len = mem::size_of::<Qhandle>() +  mem::size_of::<u64>()  + 4;
        let vec:Vec<u8> = Vec::with_capacity(len);
        // if we cannot send err msg we abort when unwrap!
        let mut out_cursor = Cursor::new(vec);
        err.0.pack(&mut out_cursor).unwrap();
        (err.1 as u64).pack(&mut out_cursor).unwrap();
        //
        // now we send cursor to the socket
        //
        let rawbuf = out_cursor.get_ref().as_slice();
        let msg = Message::from_slice(rawbuf).unwrap();
        info!("reply_error : sending err msg ");
        let res = self.socket.send_msg(msg, 0);
        if res.is_err() {
            println!(" could not even send err msg") 
        }
    } // reply_err

    // 
    fn send_quality_slice(&self, request_id : RequestId, quals:&[u8]) -> io::Result<()> {
        debug!("send_quality_slice : entering");
        let mut _bytes_written : usize = 0;
        // make a cursor for response. It must contains request id and quals.len() and then quals
        let len = mem::size_of::<Qhandle>() +  mem::size_of::<u64>()  + quals.len() + 8 + 4;
        let vec:Vec<u8> = Vec::with_capacity(len);
        let mut out_cursor = Cursor::new(vec);
        request_id.0.pack(&mut out_cursor).unwrap();
        (StatusCode::Ok as u64).pack(&mut out_cursor).unwrap();
        // now u8 do not need to be xdr encoded!, just write in cursor, preceded by len!
        (quals.len() as u64).pack(&mut out_cursor).unwrap();
        out_cursor.write_all(quals)?;
        out_cursor.flush().unwrap();
        trace!(" nb quals bytes written in cursor {:?} ", out_cursor.position());
        //
        // now we send cursor to the socket
        //
        let rawbuf = out_cursor.get_ref().as_slice();
        let msg = Message::from_slice(rawbuf).unwrap();
        debug!("send_quality_slice: sending message nb quals {}", quals.len());
        let res = self.socket.send_msg(msg, 0);
        if res.is_err() {
            println!(" get_initial_handle : error sending message") 
        }
        debug!("send_quality_slice : message sent, exiting");
        Ok(())
    } // end send_quality_slice



    
    // ============== treatment of requests =====================

    // response to get a block of a Qsequence
    fn quality_sequence_block_request(&self, request_id : RequestId, cursor_in : &mut Cursor<Vec<u8>> ) {
        info!("server entering block request_type");
        // must decode the num of sequence and send response to client
        let seqnum;
        //
        let mut resu64_xdr: xdr_codec::Result<(u64, usize)>  = xdr_codec::Unpack::unpack(cursor_in);
        match resu64_xdr {
            Ok(decode) => seqnum = decode.0,
            Err(_) => {
                info!("\n server getting block request : num of sequence decoding failed");
                self.reply_error(ResponseStatus(request_id.0, StatusCode::ErrGen));
                return;
            },
        }
        let begin_block;
        let end_block;
        //
        resu64_xdr  = xdr_codec::Unpack::unpack(cursor_in);
        match resu64_xdr {
            Ok(decode) => {begin_block = decode.0 as usize},
            Err(_) => {
                info!("\n server getting block request : begin of block sequence decoding failed");
                self.reply_error(ResponseStatus(request_id.0, StatusCode::ErrGen));
                return;
            },
        }        
        //
        resu64_xdr= xdr_codec::Unpack::unpack(cursor_in);
        match  resu64_xdr {
            Ok(decode) => {end_block = decode.0 as usize},
            Err(_) => {
                warn!("\n server getting request num of sequence  decoding failed");
                self.reply_error(ResponseStatus(request_id.0, StatusCode::ErrGen));
                return;
            },
        }
        //
        // we can check for the read_num ...
        //
        if seqnum as usize >= self.qualities.len() {
            self.reply_error(ResponseStatus(request_id.0, StatusCode::ErrArg));
            return;
        }
        // get the QSequenceRaw, we need to decompress, cause a copy.
        let seqraw = self.qualities[seqnum as usize].decompress();
        if DEBUG > 0 {
            for i in 0..cmp::min(10, seqraw.len()) {
                println!("send qual {} {}" , i, seqraw.qseq[i]);
            }
        }
        let _res = self.send_quality_slice(request_id, &seqraw.qseq.as_slice()[begin_block..end_block]);
    } // end of quality_sequence_block_request


    
    // response to a quality sequence request.
    fn quality_sequence_request(&self, request_id : RequestId, cursor_in : &mut Cursor<Vec<u8>> ) {
        info!("server entering quality_sequence_request");
        // must decode the num of sequence and send response to client
        let seqnum;
        //
        let resu64_xdr: xdr_codec::Result<(u64, usize)> = xdr_codec::Unpack::unpack(cursor_in);
        match resu64_xdr {
            Ok(decode) => { seqnum = decode.0;
            },
            Err(_e) => {
                info!("server getting request num of sequence decoding failed");
                self.reply_error(ResponseStatus(request_id.0, StatusCode::ErrGen));
                return;
            },
        }
        //
        if seqnum as usize >= self.qualities.len() {
            warn!("server got request for sequence num {} greater than number of qualities {} ", seqnum, self.qualities.len());
            self.reply_error(ResponseStatus(request_id.0, StatusCode::ErrArg));
            return;
        }
        // get the QSequenceRaw, we need to decompress, cause a copy.
        let seqraw = self.qualities[seqnum as usize].decompress();
        if DEBUG > 0 {
            for i in 0..cmp::min(10, seqraw.len()) {
                println!("send qual {} {}" , i, seqraw.qseq[i]);
            }
        }        // we can check for the read_num ...
        let _res = self.send_quality_slice(request_id, seqraw.qseq.as_slice());
    } // end of quality_sequence_request



    
    // We must get first a RequestId
    fn decode_and_treat_msg(&self, cursor_in : &mut Cursor<Vec<u8>> ) {
        info!("server entering decode_and_treat_msg");
        let request_id; 
        // get code for request
        let request_t = self.get_request_type(cursor_in);
        match request_t {
            Ok(request_tmp) => {
                request_id = RequestId(request_tmp.0, request_tmp.1);
            },
            Err(e) => {
                self.reply_error(ResponseStatus(0,e));
                return;
            }
        }
        // 
        match request_id.1 {
            RequestCode::GetQRead => {
                let _res = self.quality_sequence_request(request_id, cursor_in);
            }
            RequestCode::GetQBlock => {
                self.quality_sequence_block_request(request_id, cursor_in);
            }
            RequestCode::Exit => {
                println!("\n server received order to quit ");               
                process::exit(0);
            }
            _ => { // should log bad request, send sthing to client?
                println!("server get_opcode, got unknown opcode");
            }           
        }
    } // end of function decode_and_treat_msg


   
    fn service_loop(&self) -> usize {
        info!("entering service loop");
        loop {
            // we treat requests in sequence to achieve data/request consistency
            // get the message as a Result<Vec<u8> > , do not expect more
            // 0 means no more , msg has only one part!
            // recv_bytes returns Result<Vec<u8>>
            let res_vec = self.socket.recv_bytes(0);
            match res_vec {
                Ok(v) => self.decode_and_treat_msg(&mut Cursor::new(v)),
                Err(e) => println!("\n reception of message failed {} ", e),          
            }
        } // end of loop block
        //
    }
} // end of impl block for QualityServer



//////////////////////////////////////////////////////////////////////////////////////////

fn main() {
    
    if cfg!(verbose_1 = "1") {
        println!(" appel main : ");
    }

    
    // the reference to LOG will force the call to lazy_static! call to init_log to get LOG initialized.
    if *LOG != 1 {
        println!(" LOG = {:?}", *LOG);
    }
    let fname;


    let matches = App::new("qualityloader")
        .arg(Arg::with_name("file")
             .long("file")
             .short("f")
             .takes_value(true)
             .help("expecting a fastq file"))
        .arg(Arg::with_name("port")
             .long("port")
             .short("p")
             .takes_value(true)
             .help("expecting a port number (>= 1024 and < 65536 "))
        .arg(Arg::with_name("wavelet")
             .long("wavelet")
             .short("w")
             .help("wavelet compression for long reads"))
        .get_matches();


    if matches.is_present("file") {
        fname = matches.value_of("file").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got filename , {}", fname);
    }
    else {
        println!("-f filename is mandatory");
        println!(" usage qualityloade -f name --wawelet (or -w) --port (-p) num");
        process::exit(1);
    }

    let qseqvec : Vec<QSequenceWM>;
    
    // if present get port else default affectation
    let mut portnum = 4766;
    if matches.is_present("port") {
        portnum = matches.value_of("port").ok_or("bad value").unwrap().parse::<u64>().unwrap_or(portnum);
        info!(" will try port number {}", portnum);
    }
    

    if matches.is_present("wavelet") {
        if let Ok(qseqvec_r) = load_quality_wm(&fname) {
            qseqvec = qseqvec_r;
            println!("got nb read {} : ", qseqvec.len());
        }
        else {
            println!("could not read qualities");
            process::exit(1); 
        }
    }
    else {
        println!("only wavelet matrix option implemented : pass --wavelet or -w option");
        process::exit(1);
    }
    //
    let context = zmq::Context::new();
    let responder = context.socket(zmq::REP).unwrap();

    let mut connection_str = String::from("tcp://*:");
    connection_str.push_str(&portnum.to_string());
    println!("server binding to : {} ", connection_str);
    if responder.bind(&connection_str).is_err() {
        println!("server cannot bind to socket");
        process::exit(1);
    }
    //
    // quality are loaded, port is bind...
    //
    let server = QualityServer::new(portnum, &responder, &qseqvec);
    let _status = server.service_loop();
}   // end of main
