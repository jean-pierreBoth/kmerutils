//! This module contains code implemeting a network client
//! to quality requests for a server responsible for loading and servicing
//! qualities for a data file



extern crate zmq;
extern crate xdr_codec;
extern crate rand;


use std::mem;
use std::cmp;
use std::process;

use self::xdr_codec::{Pack, Unpack};

use std::io::{Cursor,Read};
use self::zmq::{Message};

use rand::{Rng,thread_rng};

// our mod

use crate::quality::*;
use crate::qserverclient::*;


const DEBUG : u64 = 0;




pub struct QualityClient {
    /// name of server
    _server: String,
    /// port for contacting server
    _port: u64,
    // in fact a REQ socket. Note Socket keeps an Arc to Context which it comes from.
    socket: zmq::Socket,
    //
    _filename: String,
    //
}  // end of struct QualityClient




impl QualityClient {
    pub fn new(server : String, port : u64, filename: String) -> QualityClient {
        let context = zmq::Context::new();
        let socket = context.socket(zmq::REQ).unwrap();
        let mut connection_str = String::from("tcp://");
        // concatenate servername and port converted to a string
        connection_str.push_str(&server);
        connection_str.push(':');
        connection_str.push_str(&port.to_string());
        log::info!(" connecting to {} ... ", connection_str);
        let resconn = socket.connect(&connection_str);
        assert!(resconn.is_ok());
        match resconn {
            Ok(_) => {
                println!(" connect Ok");
            }
            Err(e) => {
                println!(" got an error in connect : {}", e);
                process::exit(1);
            }
        }
        QualityClient{_server:server, _port:port, socket:socket, _filename:filename}
    }  // end of new


    // CAVEAT  this is not thread safe as thread_rng is not Sync!!
    fn get_request_handle(&self) -> Qhandle {
        thread_rng().gen::<u64>() as Qhandle
    }
    
    /// return a Quality sequence given its number
    pub fn get_quality_sequence(&self, numseq:u64) -> Result<QSequenceRaw, StatusCode> {
        // we must first send a requestId
        log::debug!("get_quality requesting qsequence {}", numseq);
        let len = mem::size_of::<RequestId>()  + 8 + 4;
        let vec:Vec<u8> = Vec::with_capacity(len);
        let mut buffer  = Cursor::new(vec);
        // construct a request with Qhandle and request code
        let handle = self.get_request_handle();
        let req_code = RequestCode::GetQRead;
        log::debug!("get_quality request handle  {}", handle);
        if handle.pack(&mut buffer).is_err() {
            return Err(StatusCode::ErrXdr);
        }
        if (req_code as u64).pack(&mut buffer).is_err() {
            return Err(StatusCode::ErrXdr);
        }
        // now send numseq
        if numseq.pack(&mut buffer).is_err() {
            return Err(StatusCode::ErrXdr);
        }
        // check how many have been packed.
        
        // now have a full buffer to transmit, get a &[u8] slice from cursor
        let rawbuf = buffer.get_ref().as_slice();
        // make up a zmq message
        let resmsg = Message::from_slice(rawbuf);
        if resmsg.is_err() {
            println!(" construction of msg from slice failed ... ");
        }        
        let msg = resmsg.unwrap();
        log::debug!(" sending message ");
        let _res = self.socket.send_msg(msg, 0).unwrap();
        //
        // must wait for handle back from server
        //
        if DEBUG > 0 {
            println!(" waiting for server ... "); 
        }
        // the 0 flag is for blocks...
        if let Ok(v) = self.socket.recv_bytes(0) {
            log::info!(" client receiving a response .. nb_bytes : {} ", v.len());
            let mut qualv: Vec<u8>;
            let mut cursor = Cursor::new(v);
            // get handle (Qhandle  = u64) and req_code (as u 64) and after data for a QSequenceRaw
            let (handle_got, _) = u64::unpack(&mut cursor).unwrap();
            assert_eq!(handle_got, handle);
            // now get response status
            let (ret_code , _) = u64::unpack(&mut cursor).unwrap();
            if ret_code != StatusCode::Ok as u64 {
                log::info!("get_quality got an error status from server");
                return Err(StatusCode::ErrGen); 
            }
            // now get len of quality vector
            if let Ok((nb_qual, _)) = u64::unpack(&mut cursor) {            
                // get QsequenceRaw. The u8 slice is not xdr encoded. We get it directly from cursor
                // after allocating a vector of the right size
                log::info!("got a qseq size: {} ", nb_qual);
                qualv = Vec::with_capacity(nb_qual as usize);
                // the following line is crucial as arg is supposed to be an initialized slice and the only
                // function cursor can use is len() !!
                unsafe { qualv.set_len(nb_qual as usize) };
                cursor.read_exact(qualv.as_mut_slice()).unwrap();
                if DEBUG > 0 {
                    for i in 0..cmp::min(10, qualv.len()) {
                        println!("client got qual {} {}" , i, qualv[i]);
                    }
                }                
                // we return our quality vector
                log::info!("returning a QSequenceRaw with size : {:?} ", cursor.position());
                return Ok(QSequenceRaw{read_num: numseq as usize, qseq: qualv});
            }
            else {
                return Err(StatusCode::ErrGen); 
            }
        }
        else {
            return Err(StatusCode::ErrGen); 
        }
    } // end of get_quality_sequence

    
    
} // end of implementation QualityClient
