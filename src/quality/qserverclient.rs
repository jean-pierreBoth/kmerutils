//! This file contains some struct that encode client server interactions

///
/// a request must have an Id
///
pub type Qhandle = u64;

///
///  RequestType
///
/// The basic request are
///       . to get quality sequence for a read number >= 0
///       . get the quality for a base in a sequence
///       . shutdown EXIT
///  But wavelet matrix provides for counting , quantile estimator, median and so on
///  which can be useful to assess global quality for a read.
///
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RequestCode {
    Zero = 0,
    GetQRead = 1,  // basic request for a whole read
    GetQBlock = 2, // get a part of a sequence from begin to end (excluded)
    GetQBase = 3,  // get just a base (must know readnum and base pos)
    Exit = 9,      // exit server
    Unknown = 10,
}

pub struct RequestId(pub Qhandle, pub RequestCode);

/// The enum representing an enum of arguments for various client requests.
pub enum RequestArg {
    /// request is a sequence described by num
    Seq(u64),
    /// request is a block we expect as args , a sequence num , begin and end of seq
    Block(u64, u64, u64),
    /// request is a base (described by a read num and a base pos)
    Base(u64, u64),
}

#[derive(Clone, Debug, PartialEq)]
pub enum StatusCode {
    Ok = 0,
    ErrGen = 1, // generic error
    ErrXdr = 2,
    ErrArg = 3, // bad argument for request
}

pub struct ResponseStatus(pub Qhandle, pub StatusCode);
