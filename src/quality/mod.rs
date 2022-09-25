//! This module gathers quality mapping, compression and service 

pub mod quality;
pub mod qserverclient;

#[cfg(withzmq)]
pub mod qualclient;
