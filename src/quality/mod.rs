//! This module gathers quality mapping, compression and service

pub mod qserverclient;
pub mod quality;

#[cfg(feature = "withzmq")]
pub mod qualclient;
