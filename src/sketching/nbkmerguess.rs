/// try to guess an upper bound for nb kmer encountered in large seq or files (tens of gb bases)
/// 


use crate::base::sequence::Sequence;

// We need a guess to allocate HashMap used with Kmer Generation
// for very long sequence we must avoid nb_kmer to sequence length! Find a  good heuristic
pub(crate) fn get_nbkmer_guess(seq : &Sequence) -> usize {
    let nb = 100_000_000 * (1usize + seq.size().ilog2() as usize);
    let nb_kmer = seq.size().min(nb);
    return nb_kmer;
} // end of get_nbkmer_guess



// We need a guess to allocate HashMap used with Kmer Generation
// for vector of sequenc coming from a non concatnated file, we must avoid nb_kmer to sequence length! Find a  good heuristic
pub(crate) fn get_nbkmer_guess_seqs(vseq : &Vec<&Sequence>) -> usize {
    let total_nb_base = vseq.iter().fold(0, |acc, seq | acc+seq.size());
    let nb_kmer = 10_000_000 * (1usize + total_nb_base.ilog2() as usize);
    nb_kmer
}  // end of get_nbkmer_guess_seqs
