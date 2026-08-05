use bio::alphabets::dna;
pub fn clip_interval(start: u64, end: u64, length: usize, flank: usize) -> (usize, usize) {
    let s = start.saturating_sub(flank as u64) as usize;
    let e = (end as usize + flank).min(length);
    (s.min(e), e)
}
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    dna::revcomp(seq)
}
