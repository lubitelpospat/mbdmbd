extern crate bio;
mod datatypes;
use datatypes::*;
mod import_functions;
use import_functions::*;

fn main() {
    // first we need paths to the nanopolish and genome files
    let npa_path: FilePath = FilePath::new("./data-files/mouse_antiadar_eventalign.txt",  1);
    let npb_path: FilePath = FilePath::new("./data-files/Mus_musculus.GRCm39.cdna.all.fa", 1);
    let genome_path: FilePath = FilePath::new("./data-files/Mus_musculus.GRCm39.dna.toplevel.fa", 0);
    let gtf_path: FilePath = FilePath::new("./data-files/Mus_musculus.GRCm39.107.gtf", 5);

    // then we need to read those into memory
    let np_a = read_nanopolish_file(npa_path);
    let np_b = read_nanopolish_file(npb_path);
    // then we need to write TWO functions for the gtf
    // tID_to_gID() and tPos_to_gPos()
    let (tid_gid, gtf_exons) = read_gtf_file(gtf_path);
    // we need to be sure that tPos_to_gPos works properly by checking the kmer
    // this means we need also to load the transcriptome fa and genome fa to compare kmer
    //let (transcriptome_fa, genome_fa) = read_sequence_file();
    // then, we iterate over all the nanopolish lines from each file
    // we get genomic coords
}

