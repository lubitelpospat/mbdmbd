use std::fs;
use super::{FilePath, Nanopolish, Sequence, TidGid, GtfExons};

// pub fn read_fa_file(fa: FilePath) -> Sequence {
// }

pub fn read_gtf_file(gtf: FilePath) -> (TidGid, GtfExons) {
    let data = fs::read_to_string(gtf.path).expect("Unable to read file");
    let lines = data.split("\n");
    // vecs
    let mut transcripts: Vec<String> = Vec::new();
    let mut transcripts2: Vec<String> = Vec::new();
    let mut genes: Vec<String> = Vec::new();
    let mut exon_number: Vec<u8> = Vec::new();
    let mut chr: Vec<u8> = Vec::new();
    let mut start: Vec<u64> = Vec::new();
    let mut end: Vec<u64> = Vec::new();
    let mut strand: Vec<char> = Vec::new();
    /////////////////// iter
    for (idx, line) in lines.enumerate() {
        if idx > gtf.headsize.into() {
            let cols: Vec<&str> = line.split("\t").collect();
            let feature_type = cols[2].to_string();
            if feature_type == "exon" {
                chr.push(cols[0].parse::<u8>().unwrap());
                start.push(cols[3].parse::<u64>().unwrap());
                end.push(cols[4].parse::<u64>().unwrap());
                strand.push(cols[6].parse::<char>().unwrap());
                let fields = cols[8].split("; ");
                for field in fields {
                    let f: Vec<&str> = field.split(" ").collect();
                    if f[0] == "gene_id" {
                        genes.push(f[1].to_string());
                    }
                    if f[0] == "exon_number" {
                        let en: String = f[1].to_string();
                        exon_number.push(en.parse::<u8>().unwrap());
                    }
                    if f[0] == "transcript_id" {
                        transcripts.push(f[1].to_string());
                        transcripts2.push(f[1].to_string());
                    }
                }
            }
        }
    }
    return (
        TidGid{transcripts, genes},
        GtfExons{transcripts: transcripts2, exon_number, chr, start, end, strand}
    )
}

pub fn read_nanopolish_file(np: FilePath) -> Nanopolish {
    let data = fs::read_to_string(np.path).expect("Unable to read file");
    let lines = data.split('\n');
    // np vecs
    let mut contig: Vec<String> = Vec::new();
    let mut position: Vec<u64> = Vec::new();
    let mut reference_kmer: Vec<String> = Vec::new();
    let mut read_index: Vec<usize> = Vec::new();
    let mut strand: Vec<char> = Vec::new();
    let mut event_index: Vec<f32> = Vec::new();
    let mut event_level_mean: Vec<f32> = Vec::new();
    let mut event_stdv: Vec<f32> = Vec::new();
    let mut event_length: Vec<f32> = Vec::new();
    let mut model_kmer: Vec<String> = Vec::new();
    let mut model_mean: Vec<f32> = Vec::new();
    let mut model_stdv: Vec<f32> = Vec::new();
    let mut standardized_level: Vec<f32> = Vec::new();
    let mut start_idx: Vec<u64> = Vec::new();
    let mut end_idx: Vec<u64> = Vec::new();
    ///////////////////////////// ITERATE
    for (idx, line) in lines.enumerate() {
        if idx > np.headsize.into(){
            let cols: Vec<&str> = line.split("\t").collect();
            contig.push(cols[0].to_string());
            position.push(cols[1].parse::<u64>().unwrap());
            reference_kmer.push(cols[2].to_string());
            read_index.push(cols[3].parse::<usize>().unwrap());
            strand.push(cols[4].parse::<char>().unwrap());
            event_index.push(cols[5].parse::<f32>().unwrap());
            event_level_mean.push(cols[6].parse::<f32>().unwrap());
            event_stdv.push(cols[7].parse::<f32>().unwrap());
            event_length.push(cols[8].parse::<f32>().unwrap());
            model_kmer.push(cols[9].to_string());
            model_mean.push(cols[10].parse::<f32>().unwrap());
            model_stdv.push(cols[11].parse::<f32>().unwrap());
            standardized_level.push(cols[12].parse::<f32>().unwrap());
            start_idx.push(cols[13].parse::<u64>().unwrap());
            end_idx.push(cols[14].parse::<u64>().unwrap());
        }
    }
    //// return NP object
    return Nanopolish {
        contig,
        position,
        reference_kmer,
        read_index,
        strand,
        event_index,
        event_level_mean,
        event_stdv,
        event_length,
        model_kmer,
        model_mean,
        model_stdv,
        standardized_level,
        start_idx,
        end_idx,
    }
}