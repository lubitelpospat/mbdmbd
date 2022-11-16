use std::fs;
use regex::Regex;
use super::{FilePath, Nanopolish, Sequence, GtfExons};

pub fn read_sequence_file(fa: FilePath) -> Sequence {
	let data = fs::read_to_string(fa.path).expect("Unable to read file");
	let lines = data.split('\n');
	// to delete this symbol: > from id
	let rg = Regex::new(r">").unwrap();
	// also remove version
	let rg2 = Regex::new(r"\..*$").unwrap();
	// vecs
	let mut chr: Vec<String> = Vec::new();
	let mut start: Vec<u64> = Vec::new();
	let mut end: Vec<u64> = Vec::new();
	let mut seq: Vec<String> = Vec::new();
	// seq accumulator
	let mut sequence: String = "".to_string();
	for line in lines {
		if !line.is_empty() {
			let first = line.chars().next().expect("problem");
			if first == '>' {
				let f: Vec<&str> = line.split(' ').collect();
				if !f[0].is_empty() {
					let trns: &str = f[0];
					let tso = rg.replace(trns, "");
					let ts = rg2.replace(&tso, "");
					chr.push(ts.to_string());
					let posf: Vec<&str> = f[2].split(':').collect();
					// parent chr posf[2]
					start.push(posf[3].parse::<u64>().unwrap());
					end.push(posf[4].parse::<u64>().unwrap());
					if !sequence.is_empty() {
						seq.push(sequence.clone());
						sequence = "".to_string();
					}
				}
			} else {
				sequence += line;
			}
		}
	}
	Sequence{chr,start,end,seq}
}

pub fn read_gtf_file(gtf: FilePath) -> GtfExons {
    let data = fs::read_to_string(gtf.path).expect("Unable to read file");
    let lines = data.split('\n');
    // to delete quotes from field values 
    let rg = Regex::new(r#"[\\"]"#).unwrap();
    // vecs
    let mut transcripts: Vec<String> = Vec::new();
    let mut transcripts2: Vec<String> = Vec::new();
    let mut genes: Vec<String> = Vec::new();
    let mut exon_number: Vec<u32> = Vec::new();
    let mut chr: Vec<String> = Vec::new();
    let mut start: Vec<u64> = Vec::new();
    let mut end: Vec<u64> = Vec::new();
    let mut strand: Vec<char> = Vec::new();
    /////////////////// iter
    for (idx, line) in lines.enumerate() {
        if idx > gtf.headsize.into() {
            let cols: Vec<&str> = line.split('\t').collect();
            if !cols[0].is_empty() {
	            let feature_type = cols[2].to_string();
	            if feature_type == "exon" {
	                chr.push(cols[0].to_string());
	                start.push(cols[3].parse::<u64>().unwrap());
	                end.push(cols[4].parse::<u64>().unwrap());
	                strand.push(cols[6].parse::<char>().unwrap());
	                let fields = cols[8].split("; ");
	                for field in fields {
	                    let f: Vec<&str> = field.split(' ').collect();
	                    if f[0] == "gene_id" {
	                    	let gid: &str = f[1];
	                    	let g = rg.replace_all(gid,"");
	                        genes.push(g.to_string());
	                    }
	                    if f[0] == "exon_number" {
	                        let en: &str = f[1];
	                        let exn = rg.replace_all(en,"");
							if let Ok(exnum) = exn.parse::<u32>() {
	                        	exon_number.push(exnum);
							} else {
								println!("{}",exn);
							}
	                    }
	                    if f[0] == "transcript_id" {
	                    	let trns: &str = f[1];
	                    	let ts = rg.replace_all(trns,"").to_string();
     		                let ts2 = rg.replace_all(trns,"").to_string();
	                        transcripts.push(ts);
	                        transcripts2.push(ts2);
	                    }
	                }
	            }
            }
        }
    }
    GtfExons{transcripts, genes, exon_number, chr, start, end, strand}
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
        if idx > np.headsize.into() {
            let cols: Vec<&str> = line.split('\t').collect();
            if !cols[0].is_empty(){
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
    }
    //// return NP object
    Nanopolish {
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
