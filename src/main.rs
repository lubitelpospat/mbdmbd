mod datatypes;
use datatypes::*;
mod import_functions;
use import_functions::*;
use std::process::exit;
use std::fs;
use regex::Regex;
use substring::Substring;
use std::collections::HashMap;

fn main() -> std::io::Result<()> {
    // first we need paths to the nanopolish and genome files
    let npa_path: FilePath = FilePath::new("./data-files/mouse_antiadar_eventalign-test.txt", 1);
    let npb_path: FilePath = FilePath::new("./data-files/mouse_scramble_eventalign-test.txt", 1);
    let gtf_path: FilePath = FilePath::new("./data-files/Mus_musculus.GRCm39.107.gtf", 5);
    let traome_path: FilePath = FilePath::new("./data-files/Mus_musculus.GRCm39.cdna.all.fa", 0);
	let genome_path: FilePath = FilePath::new("./data-files/Mus_musculus.GRCm39.dna.toplevel.fa", 0);
	let output_path: &str = "./hashmap-test.txt";

    // then we need to read those into memory

    //let np_a = read_nanopolish_file(npa_path);
    //println!("{:?}", np_a.contig[6]);
    //let np_b = read_nanopolish_file(npb_path);
    //println!("{:?}", np_b.contig[22]);

    // then we need to write TWO functions for the gtf
    // tID_to_gID() and tPos_to_gPos()

    let gtf_exons = read_gtf_file(gtf_path);
    //println!("{:?}",tid_gid.transcripts[0]);
    //println!("{:?}",gtf_exons.start[90]);

    // we need to be sure that tPos_to_gPos works properly by checking the kmer
    // this means we need also to load the transcriptome fa and genome fa to compare kmer

    let traome_fa = read_sequence_file(traome_path);
    //println!("{:?}", traome_fa.chr[64]);
    //println!("{:?}", traome_fa.seq[64]);
    let genome_fa = read_sequence_file(genome_path);
	//println!("{:?}", genome_fa.chr[5]);
	// store in (geneID, chrom, pos) format
	let mut genpos_storage: HashMap<(String,String,u64), Vec<Vec<f32>>> = HashMap::new();

    // then, we iterate over all the nanopolish lines from each file
	let data = fs::read_to_string(npb_path.path).expect("Unable to read file");
	let lines = data.split('\n');
	// regex to remove .num from after tid
	let rg = Regex::new(r"\..*$").unwrap();
	// iter
	for (idx, line) in lines.enumerate() {
		if idx > npb_path.headsize.into() {
			let cols: Vec<&str> = line.split('\t').collect();
			if !cols[0].is_empty(){		
				// check that it's a legal nanopolish aligment, model and ref kmer match
			    // ref kmer
			    let refkmer = cols[2].to_string();
				let modelkmer = cols[9].to_string();
				if refkmer == modelkmer {
					// here I make a hashmap to store the voltage distribution from each kmer
					//let mut kmer_means: HashMap<String, Vec<f32>> = HashMap::new();
					//let mut kmer_stdvs: HashMap<String, Vec<f32>> = HashMap::new();
					//let mut kmer_lengths: HashMap<String, Vec<f32>> = HashMap::new();
					//let mut kmer_dwells: HashMap<String, Vec<f32>> = HashMap::new();
				    // get these: 
				    // transcriptomic id
				    let trns: &str = cols[0];
				    let ts = rg.replace(trns, "");
				    let tid = ts.to_string();
				    // transcriptomic position
					let tpos = cols[1].parse::<u64>().unwrap();
				    // then convert to genomic pos
				   	let gid = gtf_exons.tid_to_gid(&tid)[0];
				   	let chrom = gtf_exons.tid_to_chrom(&tid)[0];
				   	let (chr, gpos) = gtf_exons.tpos_to_gpos(&tid, tpos);
				    //// transcript kmer
				    //let t_idx = traome_fa.chr.iter().position(|r| r == &tid).unwrap();
				    //let ks = tpos as usize;
				    //let ke = tpos as usize + 5;
				    //let tkmer = traome_fa.seq[t_idx].substring(ks,ke);
				    //// genome kmer
				    let g_idx = genome_fa.chr.iter().position(|r| r == chrom).unwrap();
				    //let gs = gpos as usize - 1;
				    //let ge = gpos as usize + 4;
				    //let gkmer = genome_fa.seq[g_idx].substring(gs,ge);
					//////// event data
				    let evend = cols[13].parse::<f32>().unwrap();
				    let evstart = cols[14].parse::<f32>().unwrap();
					let evlength = evstart-evend;
				    let evdwell = cols[8].parse::<f32>().unwrap() * evlength;
				    let evstdv = cols[7].parse::<f32>().unwrap() * evlength;
				    let evmean = cols[6].parse::<f32>().unwrap() * evlength;

				    // try to get array of means from distros
				    // only use new v if you can't find another one
			    	let evdata: Vec<f32> = vec![evmean, evstdv, evlength, evdwell];
			    	let new_data: Vec<Vec<f32>> = vec![vec![evmean, evstdv, evlength, evdwell]];		    	
					// put them into pos
					genpos_storage
						.entry((gid.to_string(), chr, gpos))
						.and_modify(|vecvec| vecvec.push(evdata))
						.or_insert(new_data);
				}
			}			
		}
	}
	println!("{:?}", genpos_storage);
	//fs::write(output_path, format!("{:?}", kmer_distros))?;
	Ok(())
}
