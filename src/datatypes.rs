

pub struct FilePath {
    pub path: String,
    pub headsize: u8,
}
impl FilePath {
    pub fn new(path: &str, headsize: u8) -> FilePath {
        let outpath: String = path.to_string();
        return FilePath{path: outpath, headsize}
    }
}

pub struct Nanopolish {
    pub contig: Vec<String>,
    pub position: Vec<u64>,
    pub reference_kmer: Vec<String>,
    pub read_index: Vec<usize>,
    pub strand: Vec<char>,
    pub event_index: Vec<f32>,
    pub event_level_mean: Vec<f32>,
    pub event_stdv: Vec<f32>,
    pub event_length: Vec<f32>,
    pub model_kmer: Vec<String>,
    pub model_mean: Vec<f32>,
    pub model_stdv: Vec<f32>,
    pub standardized_level: Vec<f32>,
    pub start_idx: Vec<u64>,
    pub end_idx: Vec<u64>,
}

pub struct GtfExons {
    pub transcripts: Vec<String>,
    pub genes: Vec<String>,
    pub exon_number: Vec<u32>,
    pub chr: Vec<String>,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub strand: Vec<char>,
}
impl GtfExons {
    pub fn tid_to_gid(&self, t_id: &String) -> Vec<&String> {
        let indices = self.transcripts
            .iter()
            .enumerate()
            .filter_map(|(index, r)| if r == t_id { Some(index) } else { None })
            .collect::<Vec<_>>();
        return indices.iter().map(|i| &self.genes[*i]).collect::<Vec<_>>();
    }

    pub fn tid_to_chrom(&self, t_id: &String) -> Vec<&String> {
    	let indices = self.transcripts
    		.iter()
    		.enumerate()
    		.filter_map(|(index, r)| if r == t_id { Some(index) } else { None })
    		.collect::<Vec<_>>();
    	return indices.iter().map(|i| &self.chr[*i]).collect::<Vec<_>>();
    }

    pub fn tpos_to_gpos(&self, transcript_id: &String, t_pos: u64) -> (String, u64) {
    	let mut chr: String = "".to_string();
        let mut g_pos: u64 = 0;
        // get all indices which match transcript
        let t_ind = self.transcripts
            .iter()
            .enumerate()
            .filter_map(|(index, r)| if r == transcript_id { Some(index) } else { None })
            .collect::<Vec<_>>();
        let mut cumsum: u64 = 0;
        for t in t_ind{
            // there should only be one strand per transcript... I hope
            // I also hope exons within same gene cannot be overlapping :3
            let strand = self.strand[t];
            // check if tpos is in this exon, exons count up
            // equal to incl 0, it is first exon
            if t_pos >= cumsum {
            	// if we are IN this exon
	            chr = self.chr[t].clone();
	            let start: u64 = self.start[t];
	            let end: u64 = self.end[t];
	            // send output
	            if strand == '+' {
	                g_pos = t_pos - cumsum + start;
	            } else if strand == '-' {
	                g_pos = t_pos - cumsum + start;
	            }        	
            } else {
            	// if we are not, cumulate lengths
            	let exon_len = self.end[t] - self.start[t];
            	cumsum += exon_len;
            }
        }
        (chr, g_pos)
    }
}

pub struct Sequence {
    pub chr: Vec<String>,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub seq: Vec<String>,
}
