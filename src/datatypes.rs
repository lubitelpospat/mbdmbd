

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

pub struct TidGid {
    pub transcripts: Vec<String>,
    pub genes: Vec<String>,
}
impl TidGid {
    pub fn tid_to_gid(&self, tID: String) -> Vec<&String> {
        let indices = self.transcripts
            .iter()
            .enumerate()
            .filter_map(|(index, r)| if r == &tID { Some(index) } else { None })
            .collect::<Vec<_>>();
        return indices.iter().map(|i| &self.genes[*i]).collect::<Vec<_>>();
    }
}

pub struct GtfExons {
    pub transcripts: Vec<String>,
    pub exon_number: Vec<u8>,
    pub chr: Vec<u8>,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub strand: Vec<char>,
}
impl GtfExons {
    pub fn t_pos_to_g_pos(&self, transcript_id: String, t_pos: u64) -> u64 {
        let mut g_pos: u64 = 0;
        // get all indices which match transcript
        let t_ind = self.transcripts
            .iter()
            .enumerate()
            .filter_map(|(index, r)| if r == &transcript_id { Some(index) } else { None })
            .collect::<Vec<_>>();
        let mut cumsum: u64 = 0;
        for t in t_ind{
            // there should only be one strand per transcript... I hope
            // I also hope exons within same gene cannot be overlapping :3
            let strand = self.strand[t];
            // this should start at the lowest exon number...
            let exon_len = self.end[t] - self.start[t];
            cumsum += exon_len;
            // if THIS POSITION is WITHIN this transcript
            if t_pos > self.start[t] && t_pos < self.end[t] {
                // send output
                if strand == '+' {
                    g_pos = &self.start[t] + &t_pos;
                } else if strand == '-' {
                    g_pos = &self.end[t] - &t_pos;
                }
            }
        }
        return g_pos;
    }
}

pub struct Sequence {
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub seq: Vec<char>,
}