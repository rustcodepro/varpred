use crate::expr::FastaRead;
use crate::expr::FastaStruct;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

/*
 Gaurav Sablok
 codeprog@icloud.com
*/

impl FastaRead {
    pub fn fastaread(&self) -> Result<Vec<FastaStruct>, Box<dyn Error>> {
        let fileopen = File::open(self.pathname.clone()).expect("file not present");
        let fileread = BufReader::new(fileopen);
        let mut header: Vec<String> = Vec::new();
        let mut seq: Vec<String> = Vec::new();
        let mut finalvec: Vec<FastaStruct> = Vec::new();
        for i in fileread.lines() {
            let line = i.expect("line not present");
            if line.starts_with("") {
                header.push(line.clone());
            }
            if !line.starts_with(">") {
                seq.push(line);
            }
        }
        for i in 0..header.len() {
            finalvec.push(FastaStruct {
                header: header[i].clone(),
                seq: seq[i].clone(),
            })
        }
        Ok(finalvec)
    }
}
