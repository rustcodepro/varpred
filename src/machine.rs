use crate::expr::FastaRead;
use crate::expr::FinalVariant;
use crate::expr::MachineVariant;
use crate::expr::VariantStore;
use crate::gc::gcontent;
use smartcore::linalg::basic::matrix::DenseMatrix;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
codeprog@icloud.com
*/

type GTFREAD = (String, usize, usize);
type DENSEMATRIX = (DenseMatrix<f32>, Vec<i32>, DenseMatrix<f32>);

pub fn variantpred(
    pathgtf: &str,
    pathvcf: &str,
    pathexpr: &str,
    genefastafile: &str,
    qualitythreshold: &str,
    expressionpred: &str,
) -> Result<DENSEMATRIX, Box<dyn Error>> {
    let pathgtfopen = File::open(pathgtf).expect("file not present");
    let pathvcf = File::open(pathvcf).expect("file not present");
    let pathexpression = File::open(pathexpr).expect("file not present");

    let pathvcfread = BufReader::new(pathvcf);
    let pathgtfread = BufReader::new(pathgtfopen);
    let pathexpressionread = BufReader::new(pathexpression);

    let mut vcfstruct: Vec<VariantStore> = Vec::new();
    for i in pathvcfread.lines() {
        let line = i.expect("line not present");
        let vcfvec = line.split("\t").collect::<Vec<_>>();
        vcfstruct.push(VariantStore {
            name: vcfvec[0].to_string(),
            pos: vcfvec[1].to_string(),
            id: vcfvec[2].to_string(),
            refalt: vcfvec[3].to_string(),
            altallele: vcfvec[4].to_string(),
            quality: vcfvec[5].to_string(),
            filter: vcfvec[6].to_string(),
        });
    }

    let mut filteredvariant: Vec<VariantStore> = Vec::new();
    for i in vcfstruct.iter() {
        if i.filter == "PASS" {
            filteredvariant.push(VariantStore {
                name: i.name.clone(),
                pos: i.pos.clone(),
                id: i.id.clone(),
                refalt: i.refalt.clone(),
                altallele: i.altallele.clone(),
                quality: i.quality.clone(),
                filter: i.filter.clone(),
            })
        }
    }

    let mut expressionvec: Vec<(String, f32)> = Vec::new();
    for i in pathexpressionread.lines() {
        let line = i.expect("file not present");
        let linevec = line.split("\t").collect::<Vec<_>>();
        expressionvec.push((linevec[0].to_string(), linevec[1].parse::<f32>().unwrap()));
    }

    let mut pathgtf: Vec<GTFREAD> = Vec::new();
    for i in pathgtfread.lines() {
        let line = i.expect("line not present");
        if line.starts_with("#") {
            continue;
        }
        if !line.starts_with("#") {
            let linevec = line.split("\t").collect::<Vec<_>>();
            if linevec[2] == "gene" {
                let genename = linevec[8].split(";").collect::<Vec<_>>()[0].replace("gene_id ", "");
                let inserttuple: (String, usize, usize) = (
                    genename,
                    linevec[3].parse::<usize>().unwrap(),
                    linevec[4].parse::<usize>().unwrap(),
                );
                pathgtf.push(inserttuple);
            }
        }
    }

    let pathfasta = FastaRead {
        pathname: genefastafile.to_string(),
    };
    let fastaunpack = pathfasta.fastaread().unwrap();

    let mut finalvariant: Vec<FinalVariant> = Vec::new();

    for i in filteredvariant.iter() {
        for val in pathgtf.iter() {
            for genename in fastaunpack.iter() {
                if i.pos.parse::<usize>().unwrap() <= val.1
                    || i.pos.parse::<usize>().unwrap() <= val.2 && val.0 == genename.header
                {
                    finalvariant.push(FinalVariant {
                        name: i.name.clone(),
                        pos: i.pos.clone(),
                        id: i.id.clone(),
                        refalt: i.refalt.clone(),
                        altallele: i.altallele.clone(),
                        quality: i.quality.clone(),
                        filter: i.filter.clone(),
                        gene: val.0.clone(),
                        geneseq: genename.seq.clone(),
                        genelength: val.2 - val.1,
                        genecontent: gcontent(&genename.seq).unwrap(),
                    });
                }
            }
        }
    }

    let mut machinelearning: Vec<MachineVariant> = Vec::new();
    for i in expressionvec.iter() {
        for val in finalvariant.iter() {
            if i.0 == val.gene {
                let valconvert = val.refalt.chars().collect::<Vec<char>>();
                machinelearning.push(MachineVariant {
                    name: val.name.clone(),
                    pos: val.pos.clone(),
                    id: val.pos.clone(),
                    refalt: valconvert[0],
                    altallele: val.altallele.clone(),
                    quality: val.quality.clone(),
                    filter: val.filter.clone(),
                    gene: val.gene.clone(),
                    geneseq: val.geneseq.clone(),
                    genelength: val.genelength,
                    genecontent: val.genecontent,
                    expression: i.1,
                })
            }
        }
    }

    let mut matrixquality: Vec<Vec<Vec<f32>>> = Vec::new();
    let mut class_labels: Vec<i32> = Vec::new();

    for i in machinelearning.iter() {
        if i.quality.parse::<usize>().unwrap() < qualitythreshold.parse::<usize>().unwrap() {
            let mut matrixinsert: Vec<Vec<f32>> = Vec::new();
            match i.refalt {
                'A' => matrixinsert.push(vec![
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                'T' => matrixinsert.push(vec![
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                'G' => matrixinsert.push(vec![
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                'C' => matrixinsert.push(vec![
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                _ => continue,
            }
            matrixquality.push(matrixinsert);
            class_labels.push(0);
        }
        if i.quality.parse::<usize>().unwrap() > qualitythreshold.parse::<usize>().unwrap() {
            let mut matrixinsert: Vec<Vec<f32>> = Vec::new();
            match i.refalt {
                'A' => matrixinsert.push(vec![
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                'T' => matrixinsert.push(vec![
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                'G' => matrixinsert.push(vec![
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                'C' => matrixinsert.push(vec![
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    i.quality.parse::<f32>().unwrap(),
                    i.altallele.split(",").collect::<Vec<_>>().len() as f32,
                    i.expression,
                ]),
                _ => continue,
            }
            matrixquality.push(matrixinsert);
            class_labels.push(1);
        }
    }

    let mut vectoraddexpr: Vec<Vec<f32>> = Vec::new();
    let expresspredopen = File::open(expressionpred).expect("file not found");
    let expressionopen_read = BufReader::new(expresspredopen);
    for i in expressionopen_read.lines() {
        let line = i.expect("line not present");
        vectoraddexpr.push(vec![line.parse::<f32>().unwrap()]);
    }

    let newmatrix = matrixquality.iter().flatten().cloned().collect::<Vec<_>>();

    let densematrix_quality = DenseMatrix::from_2d_vec(&newmatrix).unwrap();
    let predictionmatrix: DenseMatrix<f32> = DenseMatrix::from_2d_vec(&vectoraddexpr).unwrap();

    Ok((densematrix_quality, class_labels, predictionmatrix))
}
