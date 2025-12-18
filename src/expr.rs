#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct VariantStore {
    pub name: String,
    pub pos: String,
    pub id: String,
    pub refalt: String,
    pub altallele: String,
    pub quality: String,
    pub filter: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct FinalVariant {
    pub name: String,
    pub pos: String,
    pub id: String,
    pub refalt: String,
    pub altallele: String,
    pub quality: String,
    pub filter: String,
    pub gene: String,
    pub geneseq: String,
    pub genelength: usize,
    pub genecontent: f32,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct MachineVariant {
    pub name: String,
    pub pos: String,
    pub id: String,
    pub refalt: char,
    pub altallele: String,
    pub quality: String,
    pub filter: String,
    pub gene: String,
    pub geneseq: String,
    pub genelength: usize,
    pub genecontent: f32,
    pub expression: f32,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct FastaRead {
    pub pathname: String,
}

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct FastaStruct {
    pub header: String,
    pub seq: String,
}
