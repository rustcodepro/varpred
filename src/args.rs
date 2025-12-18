use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "varpred",
    version = "1.0",
    about = "variant logistic classifier
       ************************************************
       Gaurav Sablok,
       Email: codeprog@icloud.com
      ************************************************"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Classify a entire population
    Classifier {
        /// input path to the gtf annotation
        gtffile: String,
        /// input path to the vcf annotation
        vcffile: String,
        /// input path to the expression values
        expressionvalues: String,
        /// genes fasta file
        genesfastainput: String,
        /// quality threshold
        qualitythresholdinput: String,
        /// expression prediction file
        expressionpredinput: String,
        /// threads for the analysis
        thread: String,
    },
}
