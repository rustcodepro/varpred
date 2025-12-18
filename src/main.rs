mod args;
mod expr;
mod fasta;
mod gc;
mod machine;
use crate::args::CommandParse;
use crate::args::Commands;
use crate::machine::variantpred;
use clap::Parser;
use smartcore::linear::logistic_regression::LogisticRegression;
use smartcore::metrics::accuracy;

/*
Gaurav Sablok
codeprog@icloud.com
*/

fn main() {
    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::Classifier {
            gtffile,
            vcffile,
            expressionvalues,
            genesfastainput,
            qualitythresholdinput,
            expressionpredinput,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let commandnew = variantpred(
                    gtffile,
                    vcffile,
                    expressionvalues,
                    genesfastainput,
                    qualitythresholdinput,
                    expressionpredinput,
                )
                .unwrap();
                let machinefile =
                    LogisticRegression::fit(&commandnew.0, &commandnew.1, Default::default())
                        .unwrap();
                let machinepredict = machinefile.predict(&commandnew.2).unwrap();
                let accuracyvalue = accuracy(&commandnew.1, &machinepredict);
                println!("The accuracy of the predicted model is {:?}", accuracyvalue);
            });
        }
    }
}
