# varpred

- a variant count approach to predict variant. 
- based on expression values, number of the predicted alt alleles and sequence content. 
- Why this: 
    1. Takes into account the variant associated gene sequence features for machine learning.
    2. Any ref allele can have multiple alt allele so convert the alt allele count into a number of the predicted alt allele for feature engineering. 
    3. Quality value as the input vector and the classification. 
    4. The ref allele is always a single allele so the frequency of the same. 
- This allows you to link the expression of the pathogenic experiment to the predicted variant and classify them to others as pathogenic or not. 
- Taking into account all the features including the sequence features ones. 

```
cargo build
```
```
variant logistic classifier
       ************************************************
       Gaurav Sablok,
       Email: codeprog@icloud.com
      ************************************************

Usage: varpred <COMMAND>

Commands:
  classifier  Classify a entire population
  help        Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
varpred on  main is  0.1.0 
❯ ./target/debug/varpred classifier -h
Classify a entire population

Usage: varpred classifier <GTFFILE> <VCFFILE> <EXPRESSIONVALUES> <GENESFASTAINPUT> <QUALITYTHRESHOLDINPUT> <EXPRESSIONPREDINPUT> <THREAD>

Arguments:
  <GTFFILE>                input path to the gtf annotation
  <VCFFILE>                input path to the vcf annotation
  <EXPRESSIONVALUES>       input path to the expression values
  <GENESFASTAINPUT>        genes fasta file
  <QUALITYTHRESHOLDINPUT>  quality threshold
  <EXPRESSIONPREDINPUT>    expression prediction file
  <THREAD>                 threads for the analysis

Options:
  -h, --help  Print help
```

Gaurav Sablok \
codeprog@icloud.com
