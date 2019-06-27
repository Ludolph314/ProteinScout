# What ProteinScout is
My bachelor thesis project, a collection of scripts meant for partial automation and streamlining of annotation of unidentified protein sequences acquired through LC/MS-MS method of proteomic analysis. Input data should be two files: 
- a FASTA file with sequences and their IDs
- a TSV file with LC/MS-MS data processed by Peptide Shaker

ProteinScout is divided into three scripts:
- Module 1 formats the input data for better access and further analysis, and performs a search for homologous proteins via protein BLAST in the SwissProt database. At this point, user has to manually input the processed file with protein sequences into the prediction tools web interface (WegoLoc, LOCALIZER, NucPred, CELLO2GO, Protein Cutter), and manually download their results. These results are used by ProteinScout as input files for Module 2 and 3.

- Module 2 merges partial results from the prediction tools together.

- Module 3 merges all available data into a final table, calculates SAF and NSAF values and makes a decision (based on the available data) wheter a protein is localized in cell nucleus or not.

# What ProteinScout isn't (yet)
A program that could fully automate the annotation process from start to end. Due to the fact that the used prediction tools only have a web interface, this would require to incorporate browser automation and data mining, which would take some time (which I currently don't have) and experience in these fields (which I don't have either).

# What ProteinScout isn't (and most definitely will never be)
An all-purpose universal tool for protein annotation. ProteinScout is very particular about the input file text structure, and won't process data from any other prediction tool. It has been developed for a specific task according to the needs of the Department of Protein Biochemistry and Proteomics at the Center of the region Haná (CRH).

More info about the ProteinScout is available in my bachelor thesis (Czech language only), or you can contact me via e-mail.
