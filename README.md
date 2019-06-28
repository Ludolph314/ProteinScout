# What ProteinScout is
My bachelor thesis project, a collection of Biopython and Python scripts meant for partial automation and streamlining of annotation process for unidentified protein sequences acquired through LC-MS/MS method of proteomic analysis. 

Main focus of ProteinScout is to accumulate information necessary to estimate the cellular localization of unknown proteins - specifically whether the protein is localized within the nucleus and whether it has a structural function in plant chromosomes. This is achieved by using the protein localization prediction tools (specified below) and also by BLASTing the unknown sequences against the SwissProt database, finding homologous proteins and checking their Gene Ontology properties for specific GO terms.

Input data should be two files: 
- a standard FASTA file with protein sequences and their IDs
- a TSV file with LC-MS/MS data processed by Peptide Shaker

ProteinScout is divided into three parts:
- Module 1 formats the input data for better access and further analysis, and performs a search for homologous proteins via protein BLAST in the SwissProt database. At this point, user has to manually input the processed file with protein sequences into the prediction tools web interface (WegoLoc, LOCALIZER, NucPred, CELLO2GO, Protein Cutter), and manually download their results. These results are used by ProteinScout as input files for Module 2 and 3.

- Module 2 merges partial results from the prediction tools together.

- Module 3 merges all available data into a final table, calculates SAF and NSAF values and estimates (based on the available data) whether a protein is localized within the cell nucleus or not.

# What ProteinScout isn't (yet)
A program which could fully automate the annotation process from beginning to end. Due to the fact that all used prediction tools only have a web interface, this would require me to incorporate browser automation and data mining, which would take some time (which I currently don't have) and experience in those fields (which I don't have either).

# What ProteinScout isn't (and most probably never will be)
An all-purpose universal tool for protein annotation. ProteinScout is very particular about the text structure of its input files and won't process data from any other prediction tool. It has been developed for a specific task according to the needs of the Department of Protein Biochemistry and Proteomics at the Center of the region Haná (CRH).

# Quick guide
- Start Module 1
- Specify the input files and output folder
- Specify BLAST parameters
- While the BLAST is running, use the output FASTA file to manually start prediction webtools
- Manually save the results of the prediction webtools to the specified folder. All results should be either in the same format as downloaded, or converted to a .txt file (for Protein Cutter output, it is currently necessary to manually copy a part of HTML code containing the table with results into a text file)
- (Optional) If the original FASTA file was split into several smaller ones, merge the partial results with Module 2
- Run Module 3 to finalize the annotation

More info about the ProteinScout is available in my bachelor thesis (Czech language only), or you can contact me via e-mail.

# Useful links
- http://cello.life.nctu.edu.tw/cello2go/
- http://localizer.csiro.au
- https://nucpred.bioinfo.se/nucpred/
- https://software.cr-hana.upol.cz/proteincutter/
- http://www.btool.org/WegoLoc
