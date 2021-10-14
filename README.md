# PyBLAST

Tools for protein sequence data analytics and genetic analysis.

Useful Resources:
- BLAST Web Interace: https://www.uniprot.org/blast/
- NCBI Web Interface: https://www.ncbi.nlm.nih.gov/
- Phylogenetic Software: https://cme.h-its.org/exelixis/software.html
- Muscle Protein Sequence Alignment: https://www.ebi.ac.uk/Tools/msa/muscle/
    - Documentation: https://petrov.stanford.edu/software/src/muscle3.6/muscle3.6.html
    - Python CLI API Tool: https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/MUSCLE+Help+and+Documentation
    - Binary: http://www.drive5.com/muscle/downloads.htm
- NCBI Multiple Sequence Alignment (MSA) Viewer: https://www.ncbi.nlm.nih.gov/tools/msaviewer/
- Protein BLAST Sequence Alignment (blastp): https://blast.ncbi.nlm.nih.gov/Blast.cgi

Usefule Commands:
Running MUSCLE MSA binary:
```
$ ./muscle<version> -in <path_to_FASTA_faa> -out <aligned_FASTA_afa> -maxiters 2
```
Use `maxiters` when `.faa` FASTA file contains more than 50 sequences, only runs 2 iterations of the algorithm, default is 16