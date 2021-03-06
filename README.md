# Gene_Extract

The program here can extract genes from batches of whole genome sequences using a list of HMMs corresponding with the genes of interest.

For example, if you want to extract all of the luxI and luxR genes from a bunch of pseudomonas strains, what you do is put all of the pseudomonas genome fastas into one folder, and download HMMs which match conserved domains in luxI and luxR, respectively.

Then, you just tell the program what folder the genomes are stored in, how many HMMs are being used, and the paths to each HMM file.

A new folder will be created for each genome called <genome_name>_data_, and these folders will contain the extracted hit genes as well as some data files.

Compile using the command g++ -std=c++11 ming_pipeline ming_pipeline.cpp

Run using ./ming_pipeline genomedir numberofHMMs HMM1path HMM2path ... HMM1name HMM2name ...
