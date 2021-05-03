// @Author Peter Wengert
// @Email  pcw4943@g.rit.edu
// @Date   1/22/2020
// @Description:
// The purpose of this program is to implement a complete pipeline 
// for genome annotation and figure generation based on the paper
// "Whole genome sequencing and analysis reveal insights into the 
//  genetic structure, diversity and evolutionary relatedness of 
//  luxI and luxR homologs in bacteria belonging to the
//  Sphingomonadaceae family" (Gan et al 2015)

// Steps: 
// 1) Read in a genome. 
// 2) Using prodigal, generate protein and gene files.
// 3) Use HMMsearch on the protein file using HMMs for the genes to be found.
// 5) For each HMM hit, grab the corresponding protein sequence from the protein file
// 6) Use InterProscan to find domains in each of the hits. 
// 7) If all the domains for the protein are present, store it as a true homolog.
// 8) For each genome with a true homolog, split it into contigs.
// 9) use prokka on the contigs containing true homologs to make a gbk file
// 10) Use easyfig to build a figure from the gbk file and the subregions identified in the hmmscan

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <dirent.h>
//#include </usr/local/Cellar/libomp/9.0.1/include/omp.h>

using namespace std;

// Cleans up a line in a fasta file by removing asterisks 
// Asterisks are not accepted by InterProScan5
// @line the line to be cleaned 
// @return the line cleaned of unwanted characters
std::string lineClean(std::string line)
{
    line.erase(std::remove(line.begin(), line.end(), '*'), line.end());
    return line;
}

// Extracts a single protein, line by line, from a protein faa file, and writes it to an output fasta.
// @infile is a filestream which has reached the starting line of a desired protein
// @ofilepath is the path (and name) of the new output file
// @header is the first line (the header) of the new protein fasta file.
int writeProt(std::ifstream &infile, std::string ofilePath, std::string header)
{
    std::ofstream ofile;
    ofile.open(ofilePath);
    ofile << header << '\n';
    for(std::string line; getline(infile, line) && line[0] != '>';)
        ofile << lineClean(line) << '\n';
    ofile.close();
    return 0;
}

// Searches a protein faa file for a given protein using its protein ID (<contig#_protein#>)
// Writes the isolated protein out to a new file line by line using the writeProt() method
// @protID is the identifier of the protein being searched for (<contig#_protein#>)
// @protPath is the location of the protein file
// @protName is used to name the output files for a given protein
// @o_dir is the output directory the generated files will be put into
// @return the path of the file produced
std::string findProtByID(std::string protID, std::string protPath, std::string protName, std::string o_dir)
{
    std::ifstream infile; 
	infile.open(protPath);
    int searching = 1;
    std::stringstream path;
    path << o_dir + "/" + protName + "_" + protID + ".fasta";
	for(std::string line; searching && getline(infile, line);)
		if(line.substr(0, line.find(" ")) == ">" + protID)
            searching = writeProt(infile, path.str(), line);
	infile.close();
    return path.str();
}

// Generates fasta files for each of the identified homologs to a given protein
// @protIDs is the list of IDs for proteins which were determined to be homologs to an HMM 
// @protName is used to name the output files for a given protein
// @protPath is the location of the protein file
// @return a list of the file names for each output
std::vector<std::string> generateProtFiles(std::vector<std::string> protIDs, std::string prefix, std::string protName, std::string protPath)
{
    std::stringstream ss;
    ss << prefix << "_" << protName;
    std::string o_dir = ss.str();
    std::vector<std::string> paths;
    system(("mkdir -p " + o_dir).c_str());
    for(unsigned int i = 0; i < protIDs.size(); i++)
        paths.push_back(findProtByID(protIDs[i], protPath, protName, o_dir));
    return paths;
}

// Returns a list containing the IDs of all of the protein hits to a given HMM
// @hmmOutPath is the location of the HMM file to be parsed
std::vector<std::string> parseHMMoutput(std::string hmmOutPath)
{
	std::ifstream infile; 
	infile.open(hmmOutPath);
    std::vector<std::string> protIDs;
	for(std::string line; getline(infile, line);)
		if(line[0] != '#')
            protIDs.push_back(line.substr(0, line.find(" ")));
	infile.close();
	return protIDs;
}

// A function which runs the prodigal gene finder program on a genome
// Putative proteins and genes will be written out to two files
// @path is the path to the genome being looked at
// @prefix is the prefix of the gene and protein file names
// @return the name of the file containing the proteins
std::string prodigal(std::string path, std::string prefix)
{
    std::stringstream ss;
    std::stringstream ofname;
    ofname << prefix << ".proteins.faa";
    std::string ofstr = ofname.str();
    ss << "prodigal -i " << path << " -o " << prefix << ".genes -a " << ofstr;
    system(ss.str().c_str());
    return ofstr;
}

// A function which scans a protein file for particular proteins using a list of HMMs
// HMM hits are written out to files named based on which HMM was used in the search
// ex: HMM1.hmm hits go to Protein1.txt, HMM2 hits go to Protein2.txt and so on
// @HMMvec is a vector containing the paths to several HMM files to be used
// @names is a vector containing the names of each protein being searched for
// @prefix is the name of the file containing the predicted protein sequences
// @return an array containing the names of the created files
std::vector<std::string> HMMsearch(std::vector<std::string> HMMvec, std::vector<std::string> names, int numHMMs, std::string prefix)
{
    std::vector<std::string> ofstrings;
    for(int i = 0; i < numHMMs; i++)
    {
        std::stringstream ofname;
        std::stringstream ss;
        ofname << prefix << "_" << names[i] << ".txt";
        ofstrings.push_back(ofname.str());
        ss << "HMMsearch --tblout " << ofname.str() << " " << HMMvec.at(i) << " " << prefix << ".proteins.faa";
        system(ss.str().c_str());
    }
    return ofstrings;
}

// int parseIPRscanTSV(std::string fpath)
// {
//     return 0;
// }

// Runs an InterProScan on one file
// Produces a TSV file 
// @file is the fasta file containing the sequence to be scanned
/* void InterProScan(std::string file)
{
    std::stringstream ss;
    ss << "python ~/iprscan5.py --outformat tsv --email pcw4943@g.rit.edu --sequence " << file;
    system(ss.str().c_str());
}

// A method which runs an IPRScan on each file (no more than 30 scans run concurrently)
// @protFiles is a list of file paths containing putative HMM matches
void batchIPRScan(std::vector<std::string> protFiles)
{
    omp_set_num_threads(30);
    #pragma omp parallel for
    for(unsigned int i = 0; i < protFiles.size(); i++)
        InterProScan(protFiles.at(i));
} */

// @argv[1] is the path to the genome to be read
// @argv[2] is a prefix to be use to name the gene and protein files
// @argv[3] is the number of HMM files to be used
// @argv[4] through @argv[4+argv[3]] is the list of HMM files to be used
// @argv[5+argv[3]] through @argv[5+2*argv[3]] are names
int main(int argc, char *argv[])
{
    // This should catch a lot of input errors
    if(argc < 4)
    {
        cout << "Usage: <genome dir> <number of HMMs> <HMM 1 path> <HMM 2 path> ... <HMM 1 name> <HMM 2 name> ...\n";
        return 1;
    }

    char *p_genomeDir = argv[1];
    std::stringstream ss;
    ss << p_genomeDir;
    std::string genomeDir = ss.str();
    ss.str(std::string());
    
    // The number of HMMs can vary, so it is necessary to check how many are listed
    int numHMMs = atoi(argv[2]);
    std::vector<std::string> HMMvec (numHMMs);
    std::vector<std::string> names (numHMMs);
    std::copy(argv + 3, argv + numHMMs + 3, HMMvec.begin());
    std::copy(argv + numHMMs + 3, argv + numHMMs*2 + 3, names.begin());

    struct dirent *de;
    DIR *dr = opendir(p_genomeDir);
    if(dr == NULL)
    {
        std::cout << "Unable to open directory: " << genomeDir << '\n';
        return 0;
    }

    while((de = readdir(dr)) != NULL)
    {
        // The purpose of this is to put all the files for a given run into one folder
        std::string prefix = de->d_name;
        std::string genome = genomeDir + "/" + prefix;

        system(("mkdir -p " + genome + "Data/").c_str());
        prefix = genome + "Data/" + prefix;

        std::string proteins = prodigal(genome, prefix);
        std::vector<std::string> hitFiles = HMMsearch(HMMvec, names, numHMMs, prefix);
        std::vector<std::vector<std::string> > protIDs;
        std::vector<std::string> protFiles;
        for(unsigned int i = 0; i < hitFiles.size(); i++)
        {
            protIDs.push_back(parseHMMoutput(hitFiles[i]));
            std::vector<std::string> protFile = generateProtFiles(protIDs[i], prefix, names[i], proteins);
            protFiles.insert(protFiles.end(), protFile.begin(), protFile.end());
        }

        //batchIPRScan(protFiles);
    }
    return 0;
}

