/*****************************************************************************
  cuffToTransMain.cpp

  (c) 2011 - Aaron Quinlan
  Quinlan Laboratory
  Center for Public Health Genomics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "cuffToTrans.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "cuffToTrans"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string fastaDbFile;
    string bedFile;

    // output files
    string fastaOutFile;

    // checks for existence of parameters
    bool haveFastaDb = false;
    bool haveBed = false;
    bool haveFastaOut = false;
    bool useNameOnly = false;
    bool useFasta = true;
    bool useStrand = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) ShowHelp();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-fi", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveFastaDb = true;
                fastaDbFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-fo", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveFastaOut = true;
                fastaOutFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-bed", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveBed = true;
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-name", 5, parameterLength)) {
            useNameOnly = true;
        }
        else if(PARAMETER_CHECK("-tab", 4, parameterLength)) {
            useFasta = false;
        }
        else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
            useStrand = true;
        }
        else {
            cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    if (!haveFastaDb || !haveFastaOut || !haveBed) {
        showHelp = true;
    }

    if (!showHelp) {

        CuffToTrans *b2f = new CuffToTrans(useNameOnly, fastaDbFile, bedFile, fastaOutFile, useFasta, useStrand);
        delete b2f;

        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Extract DNA sequences into a fasta file based on feature coordinates." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta> " << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-fi\tInput FASTA file" << endl;
    cerr << "\t-bed\tBED/GFF/VCF file of ranges to extract from -fi" << endl;
    cerr << "\t-fo\tOutput file (can be FASTA or TAB-delimited)" << endl;
    cerr << "\t-name\tUse the name field for the FASTA header" << endl;

    cerr << "\t-tab\tWrite output in TAB delimited format." << endl;
    cerr << "\t\t- Default is FASTA format." << endl << endl;

    cerr << "\t-s\tForce strandedness. If the feature occupies the antisense strand," << endl;
    cerr << "\t\tthe sequence will be reverse complemented." << endl;
    cerr << "\t\t- By default, strand information is ignored." << endl << endl;



    // end the program here
    exit(1);

}
