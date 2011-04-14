/*****************************************************************************
  cuffToTrans.h

  (c) 2011 - Aaron Quinlan
  Quinlan Laboratory
  Center for Public Health Genomics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef CUFFTOTRANS_H
#define CUFFTOTRANS_H

#include "bedFile.h"
#include "sequenceUtils.h"
#include "Fasta.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class CuffToTrans {

public:

    // constructor
    CuffToTrans(bool &useName, string &dbFile, string &bedFile, string &fastaOutFile,
        bool &useFasta, bool &useStrand);

    // destructor
    ~CuffToTrans(void);

    void ExtractDNA();
    void ReportDNA(const BED &bed, string &dna);


private:

    bool _useName;
    string _dbFile;
    string _bedFile;
    string _fastaOutFile;
    bool _useFasta;
    bool _useStrand;

    // instance of a bed file class.
    BedFile  *_bed;
    ostream *_faOut;
};

#endif
