/*****************************************************************************
  intersectBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "fjoin.h"
#include <queue>

bool leftOf(const BED &a, const BED &b);


bool BedIntersect::processHits(const BED &a, const vector<BED> &hits, bool printable) {

    // how many overlaps are there b/w the bed and the set of hits?
    int s, e, overlapBases;
	int  numOverlaps = 0;		
    bool hitsFound   = false;
	int aLength      = (a.end - a.start);   // the length of a in b.p.
    
	// loop through the hits and report those that meet the user's criteria
	vector<BED>::const_iterator h       = hits.begin();
	vector<BED>::const_iterator hitsEnd = hits.end();
	for (; h != hitsEnd; ++h) {
		s            = max(a.start, h->start);
		e            = min(a.end, h->end);
		overlapBases = (e - s);				// the number of overlapping bases b/w a and b

		// is there enough overlap relative to the user's request? (default ~ 1bp)
		if ( ( (float) overlapBases / (float) aLength ) >= _overlapFraction ) { 
			// Report the hit if the user doesn't care about reciprocal overlap between A and B.
			if (_reciprocal == false) {
				hitsFound = true;
				numOverlaps++;
				if (printable == true)
    				ReportOverlapDetail(overlapBases, a, *h, s, e);
			}
			// we require there to be sufficient __reciprocal__ overlap
			else {			
				int bLength    = (h->end - h->start);
				float bOverlap = ( (float) overlapBases / (float) bLength );
				if (bOverlap >= _overlapFraction) {
					hitsFound = true;
					numOverlaps++;
					if (printable == true)
        				ReportOverlapDetail(overlapBases, a, *h, s, e);
				}
			}
		}
	}
	// report the summary of the overlaps if requested.
	ReportOverlapSummary(a, numOverlaps);
	// were hits found for this BED feature?
	return hitsFound;
}

/*
	Constructor
*/
BedIntersect::BedIntersect(string bedAFile, string bedBFile, bool anyHit, 
						   bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
						   float overlapFraction, bool noHit, bool writeCount, bool forceStrand, 
						   bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput) {

	_bedAFile            = bedAFile;
	_bedBFile            = bedBFile;
	_anyHit              = anyHit;
	_noHit               = noHit;
	_writeA              = writeA;	
	_writeB              = writeB;
	_writeOverlap        = writeOverlap;
	_writeAllOverlap     = writeAllOverlap;		
	_writeCount          = writeCount;
	_overlapFraction     = overlapFraction;
	_forceStrand         = forceStrand;
	_reciprocal          = reciprocal;
    _obeySplits          = obeySplits;
	_bamInput            = bamInput;
	_bamOutput           = bamOutput;
	
	// create new BED file objects for A and B
	_bedA = new BedFile(bedAFile);
	_bedB = new BedFile(bedBFile);
	
	IntersectBed();
}


/*
	Destructor
*/
BedIntersect::~BedIntersect(void) {
}
 
 
bool leftOf(const BED &a, const BED &b) {
    return a.end <= b.start;
}

void BedIntersect::ReportOverlapDetail(const int &overlapBases, const BED &a, const BED &b,
									   const CHRPOS &s, const CHRPOS &e) {
	// default. simple intersection only
	if (_writeA == false && _writeB == false && _writeOverlap == false) {
		_bedA->reportBedRangeNewLine(a,s,e);
	}
	//  -wa -wbwrite the original A and B 
	else if (_writeA == true && _writeB == true) {
		_bedA->reportBedTab(a);
		_bedB->reportBedNewLine(b);
	}
	// -wa write just the original A
	else if (_writeA == true) {
		_bedA->reportBedNewLine(a);
	}
	// -wb write the intersected portion of A and the original B 
	else if (_writeB == true) {
		_bedA->reportBedRangeTab(a,s,e);
		_bedB->reportBedNewLine(b);
	}
	// -wo write the original A and B plus the no. of overlapping bases.
	else if (_writeOverlap == true) {
		_bedA->reportBedTab(a);
		_bedB->reportBedTab(b);
		printf("%d\n", overlapBases);
	}
}


void BedIntersect::ReportOverlapSummary(const BED &a, const int &numOverlapsFound) {
	// -u  just report the fact that there was >= 1 overlaps
	if (_anyHit && (numOverlapsFound >= 1)) {
		_bedA->reportBedNewLine(a);
	}
	// -c  report the total number of features overlapped in B
	else if (_writeCount) {
		_bedA->reportBedTab(a); 
		printf("%d\n", numOverlapsFound);
	}
	// -v  report iff there were no overlaps
	else if (_noHit && (numOverlapsFound == 0)) {
		_bedA->reportBedNewLine(a);
	}
	// -wao the user wants to force the reporting of 0 overlap
	else if (_writeAllOverlap && (numOverlapsFound == 0)) {
		_bedA->reportBedTab(a);
		_bedB->reportNullBedTab();
		printf("0\n");
	}
}



void BedIntersect::Scan(const BED &x, vector<BED> &wX, BedLineStatus xStatus, 
          const BED &y, vector<BED> &wY, BedLineStatus yStatus, bool fromA) {

    if (xStatus != BED_VALID) {
        return;
    }

    std::vector<BED>::iterator wYIter = wY.begin();
    std::vector<BED>::iterator wYEnd  = wY.end();

    for(; wYIter != wYEnd; ++wYIter) {
        if (leftOf(*wYIter, x) == true) {
            wY.erase(wYIter);
        }
        else if (overlaps(wYIter->start, wYIter->end, x.start, x.end) > 0) {
            if (fromA == true) {
                _bedA->reportBedTab(x);
                _bedB->reportBedNewLine(*wYIter);
            }
            else {
                _bedA->reportBedTab(*wYIter);
                _bedB->reportBedNewLine(x);                
            }
        }
    }
    if (leftOf(x,y) == false) {
        wX.push_back(x);
    }
}


void cullHits(const BED &a, vector<BED> &hits) {
    vector<BED>::iterator hitIter = hits.begin();
    vector<BED>::iterator hitEnd  = hits.end();
    for(; hitIter != hitEnd; ++hitIter) {
        if (leftOf(*hitIter, a) == true)
            hits.erase(hitIter);
    }
}


void BedIntersect::IntersectBed() {                                                                                                             
	
	int aLineNum = 0;
	int bLineNum = 0;
	BED a, b, nullBed;	
	BedLineStatus aStatus, bStatus;

    
    //  This is the fjoin algorithm.
    
    bool bSentinel = false;
    bool aSentinel = false;

    vector<BED> windowA, windowB;
     
     // open the files; get the first line from each
	_bedA->Open();
    _bedB->Open();
    aStatus = _bedA->GetNextBed(a, aLineNum);
    bStatus = _bedB->GetNextBed(b, bLineNum); 
    
    while (aSentinel == false || bSentinel == false) {
	    if (aStatus != BED_INVALID || bStatus != BED_INVALID) {
            if (a.start <= b.start) {
                Scan(a, windowA, aStatus, b, windowB, bStatus, true);
                aStatus = _bedA->GetNextBed(a, aLineNum);
    	    }
    	    else {
                Scan(b, windowB, bStatus, a, windowA, aStatus, false);
                bStatus = _bedB->GetNextBed(b, bLineNum);
    	    }
	    }
	    if (aStatus == BED_INVALID) {
            aSentinel = true; 
            a.start   = INT_MAX;
        }
	    if (bStatus == BED_INVALID) {
            bSentinel = true; 
            b.start   = INT_MAX;
        }
	}
	// close the files
	_bedA->Close();
	_bedB->Close();
	
	
	
	/* My modified algorithm
    vector<BED> hits;
    bool newBFound = false;
	_bedA->Open();
    _bedB->Open();
    
    aStatus = _bedA->GetNextBed(a, aLineNum);
    bStatus = _bedB->GetNextBed(b, bLineNum); 
    while (aStatus != BED_INVALID || bStatus != BED_INVALID) 
    {   
        newBFound = false;
        if (aStatus == BED_VALID) 
        {    
            // skip all of the B features that are left of A
            while ((leftOf(b,a) == true) && (bStatus == BED_VALID))
            {
                bStatus = _bedB->GetNextBed(b, bLineNum);
                newBFound = true;
            }
        
            // create a list of all the features overlapping A
            while ((overlaps(a.start, a.end, b.start, b.end) > 0) && (bStatus == BED_VALID))
            {
                hits.push_back(b);
                bStatus = _bedB->GetNextBed(b, bLineNum);
                newBFound = true;
            }
            // report the hits and move onto the next A feature.
            processHits(a, hits, true);
            // add the current b feature (the one breaking the above while loop) 
            // to hits and then cull
            // the hits that are no longer needed
            if (newBFound == true) {
                hits.push_back(b);
            }
            cullHits(a, hits);
        }
        // move on to the next A featrure
        aStatus = _bedA->GetNextBed(a, aLineNum);
        bStatus = _bedB->GetNextBed(b, bLineNum);
	}
	_bedA->Close();
	_bedB->Close();	
	*/
}


