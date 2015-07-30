//
// Class GAMAObject:
//   A catalog entry for an GAMA Object which keeps track of the various
//   redshift information of the object
//

#ifndef GAMAOBJECT_H
#define GAMAOBJECT_H
/**
 * @file GAMAObjects.h @brief Contains a class definition for GAMA objects.
 *
 * The GAMAObject include redshift, magnitude, stellar mass, absolute luminosity,
 * and various other galaxy information available from GAMA.
 */

#include "GalaxyObjects.h"
#include <iostream>
using std::istream;
using std::ostream;
#include <list>
using std::list;
#include <map>
using std::map;
using std::multimap;
#include <string>
using std::istringstream;
#include "Std.h"
#include <CCfits/CCfits>


class GAMAObjectError : public MyException {
 public:
 GAMAObjectError(const string &m="") : MyException("GAMAObjectError: " + m) {} 
};


/*/
// GAMAObject: a single entry for the GAMAObjectList class.

Edo's text file, split in KiDS fields (GAMA_[field]_V0.5_merged.asc)

 CATAID, GROUPIDA, RA, DEC, Z,
 ABSMAG_R, DELABSMAG_R, LOGMSTAR, DELLOGMSTAR, UMINUSR, 
 DELUMINUSR, RPETRO, LOGMOVERL, DELLOGMOVERL, LOGLWAGE,
 DELLOGLWAGE, METAL, DELMETAL, LOGTAU, DELLOGTAU,
 LOGMREMNANTS, DELLOGMREMNATS, RANKBCG, NFOF, ZMAX_19P8,
 ZMAX_19P4
/*/
class GAMAObject : public GalaxyObject {
 public:

  // constructor, reads a line
  GAMAObject(const string buffer);
  // simple constructor, covers only the basics
  GAMAObject(double _ra, double _dec, double _redshift) : GalaxyObject(_ra, _dec, _redshift) {}
  // destructor
  ~GAMAObject() {}

  string getGAMAId() const { return GAMAid; }
  string getGroupID() const { return groupID; }
  
  float getAbsMagR() const { return absMagR; }
  float getAbsMagRErr() const { return absMagRErr; }
  float getLogMStar() const { return logMStar; }
  float getLogMStarErr() const { return logMStarErr; }
  float getRPetro() const { return rPetro; }
  float getMetal() const { return metal; }
  float getMetalErr() const { return metalErr; }
  int   getRankBCG() const { return rankBCG; }   // CHECK !!!
  int   getNFOF() const { return nFOF; }   // CHECK !!!

 private:

  GAMAObject(const GAMAObject& rhs) {} // hide copy constructor

  // GAMA Object data
  string GAMAid;
  string groupID;

  float  absMagR, absMagRErr;
  float  logMStar, logMStarErr;
  float  rPetro;
  float  metal, metalErr;
  int    rankBCG;   // CHECK !!!
  int    nFOF;   // CHECK !!!

  float  uminusr, uminusrErr;
  float  logMOverL, logMOverLErr;
  float  loglwage, loglwageErr;
  float  logTau, logTauErr;
  float  logRemnants, logRemnantsErr;
  float  zmax19p8, zmax19p4;
};


// GAMAObject Ptr comparison class
class Compare_GAMARedshift {
 public:
  bool operator() (const GAMAObject* lhs, const GAMAObject* rhs) const {
    return lhs->getRedshift() < rhs->getRedshift();
  }
};


//
// GAMAObjectList: 
//
class GAMAObjectList : public GalaxyObjectList {
 public:
  // constructor, creates empty list
  GAMAObjectList() {};
  // constructor, reads input stream (assumes ASCII (Edo's) file)
  GAMAObjectList(istream& is);
  // constructor, read FITS file input  (currently only read ra/dec/z)
  GAMAObjectList(const string fits_filename);

  // read istream or FITS file.  returns size of newly read objects
  int read(istream& is);                 // Edo's ASCII file
  int read(const string fits_filename);  // FITS file

  // get a subsample in redshift bins
  GAMAObjectList getZBinSubsample(double zmin, double zmax) const ;
  // sort by spectroscopic redshifts
  void sortByRedshift();
  // sort by GAMA ID for matching purposes
  void sortByGAMAId();

 private:
  void _readFITSFile(const string fits_filename);
};



#endif // GAMAOBJECT_H
