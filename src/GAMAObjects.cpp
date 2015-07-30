#include "GAMAObjects.h"
#include "StringStuff.h"
#include "Bounds.h"
using namespace std;


//
// GAMAObject : constructor
//
GAMAObject::GAMAObject(const string buffer) {

  istringstream iss(buffer);

  // read basic object info
  iss >> GAMAid >> groupID
      >> GalaxyObject::ra >> GalaxyObject::dec >> GalaxyObject::redshift
      >> absMagR >> absMagRErr
      >> logMStar >> logMStarErr
      >> uminusr >> uminusrErr
      >> rPetro
      >> logMOverL >> logMOverLErr
      >> loglwage >> loglwageErr
      >> metal >> metalErr
      >> logTau >> logTauErr
      >> logRemnants >> logRemnantsErr
      >> rankBCG >> nFOF
      >> zmax19p8 >> zmax19p4;
  
  if (!iss)
    cerr << "#** invalid GAMAObject Entry:" << endl << buffer << endl;
  
  return;
}


//
// GAMAObjectList() : constructor for ASCII input (Edo's files)
//
GAMAObjectList::GAMAObjectList(istream& is)
{
  string buffer;
  while (getlineNoComment(is, buffer)) 
    GalaxyObjectList::objPtrList.push_back(new GAMAObject(buffer));
  return;
}


//
// _readFITSFile() : read in FITS file
//
void
GAMAObjectList::_readFITSFile(const string fits_filename)
{
  // FITS file pointer
  auto_ptr<CCfits::FITS> pInfile(0);

  // open the fits table and go to the right extension
  try {
    pInfile.reset(new CCfits::FITS(fits_filename, CCfits::Read));
  } catch(MyException& m) {
    m.dump(cerr);
  } catch (CCfits::FITS::CantOpen &fitsError) {
    throw MyException(fitsError.message());
  } catch (CCfits::FitsException &fitsError) {
    throw MyException(fitsError.message());
  }
  CCfits::ExtHDU& table = pInfile->extension(1);

  // read data from FITS table
  long n_row = table.rows();
  std::vector <double> ra;
  std::vector <double> dec;
  std::vector <double> redshift;
  try {
    table.column("ALPHA_J2000").read(ra, 0, n_row-1);
    table.column("DELTA_J2000").read(dec, 0, n_row-1);
    table.column("Z").read(redshift, 0, n_row-1);
  } catch (CCfits::Table::NoSuchColumn &fitsError) {
    throw MyException(fitsError.message());
  } catch (CCfits::Column::InvalidRowNumber &fitsError) {
    throw MyException(fitsError.message());
  } catch (CCfits::FitsException &fitsError) {
    throw MyException(fitsError.message());
  }

  // things that the random catalog doesn't have...
  std::vector <double> logmstar;
  std::vector <double> absmag_r;
  try {
    table.column("ABSMAG_R").read(absmag_r, 0, n_row-1);
    table.column("LOGMSTAR").read(logmstar, 0, n_row-1);
    // add on to the object pointer list if these exist, and return  (TEMPORARY HACK)
    for (long i=0; i<n_row; ++i) {
      GalaxyObjectList::objPtrList.push_back(new GAMAObject(ra[i], dec[i], redshift[i],
							    absmag_r[i], logmstar[i]));
    }
    return;
  } catch (CCfits::Table::NoSuchColumn &fitsError) {
    cerr << fitsError.message() << endl;
  }

  // no auxilary columns found, push_back each item and return
  for (long i=0; i<n_row; ++i) {
    GalaxyObjectList::objPtrList.push_back(new GAMAObject(ra[i], dec[i], redshift[i]));
  }

}


//
// GAMAObjectList() : constructor for FITS file input
//
GAMAObjectList::GAMAObjectList(const string fits_filename)
{
  _readFITSFile(fits_filename);
  return;
}



int
GAMAObjectList::read(istream& is)
{
  resetBounds();  // old ra/dec bounds is invalid

  int size_before = GalaxyObjectList::objPtrList.size();
  string buffer;
  while (getlineNoComment(is, buffer)) 
    GalaxyObjectList::objPtrList.push_back(new GAMAObject(buffer));
  return (GalaxyObjectList::objPtrList.size() - size_before);
}


int
GAMAObjectList::read(const string fits_filename)
{
  resetBounds();  // old ra/dec bounds is invalid

  int size_before = GalaxyObjectList::objPtrList.size();
  _readFITSFile(fits_filename);
  return (GalaxyObjectList::objPtrList.size() - size_before);
}


bool Compare_GAMAID(GalaxyObject* rhs, GalaxyObject* lhs) {
  GAMAObject* gama_ptr_rhs = static_cast<GAMAObject*>(rhs);
  GAMAObject* gama_ptr_lhs = static_cast<GAMAObject*>(lhs);
  return gama_ptr_rhs->getGAMAId() < gama_ptr_lhs->getGAMAId(); // sort in increasing order
}

void 
GAMAObjectList::sortByGAMAId() {
  GalaxyObjectList::objPtrList.sort(Compare_GAMAID);  // not perfect, but min/max is ok
  return;

}


//
// cullByLogMStar() : get a subsample in log10 stellar mass
//
GAMAObjectList
GAMAObjectList::cullByLogMStar(float min_log10_mstar, float max_log10_mstar) const {

  GAMAObjectList subsample;

  list<GalaxyObject*>::const_iterator i = GalaxyObjectList::objPtrList.begin();
  for (; i != GalaxyObjectList::objPtrList.end(); ++i) {
    GAMAObject* gama_ptr = static_cast<GAMAObject*>(*i);
    float logmstar = gama_ptr->getLogMStar();
    if (logmstar >= min_log10_mstar && logmstar < max_log10_mstar)
      subsample.GalaxyObjectList::objPtrList.push_back(*i);
  }
  return subsample;
}



//
// getZBinSubsample() : get a subsample in redshift bins
//
GAMAObjectList
GAMAObjectList::getZBinSubsample(double zmin, double zmax) const {

  GAMAObjectList subsample; // initialize return list

  list<GalaxyObject*>::const_iterator i = GalaxyObjectList::objPtrList.begin();
  for (; i != GalaxyObjectList::objPtrList.end(); ++i) {
    GAMAObject* gama_ptr = static_cast<GAMAObject*>(*i);
    double z = gama_ptr->getRedshift();
    if (z >= zmin && z < zmax)
      subsample.GalaxyObjectList::objPtrList.push_back(*i);
  }

  return subsample;
}


