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
      >> GalaxyObject::ra >> GalaxyObject::dec >> redshift
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
// GAMAObjectList() : constructor
//
GAMAObjectList::GAMAObjectList(istream& is)
{
  string buffer;
  while (getlineNoComment(is, buffer)) 
    GalaxyObjectList::objPtrList.push_back(new GAMAObject(buffer));
  return;
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


bool Compare_GAMA_Redshift(GalaxyObject* rhs, GalaxyObject* lhs) {
  GAMAObject* gama_ptr_rhs = static_cast<GAMAObject*>(rhs);
  GAMAObject* gama_ptr_lhs = static_cast<GAMAObject*>(lhs);
  return gama_ptr_rhs->getRedshift() < gama_ptr_lhs->getRedshift(); // sort in increasing order
}

void 
GAMAObjectList::sortByRedshift() {

  GalaxyObjectList::objPtrList.sort(Compare_GAMA_Redshift);  // check?

  return;
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

