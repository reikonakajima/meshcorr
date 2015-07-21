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
    objPtrList.push_back(new GAMAObject(buffer));
  return;
}


bool Compare_GAMA_GAMAID(GAMAObject* rhs, GAMAObject* lhs) {
  return rhs->getGAMAId() < lhs->getGAMAId(); // sort in increasing order
}

void 
GAMAObjectList::sortByGAMAId() {
  objPtrList.sort(Compare_GAMA_GAMAID);  // not perfect, but min/max is ok
  return;

}


bool Compare_GAMA_Redshift(GAMAObject* rhs, GAMAObject* lhs) {
  return rhs->getRedshift() < lhs->getRedshift(); // sort in increasing order
}

void 
GAMAObjectList::sortByRedshift() {

  objPtrList.sort(Compare_GAMA_Redshift);  // check?

  return;
}


//
// getZBinSubsample() : get a subsample in redshift bins
//
GAMAObjectList
GAMAObjectList::getZBinSubsample(double zmin, double zmax) const {

  GAMAObjectList subsample; // initialize return list

  list<GAMAObject*>::const_iterator i = objPtrList.begin();
  for (; i != objPtrList.end(); ++i) {
    double z = (*i)->getRedshift();
    if (z >= zmin && z < zmax)
      subsample.objPtrList.push_back(*i);
  }

  return subsample;
}


