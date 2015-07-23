#include "GalaxyObjects.h"
#include "StringStuff.h"
using namespace std;


ComovingCoord::ComovingCoord(Cosmology cosmo, double ra, double dec, double redshift) {

  double comoving_r = cosmo.Dc(redshift);  // comoving distance

  double rarad = ra*DEGREE;
  double decrad = dec*DEGREE;

  double sindec = sin(decrad);
  x = comoving_r * cos(rarad) * sindec;
  y = comoving_r * sin(rarad) * sindec;
  z = comoving_r * cos(decrad);

  return;
}


bool Compare_Source_RA(GalaxyObject* rhs, GalaxyObject* lhs) {
  return rhs->getRA() < lhs->getRA(); // sort in increasing order
}
bool Compare_Source_Dec(GalaxyObject* lhs, GalaxyObject* rhs) {
  return lhs->getDec() < rhs->getDec(); // sort in increasing order
}

void 
GalaxyObjectList::sortByRA() {
  objPtrList.sort(Compare_Source_RA);
  return;
}

void 
GalaxyObjectList::sortByDec() {
  objPtrList.sort(Compare_Source_Dec);
  return;
}

GalaxyObjectList
GalaxyObjectList::cullByRA(double minra, double maxra) {
  GalaxyObjectList culledlist;
  this->sortByRA();
  list<GalaxyObject*>::iterator i0 = searchRA(objPtrList.begin(), objPtrList.end(), minra);
  list<GalaxyObject*>::iterator i1 = searchRA(i0, objPtrList.end(), maxra);
  culledlist.objPtrList.assign(i0, i1);
  return culledlist;
}

GalaxyObjectList
GalaxyObjectList::cullByDec(double mindec, double maxdec) {
  GalaxyObjectList culledlist;
  this->sortByDec();
  list<GalaxyObject*>::iterator i0 = searchDec(objPtrList.begin(), objPtrList.end(), mindec);
  list<GalaxyObject*>::iterator i1 = searchDec(i0, objPtrList.end(), maxdec);
  culledlist.objPtrList.assign(i0, i1);
  return culledlist;
}

list<GalaxyObject*>::iterator 
GalaxyObjectList::searchRA(list<GalaxyObject*>::iterator first, list<GalaxyObject*>::iterator last, 
			 const double ra) {
  list<GalaxyObject*>::iterator i = first;
  while ( ra > (*i)->getRA() && i != last)
    ++i;
  return i;
}

list<GalaxyObject*>::iterator 
GalaxyObjectList::searchDec(list<GalaxyObject*>::iterator first, list<GalaxyObject*>::iterator last, 
			  const double dec) {
  list<GalaxyObject*>::iterator i = first;
  while ( dec > (*i)->getDec() && i != last)
    ++i;
  return i;
}



bool Compare_Redshifts(GalaxyObject* rhs, GalaxyObject* lhs) {
  return rhs->getRedshift() < lhs->getRedshift(); // sort in increasing order
}

void
GalaxyObjectList::sortByRedshift() {

  objPtrList.sort(Compare_Redshifts);  // check?

  return;
}


void 
GalaxyObjectList::setBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  list<GalaxyObject*>::const_iterator i = objPtrList.begin();
  minra = (*i)->getRA();
  i = objPtrList.end();  --i;
  maxra = (*i)->getRA();

  this->sortByDec();
  i = objPtrList.begin();
  mindec = (*i)->getDec();
  i = objPtrList.end();  --i;
  maxdec = (*i)->getDec();

  bounds.setXMin(minra);
  bounds.setXMax(maxra);
  bounds.setYMin(mindec);
  bounds.setYMax(maxdec);

  return;
}


void
GalaxyObjectList::resetBounds() {
  bounds = Bounds<double>();
}


void
GalaxyObjectList::setComovingCoords(Cosmology c) {

  list<GalaxyObject*>::iterator i = objPtrList.begin();
  for (; i != objPtrList.end(); ++i) {
    coordPtrList.push_back(new ComovingCoord(c, (*i)->getRA(), (*i)->getDec(), (*i)->getRedshift()));
  }
  return;
}


vector<GalaxyObject*> 
GalaxyObjectList::getVectorForm() {

  vector<GalaxyObject*> vectorform;
  vectorform.reserve(objPtrList.size());
  list<GalaxyObject*>::const_iterator i = objPtrList.begin();
  for (; i != objPtrList.end(); ++i) {
    vectorform.push_back(*i);
  }
  return vectorform;
}


