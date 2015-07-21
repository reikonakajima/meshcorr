#include "GalaxyObjects.h"
#include "StringStuff.h"
using namespace std;



bool Compare_Source_RA(GalaxyObject* rhs, GalaxyObject* lhs) {
  return rhs->getRA() < lhs->getRA(); // sort in increasing order
}
bool Compare_Source_Dec(GalaxyObject* lhs, GalaxyObject* rhs) {
  return lhs->getDec() < rhs->getDec(); // sort in increasing order
}

void 
GalaxyObjectList::sortByRA() {
  this->sort(Compare_Source_RA);
  return;
}

void 
GalaxyObjectList::sortByDec() {
  this->sort(Compare_Source_Dec);
  return;
}

GalaxyObjectList
GalaxyObjectList::cullByRA(double minra, double maxra) {
  GalaxyObjectList culledlist;
  this->sortByRA();
  GalaxyObjectList::iterator i0 = searchRA(this->begin(), this->end(), minra);
  GalaxyObjectList::iterator i1 = searchRA(i0, this->end(), maxra);
  culledlist.assign(i0, i1);
  return culledlist;
}

GalaxyObjectList
GalaxyObjectList::cullByDec(double mindec, double maxdec) {
  GalaxyObjectList culledlist;
  this->sortByDec();
  GalaxyObjectList::iterator i0 = searchDec(this->begin(), this->end(), mindec);
  GalaxyObjectList::iterator i1 = searchDec(i0, this->end(), maxdec);
  culledlist.assign(i0, i1);
  return culledlist;
}

GalaxyObjectList::iterator 
GalaxyObjectList::searchRA(GalaxyObjectList::iterator first, GalaxyObjectList::iterator last, 
			 const double ra) {
  GalaxyObjectList::iterator i = first;
  while ( ra > (*i)->getRA() && i != last)
    ++i;
  return i;
}

GalaxyObjectList::iterator 
GalaxyObjectList::searchDec(GalaxyObjectList::iterator first, GalaxyObjectList::iterator last, 
			  const double dec) {
  GalaxyObjectList::iterator i = first;
  while ( dec > (*i)->getDec() && i != last)
    ++i;
  return i;
}


void 
GalaxyObjectList::setBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  GalaxyObjectList::const_iterator i = this->begin();
  minra = (*i)->getRA();
  i = this->end();  --i;
  maxra = (*i)->getRA();

  this->sortByDec();
  i = this->begin();
  mindec = (*i)->getDec();
  i = this->end();  --i;
  maxdec = (*i)->getDec();

  bounds.setXMin(minra);
  bounds.setXMax(maxra);
  bounds.setYMin(mindec);
  bounds.setYMax(maxdec);

  return;
}


vector<GalaxyObject*> 
GalaxyObjectList::getVectorForm() {

  vector<GalaxyObject*> vectorform;
  vectorform.reserve(this->size());
  GalaxyObjectList::const_iterator i = this->begin();
  for (; i != this->end(); ++i) {
    vectorform.push_back(*i);
  }
  return vectorform;
}


