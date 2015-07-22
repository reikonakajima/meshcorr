#include "GalaxyObjects.h"
#include "StringStuff.h"
using namespace std;



bool Compare_Source_RA(GalaxyObject* rhs, GalaxyObject* lhs) {
  return rhs->getRA() < lhs->getRA(); // sort in increasing order
}
bool Compare_Source_Dec(GalaxyObject* lhs, GalaxyObject* rhs) {
  return lhs->getDec() < rhs->getDec(); // sort in increasing order
}

template <class Tgalobj>
void 
GalaxyObjectList<Tgalobj>::sortByRA() {
  objPtrList.sort(Compare_Source_RA);
  return;
}

template <class Tgalobj>
void 
GalaxyObjectList<Tgalobj>::sortByDec() {
  objPtrList.sort(Compare_Source_Dec);
  return;
}

template <class Tgalobj>
GalaxyObjectList<Tgalobj>
GalaxyObjectList<Tgalobj>::cullByRA(double minra, double maxra) {
  GalaxyObjectList<Tgalobj> culledlist;
  this->sortByRA();
  typename list<Tgalobj*>::iterator i0 = searchRA(objPtrList.begin(), objPtrList.end(), minra);
  typename list<Tgalobj*>::iterator i1 = searchRA(i0, objPtrList.end(), maxra);
  culledlist.objPtrList.assign(i0, i1);
  return culledlist;
}

template <class Tgalobj>
GalaxyObjectList<Tgalobj>
GalaxyObjectList<Tgalobj>::cullByDec(double mindec, double maxdec) {
  GalaxyObjectList<Tgalobj> culledlist;
  this->sortByDec();
  typename list<Tgalobj*>::iterator i0 = searchDec(objPtrList.begin(), objPtrList.end(), mindec);
  typename list<Tgalobj*>::iterator i1 = searchDec(i0, objPtrList.end(), maxdec);
  culledlist.objPtrList.assign(i0, i1);
  return culledlist;
}

template <class Tgalobj>
typename list<Tgalobj*>::iterator 
GalaxyObjectList<Tgalobj>::searchRA(typename list<Tgalobj*>::iterator first,
				    typename list<Tgalobj*>::iterator last, 
				    const double ra) {
  typename list<Tgalobj*>::iterator i = first;
  while ( ra > (*i)->getRA() && i != last)
    ++i;
  return i;
}

template <class Tgalobj>
typename list<Tgalobj*>::iterator 
GalaxyObjectList<Tgalobj>::searchDec(typename list<Tgalobj*>::iterator first,
				     typename list<Tgalobj*>::iterator last, 
				     const double dec) {
  typename list<Tgalobj*>::iterator i = first;
  while ( dec > (*i)->getDec() && i != last)
    ++i;
  return i;
}


template <class Tgalobj>
void 
GalaxyObjectList<Tgalobj>::setBounds() {

  double minra, maxra, mindec, maxdec;

  this->sortByRA();
  typename list<Tgalobj*>::const_iterator i = objPtrList.begin();
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


template <class Tgalobj>
vector<Tgalobj*> 
GalaxyObjectList<Tgalobj>::getVectorForm() {

  vector<Tgalobj*> vectorform;
  vectorform.reserve(objPtrList.size());
  typename list<Tgalobj*>::const_iterator i = objPtrList.begin();
  for (; i != objPtrList.end(); ++i) {
    vectorform.push_back(*i);
  }
  return vectorform;
}


