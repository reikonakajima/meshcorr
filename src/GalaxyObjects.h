//
// GalaxyObjects.h
//

#ifndef GALAXYOBJECTS_H
#define GALAXYOBJECTS_H
/**
 * @file GalaxyObjects.h @brief Contains a class definition for Galaxy objects.
 *
 * The GalaxyObject contain the minimum ra/dec information,
 * and is a base class for other galaxy survey datasets.
 */

#include <iostream>
#include <list>
#include <vector>
#include "Cosmology.h"
#include "Bounds.h"
#include "AstronomicalConstants.h"
using namespace cosmology;
using std::list;
using std::istringstream;


class GalaxyObjectsError : public MyException {
 public:
 GalaxyObjectsError(const string& m="") :
  MyException("GalaxyObjectsError: " +m) {}
};


//
// GalaxyObject class
//   Base class for LensObject and SourceObject
//
class GalaxyObject {
 public:
  GalaxyObject(double _ra, double _dec) {ra=_ra, dec=_dec;}
  GalaxyObject() {}

  // RA/Dec are stored in radians
  double getRA() const {return ra;}
  double getDec() const {return dec;}
  Position<double> getRADec() const {return Position<double>(ra,dec);}

 protected:
  double ra, dec;
};


class GalaxyObjectList {
 public:
  GalaxyObjectList() {}  // empty list
  void sortByRA();
  void sortByDec();
  GalaxyObjectList cullByRA(double minra, double maxra);
  GalaxyObjectList cullByDec(double mindec, double maxdec);
  
  int size() const { return objPtrList.size(); }
  Bounds<double> getBounds() { if (!bounds) setBounds(); return bounds;}
  void setBounds();   // defines bounds
  void resetBounds(); // undefines bounds
  vector<GalaxyObject*> getVectorForm();

 protected:
  list<GalaxyObject*> objPtrList;

 private:
  Bounds<double> bounds;
  list<GalaxyObject*>::iterator searchRA(list<GalaxyObject*>::iterator first, 
				    list<GalaxyObject*>::iterator last,
				    const double ra);
  list<GalaxyObject*>::iterator searchDec(list<GalaxyObject*>::iterator first, 
				     list<GalaxyObject*>::iterator last,
				     const double dec);

};


#endif // GALAXYOBJECTS_H
