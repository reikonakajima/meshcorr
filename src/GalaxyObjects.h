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
#include "Std.h"
#include "Cosmology.h"
#include "Bounds.h"
#include "AstronomicalConstants.h"
#include "Histogram.h"
#include "Mesh.h"
using namespace cosmology;
using std::list;
using std::istringstream;


class GalaxyObjectsError : public MyException {
 public:
 GalaxyObjectsError(const string& m="") :
  MyException("GalaxyObjectsError: " +m) {}
};


class ComovingCoord {
 public:
  ComovingCoord(Cosmology cosmo, double ra, double dec, double redshift);
  double getX() { return x; }
  double getY() { return y; }
  double getZ() { return z; }

 protected:
  double x, y, z;
};


//
// GalaxyObject class
//   Base class for LensObject and SourceObject
//
class GalaxyObject {
 public:
  GalaxyObject(double _ra, double _dec, double _redshift) {ra=_ra, dec=_dec, redshift=_redshift;}
  GalaxyObject() {}

  // RA/Dec are stored in radians
  double getRA() const {return ra;}
  double getDec() const {return dec;}
  Position<double> getRADec() const {return Position<double>(ra,dec);}

  double getRedshift() const { return redshift; }

  double getX() const {return coordPtr->getX();}
  double getY() const {return coordPtr->getY();}
  double getZ() const {return coordPtr->getZ();}

  void setCoordPtr(ComovingCoord* p) { coordPtr = p; return; }

 protected:
  double ra, dec;
  double redshift;
  mutable ComovingCoord* coordPtr;
};


class GalaxyObjectList {
 public:
  GalaxyObjectList() {}  // empty list
  void sortByRA();
  void sortByDec();
  GalaxyObjectList cullByRA(double minra, double maxra);
  GalaxyObjectList cullByDec(double mindec, double maxdec);
  
  void sortByRedshift();

  int size() const { return objPtrList.size(); }
  Bounds<double> getBounds() { if (!bounds) setBounds(); return bounds;}
  void setBounds();   // defines bounds
  void resetBounds(); // undefines bounds

  vector<GalaxyObject*> getVectorForm() const;
  void getXYZMinMax(double& xmin, double& xmax,
		    double& ymin, double& ymax,
		    double& zmin, double& zmax,
		    bool addEpsilon=false);

  void setComovingCoords(Cosmology cosmo);

  list<GalaxyObject*>::iterator objListBegin() { return objPtrList.begin(); }
  list<GalaxyObject*>::iterator objListEnd() { return objPtrList.end(); }

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


// Calculates Lande-Szalay correlation estimator (DD - DR - RD + RR) / RR
vector<double> LandeSzalay(GalaxyObjectList& data, GalaxyObjectList& randoms,
			   const HistogramLogBin& rbin, const double mesh_dx,
			   vector<double>& DD, vector<double>& DR,
			   vector<double>& RD, vector<double>& RR,
			   vector<double>& mean_r,
			   bool isAutoCorr=false);



#endif // GALAXYOBJECTS_H
