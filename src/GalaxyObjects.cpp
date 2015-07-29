#include "GalaxyObjects.h"
#include "StringStuff.h"
using namespace std;


ComovingCoord::ComovingCoord(Cosmology cosmo, double ra, double dec, double redshift) {

  double comoving_r = cosmo.Dc(redshift) * HubbleLengthMpc;  // comoving distance in Mpc/h (c/H0)

  double rarad = ra*DEGREE;    // convert to radian units
  double decrad = dec*DEGREE;
  double cosdec = cos(decrad);

  x = comoving_r * cos(rarad) * cosdec;
  y = comoving_r * sin(rarad) * cosdec;
  z = comoving_r * sin(decrad);

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

  if (objPtrList.size()==0) {
    return;
  }

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
    ComovingCoord* cc = new ComovingCoord(c, (*i)->getRA(), (*i)->getDec(), (*i)->getRedshift());
    (*i)->setCoordPtr(cc);
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


void
GalaxyObjectList::getXYZMinMax(double& xmin, double& xmax,
			       double& ymin, double& ymax,
			       double& zmin, double& zmax,
			       bool addEpsilon) {

  list<GalaxyObject*>::const_iterator i = objPtrList.begin();
  xmin = xmax = (*i)->getX();
  ymin = ymax = (*i)->getY();
  zmin = zmax = (*i)->getZ();
  ++i;
  double x, y, z;
  for (; i != objPtrList.end(); ++i) {
    x = (*i)->getX();
    y = (*i)->getY();
    z = (*i)->getZ();

    if (xmin > x)  xmin = x;
    if (xmax < x)  xmax = x;
    if (ymin > y)  ymin = y;
    if (ymax < y)  ymax = y;
    if (zmin > z)  zmin = z;
    if (zmax < z)  zmax = z;
  }
  if (addEpsilon) {
    const double tinyfrac = 1e-6;
    double epsilon = (xmax - xmin) * tinyfrac;
    xmin -= epsilon;    xmax += epsilon;
    epsilon = (ymax - ymin) * tinyfrac;
    ymin -= epsilon;    ymax += epsilon;
    epsilon = (zmax - zmin) * tinyfrac;
    zmin -= epsilon;    zmax += epsilon;
  }
}


//
// accumulates pair counts into hist[], according to rbin separation
//
void
_accumulatePairCounts (std::multimap<int,double> nbr_index,
		      const HistogramLogBin& rbin, vector<double>& hist){

    for (std::multimap<int,double>::iterator ii=nbr_index.begin(); ii!=nbr_index.end(); ii++) {
      if ((*ii).second < rbin.rsq(0)) continue;
      for (int i=0; i<rbin.size(); ++i) {
	if ((*ii).second < rbin.rsq(i+1)) {
	  ++hist[i];
	  break;
	}
      }
    }

  return;
}


//
// Calculate correlation of "obj_mesh galaxies" around "obj_list galaxies" at separations in rbin[]
//
// Used in calculating DD, DR, RD, and RR for the Lande-Szalay correlation function estimator.
// Assumes the use of "randoms" which has similar volume coverage as the "data";
// hence no knowledge of the volume is necessary, as we take ratios of the
//  (number-normalized) correlations and the volume cancels out for the quantities DD/RR, DR/RR, etc.
//
vector<double>
_calculateCorrelation(const HistogramLogBin& rbin,
		      vector<GalaxyObject*> obj_list, Mesh<GalaxyObject*> obj_mesh) {

  vector<double> corr(rbin.size(), 0.);  // initialize return value

  int n_bar_1 = obj_list.size();      // n_bar = mean obj density (modified by Vol fraction)
  int n_bar_2 = obj_mesh.obj_count();

  // for each object in obj_list, find nearest neighbor in obj_mesh
  for (int i = 0; i < obj_list.size(); ++i) {
    double x0 = obj_list[i]->getX();
    double y0 = obj_list[i]->getY();
    double z0 = obj_list[i]->getZ();
    multimap<int, double> nbr_index = obj_mesh.getNearMeshMap(x0, y0, z0,
							      rbin.getRangeMax(),
							      rbin.getRangeMin());
    // accumulate pair counts into corr[]
    _accumulatePairCounts(nbr_index, rbin, corr);
  }

  // normalize to the mean average density
  for (int i = 0; i < rbin.size(); ++i) {
    double r1cubed = rbin[i] * rbin[i] * rbin[i];
    double r2cubed = rbin[i+1] * rbin[i+1] * rbin[i+1];
    // dn = number of expected obj_mesh galaxies in a spherical shell (modified by Vol fraction)
    double dn = 4.*PI/3. * n_bar_2 * (r2cubed - r1cubed);
    corr[i] = corr[i] / dn / n_bar_1;
  }
  return corr;
}



vector<double>
LandeSzalay(GalaxyObjectList data, GalaxyObjectList random,
	    const HistogramLogBin& rbin, double mesh_dx,
	    vector<double>& DD, vector<double>& DR, vector<double>& RD, vector<double>& RR) {

  // prepare inputs
  bool periodic = false;
  vector<GalaxyObject*> data_vector = data.getVectorForm();
  vector<GalaxyObject*> random_vector = random.getVectorForm();
  cerr << "data size: "  << data.size() << endl;
  cerr << "random size: "  << random.size() << endl;
  // use a common mesh with the same dimension
  double xmin, xmax, ymin, ymax, zmin, zmax;
  bool addEpsilon = true;
  random.getXYZMinMax(xmin, xmax, ymin, ymax, zmin, zmax, addEpsilon); // randoms are denser
  cerr << xmin << " " << xmax << " "
       << ymin << " " << ymax << " "
       << zmin << " " << zmax << endl;
  cerr << "dx=dy=dz: " << mesh_dx << endl;
  Mesh<GalaxyObject*, double> data_mesh(mesh_dx, mesh_dx, mesh_dx, data_vector, periodic,
					xmin, xmax, ymin, ymax, zmin, zmax);
  Mesh<GalaxyObject*, double> random_mesh(mesh_dx, mesh_dx, mesh_dx, random_vector, periodic,
					  xmin, xmax, ymin, ymax, zmin, zmax);

  // prepare output
  vector<double> xi(rbin.size());
  DD = _calculateCorrelation(rbin, data_vector, data_mesh);
  DR = _calculateCorrelation(rbin, data_vector, random_mesh);
  RD = _calculateCorrelation(rbin, random_vector, data_mesh);
  RR = _calculateCorrelation(rbin, random_vector, random_mesh);

  // calculate the Lande-Szalay estimator
  for (int i = 0; i < rbin.size(); ++i) {
    xi[i] = (DD[i] - DR[i] - RD[i]) / RR[i] + 1.;
  }  

  return xi;
}
