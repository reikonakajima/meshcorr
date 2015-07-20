#ifndef NBODYSIMULATIONDATA_H
#define NBODYSIMULATIONDATA_H

#include <iostream>
#include "Std.h"
#include "Histogram.h"

namespace nbody{
  
  class NBodySimError : public MyException {
  public:
    NBodySimError(const string &m=""):
      MyException("NBodySimError: " + m) {}
  };

  // wraps "x" periodically into the range [0,1)
  float periodic(float x);   

  // wraps "x" periodically into the range [-0.5,0.5)
  float diffperiodic(float x); 


  class GalaxyObject {
  public:
    GalaxyObject() {}
    GalaxyObject(const string buffer);
    GalaxyObject(const GalaxyObject& rhs);
    GalaxyObject& operator=(const GalaxyObject& rhs);
    void read(const string buffer);
    float pos[3];  // in units of box size, contents should be [0,1)
    float vel[3];  // in ??? units
    float lum;
  private:
    bool valid;
  };


  class GalaxyList {
  public:
    GalaxyList(istream& is, const int np_);  //  note: (np_ >= galaxy count)
    //~GalaxyList() {}
    GalaxyObject& operator[](int gid) const { return glist[gid]; }
    int size() const { return np; }
    void printPositions(ostream& os) const;
  private:
    GalaxyObject* glist;
    int np;  // galaxy count 
  };
  

  class FileHeader {
  public:
    int   npart;          /* Total number of particles. */
    int   nsph;           /* Number of gas particles.   */
    int   nstar;          /* Number of star particles.  */
    float aa;             /* Scale factor. */
    float softlen;        /* Gravitational softening    */
  };


  class DMParticle {
  public:
    float pos[3];
  };


  class DMParticleList {
  public:
    DMParticleList(char* fname);
  private:
    int pcount;
    DMParticle* plist;
  };


  class ParticleDistancePair {
  public:
    int   pid;     // particle id
    float distsq;  // distance squared
  };


  class GalaxyChainingMesh {
  public:
    GalaxyChainingMesh(GalaxyList glist, int meshside);
    ~GalaxyChainingMesh();
    void printCounts(ostream& os) const;
    void calculateAutoCorrelation(const HistogramLogBin& rbin, 
				  const GalaxyList& glist,
				  double* hist) const;
    int  getNeighborList(const GalaxyList& glist, int gindex,
			 float rmin, float rmax, 
			 ParticleDistancePair* returnlist, 
			 int returnlistsize) const;
    //void printAutoCorrelation(ostream& os, float rmin, float rmax);
  private:
    void initialize();  // reserves space for Mhead, Mnext, Mtail
    int  np;    // number of particles (galaxies)
    int  Nm;    // number of mesh on a side 
    int* Mhead;
    int* Mnext;
    int* Mtail;
  };


  class ParticleChainingMesh {
  public:
    ParticleChainingMesh(DMParticleList plist, int meshside);
    ~ParticleChainingMesh();
    void printCounts(ostream& os) const;
    void calculateAutoCorrelation(const HistogramLogBin& rbin, 
				  const DMParticleList& plist,
				  double* hist) const;
    void calculateGalaxyCorrelation(const HistogramLogBin& rbin, 
				    const GalaxyList& glist,
				    double* hist) const;
    int  getNeighborList(const DMParticleList& plist, int gindex,
			 float rmin, float rmax, 
			 ParticleDistancePair* returnlist, 
			 int returnlistsize) const;
    //void printAutoCorrelation(ostream& os, float rmin, float rmax);
  private:
    void initialize();  // reserves space for Mhead, Mnext, Mtail
    int  np;    // number of particles (galaxies)
    int  Nm;    // number of mesh on a side 
    int* Mhead;
    int* Mnext;
    int* Mtail;
  };

}




#endif // ifndef NBODYSIMULATIONDATA_H
