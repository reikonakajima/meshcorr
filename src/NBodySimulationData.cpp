#include "NBodySimulationData.h"
#include "StringStuff.h"
using std::istringstream;
using std::setw;

namespace nbody {
  

  float periodic (float x) {
    float tmp = x;
    if (tmp >=1.0) tmp = x-floor(x);
    if (tmp < 0.0) tmp = 1.0 -(-x-floor(-x));
    return tmp;
  }


  float diffperiodic (float x) {
    float tmp = x;
    do {
      if (tmp >= 0.5) tmp -= 1.0;
      if (tmp < -0.5) tmp += 1.0;
    } while (tmp >=0.5 || tmp < -0.5);
    return tmp;
  }


  GalaxyObject::GalaxyObject(const string buffer) {
    valid = false;
    istringstream iss(buffer);
    if ((iss >> pos[0] >> pos[1] >> pos[2]
	     >> vel[0] >> vel[1] >> vel[2] >> lum)) {
      valid = true;
      pos[0] = periodic(pos[0]);
      pos[1] = periodic(pos[1]);
      pos[2] = periodic(pos[2]);
    } 
  }


  GalaxyObject::GalaxyObject(const GalaxyObject& rhs) {
    valid = rhs.valid;
    pos[0] = rhs.pos[0];
    pos[1] = rhs.pos[1];
    pos[2] = rhs.pos[2];
    vel[0] = rhs.vel[0];
    vel[1] = rhs.vel[1];
    vel[2] = rhs.vel[2];
    lum    = rhs.lum;
  }


  GalaxyObject&
  GalaxyObject::operator=(const GalaxyObject& rhs) {
    valid = rhs.valid;
    pos[0] = rhs.pos[0];
    pos[1] = rhs.pos[1];
    pos[2] = rhs.pos[2];
    vel[0] = rhs.vel[0];
    vel[1] = rhs.vel[1];
    vel[2] = rhs.vel[2];
    lum    = rhs.lum;
  }


  void
  GalaxyObject::read(const string buffer) {
    valid = false;
    istringstream iss(buffer);
    if ((iss >> pos[0] >> pos[1] >> pos[2]
	     >> vel[0] >> vel[1] >> vel[2] >> lum)) {
      valid = true;
      pos[0] = periodic(pos[0]);
      pos[1] = periodic(pos[1]);
      pos[2] = periodic(pos[2]);
    } 
  }


  GalaxyList::GalaxyList(istream& is, int np_) {

    glist = new GalaxyObject[np_+1];
    
    np = 0;
    string buffer;
    while (getlineNoComment(is, buffer)) {

      glist[np].read(buffer);
      if (++np > (np_)) {
	cerr << "suggested gal count was " << np_ << endl;
	break;
      }
    }
    cerr << "galaxy count: " << np << endl;
  }
    

  void 
  GalaxyList::printPositions(ostream& os) const {

    for (int i=0; i<np; ++i) {
      os << setw(8) << glist[i].pos[0] << " " 
	 << setw(8) << glist[i].pos[1] << " " 
	 << setw(8) << glist[i].pos[2] << endl;
    }

    return;
  }


  DMParticleList::DMParticleList(char* fname) {

    int junk, hsize;
    FileHeader header;
    
    FILE *fp;
    fp = fopen(fname, "rb");

    if (!fp) 
      throw NBodySimError("cannot open DM particle file");

    int n = fread(&junk, 1, sizeof(int), fp);
    //cerr << "Format test (should be \"1\"): " << junk << endl;
    n = fread(&hsize, 1, sizeof(int), fp);
    n = fread(&header, 1, sizeof(FileHeader), fp);
    cerr << "Number of dm particles ..... " << header.npart << endl;
    //cerr << "Number of gas particles .... " << header.nsph << endl;
    //cerr << "Number of star particles ... " << header.nstar << endl;
    cerr << "Scale factor ............... " << header.aa << endl;
    //cerr << "Softening length ........... " << header.softlen << endl;

    pcount = header.npart;

    plist = (DMParticle*) malloc(sizeof(DMParticle)*pcount);
    size_t readsize;
    readsize = fread(plist, pcount, sizeof(DMParticle), fp);
  }


  GalaxyChainingMesh::GalaxyChainingMesh(GalaxyList gallist, int Nmesh) :
    Nm(Nmesh), np(gallist.size()) 
  {
    initialize();
    for (int nn = 0; nn<np; nn++) {                 // nn = galaxy index
      const int Nmesh = Nm*Nm*Nm;
      int mx = static_cast<int>(periodic(gallist[nn].pos[0]) * Nm);
      int my = static_cast<int>(periodic(gallist[nn].pos[1]) * Nm);
      int mz = static_cast<int>(periodic(gallist[nn].pos[2]) * Nm);
      int mm = (Nm*Nm*mx + Nm*my + mz) % (Nmesh);  // mm = mesh index
      if (Mhead[mm]<0) { // this cell is empty
	Mhead[mm] = Mtail[mm] = nn;
      } else {
	Mnext[Mtail[mm]] = nn;
	Mtail[mm] = nn;
      }
    }
  }


  GalaxyChainingMesh::~GalaxyChainingMesh() {
    if (Mhead) delete Mhead;
    if (Mnext) delete Mnext;
    if (Mtail) delete Mtail;
  }


  void 
  GalaxyChainingMesh::initialize() {

    int Nmesh = Nm*Nm*Nm;
    cerr << "Nmesh: " << Nmesh << endl;
    Mhead = new int[Nmesh];
    Mtail = new int[Nmesh];
    Mnext = new int[np];   
    for (int nn = 0; nn < Nmesh; nn++)
      Mhead[nn] = Mtail[nn] = -1;
    for (int ii = 0; ii < np; ii++)
      Mnext[ii] = -1;
    return;
  }


  void 
  GalaxyChainingMesh::printCounts(ostream& os) const {
    const int Nmesh = Nm*Nm*Nm;
    int total = 0;
    for (int nn = 0; nn < Nmesh; ++nn) {
      int count = 0;
      int nextgal = Mhead[nn];
      if (nextgal >= 0) {
	do {
	  count++;
	  nextgal = Mnext[nextgal];
	} while (nextgal >= 0);
      }
      os << setw(8) << nn << " "
	 << setw(10) << count << endl;
      total += count;
    }
    cerr << "total count in mesh: " << total << endl;
    return;
  }


  int  
  GalaxyChainingMesh::getNeighborList(const GalaxyList& glist, int gindex,
				      float rmin, float rmax, 
				      ParticleDistancePair* returnlist, 
				      int maxlistsize) const {
    // keep track of neighbor counts
    int count = 0;
    const int Ndim = 3;
    const float rminsq = rmin*rmin;
    const float rmaxsq = rmax*rmax;

    // find mesh position of current object
    int gx = static_cast<int>(Nm * glist[gindex].pos[0]);
    int gy = static_cast<int>(Nm * glist[gindex].pos[1]);
    int gz = static_cast<int>(Nm * glist[gindex].pos[2]);
    //cerr << mx << " " << my << " " << mz << endl;

    // set maximum mesh search radius
    int sr = static_cast<int>(Nm * rmax) + 1;
    int sc2 = (sr + 2) * (sr + 2);
    
    // identify all mesh within search radius
    for (int ix = -sr; ix <= sr; ix++)
      for (int iy = -sr; iy <= sr; iy++)
	for (int iz = -sr; iz <= sr; iz++) {
	  int ii;  
	  if (ix*ix + iy*iy + iz*iz < sc2) {
	    // calculate mesh id
	    ii = Nm*Nm*( (gx+ix+Nm)%Nm ) +
	            Nm*( (gy+iy+Nm)%Nm ) +
	               ( (gz+iz+Nm)%Nm );
	    // calculate distances to all galaxies in this mesh box
	    int p;  // particle index
	    if ( (p=Mhead[ii]) >= 0 ) {
	      //cerr << endl << "mesh id = " << ii << endl;
	      do {
		float distsq = 0;
		for (int j = 0; j < Ndim; ++j) {
		  float dr = 
		    diffperiodic( glist[p].pos[j] - glist[gindex].pos[j] );
		  distsq += dr*dr;
		}
		// include in neighbor list if within range
		if (distsq < rmaxsq && distsq > rminsq) {
		  returnlist[count].pid = p;
		  returnlist[count].distsq = distsq;
		  count++;
		  if (count >= maxlistsize) {
		    throw NBodySimError("List length exceeded max. size");
		  }
		}
	      } while ( (p=Mnext[p]) >= 0 );
	    }
	  }
	}

    return count;
  }


  void
  AccumulateDistanceCounts(ParticleDistancePair pdpair,
			   const HistogramLogBin& rbin, 
			   double hist[]) {
    if (pdpair.distsq < rbin.rsq(0))
      return;
    for (int i = 0; i < rbin.size(); ++i) {
      if (pdpair.distsq < rbin.rsq(i+1)) {
	++hist[i];
	/*
	cerr << "pairdist = " << sqrt(pdpair.distsq) << endl;
	cerr << rbin[i] << " < rbin < " << rbin[i+1] << endl;
	cerr << "i = " << i << " incremented" << endl;
	*/
	break;
      }
    }
    return;
  }


  void
  GalaxyChainingMesh::calculateAutoCorrelation(const HistogramLogBin& rbin,
					       const GalaxyList& glist,
					       double hist[]) const {
    double rmax = rbin.getRangeMax();
    double rmin = rbin.getRangeMin();
    if (rmax > 0.5) 
      throw NBodySimError("normalized rmax > 0.5 w/ periodic boundary cond.");
    
    int Nbin = rbin.size();
    for (int i = 0; i < Nbin; ++i) {
      hist[i] = 0;
    }
    
    // calculate summation
    int maxN = np/4;
    if (maxN < 1000) 
      maxN = 1000;
    for (int i = 0; i < np; ++i) {
      /*
      cerr << "gal id = " << i << ", "
	   << "pos = (" << glist[i].pos[0] << ","
	   << glist[i].pos[1] << ","
	   << glist[i].pos[2] << ")"
	   << endl;
      */
      ParticleDistancePair* templist = new ParticleDistancePair[maxN];
      int listsize = this->getNeighborList(glist,i,rmin,rmax,templist,maxN);
      for (int j = 0; j < listsize; ++j) {
	AccumulateDistanceCounts(templist[j],rbin,hist);
      }
      delete[] templist;
    }
    
    // calculate xi(r)
    double nbar = static_cast<double>(np);  // Volume = 1.0
    double xi[rbin.size()];
    for (int i = 0; i < rbin.size(); ++i) {
      //cerr << "hist[" << i << "] = " << hist[i] << endl;
      double r1cubed = rbin[i] * rbin[i] * rbin[i];
      double r2cubed = rbin[i+1] * rbin[i+1] * rbin[i+1];
      double dV = 4.*PI/3. * nbar * (r2cubed - r1cubed);
      hist[i] = hist[i] / dV / nbar - 1.0;
    }
      
    return;
  }




} // end namespace "nbody"

