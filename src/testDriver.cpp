//
// meshDriver.cpp   : For testing the Mesh C++ module
//
#include <iostream>
#include "Std.h"
#include "StringStuff.h"
#include "Bounds.h"
#include "GAMAObjects.h"
using std::ostringstream;
using std::setw;
using std::setfill;
using std::fixed;
using std::setprecision;


const string usage =
  "\n"
  "testDriver: calculate tangential shear around a point (galaxy-galaxy lensing signal)\n"
  "\n"
  " usage: testDriver <GAMA_catalog>\n"
  "  GAMA_catalog:  lens catalog\n"
  "  \n"
  //  " output #1: file name:\" "+outfprefix+suffix+"\"\n"
  " stdin:  (none)\n"
  " stdout: (none)\n";


int
main(int argc, char* argv[]) {

  try {

    //
    // process arguments
    //
    if (argc < 2) {
      cerr << usage;
      exit(2);
    }

    // read in a series of gama files
    GAMAObjectList gama_list;
    for (int i_arg=1; i_arg<argc; i_arg++) {
      const string gama_filename = argv[i_arg];
    
      /// open gama file
      ifstream gamaf(gama_filename.c_str());
      if (!gamaf) 
	throw MyException("GAMA catalog file " + gama_filename + " not found");
      gama_list.read(gamaf);
    }


    //
    // diagnostic error messages
    //
    cerr << "=== testDriver ===" << endl;
    cerr << "GAMA catalog ...... " << argv[1] << endl;
    cerr << "     count ........ " << gama_list.size() << endl;
    cerr << "     bounds ....... " << gama_list.getBounds() << endl;

    if (gama_list.size() == 0) {
      cerr << "no gama objects, exiting" << endl;
      return(9);
    }


  } catch (MyException& m) {
    m.dump(cerr);
    exit(1);
  }

}
