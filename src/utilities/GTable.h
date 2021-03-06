// Classes to represent lookup tables.
// A is the argument class, which must have ordering
// operations, and +-*/ for interpolation.
// D is the value class, which must have + and * operations
// to permit interpolation.
// 	$Id: Table.h,v 1.4 2004-03-03 14:26:58 garyb Exp $
#ifndef GTABLE_H
#define GTABLE_H

#include "Function1d.h"
#include "Std.h"

#include <vector>
using std::vector;
#include <algorithm>
#include <string>
#include <sstream>

// Exception classes:
class GTableError: public MyException {
 public:
  GTableError(const string &m=""): MyException("GTable Error: " +m) {}
};
class GTableOutOfRange: public GTableError {
public:
  GTableOutOfRange(): GTableError("Argument out of range") {}
};
class GTableReadError: public GTableError {
public:
  GTableReadError(const string &c): 
    GTableError("Data read error for line ->"+c) {}
};

// GTable element:
template<class V=double, class A=double>
class GTableEntry {
public:
  GTableEntry(A a, V v): arg(a), val(v) {}
  A arg;
  V val;
  bool operator==(const GTableEntry rhs) const {return arg==rhs.arg;}
  bool operator==(const A rhs) const {return arg==rhs;}
  bool operator!=(const GTableEntry rhs) const {return arg!=rhs.arg;}
  bool operator!=(const A rhs) const {return arg!=rhs;}
  bool operator>(const GTableEntry rhs) const {return arg>rhs.arg;}
  bool operator>(const A rhs) const {return arg>rhs;}
  bool operator<(const GTableEntry rhs) const {return arg<rhs.arg;}
  bool operator<(const A rhs) const {return arg<rhs;}
  bool operator<=(const GTableEntry rhs) const {return arg<=rhs.arg;}
  bool operator<=(const A rhs) const {return arg>=rhs;}
  bool operator>=(const GTableEntry rhs) const {return arg>=rhs.arg;}
  bool operator>=(const A rhs) const {return arg>=rhs;}
};

// The GTable itself:
template<class V=double, class A=double>
class GTable: public Function1d<V,A> {
public:
  enum interpolant {linear, spline, floor, ceil};
  //Construct empty table
  GTable(interpolant i=linear): v(), iType(i), isReady(false), y2() {} 
  //GTable from two arrays:
  GTable(const A* argvec, const V* valvec, int N, interpolant in=linear) ;
  GTable(const vector<A> &a, const vector<V> &v, interpolant in=linear) ;
  GTable(istream &is, interpolant in=linear): v(), iType(in), isReady(),
					     y2() {read(is);}
  void clear() {v.clear(); isReady=false;}
  void read(istream &is);
  void addEntry(const A a, const V v) ; //new element for table.
  V operator() (const A a) const ;	//lookup & interp. function value.
  V lookup(const A a) const ;	//interp, but exception if beyond bounds
  int size() const {return v.size();}	//size of table
  A argMin() const { setup(); if (v.size()>0) return v.front().arg;
   else throw GTableError("argMin for null GTable");}	//Smallest argument
  A argMax() const { setup(); if (v.size()>0) return v.back().arg;
   else throw GTableError("argMax for null GTable");} 	//Largest argument

  template <class T>
  void TransformVal(T &xfrm) {
    for (iter p=v.begin(); p!=v.end(); ++p)
      p->val = xfrm(p->arg, p->val);
    isReady=false; setup();
  }
  template <class T>
  void TransformArg(T &xfrm) {
    for (iter p=v.begin(); p!=v.end(); ++p)
      p->arg = xfrm(p->arg, p->val);
    isReady=false; setup();
  }

/**/void dump() const {setup(); for (citer p=v.begin(); p!=v.end(); ++p) 
		 cout << p->arg << " " << p->val << endl; }
private:
  typedef GTableEntry<V,A> Entry;
  typedef typename vector<Entry>::const_iterator citer;
  typedef typename vector<Entry>::iterator iter;
  interpolant iType;

  mutable vector<Entry> v;
  mutable int	lastIndex;	//Index for last lookup into table.
  mutable bool  isReady;	//Flag if table has been prepped.
  mutable vector<V> y2;		//vector of 2nd derivs for spline
  //get index to 1st element >= argument.  Can throw the exception here.
  iter upperIndex(const A a) const;
  void sortIt() const {std::sort(v.begin(), v.end());}
  void setup() const;	//Do any necessary preparation;
//Interpolate value btwn p & --p:
  V interpolate(const A a, const citer p) const; 
};
#endif
