// 	$Id: BinomFact.cpp,v 1.2 2002-08-16 15:49:44 garyb Exp $	

#include "BinomFact.h"

double fact(int i)
// return i!
{
  Assert(i>=0);
  static myvector<double> f(10);
  static bool first=true;
  if (first) {
    f[0] = f[1] = 1.;
    for(uint j=2;j<10;j++) f[j] = f[j-1]*(double)j;
    first = false;
  }
  if (i>=(int)f.size()) {
    for(int j=f.size();j<=i;j++)
      f.push_back(f[j-1]*(double)j);
    Assert(i==(int)f.size()-1);
  }
  Assert(i<(int)f.size());
  return f[i];
}

double sqrtfact(int i)
// return sqrt(i!)
{
  static myvector<double> f(10);
  static bool first=true;
  if (first) {
    f[0] = f[1] = 1.;
    for(uint j=2;j<10;j++) f[j] = f[j-1]*sqrt((double)j);
    first = false;
  }
  if (i>=(int)f.size())
    for(int j=f.size();j<=i;j++)
      f.push_back(f[j-1]*sqrt((double)j));
  Assert(i<(int)f.size());
  return f[i];
}

double binom(int i,int j)
// return iCj, i!/(j!(i-j)!)
{
  static myvector<myvector<double> > f(10);
  static bool first=true;
  if (first) {
    f[0] = myvector<double>(1,1.);
    f[1] = myvector<double>(2,1.);
    for(int i1=2;i1<10;i1++) {
      f[i1] = myvector<double>(i1+1);
      f[i1][0] = f[i1][i1] = 1.;
      for(int j1=1;j1<i1;j1++) f[i1][j1] = f[i1-1][j1-1] + f[i1-1][j1];
    }
    first = false;
  }
  if (j<0 || j>i) return 0.;
  if (i>=(int)f.size()) {
    for(int i1=f.size();i1<=i;i1++) {
      f.push_back(myvector<double>(i1+1,1.));
      for(int j1=1;j1<i1;j1++) f[i1][j1] = f[i1-1][j1-1] + f[i1-1][j1];
    }
    Assert(i==(int)f.size()-1);
  }
  Assert(i<(int)f.size());
  Assert(j<(int)f[i].size());
  return f[i][j];
}

double sqrtn(int i)
// return sqrt(i)
{
  static myvector<double> f(10);
  static bool first=true;
  if (first) {
    for(uint j=0;j<10;j++) f[j] = sqrt((double)j);
    first = false;
  }
  if (i>=(int)f.size())
    for(int j=f.size();j<=i;j++)
      f.push_back(sqrt((double)j));
  Assert(i<(int)f.size());
  return f[i];
}

