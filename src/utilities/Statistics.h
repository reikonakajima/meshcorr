// $Id: Statistics.h,v 2.4 2007-04-23 18:48:57 reiko3 Exp $
// Some statistical operations on vectors
#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
using std::vector;
#include <algorithm>
#include <numeric>

namespace stats {

  template <class T>
  int getPctlIndex(vector<T>& pdist, float pct, float pnorm=0.) {
    if (pdist.empty()) return -1.;
    int i_marker = -1.;
    if (pnorm == 0.)
      pnorm = accumulate(pdist.begin(), pdist.end(), 0.);
    float pctlval = pnorm * pct;
    float lastpsum = 0.;
    float currentsum = 0.;
    for (int i=0; i < pdist.size(); ++i) {
      lastpsum = currentsum;
      currentsum += pdist[i];
      if (currentsum >= pctlval) {
	i_marker = i;
	break;
      }
    }
    return i_marker;
  }

  template <class T>
  int getLo2SigmaIndex(vector<T>& pdist, float pnorm=0.) {
    return getPctlIndex(pdist, 0.0275, pnorm);
  }
  template <class T>
  int getLo1SigmaIndex(vector<T>& pdist, float pnorm=0.) {
    return getPctlIndex(pdist, 0.15865, pnorm);
  }
  template <class T>
  int getMedianIndex(vector<T>& pdist, float pnorm=0.) {
    return getPctlIndex(pdist, 0.5, pnorm);
  }
  template <class T>
  int getHi1SigmaIndex(vector<T>& pdist, float pnorm=0.) {
    return getPctlIndex(pdist, 0.84135, pnorm);
  }
  template <class T>
  int getHi2SigmaIndex(vector<T>& pdist, float pnorm=0.) {
    return getPctlIndex(pdist, 0.9725, pnorm);
  }
  
  template <class T>
  T percentile(vector<T>& v, float pct) {
    if (v.empty()) return 0;
    int pctileIndex= static_cast<int> ((v.size()-1)*pct);
    if (pctileIndex>=v.size()) pctileIndex=v.size()-1;
    std::nth_element(v.begin(), v.begin()+pctileIndex, v.end());
    return v[pctileIndex];
  }

  template <class T>
  T median(vector<T>& v) {return percentile(v,0.5);}

  template <class T>
  T getLo2Sigma(vector<T>& v) {return percentile(v,0.0275);}
  template <class T>
  T getLo1Sigma(vector<T>& v) {return percentile(v,0.15865);}
  template <class T>
  T getHi1Sigma(vector<T>& v) {return percentile(v,0.84135);}
  template <class T>
  T getHi2Sigma(vector<T>& v) {return percentile(v,0.9725);}

  template <class T>
  T mean(vector<T>& v) {
    if (v.empty()) return 0;
    return std::accumulate(v.begin(),v.end(),0.) / v.size();
  }

  // Use median-16th pctile as sigma to clip a vector
  template <class T>
  void pctileClip(vector<T>& v, float nSigma) {
    T med = percentile(v,0.5);
    T lsig = percentile(v,0.16);
    T llim=med - nSigma*(med-lsig);
    T ulim=med + nSigma*(med-lsig);
    typename vector<T>::iterator dest=v.begin();
    for (typename vector<T>::iterator src=v.begin();
	 src != v.end();
	 ++src)
      if ( *src > llim && *src<ulim) *(dest++)=*src;
    v.erase(dest, v.end());
  }

  template <class T>
  T max(vector<T>& v) {
    T max = v[0];
    for (int i=1; i<v.size(); i++) if (v[i] > max) max = v[i];
    return max;
  }

  template <class T>
  T min(vector<T>& v) {
    T min = v[0];
    for (int i=1; i<v.size(); i++) if (v[i] < min) min = v[i];
    return min;
  }

  template <class T>
  T covariance(vector<T>& v1, vector<T>& v2) {
    int N = v1.size();
    if (N != v2.size()) {
      cerr << "covariance vector size incompatible";
      return 0;
    }
    if (N == 0) return 0;
    T cov = 0, mean1 = mean(v1), mean2 = mean(v2);
    for (int i=0; i<N; i++) cov += v1[i] * v2[i];
    cov = (cov / N) - (mean1 * mean2);
    return cov;
  }

  template <class T>
  T variance(vector<T>& v) {
    return covariance(v, v);
  }

  // sample covariance
  template <class T>
  T scovar(vector<T>& v1, vector<T>& v2) {
    int N = v1.size();
    if (N == 0) return 0;
    return covariance(v1,v2) * N / (N-1);
  }

  // sample variance
  template <class T>
  T svar(vector<T>& v) {
    return scovar(v, v);
  }

  template <class T>
  T rms(vector<T>& v) {
    int N = v.size();
    T rms = 0;
    for (int i=0; i<N; i++) rms += v[i] * v[i];
    rms /= N;
    return sqrt(rms);
  }

} // namespace
#endif 
