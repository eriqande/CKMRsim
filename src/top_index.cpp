#include <Rcpp.h>
#include <queue>

using namespace Rcpp;

// Top index stuff from Roman:  http://gallery.rcpp.org/articles/top-elements-from-vectors-using-priority-queue/

template <int RTYPE>
class IndexComparator {
public:
  typedef typename Rcpp::traits::storage_type<RTYPE>::type STORAGE;

  IndexComparator(const Vector<RTYPE>& data_) : data(data_.begin()) {}

  inline bool operator()(int i, int j) const {
    return data[i] > data[j] || (data[i] == data[j] && j > i);
  }

private:
  const STORAGE* data;
};

template <>
class IndexComparator<STRSXP> {
public:
  IndexComparator( const CharacterVector& data_ ) : data(data_.begin()) {}

  inline bool operator()(int i, int j) const {
    return (String)data[i] > (String)data[j] || (data[i] == data[j] && j > i );
  }

private:
  const Vector<STRSXP>::const_iterator data;
};

template <int RTYPE>
class IndexQueue {
public:
  typedef std::priority_queue<int, std::vector<int>, IndexComparator<RTYPE> > Queue;

  IndexQueue(const Vector<RTYPE>& data_) : comparator(data_), q(comparator), data(data_) {}

  inline operator IntegerVector() {
    int n = q.size();
    IntegerVector res(n);
    for( int i=0; i<n; i++) {
      // +1 for 1-based R indexing
      res[i] = q.top() + 1;
      q.pop();
    }
    return res;
  }
  inline void input( int i) {
    // if( data[ q.top() ] < data[i] ) {
    if (comparator(i, q.top())) {
      q.pop();
      q.push(i);
    }
  }
  inline void pop() { q.pop(); }
  inline void push(int i) { q.push(i); }

private:
  IndexComparator<RTYPE> comparator;
  Queue q;
  const Vector<RTYPE>& data;
};


template <int RTYPE>
IntegerVector top_index(Vector<RTYPE> v, int n) {
  int size = v.size();

  // not interesting case. Less data than n
  if( size < n){
    return seq( 0, n-1 );
  }

  IndexQueue<RTYPE> q( v );
  for( int i=0; i<n; i++) q.push(i);
  for( int i=n; i<size; i++) q.input(i);
  return q;
}

// [[Rcpp::export]]
IntegerVector top_index(SEXP x, int n) {
  switch (TYPEOF(x)) {
  case INTSXP: return top_index<INTSXP>(x, n);
  case REALSXP: return top_index<REALSXP>(x, n);
  case STRSXP: return top_index<STRSXP>(x, n);
  default: stop("type not handled");
  }
  return IntegerVector() ; // not used
}
