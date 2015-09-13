/* code by: iperetta@ieee.org */
#ifndef MYFUN_H_
#define MYFUN_H_

#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <limits>
#include <vector>
#if defined(_WIN32) || defined(_WIN64)
 #define _USE_MATH_DEFINES
#endif
#include <cmath>
using namespace std;

namespace my {
	const double EPS = numeric_limits<double>::epsilon();
	bool isLikelyZero(const double &a, const double &k = 1.) 
	{
		return (std::isnan(a))? false : (abs(a) <= k*EPS);
	}
	bool isinf(const double &number, const double 
		&ref = 4503599627370496.0) // reference 2^52 = 4503599627370496
	{
		return (number > ref);
	}
	template<typename T>
	string num2str(const T &number) {
		return static_cast<ostringstream*>
			(&(ostringstream() << number))->str();
	}
	template<typename T>
	T str2num(const string &text) {
		T number;
		stringstream(text) >> number;
		return number;
	}
	template<typename T>
	void show_numvec(const vector<T>& vec) {
		string ST = "Size [" + my::num2str<size_t>(vec.size())+ "] => ";
		for(size_t i = 0; i < vec.size(); i++)
			ST += my::num2str<T>(vec.at(i)) + 
				((i == vec.size()-1) ? '\0' : '\t');
		cout << ST << endl;
	}
	template<typename T>
	void show_vec(const vector<T>& vec) {
		cout << "Size [" + my::num2str<size_t>(vec.size()) + "] => ";
		for(size_t i = 0; i < vec.size(); i++)
			cout << vec.at(i) << ((i == vec.size()-1) ? '\0' : '\t');
		cout << endl;
	}
	template<typename T>
	bool is_in(T data, const vector<T> &vec) {
		bool chk = false;
		if(!vec.empty()) {
			size_t i = 0;
			while(i < vec.size() && !chk) {
				if(data == vec.at(i)) 
					chk = true;
				i++;
			}
		}
		return chk;			
	}
	template<typename T>
	size_t find_in(T data, const vector<T> &vec) {
		size_t N = vec.size(), pos = N;
		size_t i = 0; 
		while(i < vec.size() && pos == N) {
			if(data == vec.at(i)) pos = i;
			i++;
		}
		if(pos == N)
			cerr << "[Warning] Data not found in vector!" << endl;
		return pos;			
	}
	template<typename T>
	vector<T> catsortvec(const vector<T> &veca, const vector<T> &vecb) {
		vector<T> idx,
				  idxa(veca),
			 	  idxb(vecb);
		typename vector<T>::iterator ia = idxa.begin(),
									 ib = idxb.begin();
		while(ia < idxa.end() || ib < idxb.end()) {
			if((ia < idxa.end()) && (ib < idxb.end())) {
				if(*ia < *ib) {
					idx.push_back(*ia); 
					ia++;
				} else {
					idx.push_back(*ib); 
					ib++;
				}
			}
			if((ia < idxa.end()) && !(ib < idxb.end())) {
				idx.push_back(*ia); 
				ia++;
			}
			if(!(ia < idxa.end()) && (ib < idxb.end())) {
				idx.push_back(*ib); 
				ib++;
			}
		}
		return idx;
	}
	double Sum(const vector<double> &vec)
	{
		double s = 0.;
		for(vector<double>::const_iterator i = vec.begin();
			i != vec.end(); i++)
				s += *i;
		return s;
	}
	double SumAbs(const vector<double> &vec)
	{
		double s = 0.;
		for(vector<double>::const_iterator i = vec.begin();
			i != vec.end(); i++)
				s += abs(*i);
		return s;
	}
	double Prod(const vector<double> &vec)
	{
		double p = 1.;
		for(vector<double>::const_iterator i = vec.begin();
			i != vec.end(); i++)
				p *= *i;
		return p;
	}
	double combination(const double &n, const double &k) {
		return tgamma(n+1)/(tgamma(k+1)*tgamma(n-k+1));
	}
	double pochhammer(const double &a, const double &n) {
		if(n < 0. && (my::isLikelyZero(n + a) || -n > a))
			return NAN;
		if(my::isLikelyZero(n)) // n == 0
			return 1.;
		if(my::isLikelyZero(n - 1.)) // n == 1
			return a;
		if(my::isLikelyZero(n + 1.)) // n == -1
			return 1./(a - 1.);
		if(my::isLikelyZero(a)) { // a == 0
			if(n < 0.)
				return pow(-1,n)/tgamma(n + 1.);
			else
				return 0.;
		}
		if(my::isLikelyZero(a - 1.)) // a == 1
			return tgamma(n + 1.);
		if(my::isLikelyZero(a - 2.)) // a == 2
			return tgamma(n + 2.);
		if(n < 0.) // n < 0, 
			return 1./pochhammer(a + n, -n);
		if(a < 0.) { // as a == 0 is already covered
			if(n < -a || my::isLikelyZero(n + a)) 
				return pow(-1,n)*tgamma(1. - a)/tgamma(1. - a - n);
			else //	if(n > -a)
				return 0.;
		}
		return tgamma(a + n)/tgamma(a);
	}
	double dist(const double &a, const double &b) {
		return abs(b-a);
	}
	double dist(const vector<double> &a, const vector<double> &b) {
		if(a.size() != b.size()) {
			cerr << "[Error] In " << __FILE__ << " at line " << 
				__LINE__ << endl;
			string msg = "Vectors of different sizes";
			throw std::out_of_range(msg.c_str());
		}
		double sqdiff = 0.;
		for(size_t i = 0; i < a.size(); i++)
			sqdiff += pow(b.at(i) - a.at(i),2.);
		return sqrt(sqdiff);
	}
	vector<double> linspace(const double &from, const double &to, const 
		size_t &howmany)
	{
		vector<double> markers(howmany);
		if(howmany > 2)
			for(size_t i = 0; i < howmany; i++)
				markers.at(i) = from + i*(to-from)/((double)howmany-1.);
		return markers;
	}
	int which_min(const vector<double> &a) {
		if(!a.empty())
		{
			double mina = a.at(0);
			int idx = 0;
			for(size_t i = 1; i < a.size(); i++)
				if(mina > a.at(i))
				{
					mina = a.at(i);
					idx = (int)i;
				}
			return idx;
		}
		else return -1;
	}
	double min(const vector<double> &a) {
		if(!a.empty())
		{
			double mina = a.at(0);
			for(size_t i = 1; i < a.size(); i++)
				if(mina > a.at(i))
					mina = a.at(i);
			return mina;
		}
		else return NAN;
	}
	int which_max(const vector<double> &a) {
		if(!a.empty())
		{
			double maxa = a.at(0);
			int idx = 0;
			for(size_t i = 1; i < a.size(); i++)
				if(maxa < a.at(i))
				{
					maxa = a.at(i);
					idx = (int)i;
				}
			return idx;
		}
		else return -1;
	}
	double max(const vector<double> &a) {
		if(!a.empty())
		{
			double maxa = a.at(0);
			for(size_t i = 1; i < a.size(); i++)
				if(maxa < a.at(i))
					maxa = a.at(i);
			return maxa;
		}
		else return NAN;
	}
	double mean(const double &a, const double &b) {
		return (b+a)/2.;
	}
	double mean(const vector<double> &a) {
		double a_sum = 0.;
		for(size_t i = 0; i < a.size(); i++)
			a_sum += a.at(i);
		return a_sum/a.size();
	}
	double protected_mean(const vector<double> &a) {
		double a_sum = 0.;
		int N = a.size();
		for(size_t i = 0; i < a.size(); i++)
			if(std::isnan(a.at(i)) || isinf(a.at(i)))
				N--;
			else
				a_sum += a.at(i);
		return a_sum/N;
	}
	double variance(const vector<double> &a) {
		double mu = mean(a);
		double a_sum = 0.;
		for(size_t i = 0; i < a.size(); i++)
			a_sum += pow(a.at(i) - mu,2.);
		return a_sum/(a.size() - 1.);
	}
	double protected_variance(const vector<double> &a) {
		double mu = protected_mean(a);
		int N = a.size();
		double a_sum = 0.;
		for(size_t i = 0; i < a.size(); i++)
			if(std::isnan(a.at(i)) || isinf(a.at(i)))
				N--;
			else
				a_sum += pow(a.at(i) - mu,2.);
		return a_sum/(N - 1.);
	}
	double std(const vector<double> &a) {
		return sqrt(variance(a));
	}
	double protected_std(const vector<double> &a) {
		return sqrt(protected_variance(a));
	}
	vector<double> mean(const vector<double> &a, 
			const vector<double> &b) {
		if(a.size() != b.size()) {
			cerr << "[Error] In " << __FILE__ << " at line " << 
				__LINE__ << endl;
			string msg = "Vectors of different sizes";
			throw std::out_of_range(msg.c_str());
		}
		vector<double> a_sum;
		for(size_t i = 0; i < a.size(); i++)
			a_sum.push_back((b.at(i)+a.at(i))/2.);
		return a_sum;
	}
	vector<double> mean(const vector<vector<double> > &a) {
		for(size_t i = 1; i < a.size(); i++) {
			if(a.at(0).size() != a.at(i).size()) {
				cerr << "[Error] In " << __FILE__ << " at line " << 
					__LINE__ << endl;
				string msg = "Vectors of different sizes";
				throw std::out_of_range(msg.c_str());
			}
		}
		vector<double> a_sum(a.at(0).size(),0.);
		for(size_t i = 0; i < a.at(0).size(); i++) {
			for(size_t j = 0; j < a.size(); j++)
				a_sum.at(i) += a.at(j).at(i);
			a_sum.at(i) /= a.size();
		}
		return a_sum;
	}
	inline int factorial(int x) 
	{
	  return ((x == 0 || x == 1)? 1 : x * factorial(x - 1));
	}
	void saveNumVecAs(const vector<double> &a, const string &namefile)
	{
		ofstream myfile;
		myfile.open(namefile);
		if (myfile.is_open())
		{
			for(size_t i = 0; i < a.size(); i++)
				myfile << a.at(i) << endl;
		}
		else
		{
			cerr << "Unable to open " << namefile << " file!" << endl;
		}
		myfile.close();
	}
	template<typename T>
	class STACK {
	  private:
		vector<T> V;
	  public:
	    STACK() {};
		void flush()
		{
			V.clear();
		}
		void push(const T &A)
		{
			V.push_back(A);
		}
		T pop()
		{
			T A = V.back();
			V.pop_back();
			return A;
		}
		~STACK() {};
	};
	vector<vector<int> > integer_partition(const int &n, const int 
		&nvar_ = 0)
	{
		// Adapted from ZS1
		// FAST ALGORITHMS FOR GENERATING INTEGER PARTITIONS 
		// ANTOINE ZOGHBIU and IVAN STOJMENOVIC (1998)
		int nvar = (nvar_ == 0)? n : nvar_;
		vector<int> x(n,1), aux;
		vector<vector<int> > output, partitions;
		x.at(0) = n;
		int m = 1, h = 1, r, t;
		aux = x;
		aux.resize(m);
		aux.resize(nvar, 0);
		output.push_back(aux);
		while(x.at(0) != 1)
		{
			if(x.at(h-1) == 2)
			{
				m++;
				x.at(h-1) = 1;
				h--;
			}
			else
			{
				r = x.at(h-1) - 1;
				t = m - h + 1;
				x.at(h-1) = r;
				while(t >= r)
				{
					h++;
					x.at(h-1) = r;
					t -= r;
				}
				if(t == 0)
					m = h;
				else
				{
					m = h + 1;
					if(t > 1)
					{
						h++;
						x.at(h-1) = t;
					}
				}
			}
			aux = x;
			aux.resize(m);
			aux.resize(nvar, 0);
			output.push_back(aux);
		}
		size_t i, j;
		int sum;
		for(i = 0; i < output.size(); i++)
		{
			sum = 0;
			for(j = 0; j < output.at(i).size(); j++)
				sum += output.at(i).at(j);
			if(sum == n)
				partitions.push_back(output.at(i));
		}
		return partitions;
	}
}

/*
//-------------- counter ------------------------------------\

size_t dimension = n_vars;
size_t min_count = 0, 
	   max_count = up2degree;
vector<size_t> counter(dimension,min_count);
while(counter.back() <= max_count) {
	//+++++++++++++++++++++++++++++++++++++++++++++++++++\
	// HERE GOES YOUR CODE!!!!
	// VECTOR *counter* HAS THE COUNTINGS FROM  
	// *min_count* TO *max_count* THAT YOU NEED
	// EX.: [ i,j,k,...n ] OR [ X(i),X(j),X(k),...X(n)])		
	my::show_numvec<size_t>(counter);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++/
	counter.at(0)++;
	for(size_t i = 0; i < counter.size()-1; i++) {
		if(counter.at(i) > max_count) {
			counter.at(i) = min_count;
			counter.at(i+1)++;
		}
	}
}
//-------------- counter ------------------------------------/
*/

namespace alecJ {

	template <class T>
	void sort(
		vector<T> &unsorted,
		vector<T> &sorted,
		vector<size_t> &index_map);
		
	template <class T>
	void reverse_sort(
		vector<T> &unsorted,
		vector<T> &sorted,
		vector<size_t> &index_map);

	// Act like matlab's Y = X[I]
	// where I contains a vector of indices so that after,
	// Y[j] = X[I[j]] for index j
	// this implies that Y.size() == I.size()
	// X and Y are allowed to be the same reference
	template< class T >
	void reorder(
	  vector<T> & unordered, 
	  vector<size_t> const & index_map, 
	  vector<T> & ordered);

	////////////////////////////////////////////////////////////////////
	// Implementation
	////////////////////////////////////////////////////////////////////


	// Comparison struct used by sort
	// http://bytes.com/topic/c/answers/132045-sort-get-index
	
	template<class T> struct index_cmp 
	{
	  index_cmp(const T arr) : arr(arr) {}
	  bool operator()(const size_t a, const size_t b) const
	  { 
		return arr[a] < arr[b];
	  }
	  const T arr;
	};
	//sort ascend
	template <class T>
	void sort(
	  vector<T> & unsorted,
	  vector<T> & sorted,
	  vector<size_t> & index_map)
	{
	  // Original unsorted index map
	  index_map.resize(unsorted.size());
	  for(size_t i=0;i<unsorted.size();i++)
	  {
		index_map[i] = i;
	  }
	  // Sort the index map, using unsorted for comparison
	  std::sort(
		index_map.begin(), 
		index_map.end(), 
		index_cmp<vector<T>& >(unsorted));

	  sorted.resize(unsorted.size());
	  reorder(unsorted,index_map,sorted);
	}
	
	// =============== by IPERETTA =============================
	
	template<class T> struct rev_index_cmp 
	{
	  rev_index_cmp(const T arr) : arr(arr) {}
	  bool operator()(const size_t a, const size_t b) const
	  { 
		return arr[a] > arr[b];
	  }
	  const T arr;
	};
	// sort descend
	template <class T>
	void reverse_sort(
	  vector<T> & unsorted,
	  vector<T> & sorted,
	  vector<size_t> & index_map)
	{
	  // Original unsorted index map
	  index_map.resize(unsorted.size());
	  for(size_t i=0;i<unsorted.size();i++)
	  {
		index_map[i] = i;
	  }
	  // Sort the index map, using unsorted for comparison
	  std::sort(
		index_map.begin(), 
		index_map.end(), 
		rev_index_cmp<vector<T>& >(unsorted));

	  sorted.resize(unsorted.size());
	  reorder(unsorted,index_map,sorted);
	}
	
	// =============== by IPERETTA =============================
	
	// This implementation is O(n), but also uses O(n) extra memory
	template< class T >
	void reorder(
	  vector<T> & unordered, 
	  vector<size_t> const & index_map, 
	  vector<T> & ordered)
	{
	  // copy for the reorder according to index_map, because unsorted 
	  // may also be sorted
	  vector<T> copy = unordered;
	  ordered.resize(index_map.size());
	  for(int i = 0; i<index_map.size();i++)
	  {
		ordered[i] = copy[index_map[i]];
	  }
	}

}

#endif
