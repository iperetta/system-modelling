/* code by: iperetta@ieee.org */
#ifndef PREP_MODEL_H_
#define PREP_MODEL_H_

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include "myfun.h"

using std::vector;
using std::string;
using std::ifstream;
using std::istringstream;
using std::cerr;
using std::endl;
using std::ostringstream;

namespace model {

	struct _Point {
		size_t dimension;
		vector<double> location;
		_Point(const size_t &dim = 1) : dimension(dim) {
			location.assign(dim,0.);
		}
		_Point(const vector<double> &coord) : dimension(coord.size()),
			location(coord) { }
		_Point(const double &x) : dimension(1) {
			location.assign(1,x);
		}
		_Point(const _Point &p) : dimension(p.dimension), 
			location(p.location) {}
		void set(const size_t &dim = 1) {
			dimension = dim;
			location.assign(dim,0.);
		}
		void operator= (const _Point &p) {
			if(dimension == p.dimension)
				location = p.location;
			else {
				dimension = p.dimension;
				location.clear();
				location.resize(dimension);
				location = p.location;
			}
		}
		void operator= (const vector<double> &v) {
			if(dimension == v.size())
				location = v;
			else {
				dimension = v.size();
				location.clear();
				location.resize(v.size());
				location = v;
			}
		}
		// Lexicographical order
		bool operator== (const _Point &p) const {
			if(dimension == p.dimension) {
				bool isequal = true;
				for(size_t i = 0; i < dimension; i++) {
					isequal &= (my::isLikelyZero(location.at(i)
						- p.location.at(i)));
					if(!isequal) break;
				}
				return isequal;
			} else { 
				cerr << "Points with different dimensions!" << endl; 
				return false;
			}
		}
		bool operator!= (const _Point &p) const {
			return !this->operator==(p);
		}
		bool operator< (const _Point &p) const {
			if(dimension == p.dimension) {
				bool islesser = location.at(0) < p.location.at(0);
				size_t i = 0;
				while(my::isLikelyZero(location.at(i)-p.location.at(i))
					&& i < dimension-1)
				{
					i++;
					if(location.at(i) < p.location.at(i)) 
					{
						islesser = true;
						break;
					}
				}
				return islesser;
			} else { 
				cerr << "Points with different dimensions!" << endl; 
				return false;
			}
		}
		bool operator<= (const _Point &p) const {
			return (this->operator==(p))? true : this->operator<(p);
		}
		bool operator> (const _Point &p) const {
			if(dimension == p.dimension) {
				bool isgreater = location.at(0) > p.location.at(0);;
				size_t i = 0;
				while(my::isLikelyZero(location.at(i)-p.location.at(i))
					&& i < dimension-1)
				{
					i++;
					if(location.at(i) > p.location.at(i)) 
					{
						isgreater = true;
						break;
					}
				}
				return isgreater;
			} else { 
				cerr << "Points with different dimensions!" << endl; 
				return false;
			}
		}
		bool operator>= (const _Point &p) const {
			return (this->operator==(p))? true : this->operator>(p);
		}
		void assign(const vector<double>::iterator &first, 
				const vector<double>::iterator &last) {
			location.assign(first, last);
			dimension = location.size();
		}
		void assign(const size_t &qty, const double &x) {
			location.assign(qty, x);
			dimension = qty;
		}
		double& at(const size_t &pos) { 
			return location.at(pos); 
		}
		const double& at(const size_t &pos) const {
			return location.at(pos);
		}
		string str() const {
			string crd = "( ";
			for(size_t i = 0; i < dimension; i++)
				crd += static_cast<ostringstream*>(&(ostringstream() 
					<< location.at(i)))->str() + 
					((i == dimension-1)? " " : ", ");
			crd += ")";
			return crd;
		}
	};

	struct _RawDataSet {
		size_t n_variables;
		vector<_Point> coordinates;
		vector<double> values;
		void loadfile(const string &filename) {
			string line;
			double dbl;
			size_t count_dbl = 0, check;
			vector<double> aux;
			vector<double>::iterator it;
			ifstream loadFile;
			// *** clean old RAW_DATA ***
			coordinates.clear();
			values.clear();
			// *** open file for reading ***
			loadFile.open(filename.c_str(), ifstream::in); 
			if(loadFile.is_open()) {
				getline(loadFile,line); // first line of file
				istringstream iss(line);
				while (iss >> dbl) {
					// update # of doubles in the first line
					aux.push_back(dbl);
					count_dbl++; 
				}
				// *** verify consistency ***
				if(count_dbl == 0 || count_dbl == 1) {
					cerr << "No useful RAW_DATA in " << filename << "!" 
						 << endl;
					loadFile.close(); // close file
					exit(1);
				}
				// *** go on ***
				while (!loadFile.eof()) { // next lines
					check = 0; 
					getline(loadFile,line);
					istringstream iss(line);
					while (iss >> dbl) {
						aux.push_back(dbl);
						check++;
					}
					// *** verify consistency ***
					if(check != count_dbl && !loadFile.eof()) {
						cerr << filename << ": file corrupted!" << endl;
						loadFile.close(); // close file
						exit(1);
					}
				}
				loadFile.close(); // close file
			} else {
				cerr << "Unable to open " << filename << "!" << endl;
				exit(1);
			}
			// *** organize RAW_DATA ***
			n_variables = count_dbl - 1;
			_Point buildcoord(n_variables);
			it = aux.begin();
			while(it != aux.end()) {
				buildcoord.assign(it, it + n_variables);
				coordinates.push_back(buildcoord);
				it += n_variables;
				values.push_back(*it);
				it++;
			}
		}
	} RAW_DATA;

	struct _SetOfDataPoints {
		vector<_Point> coordinates;
		vector<double> values;
		void clear() {
			coordinates.clear();
			values.clear();
		}
		double getInferiorBoundaryAt(const size_t &pos) {
			double minat = coordinates.at(0).at(pos);
			for(size_t i = 1; i < coordinates.size(); i++)
				if(minat > coordinates.at(i).at(pos))
					minat = coordinates.at(i).at(pos);
			return minat;
		}
		double getSuperiorBoundaryAt(const size_t &pos) {
			double maxat = coordinates.at(0).at(pos);
			for(size_t i = 1; i < coordinates.size(); i++)
				if(maxat < coordinates.at(i).at(pos))
					maxat = coordinates.at(i).at(pos);
			return maxat;
		}
		vector<vector<double> > consolidated() const {
			vector<vector<double> > cons;
			for(size_t i = 0; i < coordinates.size(); i++)
				cons.push_back(coordinates.at(i).location);
			return cons;
		}
		void show() const {
			for(size_t i = 0; i < coordinates.size(); i++)
				cout << "point " << i << " @" << coordinates.at(i).str() 
					<< " value of " << values.at(i) << endl;
		}
	};

	double dist(const _Point &a, const _Point &b) {
		return my::dist(a.location, b.location);
	}
	
	double hypervolume(const _Point &a, const _Point &b) {
		double hv = 1.;
		for(size_t i = 0; i < a.dimension; i++)
			hv *= abs(b.at(i) - a.at(i));
		return hv;
	}
	
	_Point centroid(const _SetOfDataPoints &auxset) {
		_Point cent(my::mean(auxset.consolidated()));
		return cent;
	}
	
	_Point min_coords(const vector<_Point> &set)
	{
		_Point minimal(set.at(0).dimension);
		for(size_t i = 0; i < minimal.dimension; i++)
		{
			minimal.at(i) = set.at(0).at(i);
			for(size_t j = 1; j < set.size(); j++)
				if(set.at(j).at(i) < minimal.at(i))
					minimal.at(i) = set.at(j).at(i);
		}
		return minimal;
	}
	
	_Point max_coords(const vector<_Point> &set)
	{
		_Point maximal(set.at(0).dimension);
		for(size_t i = 0; i < maximal.dimension; i++)
		{
			maximal.at(i) = set.at(0).at(i);
			for(size_t j = 1; j < set.size(); j++)
				if(set.at(j).at(i) > maximal.at(i))
					maximal.at(i) = set.at(j).at(i);
		}
		return maximal;
	}
	
	struct _MapDomainData {
		size_t n_variables;
		size_t points_inset;
		size_t up2next;
		vector<_Point> bound_inf;
		vector<_Point> bound_sup;
		vector<_SetOfDataPoints> domain;
		vector<_SetOfDataPoints> extras;
		vector<_SetOfDataPoints> reference;
		void init() {
			srand (time(NULL));
			vector<size_t> index;
			_SetOfDataPoints auxset;
			_SetOfDataPoints allpoints;
			vector<_Point> coord_data(RAW_DATA.coordinates);
			size_t npoints = coord_data.size();
			vector<double> values(RAW_DATA.values);
			// *** calculate number of points needed 2^(d+1) ***
			_Point binf(RAW_DATA.n_variables);
			_Point bsup(RAW_DATA.n_variables);
			cout << "RAW DATA contains " << RAW_DATA.n_variables << 
				" variables." << endl;
			up2next = 1;
			for(size_t i = 0; i < RAW_DATA.n_variables; i++) 
				up2next *= 2;
			points_inset = 3*up2next;
			// *** verify consistency of RAW_DATA ***
			if(RAW_DATA.coordinates.size() < points_inset) {
				cerr << "Not enough useful RAW_DATA to map!" << endl;
				exit(1);
			}
			domain.clear();
			extras.clear();
			reference.clear();
			bound_inf.clear();
			bound_sup.clear();
			// *** organize sets ***
			n_variables = RAW_DATA.n_variables;
			if(n_variables == 1) {
				alecJ::sort(coord_data, coord_data, index);
				alecJ::reorder(values, index, values);
				for(size_t i = 0; i < coord_data.size() - points_inset 
					+ 1; i += up2next)
				{
					// 1 variable: 4 points, 2 extrema
					allpoints.clear();
					for(size_t j = 0; j < 6; j++)
						allpoints.coordinates.push_back(coord_data.
							at(i+j));
					for(size_t j = 0; j < n_variables; j++)
					{
						binf.at(j) = allpoints.getInferiorBoundaryAt(j);
						bsup.at(j) = allpoints.getSuperiorBoundaryAt(j);
					}	
					bound_inf.push_back(binf);
					bound_sup.push_back(bsup);
					auxset.clear();
					auxset.coordinates.push_back(coord_data.at(i));
					auxset.coordinates.push_back(coord_data.at(i+5));
					auxset.values.push_back(values.at(i));
					auxset.values.push_back(values.at(i+5));
					domain.push_back(auxset);
					auxset.clear();
					auxset.coordinates.push_back(coord_data.at(i+1));
					auxset.coordinates.push_back(coord_data.at(i+3));
					auxset.values.push_back(values.at(i+1));
					auxset.values.push_back(values.at(i+3));
					extras.push_back(auxset);
					auxset.clear();
					auxset.coordinates.push_back(coord_data.at(i+2));
					auxset.coordinates.push_back(coord_data.at(i+4));
					auxset.values.push_back(values.at(i+2));
					auxset.values.push_back(values.at(i+4));
					reference.push_back(auxset);
				}
			} else { 
				binf = min_coords(coord_data);
				bsup = max_coords(coord_data);
				int nclustersperdim = (int) ceil(pow((double)npoints/
					(double)points_inset,1./(double)n_variables));
				if(nclustersperdim < 3)
					nclustersperdim = 3;
				if(pow(nclustersperdim,n_variables)*points_inset > 
					npoints)
				{
					cerr << "Low number of avaliable points: " << 
						npoints << "; please consider about " <<
						"increasing it to " << pow(nclustersperdim,
						n_variables)*points_inset << "."  << endl;
				}
				vector<vector<double> > markers(n_variables);
				for(size_t i = 0; i < n_variables; i++)
					markers.at(i) = my::linspace(binf.at(i),bsup.at(i),
						nclustersperdim);
				vector<_Point> pseudo_cluster;
				_Point aux_pt(n_variables);
				//-------------- counter -------------------------------
				size_t dimension = n_variables;
				size_t min_count = 0, 
					   max_count = nclustersperdim-1;
				vector<size_t> counter(dimension,min_count);
				while(counter.back() <= max_count) {
					//++++++++++++++++++++++++++++++++++++++++++++++++++
					for(size_t j = 0; j < n_variables; j++)
						aux_pt.at(j) = markers.at(j).at(counter.at(j));
					pseudo_cluster.push_back(aux_pt);
					//++++++++++++++++++++++++++++++++++++++++++++++++++
					counter.at(0)++;
					for(size_t c = 0; c < counter.size()-1; c++) {
						if(counter.at(c) > max_count) {
							counter.at(c) = min_count;
							counter.at(c+1)++;
						}
					}
				}
				//-------------- counter -------------------------------
				vector<double> distances;
				vector<_Point> ingroup;
				vector<double> ingroup_val;
				vector<_Point> outgroup;
				vector<double> outgroup_val;
				vector<_Point> hypercage;
				vector<bool> mychoice;
				for(size_t i = 0; i < pseudo_cluster.size(); i++)
				{
					ingroup.resize(points_inset);
					ingroup_val.resize(points_inset);
					distances.assign(npoints,0.);
					for(size_t j = 0; j < npoints; j++)
						distances.at(j) = dist(coord_data.at(j),
							pseudo_cluster.at(i));
					alecJ::sort(distances, distances, index);
					for(size_t j = 0; j < points_inset; j++)
					{
						ingroup.at(j) = coord_data.at(index.at(j));
						ingroup_val.at(j) = values.at(index.at(j));
					}
					binf = min_coords(ingroup); 
					bsup = max_coords(ingroup);
					bound_inf.push_back(binf);
					bound_sup.push_back(bsup);
					hypercage.clear();
					//-------------- counter ---------------------------
					dimension = n_variables;
					min_count = 0; 
					max_count = 1;
					counter.assign(dimension,min_count);
					while(counter.back() <= max_count) {
						//++++++++++++++++++++++++++++++++++++++++++++++
						for(size_t j = 0; j < n_variables; j++)
							aux_pt.at(j) = (counter.at(j) == 0) ? 
								binf.at(j) : bsup.at(j);
						hypercage.push_back(aux_pt);
						//++++++++++++++++++++++++++++++++++++++++++++++
						counter.at(0)++;
						for(size_t c = 0; c < counter.size()-1; c++) {
							if(counter.at(c) > max_count) {
								counter.at(c) = min_count;
								counter.at(c+1)++;
							}
						}
					}
					//-------------- counter ---------------------------
					auxset.clear(); // Building domain
					distances.assign(points_inset,0.);
					for(size_t j = 0; j < up2next; j++) 
					{
						outgroup.clear();
						outgroup_val.clear();
						for(size_t k = 0; k < ingroup.size(); k++)
							distances.at(k) = 0.5*(dist(ingroup.at(k),
								hypercage.at(j)) + pow(hypervolume(
								ingroup.at(k),hypercage.at(j)),1./
								n_variables)); //important metric, should be analysed afterwards
						alecJ::sort(distances, distances, index);
						int igs = ingroup.size()-1;
						for(int k = igs; k >= 0; k--)
						{
							if(k == (int) index.at(0))
							{
								auxset.coordinates.push_back(ingroup.
									back());
								auxset.values.push_back(ingroup_val.
									back());
								ingroup.pop_back();
								ingroup_val.pop_back();
								break;
							} 
							else
							{
								outgroup.push_back(ingroup.back());
								outgroup_val.push_back(ingroup_val.
									back());
								ingroup.pop_back();
								ingroup_val.pop_back();
							}
						}
						while(!outgroup.empty())
						{
							ingroup.push_back(outgroup.back());
							ingroup_val.push_back(outgroup_val.back());
							outgroup.pop_back();
							outgroup_val.pop_back();
						}
					}
					domain.push_back(auxset);
					auxset.clear(); // Building extra
					size_t count = 0, walk = 0;
					mychoice.assign(ingroup.size(),false);
					while(count < ingroup.size()/2)
					{
						if(((double)rand()/RAND_MAX) < 0.5)
						{
							mychoice.at(walk) = true;
							count = 0;
							for(size_t k = 0; k < mychoice.size(); k++)
								if(mychoice.at(k)) count++;
						}
						walk++;
						if(walk == mychoice.size()) walk = 0;
					}
					for(size_t k = 0; k < mychoice.size(); k++)
					{
						if(mychoice.at(k)) 
						{
							auxset.coordinates.push_back(ingroup.at(k));
							auxset.values.push_back(ingroup_val.at(k));
						}
					}
					extras.push_back(auxset);
					auxset.clear(); // Building reference
					for(size_t k = 0; k < mychoice.size(); k++)
					{
						if(!mychoice.at(k)) 
						{
							auxset.coordinates.push_back(ingroup.at(k));
							auxset.values.push_back(ingroup_val.at(k));
						}
					}
					reference.push_back(auxset);
				}
			}
		}
	} DATA;

}
#endif
