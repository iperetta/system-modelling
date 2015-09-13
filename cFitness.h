/* code by: iperetta@ieee.org */
#ifndef GMFITNESS_H_
#define GMFITNESS_H_
// =====================================================================
#include <iostream>
#ifndef CRANDOMXS_INIT_
	#define CRANDOMXS_INIT_
	#if defined(__linux__) || defined(__APPLE__)
		#include <unistd.h>
	long PID = getpid();
	#elif defined(_WIN32) || defined(_WIN64)
		#include <windows.h>
		long PID = GetCurrentProcessId();
	#else 
		#error "Platform not supported"
	#endif
	#include <ctime>
	#include "cPRNGXS.h"
	static cPRNGXS_1024 PRNG(time(NULL)+PID);
#endif
#include "prepareData.h"
#include "cexpAAST.h"
#include "cMatrix.h"
using namespace std;

class myFitness
{
  private:
    size_t nvar;
    size_t ndomains;
    size_t N;
    size_t Q;
    unsigned long nmcpts;
    string namefile;
    model::_MapDomainData dataset;
    vector<vector<int> > pdeg_vec;
    vector<vector<int> > dord_vec;
    vector<vector<double > > xis; //per xi and per var
    vector<Matrix> matrixPn; //per xi
    vector<Matrix> matrixPqm; //per xi
    vector<Matrix> CoV_diff; //per domain
    vector<Matrix> U; //per domain
  public:
    myFitness() : nvar(0), ndomains(0), N(0), Q(0), nmcpts(0), 
		namefile("") {};
    myFitness(const string &namefile, const size_t &polydeg, const 
		size_t &difford, const int &powerof2);
	void init(const string &namefile_, const size_t &polydeg, const 
		size_t &difford, const int &powerof2);
	void saveAs(const Matrix &U, const vector<double> &a, const 
		vector<double> &b, const bool &append) const;
	size_t getNVars() const { return this->nvar; };
	vector<vector<int> > getPolDegs() const { return this->pdeg_vec; };
	vector<vector<int> > getDifOrds() const { return this->dord_vec; };
    double jacobiP(const int &n, const double &alpha, const double 
		&beta, const double &xi) const;
    double P(const int &n, const double &alpha, const double &beta, 
		const double &xi) const;
    double dP(const int &k, const int &n, const double &alpha, const 
		double &beta, const double &xi) const;
    double xi_mapto_x(const double &xi, const double &a, const double 
		&b)	const;
    double x_mapto_xi(const double &x, const double &a, const double &b) 
		const;
    double u_hat(const vector<double> &uk, const vector<double>	&x, 
		const vector<double> &a, const vector<double> &b ) const;
	double mcint(const vector<double> &f, const double &a, const double 
		&b) const;
	double mcint(const vector<double> &f, const vector<double> 	&a, 
		const vector<double> &b) const;
    double solve_for(const vector<AAST> &K_, const AAST &s_, const bool 
		&saveflag = false);
};

myFitness::myFitness(const string &namefile, const size_t &polydeg, 
	const size_t &difford, const int &powerof2) 
{	
	this->init(namefile, polydeg, difford, powerof2);   
}

void myFitness::init(const string &namefile_, const size_t &polydeg, 
	const size_t &difford, const int &powerof2)
{
	// Load DATA =======================================================
	model::RAW_DATA.loadfile(namefile_);
	this->namefile =  "TGE/TGE_" + my::num2str(PID) + "_" + 
		namefile_;
	model::DATA.init();
	this->dataset = model::DATA;
	this->nvar = this->dataset.n_variables;
	this->ndomains = this->dataset.domain.size();
	// Additional Information ==========================================
	this->pdeg_vec = getPowers(this->nvar, polydeg);
	this->N = pdeg_vec.size();
	this->dord_vec = getPowers(this->nvar, (difford > polydeg) ? 
		polydeg : difford);
	this->Q = dord_vec.size();
	// ** BASE FOR ALL ** ==============================================
	// Galerkin truncate expansion coeffcients
	this->U.resize(this->ndomains);
	// Generate random numbers to be used within interval [-1,1]
	this->nmcpts = 1L << powerof2;
	this->xis.resize(this->nmcpts);
	vector<double> auxis(this->nvar);
	for(vector<vector<double> >::iterator it = this->xis.begin(); it != 
			this->xis.end(); it++)
	{
		for(size_t v = 0; v < this->nvar; v++)
			auxis.at(v) = PRNG.real(-1.,1.);
		*it = auxis;
	}
	// Calculate matrices of polynomials
	Matrix mat_Pn(N,1);
	this->matrixPn.resize(this->nmcpts);
	Matrix mat_Pqm(Q,N);
	this->matrixPqm.resize(this->nmcpts);
	double prod_;
	for(unsigned long w = 0; w < this->nmcpts; w++)
	{
		for(size_t i = 0; i < this->N; i++)
		{
			prod_ = 1.;
			for(size_t v = 0; v < this->nvar; v++)
				prod_ *= P(this->pdeg_vec.at(i).at(v),0.,0.,this->
					xis.at(w).at(v));
			mat_Pn.at(i) = prod_;
			for(size_t k = 0; k < this->Q; k++)
			{
				prod_ = 1.;
				for(size_t v = 0; v < this->nvar; v++)
					prod_ *= dP(this->dord_vec.at(k).at(v),
						this->pdeg_vec.at(i).at(v),0.,0.,this->
						xis.at(w).at(v));
				mat_Pqm.at(k,i) = prod_;
			}
		}
		this->matrixPn.at(w) = mat_Pn;
		this->matrixPqm.at(w) = mat_Pqm;
	}
	// ** DOMAIN-RELATED ** ============================================
	// For every differential, for each domain
	Matrix CoVvec(1,this->Q);
	this->CoV_diff.resize(this->ndomains);
	for(size_t d = 0; d < this->ndomains; d++)
	{
		for(size_t k = 0; k < this->Q; k++)
		{
			// Vector for change of variables in derivatives
			prod_ = 1.;
			for(size_t v = 0; v < this->nvar; v++)
				prod_ *= pow(2./(this->dataset.bound_sup.at(d).at(v) - 
					this->dataset.bound_inf.at(d).at(v)), 
					(double)this->dord_vec.at(k).at(v));
			CoVvec.at(k) = prod_;
		}
		this->CoV_diff.at(d) = CoVvec;
	}
	// Note that the determinant of the Jacobian matrix is not necessary
}

void myFitness::saveAs(const Matrix &U, const vector<double> &a, 
	const vector<double> &b, const bool &append) const
{
	if(a.size() != this->nvar || b.size() != this->nvar)
		cerr << "Limits for saving information are not consistent!"
			<< endl;
	else
	{
		vector<double> numvec(U.as_vector());
		ofstream myfile;
		myfile.open(this->namefile,((append)? ios::app : ios::trunc));
		if (myfile.is_open())
		{
			for(size_t v = 0; v < this->nvar; v++)
			{
				myfile << a.at(v) << '\t';
				myfile << b.at(v) << '\t';
			}
			for(size_t i = 0; i < numvec.size(); i++)
				myfile << numvec.at(i) << '\t';
			myfile << endl;
		}
		else
			cerr << "Unable to open " << namefile << " file!" << endl;
		myfile.close();
	}
}

double myFitness::jacobiP(const int &n, const double &alpha, const 
	double &beta, const double &xi) const
{
	double expr; 
	if(n <= 0) {
		if (n == 0)
			expr = 1.;
		else
			expr = 0.;
	} else {
		expr = 0.;
		double coeff = my::combination(n + alpha, n);
		double sumcoeff = 0.;
		for(int k = 0; k <= n; k++) {
			sumcoeff = my::pochhammer(-(n), (double)k)*
				my::pochhammer(n+alpha+beta+1., 
				(double)k)/(my::pochhammer(alpha+1.,(double)k)*
				tgamma((double)k+1.));
			expr += sumcoeff*pow(.5*(1. - xi), (double)k);
		}
		expr *= coeff;
	}
	return expr;
}

// Jacobi polynomials:
double myFitness::P(const int &n, const double &alpha, const double 
	&beta, const double &xi) const
{
	return jacobiP(n,alpha,beta,xi);
}

// differentials for Jacobi polynomials:
double myFitness::dP(const int &k, const int &n, const double &alpha, 
	const double &beta, const double &xi) const
{
	return P(n-k, alpha+k, beta+k, xi)*
		tgamma(n+alpha+beta+k+1.)/(pow(2., k)*tgamma(n+alpha+beta+1.));
}

double myFitness::xi_mapto_x(const double &xi, const double &a, const 
	double &b) const
{
	return ((b - a)*xi + (b + a))*0.5;
}

double myFitness::x_mapto_xi(const double &x, const double &a, const 
	double &b) const
{
	return 2.*(x - a)/(b - a) - 1.;
}

double myFitness::u_hat(const vector<double> &uk, const vector<double>
	&x, const vector<double> &a, const vector<double> &b ) const
{
	double uprod;
	double usum = 0.;
	if(a.size() != this->nvar || b.size() != this->nvar || uk.size() 
		!= this->N)
		cerr << "Information for uhat are not consistent!"
			<< endl;
	else
	{
		for(size_t i = 0; i < uk.size(); i++)
		{
			uprod = 1.;
			for(size_t v = 0; v < this->nvar; v++)
				uprod *= jacobiP(this->pdeg_vec.at(i).at(v),0.,0.,
					x_mapto_xi(x.at(v),a.at(v),b.at(v)));
			usum += uk.at(i)*uprod;
		}
	}
	return usum;
}

// Monte Carlo Integration
double myFitness::mcint(const vector<double> &f, const double &a, const 
	double &b) const
{
	vector<double> A(this->nvar,a);
	vector<double> B(this->nvar,b);
	return mcint(f, A, B);
}

double myFitness::mcint(const vector<double> &f, const vector<double> 
	&a, const vector<double> &b) const
{
	long long qty = f.size();
	double prodf = 1.;
	double sumf = 0.;
	if(a.size() != this->nvar || b.size() != this->nvar)
	{
		cerr << "Limits for monte carlo integration are not consistent!"
			<< endl;
		return NAN;
	}
	for(size_t v = 0; v < this->nvar; v++)
		prodf *= b.at(v) - a.at(v);
	for(vector<double>::const_iterator it = f.begin(); it != f.end(); 
		it++)
	{
		if(std::isnan(*it) || my::isinf(*it))
			qty--;
		else
		{
			sumf += *it;
		}
	}
	sumf *= (((double)qty < 0.95*(double)f.size())? NAN : prodf/qty);
	return (my::isinf(sumf)? NAN : sumf);
}

double myFitness::solve_for(const vector<AAST> &K_, const AAST &s_, 
	const bool &saveflag)
{
	AAST temp;
	if(K_.size() != this->Q)
	{
		cerr << "Missing terms of K! " << K_.size() << " != " << 
			this->Q << endl;
		return NAN;
	}
	// Not trivial solutions ................
	bool trivial = true;
	for(size_t i = 0; i < this->Q; i++)
	{
		temp = K_.at(i);
		temp.simplify();
		trivial &= my::isLikelyZero(temp.getRealValue());
	}
	temp = s_;
	temp.simplify();
	trivial &= my::isLikelyZero(temp.getRealValue());
	if(trivial)
		return HUGE_VAL;
	//.......................................
	vector<vector<RPN> > K(this->Q);
	for(size_t i = 0; i < this->Q; i++) 
		K.at(i) = convert_AASTtoRPN(K_.at(i));
	vector<RPN> s = convert_AASTtoRPN(s_);
	vector<double> auxvar(this->nvar);
	vector<vector<double> > grandauxvar(this->nmcpts);                    
	vector<vector<double> > mector_K(this->Q); //per xi
	vector<double> vector_S; //per xi
	vector<double> partialfit(this->ndomains);
	for(size_t d = 0; d < this->ndomains; d++)
	{
		for(unsigned long w = 0; w < this->nmcpts; w++)
		{
			for(size_t v = 0; v < this->nvar; v++)
			{
				auxvar.at(v) = xi_mapto_x(this->xis.at(w).at(v), 
					this->dataset.bound_inf.at(d).at(v), 
					this->dataset.bound_sup.at(d).at(v));
			}
			grandauxvar.at(w) = auxvar;
		}
		vector_S = evalRPNexpr(s, grandauxvar, 'n');
		for(size_t k = 0; k < this->Q; k++) 
			mector_K.at(k) = evalRPNexpr(K.at(k), grandauxvar, 'n');
		// Build GSE ---------------------------------------------------
		Matrix MLHS(this->N,this->N);
		Matrix MRHS(this->N,1);
		Matrix Maux(this->Q,this->N);
		vector<double> integrand1(this->nmcpts);
		vector<double> integrand2(this->nmcpts);
		size_t Nc = this->dataset.domain.at(d).coordinates.size();
		// Galerkin
		for(size_t i = 0; i < this->N - Nc; i++)
		{
			for(size_t j = 0; j < this->N; j++)
			{
				for(size_t k = 0; k < this->Q; k++)
				{
					for(unsigned long w = 0; w < this->nmcpts; w++)
					{
						integrand1.at(w) = mector_K.at(k).at(w)*this->
							matrixPn.at(w).at(i)*this->matrixPqm.
							at(w).at(k,j);
					}
					Maux.at(k,j) = mcint(integrand1,-1.,1.);
				}
			}
			MLHS.assign2row(i,((this->CoV_diff.at(d))*Maux).
				as_vector());
			for(unsigned long w = 0; w < this->nmcpts; w++)
				integrand2.at(w) = this->matrixPn.at(w).at(i)*
					vector_S.at(w);
			MRHS.at(i) = mcint(integrand2,-1.,1.);
		}
		// Boundary
		for(size_t i = 0; i < Nc; i++)
		{
			for(size_t j = 0; j < this->N; j++)
			{
				double prod_ = 1.;
				for(size_t v = 0; v < this->nvar; v++)
					prod_ *= P(this->pdeg_vec.at(j).at(v), 0., 0.,
						x_mapto_xi(this->dataset.domain.at(d).
						coordinates.at(i).at(v), this->dataset.
						bound_inf.at(d).at(v), this->dataset.bound_sup.
						at(d).at(v)));
					MLHS.at(this->N-i-1,j) = prod_;
			}
			MRHS.at(this->N-i-1) = this->dataset.domain.at(d).values.
				at(i);
		}
		// Enough?
		if(!MLHS.isValid() || !MRHS.isValid())
		{
			partialfit.at(d) = NAN;
			break;
		}
		else
		{
			size_t q = 0;
			while (!MLHS.isFullRank() && q < Nc)
			{ // same extras as in domain
				for(size_t j = 0; j < this->N; j++)
				{
					double prod_ = 1.;
					for(size_t v = 0; v < this->nvar; v++)
						prod_ *= P(this->pdeg_vec.at(j).at(v), 0., 0.,
							x_mapto_xi(this->dataset.extras.at(d).
							coordinates.at(q).at(v), this->dataset.
							bound_inf.at(d).at(v), this->dataset.
							bound_sup.at(d).at(v)));
						MLHS.at(this->N-Nc-q-1,j) = prod_;
				}
				MRHS.at(this->N-Nc-q-1) = this->dataset.extras.at(d).
					values.at(q);
				q++;
			}
		}
		if(!MLHS.isFullRank()) // still not full rank
			cerr << "still not full rank ..." << endl;
		this->U.at(d) = MLHS.pseudoinverse()*MRHS;
		double cond_err = my::SumAbs((MRHS - MLHS*this->U.at(d)).
			as_vector());
		if(cond_err/this->N > 1.) // ill conditioned 
			// difference between them is greater than the unity
			partialfit.at(d) = NAN;
		else
		{
			if(saveflag) //=================================
			{
				this->saveAs(U.at(d), this->dataset.bound_inf.at(d).
					location, this->dataset.bound_sup.at(d).location, 
					((d == 0)? false : true));
			} //============================================
			partialfit.at(d) = 0.;
			for(size_t i = 0; i < this->dataset.reference.at(d).
				coordinates.size(); i++)
			{
				partialfit.at(d) += 
					abs(this->dataset.reference.at(d).
					values.at(i) - u_hat(this->U.at(d).as_vector(), 
					this->dataset.reference.at(d).coordinates.at(i).
					location, this->dataset.bound_inf.at(d).
					location, this->dataset.bound_sup.at(d).
					location));
			}
			partialfit.at(d) /= this->dataset.reference.at(d).
				coordinates.size();
		}
	}
	// Fitness
	int nn = partialfit.size();
	double overall_fitness = 0;
	for(vector<double>::const_iterator it = partialfit.begin(); 
		it != partialfit.end(); it++)
	{
		if(std::isnan(*it) || my::isinf(*it) || my::isLikelyZero(*it)) 
		{
			nn--;
		}
		else
		{
			overall_fitness += *it;
		}
	}
	int hmb = partialfit.size() - nn;
	return pow(2,hmb)*(overall_fitness/nn+hmb);
}
// =====================================================================
#endif
