/* code by: iperetta@ieee.org */
#ifndef CEXPAAST_H_
#define CEXPAAST_H_

#include "cAAST.h"
using namespace std;
// =====================================================================
struct UNK : Node { // Unknown
	size_t index;
	size_t reference;
	UNK(const size_t &i = 0, const size_t &ref = 0) : 
		index(i), reference(ref) {}
	char isType() const { return 'k'; }
	string getSymbol() const { return "k"+my::num2str<size_t>(index); }
	size_t getArity() const { return 0; }
	bool operator== (const UNK &other) const {
		return (index == other.index && reference == other.reference);
	}
	bool operator!= (const UNK &other) const {
		return !(operator==(other));
	}
};
// =====================================================================
struct U_X : Node  { // Unknown function - solution of PDE
	size_t reference;
	U_X(const size_t &ref = 0) : reference(ref) {}
	char isType() const { return 'u'; }
	string getSymbol() const { return "UX"; }
	size_t getArity() const { return 0; }
	bool operator== (const U_X &other) const {
		return (reference == other.reference);
	}
	bool operator!= (const U_X &other) const {
		return !(operator==(other));
	}
};
// =====================================================================
struct DDX : Node { // Derivative operator
	size_t wrt_id;
	size_t order;
	DDX(const size_t &wrt = 0, const size_t &ord = 1) : 
		wrt_id(wrt), order(ord) {}
	char isType() const { return 'd'; }
	string getSymbol() const { 
		return "D_x" + my::num2str<size_t>(wrt_id) + ((order == 1) ? "" 
			: "_" + my::num2str<size_t>(order)); }
	size_t getArity() const { return 1; }
	bool operator== (const DDX &other) const {
		return (wrt_id == other.wrt_id && order == other.order);
	}
	bool operator!= (const DDX &other) const {
		return !(operator==(other));
	}
};
// =====================================================================
struct POL : Node { // Polynomial representation
	int degree;
	size_t wrt_id;
	double alpha;
	double beta;
	POL(const int &deg, const size_t &wrt = 0, const double &alp = 0., 
		const double &bet = 0.) : degree(deg), wrt_id(wrt), 
		alpha(alp), beta(bet) {}
	char isType() const { return 'p'; }
	string getSymbol() const { 
		return "P_" + my::num2str<int>(degree) + ( 
			(my::isLikelyZero(alpha) && my::isLikelyZero(beta)) ? "" : 
			"^{" + my::num2str<double>(alpha) + "," + 
			my::num2str<double>(beta) + "}" ) + "(x" + 
			my::num2str<size_t>(wrt_id) + ")"; 
	}
	size_t getArity() const { return 0; }
	bool operator== (const POL &other) const {
		return (degree == other.degree && wrt_id == other.wrt_id &&
			alpha == other.alpha && beta == other.beta);
	}
	bool operator!= (const POL &other) const {
		return !(operator==(other));
	}
};
// =====================================================================
// Jacobi Polynomials  -------------------------------------------------
struct JacobiP {
	int n;
	double alpha;
	double beta;
	size_t wrt_id;
	double bound_inf;
	double bound_sup;
    JacobiP() : n(0), alpha(0.), beta(0.),	wrt_id(0), bound_inf(-1.),
		bound_sup(1.) {};
    JacobiP(const int &n, const double &_alpha, const double &_beta, 
		const size_t &varid = 0, const double &inf = -1., 
		const double &sup = 1.) : n(n), alpha(_alpha), beta(_beta), 
		wrt_id(varid), bound_inf((inf < sup)? inf : sup), 
		bound_sup((inf < sup)? sup : inf) {};
	void set(const int &n, const double &_alpha, const double &_beta, 
		const size_t &varid = 0, const double &inf = -1., 
		const double &sup = 1.);
	AAST T() const;
};

void JacobiP::set(const int &n, const double &_alpha, const double 
		&_beta, const size_t &varid, const double &inf, const double 
		&sup) {
	this->alpha = _alpha;
	this->beta = _beta;
	this->n = n;
	this->wrt_id = varid;
	this->bound_inf = inf;
	this->bound_sup = sup;
}
AAST JacobiP::T() const {
	// Orthogonal Polynomials - Gabor Szegö (1959) - pp.62
	/*
	 * MAXIMA Implementation
	 * 
	 * P(n,a,b,x):=ratsimp(1/(n!)*sum((if v = 0 then prod(a+i,i,1,n) 
	 * else if v = n then prod(n+a+b+i,i,1,n) else 
	 * (n!/(v!*(n-v)!))*prod(n+a+b+i,i,1,v)*prod(a+i,i,v+1,n))*
	 * ((x-1)/2)^v,v,0,n));
	 * 
	 * same as:
	 * 
	 * expand_hypergeometric: true$
	 * P(n,alpha,beta,x):= if n < 0 then 0 else if n+alpha+beta+1 = 0 
	 * then 1 else ratsimp(gamma(n+alpha+1)/(gamma(alpha+1)*gamma(n+1))*
	 * hypergeometric([-n,1+n+alpha+beta],[alpha+1],1/2*(1-x)));
	 *  
     */
    AAST expr; 
	if(this->n <= 0) {
		if (this->n == 0)
			expr = 1.;
		else
			expr = 0.;
	} else {
		FNC pwr("POW",2);
		VAR wrt(this->wrt_id);
		AAST wrt_term = (this->bound_sup - wrt.T())/(this->bound_sup 
			- this->bound_inf);
		AAST aux;
		expr = 0.;
		double coeff = my::combination(this->n + this->alpha, this->n);
		double sumcoeff = 0.;
		for(int k = 0; k <= this->n; k++) {
			sumcoeff = my::pochhammer(-(this->n), (double)k)*
				my::pochhammer(this->n+this->alpha+this->beta+1., 
				(double)k)/(my::pochhammer(this->alpha+1.,(double)k)*
				tgamma((double)k+1.));
			aux = pwr.T(wrt_term << (double)k);
			expr += sumcoeff*aux;
		}
		expr *= coeff;
		expr.simplify();
	}
	return expr;
}
// =====================================================================
vector<vector<int> > getPowers(const int &n_var, const int &p_deg) 
{
	vector<vector<int> > int_part, int_pwrs_all; 
	vector<int> int_pwrs;
	size_t i, j;
	vector<int> ith_powers(n_var,0.);
	vector<vector<int> > powers;
	powers.push_back(ith_powers);
	for(int p = 1; p <= p_deg; p++)
	{
		int_pwrs_all.clear();
		int_part = my::integer_partition(p,n_var);
		for(i = 0; i < int_part.size(); i++) 
		{
			int_pwrs = int_part.at(i);
			int_pwrs_all.push_back(int_pwrs);
			while(prev_permutation(int_pwrs.begin(),int_pwrs.end()))
				int_pwrs_all.push_back(int_pwrs);
		}
		sort(int_pwrs_all.begin(), int_pwrs_all.end());
		reverse(int_pwrs_all.begin(), int_pwrs_all.end());
		for(i = 0; i < int_pwrs_all.size(); i++) 
		{
			ith_powers.clear();
			for(j = 0; j < (size_t)n_var; j++)
				ith_powers.push_back(int_pwrs_all.at(i).at(j));
			powers.push_back(ith_powers);
		}
	}
	return powers;
}
// =====================================================================
class Basis {
  public:
	size_t n_variables, 
		   span_degree;
	double alpha, 
		   beta;
	vector<AAST> span;
	vector<AAST> span_mapped;
	Basis(const size_t &updeg, const double &alpha = 0., 
		const double &beta = 0.);
	Basis(const size_t &nvar, const size_t &updeg, const double 
		&alpha, const double &beta);
	~Basis() {};
	void map_to(const double &bd_inf, const double &bd_sup);
	void map_to(const vector<double> &bd_inf, const 
		vector<double> &bd_sup);
};

Basis::Basis(const size_t &updeg, const double 
	&alpha, const double &beta) : n_variables(1), span_degree(updeg), 
	alpha(alpha), beta(beta)
{
	AAST aux;
	JacobiP polynomial;
	size_t N = updeg+1;
	for(size_t n = 0; n < N; n++)
	{
		polynomial.set(n, alpha, beta, 0, -1., 1.);
		aux = polynomial.T();
		aux.simplify();
		this->span.push_back(aux);
	}
}
Basis::Basis(const size_t &nvar, const size_t &updeg, const double 
	&alpha, const double &beta) : n_variables(nvar), span_degree(updeg), 
	alpha(alpha), beta(beta)
{
	AAST aux;
	JacobiP polynomial;
	vector<vector<int> > powers = getPowers(nvar, updeg);
	size_t N = powers.size();
	for(size_t n = 0; n < N; n++)
	{
		aux = 1.;
		for(size_t i = 0; i < nvar; i++)
		{
			polynomial.set(powers.at(n).at(i), alpha, beta, i);
			aux *= polynomial.T();
		}
		aux.simplify();
		this->span.push_back(aux);
	}
}
void Basis::map_to(const double &bd_inf, const double &bd_sup)
{
	vector<double> bd_infv(1,bd_inf), bd_supv(1,bd_sup);
	this->map_to(bd_infv, bd_supv);
}
void Basis::map_to(const vector<double> &bd_inf, const 
	vector<double> &bd_sup)
{
	if(bd_inf.size() != this->n_variables || bd_sup.size() != 
		this->n_variables)
	{
		cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
			<< endl;
		string msg = "Invalid number of boundary values to map";
		throw std::invalid_argument(msg.c_str());
	}
	AAST aux;
	vector<AAST> mapvar;
	for(size_t j = 0; j < this->n_variables; j++) 
	{
		aux = 2.*(static_cast<VAR>(j).T() - bd_inf.at(j))/
			(bd_sup.at(j) - bd_inf.at(j)) - 1.;
		aux.simplify();
		mapvar.push_back(aux);
	}
	span_mapped.clear();
	for(size_t i = 0; i < this->span.size(); i++)
	{
		aux = span.at(i);
		for(size_t j = 0; j < this->n_variables; j++)
		{
			aux.map_variable_to(j,mapvar.at(j));
		}
		aux.simplify();
		span_mapped.push_back(aux);
	}
}
// =====================================================================
AAST::AAST(const AAST &source) { //copy constructor
	this->parent = NULL;
	this->children = source.children;
	if(!this->children.empty())
		for(size_t i = 0; i < this->children.size(); i++)
			this->children.at(i).parent = this;
    switch(source.node->isType()) {
	  case 'c':
		this->node = new CNS(static_cast<CNS*>(source.node)->value);
		break;
	  case 'v':
		this->node = new VAR(static_cast<VAR*>(source.node)->index);
		break;
	  case 'o':
		this->node = new OPR(static_cast<OPR*>(source.node)->symbol,
			static_cast<OPR*>(source.node)->arity);
		break;
	  case 'f':
		this->node = new FNC(static_cast<FNC*>(source.node)->symbol,
			static_cast<FNC*>(source.node)->arity);
		break;
	  case 'k':
		this->node = new UNK(static_cast<UNK*>(source.node)->index, 
			static_cast<UNK*>(source.node)->reference);
		break;
	  case 'u':
		this->node = new U_X(static_cast<U_X*>(source.node)->reference);
		break;
	  case 'd':
		this->node = new DDX(static_cast<DDX*>(source.node)->wrt_id,
			static_cast<DDX*>(source.node)->order);
		break;
	  case 'p':
		this->node = new POL(static_cast<POL*>(source.node)->degree,
			static_cast<POL*>(source.node)->wrt_id, 
			static_cast<POL*>(source.node)->alpha, 
			static_cast<POL*>(source.node)->beta);
		break;
	}
}
bool AAST::operator==(const AAST &other) const {
	bool chk;
	switch(this->node->isType())
	{
		case 'c':
			chk = ((other.node->isType() == 'c') ? static_cast<CNS*>
				(this->node)->operator==(*static_cast<CNS*>(other.node))
				 : false);
		break;
		case 'v':
			chk = ((other.node->isType() == 'v') ? static_cast<VAR*>
				(this->node)->operator==(*static_cast<VAR*>(other.node))
				 : false);
		break;
		case 'o':
			chk = ((other.node->isType() == 'o') ? static_cast<OPR*>
				(this->node)->operator==(*static_cast<OPR*>(other.node))
				 : false);
		break;
		case 'f':
			chk = ((other.node->isType() == 'f') ? static_cast<FNC*>
				(this->node)->operator==(*static_cast<FNC*>(other.node))
				 : false);
		break;
		case 'k':
			chk = ((other.node->isType() == 'c') ? static_cast<UNK*>
				(this->node)->operator==(*static_cast<UNK*>(other.node))
				 : false);
		break;
		case 'u':
			chk = ((other.node->isType() == 'v') ? static_cast<U_X*>
				(this->node)->operator==(*static_cast<U_X*>(other.node))
				 : false);
		break;
		case 'd':
			chk = ((other.node->isType() == 'o') ? static_cast<DDX*>
				(this->node)->operator==(*static_cast<DDX*>(other.node))
				 : false);
		break;
		case 'p':
			chk = ((other.node->isType() == 'f') ? static_cast<POL*>
				(this->node)->operator==(*static_cast<POL*>(other.node))
				 : false);
		break;
	}
	if(!chk || this->children.size() != other.children.size())
		return false;
	else
	{
		if(!this->children.empty())
		{
			vector<size_t> indices;
			bool found;
			for(size_t i = 0; i < this->children.size(); i++)
			{
				found = false;
				for(size_t j = 0; j < other.children.size(); j++)
				{
					if(!my::is_in(j,indices))
						if(this->children.at(i) == other.children.at(j))
						{
							found = true;
							indices.push_back(j);
							break;
						}
				}
				if(!found)
				{
					chk = false;
					break;
				}
			}
		}
		return chk;
	}
}
void AAST::assign2unknowns(const size_t &qty, const double &value) {
	vector<double> vvalue(qty,value);
	this->assign2unknowns(vvalue);
}
void AAST::assign2unknowns(const double &value) {
	vector<double> vvalue(1,value);
	this->assign2unknowns(vvalue);
}
void AAST::assign2unknowns(const vector<double> &values) {
	if(this->node->isType() == 'k')
	{
		size_t index = static_cast<UNK*>(this->node)->index;
		AAST newCNS = values.at(index);
		this->safe_swap(newCNS,this->parent);
	}
	else
		if(!this->children.empty())
			for(size_t i = 0; i < this->children.size(); i++)
				this->children.at(i).assign2unknowns(values);
}
double AAST::eval_with_unknowns(const vector<double> &values) const {
	AAST tree(*this);
	tree.assign2variables(values);
	tree.simplify();
	return tree.getRealValue();
}
vector<size_t> AAST::getTypeIdx(const string &type) const {
	vector<size_t> idx;
	vector<const AAST*> mapped = this->map();
	for(size_t i = 0; i < mapped.size(); i++) {
		for(size_t j = 0; j < type.size(); j++) {
			if(mapped.at(i)->node->isType() == type.at(j))
				idx.push_back(i);
		}
	}
	return idx;
}
void AAST::shrink_derivative_operator(const size_t &n_vars) {
	vector<size_t> idx, opers;
	vector<int> wrt;
	vector<vector<size_t> > clusters;
	int i, j, k; 
	size_t localroot;
	AAST newopr, aux;
	DDX dif;
	idx = this->getTypeIdx("d");
	while(idx.size() > 0) {
		i = idx.back();
		opers.push_back(i);
		idx.pop_back();
		if((idx.size() == 0) ? true : !(&this->at(idx.back()) == 
				this->at(i).parent)) {
			if(opers.size() > 1) {
				clusters.push_back(opers);
			}
			opers.clear();
		}
	}
	for(i = 0; i < (int)clusters.size(); i++) {
		wrt.assign(n_vars,0);
		localroot = clusters.at(i).back();
		for(j = 0; j < (int)clusters.at(i).size(); j++) {
			DDX* ptddx = static_cast<DDX*>(this->at(clusters.at(i).
				at(j)).node);
			k = ptddx->wrt_id;
			if(k >= (int)n_vars) {
				cerr << "[Error] In " << __FILE__ << " at line " << 
					__LINE__ << endl;
				string msg = "Derivative wrt invalid variable; review.";
				throw std::out_of_range(msg.c_str());
			}
			wrt.at(k) += ptddx->order;
		}
		newopr = this->at(clusters.at(i).front()).children.at(0);
		for(j = wrt.size()-1; j >= 0; j--)
			if(wrt.at(j) != 0) {
				dif.wrt_id = j;
				dif.order = wrt.at(j);
				newopr = dif.T(newopr);
			}
		this->at(localroot).safe_swap(newopr,
			this->at(localroot).parent);
	}
}

void AAST::expand_unknown_function(const size_t &n_vars, 
		const size_t &up2degree) {
	vector<vector<int> > basispowers = getPowers(n_vars, up2degree);
	UNK un;
	POL polynomial(0,0,0.,0.);
	AAST aux, newopr;
	AAST *opr;
	newopr = 0.;
	for(size_t n = 0; n < basispowers.size(); n++) {
		un.index = n;
		aux = un.T();
		for(size_t i = 0; i < basispowers.at(n).size(); i++) {
			polynomial.wrt_id = i;
			polynomial.degree = basispowers.at(n).at(i);
			aux *= polynomial.T();
		}
		newopr += aux;
	}
	newopr.simplify();
	vector<size_t> idx = this->getTypeIdx("u");
	for(int i = idx.size()-1; i >= 0; i--) {
		opr = &(this->at(idx.at(i)));
		aux = newopr;
		opr->safe_swap(aux,opr->parent);
	}
}

void AAST::expand_derivative_argument(const vector<double> &bound_inf, 
		const vector<double> &bound_sup) {
	double coeff;
	size_t ord, varid;
	POL poly(0);
	AAST aux, newopr;
	AAST *opr;
	vector<size_t> idx = this->getTypeIdx("d");
	vector<const AAST*> mapped;
	for(int i = idx.size()-1; i >= 0; i--) {
		opr = &(this->at(idx.at(i)));
		newopr = opr->children.at(0);
		mapped = newopr.map();
		for(int j = mapped.size()-1; j >= 0; j--) {
			switch(mapped.at(j)->node->isType()) {
			  case 'c':
				aux = *mapped.at(j);
				break;
			  case 'k':
				aux = *mapped.at(j);
				break;
			  case 'p':
				poly = *(static_cast<POL*>(mapped.at(j)->node));
				varid = static_cast<DDX*>(opr->node)->wrt_id;
				if(varid > bound_inf.size() || varid > bound_sup.size())
				{
					cerr << "[Error] In " << __FILE__ << " at line " 
						<< __LINE__ << endl;
					string msg = "Invalid attempt to expand polynomial";
					msg += "; review inf and/or sup";
					throw std::out_of_range(msg.c_str());
				}
				if(varid == poly.wrt_id){
					ord = static_cast<DDX*>(opr->node)->order;
					coeff = tgamma((double)poly.degree + poly.alpha + 
						poly.beta + (double)ord + 1.)/(pow(
						max(bound_sup.at(varid),bound_inf.at(varid)) - 
						min(bound_sup.at(varid),bound_inf.at(varid)),
						(double)ord)*tgamma((double)poly.degree +
						poly.alpha + poly.beta + 1.));
					poly.degree -= (int)ord;
					poly.alpha += (double)ord;
					poly.beta += (double)ord;
					aux = coeff*poly.T();
				} else {
					aux = poly.T();
				}
				break;
			  case 'o':
				aux = *mapped.at(j);
				if(mapped.at(j)->node->getSymbol().compare("*") == 0 ||
					mapped.at(j)->node->getSymbol().compare("+") == 0)
						break;
			  default:
				cerr << "[Error] In " << __FILE__ << " at line " << 
					__LINE__ << endl;
				string msg(1,mapped.at(j)->node->isType());
				msg += ": " + mapped.at(j)->node->getSymbol();
				msg += ": unknown/invalid sub-node for derivative.";
				throw std::invalid_argument(msg.c_str());
			}
			newopr.at(j).safe_swap(aux,newopr.at(j).parent);
		}
		opr->safe_swap(newopr,opr->parent);
	}
}

void AAST::expand_polynomials(const vector<double> &inf, const 
	vector<double> &sup) 
{
	vector<size_t> idx = this->getTypeIdx("p");
	AAST newpoly, newvar;
	AAST *polyaddr;
	int n;
	size_t varid;
	double alpha, beta;
	JacobiP jpoly;
	for(int i = idx.size()-1; i >= 0; i--) 
	{
		polyaddr = &(this->at(idx.at(i)));
		n = static_cast<POL*>(polyaddr->node)->degree;
		if(n < 0)
			newpoly = 0.;
		else
		{
			if(n == 0)
				newpoly = 1.;
			else
			{
				varid = static_cast<POL*>(polyaddr->node)->wrt_id;
				alpha = static_cast<POL*>(polyaddr->node)->alpha;
				beta = static_cast<POL*>(polyaddr->node)->beta;
				jpoly.set(n,alpha,beta,varid,inf.at(varid),
					sup.at(varid));
				newpoly = jpoly.T();
			}
		}
		polyaddr->safe_swap(newpoly,polyaddr->parent);
	}
}

// =====================================================================
// Node AASTree representation  ----------------------------------------
AAST Node::T() const {
	AAST tree;
	switch(this->isType()) {
	  case 'c':
		tree.node = new CNS(static_cast<const CNS*>(this)->value);
		break;
	  case 'v':
		tree.node = new VAR(static_cast<const VAR*>(this)->index);
		break;
	  case 'k':
		tree.node = new UNK(static_cast<const UNK*>(this)->index,
			static_cast<const UNK*>(this)->reference);
		break;
	  case 'u':
		tree.node = new U_X(static_cast<const U_X*>(this)->reference);
		break;
	  case 'p':
		tree.node = new POL(static_cast<const POL*>(this)->degree,
			static_cast<const POL*>(this)->wrt_id,
			static_cast<const POL*>(this)->alpha,
			static_cast<const POL*>(this)->beta);
		break;
	  default:
		cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
			<< endl;
		string msg = this->getSymbol() + " need arguments!";
		throw std::invalid_argument(msg.c_str());
	}
	return tree;
}
AAST Node::T(const AAST &arg) const {
	AAST tree;
	char type = this->isType();
	if(this->getArity() != 1)
		type = '\0';
	switch(type) {
	  case 'f':
		tree.node = new FNC(static_cast<const FNC*>(this)->symbol,
			static_cast<const FNC*>(this)->arity);
		break;
	  case 'd':
		tree.node = new DDX(static_cast<const DDX*>(this)->wrt_id, 
			static_cast<const DDX*>(this)->order);
		if(arg.node->isType() == 'd' || arg.node->isType() == 'u')
			break;
	  default:
		cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
			<< endl;
		string msg = "Review arguments for: " + this->getSymbol();
		throw std::invalid_argument(msg.c_str());
	}
	tree.children.push_back(arg);
	tree.children.at(0).parent = &tree;
	return tree;
}
// =====================================================================
// =====================================================================
#define nchar_RPN_label 5

struct RPN {
	char type;
	double value;
	size_t index_arity;
	char label[nchar_RPN_label];
	RPN() : type('?'), value(NAN), index_arity(0)
	{
		label[0] = '\0';
	}
	RPN(const double &vl) : type('c'), value(vl), index_arity(0)
	{
		label[0] = '\0';
	} 
	RPN(const size_t &id) : type('v'), value(NAN), index_arity(id)
	{
		label[0] = '\0';
	} 
	RPN(const char *lbl, const size_t &ar) : type('f'), 
		value(NAN), index_arity(ar) 
	{
		strncpy(label, lbl, sizeof(label)/sizeof(char));
	}
	void CNS(const double &vl)
	{
		type = 'c';
		value = vl;
		index_arity = 0;
		label[0] = '\0';
	}
	void VAR(const size_t &id)
	{
		type = 'v';
		value = NAN;
		index_arity = id;
		label[0] = '\0';
	}
	void FOP(const char *lbl, const size_t &ar)
	{
		type = 'f';
		value = NAN;
		index_arity = ar;
		strncpy(label, lbl, nchar_RPN_label-1);
		label[nchar_RPN_label-1] = '\0';
	}
	char isType() const { return type; };
};

vector<RPN> convert_AASTtoRPN(const AAST &tree) 
{
	RPN aux;
	vector<RPN> rpnexpr;
	vector<const AAST*> mapped = tree.map();
	for(int i = mapped.size()-1, j = 0; i >= 0; i--, j++) {
		switch(mapped.at(i)->getNode()->isType()) {
		  case 'c':
			aux.CNS(static_cast<CNS*>(mapped.at(i)->getNode())->value);
			break;
		  case 'v':
			aux.VAR(static_cast<VAR*>(mapped.at(i)->getNode())->index);
			break;
		  case 'k':
			aux.VAR(static_cast<UNK*>(mapped.at(i)->getNode())->index);
			break;
		  case 'f':
		  case 'o':
			aux.FOP(mapped.at(i)->getNode()->getSymbol().c_str(), 
				mapped.at(i)->getNode()->getArity());
			break;
		  default:
			cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
				<< endl;
			throw std::invalid_argument("Conversion not allowed");
		}
		rpnexpr.push_back(aux);
	}
	return rpnexpr;
}

void showRPNexpr(const vector<RPN> &rpn_expr)
{
	for(vector<RPN>::const_iterator i = rpn_expr.begin(); i != 
		rpn_expr.end(); i++)
	{
		cout << (((*i).type == 'c') ? my::num2str((*i).value) : 
			(((*i).type == 'v') ?  "x"+my::num2str((*i).index_arity) : 
			(*i).label)) << " | ";
	}
	cout << endl;
}

double evalRPNexpr(const vector<RPN> &rpn_expr, 
	const vector<double> &var_values, const char &val_of_nan = 'n') 
{
	size_t N_expr = rpn_expr.size();
	my::STACK<double> stack;
	double A, B;
	double auxv;
	string labels;
	double results;
	for(size_t i = 0; i < N_expr; i++) 
	{
		switch(rpn_expr.at(i).type) 
		{
		  case 'c':
			stack.push(rpn_expr.at(i).value);
			break;
		  case 'v':
			stack.push(var_values.at(rpn_expr.at(i).index_arity));
			break;
		  case 'f':
			labels = rpn_expr.at(i).label;
			if(rpn_expr.at(i).index_arity == 1) 
			{
				A = stack.pop();
				if(labels.compare("EXP") == 0)
					auxv = exp(A);
				else
				if(labels.compare("LOG") == 0)
					auxv = (A <= 0) ? 1. : log(A);
				else
				if(labels.compare("LOG2") == 0)
					auxv = (A <= 0) ? 1. : log(A)/log(2.);
				else
				if(labels.compare("LOG10") == 0)
					auxv = (A <= 0) ? 1. : log10(A);
				else
				if(labels.compare("SQRT") == 0)
					auxv = (A <= 0) ? 1. : sqrt(A);
				else
				if(labels.compare("COS") == 0)
					auxv = cos(A);
				else
				if(labels.compare("SIN") == 0)
					auxv = sin(A);
				else
				if(labels.compare("TAN") == 0)
					auxv = tan(A);
				else
				if(labels.compare("ACOS") == 0)
					auxv = acos(A);
				else
				if(labels.compare("ASIN") == 0)
					auxv = asin(A);
				else
				if(labels.compare("ATAN") == 0)
					auxv = atan(A);
				else
				if(labels.compare("COSH") == 0)
					auxv = cosh(A);
				else
				if(labels.compare("SINH") == 0)
					auxv = sinh(A);
				else
				if(labels.compare("TANH") == 0)
					auxv = tanh(A);
				else
				if(labels.compare("POW2") == 0)
					auxv = A*A;
				else
				if(labels.compare("POW3") == 0)
					auxv = A*A*A;
				else
				if(labels.compare("CBRT") == 0)
					auxv = pow(A,1./3.);
				else
				if(labels.compare("EXPn") == 0)
					auxv = exp(-A);
				else
				if(labels.compare("INV") == 0)
					auxv = 1./A;
				else
				if(labels.compare("NEG") == 0)
					auxv = -A;
			}
			else 
			{
			if(rpn_expr.at(i).index_arity == 2) {
				A = stack.pop();
				B = stack.pop();
				if(labels.compare("+") == 0)
					auxv = A+B; 
				else
				if(labels.compare("-") == 0)
					auxv = A-B; 
				else
				if(labels.compare("*") == 0)
					auxv = A*B; 
				else
				if(labels.compare("/") == 0)
					auxv = A/B; 
				else
				if(labels.compare("POW") == 0)
					auxv = pow(A,B);
				else
				if(labels.compare("ATAN2") == 0)
					auxv = atan2(A,B);
			}}
			stack.push((std::isnan(auxv) || my::isinf(auxv)) ?	
				1. : ((my::isLikelyZero(auxv)) ? 0. : auxv)); // protection
			break;
		}
	}
	A = stack.pop();
	if(!std::isnan(A) && !my::isinf(A)) 
		results = A;
	else
	{
		switch(val_of_nan)
		{
			case '0':
				results = 0.;
				break;
			case '1':
				results = 1.;
				break;
			case 'i':
				results = HUGE_VAL;
				break;
			default: //case 'n'
				results = NAN;
		}
	}
	return results;
}

vector<double> evalRPNexpr(const vector<RPN> &rpn_expr, 
	const vector<vector<double> > &var_values, 
	const char &val_of_nan = 'n') 
{
	size_t N_values = var_values.size();
	vector<double> results;
	for(size_t n = 0; n < N_values; n++)
	{
		results.push_back(evalRPNexpr(rpn_expr, var_values.at(n), 
			val_of_nan));
	}
	return results;
}
// =====================================================================
#endif
