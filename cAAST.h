/* code by: iperetta@ieee.org */
#ifndef CAAST_H_
#define CAAST_H_

#include <stdexcept>
#include <cmath>
#include <cstring>
#include "myfun.h"
using namespace std;

class Double {
  friend ostream& operator<<(ostream &os, const Double& v);
  public:
	double value;
	Double() : value(0.) {};
	Double(const double &v) : value(v) {};
	~Double(){};
};

ostream& operator<<(ostream &os, const Double& v) {
	return os << v.value;
}

class AAST;
#ifdef CEXPAAST_H_ 
class Basis;
#endif

struct Node { // Generic representation
	virtual char isType() const = 0;
	virtual string getSymbol() const = 0;
	virtual size_t getArity() const = 0;
	AAST T() const;
	AAST T(const AAST &arg) const;
	AAST T(const vector<AAST> &arg) const;
	virtual ~Node() {};
};

struct CNS : Node { // Constant
	double value;
	CNS(const double &v = NAN) : value(v) {}
	char isType() const { return 'c'; }
	string getSymbol() const { return my::num2str<double>(value); }
	size_t getArity() const { return 0; }
	bool operator== (const CNS &other) const {
		return (my::isLikelyZero(value - other.value));
	}
	bool operator!= (const CNS &other) const {
		return !(operator==(other));
	}
};

struct VAR : Node { // Variable
	size_t index;
	VAR(const size_t &i = 0) : index(i) {}
	char isType() const { return 'v'; }
	string getSymbol() const { return "x"+my::num2str<size_t>(index); }
	size_t getArity() const { return 0; }
	bool operator== (const VAR &other) const {
		return (index == other.index);
	}
	bool operator!= (const VAR &other) const {
		return !(operator==(other));
	}
};

struct OPR : Node { // Operator
	string symbol;
	size_t arity;
	OPR(const char &sym, const size_t &ar = 2) : 
		symbol(1,sym), arity(ar) {}
	OPR(const string &sym, const size_t &ar = 2) : 
		symbol(sym), arity(ar) {}
	char isType() const { return 'o'; }
	string getSymbol() const { return symbol; }
	size_t getArity() const { return arity; }
	bool operator== (const OPR &other) const {
		return (symbol.compare(other.symbol) == 0 && 
			arity == other.arity);
	}
	bool operator!= (const OPR &other) const {
		return !(operator==(other));
	}
};

struct FNC : Node  { // Function
	string symbol;
	size_t arity;
	FNC(const string &sym, const size_t &ar) : 
		symbol(sym), arity(ar) {}
	char isType() const { return 'f'; }
	string getSymbol() const { return symbol; }
	size_t getArity() const { return arity; }
	bool operator== (const FNC &other) const {
		return (symbol.compare(other.symbol) == 0 && 
			arity == other.arity);
	}
	bool operator!= (const FNC &other) const {
		return !(operator==(other));
	}
};

class AAST { // Algebraic Abstract Syntax Tree
	friend struct Node;
	friend struct CNS;
	friend struct VAR;
	friend struct OPR;
	friend struct FNC;
#ifdef CEXPAAST_H_ 
	friend class Basis;
#endif
	friend AAST operator+(const double &value, const AAST &rhs);
	friend AAST operator-(const double &value, const AAST &rhs);
	friend AAST operator*(const double &value, const AAST &rhs);
	friend AAST operator/(const double &value, const AAST &rhs);
	friend vector<AAST> operator<<(const vector<AAST> &arg1, 
		const AAST &arg2);
	friend vector<AAST> operator<<(const double &arg1, 
		const AAST &arg2);
	friend vector<AAST> operator<<(const double &arg1, 
		const vector<AAST> &arg2);
	friend vector<AAST> operator<<(const Double &arg1, 
		const Double &arg2);
	friend vector<AAST> operator<<(const vector<AAST> &arg1, 
		const vector<AAST> &arg2);
	friend void operator<<=(vector<AAST> &arg1, const AAST &arg2);
	friend void operator<<=(vector<AAST> &arg1, 
		const vector<AAST> &arg2);
  protected:
    Node* node;
    AAST* parent;
    vector<AAST> children;
  public:
    AAST(): node(NULL), parent(NULL) {};
    AAST(const double& constant);
    AAST(const AAST &source);    
    AAST& operator= (AAST source); 
    void swap(AAST &other);
    void safe_swap(AAST &other, AAST *root);
	void fixParents(AAST *root);
    ~AAST();
    Node* getNode() const { return node; };
    AAST* getParent() const { return parent; };
    vector<AAST> getChildren() const {return children; };
    vector<const AAST*> map() const;
    size_t size() const;
	AAST& at(const size_t &i_);
    const AAST& at(const size_t &i_) const;
    string str(const bool &infix = false) const;
    void show(const bool &infix = false) const;
    bool operator==(const AAST &other) const;
    bool operator!=(const AAST &other) const;
    vector<AAST> operator<<(const double &other);
    vector<AAST> operator<<(const AAST &other);
    vector<AAST> operator<<(const vector<AAST> &other);
	void operator+=(const double &value);
	void operator+=(const AAST &operand);
	void operator-=(const double &value);
	void operator-=(const AAST &operand);
	void operator*=(const double &value);
	void operator*=(const AAST &operand);
	void operator/=(const double &value);
	void operator/=(const AAST &operand);
	AAST operator+(const double &value);
	AAST operator+(const AAST &operand);
	AAST operator-(const double &value);
	AAST operator-(const AAST &operand);
	AAST operator*(const double &value);
	AAST operator*(const AAST &operand);
	AAST operator/(const double &value);
	AAST operator/(const AAST &operand);
	// simplify
	void simp_DivisionToProduct();
	void simp_DifferenceToSum();
	bool simp_JoinProducts();
	bool simp_JoinSums();
	bool simp_Constants();
	bool simp_Distributive();
	bool simp_Identities();
	void simplify();
	void assign2variables(const size_t &qty, const double &value);
	void assign2variables(const double &value);
	void assign2variables(const vector<double> &values);
	double getRealValue() const;
	double eval_with_variables(const vector<double> &values) const;
	void map_variable_to(const size_t &index, const AAST &newvar);
#ifdef CEXPAAST_H_
	void assign2unknowns(const size_t &qty, const double &value);
	void assign2unknowns(const double &value);
	void assign2unknowns(const vector<double> &values);
	double eval_with_unknowns(const vector<double> &values) const;
	vector<size_t> getTypeIdx(const string &type) const;
	void shrink_derivative_operator(const size_t &n_vars);
	void expand_unknown_function(const size_t &n_vars, 
		const size_t &up2degree);
	void expand_derivative_argument(const vector<double> &b_inf, 
		const vector<double> &b_sup);
	void expand_polynomials(const vector<double> &inf, const 
		vector<double> &sup);
#endif	
};

AAST::AAST(const double& constant) : parent(NULL) {
	this->node = new CNS(constant);
	this->children.clear();
}
#ifndef CEXPAAST_H_
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
	}
}
#endif
AAST& AAST::operator= (AAST source) { //assignement operator
    if(this != &source) {
        this->swap(source);
    }
    return *this;
}
void AAST::swap(AAST &other) {
	std::swap(this->node, other.node); 
    std::swap(this->parent, other.parent); 
    std::swap(this->children, other.children);
}
void AAST::safe_swap(AAST &other, AAST *root) {
	this->swap(other);
	this->fixParents(root);
}
void AAST::fixParents(AAST *root) {
	this->parent = root;
	for(size_t i = 0; i < this->children.size(); i++)
		this->children.at(i).fixParents(this);
}
AAST::~AAST() {	//destructor
	// Never try to delete this->parent, 
	// it is not malloc'ed neither new'ed!!
	if(this->node != NULL) 
		delete this->node;
	this->children.clear();
}
vector<const AAST*> AAST::map() const {
	vector<const AAST*> mapped;
	vector<const AAST*> visited_nodes;
	visited_nodes.push_back(this);
	while(visited_nodes.size() > 0) {
		mapped.push_back(visited_nodes.back());
		visited_nodes.pop_back();
		for(int i = mapped.back()->children.size()-1; i >= 0; i--)
			visited_nodes.push_back(&(mapped.back()->children.at(i)));
	}
	return mapped;
}
size_t AAST::size() const {
	size_t s = 0;
	if(this->node != NULL)
	{
		s++;
		if(!this->children.empty())
			for(size_t i = 0; i < this->children.size(); i++)
				s += this->children.at(i).size();
	}
	return s;
}
AAST& AAST::at(const size_t &pos) {
	int pos_ = pos;
	vector<AAST*> visited_nodes;
	if(this->node != NULL) {
		AAST* current_node;
		visited_nodes.push_back(this);
		while(visited_nodes.size() > 0 && pos_ > 0) {
			current_node = visited_nodes.back();
			visited_nodes.pop_back();
			for(int i = current_node->children.size()-1; i >= 0; i--)
				visited_nodes.push_back(&(current_node->
					children.at(i)));
			pos_--;
		}
	}
	if(visited_nodes.empty() || pos_ > 0) {
		cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
			<< endl;
		string msg = "Index exceeded tree size.";
		throw std::out_of_range(msg.c_str());
	}
	return *(visited_nodes.back());
}
const AAST& AAST::at(const size_t &pos) const {
	int pos_ = pos;
	vector<const AAST*> visited_nodes;
	if(this->node != NULL) {
		const AAST* current_node;
		visited_nodes.push_back(this);
		while(visited_nodes.size() > 0 && pos_ > 0) {
			current_node = visited_nodes.back();
			visited_nodes.pop_back();
			for(int i = current_node->children.size()-1; i >= 0; i--)
				visited_nodes.push_back(&(current_node->
					children.at(i)));
			pos_--;
		}
	}
	if(visited_nodes.empty() || pos_ > 0) {
		cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
			<< endl;
		string msg = "Index exceeded tree size.";
		throw std::out_of_range(msg.c_str());
	}
	return *(visited_nodes.back());
}
string AAST::str(const bool &infix) const {
	string ST("< empty > ");
	if(this->node != NULL) {
		if(infix) {
			switch(this->node->getArity()) {
			  case 0:
				ST = this->node->getSymbol();
				break;
			  case 1:
			    if(this->node->getSymbol().compare("POW2") == 0) {
					ST = "( ";
					ST += this->children.at(0).str(true);
					ST += " ^ 2 )";
				} else {
				if(this->node->getSymbol().compare("POW3") == 0) {
					ST = "( ";
					ST += this->children.at(0).str(true);
					ST += " ^ 3 )";
				} else {
					ST = this->node->getSymbol();
					ST += "( ";
					ST += this->children.at(0).str(true);
					ST += " )";
				}}
				break;
			  case 2:
				if(this->node->isType() == 'o' || 
						this->node->getSymbol().compare("POW") == 0) {
					ST = "( ";
					ST += this->children.at(0).str(true);
					ST += " ";
					ST += (this->node->getSymbol().compare("POW") == 0) 
						? "^" : this->node->getSymbol();
					ST += " ";
					ST += this->children.at(1).str(true);
					ST += " )";
				} else {
					ST = this->node->getSymbol();
					ST += "( ";
					ST += this->children.at(0).str(true);
					ST += ", ";
					ST += this->children.at(1).str(true);
					ST += " )";
				}
				break;
			  default:
				ST = "( ";
				for(size_t i = 0; i < this->node->getArity(); i++) {
					ST += this->children.at(i).str(true);
					ST += (i == this->node->getArity() - 1) ? "" : 
						" " + this->node->getSymbol() + " ";
				}
				ST += " )";
			}
		} else {
			ST = (this->node->getArity() == 0) ? "" : "( ";
			ST += this->node->getSymbol() + " ";
			for(size_t i = 0; i < this->node->getArity(); i++)
				ST += this->children.at(i).str(false);
			ST += (this->node->getArity() == 0) ? "" : ") ";
		}
	}
	return ST;
}
void AAST::show(const bool &infix) const {
	cout << this->str(infix) << endl;
}
#ifndef CEXPAAST_H_
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
		default:
			cerr << "AAST == does not recognize this type!" << endl;
			chk = false;
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
#endif
bool AAST::operator!=(const AAST &other) const {
	return !(this->operator==(other));
}
vector<AAST> AAST::operator<<(const double &other) {
	return this->operator<<(static_cast<AAST>(other));
}

vector<AAST> AAST::operator<<(const AAST &other) {
	vector<AAST> aux;
	aux.push_back(*this);
	aux.push_back(other);
	return aux;
}
vector<AAST> AAST::operator<<(const vector<AAST> &other) {
	vector<AAST> aux;
	aux.push_back(*this);
	for(size_t i = 0; i < other.size(); i++)
		aux.push_back(other.at(i));
	return aux;
}
void AAST::operator+=(const double &value) {
	CNS operand(value);
	this->operator+=(operand.T());
}
void AAST::operator+=(const AAST &operand) {
	OPR operation('+');
	vector<AAST> arg;
	arg.push_back(*this);
	arg.push_back(operand);
	AAST result = operation.T(arg);
	this->safe_swap(result,NULL);
}
AAST AAST::operator+(const double &value) {
	AAST result(*this);
	CNS operand(value);
	result += operand.T();
	return result;
}
AAST AAST::operator+(const AAST &operand) {
	AAST result(*this);
	result += operand;
	return result;
}

void AAST::operator-=(const double &value) {
	CNS operand(value);
	this->operator-=(operand.T());
}
void AAST::operator-=(const AAST &operand) {
	OPR operation('-');
	vector<AAST> arg;
	arg.push_back(*this);
	arg.push_back(operand);
	AAST result = operation.T(arg);
	this->safe_swap(result,NULL);
}
AAST AAST::operator-(const double &value) {
	AAST result(*this);
	CNS operand(value);
	result -= operand.T();
	return result;
}
AAST AAST::operator-(const AAST &operand) {
	AAST result(*this);
	result -= operand;
	return result;
}

void AAST::operator*=(const double &value) {
	CNS operand(value);
	this->operator*=(operand.T());
}
void AAST::operator*=(const AAST &operand) {
	OPR operation('*');
	vector<AAST> arg;
	arg.push_back(*this);
	arg.push_back(operand);
	AAST result = operation.T(arg);
	this->safe_swap(result,NULL);
}
AAST AAST::operator*(const double &value) {
	AAST result(*this);
	CNS operand(value);
	result *= operand.T();
	return result;
}
AAST AAST::operator*(const AAST &operand) {
	AAST result(*this);
	result *= operand;
	return result;
}
void AAST::operator/=(const double &value) {
	CNS operand(value);
	this->operator/=(operand.T());
}
void AAST::operator/=(const AAST &operand) {
	OPR operation('/');
	vector<AAST> arg;
	arg.push_back(*this);
	arg.push_back(operand);
	AAST result = operation.T(arg);
	this->safe_swap(result,NULL);
}
AAST AAST::operator/(const double &value) {
	AAST result(*this);
	CNS operand(value);
	result /= operand.T();
	return result;
}
AAST AAST::operator/(const AAST &operand) {
	AAST result(*this);
	result /= operand;
	return result;
}
//======================================================================
// simplify
//======================================================================
void AAST::simp_DivisionToProduct() {
	if((this->node->isType() == 'o') ? 
		this->node->getSymbol().compare("/") == 0 : false)
	{
		vector<AAST> args_prod, arg_power;
		FNC power("POW",2);
		args_prod <<= this->children.front();
		for(size_t i = 1; i < this->children.size(); i++)
			arg_power <<= this->children.at(i);
		args_prod <<= power.T(arg_power << -1.);
		OPR newop("*",args_prod.size());
		AAST aux = newop.T(args_prod);
		this->safe_swap(aux,this->parent);
	}
	for(size_t i = 0; i < this->children.size(); i++)
		this->children.at(i).simp_DivisionToProduct();
}
void AAST::simp_DifferenceToSum() {
	if((this->node->isType() == 'o') ? 
		this->node->getSymbol().compare("-") == 0 : false)
	{
		vector<AAST> args_sum, arg_minus;
		args_sum <<= this->children.front();
		arg_minus <<= -1.;
		for(size_t i = 1; i < this->children.size(); i++)
			arg_minus <<= this->children.at(i);
		OPR opaux("*",arg_minus.size());
		args_sum <<= opaux.T(arg_minus);
		OPR newop("+",args_sum.size());
		AAST aux = newop.T(args_sum);
		this->safe_swap(aux,this->parent);
	}
	for(size_t i = 0; i < this->children.size(); i++)
		this->children.at(i).simp_DifferenceToSum();
}
bool AAST::simp_JoinProducts() {
	bool perform = false;
	if((this->node->isType() == 'o') ? 
		this->node->getSymbol().compare("*") == 0 : false)
	{
		bool joined = false;
		vector<AAST> args_prod;
		for(size_t i = 0; i < this->children.size(); i++)
		{
			if((this->children.at(i).node->isType() == 'o') ? 
				this->children.at(i).node->getSymbol().compare("*") 
				== 0 : false)
			{
				joined = true;
				args_prod <<= this->children.at(i).children;
			}
			else
				args_prod <<= this->children.at(i);
		}
		if(joined)
		{
			OPR newop("*",args_prod.size());
			AAST aux = newop.T(args_prod);
			this->safe_swap(aux,this->parent);
			perform = true;
		}
	}
	for(size_t i = 0; i < this->children.size(); i++)
		perform |= this->children.at(i).simp_JoinProducts();
	return perform;
}
bool AAST::simp_JoinSums() {
	bool perform = false;
	if((this->node->isType() == 'o') ? 
		this->node->getSymbol().compare("+") == 0 : false)
	{
		bool joined = false;
		vector<AAST> args_sum;
		for(size_t i = 0; i < this->children.size(); i++)
		{
			if((this->children.at(i).node->isType() == 'o') ? 
				this->children.at(i).node->getSymbol().compare("+") 
				== 0 : false)
			{
				joined = true;
				args_sum <<= this->children.at(i).children;
			}
			else
				args_sum <<= this->children.at(i);
		}
		if(joined)
		{
			OPR newop("+",args_sum.size());
			AAST aux = newop.T(args_sum);
			this->safe_swap(aux,this->parent);
			perform = true;
		}
	}
	for(size_t i = 0; i < this->children.size(); i++)
		perform |= this->children.at(i).simp_JoinProducts();
	return perform;
}
bool AAST::simp_Constants() {
	bool perform = false;
	if(this->node->isType() == 'o' || this->node->isType() == 'f')
	{
		double result;
		AAST newnode;
		string symbol = this->node->getSymbol();
		if(this->node->getArity() == 1 || this->node->getArity() == 2)
		{
			bool allCNS = (this->children.front().node->isType() == 'c')
				&& (this->children.back().node->isType() == 'c');
			if(allCNS) {
				double arg1 = static_cast<CNS*>(this->children.front().
						      node)->value, 
					   arg2 = static_cast<CNS*>(this->children.back().
						      node)->value;
				if(symbol.compare("EXP") == 0) { 
					result = exp(arg1); 
				} else {
				if(symbol.compare("LOG") == 0) {
					result = (arg1 <= 0)? 1. : log(arg1); 
				} else {
				if(symbol.compare("LOG10") == 0) {
					result = (arg1 <= 0)? 1. : log10(arg1); 
				} else {
				if(symbol.compare("SQRT") == 0) {
					result = (arg1 < 0)? 1. : sqrt(arg1); 
				} else {
				if(symbol.compare("COS") == 0) {
					result = cos(arg1);
				} else {
				if(symbol.compare("SIN") == 0) {
					result = sin(arg1);
				} else {
				if(symbol.compare("TAN") == 0) {
					result = tan(arg1);
				} else {
				if(symbol.compare("ACOS") == 0) {
					result = (arg1 < -1. || arg1 > 1.)? 1. : acos(arg1); // protected 
				} else {
				if(symbol.compare("ASIN") == 0) {
					result = (arg1 < -1. || arg1 > 1.)? 1. : asin(arg1); // protected 
				} else {
				if(symbol.compare("ATAN") == 0) {
					result = atan(arg1);
				} else {
				if(symbol.compare("COSH") == 0) {
					result = cosh(arg1);
				} else {
				if(symbol.compare("SINH") == 0) {
					result = sinh(arg1);
				} else {
				if(symbol.compare("TANH") == 0) {
					result = tanh(arg1);
				} else {
				if(symbol.compare("POW") == 0) {
					result = pow(arg1,arg2); 
				} else {
				if(symbol.compare("LOGN") == 0) {
					result = (arg1 <= 0 || arg2 <= 0)? 1. : log(arg2)/
					log(arg1); // protected 
				} else {
				if(symbol.compare("ATAN2") == 0) {
					result = atan2(arg1,arg2);
				} else {
				if(symbol.compare("+") == 0) {
					result = arg1 + arg2;
				} else {
				if(symbol.compare("-") == 0) {
					result = arg1 - arg2;
				} else {
				if(symbol.compare("*") == 0) {
					result = arg1 * arg2;
				} else {
				if(symbol.compare("/") == 0) {
					result = (arg2 == 0)? 1. : arg1 / arg2; // protected 
				} else {
				if(symbol.compare("POW2") == 0) {
					result = arg1*arg1;
				} else {
				if(symbol.compare("POW3") == 0) {
					result = arg1*arg1*arg1;
				} else {
				if(symbol.compare("CBRT") == 0) {
					result =  pow(arg1, 1./3.);
				} else {
				if(symbol.compare("EXPn") == 0) {
					result = exp(-arg1);
				} else {
				if(symbol.compare("INV") == 0) {
					result = 1./arg1;
				}  else {
				if(symbol.compare("NEG") == 0) {
					result = -arg1;
				} else {
					cerr << "[Error] In " << __FILE__ << " at line " << 
						__LINE__ << endl;
					string msg = "Unknown function/operator: " + symbol;
					throw std::invalid_argument(msg.c_str());
				}}}}}}}}}}}}}}}}}}}}}}}}}}
				newnode = (std::isnan(result) || my::isinf(result)) ?	
					1. : ((my::isLikelyZero(result)) ? 0. : result); // protection
				//newnode = result;
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		else
		{
			vector<AAST> yesCNS, nonCNS;
			for(size_t i = 0; i < this->children.size(); i++)
			{
				if(this->children.at(i).node->isType() == 'c')
					yesCNS.push_back(this->children.at(i));
				else
					nonCNS.push_back(this->children.at(i));
			}
			if(yesCNS.size() >= 2)
			{
				if(symbol.compare("+") == 0) {
					result = 0.;
					for(size_t i = 0; i < yesCNS.size(); i++)
						result += static_cast<CNS*>(yesCNS.at(i).node)
							->value;
				} else {
				if(symbol.compare("*") == 0) {
					result = 1.;
					for(size_t i = 0; i < yesCNS.size(); i++)
						result *= static_cast<CNS*>(yesCNS.at(i).node)
							->value;
				} else {
					cerr << "[Error] In " << __FILE__ << " at line " << 
						__LINE__ << endl;
					string msg = "Unknown function/operator: " + symbol;
					throw std::invalid_argument(msg.c_str());
				}}
				if(nonCNS.empty())
					newnode = (std::isnan(result) || my::isinf(result)) ?	
					1. : ((my::isLikelyZero(result)) ? 0. : result); // protection
				else
				{
					OPR newop(symbol,nonCNS.size()+1);
					newnode = newop.T(result << nonCNS);
				}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
	}
	if(!this->children.empty())
		for(size_t i = 0; i < this->children.size(); i++)
			perform |= this->children.at(i).simp_Constants();
	return perform;
}
bool AAST::simp_Distributive() {
	bool perform = false;
	if((this->node->isType() == 'o') ? this->node->getSymbol().
		compare("*") == 0 : false)
	{
		AAST aux;
		bool found_factor = false, found_addition = false;
		size_t factor, addition;
		for(size_t i = 0; i < this->children.size(); i++)
		{
			if(!found_factor && (this->children.at(i).node->isType() == 
				'c' || this->children.at(i).node->isType() == 'v'))
			{
				factor = i;
				found_factor = true;
			}
			if((this->children.at(i).node->isType() == 'o') ? 
				this->children.at(i).node->getSymbol().compare("+") 
				== 0 : false)
			{
				addition = i;
				found_addition = true;
			}
			if(found_factor && found_addition)
			{
				for(size_t j = 0; j < this->children.at(addition).
					children.size(); j++)
				{
					aux = this->children.at(factor)*this->children.
						at(addition).children.at(j);
					this->children.at(addition).children.at(j).
						safe_swap(aux,this);
				}
				aux = 1.;
				this->children.at(factor).safe_swap(aux,this);
				perform = true;
				found_factor = false;
				found_addition = false;
			}
		}
	}
	if(!this->children.empty())
		for(size_t i = 0; i < this->children.size(); i++)
			perform |= this->children.at(i).simp_Distributive();
	return perform;
}
bool AAST::simp_Identities() {
	bool perform = false;
	if(this->node->isType() == 'o' || ((this->node->isType() == 'f') ? 
		this->node->getSymbol().compare("POW") == 0 : false))
	{
		AAST newnode, pattern;
		vector<AAST> yesvec, novec, auxvec, auxvec2;
		vector<size_t> indices;
		FNC newpow("POW",2);
		size_t count;
		if(this->node->getSymbol().compare("+") == 0) { //::::::::::::::
			// Additive identity (x + 0 = x) ...........................
			yesvec.clear();
			novec.clear();
			for(size_t i = 0; i < this->children.size(); i++)
				if((this->children.at(i).node->isType() == 'c') ?
						my::isLikelyZero(static_cast<CNS*>(
						this->children.at(i).node)->value) : false)
					yesvec <<= static_cast<CNS*>(this->children.at(i).
						node)->value;
				else
					novec <<= this->children.at(i);
			if(!yesvec.empty())
			{
				if(novec.empty()) {
					newnode = 0.;
				} else {
				if(novec.size() == 1) {
					newnode = novec.front();
				} else {
					OPR newop(this->node->getSymbol(),novec.size());
					newnode = newop.T(novec); 
				}}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("+") == 0) { //::::::::::::::
			// Additive inverse (x - x = 0) ............................
			indices.clear();
			for(size_t i = 0; i < this->children.size()-1; i++)
			{
				if (!my::is_in(i,indices))
				{
					pattern = -1.*this->children.at(i);
					for(size_t j = i+1; j < this->children.size(); j++)
					{
						if(this->children.at(j) == pattern)
						{
							indices.push_back(i);
							indices.push_back(j);
							break;
						}
					}
				}
			}
			if(!indices.empty())
			{
				novec.clear();
				for(size_t i = 0; i < this->children.size(); i++)
					if (!my::is_in(i,indices))
						novec <<= this->children.at(i);
				if(novec.empty()) {
					newnode = 0.;
				} else {
				if(novec.size() == 1) {
					newnode = novec.front();
				} else {
					OPR newop(this->node->getSymbol(),novec.size());
					newnode = newop.T(novec); 
				}}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("+") == 0) { //::::::::::::::
			// Distributive type I (x+x+...x = n*x) ....................
			yesvec.clear();
			indices.clear();
			for(size_t i = 0; i < this->children.size()-1; i++)
			{
				count = 1;
				indices.push_back(i);
				pattern = this->children.at(i);
				for(size_t j = i+1; j < this->children.size(); j++)
				{
					if (!my::is_in(j,indices))
						if(this->children.at(j) == pattern)
						{
							count++;
							indices.push_back(j);
						}
				}
				if(count > 1)
					yesvec <<= (double)count*pattern;
				else
					indices.pop_back();
			}
			if(!yesvec.empty())
			{
				novec.clear();
				for(size_t i = 0; i < this->children.size(); i++)
					if (!my::is_in(i,indices))
						novec <<= this->children.at(i);
				if(novec.empty()) {
					if(yesvec.size() > 1)
					{
						OPR newop(this->node->getSymbol(),yesvec.
							size());
						newnode = newop.T(yesvec);
					}
					else
						newnode = yesvec.front();
				} else {
					OPR newop(this->node->getSymbol(),yesvec.size()+
						novec.size());
					newnode = newop.T(yesvec << novec); 
				}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		} 
		////============================================================
		if(this->node->getSymbol().compare("*") == 0) { //::::::::::::::
			// Zero element (x*0 = 0) ..................................
			yesvec.clear();
			novec.clear();
			for(size_t i = 0; i < this->children.size(); i++)
				if((this->children.at(i).node->isType() == 'c') ?
						my::isLikelyZero(static_cast<CNS*>(
						this->children.at(i).node)->value) : false)
					yesvec <<= static_cast<CNS*>(this->children.at(i).
						node)->value;
				else
					novec <<= this->children.at(i);
			if(!yesvec.empty())
			{
				newnode = 0.;
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("*") == 0) { //::::::::::::::
			// Multiplicative identity (x*1 = x) .......................
			yesvec.clear();
			novec.clear();
			for(size_t i = 0; i < this->children.size(); i++)
				if((this->children.at(i).node->isType() == 'c') ?
						my::isLikelyZero(1.-static_cast<CNS*>(
						this->children.at(i).node)->value) : false)
					yesvec <<= static_cast<CNS*>(this->children.at(i).
						node)->value;
				else
					novec <<= this->children.at(i);
			if(!yesvec.empty())
			{
				if(novec.empty()) {
					newnode = 1.;
				} else {
				if(novec.size() == 1) {
					newnode = novec.front();
				} else {
					OPR newop(this->node->getSymbol(),novec.size());
					newnode = newop.T(novec); 
				}}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("*") == 0) { //::::::::::::::
			// Multiplicative inverse (x*x^(-1) = 1) ...................
			indices.clear();
			for(size_t i = 0; i < this->children.size()-1; i++)
			{
				if (!my::is_in(i,indices))
				{
					pattern = newpow.T(this->children.at(i) << -1.);
					for(size_t j = i+1; j < this->children.size(); j++)
					{
						if(this->children.at(j) == pattern)
						{
							indices.push_back(i);
							indices.push_back(j);
							break;
						}
					}
				}
			}
			if(!indices.empty())
			{
				novec.clear();
				for(size_t i = 0; i < this->children.size(); i++)
					if (!my::is_in(i,indices))
						novec <<= this->children.at(i);
				if(novec.empty()) {
					newnode = 1.;
				} else {
				if(novec.size() == 1) {
					newnode = novec.front();
				} else {
					OPR newop(this->node->getSymbol(),novec.size());
					newnode = newop.T(novec); 
				}}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("*") == 0) { //::::::::::::::
			// Exponentiation type I (x*x*...x = x^n) ..................
			yesvec.clear();
			for(size_t i = 0; i < this->children.size()-1; i++)
			{
				count = 1;
				indices.push_back(i);
				pattern = this->children.at(i);
				for(size_t j = i+1; j < this->children.size(); j++)
				{
					if (!my::is_in(j,indices))
						if(this->children.at(j) == pattern) 
						{
							count++;
							indices.push_back(j);
						}
				}
				if(count > 1)
					yesvec <<= newpow.T(pattern << (double)count);
				else
					indices.pop_back();
			}
			if(!yesvec.empty())
			{
				novec.clear();
				for(size_t i = 0; i < this->children.size(); i++)
					if (!my::is_in(i,indices))
						novec <<= this->children.at(i);
				if(novec.empty()) {
					if(yesvec.size() > 1)
					{
						OPR newop(this->node->getSymbol(),
							yesvec.size());
						newnode = newop.T(yesvec);
					}
					else
						newnode = yesvec.front();
				} else {
					OPR newop(this->node->getSymbol(),yesvec.size()+
						novec.size());
					newnode = newop.T(yesvec << novec); 
				}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		} 
		if(this->node->getSymbol().compare("*") == 0) { //::::::::::::::
			// Exponentiation type II (x^a*x^b*...x^c = x^(a+b+...c)) ..
			yesvec.clear();
			indices.clear();
			for(size_t i = 0; i < this->children.size()-1; i++)
			{
				if((this->children.at(i).node->isType() == 'f') ? 
					this->children.at(i).node->getSymbol().
					compare("POW") == 0 : false)
				{
					auxvec.clear();
					indices.push_back(i);
					pattern = this->children.at(i).children.front();
					auxvec <<= this->children.at(i).children.back();
					for(size_t j = i+1; j < this->children.size(); j++)
					{
						if((this->children.at(j).node->isType() == 'f') 
							? this->children.at(j).node->getSymbol().
							compare("POW") == 0 : false)
						{
							if (!my::is_in(j,indices))
								if(this->children.at(j).children.front()
									== pattern)
								{
									indices.push_back(j);
									auxvec <<= this->children.at(j).
										children.back();
								}
						}
					}
					if(auxvec.size() > 1)
					{
						OPR newsum("+",auxvec.size());
						yesvec <<= newpow.T(pattern << 
							newsum.T(auxvec));
					}
					else
						indices.pop_back();
				}
			}
			if(!yesvec.empty())
			{
				novec.clear();
				for(size_t i = 0; i < this->children.size(); i++)
					if (!my::is_in(i,indices))
						novec <<= this->children.at(i);
				if(novec.empty()) {
					if(yesvec.size() > 1)
					{
						OPR newop(this->node->getSymbol(),
							yesvec.size());
						newnode = newop.T(yesvec);
					}
					else
						newnode = yesvec.front();
				} else {
					OPR newop(this->node->getSymbol(),yesvec.size()+
						novec.size());
					newnode = newop.T(yesvec << novec); 
				}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		} 
		if(this->node->getSymbol().compare("*") == 0) { //::::::::::::::
			// Exponentiation type III (x*x^a = x^(a+1)) ...............
			yesvec.clear();
			indices.clear();
			//////////////////////////////// SIZE() INSTEAD OF SIZE()-1
			for(size_t i = 0; i < this->children.size(); i++) 
			{
				auxvec.clear();
				indices.push_back(i);
				pattern = this->children.at(i);
				auxvec <<= 1.;
				//////////////////////////////////// ZERO INSTEAD OF I+1
				for(size_t j = 0; j < this->children.size(); j++) 
				{
					if (!my::is_in(j,indices))
						if((this->children.at(j).node->isType() == 'f')  
							? this->children.at(j).node->getSymbol().
							compare("POW") == 0 : false)
						{
							if(this->children.at(j).children.front()
								== pattern)
							{
								indices.push_back(j);
								auxvec <<= this->children.at(j).
									children.back();
							}
						}
				}
				if(auxvec.size() > 1)
				{
					OPR newsum("+",auxvec.size());
					yesvec <<= newpow.T(pattern << newsum.T(auxvec));
				}
				else
					indices.pop_back();
			}
			if(!yesvec.empty())
			{
				novec.clear();
				for(size_t i = 0; i < this->children.size(); i++)
					if (!my::is_in(i,indices))
						novec <<= this->children.at(i);
				if(novec.empty()) {
					if(yesvec.size() > 1)
					{
						OPR newop(this->node->getSymbol(),
							yesvec.size());
						newnode = newop.T(yesvec);
					}
					else
						newnode = yesvec.front();
				} else {
					OPR newop(this->node->getSymbol(),yesvec.size()+
						novec.size());
					newnode = newop.T(yesvec << novec); 
				}
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("POW") == 0) { //::::::::::::
			// Powers to zero (x^0 = 1, undef 0^0 here defined as 1) ...
			if((this->children.back().node->isType() == 'c') ?
					my::isLikelyZero(static_cast<CNS*>(this->
					children.back().node)->value) : false)
			{
				newnode = 1.;
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("POW") == 0) { //::::::::::::
			// Powers of zero (0^x = 0) ................................
			if((this->children.front().node->isType() == 'c') ?
					my::isLikelyZero(static_cast<CNS*>(this->
					children.front().node)->value) : false)
			{
				newnode = 0.;
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("POW") == 0) { //::::::::::::
			// Powers to one (x^1 = x) .................................
			if((this->children.back().node->isType() == 'c') ?
					my::isLikelyZero(1.-static_cast<CNS*>(this->
					children.back().node)->value) : false)
			{
				newnode = this->children.front();
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("POW") == 0) { //::::::::::::
			// Power of one (1^x = 1) ..................................
			if((this->children.front().node->isType() == 'c') ?
					my::isLikelyZero(1.-static_cast<CNS*>(this->
					children.front().node)->value) : false)
			{
				newnode = 1.;
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
		if(this->node->getSymbol().compare("POW") == 0) { //::::::::::::
			// Exponentiation power to power ((x^a)^b = x^(a*b)) .......
			if((this->children.front().node->isType() == 'f') ?
					this->children.front().node->getSymbol().
					compare("POW") == 0 : false)
			{
				OPR newprod("*",2);
				newnode = newpow.T(this->children.front().children.
					front() << newprod.T(this->children.front().
					children.back() << this->children.back()));
				this->safe_swap(newnode,this->parent);
				perform = true;
			}
		}
	}
	for(size_t i = 0; i < this->children.size(); i++)
		perform |= this->children.at(i).simp_Identities();
	return perform;
}
void AAST::simplify() {
	bool performed, performed2;
	//int count = 0;
	this->simp_DivisionToProduct();
	this->simp_DifferenceToSum();
	do {
		performed = false;
		performed |= this->simp_JoinProducts();
		performed |= this->simp_JoinSums();
		performed |= this->simp_Distributive();
		performed |= this->simp_Identities();
		do {
			performed2 = false;
			performed2 |= this->simp_Constants();
		} while(performed2);
		performed |= performed2;
	} while(performed);
}
void AAST::assign2variables(const size_t &qty, const double &value) {
	vector<double> vvalue(qty,value);
	this->assign2variables(vvalue);
}
void AAST::assign2variables(const double &value) {
	vector<double> vvalue(1,value);
	this->assign2variables(vvalue);
}
void AAST::assign2variables(const vector<double> &values) {
	if(this->node->isType() == 'v')
	{
		size_t index = static_cast<VAR*>(this->node)->index;
		CNS value(values.at(index));
		AAST newCNS = value.T();
		this->safe_swap(newCNS,this->parent);
	}
	else
		if(!this->children.empty())
			for(size_t i = 0; i < this->children.size(); i++)
				this->children.at(i).assign2variables(values);
}
double AAST::getRealValue() const {
	if(this->node->isType() == 'c')
		return static_cast<CNS*>(this->node)->value;
	else
		return NAN;
}
double AAST::eval_with_variables(const vector<double> &values) const {
	AAST tree(*this);
	tree.assign2variables(values);
	tree.simplify();
	tree.simplify();
	tree.simplify();
	cout << " evalw/var ";
	tree.show();
	return tree.getRealValue();
}
void AAST::map_variable_to(const size_t &index, const AAST &newvar) 
{
	if((this->node->isType() == 'v') ? index == static_cast<VAR*>
		(this->node)->index : false)
	{
		AAST volvar(newvar);
		this->safe_swap(volvar,this->parent);
	}
	else
	{
		if(!this->children.empty())
		{
			for(size_t i = 0; i < this->children.size(); i++)
				this->children.at(i).map_variable_to(index, newvar);
		}
	}
}
// =====================================================================
// Node AASTree representation  ----------------------------------------
#ifndef CEXPAAST_H_
AAST Node::T() const {
	AAST tree;
	switch(this->isType()) {
	  case 'c':
		tree.node = new CNS(static_cast<const CNS*>(this)->value);
		break;
	  case 'v':
		tree.node = new VAR(static_cast<const VAR*>(this)->index);
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
#endif
AAST Node::T(const vector<AAST> &arg) const {
	AAST tree;
	char type = this->isType();
	if(this->getArity() != arg.size()) {
		cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
			<< endl;
		string msg = "Review argument size for: " + this->getSymbol();
		throw std::invalid_argument(msg.c_str());
	}
	switch(type) {
	  case 'f':
		tree.node = new FNC(static_cast<const FNC*>(this)->symbol,
			static_cast<const FNC*>(this)->arity);
		break;
	  case 'o':
		tree.node = new OPR(static_cast<const OPR*>(this)->symbol,
			static_cast<const OPR*>(this)->arity);
		break;
	  default:
		cerr << "[Error] In " << __FILE__ << " at line " << __LINE__ 
			<< endl;
		string msg = "Review arguments for: " + this->getSymbol();
		throw std::invalid_argument(msg.c_str());
	}
	tree.children = arg;
	for(size_t i = 0; i < arg.size(); i++)
		tree.children.at(i).parent = &tree;
	return tree;
}
// Friends -------------------------------------------------------------
AAST operator+(const double &value, const AAST &rhs) {
	CNS lhs(value);
	return lhs.T() + rhs;
}
AAST operator-(const double &value, const AAST &rhs) {
	CNS lhs(value);
	return lhs.T() - rhs;
}
AAST operator*(const double &value, const AAST &rhs) {
	CNS lhs(value);
	return lhs.T() * rhs;
}
AAST operator/(const double &value, const AAST &rhs) {
	CNS lhs(value);
	return lhs.T() / rhs;
}
vector<AAST> operator<<(const vector<AAST> &arg1, const AAST &arg2) {
	vector<AAST> aux(arg1);
	aux.push_back(arg2);
	return aux;
}
vector<AAST> operator<<(const double &arg1, const AAST &arg2) {
	vector<AAST> aux;
	aux.push_back(static_cast<AAST>(arg1));
	aux.push_back(arg2);
	return aux;
}
vector<AAST> operator<<(const double &arg1, const vector<AAST> &arg2) {
	vector<AAST> aux;
	aux.push_back(static_cast<AAST>(arg1));
	for(size_t i = 0; i < arg2.size(); i++)
		aux.push_back(arg2.at(i));
	return aux;
}
vector<AAST> operator<<(const Double &arg1, const Double &arg2) {
	AAST aux1(arg1.value);
	AAST aux2(arg2.value);
	vector<AAST> auxv;
	auxv.push_back(aux1);
	auxv.push_back(aux2);
	return auxv;
}
vector<AAST> operator<<(const vector<AAST> &arg1, 
		const vector<AAST> &arg2) {
	vector<AAST> aux(arg1);
	for(size_t i = 0; i < arg2.size(); i++)
		aux.push_back(arg2.at(i));
	return aux;
}
void operator<<=(vector<AAST> &arg1, const AAST &arg2) {
	arg1.push_back(arg2);
}
void operator<<=(vector<AAST> &arg1, const vector<AAST> &arg2) {
	for(size_t i = 0; i < arg2.size(); i++)
		arg1.push_back(arg2.at(i));
}
// =====================================================================
#endif
