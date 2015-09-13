/* code by: iperetta@ieee.org */
#ifndef GENPROG_H_
#define GENPROG_H_
// =====================================================================
#include "cFitness.h"
#include "cexpAAST.h"
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
using namespace std;

class GenProg;

class stTerminalSet 
{
	friend class GenProg;
  private:
	void include_vars(const size_t &upto);
  public:
	vector<CNS> constants;
	vector<VAR> variables;
	stTerminalSet() {};
	~stTerminalSet() {};
	void clear() { constants.clear(); variables.clear(); };
	size_t size() { return constants.size()+variables.size(); };
	bool isempty() { return constants.empty() && variables.empty(); };
	void include(const double &newone);
	void include(const double *newone, const size_t &n);
	void include(const vector<double> &newones);
	void include_basic();
	void include_common();
	string get_symbol(const size_t &i_);
	AAST get(const size_t &i_);
};

void stTerminalSet::include(const double &newone)
{
	constants.push_back(static_cast<CNS>(newone)); 
}

void stTerminalSet::include(const double *newone, const size_t &n)
{
	for(size_t i = 0; i < n; i++)
		this->include(newone[i]);
}

void stTerminalSet::include(const vector<double> &newones)
{
	for(size_t i = 0; i < newones.size(); i++)
		this->include(newones.at(i));
}

void stTerminalSet::include_basic()
{
	this->include(-1.);
	this->include(0.);
	this->include(1.);
	this->include(2.);
	this->include(3.);
}

void stTerminalSet::include_common()
{
	this->include_basic();
	this->include_basic();
	this->include_basic();
	this->include(M_PI);
}

void stTerminalSet::include_vars(const size_t &upto)
{
	for(size_t i = 0; i < upto; i++)
		variables.push_back(static_cast<VAR>(i));
}

string stTerminalSet::get_symbol(const size_t &i_)
{
	if(i_ >= this->size()) 
		throw std::out_of_range("by: TerminalSet get_symbol()");
	if(i_ >= this->constants.size())
		return variables.at(i_ - this->constants.size()).getSymbol();
	else
		return constants.at(i_).getSymbol();
}

AAST stTerminalSet::get(const size_t &i_)
{
	if(i_ >= this->size()) 
		throw std::out_of_range("by: TerminalSet get()");
	if(i_ >= this->constants.size())
		return variables.at(i_ - this->constants.size()).T();
	else
		return constants.at(i_).T();
}

class stFuncionSet 
{
	friend class GenProg;
  public:
	vector<OPR> operators;
	vector<FNC> functions;
	stFuncionSet() {};
	~stFuncionSet() {};
	void clear() { operators.clear(); functions.clear(); };
	size_t size(const bool &allowfunc = true);
	bool isempty() { return operators.empty() && functions.empty(); };
	void include(const OPR &newone);
	void include_ops();
	void include(const FNC &newone);
	void include_basic();
	void include_common();
	void include_trigonometric();
	char type(const size_t &i_);
	size_t arity(const size_t &i_);
	string get_symbol(const size_t &i_);
	AAST get(const size_t &i_, const AAST &arg);
	AAST get(const size_t &i_, const vector<AAST> &arg);
};

size_t stFuncionSet::size(const bool &allowfunc) 
{ 
	return operators.size() + ((allowfunc) ? functions.size() : 0); 
}

void stFuncionSet::include(const OPR &newone)
{
	operators.push_back(newone);
}

void stFuncionSet::include_ops()
{
	operators.push_back(static_cast<OPR>("+"));
	operators.push_back(static_cast<OPR>("-"));
	operators.push_back(static_cast<OPR>("*"));
	operators.push_back(static_cast<OPR>("/"));
}
	
void stFuncionSet::include(const FNC &newone)
{
	functions.push_back(newone);
}

void stFuncionSet::include_basic()
{
	this->include_ops();
}

void stFuncionSet::include_common()
{
	this->include_basic();
	this->include_basic();
	this->include_basic();
	FNC f1("INV",1), f2("POW3",1), 
		f3("EXPn",1), f4("SQRT",1), 
		f5("CBRT",1), f6("POW2",1), 
		f7("NEG",1);	
	this->include(f1);
	this->include(f2);
	this->include(f3);
	this->include(f4);
	this->include(f5);
	this->include(f6);
	this->include(f7);
}

void stFuncionSet::include_trigonometric()
{
	FNC f1("SIN",1), f2("COS",1), 
		f3("TAN",1);	
	this->include(f1);
	this->include(f2);
	this->include(f3);
}
char stFuncionSet::type(const size_t &i_)
{
	if(i_ >= this->size()) 
		throw std::out_of_range("by: FunctionSet type()");
	if(i_ >= this->operators.size())
		return functions.at(i_ - this->operators.size()).isType();
	else
		return operators.at(i_).isType();
}

size_t stFuncionSet::arity(const size_t &i_)
{
	if(i_ >= this->size()) 
		throw std::out_of_range("by: FunctionSet arity()");
	if(i_ >= this->operators.size())
		return functions.at(i_ - this->operators.size()).getArity();
	else
		return operators.at(i_).getArity();
}

string stFuncionSet::get_symbol(const size_t &i_)
{
	if(i_ >= this->size()) 
		throw std::out_of_range("by: FunctionSet get_symbol()");
	if(i_ >= this->operators.size())
		return functions.at(i_ - this->operators.size()).getSymbol();
	else
		return operators.at(i_).getSymbol();
}

AAST stFuncionSet::get(const size_t &i_, const AAST &arg)
{
	if(i_ >= this->size()) 
		throw std::out_of_range("by: FunctionSet get()");
	if(i_ >= this->operators.size())
		return functions.at(i_ - this->operators.size()).T(arg);
	else
		return operators.at(i_).T(arg);
}

AAST stFuncionSet::get(const size_t &i_, const vector<AAST> &arg)
{
	if(i_ >= this->size()) 
		throw std::out_of_range("by: FunctionSet get()");
	if(i_ >= this->operators.size())
		return functions.at(i_ - this->operators.size()).T(arg);
	else
		return operators.at(i_).T(arg);
}

class GenProg 
{
	friend class stTerminalSet;
	friend class stFuncionSet;
  public:
	size_t nPop;
	size_t nGen;
	size_t nTourn;
	double prXO;
	double prMT;
	size_t bestID;
	size_t nVar;
	size_t nAlleles;
	size_t tree_depth;
	double current_stdFit;
	myFitness FitControl;
	stTerminalSet TerminalSet;
	stFuncionSet FunctionSet;
	vector<vector<AAST> > Population;
	vector<double> Fitness;
	vector<vector<AAST> > Offspring;
	vector<double> OffFitness;
	vector<double> meanFit_vec;
	vector<double> bestFit_vec;
	GenProg(const size_t &nPop_ = 100, const size_t &nGen_ = 70, const 
		size_t &nTourn_ = 3, const double &prXO_ = 0.9, const double 
		&prMT_ = 0.15) : nPop(nPop_), nGen(nGen_), nTourn(nTourn_),
		prXO(prXO_), prMT(prMT_), bestID(0), nVar(0), nAlleles(0), 
		tree_depth(3), current_stdFit(HUGE_VAL)
		{ 
			this->Population.resize(nPop_);
			this->Fitness.resize(nPop_);
			this->Offspring.resize(nPop_);
			this->OffFitness.resize(nPop_);
		}
	~GenProg() {};
	void init(const string &namefile, const size_t &poly_deg, const 
		size_t &diff_ord, const size_t &power_of_2 = 9);
	void sets_basic();
	void sets_common();
	AAST buildTree(int max_d, const bool &method, bool allowfunc = 
		true);
	void initializePopulation(const size_t &tree_depth = 3);
	vector<double> evalFitness(const vector<vector<AAST> > &pop);
	size_t tournament();
	size_t selectNode(const AAST &tree, const bool &onlyfunc_ = true);
	vector<AAST> recombination(const size_t &idPar1, const size_t 
		&idPar2);
	vector<AAST> crossover(const size_t &idPar1, const size_t &idPar2);
	vector<AAST> bartering(const size_t &idPar1, const size_t &idPar2);
	void mutation(const size_t &ID, vector<AAST> &individual);
	void permutation(const size_t &ID1, vector<AAST> &individual);
	void pruning(const size_t &ID, vector<AAST> &individual);
	void generateOffspring();
	void fitnessPopulation();	
	void fitnessOffspring();
	void gatherData();
	void reducePopulation(const bool &keep_best = false);
	void show_sets();
	void show(const size_t &which);
	void getTGEC(const size_t &which);
};

void GenProg::init(const string &namefile, const size_t &poly_deg, const 
	size_t &diff_ord, const size_t &power_of_2)
{
	if(TerminalSet.isempty() || FunctionSet.isempty())
		throw std::invalid_argument(
			"Either terminal set or function set is empty!");
	this->FitControl.init(namefile, poly_deg, diff_ord, power_of_2);
	this->nVar = this->FitControl.getNVars();
	this->nAlleles = this->FitControl.getDifOrds().size() + 1;
	meanFit_vec.clear();
	bestFit_vec.clear();
	this->TerminalSet.variables.clear();
	this->TerminalSet.include_vars(this->nVar);
}

void GenProg::sets_basic()
{
	this->TerminalSet.include_basic();
	this->FunctionSet.include_basic();
}

void GenProg::sets_common()
{
	this->TerminalSet.include_common();
	this->FunctionSet.include_common();
}

AAST GenProg::buildTree(int max_d, const bool &method, bool allowfunc_)
{
	bool allowfunc = true;
	size_t TSs = TerminalSet.size(),
		   FSs = FunctionSet.size();
	if(max_d  == 0 || ( method && PRNG.real() < 
		(double)TSs/((double)TSs+FSs) ))
	{
		return TerminalSet.get(PRNG.integer(0,TSs-1));
	}
	else
	{
		size_t chosen = PRNG.integer(0,FunctionSet.size(allowfunc)-1),
		       arity = FunctionSet.arity(chosen);
		if(FunctionSet.type(chosen) == 'f')
			allowfunc = false;
		if(arity == 1)
			return FunctionSet.get(chosen, buildTree(max_d-1, method, 
				allowfunc));
		else
		{
			vector<AAST> arg(arity);
			for(size_t i = 0; i < arity; i++)
				arg.at(i) = buildTree(max_d-1, method, allowfunc);
			return FunctionSet.get(chosen, arg);
		}
	}
}

void GenProg::initializePopulation(const size_t &tree_depth_)
{
	this->tree_depth = tree_depth_;
	vector<AAST> individual(this->nAlleles);
	for(size_t i = 0; i < this->nPop; i++)
	{
		for(size_t j = 0; j < this->nAlleles; j++)
			individual.at(j) = buildTree(this->tree_depth, (PRNG.real() 
				< 0.5)? true : false);
		Population.at(i) = individual;
	}
}

vector<double> GenProg::evalFitness(const vector<vector<AAST> > &pop)
{
	vector<double> fitness(this->nPop);
	AAST s;
	vector<AAST> K(this->nAlleles - 1);
	for(size_t i = 0; i < this->nPop; i++)
	{
		//cout << ">> =============================== " << endl;
		//cout << ">> .:: Individual # " << i+1 << " ::." << endl;
		s = pop.at(i).at(0);
		s.simplify();
		//cout << "allele # 0 (sx): " << s.str() << endl;
		for(size_t j = 0; j < this->nAlleles - 1; j++)
		{
			K.at(j) = pop.at(i).at(j+1);
			K.at(j).simplify();
			//cout << "allele # " << j+1 << ": " << K.at(j).str() << endl;
		}
		//this->show(i);
		fitness.at(i) = FitControl.solve_for(K, s, false);
		//cout << ">> Fitness: " << fitness.at(i) << endl;
		//cout << ">> =============================== " << endl;
	}
	
	return fitness;
}

size_t GenProg::tournament()
{
	size_t idwin;
	vector<size_t> idP = PRNG.randperm(this->nPop);
	idwin = idP.at(0);
	for(size_t i = 1; i < this->nTourn; i++)
		if(this->Fitness.at(idP.at(i)) < this->Fitness.at(idwin))
			idwin = idP.at(i);
	return idwin;
}

size_t GenProg::selectNode(const AAST &tree, const bool &onlyfunc_)
{
	bool onlyfunc = false;
	vector<size_t> index_f = my::catsortvec<size_t>(
		tree.getTypeIdx("o"), tree.getTypeIdx("f"));
	vector<size_t> index_t = my::catsortvec<size_t>(tree.getTypeIdx("c")
		, tree.getTypeIdx("v"));
	if((PRNG.real() < 0.9 || onlyfunc) && !index_f.empty())
		return index_f.at(PRNG.integer(0,index_f.size()-1));
	else
		return index_t.at(PRNG.integer(0,index_t.size()-1));
}

vector<AAST> GenProg::recombination(const size_t &idPar1, const size_t 
	&idPar2)
{
	vector<AAST> offspr1(this->Population.at(idPar1));
	vector<AAST> offspr2(this->Population.at(idPar2));
	AAST temp;
	size_t oneID = PRNG.integer(0,this->nAlleles-1);
	size_t id1, id2;
	for(size_t i = 0; i < this->nAlleles; i++)
		if(PRNG.real() < 1./this->nAlleles || i == oneID)
		{
			id2 = this->selectNode(offspr2.at(i));
			temp = offspr2.at(i).at(id2);
			id1 = this->selectNode(offspr1.at(i), true);
			offspr1.at(i).at(id1).safe_swap(temp, offspr1.at(i).at(id1).
				getParent());
			offspr2.at(i).at(id2).safe_swap(temp, offspr2.at(i).at(id2).
				getParent());
		}
	return (PRNG.real() < 0.5)? offspr1 : offspr2;
}

vector<AAST> GenProg::crossover(const size_t &idPar1, const size_t 
	&idPar2)
{
	vector<AAST> offspr1(this->Population.at(idPar1));
	vector<AAST> offspr2(this->Population.at(idPar2));
	AAST temp;
	size_t oneID = PRNG.integer(0,this->nAlleles-1);
	for(size_t i = 0; i < this->nAlleles; i++)
		if(PRNG.real() < 1./this->nAlleles || i == oneID)
		{
			temp = offspr1.at(i);
			offspr1.at(i) = offspr2.at(i);
			offspr2.at(i) = temp;
		}
	return (PRNG.real() < 0.5)? offspr1 : offspr2;
}


vector<AAST> GenProg::bartering(const size_t &idPar1, const size_t 
	&idPar2)
{
	vector<AAST> offspr1(this->Population.at(idPar1));
	vector<AAST> offspr2(this->Population.at(idPar2));
	size_t ID1 = PRNG.integer(0,this->nAlleles-1);
	size_t ID2 = PRNG.integer(0,this->nAlleles-1);
	AAST aux = offspr1.at(ID1);
	offspr1.at(ID1) = offspr2.at(ID2);
	offspr2.at(ID2) = aux;
	return (PRNG.real() < 0.5)? offspr1 : offspr2;
}


void GenProg::mutation(const size_t &ID, vector<AAST> &individual)
{
	AAST temp;
	temp = buildTree(PRNG.integer(0,this->tree_depth), (PRNG.real() < 
		0.5)? true : false);
	size_t id = this->selectNode(individual.at(ID));
	individual.at(ID).at(id).safe_swap(temp, individual.at(ID).at(id).
		getParent());
}

void GenProg::permutation(const size_t &ID1, vector<AAST> &individual)
{
	size_t ID2 = PRNG.integer(0,this->nAlleles-1);
	AAST aux = individual.at(ID1);
	individual.at(ID1) = individual.at(ID2);
	individual.at(ID2) = aux;
}

void GenProg::pruning(const size_t &ID, vector<AAST> &individual)
{
	bool to_one = false;
	if(ID >= 2) // is a coefficient of a derivative ?
	{
		size_t count = 0;
		for(size_t i = 2; i < this->nAlleles; i++) // how many are zeros?
			if(my::isLikelyZero(individual.at(i).getRealValue()))
				count++;
		if (count == this->nAlleles-2)
			to_one = true;
	}
	AAST temp = (to_one)? 1. : 0.;
	size_t id = this->selectNode(individual.at(ID));
	individual.at(ID).at(id).safe_swap(temp, individual.at(ID).at(id).
		getParent());
}

void GenProg::generateOffspring()
{
	double chance;
	for(size_t i = 0; i < this->nPop; i++)
	{
		size_t idPar1, idPar2;
		idPar1 = tournament();
		// Recombination (xover) or Reproduction (copy)
		if(PRNG.real() < this->prXO)
		{
			do { 
				idPar2 = tournament();	
			} while(idPar2 == idPar1);
			chance = PRNG.real();
			if(chance < 1./3.)
				this->Offspring.at(i) = this->recombination(idPar1, 
					idPar2);
			else 
			if(chance < 2./3.)
				this->Offspring.at(i) = this->crossover(idPar1, idPar2);
			else
				this->Offspring.at(i) = this->bartering(idPar1, idPar2);
		}
		else
		{
			this->Offspring.at(i) = this->Population.at(idPar1);
		}
		// One parent operators
		if(PRNG.real() < this->prMT)
		{
			for(size_t j = 0; j < this->nAlleles; j++)
			{
				if(PRNG.real() < 1./this->nAlleles)
				{
					chance = PRNG.real();
					if(chance < 1./3.)
						this->mutation(j,this->Offspring.at(i));
					else
					if(chance <2./3.)
						this->permutation(j,this->Offspring.at(i));
					else
						this->pruning(j,this->Offspring.at(i));
				}
			}
		}
	}
}

void GenProg::fitnessPopulation()
{
	this->Fitness = this->evalFitness(this->Population);
	this->bestID = my::which_min(this->Fitness);
}

void GenProg::fitnessOffspring() 
{
	this->OffFitness = this->evalFitness(this->Offspring);
}

void GenProg::gatherData()
{
	this->meanFit_vec.push_back(my::protected_mean(this->Fitness));
	this->bestFit_vec.push_back(this->Fitness.at(this->bestID));
	this->current_stdFit = my::protected_std(this->Fitness);
}

void GenProg::reducePopulation(const bool &keep_best)
{
	for(size_t i = 0; i < this->nPop; i++)
	{
		if(this->OffFitness.at(i) < this->Fitness.at(i))
		{
			this->Population.at(i) = this->Offspring.at(i);
			this->Fitness.at(i) = this->OffFitness.at(i);
		}
	}
	this->bestID = my::which_min(this->Fitness);
	this->Offspring.clear();
	this->OffFitness.clear();
	this->Offspring.resize(this->nPop);
	this->OffFitness.resize(this->nPop);
}

void GenProg::show_sets()
{
	cout << "Terminal Set: ";
	for(size_t i = 0; i < this->TerminalSet.size(); i++)
		cout << this->TerminalSet.get_symbol(i) << " ";
	cout << endl;
	cout << "Function Set: ";
	for(size_t i = 0; i < this->FunctionSet.size(); i++)
		cout << this->FunctionSet.get_symbol(i) << " ";
	cout << endl;
}

void GenProg::show(const size_t &which)
{
	vector<vector<int> > dord_vec = this->FitControl.getDifOrds();
	DDX diff_;
	U_X u_;
	AAST diffterm, expression = this->Population.at(which).at(0);
	for(size_t i = 1; i < this->nAlleles; i++)
	{
		diffterm = this->Population.at(which).at(i);
		for(size_t j = 0; j < this->nVar; j++)
		{
			diff_.wrt_id = j;
			diff_.order = dord_vec.at(i-1).at(j);
			diffterm *= (diff_.order == 0) ? 1. : diff_.T(u_.T());
		}
		diffterm *= (i == 1) ? u_.T() : 1.;
		expression += diffterm;
	}
	expression.simplify();
	expression.show(true);
}

void GenProg::getTGEC(const size_t &which)
{
	AAST s = this->Population.at(which).at(0);
	s.simplify();
	vector<AAST> K(this->nAlleles - 1);
	for(size_t j = 0; j < this->nAlleles - 1; j++)
	{
		K.at(j) = this->Population.at(which).at(j+1);
		K.at(j).simplify();
	}
	FitControl.solve_for(K, s, true);
}
// =====================================================================
#endif
