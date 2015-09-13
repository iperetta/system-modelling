/* code by: iperetta@ieee.org */
#include <iostream>
#include <stdexcept>
#if defined(_WIN32) || defined(_WIN64)
 #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include "myfun.h"
#include "GenProg.h"
using namespace std;

int main(int argc, char **argv)
{
	if(argc != 8)
	{
		cerr << "Arguments == 0: set basic, 1: set common / namefile"
			<< " / poly_degree / max_diff_order / mcint powerof2"
			<< " / population / maxgen !!" << endl;
		return -1;
	}
	else
	{
		// GP: nPop_, nGen_, nTourn_, prXO_, prMT_
		GenProg GP(my::str2num<size_t>(argv[6]), 
			my::str2num<size_t>(argv[7]), 3, 0.95, 0.15);
		// Terminal and Function sets
		if(my::str2num<bool>(argv[1]))
		{
			GP.sets_common(); // set_basic , other functions included
			cout << "Using set_common" << endl;
		}
		else
		{
			GP.sets_basic(); // ARITHMETIC OPERATORS (+,-,*,/) 
			cout << "Using set_basic" << endl;
		} 
		// problem-derived constants (mass, gravity, etc)
		vector<double> numbers = {1.,2.};//-4.*M_PI*sqrt(1./(8.*pow(M_PI,3)))};
		GP.TerminalSet.include(numbers);
		// GP.INIT: namefile, poly_deg, diff_ord, MC power_of_2
		GP.init(argv[2], my::str2num<size_t>(argv[3]), 
			my::str2num<size_t>(argv[4]), 
			my::str2num<double>(argv[5]));
		GP.show_sets();
		GP.initializePopulation(7);
		GP.fitnessPopulation();
		GP.gatherData();
		cout << "best from random initialized population: " << endl;
		GP.show(GP.bestID);
		cout << "with fitness [ " << GP.Fitness.at(GP.bestID) << 
			" ]" << endl;
		int ctrl = 0;
		double oldSTD = 0.;
		for(size_t g = 0; g < GP.nGen; g++)
		{
			cout << "*** Generation " << g+1 << " ***" << endl;
			GP.generateOffspring();
			GP.fitnessOffspring();
			GP.reducePopulation(true);
			GP.gatherData();
			cout << "best so far: " << endl;
			GP.show(GP.bestID);
			cout << "with fitness [ " << GP.Fitness.at(GP.bestID) << 
				" ]" << endl;
			cout << "------------" << endl;
			if(my::isLikelyZero(oldSTD - GP.current_stdFit))
				ctrl++;
			else
				oldSTD = GP.current_stdFit;
			if(ctrl >= 7)
				break;
		}
		cout << "*** END ***" << endl;
		GP.getTGEC(GP.bestID);
		cout << "Best Individual: " << endl;
		GP.show(GP.bestID);
		cout << "Best Fitness: " << GP.Fitness.at(GP.bestID) << endl;
		cout << "PID: " << PID << endl;
		cout << "------------" << endl;
		my::saveNumVecAs(GP.meanFit_vec, "graph/GP_"+my::num2str(PID)+
			"meanfit_vec.txt");
		my::saveNumVecAs(GP.bestFit_vec, "graph/GP_"+my::num2str(PID)+
			"bestfit_vec.txt");
		return 0;
	}
}

