#include "NEGA2.hpp"

#define DTLZ1
#include "TESTFUNCTION.hpp"

extern double limits_array[10][2];


int main()
{
	NSGA2 nsga(10, 3, 500, 0.01, 0.9, 15.0, 20.0,test_problem);
	nsga.set_limits(limits_array, 10);
	nsga.initial_population();


	
	nsga.run(100);
	nsga.display();

	return 0;
}