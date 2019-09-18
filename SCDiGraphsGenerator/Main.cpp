#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include "StronglyConnectedGraphGenerator.h"

using namespace SCDiGraphs;

int main()
{
	while (true)
	{
		auto res = Generator::generateGraphAndInducedSubgraphs(100, 100);
		//Generator::print(res);
		if (!Generator::isStronglyConnected(res)) { /*Generator::print(res);*/ break; }
	}
	return 0;
}
