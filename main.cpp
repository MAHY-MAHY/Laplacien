#include <iostream>

#include "laplacien.h"
#include "laplacienanalyse.h"

int main(int , char *[] )
{
	std::cout << "Hello World!\n";

	LaplacienAnalyse LA;

	LA.convergence("erreur");
	LA.temps("erreur");

	return 0;
}
