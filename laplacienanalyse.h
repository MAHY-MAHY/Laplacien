#ifndef LAPLACIENANALYSE_H
#define LAPLACIENANALYSE_H

#include <vector>

#include "laplacien.h"

class LaplacienAnalyse
{
public:
	LaplacienAnalyse();
	void convergence(const std::string &fichier) const;
	void graph() const;
	void temps(const std::string &fichier) const;


private:

	const int N_erreur[6] = {10,50,90,200,600,1000};
};

#endif // LAPLACIENANALYSE_H
