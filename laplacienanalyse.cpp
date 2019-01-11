#include <chrono>
#include <fstream>
#include <string>

#include "laplacienanalyse.h"

LaplacienAnalyse::LaplacienAnalyse()
{

}

//Analyse de la convergence des méthodes de résolution
void LaplacienAnalyse::convergence(const std::string &fichier) const
{
	std::ofstream myfile;
	myfile.open(fichier, std::ios::in | std::ios::out | std::ios::trunc);

	for(int i=0; i<6; ++i)
	{
		std::vector<double> tmp;
		tmp.resize(7);
		Laplacien *L = new Laplacien(0.,1, N_erreur[i]);

		L->SolveurDirect_LU();
		tmp.at(0) = L->CalcErreur();

		L->SolveurDirect_QR();
		tmp.at(1) = L->CalcErreur();

		L->SolveurDirect_LLT();
		tmp.at(2) = L->CalcErreur();

		L->SolveurDirect_LDLT();
		tmp.at(3) = L->CalcErreur();

		L->SolveurIteratif_BiCGSTAB();
		tmp.at(4) = L->CalcErreur();

		L->SolveurIteratif_ConjugateGradient();
		tmp.at(5) = L->CalcErreur();

		L->SolveurIteratif_LeastSquaresConjugateGradient();
		tmp.at(6) = L->CalcErreur();

		//Sauvegarde dans un fichier
		myfile << std::log(N_erreur[i]);
		for(unsigned int j=0; j<7; ++j)
		{
			myfile << " " << std::log(tmp.at(j));
		}
		myfile << "\n";

		delete L;
	}
}

//Crée differents graph pour les deux méthodes choisies
void LaplacienAnalyse::graph() const
{
	Laplacien *L = new Laplacien(0.,1, 10);
	L->SolveurDirect_LU();
	L->Save("Directe_10");
	L->SolveurIteratif_BiCGSTAB();
	L->Save("Indirecte_10");
	delete L;

	L = new Laplacien(0.,1, 100);
	L->SolveurDirect_LU();
	L->Save("Directe_100");
	L->SolveurIteratif_BiCGSTAB();
	L->Save("Indirecte_100");
	delete L;

	L = new Laplacien(0.,1, 500);
	L->SolveurDirect_LU();
	L->Save("Directe_500");
	L->SolveurIteratif_BiCGSTAB();
	L->Save("Indirecte_500");
	delete  L;
}

//Analyse du temps utilisés par les méthodes de résolution
void LaplacienAnalyse::temps(const std::string &fichier) const
{
	std::ofstream myfile;
	myfile.open(fichier, std::ios::in | std::ios::out | std::ios::trunc);

	for(int i=0; i<6; ++i)
	{
		std::vector<double> tmp;
		tmp.resize(7);
		Laplacien *L = new Laplacien(0.,1, N_erreur[i]);

		myfile << N_erreur[i];

		auto start = std::chrono::high_resolution_clock::now();//Temps de réference
		L->SolveurDirect_LU();
		auto diff = std::chrono::high_resolution_clock::now() - start;//Temps passé
		myfile << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();

		start = std::chrono::high_resolution_clock::now();
		L->SolveurDirect_QR();
		diff = std::chrono::high_resolution_clock::now() - start;
		myfile << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();

		start = std::chrono::high_resolution_clock::now();
		L->SolveurDirect_LLT();
		diff = std::chrono::high_resolution_clock::now() - start;
		myfile << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();

		start = std::chrono::high_resolution_clock::now();
		L->SolveurDirect_LDLT();
		diff = std::chrono::high_resolution_clock::now() - start;
		myfile << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();

		start = std::chrono::high_resolution_clock::now();
		L->SolveurIteratif_BiCGSTAB();
		diff = std::chrono::high_resolution_clock::now() - start;
		myfile << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();

		start = std::chrono::high_resolution_clock::now();
		L->SolveurIteratif_ConjugateGradient();
		diff = std::chrono::high_resolution_clock::now() - start;
		myfile << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();

		start = std::chrono::high_resolution_clock::now();
		L->SolveurIteratif_LeastSquaresConjugateGradient();
		diff = std::chrono::high_resolution_clock::now() - start;
		myfile << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();

		myfile << "\n";

		delete L;
	}
}
