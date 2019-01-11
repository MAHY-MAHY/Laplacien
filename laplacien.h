#ifndef egein_hpp
#define egein_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Laplacien
{
public:
	Laplacien();
	Laplacien(const double a, const double b,const int N);

	double f(double x) const;
	void Solexact();

	void MatLaplacien();
	void TermeSource();

	void SolveurDirect_QR();
	void SolveurDirect_LU();
	void SolveurDirect_LLT();
	void SolveurDirect_LDLT();

	void SolveurIteratif_ConjugateGradient();
	void SolveurIteratif_LeastSquaresConjugateGradient();
	void SolveurIteratif_BiCGSTAB();

	double CalcErreur();

	void Save(const std::string &fichier);


private:
	const double m_a;
	const double m_b;
	const int m_N;
	const double m_h;

	Eigen::SparseMatrix<double> MAT_A;

	Eigen::VectorXd Vec_X;
	Eigen::VectorXd Vec_F;

	Eigen::VectorXd Vec_Uapp;
	Eigen::VectorXd Vec_Uex;
	double m_er;
};
#endif
