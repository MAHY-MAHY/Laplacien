#include <fstream>
#include <iostream>
#include<Eigen/SparseQR>
#include<Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>

#include "laplacien.h"

Laplacien::Laplacien(): Laplacien(0., 1., 1)
{
}

Laplacien::Laplacien(const double a, const double b, const int N ): m_a(a), m_b(b), m_N(N), m_h( (m_b-m_a)/(m_N+1)), m_er(0.)
{
	Vec_X.resize(N+2);
	for(int i=0; i<=N+1; ++i)
	{
		Vec_X(i) = m_a + i*m_h;
	}

	// Calcul du terme source
	TermeSource();
	// Calcul de la matrice
	MatLaplacien();
	// Calcule de la solution exacte
	Solexact();
}

double Laplacien::f(double x) const
{
	return -std::exp(-x) * (2.*(1. - 2.*x + m_a + m_b) +(x-m_a)*(x-m_b));
}

void Laplacien:: Solexact()
{
	Vec_Uex.resize(m_N);
	for(int i=1; i<=m_N; i++)
	{
		double x = Vec_X(i);
		Vec_Uex(i-1) = (x - m_a)*(x - m_b)*std::exp( -x);
	}
}

void Laplacien::MatLaplacien()
{
	int i;
	MAT_A.resize(m_N,m_N);
	for(i=0;i<m_N;i++){
		MAT_A.coeffRef(i, i)=2./m_h/m_h;
	}
	for(i=0;i<m_N-1;i++){
		MAT_A.coeffRef(i, i+1)=-1./m_h/m_h;
		MAT_A.coeffRef(i+1,i)=-1./m_h/m_h;
	}
}

void Laplacien::TermeSource()
{
	Vec_F.resize(m_N);
	for(int i=1; i<=m_N; ++i)
	{
		Vec_F(i-1) = f(Vec_X(i));
	}
}

void Laplacien::SolveurDirect_QR()
{
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	MAT_A.makeCompressed();
	solver.analyzePattern(MAT_A);
	// Compute the numerical factorization
	solver.factorize(MAT_A);

	//Use the factors to solve the linear system
	Vec_Uapp = solver.solve(Vec_F);

}

void Laplacien::SolveurDirect_LU()
{
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(MAT_A);
	// Compute the numerical factorization
	solver.factorize(MAT_A);

	//Use the factors to solve the linear system
	Vec_Uapp = solver.solve(Vec_F);
}

void Laplacien::SolveurDirect_LLT()
{
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(MAT_A);
	// Compute the numerical factorization
	solver.factorize(MAT_A);

	//Use the factors to solve the linear system
	Vec_Uapp = solver.solve(Vec_F);
}

void Laplacien::SolveurDirect_LDLT()
{
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(MAT_A);
	// Compute the numerical factorization
	solver.factorize(MAT_A);

	//Use the factors to solve the linear system
	Vec_Uapp = solver.solve(Vec_F);
}

void Laplacien::SolveurIteratif_ConjugateGradient()
{
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > C_Grad;
	C_Grad.compute(MAT_A);
	Vec_Uapp=C_Grad.solve(Vec_F);
}

void Laplacien::SolveurIteratif_LeastSquaresConjugateGradient()
{
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > LsC_Grad;
	LsC_Grad.compute(MAT_A);
	Vec_Uapp=LsC_Grad.solve(Vec_F);

}

void Laplacien::SolveurIteratif_BiCGSTAB()
{
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > BCG;
	BCG.compute(MAT_A);
	Vec_Uapp=BCG.solve(Vec_F);

}

double Laplacien::CalcErreur()
{
	int i;
	m_er=0;
	for(i=0;i<m_N;i++){
		if(std::abs(Vec_Uex(i)-Vec_Uapp(i))>m_er){
			m_er=std::abs(Vec_Uex(i)-Vec_Uapp(i));
		}
	}
	return m_er;
}

void Laplacien::Save(const std::string &fichier)
{
	std::ofstream myfile;
	myfile.open(fichier, std::ios::in | std::ios::out | std::ios::trunc);

	myfile << Vec_X(0) << " " << 0. << " " << 0. << "\n";

	for(int i=1;i<=m_N;++i)
	{
		myfile << Vec_X(i) << " " << Vec_Uex(i-1) << " " << Vec_Uapp(i-1) << "\n";
	}

	myfile << Vec_X(m_N+1) << " " << 0. << " " << 0. << "\n";
	myfile.close();
}
