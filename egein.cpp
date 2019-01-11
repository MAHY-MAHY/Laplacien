#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "egein.hpp"

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
}

double Laplacien::f(double x) const
{
	return -std::exp(-x) * (2. + 2.*x - m_a - m_b - (x-m_a)*(x-m_b));
}

void Laplacien::Matlaplacien(){}

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
	
}

void Laplacien::SolveurDirect_LU()
{
	
}

void Laplacien::SolveurDirect_Cholesky()
{
	
}

void Laplacien::SolveurIteratif_ConjugateGradient()
{
	
}

void Laplacien::SolveurIteratif_LeastSquaresConjugateGradient()
{
	
}

void Laplacien::SolveurIteratif_BiCGSTAB()
{
	
}

double Laplacien::CalcErreur(){}

void Laplacien::Save(const std::string &fichier){}
