#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "egein.hpp"

Laplacien::Laplacien(): m_a(0.), m_b(1.), m_N(2), m_h((m_b-m_a)/(m_N+1))
{
}

Laplacien::Laplacien(const double a, const double b, const int N ): m_a(a), m_b(b), m_N(N), m_h( (m_b-m_a)/(m_N+1))
{
}

double Laplacien::f(double x) const
{
	return -std::exp(-x) * (2. + 2.*x - m_a - m_b - (x-m_a)*(x-m_b));
}

void Laplacien::Matlaplacien(){}

void Laplacien::TermeSource(){}

void Laplacien::SolveurDirect(){}

double Laplacien::CalcErreur(){}

void Laplacien::Save(std::string fichier){}
