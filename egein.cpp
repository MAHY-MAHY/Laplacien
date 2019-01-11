#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "egein.hpp"

Laplacien(){
  m_a=0;
  m_b=1;
  m_N=2;
  m_h=(m_b-m_a)/(m_N+1);
}

Laplacien(const double a, const double b, const int N ){
  m_a=a;
  m_b=b;
  m_N=N;
  m_h=(b-a)/(N+1);
}

void Matlaplacien(){}

void TermeSource(){}

void SolveurDirect(){}

double CalcErreur(){}

void Save(std::string fichier){}
