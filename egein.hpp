#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifndef egein_hpp
#define egein_hpp

class Laplacien{
private:
  const double a;
  const double b;
  const int N;
  const double h;
  Eigen::SparseMatrix<double> MAT_A;
  Eigen::VectorXd Vec_X;
  Eigen::VectorXd Vec_F;
  Eigen::VectorXd Vec_Uapp;
  Eigen::VectorXd Vec_Uex;
  double er;

public:
  Laplacien();
  Laplacien(const double a, const double b,const int N);
  void MatLaplacien();
  void TermeSource();
  void SolveurDirect();
  void SolveurIteratif();
  double CalcErreur();
  void Save(std::string fichier);
  
  
  
};
#endif

  
