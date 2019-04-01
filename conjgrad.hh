#pragma once
#include <vector>
#include "ai.hh"

void multiplyDiag(std::vector<double> &y,
                        //std::vector<double>& A0,
                        std::vector<double>& A1,
                        std::vector<double>& A2,
                        std::vector<double>& A3,
                        std::vector<double>& A4,
                        std::vector<double>& A5,
                        //std::vector<double>& A6,
                        std::vector<double >&x,
                        size_t Nx,
                        size_t NxNy,
                        size_t N_dof);

double MultiplyVV(std::vector<double>&a, std::vector<double>&b);

void MultiplyVN(std::vector<double>&a, double b);

void DevideVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c);

void SumVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c);

void SubstractVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c);

double NormV(std::vector<double>&x);



void conjGrad(std::vector<double> &x,                   // выход функции
                            //std::vector<double>& A0,
                            std::vector<double>& A1,
                            std::vector<double>& A2,
                            std::vector<double>& A3,
                            std::vector<double>& A4,
                            std::vector<double>& A5,
                            //std::vector<double>& A6,// семидиагональная матрица
                            std::vector<double>& b,       // вектор раскрытий
                            std::vector<double>& r1,
                            std::vector<double>& r2,
                            std::vector<double>& p,
                            std::vector<double>& A_p,
                            size_t Nx , size_t NxNy , size_t N_dof);
