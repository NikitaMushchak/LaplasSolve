#pragma once
#include <vector>

void multiplyDiag(std::vector<double> &y, std::vector<double >&x,
                        std::vector<std::vector<double> >&A,size_t Nx , size_t NxNy , size_t N_dof);

void MultiplyVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c);

void DevideVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c);

void SumVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c);

void SubstractVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c);

double NormV(std::vector<double>&x);



void conjGrad(std::vector<std::vector<double> >&A,
                            std::vector<double>&b,
                            size_t Nx , size_t NxNy , size_t N_dof);
