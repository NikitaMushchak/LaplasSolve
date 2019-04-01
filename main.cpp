#include <iostream>
#include <algorithm>
#include <cmath>

#include <vector>

#include "ai.hh"
#include "createMatrix.hh"
#include "conjgrad.hh"

#define M_PI           3.141592653589793238463

int main(){
    size_t Nx = 15;
    size_t Ny = 2 * Nx - 1;
    size_t Nz = Nx;


    double   dx = 1.;   //  mesh size in x, y, and z directions
    double    E = 1.;   //Young modulus
    double    R = 10.;
    double   nu = 0.25; // Poisson


    size_t N_dof = Nx * Ny * Nz;    // number of unknowns in Laplace equation
    std::vector<double> T(N_dof, 0.);  // unknowns in  finite difference discretization of Laplace equation (AT = b)

    //std::vector<double> A0(N_dof, 0.);
    std::vector<double> A1(N_dof - Nx, 0.);
    std::vector<double> A2(N_dof, 0.);
    std::vector<double> A3(N_dof, 0.);
    std::vector<double> A4(N_dof, 0.);
    std::vector<double> A5(N_dof - Nx, 0.);
    //std::vector<double> A6(N_dof, 0.);


    std::vector<double> b(N_dof, 0.);  // right side of finite difference discretization of Laplace equation
    std::vector<double> r1(N_dof, 0.);
    std::vector<double> r2(N_dof, 0.);
    std::vector<double> p(N_dof, 0.);
    std::vector<double> A_p(N_dof, 0.);


	std::cout<<"Nx = "<<Nx<<" Ny = "<<Ny<<"  Nz= "<<Nz<<std::endl;;
    std::cout<<"A3 size : "<<A3.size()<<std::endl;
//

    std::vector<std::vector<double> > Tz1; //раскрытие
    Tz1.resize(Nx);
    for(size_t i = 0 ; i < Tz1.size();++i){
        Tz1[i].resize(Ny);
    }
  double r;

    for(size_t j = 0 ; j < Ny; ++j){
        for(size_t i = 0 ; i < Nx; ++i){
                r = std::sqrt( i*i+std::pow(j-0.5*(Ny)-0.5, 2) );
                if(r < R){
                    Tz1[i][j-1] = 0.5*std::cos((0.5 * M_PI *r)/R);
                }
        }
    }
auto start = ai::time();
    // ai::saveMatrix("Tz1",Tz1);

    // create matrix corresponding to finite difference discretization of Laplace equation (AT = b)

    createMatrixDiag(   A1,
                        A2,
                        A3,
                        A4,
                        A5,
                        N_dof,
                        Nx,
                        Ny,
                        Nz);
        //ai::printMarker();
    //ai::saveVector("A0",A0);
    // ai::saveVector("A1",A1);
    // ai::saveVector("A2",A2);
    // ai::saveVector("A3",A3);
    // ai::saveVector("A1",A4);
    // ai::saveVector("A5",A5);
    // ai::saveVector("A6",A6);
//     % create right side using boundary condition T|_{z=0} = Tz1(x,y)
// for jj=1:Ny
//     for ii=1:Nx
//         b(ii + (jj-1)*Nx) = -Tz1(ii,jj);
//     end
// end
    for(size_t j = 1 ; j < Ny; ++j){
        for(size_t i = 0 ; i < Nx; ++i){
            b[i + j*Nx] = - Tz1[i][j];
        }
    }
    // ai::saveVector("b",b);

    // ai::saveMatrix("A", A);

    //     %  solve A*T=b using conjugate gradient method
    // [T n_iter] = conj_grad(A, b, Nx, Nx*Ny, N_DOF)
    auto t1 = ai::time();
    // conjGrad(T ,A, b, Nx , Nx*Ny , N_dof);
    conjGrad(T,                   // выход функции
                A1,
                A2,
                A3,
                A4,
                A5,
                b,       // вектор раскрытий
                r1,
                r2,
                p,
                A_p,
                Nx,
                Nx*Ny,
                N_dof);
    auto t2 = ai::time();
    std::cout<<"conj grad time = "<<ai::duration(t1,t2, "us")<<" us"<<std::endl;
    // ai::saveVector("Pres", T);


    std::vector<std::vector<double> > press;
    press.resize(Nx);
    for(size_t i = 0;i < press.size();++i)
        press[i].resize(Ny);

    for(size_t i = 0; i < Nx;++i){
        for(size_t j = 0; j < Ny;++j){
            press[i][j] = (0.5*E/(1.- nu*nu)) * (-b[i + j*Nx] - T[i + j*Nx])/dx;
        }
    }

    auto finish = ai::time();
    std::cout<<"TIME = "<<ai::duration(start, finish , "us")<<" us "<<std::endl;
    // ai::saveMatrix("pressure", press);

    return 1;
}
