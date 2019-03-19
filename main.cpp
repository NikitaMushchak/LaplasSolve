#include <iostream>
#include <algorithm>
#include <cmath>

#include <vector>

#include "ai.hh"
#include "createMatrix.hh"
#include "conjgrad.hh"

# define M_PI           3.14159265358979323846

int main(){
    size_t Nx = 15;
    size_t Ny = 2 * Nx - 1;
    size_t Nz = Nx;


    double    R = 10.;
    double   dx = 1.;   //  mesh size in x, y, and z directions
    double    E = 1.;   //Young modulus
    double   nu = 0.25; // Poisson


    double N_dof = Nx * Ny * Nz;    // number of unknowns in Laplace equation
    std::vector<double> T(N_dof, 0.);  // unknowns in  finite difference discretization of Laplace equation (AT = b)
    std::vector<std::vector<double> > A; // matrix corresponding to finite difference discretization of Laplace equation (AT = b)
    A.resize(N_dof);                     //. Nonzero elements are stored only.
    for(size_t i = 0 ; i<A.size(); ++i){
        A[i].resize(7);
    }
    std::vector<double> b(N_dof, 0.);  // right side of finite difference discretization of Laplace equation

    std::cout<<"A size : "<<A.size()<<" x "<<A[0].size()<<std::endl;
//     % set boundary conditions at z = 0 and z=(Nz-1)*dx
// Tz1 = zeros(Nx, Ny);  % at z = 0
// for jj=1:Ny
//     for ii=1:Nx
//         r = sqrt((ii-1)^2+(jj-0.5*Ny-0.5)^2);
//         if r < R
//             Tz1(ii,jj) = 0.5*cos(0.5*pi*r/R);%cos(0.5*pi*r/R); % Tz1(x,y) is equal to a HALF of crack opening 0.5*w(x,y)
//         end
//     end
// end

    std::vector<std::vector<double> > Tz1;
    Tz1.resize(Nx);
    for(size_t i = 0 ; i < Tz1.size();++i){
        Tz1[i].resize(Ny);
    }
  double r;
    for(size_t j = 0 ; j < Ny; ++j){
        for(size_t i = 0 ; i < Nx; ++i){
                r = std::sqrt(i*i+std::pow(j-0.5*Ny-0.5,2));
                if(r < R){
                    Tz1[i][j] = 0.5*cos(0.5 * M_PI *r/R);
                }
        }
    }

    ai::saveMatrix("matr",Tz1);

    // create matrix corresponding to finite difference discretization of Laplace equation (AT = b)

    createMatrixDiag( A, N_dof,  Nx, Ny, Nz);


//     % create right side using boundary condition T|_{z=0} = Tz1(x,y)
// for jj=1:Ny
//     for ii=1:Nx
//         b(ii + (jj-1)*Nx) = -Tz1(ii,jj);
//     end
// end
    for(size_t j = 1 ; j < Ny; ++j){
        for(size_t i = 0 ; i < Nx; ++i){
            b[i + (j-1)*Nx] = - Tz1[i][j];
        }
    }
    ai::saveVector("b",b);

    ai::saveMatrix("A", A);



    //     %  solve A*T=b using conjugate gradient method
    // [T n_iter] = conj_grad(A, b, Nx, Nx*Ny, N_DOF)
    conjGrad(A, b, Nx , Nx*Ny , N_dof);

    ai::saveMatrix("Pres", A);
    return 1;
}
