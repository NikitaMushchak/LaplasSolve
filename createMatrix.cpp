
#include <vector>

#include "ai.hh"

void createMatrixDiag(
                        std::vector<double>& A1,
                        std::vector<double>& A2,
                        std::vector<double>& A3,
                        std::vector<double>& A4,
                        std::vector<double>& A5,
                        //std::vector<double>& A6,
                            size_t N_dof, size_t Nx, size_t Ny, size_t Nz){

    size_t NxNy = Nx * Ny;
//     % create all diagonals
    for(size_t i = 0; i < N_dof ;++i){
        // A0[i] = 1.;
        //A1[i] = 1.;
        A2[i] = 1.;
        A3[i] = -6.;
        A4[i] = 1.;
        //A5[i] = 1.;
        //A6[i] = 1.;
    }
    for(size_t i = 0; i<N_dof - Nx; ++i){
        A1[i] = 1.;
        A5[i] = 1.;
    }
    // ai::printMarker();
    //% set some elements of A equal to zero using boundary conditions
   //     % zero flux at x=0 and x=Nx*dx

   double n;
   double p ;
   for(size_t j = 0 ; j < Ny; ++j){
       for(size_t k = 0 ; k < Nz; ++k){
           n = j * Nx + k * NxNy;
           A2[n] = 0.;
           A4[n] = 2.;

           p = n + Nx - 1;
           A4[p] = 0.;
           A3[p] = -5.;
       }
   }
// ai::printMarker();
   for(size_t i = 0; i < Nx; ++i ){
      for(size_t k = 0; k < Nz; ++k){
          n = i  + k * NxNy;
          if(n-Nx < N_dof-Nx){
            A1[n-Nx] = 0.;
          }
          A3[n] = -5.;

          p = n + (Ny-1)*Nx;
          if(p <= N_dof-Nx-1){
            A5[p] = 0.;
          }
          A3[p] = -5.;
      }
  }
// ai::printMarker();
  for(size_t i = 0; i<Nx ;++i){
      for(size_t j = 0; j < Ny; ++j){
          n = i  + j * Nx;

          if(n + (Nz-1)*NxNy < N_dof-Nx){
              A5[n + (Nz-1)*NxNy] = 0.;
          }
          A3[n + (Nz-1)*NxNy] = -5.;
      }
  }
}
