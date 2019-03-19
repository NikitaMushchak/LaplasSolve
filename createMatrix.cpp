
#include <vector>

#include "ai.hh"


// function A = create_matrix_diag(N_DOF, Nx, Ny, Nz)
//
//     Nx_Ny = Nx*Ny;
//     % create all diagonals
//     for ii = 1:N_DOF
//         A(ii,1) = 1;
//         A(ii,2) = 1;
//         A(ii,3) = 1;
//         A(ii,4) = -6;
//         A(ii,5) = 1;
//         A(ii,6) = 1;
//         A(ii,7) = 1;
//     end
//
//     % set some elements of A equal to zero using boundary conditions
//     % zero flux at x=0 and x=Nx*dx
//     for jj=1:Ny
//         for kk=1:Nz
//             n = 1 + (jj-1)*Nx + (kk-1)*Nx_Ny; % i=1
//             A(n, 3) = 0;
//             A(n, 5) = 2;
//
//             p = n + Nx -1; % i=Nx
//             A(p, 5) =  0;
//             A(p, 4) = -5;
//         end
//     end
//
//     % zero flux at y=0 and y=Ny*dy
//     for ii=1:Nx
//         for kk=1:Nz
//             n = ii + (kk-1)*Nx_Ny; % j=1
//             A(n, 2) =  0;
//             A(n, 4) = -5;
//
//             p = n + (Ny-1)*Nx; % j=Ny
//             A(p, 6) =  0;
//             A(p, 4) = -5;
//         end
//     end
//
//     % prescribed temperature at z=0 and zero flux at z=Nz*dz
//     for ii=1:Nx
//         for jj=1:Ny
//             n = ii + (jj-1)*Nx;
//             A(ii + (jj-1)*Nx, 1) = 0; % k=1
//
//             A(n + (Nz-1)*Nx_Ny, 6) =  0; % k=Nz
//             A(n + (Nz-1)*Nx_Ny, 4) = -5;
//         end
//     end
// end



void createMatrixDiag(std::vector<std::vector<double> >& A,
                            size_t N_dof, size_t Nx, size_t Ny, size_t Nz){

    size_t NxNy = Nx * Ny;
//     % create all diagonals
    for(size_t i = 0; i < N_dof ;++i){
        A[i][0] = 1.;
        A[i][1] = 1.;
        A[i][2] = 1.;
        A[i][3] = -6.;
        A[i][4] = 1.;
        A[i][5] = 1.;
        A[i][6] = 1.;
    }
    // ai::printMarker();
    //% set some elements of A equal to zero using boundary conditions
   //     % zero flux at x=0 and x=Nx*dx

   double n;
   double p ;
   for(size_t j = 0 ; j < Ny; ++j){
       for(size_t k = 0 ; k < Nz; ++k){
           n = j * Nx + k * NxNy;
           A[n][2] = 0.;
           A[n][4] = 2.;

           p = n + Nx - 1;
           A[p][4] = 0.;
           A[p][3] = -5.;
       }
   }
// ai::printMarker();
  // % zero flux at y=0 and y=Ny*dy
  //     for ii=1:Nx
  //         for kk=1:Nz
  //             n = ii + (kk-1)*Nx_Ny; % j=1
  //             A(n, 2) =  0;
  //             A(n, 4) = -5;
  //
  //             p = n + (Ny-1)*Nx; % j=Ny
  //             A(p, 6) =  0;
  //             A(p, 4) = -5;
  //         end
  //     end

   for(size_t i = 0; i < Nx; ++i ){
      for(size_t k = 0; k < Nz; ++k){
          n = i  + k * NxNy;
          A[n][1] = 0.;
          A[n][3] = -5.;

          p = n + (Ny-1)*Nx;
          A[p][5] = 0.;
          A[p][3] = -5.;
      }
  }
// ai::printMarker();
  for(size_t i = 0; i<Nx ;++i){
      for(size_t j = 0; j < Ny; ++j){
          n = i  + j * Nx;
          A[i + j*Nx][0] = 0.; //% k=1

          A[n + (Nz-1)*NxNy][5] = 0.;
          A[n + (Nz-1)*NxNy][3] = -5.;
      }
  }
// ai::printMarker();
}
