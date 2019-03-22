#include "ai.hh"
#include "conjgrad.hh"
#include <vector>
// % see Algorithm 6.17 in Y. Saad, Iterative Methods for Sparse Linear Systems
// %
// function [x n_iter] = conj_grad(A, b, Nx, Nx_Ny, N_DOF)
//    residual = zeros(30, 1);
//    x  = zeros(N_DOF,1);
//    r1 = zeros(N_DOF,1);
//    r2 = zeros(N_DOF,1);
//    p  = zeros(N_DOF,1);
//
// %    r1 = b - multiply_diag(A, x, Nx, Nx_Ny, N_DOF); % for non-zero initial gues
//    r1 = b;% initial guess is x=0 !!!
//    p = r1;
//    eps = 1;
//    n_iter = 1;
//
//    while eps > 0.01
//         A_p = multiply_diag(A, p, Nx, Nx_Ny, N_DOF);
//         alpha = r1'*r1/(A_p'*p);
//         x = x + alpha*p;
//         r2 = r1 - alpha * A_p;
//         beta = (r2'*r2)/(r1'*r1);
//         p = r2 + beta * p;
//         r1 = r2;
//         eps = norm(r2)/norm(x);
//         residual(n_iter) = eps;
//         n_iter = n_iter + 1;
//    end
//    figure(11); plot(residual); title(['number of iterations = ', num2str(n_iter)]);
// end
void conjGrad(std::vector<double> &x,
                    std::vector<std::vector<double> >&A,
                            std::vector<double>&b,
                            size_t Nx , size_t NxNy , size_t N_dof){


    std::vector<double> residual;
    //std::vector<double> x(N_dof, 0.);
    std::vector<double> r1(N_dof, 0.);
    std::vector<double> r2(N_dof, 0.);
    std::vector<double> p(N_dof, 0.);


    r1 = b;
    p = r1;
    double eps = 1.;
    size_t n_iter = 1;
    
    std::vector<double> x1(N_dof , 0.);
    std::vector<double> p1(N_dof , 0.);

    std::vector<double> A_p(N_dof , 0.);
    std::vector<double> r1r1(N_dof , 0.);
    std::vector<double> r2r2(N_dof , 0.);
    std::vector<double> A_p_p(N_dof , 0.);
    double alpha;
    double beta;
    std::vector<double> alphap(N_dof , 0.);
    std::vector<double> alphaA_p(N_dof , 0.);

    std::vector<double> betap(N_dof , 0.);
    // while(0.01 < eps){
        multiplyDiag(A_p, A, p, Nx , NxNy , N_dof);
        //ai::printMarker();
		ai::saveVector("p",p);
		ai::saveVector("A_p", A_p);
        alpha = MultiplyVV(r1,r1)/MultiplyVV(A_p,p);

        // std::cout<<"alpha = "<<std::fixed<<std::setprecision(15)<<alpha<<std::endl;
        // std::cout<<"MultiplyVV(r1,r1) = "<<MultiplyVV(r1,r1)<<std::endl;
        // std::cout<<"MultiplyVV(A_p,p) = "<<MultiplyVV(A_p,p)<<std::endl;
        // MultiplyVN(p,alpha);
        // // std::cout<<"Norm(r1) = "<< NormV(r1)<<std::endl;
        // SumVV(x, p, x1);
        //x=x1;

        //x = x + alpha*p;

        for(size_t i = 0 ; i< x.size() ;++i){
            x[i]+= alpha*p[i];
        }

        // ai::saveVector("x",x);
        // MultiplyVN(A_p,-alpha);
        //
        // SumVV(r1 ,A_p, r2);

        //r2 = r1 - alpha * A_p;

        for(size_t i = 0 ; i < r2.size(); i++){
            r2[i] = r1[i] - alpha * A_p[i];
        }
        // ai::saveVector("r2",r2);
        beta = MultiplyVV(r2,r2)/MultiplyVV(r1,r1);

         //std::cout<<"beta = "<<std::fixed<<std::setprecision(15)<<beta<<std::endl;
        // MultiplyVN(p, beta);
        //
        // SumVV(r2,p,p1);
        //p=p1;

        //p = r2 + beta * p

        for(size_t i = 0 ; i < p.size(); i++){
            p[i] = r2[i] + beta*p[i];
        }

         //ai::saveVector("p",p);
        r1 = r2;

        eps = NormV(r2)/NormV(x);
        // std::cout<<"eps = "<<eps<<std::endl;
        residual.push_back(eps);
        ++n_iter;

     //}

    std::cout<<"iter = "<<n_iter<<std::endl;

}
// function y = multiply_diag(A, x, Nx, Nx_Ny, N_DOF)
//     y = zeros(N_DOF,1);
//
//     % multiply by 3rd, 4th, 5th diagonals
//     y(1) = A(1,4) * x(1) + A(1,5) *x(2);
//     y(N_DOF) = A(N_DOF, 4) * x(N_DOF) + A(N_DOF, 3) * x(N_DOF-1);
//     for ii = 2:N_DOF-1
//         y(ii) = A(ii,3) * x(ii-1) + A(ii,4) * x(ii) + A(ii,5) * x(ii+1);
//     end
//
//     % multiply by the 1st diagonal
//     for ii = Nx_Ny+1:N_DOF
//         y(ii) = y(ii) + A(ii, 1) * x(ii-Nx_Ny);
//     end
//
//     % multiply by the 2nd diagonal
//     for ii = Nx+1:N_DOF
//         y(ii) = y(ii) + A(ii, 2) * x(ii-Nx);
//     end
//
//     % multiply by the 6th diagonal
//     for ii = 1:N_DOF-Nx
//         y(ii) = y(ii) + A(ii, 6) * x(ii+Nx);
//     end
//
//     % multiply by the 7th diagonal
//     for ii = 1:N_DOF-Nx_Ny
//         y(ii) = y(ii) + A(ii, 7) * x(ii+Nx_Ny);
//     end
// end
//
void multiplyDiag(std::vector<double> &y, std::vector<std::vector<double> >&A,std::vector<double >&x,
                        size_t Nx , size_t NxNy , size_t N_dof){
    //std::vector<double> y(N_dof, 0.);
    //% multiply by 3rd, 4th, 5th diagonals
    y[0] = A[0][3]*x[0] + A[0][4]*x[1];
    y[N_dof-1] = A[N_dof-1][3]*x[N_dof-1] + A[N_dof-1][2]*x[N_dof-2];
    for(size_t i = 1; i < N_dof - 2;++i){
        y[i] = A[i][2] * x[i-1] + A[i][3] * x[i] + A[i][4] * x[i+1];
    }

    //% multiply by the 1st diagonal
        for(size_t i = NxNy ; i < N_dof ; ++i){
            y[i] += A[i][0] * x[i-NxNy];
    }
    //% multiply by the 2nd diagonal
    for(size_t i = Nx; i < N_dof ; ++i){
        y[i] += A[i][1] * x[i-Nx];
    }
    //     % multiply by the 6th diagonal
    for(size_t i = 0; i < N_dof-Nx;++i){
        y[i] += A[i][5]*x[i+Nx];
    }

    //% multiply by the 7th diagonal
   for(size_t i = 0 ; i<N_dof - NxNy;++i ){
       y[i] += A[i][6]*x[i+NxNy];
   }
}

double MultiplyVV(std::vector<double>&a, std::vector<double>&b){
    if( a.size()!=b.size() ){
        std::cout<<"Warning!!"<<std::endl;
    }
    double c = 0.;
    for(size_t i = 0; i < a.size(); i++){
        c=c+ a[i]*b[i];
    }
    return c;
}

void MultiplyVN(std::vector<double>&a, double b){


    for(size_t i = 0; i < a.size();i++ ){
        a[i]*=b;
    }
}

void DevideVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c){
    if(a.size()!=b.size() && c.size()!=b.size()){
        std::cout<<"Warning!!"<<std::endl;
    }
    for(size_t i =0 ; i < a.size();i++ )
        c[i]=a[i]/b[i];
}

void SumVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c){
    if(a.size()!=b.size() && c.size()!=b.size()){
        std::cout<<"Warning!!"<<std::endl;
    }
    std::fill(c.begin(),c.end(), 0.);
    for(size_t i = 0 ; i < a.size();i++ )
        c[i]+=a[i]+b[i];
}

void SubstractVV(std::vector<double>&a, std::vector<double>&b, std::vector<double>&c){
    if(a.size()!=b.size() && c.size()!=b.size()){
        std::cout<<"Warning!!"<<std::endl;
    }
    for(size_t i =0 ; i < a.size();i++ )
        c[i]=a[i]-b[i];
}

double NormV(std::vector<double>&x){
    double a = 0.;
    for(size_t i = 0 ; i< x.size(); ++i){
        a = a+x[i]*x[i];
    }
    return sqrt(a);
}
