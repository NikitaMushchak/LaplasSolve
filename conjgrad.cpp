#include "ai.hh"
#include "conjgrad.hh"
#include <vector>


void conjGrad(std::vector<double> &x,                   // выход функции
                        //    std::vector<double>& A0,
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
                            size_t Nx , size_t NxNy , size_t N_dof){

for(size_t i = 0 ; i < N_dof ; ++i){
        r2[i] = 0.;
        p[i]  = 0.;
        A_p[i]= 0.;
        r1[i] = b[i];
        p[i] = r1[i];
    }
    double eps = 1.;
    double alpha;
    double beta;

    while(0.0001 < eps){
        multiplyDiag(A_p, A1, A2, A3, A4, A5, p, Nx , NxNy , N_dof);
        //ai::printMarker();
		// ai::saveVector("p",p);
		// ai::saveVector("A_p", A_p);
        alpha = MultiplyVV(r1,r1)/MultiplyVV(A_p,p);
        for(size_t i = 0 ; i < x.size() ;++i){
            x[i]+= alpha*p[i];
            r2[i] = r1[i] - alpha * A_p[i];
        }
        // ai::saveVector("x",x);
        // ai::saveVector("r2",r2);
        beta = MultiplyVV(r2,r2)/MultiplyVV(r1,r1);
         //std::cout<<"beta = "<<std::fixed<<std::setprecision(15)<<beta<<std::endl;
        for(size_t i = 0 ; i < p.size(); i++){
            p[i] = r2[i] + beta*p[i];
            r1[i] = r2[i];
        }
         //ai::saveVector("p",p);
        // r1 = r2;
        eps = NormV(r2)/NormV(x);
        // std::cout<<"eps = "<<eps<<std::endl;
     }
}
//
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
                        size_t N_dof){
    //std::vector<double> y(N_dof, 0.);
    //% multiply by 3rd, 4th, 5th diagonals
    y[0] = A3[0]*x[0] + A4[0]*x[1];
    y[N_dof-1] = A3[N_dof-1]*x[N_dof-1] + A2[N_dof-1]*x[N_dof-2];

    for(size_t i = 1; i < N_dof - 2;++i){
        y[i] = A2[i] * x[i-1] + A3[i] * x[i] + A4[i] * x[i+1];
    }
    //% multiply by the 1st diagonal
        for(size_t i = NxNy ; i < N_dof ; ++i){
            y[i] += x[i-NxNy];
    }
    //% multiply by the 2nd diagonal
    // for(size_t i = Nx; i < N_dof ; ++i){
    //     y[i] += A1[i] * x[i-Nx];
    // }
    for(size_t i = 0; i < N_dof - Nx ; ++i){
        y[i+Nx] += A1[i] * x[i];
    }
    //     % multiply by the 6th diagonal
    for(size_t i = 0; i < N_dof-Nx;++i){
        y[i] += A5[i]*x[i+Nx];
    }

    //% multiply by the 7th diagonal
   for(size_t i = 0 ; i < N_dof - NxNy;++i ){
       y[i] += x[i + NxNy];
   }
}

double MultiplyVV(std::vector<double>&a, std::vector<double>&b){

    double c = 0.;
    for(size_t i = 0; i < a.size(); i++){
        c+= a[i]*b[i];
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
    return a;
}
