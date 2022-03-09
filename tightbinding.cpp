#include <iostream>
#include<bits/stdc++.h>
#include <cmath>
#include<complex>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include<vector>
using namespace std;
using namespace Eigen;



int main()
{   
    int N=100;
    //double ep[2*N] ={};
    //double en[2*N] ={};
    double ep1[N]={};
    double en1[N]={};
    double ep2[N]={};
    double en2[N]={};
    double a = 2.46;
    double t0 = 3.033*2.46*2.46/(a*a); 
    double s0 = 0;
    double k_x1 = 0;
    double k_y1 = 0;
    double k_x = (2*M_PI/(a*sqrt(3)));
    double k_y = (2*M_PI/(3*a));
    double rad2deg =   1 ; // 180/M_PI;
    const complex<double> i(0.0,1.0);
   /* cout<<"k_x"<<k_x<<endl;                                 these are ok
    cout<<"k_y"<<k_y<<endl;   */ 
    for(int loop=0;loop<N;loop++)
    {

         //for k-to-M path  ****************
        double dx =2*M_PI/(a*sqrt(3));
        double dy =2*M_PI/(3*a)/N;
        //for k-to-M path
      //  complex<double> m= cos(k_x*a/sqrt(3)) + i*sin(k_x*a/sqrt(3)) ;
        complex<double> f = exp(rad2deg*i*k_x*a/sqrt(3)) + 2*(cos(rad2deg*k_y*a/2))*exp(-i*rad2deg*k_x*a/(2*sqrt(3)));
        
                   //std::conj(mycomplex)
      //  complex<double> f_c =
         complex<double> f_c =conj(f);
        //std::cout<< f << endl;    //
        /*std::cout<< f_c << endl; */
       // std:: cout<<'f'<< sqrt(f_c*f) <<std::endl     ;
        complex<double> sd= 1;//1/(1- (s0*s0*(abs(f)*abs(f))));
        complex<double> ss= 0;//s0*t0*(abs(f))*(abs(f));
        
        //std::cout<<"SD : "<< sd <<"SS : "<< ss <<std::endl;
        //std::cout<<"SS : "<< ss <<std::endl; 

        Eigen::MatrixXcd tbh_K_M(2,2);
        tbh_K_M(0,0)= sd*ss;
        tbh_K_M(1,0)= sd*(-t0*f);
        tbh_K_M(0,1)= sd*(-t0*f_c);
        tbh_K_M(1,1)= sd*ss;

        //std:: cout<<"f&f_c"<< tbh_K_M <<endl;

       // cout<< tbh_K_M(0,0) << tbh_K_M(0,1)<< tbh_K_M(1,0)<< tbh_K_M(1,1);
        //std:: cout<< tbh_K_M <<endl;
        Eigen::ComplexEigenSolver<MatrixXcd> ces;
        ces.compute(tbh_K_M);    // computes complex eigenvalues and eigenvectors
        //std::cout << "The eigenvalues of tbh_K_M are:" << endl << (ces.eigenvalues()).real() << endl;
        
        //std::cout << "The matrix of eigenvectors, V, is:" << endl << ces.eigenvectors() << endl << endl;
        
        ep2[loop] = (ces.eigenvalues()[0]).real();
        en2[loop] = (ces.eigenvalues()[1]).real();

        if (ep2[loop] < en2[loop]){
            ep2[loop] = (ces.eigenvalues()[1]).real();
            en2[loop] = (ces.eigenvalues()[0]).real();
        }    
       // std:: cout<< "ep2 = " << ep2[loop] << "  en2 = " << en2[loop] <<endl;
      //  std:: cout<< "en2 = " << en2[loop] <<endl;    // identified some *10 fator missing

        k_x= dx;
        k_y= k_y -dy;


        // for Gamma-to-K path  ************
        double dx1 =2*M_PI/(a*sqrt(3))/N;
        double dy1 =2*M_PI/(a*3)/N;


        complex<double> f1 = exp(rad2deg*i*k_x1*a/sqrt(3)) + 2*(cos(rad2deg*k_y1*a/2))*exp(-1*rad2deg*i*k_x1*a/(2*sqrt(3)));
        complex<double> f_c1 = conj(f1);
       // std::cout<< f1 << endl;
       // std::cout<< f_c1 << endl;
        complex<double> sd1= 1;//1/(1- (s0*s0*(abs(f1))*(abs(f1))));
        complex<double> ss1= 0;//s0*t0*(abs(f1))*(abs(f1));
        //std::cout<<"sd1:"<< sd1<< "ss1:"<<ss1 <<std::endl;    correct for both-
        //std::cout<< ss1 <<std::endl;

        Eigen::MatrixXcd tbh_G_K(2,2);  // complex dynamic matrix representing tight binding hamiltonian for monolayer graphene
        tbh_G_K(0,0)= sd1*ss1;
        tbh_G_K(1,0)= sd1*(-t0*f_c1);
        tbh_G_K(0,1)= sd1*(-t0*f1);
        tbh_G_K(1,1)= sd1*ss1;
       // std:: cout<<"fc&f_c1"<< tbh_G_K <<endl;
        //Eigen::ComplexEigenSolver<MatrixXcd> ces;
        ces.compute(tbh_G_K);
        //std::cout << "The eigenvalues of tbh_G_M are:" << endl << (ces.eigenvalues()).real() << endl;
        
        //std::cout << "The matrix of eigenvectors, V, is:" << endl << ces.eigenvectors() << endl << endl;
        ep1[loop] = (ces.eigenvalues()[0]).real();
        en1[loop] = (ces.eigenvalues()[1]).real();
        if(ep1[loop]<en1[loop]){
            ep1[loop] = (ces.eigenvalues()[1]).real();
            en1[loop] = (ces.eigenvalues()[0]).real();
        }
    //    std:: cout<< "ep1 = " << ep1[loop] << " en1 = " << en1[loop]  <<endl;  //  checked but some error may (values doesn't match), checked in ||
      //  std:: cout<< "en1 = " << en1[loop] <<endl;   // with ep2 & en2

        k_x1= k_x1 + dx1;
        k_y1= k_y1 + dy1;


    }
    
    // eigenvalues of Gamma-K  +  K-M is saved 
    FILE *fid1;
    fid1 = fopen("eigenvalues.dat","w");
  //  fprintf(fid1,"   ep1[ ]                      en1[ ]                      ep2[ ]                     en2[ ]  \n");
    for(int l=0;l<N;l++){
        
        fprintf(fid1,"  %e               %e               %e               %e\n",ep1[l], en1[l], ep2[l], en2[l]);   
          
    }
    fclose(fid1);


