#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace itensor;
using namespace std;



// function returning theta for each K
float thetaK(double k, double g) 
{
   // local variable declaration
   double result;
   result = atan(sin(k)/(g-cos(k)))/2;
 
   return result; 
}





int
main()
    {
    auto PI = 3.141592653589793238463;
    int N = 800;
    auto J = 1.0;
    auto g0 = 1.5;
    auto g1 = 0.2;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
   

    //Initial Hamiltonian before Quench
    auto ampo0 = AutoMPO(sites);
    for(int j = 1; j <N; ++j)
        {
	    ampo0 += 4.0*J,"Sz",j,"Sz",j+1;
        ampo0 += 2*g0*J,"Sx",j;
        }
    ampo0 += 2*g0*J,"Sx",N;
    auto H = toMPO(ampo0);
    


    //Find ground state before quench
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,50,100,200,300;
    sweeps.cutoff() = 1E-10;

    auto psi = randomMPS(sites);
    auto [energy,psi0] = dmrg(H,psi,sweeps,{"Quiet=",true});
    
    printfln("\n E0/N : %.12f ", energy/N);

   

    
    ////////////////////////////////// 
    //final Hamiltonian after quench
    auto ampo = AutoMPO(sites);
    for(int j = 1; j <N; ++j)
        {
        //ampo += 0.5,"S+",j,"S-",j+1;
        //ampo += 0.5,"S-",j,"S+",j+1;
        ampo += 4.0*J, "Sz",j,"Sz",j+1;
        ampo += 2*g1*J,"Sx",j;
	}
    ampo += 2*g1*J,"Sx",N;
    ////////////////////////////////////
    




    auto tau = 0.02; //Time step
    //Approx. with error O(tau^4) [thanks to Kemal]
    auto t11 = 0.10566243270259355887-0.39433756729740644113*Cplx_i;
    auto t22 = Cplx_i*t11;
    auto t33 = conj(t22);
    auto t44 = Cplx_i*t33;

    auto t1 = t11*Cplx_i*tau;
    auto t2 = t22*Cplx_i*tau;
    auto t3 = t33*Cplx_i*tau;
    auto t4 = t44*Cplx_i*tau;

    auto expH1 = toExpH(ampo,t1);
    auto expH2 = toExpH(ampo,t2);
    auto expH3 = toExpH(ampo,t3);
    auto expH4 = toExpH(ampo,t4); 
     
    
    //////////////////////////////////
    auto args = Args("Method=","Fit","Cutoff=",1E-10,"MaxDim=",1000);
    auto ttotal = 5.0;
    auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
    

 
    ofstream OUT;
    OUT.precision(16);
    OUT.open("return_probabiliy_N_800_chi_1000_Cutoff_1E-10.txt");
    OUT <<"# time       return_probability   rate_function_DMRG     rate_function_Exact"<<endl;
    
    auto RP=0.0;
    auto RF=0.0;
    psi=psi0;

    for(int n = 0; n <= nt; ++n) // loop over time
        { 

        //////////////////////
        //  Excat calculation
        ///////////////////////
        auto  gamaT=0.0;
        auto  Ek_g1=0.0;
        auto  Phi_k=0.0;
        for(int i = 0; i < N; ++i) // loop over K
                {
                 Ek_g1 = 2*J*sqrt(sqr(g1-cos(i*2*PI/N))+sqr(sin(i*2*PI/N)));
                 Phi_k = thetaK(i*2*PI/N,g0)-thetaK(i*2*PI/N,g1);
                 gamaT += log(abs(sqr(cos(Phi_k))+sqr(sin(Phi_k))*exp(-2*Cplx_i*n*tau*Ek_g1)));

                }

        gamaT = -gamaT/N; //-2*gamaT/N;
        /////////////////////////////////////////

        RP = abs(innerC(psi0,psi));
	RF = -log(sqr(RP))/N;
        OUT <<n*tau<<"		"<<RP<<"	"<<RF<<"	"<<gamaT<<endl;
        printfln("Tim : %.2f, RP :  %.8f,  RF : %.8f,  Exact : %.8f",n*tau,RP,RF,gamaT);
        
        //printfln("\n return probability : %.12f ", abs(innerC(psi0,psi)));
        //printfln("\n return probability :  %.12f ", abs(overlapC(psi0,psi)));
        //printfln("\n ====================================");


        psi = applyMPO(expH1,psi,args);
        psi.noPrime().normalize();

        psi = applyMPO(expH2,psi,args);
        psi.noPrime().normalize();

        psi = applyMPO(expH3,psi,args);
        psi.noPrime().normalize();

        psi = applyMPO(expH4,psi,args);
	psi.noPrime().normalize();
   	}

    OUT.close();
    return 0;
   } 
