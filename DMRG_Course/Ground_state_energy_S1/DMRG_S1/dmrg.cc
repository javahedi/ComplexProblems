#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <sstream>


using namespace itensor;
using namespace std;


int 
main()
    {
    int N      = 20;
    auto delta = 1.0;
    auto B     = 0.0;
    auto J     = 1.0;
    auto sites = SpinOne(N,{"ConserveQNs=",false});

    auto ampo = AutoMPO(sites);
    for(int j = 1; j <N; ++j)
        {
            //cout<<"N="<<j<<endl;
            ampo += J,"Sx",j,"Sx",j+1;
            ampo += J,"Sy",j,"Sy",j+1;
            ampo += J,"Sz",j,"Sz",j+1;
            //ampo += B,"Sx",j;
        }
    //ampo += B,"Sx",N;
        
    //for(int j = 1; j <N-1; ++j)
    //    {
    //        ampo += delta,"Sz",j,"Sz",j+2;
    //    }

    auto H = toMPO(ampo);

    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,50,200,400,600;
    sweeps.cutoff() = 1E-10;

    auto psi0 = randomMPS(sites);

    auto [energy,psi1] = dmrg(H,psi0,sweeps,{"Quiet=",true});
     
    println("Ground state energy per spin : ", energy/N);
    //
    
    
    return 0;
    }

