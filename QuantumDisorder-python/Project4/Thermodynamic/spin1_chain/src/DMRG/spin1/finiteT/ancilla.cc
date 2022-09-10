#include "itensor/all.h"
#include "TStateObserver.h"
#include "S2.h"
#include <math.h>       /* sqrt */
#include "itensor/util/print_macro.h"

using namespace std;
using namespace itensor;

int
main(int argc, char* argv[])
    {
    //Get parameter file
    if(argc != 2)
        {
        printfln("Usage: %s inputfile.",argv[0]);
        return 0;
        }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N",10);
    auto K = input.getReal("K",0.0); // single anisotropy


    auto beta = input.getReal("beta",1);
    auto tau = input.getReal("tau",0.01);

    auto maxdim = input.getInt("maxdim",1000);
    auto cutoff = input.getReal("cutoff",1E-11);

    auto verbose = input.getYesNo("verbose",false);

    Args args;
    args.add("MaxDim",maxdim);
    args.add("Cutoff",cutoff);
    args.add("Verbose",verbose);
    args.add("Method","DensityMatrix");
    //args.add("Method","Fit");

    //auto sites = SpinHalf(2*N,{"ConserveQNs=",false});
    //auto sites = SpinOne(2*N,{"SHalfEdge=",true,"ConserveQNs=",false});
    auto sites = SpinOne(2*N,{"ConserveQNs=",false});

    /*
    j : 1 , s1 : 1 , s2 : 3
    j : 2 , s1 : 3 , s2 : 5
    j : 3 , s1 : 5 , s2 : 7
    j : 4 , s1 : 7 , s2 : 9
    j : 5 , s1 : 9 , s2 : 11
    j : 6 , s1 : 11 , s2 : 13
    j : 7 , s1 : 13 , s2 : 15
    j : 8 , s1 : 15 , s2 : 17
    j : 9 , s1 : 17 , s2 : 19
    */
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        auto s1 = 2*j-1,
             s2 = 2*j+1;
        ampo += 0.5,"S+",s1,"S-",s2;
        ampo += 0.5,"S-",s1,"S+",s2;
        ampo +=     "Sz",s1,"Sz",s2;
        //ampo += -1.0*K ,"Sz",s1,"Sz",s1;
        ampo += -1.0*K ,"Sz2",s1;

        //printfln(" j : %i , s1 : %i , s2 : %i ",j,s1,s2);
        }
    auto  j = N-1;
    auto  s1 = 2*j+1;
    ampo +=  -1.0*K ,"Sz2",s1;

    // long range
    //ampo += 0.5,"S+",1,"S-",2*N-1;
    //ampo += 0.5,"S-",1,"S+",2*N-1;
    //ampo +=     "Sz",1,"Sz",2*N-1;


    auto H = toMPO(ampo);
    auto expH = toExpH(ampo,tau);
    auto S2 = makeS2(sites,{"SkipAncilla=",true});
    //auto Sz2 = makeTotSz2(sites,{"SkipAncilla=",true});


    //
    // Make initial 'wavefunction' which is a product
    // of perfect singlets between neighboring sites
    //
    auto psi = MPS(sites);
    for(int n = 1; n <= 2*N; n += 2)
        {
        auto s1 = sites(n);
        auto s2 = sites(n+1);
        auto wf = ITensor(s1,s2);
        // https://github.com/ITensor/ITensor/blob/v3/itensor/mps/sites/spinhalf.h
        // spin-1/2   singlet state , which conserve Sz=0
        // // 1/sqrt(2)*(|01>-|10>)
        //  s1:P, s2:Q
        // s1(1)= up-spin,   s1(2)=down-spin
        //wf.set(s1(1),s2(2), ISqrt2);
        //wf.set(s1(2),s2(1), -ISqrt2);
        //wf.set(s1(1),s2(2),  1.0/sqrt(2.0));
        //wf.set(s1(2),s2(1), -1.0/sqrt(2.0));


        //https://github.com/ITensor/ITensor/blob/v3/itensor/mps/sites/spinone.h
        // spin-1  maximally entagled state
        // 1/sqrt(3)*(|11>-|00>+|-1,-1>)
        wf.set(s1(1),s2(1),  1.0/sqrt(3.0));
        wf.set(s1(2),s2(2), -1.0/sqrt(3.0));
        wf.set(s1(3),s2(3),  1.0/sqrt(3.0));

        ITensor D;
        psi.ref(n) = ITensor(s1);
        psi.ref(n+1) = ITensor(s2);
        svd(wf,psi.ref(n),D,psi.ref(n+1));
        psi.ref(n) *= D;
        }



    auto obs = TStateObserver(psi);

    auto ttotal = beta/2.;
    const int nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
    if(fabs(nt*tau-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }
    printfln("Doing %d steps of tau=%f",nt,tau);

    auto targs = args;

    auto En     = Vector(nt);
    auto Sus    = Vector(nt);
    auto SzSz   = Vector(nt);
    auto Betas  = Vector(nt);
    Real tsofar = 0;

    std::ofstream enf("en_K0.5_cutoffE-10_Method_DensityMatrix.dat");
    std::ofstream susf("sus_K0.5_cutoffE-10_Method_DensityMatrix.dat");
    for(int tt = 1; tt <= nt; ++tt)
        {
        psi = applyMPO(expH,psi,args);

        //Normalize wavefunction
        psi.ref(1) /= norm(psi(1));
        //psi.noPrime().normalize();

        tsofar += tau;
        targs.add("TimeStepNum",tt);
        targs.add("Time",tsofar);
        targs.add("TotalTime",ttotal);
        obs.measure(targs);

        //Record beta value
        auto bb = (2*tsofar);
        Betas(tt-1) = bb;

        //
        // Measure Energy
        //
        auto en = inner(psi,H,psi);
        En(tt-1) = en/N;

        //
        // Measure Susceptibility
        //
        auto s2val = inner(psi,S2,psi);
        Sus(tt-1) = (s2val*bb/3.)/N;

        //
        // Measuring <SziSzj>
        //


        /*
        //http://itensor.org/docs.cgi?page=formulas/two_mps
        auto j = int(N/2);
        auto SziSzj = 0.0;
        for(int i = 2; i <N; ++i)
            {
              auto op_i = op(sites,"Sz",i);
              auto op_j = op(sites,"Sz",j);

              auto psidag = dag(psi);
              auto M = psi(1)*psidag(1);
              for(auto n : range1(2,N))
                  {
                  M *= psi(n);
                  if(n == i)
                      {
                      M *= op_i*prime(psidag(i),"Site");
                      }
                  else if(n == j)
                      {
                      M *= op_j*prime(psidag(j),"Site");
                      }
                  else
                      {
                      M *= psidag(n);
                      }
                  }
            SziSzj += elt(M);
          }
        */

         /*
        // http://itensor.org/docs.cgi?vers=cppv3&page=formulas/correlator_mps
        auto SziSzj = 0.0;
        auto i = 2;
        for(int j = 3; j <N; ++j)
            {

              //Make the operators you want to measure
              auto op_i = op(sites,"Sz",i);
              auto op_j = op(sites,"Sz",j);

              //'gauge' the MPS to site i
              //any 'position' between i and j, inclusive, would work here
              //
              // This allows for the left-normalized (purple) blocks of the MPS to give the identity
              // on contraction till site i, so we don't need to program anything special to deal
              //  with these
              psi.position(i);
              //Create the bra/dual version of the MPS psi
              auto psidag = dag(psi);


              //Prime the link indices to make them distinct from
              //the original ket links
              psidag.prime("Link");

              //index linking i-1 to i:
              auto li_1 = leftLinkIndex(psi,i);

              auto C = prime(psi(i),li_1) * op_i;
              C *= prime(psidag(i),"Site");
              for(int k = i+1; k < j; ++k)
                  {
                      C *= psi(k);
                      C *= psidag(k);
                  }
              //index linking j to j+1:
              auto lj = rightLinkIndex(psi,j);

              C *= prime(psi(j),lj) * op_j;
              C *= prime(psidag(j),"Site");

              SziSzj += elt(C);

            }
        */

        SzSz(tt-1) = 0.0;//(SziSzj*bb/3.)/N;
        printfln("\nEnergy/N %.4f %.20f, Sus/N %.20f, SzSz/N %.20f",bb,en/N,Sus(tt-1),SzSz(tt-1));
        println();
        enf  << format("%.14f %.14f %.14f\n",Betas(tt-1), 1.0/Betas(tt-1), En(tt-1));
        susf << format("%.14f %.14f %.14f\n",Betas(tt-1), 1.0/Betas(tt-1), Sus(tt-1), SzSz(tt-1));
        }

    //std::ofstream enf("en_K1.0_cutoffE-12.dat");
    //std::ofstream susf("sus_K0.5_cutoffE-12.dat");
    //for(auto n : range(Betas))
    //    {
    //    enf << format("%.14f %.14f %.14f\n",Betas(n), 1.0/Betas(n), En(n));
    //    susf << format("%.14f %.14f %.14f\n",Betas(n),1.0/Betas(n),Sus(n));
    //    }
    enf.close();
    susf.close();

    return 0;
    }
