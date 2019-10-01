#include "itensor/all.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto h = input.getReal("h");
    auto nsweeps = input.getInt("nsweeps");
    auto quiet = input.getYesNo("quiet",true);
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    auto pi = M_PI;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        auto coeff = pow(sin(pi/N*j),2.0);
        ampo += 0.5 * coeff,"S+",j,"S-",j+1;
        ampo += 0.5 * coeff,"S-",j,"S+",j+1;
        ampo +=       coeff,"Sz",j,"Sz",j+1;
        }
    for(int j = 1; j <= N; ++j)
        {
        auto coeff = pow(sin(pi/N*(j-0.5)),2.0);
        ampo += -h * coeff,"Sz",j;
        }
    auto H = toMPO(ampo);

    auto psi0 = randomMPS(sites);

    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet",quiet});

    println("\n");
    printfln("# ene %.10f",energy/N);
    println("\n");

    auto mzave = 0.0;
    for(int j = N/2-1; j <= N/2+2; ++j)
        {
        psi.position(j);
        auto wfave = psi.A(j);
        mzave += (dag(prime(wfave,"Site")) * sites.op("Sz",j) * wfave).real();
        }
    printfln("# avemag %.10f",mzave/4.0);
    println("\n");

    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        auto wf = psi.A(j);
        auto mz = (dag(prime(wf,"Site")) * sites.op("Sz",j) * wf).real();
        printfln("# mag %4d %.10f",j,mz);
        }
    println("\n");

    return 0;
    }
