#include <iomanip>
#include "itensor/all.h"
#include "Timer.h"
Timers timer;
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "MPSUtility.h"
#include "SingleParticle.h"
#include "MixedBasis.h"
#include "SystemStruct.h"
#include "ContainerUtility.h"
#include "Corr.h"
#include "TDVPObserver.h"
#include "MeaCurrent.h"
#include "tdvp.h"
#include "basisextension.h"
#include "InitState.h"
#include "Tools.h"
using namespace vectool;
using namespace itensor;
using namespace std;
using namespace iutility;

struct Para
{
    vector<tuple<string,string,int,int,Real>> hops;
    Real Ec, Ng, Delta;

    void write (ostream& s) const
    {
        iutility::write(s,hops);
        iutility::write(s,Ec);
        iutility::write(s,Ng);
        iutility::write(s,Delta);
    }

    void read (istream& s)
    {
        iutility::read(s,hops);
        iutility::read(s,Ec);
        iutility::read(s,Ng);
        iutility::read(s,Delta);
    }
};

void writeAll (const string& filename,
               const MPS& psi, const MPO& H,
               const WireSystem& system,
               const Para& para,
               int step)
{
    ofstream ofs (filename);
    itensor::write (ofs, psi);
    itensor::write (ofs, H);
    itensor::write (ofs, step);
    para.write (ofs);
    system.write (ofs);
}

void readAll (const string& filename,
              MPS& psi, MPO& H,
              WireSystem& system,
              Para& para,
              int& step)
{
    ifstream ifs = open_file (filename);
    itensor::read (ifs, psi);
    itensor::read (ifs, H);
    itensor::read (ifs, step);
    para.read (ifs);
    system.read (ifs);
}

int main(int argc, char* argv[])
{
    string infile = argv[1];
    InputGroup input (infile,"basic");

    auto L_lead   = input.getInt("L_lead");
    auto L_device   = input.getInt("L_device");
    auto t_lead     = input.getReal("t_lead");
    auto t_device   = input.getReal("t_device");
    auto t_contactL = input.getReal("t_contactL");
    auto t_contactR = input.getReal("t_contactR");
    auto mu_leadL   = input.getReal("mu_leadL");
    auto mu_leadR   = input.getReal("mu_leadR");
    auto mu_device  = input.getReal("mu_device");
    auto mu_biasL   = input.getReal("mu_biasL");
    auto mu_biasS   = input.getReal("mu_biasS");
    auto mu_biasR   = input.getReal("mu_biasR");
    auto Delta      = input.getReal("Delta");
    auto Ec         = input.getReal("Ec");
    auto Ng         = input.getReal("Ng");
    auto damp_decay_length = input.getInt("damp_decay_length",10000000);
    auto maxCharge  = input.getInt("maxCharge");

    auto dt            = input.getReal("dt");
    auto time_steps    = input.getInt("time_steps");
    auto NumCenter     = input.getInt("NumCenter");
    auto Truncate      = input.getYesNo("Truncate");
    auto mixNumCenter  = input.getYesNo("mixNumCenter",false);
    auto globExpanN          = input.getInt("globExpanN",std::numeric_limits<int>::max());
    auto globExpanItv        = input.getInt("globExpanItv",1);
    auto globExpanCutoff     = input.getReal("globExpanCutoff",1e-8);
    auto globExpanKrylovDim  = input.getInt("globExpanKrylovDim",3);
    auto globExpanHpsiCutoff = input.getReal("globExpanHpsiCutoff",1e-8);
    auto globExpanHpsiMaxDim = input.getInt("globExpanHpsiMaxDim",300);
    auto globExpanMethod     = input.getString("globExpanMethod","DensityMatrix");

    auto UseSVD        = input.getYesNo("UseSVD",true);
    auto SVDmethod     = input.getString("SVDMethod","gesdd");  // can be also "ITensor"
    auto WriteDim      = input.getInt("WriteDim");

    auto write         = input.getYesNo("write",false);
    auto write_dir     = input.getString("write_dir",".");
    auto write_file    = input.getString("write_file","");
    auto read          = input.getYesNo("read",false);
    auto read_dir      = input.getString("read_dir",".");
    auto read_file     = input.getString("read_file","");

    auto sweeps        = iut::Read_sweeps (infile, "sweeps");

    MPS psi;
    MPO H;
    // Define 
    auto system = WireSystem();
    int step = 1;
    auto sites = MixedBasis();
    Para para;
    // Initialization
    if (!read)
    {
        // Factor for exponential-decay hopping
        Real damp_fac = exp(-1./damp_decay_length);
        // Single-particle Hamiltonians
        cout << "H left lead" << endl;
        auto H_leadL = Hamilt_k (L_lead, t_lead, mu_leadL, damp_fac, true, true);
        cout << "H right lead" << endl;
        auto H_leadR = Hamilt_k (L_lead, t_lead, mu_leadR, damp_fac, false, true);
        cout << "H dev" << endl;
        auto H_dev   = Hamilt_k (L_device, t_device, mu_device);
        auto H_zero  = Matrix(1,1);

        // WireSystem
        system.add_chain ("L",H_leadL);
        system.add_chain ("R",H_leadR);
        system.add_chain ("S",H_dev);
        system.add_chain ("C",H_zero);
        //system.sort_basis (sort_by_energy_S_middle (system.part("S"), {system.part("L"), system.part("R")}));
        system.sort_basis (sort_by_energy_S_middle_charging (system.parts().at("S"), system.parts().at("C"), {system.parts().at("L"), system.parts().at("R")}));
        system.print_orbs();
        cout << "device site = " << system.idevL() << " " << system.idevR() << endl;

        // SiteSet
        int N = system.N();
        int charge_site = system.to_glob ("C",1);
        // Find sites in S
        vector<int> scatter_sites;
        for(int i = 0; i < system.orbs().size(); i++)
        {
            if (get<0>(system.orbs().at(i)) == "S")
                scatter_sites.push_back (i+1);
        }
        // Make SiteSet
        Args args_basis = {"MaxOcc",maxCharge,"SC_scatter",(Delta != 0.)};
        sites = MixedBasis (N, scatter_sites, charge_site, args_basis);
        cout << "charge site = " << charge_site << endl;

        // Make Hamiltonian MPO
        para.Ec = Ec;   para.Ng = Ng;   para.Delta = Delta;
        para.hops.clear();
        para.hops.emplace_back ("L","S",-1,1,t_contactL);
        para.hops.emplace_back ("R","S",1,-1,t_contactR);
        auto ampo = get_ampo (system, sites, para);
        H = toMPO (ampo);
        cout << "MPO dim = " << maxLinkDim(H) << endl;

        // Initialze MPS
        if (Delta != 0.)
        {
            auto DMRG_sweeps = iut::Read_sweeps (infile, "DMRG_sweeps");
            psi = get_ground_state_SC (system, sites, mu_biasL, mu_biasS, mu_biasR, para, DMRG_sweeps, args_basis);
        }
        else
            psi = get_non_inter_ground_state (system, sites, mu_biasL, mu_biasS, mu_biasR);
        psi.position(1);
    }
    else
    {
        readAll (read_dir+"/"+read_file, psi, H, system, para, step);
        sites = MixedBasis (siteInds(psi));
    }
    // ======================= Time evolution ========================
    // Observer
    auto obs = TDVPObserver (sites, psi, {"charge_site",system.to_glob ("C",1)});
    // Current MPO
    int lenL = system.parts().at("L").L();
    int lenS = system.parts().at("S").L();
    vector<int> spec_links = {lenL-1, lenL+lenS+1};
    /*for(int i = 1; i < system.N_phys(); i++)
        spec_links.push_back (i);*/
    int N = system.N();
    vector<MPO> JMPOs (N);
    for(int i : spec_links)
    {
        auto mpo = get_current_mpo (sites, system, i);
        JMPOs.at(i) = mpo;
    }
    // Correlation
    /*if (SubCorrN < 0 or SubCorrN > length(psi))
    {
        SubCorrN = length(psi);
    }
    else if (SubCorrN < idevR-idevL+1)
    {
        cout << "sub-correlation size is too small: " << SubCorrN << endl;
        cout << "device size = " << idevL-idevR+1 << endl;
        throw;
    }
    int l = (N - SubCorrN) / 2;
    int ibeg = l+1,
        iend = l+SubCorrN;
    auto sub_corr = SubCorr (system, ibeg, iend);*/

    // Time evolution
    cout << "Start time evolution" << endl;
    cout << sweeps << endl;
    psi.position(1);
    Real en, err;
    Args args_tdvp_expansion = {"Cutoff",globExpanCutoff, "Method","DensityMatrix",
                                "KrylovOrd",globExpanKrylovDim, "DoNormalize",true, "Quiet",true};
    Args args_tdvp  = {"Quiet",true,"NumCenter",NumCenter,"DoNormalize",true,"Truncate",Truncate,
                       "UseSVD",UseSVD,"SVDmethod",SVDmethod,"WriteDim",WriteDim,"mixNumCenter",mixNumCenter};
    LocalMPO PH (H, args_tdvp);
    while (step <= time_steps)
    {
        cout << "step = " << step << endl;

        // Subspace expansion
        if (maxLinkDim(psi) < sweeps.mindim(1) or (step < globExpanN and (step-1) % globExpanItv == 0))
        {
            timer["glob expan"].start();
            addBasis (psi, H, globExpanHpsiCutoff, globExpanHpsiMaxDim, args_tdvp_expansion);
            PH.reset();
            timer["glob expan"].stop();
        }
        // Time evolution
        timer["tdvp"].start();
        //tdvp (psi, H, 1_i*dt, sweeps, obs, args_tdvp);
        TDVPWorker (psi, PH, 1_i*dt, sweeps, obs, args_tdvp);
        timer["tdvp"].stop();
        auto d1 = maxLinkDim(psi);

        // Measure currents by correlations
        /*timer["current corr"].start();
        sub_corr.measure<MixedBasis> (psi, {"Cutoff",corr_cutoff});
        timer["current corr"].stop();
        for(int j = 1; j < system.N_phys(); j++)
        {
            timer["current corr"].start();
            auto J = sub_corr.get_current (j);
            timer["current corr"].stop();
            cout << "\t*current " << j << " " << j+1 << " = " << J << endl;
        }*/
        // Measure currents by MPO
        for(int j : spec_links)
        {
            timer["current mps"].start();
            auto J = get_current (JMPOs.at(j), psi);
            timer["current mps"].stop();
            cout << "\t*I " << j << " " << j+1 << " " << J << endl;
        }

        Real NN = den (sites, psi, 1, system.idevL()-1);
        NN += den (sites, psi, system.idevL()+1, length(sites));
        NN += den (sites, psi, system.to_glob ("C",1));
        cout << "tot N = " << NN << endl;

        step++;
        if (write)
        {
            timer["write"].start();
            writeAll (write_dir+"/"+write_file, psi, H, system, para, step);
            timer["write"].stop();
        }
    }
    timer.print();
    return 0;
}
