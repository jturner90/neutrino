#include <fstream>

#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

//pythia header files

#include "tauola_print_parameters.h"
using namespace std;
using namespace Tauolapp;

#include <fstream>
#include <boost/algorithm/string.hpp>
#include "Eigen/Core"
#include "opts.h"
#include "tools.h"


std::tuple<Eigen::Vector4d, Eigen::Vector4d, Eigen::Vector4d> mkParticles(Eigen::VectorXd const & DIN) {
    Eigen::Vector4d nu, tau, w;
    tau << DIN[1], DIN[2], DIN[3], DIN[4];
    nu << DIN[0], 0, 0, DIN[0];
    //w =  tau - nu;
    w = tau + nu;
    return std::make_tuple(nu,tau,w);
    
}

inline HepMC::FourVector eigenToHepMC(Eigen::Vector4d const & vin) {
    HepMC::FourVector temp(vin[1], vin[2], vin[3], vin[0]);
    return temp;
}

int main(int argc, char * argv[]) {
    int mode(0);
    int nmax(-1);
    string fin = "";
    string foutname = "testout.csv";
    using namespace opts;
    Options ops(argc, argv);
    bool norad     = ops >> Present("norad", "switch off photon radiation");
    bool nospincorr  = ops >> Present("nospincorr", "switch off spin correlation radiation");
    bool verbose     = ops >> Present('v', "verbose", "verbose output");
    bool tauplus     = ops >> Present('p', "tauplus", "Treat taus as tau+");
    ops >> Option('m', "mode",   mode,   "Tauolo decay mode");
    ops >> Option('n', "nmax",   nmax,   "Max number of events to generate");
    ops >> Option('f', "ftau",   fin,   "Max number of events to generate");
    ops >> Option('o', "fout",   foutname,   "Max number of events to generate");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    if (fin ==foutname) {
      std::cerr << "Do not overwrite input file, exiting\n";
      exit(1);
    }

    std::cout<< "Reading taus from " << fin << "\n"; 
    std::cout<< "Running Tauola mode " << mode << "\n"; 
    if (norad) std::cerr << "Photon radiation turned OFF\n";
    else std::cerr << "Photon radiation turned on\n";
    auto const data = readCSV(fin.c_str());
    auto [eventids, taus] = convert(data);

    int nevents = taus.rows();
    if ( (nmax>0) &&  (nmax < nevents)) {
      nevents = nmax;
    }
    std::cout << "Processing " << nevents << " tau events\n";


    // Enforce a single decay mode
    Tauola::setSameParticleDecayMode(mode);
    Tauola::setOppositeParticleDecayMode(mode);
    Tauola::setRadiation(!norad);
    Tauola::setEtaK0sPi(1,1,1);

    Tauola::initialize();
    Tauola::spin_correlation.setAll(!nospincorr);

    if (verbose) tauola_print_parameters(); // Prints TAUOLA  parameters (residing inside its library): e.g. to test user interface
    
    //const double tau_mass = parmas_.amtau;
    
    std::fstream fout;
    // Open the file in truncate mode if first line else in Append Mode
    fout.open(foutname, std::ios::out | (std::ios::trunc));

    fout << "EVENTID,PDGID,E,PX,PY,PZ\n";

    for (int iEvent = 0; iEvent < nevents; ++iEvent) {
        HepMC::GenEvent * event = new HepMC::GenEvent();
       
        //std::cerr <<  taus.row(iEvent) << "\n";
        auto [nu, tau, w] = mkParticles(taus.row(iEvent));
        //abort();

        auto mom_nu  = eigenToHepMC(nu);
        auto mom_tau = eigenToHepMC(tau);
        auto mom_w   = eigenToHepMC(w);
        HepMC::GenParticle * p_tau;
        HepMC::GenParticle * p_nu   = new HepMC::GenParticle(mom_nu,  16, 4); // Incident neutrino
        if (tauplus)
        {
          p_tau  = new HepMC::GenParticle(mom_tau, -15, 1); // Positively charged Tau 
        }
        else
        {
          p_tau  = new HepMC::GenParticle(mom_tau,  15, 1); // Negatively charged Tau 
        }
        HepMC::GenParticle * p_w    = new HepMC::GenParticle(mom_w,   24, 4); // The W


        HepMC::GenVertex * vertex = new HepMC::GenVertex();
        
        //vertex->add_particle_in(p_nu);
        //vertex->add_particle_in(p_w);
        //vertex->add_particle_out(p_tau);

        vertex->add_particle_out(p_nu);
        vertex->add_particle_in(p_w);
        vertex->add_particle_out(p_tau);

        event->add_vertex(vertex);
        if (verbose) {
          cout << "BEFORE:"<<endl;
          event->print();
        }
        TauolaHepMCEvent * t_event = new TauolaHepMCEvent(event);
        t_event->decayTaus();


        for (auto v = event->vertices_begin(); v != event->vertices_end(); ++v ) {
          HepMC::GenVertex::particles_out_const_iterator p1;
          for (auto p1 = (*v)->particles_out_const_begin();p1 != (*v)->particles_out_const_end(); ++p1 ) {
            auto pout=**p1;
            if (pout.status() !=1) continue; // Skip all non-final particles
            fout << eventids[iEvent] << "," << pout.pdg_id()
              << "," << pout.momentum().e() 
              << "," << pout.momentum().px()
              << "," << pout.momentum().py()
              << "," << pout.momentum().pz() << "\n";
          }
        }
        if (verbose) { 
          cout << "AFTER:"<<endl;
          event->print();
        }

        delete event;
        delete t_event;
    }
    fout.close();

    return 0;
}
