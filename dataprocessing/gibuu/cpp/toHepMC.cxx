#include <iostream>
#include <algorithm>

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/IO_GenEvent.h"

#include "opts.h"
#include "tools.h"


inline HepMC::FourVector eigenToHepMC(Eigen::Vector4d const & vin) {
    HepMC::FourVector temp(vin[1], vin[2], vin[3], vin[0]);
    return temp;
}


void processSignalEvents(std::string ftau, std::string fdec, std::string fcor, std::string _fout, bool tauplus) {
    auto const tau_data = readCSV(ftau);
    auto [tau_eventids, tau_array] = convert(tau_data);
    
    auto const decay_data = readCSV(fdec);
    auto [decay_eventids, decay_array] = convert(decay_data);
    
    auto const oth_data = readCSV(fcor);
    auto [oth_eventids, oth_array] = convert(oth_data, 8);

    auto map_decay = indexMap(decay_eventids);
    auto map_oth   = indexMap(oth_eventids);
    auto map_tau   = indexMap(tau_eventids);


    int nEvents = tau_eventids.size();//-2;

    std::ofstream fout;
    fout.open(_fout);

    HepMC::IO_GenEvent * iogen = new HepMC::IO_GenEvent(fout);
    // Event loop
    for (unsigned int i=0;i<nEvents-1;++i) { // TODO currently skipping last two events because map algorithm not stable enough

        auto id = tau_eventids[i];
        auto bla = map_decay.find(id);
        if (bla == map_decay.end()) {
          std::cerr << "No more work\n";
          break;
        }
        
        HepMC::GenVertex * vtx_prod  = new HepMC::GenVertex();
        HepMC::GenVertex * vtx_tau= new HepMC::GenVertex();

        HepMC::GenParticle * ptau;
        if (tauplus)
        {
          ptau  = new HepMC::GenParticle(eigenToHepMC(tau_array.row(i).segment(1,4)), -15, 2); // Positively charged Tau 
        }
        else
        {
          ptau  = new HepMC::GenParticle(eigenToHepMC(tau_array.row(i).segment(1,4)),  15, 2); // Negatively charged Tau 
        }

        vtx_prod->add_particle_out(ptau);
        vtx_tau->add_particle_in( ptau);

        for (auto idec : map_decay[id]) {
            HepMC::GenParticle * p = new HepMC::GenParticle(eigenToHepMC(decay_array.row(idec).segment(1,4)), decay_array.row(idec)[0], 1);
            vtx_tau->add_particle_out(p);
        }

        double weight(0), genxs(0), Enu(0);
        for (auto ioth : map_oth[id]  ) {
            HepMC::GenParticle * p = new HepMC::GenParticle(eigenToHepMC(oth_array.row(ioth).segment(1,4)), oth_array.row(ioth)[0], 1);
            vtx_prod->add_particle_out(p);
            weight= oth_array.row(ioth)[6];
            genxs = oth_array.row(ioth)[7];
            Enu = oth_array.row(ioth)[5];
        }
        HepMC::GenEvent   * event  = new HepMC::GenEvent();
        event->set_event_number(i);
        event->set_event_scale(Enu);

        event->add_vertex(vtx_prod);
        event->add_vertex(vtx_tau);

        HepMC::GenCrossSection xsec;
        xsec.set_cross_section( genxs,0);
        event->set_cross_section(xsec);

        HepMC::WeightContainer wc;
        wc["Weight"] = weight;
        event->weights()=wc;
        iogen->write_event(event);

        delete event;
    }

    delete iogen;
    fout.close();
}

void processBGEvents(std::string fcor, std::string _fout) {
    auto const oth_data = readCSV(fcor);
    auto [oth_eventids, oth_array] = convert(oth_data, 8);

    auto map_oth   = indexMap(oth_eventids);

    std::ofstream fout;
    fout.open(_fout);

    HepMC::IO_GenEvent * iogen = new HepMC::IO_GenEvent(fout);

    size_t eventNumber = 1;
    for( auto const& [evtId, evtData] : map_oth ) {

      if (evtData.size()<2) {
        continue;
      }

        
        HepMC::GenVertex * vtx_prod  = new HepMC::GenVertex();

        double weight(0), genxs(0), Enu(-1);
        for (auto ioth : evtData  ) {
            HepMC::GenParticle * p = new HepMC::GenParticle(eigenToHepMC(oth_array.row(ioth).segment(1,4)), oth_array.row(ioth)[0], 1);
            vtx_prod->add_particle_out(p);
            weight= oth_array.row(ioth)[6];
            genxs = oth_array.row(ioth)[7];
            Enu = oth_array.row(ioth)[5];
        }


        HepMC::GenEvent   * event  = new HepMC::GenEvent();
        event->set_event_number(eventNumber);
        event->set_event_scale(Enu);

        event->add_vertex(vtx_prod);

        HepMC::GenCrossSection xsec;
        xsec.set_cross_section( genxs,0);
        event->set_cross_section(xsec);

        HepMC::WeightContainer wc;
        wc["Weight"] = weight;
        event->weights()=wc;
        iogen->write_event(event);

        delete event;
        eventNumber++;
    }

    delete iogen;
    fout.close();
}


// 3 input files:
// * taus
// * tau decays
// * other
// 1 output file for hepmc
int main(int argc, char * argv[]) {
    int nmax(-1);
    std::string fcor = "";
    std::string ftau = "";
    std::string fdec = "";
    std::string _fout = "output.hepmc2g";
    using namespace opts;
    Options ops(argc, argv);
    bool verbose     = ops >> Present('v', "verbose", "verbose output");
    bool quiet       = ops >> Present('q', "quiet", "least verbose mode");
    bool tauplus     = ops >> Present('p', "tauplus", "Treat taus as tau+");
    ops >> Option('n', "nmax",   nmax,   "Max number of events");
    ops >> Option('c', "fcor",   fcor,   "Input file with hadronic remnant");
    ops >> Option('t', "ftau",   ftau,   "Input file UNDECAYED taus");
    ops >> Option('d', "fdec",   fdec,   "Input file DECAYED taus (i.e. tauola output)");
    ops >> Option('o', "fout",   _fout,   "Output file");
    if (ops >> Present('h', "help", "Show help")) {
        std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    if ( (ftau != "") && (fcor != "") && (fdec != "") ) {
      processSignalEvents(ftau, fdec, fcor, _fout, tauplus);
    }
    else if ( (ftau == "") && (fcor != "") && (fdec == "") ) {
      processBGEvents(fcor, _fout);
    }
    else {
        std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    return 0;
}
