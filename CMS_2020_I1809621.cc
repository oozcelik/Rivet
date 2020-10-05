// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2020_I1809621 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2020_I1809621);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");

      // Book histograms
      book (_h_eta1, 10, -2.4, 2.4);
      book (_h_eta2, 10, -2.4, 2.4);
      book (q2, 10, 0, 20);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here

/*      vector<double> val_dNdEta;
      vector<double> val_dNdPt;
      vector<double>q2;
      val_dNdEta.clear();
      val_dNdPt.clear();*/
      const FinalState& cfs = apply<FinalState>(event, "FS");
      for (const Particle& p1 : cfs.particles()) {
        int m1id = p1.pid();
	if (m1id != 13) continue;        
       
      const double m1 = p1.mass();
      const double eta = p1.momentum().eta();
      const double pt = p1.momentum().pT();
      _h_eta1->fill(eta);
   
      for (const Particle& p2 : cfs.particles()) {
        int m2id = p2.pid();
        if (m1id != 13 && m2id != -13 && p1 == p2) continue;

         const double m2 = p2.mass();
         const double eta2 = p2.momentum().eta();
         const double pt2 = p2.momentum().pT();
         _h_eta2->fill(eta2); 
         q2->fill(m1+m2);    
     }
   }

    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_eta1); // normalize to unity
      normalize(_h_eta2);    
      normalize(q2);
    //  scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_eta1, _h_eta2, q2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2020_I1809621);


}
