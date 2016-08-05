// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"



namespace Rivet {


  class CDF_2013_I1225278 : public Analysis {
  public:

    /// Constructor
    CDF_2013_I1225278()
      : Analysis("CDF_2013_I1225278")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // final state
      FinalState fs;
      addProjection(fs, "FS");

      /// leading photon
      LeadingParticlesFinalState photonfs(FinalState(-1.0, 1.0, 30.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      addProjection(photonfs, "LeadingPhoton");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      addProjection(vfs, "JetFS");
     
      /// jet final state  TODO: maximal eta range?
      addProjection(FastJets(vfs, FastJets::CDFJETCLU, 0.4), "Jets");

      /// booking histograms
       _h_Et_photon_bottom = bookHisto1D(1, 1, 1);
       _h_Et_photon_charm = bookHisto1D(2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      bool check = false; //flag to stop iterating over jets
      
      // get leading photon
      Particles photons = applyProjection<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size()==0)   vetoEvent;
      FourMomentum ph = photons[0].momentum();
      
      //check photon isolation
      Particles vfs = applyProjection<VetoedFinalState>(event, "JetFS").particles();
      FourMomentum mom_in_cone;
      foreach (const Particle& p, vfs){
         if (deltaR(ph, p.momentum())< 0.4) mom_in_cone +=p.momentum();
      }
      if (mom_in_cone.E() > 2.0*GeV) vetoEvent;
 
      //get jets
      Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(20*GeV && Cuts::abseta < 1.5);
			if ( jets.empty() ) vetoEvent;
      
     /// loop over jets, check if it contains bottom /charms and if the photon is outside the jet
     // break loop, if a good jet is found
      for (Jets::const_iterator jt = jets.begin(); jt!=jets.end(); ++jt){
         if (deltaR(jt->momentum(), ph)>0.4){
            if (jt->containsCharm() ){
               check=true;
               _h_Et_photon_charm->fill(ph.pT(), weight);
               }
            if (jt->containsBottom() ){
               check=true;
               _h_Et_photon_bottom->fill(ph.pT(), weight);
               }
         }
         if (check==true) break;
         
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

       scale(_h_Et_photon_charm, crossSection()/sumOfWeights()); // norm to cross section
       scale(_h_Et_photon_bottom, crossSection()/sumOfWeights()); // norm to cross section
      

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


    /// @name Histograms
    //@{
    Histo1DPtr _h_Et_photon_charm;
    Histo1DPtr _h_Et_photon_bottom;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2013_I1225278);


}
