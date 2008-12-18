/* \class ZMuMu_MCanalyzer
 * 
 * author: Davide Piccolo
 *
 * ZMuMu MC analyzer:
 * check muon reco efficiencies from MC truth, 
 *
 */

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h" 
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <vector>
using namespace edm;
using namespace std;
using namespace reco;

class ZMuMu_MCanalyzer : public edm::EDAnalyzer {
public:
  ZMuMu_MCanalyzer(const edm::ParameterSet& pset);
private:
  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup);
  bool check_ifZmumu(const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2); 
  float getParticlePt(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2); 
  float getParticleEta(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2); 
  float getParticlePhi(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2); 
  Particle::LorentzVector getParticleP4(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2); 
  virtual void endJob();

  edm::InputTag zMuMu_, zMuMuMatchMap_; 
  edm::InputTag zMuStandAlone_, zMuStandAloneMatchMap_;
  edm::InputTag zMuTrack_, zMuTrackMatchMap_; 
  edm::InputTag muons_, muonMatchMap_, muonIso_;
  edm::InputTag tracks_, trackIso_;
  edm::InputTag genParticles_;

  bool bothMuons_;

  double etamax_, ptmin_, massMin_, massMax_, isoMax_;

  reco::CandidateBaseRef globalMuonCandRef_, trackMuonCandRef_, standAloneMuonCandRef_;
  OverlapChecker overlap_;

  // general histograms
  TH1D *h_trackProbe_eta, *h_trackProbe_pt, *h_staProbe_eta, *h_staProbe_pt, *h_ProbeOk_eta, *h_ProbeOk_pt;

  // global counters
  int nGlobalMuonsMatched_passed;    // total number of global muons MC matched and passing cuts (and triggered)
  int nGlobalMuonsMatched_passedIso;    // total number of global muons MC matched and passing cuts including Iso
  int n2GlobalMuonsMatched_passedIso;    // total number of Z->2 global muons MC matched and passing cuts including Iso
  int nStaMuonsMatched_passedIso;       // total number of sta only muons MC matched and passing cuts including Iso
  int nTracksMuonsMatched_passedIso;    // total number of tracks only muons MC matched and passing cuts including Iso
  int n2GlobalMuonsMatched_passedIso2Trg;    // total number of Z->2 global muons MC matched and passing cuts including Iso and both triggered
  int nMu0onlyTriggered;               // n. of events zMuMu with mu0 only triggered
  int nMu1onlyTriggered;               // n. of events zMuMu with mu1 only triggered

  int nZMuMu_matched;               // n. of events zMuMu macthed
  int nZMuSta_matched;              // n  of events zMuSta macthed 
  int nZMuTrk_matched;              // n. of events zMuTrk mathced
};

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <iostream>
#include <iterator>
#include <cmath>
using namespace std;
using namespace reco;
using namespace edm;


typedef edm::ValueMap<float> IsolationCollection;

ZMuMu_MCanalyzer::ZMuMu_MCanalyzer(const ParameterSet& pset) : 
  zMuMu_(pset.getParameter<InputTag>("zMuMu")), 
  zMuMuMatchMap_(pset.getParameter<InputTag>("zMuMuMatchMap")), 
  zMuStandAlone_(pset.getParameter<InputTag>("zMuStandAlone")), 
  zMuStandAloneMatchMap_(pset.getParameter<InputTag>("zMuStandAloneMatchMap")), 
  zMuTrack_(pset.getParameter<InputTag>("zMuTrack")), 
  zMuTrackMatchMap_(pset.getParameter<InputTag>("zMuTrackMatchMap")), 
  muons_(pset.getParameter<InputTag>("muons")), 
  tracks_(pset.getParameter<InputTag>("tracks")), 
  genParticles_(pset.getParameter<InputTag>( "genParticles" ) ),

  bothMuons_(pset.getParameter<bool>("bothMuons")), 

  etamax_(pset.getUntrackedParameter<double>("etamax")),  
  ptmin_(pset.getUntrackedParameter<double>("ptmin")), 
  massMin_(pset.getUntrackedParameter<double>("zMassMin")), 
  massMax_(pset.getUntrackedParameter<double>("zMassMax")), 
  isoMax_(pset.getUntrackedParameter<double>("isomax")) { 
  Service<TFileService> fs;

  // binning of entries array (at moment defined by hand and not in cfg file)
  double etaRange[8] = {-2.5,-2.,-1.2,-0.8,0.8,1.2,2.,2.5};
  double ptRange[4] = {20.,40.,60.,100.};

  // general histograms
  h_trackProbe_eta = fs->make<TH1D>("trackProbeEta","Eta of tracks",7,etaRange);
  h_trackProbe_pt = fs->make<TH1D>("trackProbePt","Pt of tracks",3,ptRange);
  h_staProbe_eta = fs->make<TH1D>("standAloneProbeEta","Eta of standAlone",7,etaRange);
  h_staProbe_pt = fs->make<TH1D>("standAloneProbePt","Pt of standAlone",3,ptRange);
  h_ProbeOk_eta = fs->make<TH1D>("probeOkEta","Eta of probe Ok",7,etaRange);
  h_ProbeOk_pt = fs->make<TH1D>("probeOkPt","Pt of probe ok",3,ptRange);

  // clear global counters
  nGlobalMuonsMatched_passed = 0;
  nGlobalMuonsMatched_passedIso = 0;
  n2GlobalMuonsMatched_passedIso = 0;
  nStaMuonsMatched_passedIso = 0;
  nTracksMuonsMatched_passedIso = 0;
  n2GlobalMuonsMatched_passedIso2Trg = 0;
  nMu0onlyTriggered = 0;
  nMu1onlyTriggered = 0;
  nZMuMu_matched = 0;
  nZMuSta_matched = 0;
  nZMuTrk_matched = 0;
}

void ZMuMu_MCanalyzer::analyze(const Event& event, const EventSetup& setup) {
  Handle<CandidateView> zMuMu;  
  Handle<GenParticleMatch> zMuMuMatchMap; //Map of Z made by Mu global + Mu global 
  Handle<CandidateView> zMuStandAlone;  
  Handle<GenParticleMatch> zMuStandAloneMatchMap; //Map of Z made by Mu + StandAlone
  Handle<CandidateView> zMuTrack;  
  Handle<GenParticleMatch> zMuTrackMatchMap; //Map of Z made by Mu + Track
  Handle<CandidateView> muons; //Collection of Muons
  Handle<CandidateView> tracks; //Collection of Tracks

  Handle<GenParticleCollection> genParticles;  // Collection of Generatd Particles
  
  event.getByLabel(zMuMu_, zMuMu); 
  event.getByLabel(zMuStandAlone_, zMuStandAlone); 
  event.getByLabel(zMuTrack_, zMuTrack); 
  event.getByLabel(genParticles_, genParticles);
  event.getByLabel(muons_, muons); 
  event.getByLabel(tracks_, tracks); 

  /*
  cout << "*********  zMuMu         size : " << zMuMu->size() << endl;
  cout << "*********  zMuStandAlone size : " << zMuStandAlone->size() << endl;
  cout << "*********  zMuTrack      size : " << zMuTrack->size() << endl;
  cout << "*********  muons         size : " << muons->size() << endl; 	    
  cout << "*********  tracks        size : " << tracks->size() << endl;
  */
  //      std::cout<<"Run-> "<<event.id().run()<<std::endl;
  //      std::cout<<"Event-> "<<event.id().event()<<std::endl; 


  bool zMuMu_found = false;

  // loop on ZMuMu
  if (zMuMu->size() > 0 ) {
    event.getByLabel(zMuMuMatchMap_, zMuMuMatchMap); 
    for(size_t i = 0; i < zMuMu->size(); ++i) { //loop on candidates
      const Candidate & zMuMuCand = (*zMuMu)[i]; //the candidate
      CandidateBaseRef zMuMuCandRef = zMuMu->refAt(i);

      const Candidate * lep0 = zMuMuCand.daughter( 0 );
      const Candidate * lep1 = zMuMuCand.daughter( 1 );
      const pat::Muon & muonDau0 = dynamic_cast<const pat::Muon &>(*lep0->masterClone());
      double trkiso0 = muonDau0.trackIso();
      const pat::Muon & muonDau1 = dynamic_cast<const pat::Muon &>(*lep1->masterClone());
      double trkiso1 = muonDau1.trackIso();
      double pt0 = zMuMuCand.daughter(0)->pt();
      double pt1 = zMuMuCand.daughter(1)->pt();
      double eta0 = zMuMuCand.daughter(0)->eta();
      double eta1 = zMuMuCand.daughter(1)->eta();
      double mass = zMuMuCand.mass();

      // HLT match
      const std::vector<pat::TriggerPrimitive> & trig0 =muonDau0.triggerMatches();//vector of triggerPrimitive
      const std::vector<pat::TriggerPrimitive> & trig1 =muonDau1.triggerMatches();
      bool trig0found = false;
      bool trig1found = false;
      for (unsigned int j=0; j<trig0.size();j++) {
	if (trig0[j].filterName()=="hltSingleMuNoIsoL3PreFiltered") trig0found = true;
      }
      for (unsigned int j=0; j<trig1.size();j++) {
	if (trig1[j].filterName()=="hltSingleMuNoIsoL3PreFiltered") trig1found = true;
      }
      
      GenParticleRef zMuMuMatch = (*zMuMuMatchMap)[zMuMuCandRef];
      if(zMuMuMatch.isNonnull()) {  // ZMuMu matched
	zMuMu_found = true;
	nZMuMu_matched++;	
	if (pt0>ptmin_ && pt1>ptmin_ && abs(eta0)<etamax_ && abs(eta1) <etamax_ && mass >massMin_ && mass < massMax_ && (trig0found || trig1found)) { // kinematic and trigger cuts passed
	  nGlobalMuonsMatched_passed++; // first global Muon passed kine cuts 
	  nGlobalMuonsMatched_passed++; // second global muon passsed kine cuts
	  if (trkiso0<isoMax_) nGlobalMuonsMatched_passedIso++;       // first global muon passed the iso cut
	  if (trkiso1<isoMax_) nGlobalMuonsMatched_passedIso++;       // second global muon passed the iso cut
	  if (trkiso0<isoMax_ && trkiso1<isoMax_) {
	    n2GlobalMuonsMatched_passedIso++;  // both muons passed iso cut
	    if (trig0found && trig1found) n2GlobalMuonsMatched_passedIso2Trg++;  // both muons have HLT
	    if (trig0found && !trig1found) nMu0onlyTriggered++;
	    if (trig1found && !trig0found) nMu1onlyTriggered++;
	    // histograms vs eta and pt
	    if (trig1found) {         // check efficiency of muon0 not imposing the trigger on it 
	      h_trackProbe_eta->Fill(eta0);
	      h_trackProbe_pt->Fill(pt0);
	      h_staProbe_eta->Fill(eta0);
	      h_staProbe_pt->Fill(pt0);
	      h_ProbeOk_eta->Fill(eta0);
	      h_ProbeOk_pt->Fill(pt0);
	    }
	    if (trig0found) {         // check efficiency of muon1 not imposing the trigger on it 
	      h_trackProbe_eta->Fill(eta1);
	      h_staProbe_eta->Fill(eta1);
	      h_trackProbe_pt->Fill(pt1);
	      h_staProbe_pt->Fill(pt1);
	      h_ProbeOk_eta->Fill(eta1);
	      h_ProbeOk_pt->Fill(pt1);
	    }
	  }
	}
      } // end MC match

    }  // end loop on ZMuMu cand
  }    // end if ZMuMu size > 0

  // loop on ZMuSta
  bool zMuSta_found = false;
  if (!zMuMu_found && zMuStandAlone->size() > 0 ) {
    event.getByLabel(zMuStandAloneMatchMap_, zMuStandAloneMatchMap); 
    for(size_t i = 0; i < zMuStandAlone->size(); ++i) { //loop on candidates
      const Candidate & zMuStandAloneCand = (*zMuStandAlone)[i]; //the candidate
      CandidateBaseRef zMuStandAloneCandRef = zMuStandAlone->refAt(i);
      GenParticleRef zMuStandAloneMatch = (*zMuStandAloneMatchMap)[zMuStandAloneCandRef];

      const Candidate * lep0 = zMuStandAloneCand.daughter( 0 );
      const Candidate * lep1 = zMuStandAloneCand.daughter( 1 );
      const pat::Muon & muonDau0 = dynamic_cast<const pat::Muon &>(*lep0->masterClone());
      double trkiso0 = muonDau0.trackIso();
      const pat::Muon & muonDau1 = dynamic_cast<const pat::Muon &>(*lep1->masterClone());
      double trkiso1 = muonDau1.trackIso();
      double pt0 = zMuStandAloneCand.daughter(0)->pt();
      double pt1 = zMuStandAloneCand.daughter(1)->pt();
      double eta0 = zMuStandAloneCand.daughter(0)->eta();
      double eta1 = zMuStandAloneCand.daughter(1)->eta();
      double mass = zMuStandAloneCand.mass();

      // HLT match (check just dau0 the global)
      const std::vector<pat::TriggerPrimitive> & trig0 =muonDau0.triggerMatches();//vector of triggerPrimitive
      bool trig0found = false;
      for (unsigned int j=0; j<trig0.size();j++) {
	if (trig0[j].filterName()=="hltSingleMuNoIsoL3PreFiltered") trig0found = true;
      }
      
      if(zMuStandAloneMatch.isNonnull()) {  // ZMuStandAlone matched
	zMuSta_found = true;
	nZMuSta_matched++;	
	if (pt0>ptmin_ && pt1>ptmin_ && abs(eta0)<etamax_ && abs(eta1) <etamax_ && mass >massMin_ && 
	    mass < massMax_ && trkiso0<isoMax_ && trkiso1 < isoMax_ && trig0found) { // all cuts and trigger passed
	  nStaMuonsMatched_passedIso++;
	  // histograms vs eta and pt
	  h_staProbe_eta->Fill(eta1);
	  h_staProbe_pt->Fill(pt1);
	}
      } // end MC match
    }  // end loop on ZMuStandAlone cand
  }    // end if ZMuStandAlone size > 0


  // loop on ZMuTrack
  bool zMuTrack_found = false;
  if (!zMuMu_found && !zMuSta_found && zMuTrack->size() > 0 ) {
    event.getByLabel(zMuTrackMatchMap_, zMuTrackMatchMap); 
    for(size_t i = 0; i < zMuTrack->size(); ++i) { //loop on candidates
      const Candidate & zMuTrackCand = (*zMuTrack)[i]; //the candidate
      CandidateBaseRef zMuTrackCandRef = zMuTrack->refAt(i);
      const Candidate * lep0 = zMuTrackCand.daughter( 0 );
      const Candidate * lep1 = zMuTrackCand.daughter( 1 );
      const pat::Muon & muonDau0 = dynamic_cast<const pat::Muon &>(*lep0->masterClone());
      double trkiso0 = muonDau0.trackIso();
      const pat::GenericParticle & trackDau1 = dynamic_cast<const pat::GenericParticle &>(*lep1->masterClone());
      double trkiso1 = trackDau1.trackIso();
      double pt0 = zMuTrackCand.daughter(0)->pt();
      double pt1 = zMuTrackCand.daughter(1)->pt();
      double eta0 = zMuTrackCand.daughter(0)->eta();
      double eta1 = zMuTrackCand.daughter(1)->eta();
      double mass = zMuTrackCand.mass();

      // HLT match (check just dau0 the global)
      const std::vector<pat::TriggerPrimitive> & trig0 =muonDau0.triggerMatches();//vector of triggerPrimitive
      bool trig0found = false;
      for (unsigned int j=0; j<trig0.size();j++) {
	if (trig0[j].filterName()=="hltSingleMuNoIsoL3PreFiltered") trig0found = true;
      }

      GenParticleRef zMuTrackMatch = (*zMuTrackMatchMap)[zMuTrackCandRef];
      if(zMuTrackMatch.isNonnull()) {  // ZMuTrack matched
	zMuTrack_found = true;
	nZMuTrk_matched++;
	if (pt0>ptmin_ && pt1>ptmin_ && abs(eta0)<etamax_ && abs(eta1) <etamax_ && mass >massMin_ && 
	    mass < massMax_ && trkiso0<isoMax_ && trkiso1 < isoMax_ && trig0found) { // all cuts and trigger passed
	  nTracksMuonsMatched_passedIso++;
	  // histograms vs eta and pt
	  h_trackProbe_eta->Fill(eta1);
	  h_trackProbe_pt->Fill(pt1);
	}
      }  // end MC match
    }  // end loop on ZMuTrack cand
  }    // end if ZMuTrack size > 0

}       // end analyze

bool ZMuMu_MCanalyzer::check_ifZmumu(const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2)
{
  int partId0 = dauGen0->pdgId();
  int partId1 = dauGen1->pdgId();
  int partId2 = dauGen2->pdgId();
  bool muplusFound=false;
  bool muminusFound=false;
  bool ZFound=false;
  if (partId0==13 || partId1==13 || partId2==13) muminusFound=true;
  if (partId0==-13 || partId1==-13 || partId2==-13) muplusFound=true;
  if (partId0==23 || partId1==23 || partId2==23) ZFound=true;
  return muplusFound*muminusFound*ZFound;   
}
 
float ZMuMu_MCanalyzer::getParticlePt(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2)
{
  int partId0 = dauGen0->pdgId();
  int partId1 = dauGen1->pdgId();
  int partId2 = dauGen2->pdgId();
  float ptpart=0.;
  if (partId0 == ipart) {
    for(size_t k = 0; k < dauGen0->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen0->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	ptpart = dauMuGen->pt();
      }
    }
  }
  if (partId1 == ipart) {
    for(size_t k = 0; k < dauGen1->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen1->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	ptpart = dauMuGen->pt();
      }
    }
  }
  if (partId2 == ipart) {
    for(size_t k = 0; k < dauGen2->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen2->daughter(k);
      if(abs(dauMuGen->pdgId()) == ipart && dauMuGen->status() ==1) {
	ptpart = dauMuGen->pt();
      }
    }
  }
  return ptpart;
}
 
float ZMuMu_MCanalyzer::getParticleEta(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2)
{
  int partId0 = dauGen0->pdgId();
  int partId1 = dauGen1->pdgId();
  int partId2 = dauGen2->pdgId();
  float etapart=0.;
  if (partId0 == ipart) {
    for(size_t k = 0; k < dauGen0->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen0->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	etapart = dauMuGen->eta();
      }
    }
  }
  if (partId1 == ipart) {
    for(size_t k = 0; k < dauGen1->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen1->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	etapart = dauMuGen->eta();
      }
    }
  }
  if (partId2 == ipart) {
    for(size_t k = 0; k < dauGen2->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen2->daughter(k);
      if(abs(dauMuGen->pdgId()) == ipart && dauMuGen->status() ==1) {
	etapart = dauMuGen->eta();
      }
    }
  }
  return etapart;
}

float ZMuMu_MCanalyzer::getParticlePhi(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2)
{
  int partId0 = dauGen0->pdgId();
  int partId1 = dauGen1->pdgId();
  int partId2 = dauGen2->pdgId();
  float phipart=0.;
  if (partId0 == ipart) {
    for(size_t k = 0; k < dauGen0->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen0->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	phipart = dauMuGen->phi();
      }
    }
  }
  if (partId1 == ipart) {
    for(size_t k = 0; k < dauGen1->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen1->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	phipart = dauMuGen->phi();
      }
    }
  }
  if (partId2 == ipart) {
    for(size_t k = 0; k < dauGen2->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen2->daughter(k);
      if(abs(dauMuGen->pdgId()) == ipart && dauMuGen->status() ==1) {
	phipart = dauMuGen->phi();
      }
    }
  }
  return phipart;
}

Particle::LorentzVector ZMuMu_MCanalyzer::getParticleP4(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2)
{
  int partId0 = dauGen0->pdgId();
  int partId1 = dauGen1->pdgId();
  int partId2 = dauGen2->pdgId();
  Particle::LorentzVector p4part(0.,0.,0.,0.);
  if (partId0 == ipart) {
    for(size_t k = 0; k < dauGen0->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen0->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	p4part = dauMuGen->p4();
      }
    }
  }
  if (partId1 == ipart) {
    for(size_t k = 0; k < dauGen1->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen1->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	p4part = dauMuGen->p4();
      }
    }
  }
  if (partId2 == ipart) {
    for(size_t k = 0; k < dauGen2->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen2->daughter(k);
      if(abs(dauMuGen->pdgId()) == ipart && dauMuGen->status() ==1) {
	p4part = dauMuGen->p4();
      }
    }
  }
  return p4part;
}
 


void ZMuMu_MCanalyzer::endJob() {
  
  
  double eff_Iso = double(nGlobalMuonsMatched_passedIso)/nGlobalMuonsMatched_passed;
  double err_effIso = sqrt(eff_Iso*(1-eff_Iso)/nGlobalMuonsMatched_passed);
 
  double n1_afterIso = 2*n2GlobalMuonsMatched_passedIso2Trg+nMu0onlyTriggered+nMu1onlyTriggered+nTracksMuonsMatched_passedIso;
  double n2_afterIso = 2*n2GlobalMuonsMatched_passedIso2Trg+nMu0onlyTriggered+nMu1onlyTriggered+nStaMuonsMatched_passedIso;
  double effSta_afterIso = (2*n2GlobalMuonsMatched_passedIso2Trg+nMu0onlyTriggered+nMu1onlyTriggered)/n1_afterIso;
  double effTrk_afterIso = (2*n2GlobalMuonsMatched_passedIso2Trg+nMu0onlyTriggered+nMu1onlyTriggered)/n2_afterIso;
  double err_effsta_afterIso = sqrt(effSta_afterIso*(1-effSta_afterIso)/n1_afterIso);
  double err_efftrk_afterIso = sqrt(effTrk_afterIso*(1-effTrk_afterIso)/n2_afterIso);
 

  cout << "------------------------------------  Counters  --------------------------------" << endl;

  cout << "number of events zMuMu matched " << nZMuMu_matched << endl;
  cout << "number of events zMuSta matched " << nZMuSta_matched << endl;
  cout << "number of events zMuTk matched " << nZMuTrk_matched << endl;
  cout << "number of events zMuMu with mu0 only triggered " << nMu0onlyTriggered << endl;
  cout << "number of events zMuMu with mu1 only triggered " << nMu1onlyTriggered << endl;
  cout << "=========================================" << endl;
  cout << "n. of global muons MC matched and passing cuts:           " << nGlobalMuonsMatched_passed << endl;
  cout << "n. of global muons MC matched and passing also Iso cut:       " << nGlobalMuonsMatched_passedIso << endl;
  cout << "n. of Z -> 2 global muons MC matched and passing ALL cuts:    " << n2GlobalMuonsMatched_passedIso << endl;
  cout << "n. of ZMuSta MC matched and passing ALL cuts:    " << nStaMuonsMatched_passedIso << endl;
  cout << "n. of ZmuTrck MC matched and passing ALL cuts:   " << nTracksMuonsMatched_passedIso << endl;
  cout << "n. of Z -> 2 global muons MC matched and passing ALL cuts and both triggered: " << n2GlobalMuonsMatched_passedIso2Trg << endl;
  cout << "=================================================================================" << endl;
  cout << "Iso efficiency: " << eff_Iso << " +/- " << err_effIso << endl;
  cout << "eff StandAlone (after Isocut) : " << effSta_afterIso << "+/-" << err_effsta_afterIso << endl;
  cout << "eff Tracker (after Isocut)    : " << effTrk_afterIso << "+/-" << err_efftrk_afterIso << endl;
  
}
  
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ZMuMu_MCanalyzer);
  