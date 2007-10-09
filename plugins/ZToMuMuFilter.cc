/* \class ZToMuMuFilter
 *
 * \author Juan Alcaraz, CIEMAT
 *
 */
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class ZToMuMuFilter : public edm::EDFilter {
public:
  ZToMuMuFilter(const edm::ParameterSet &);
private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  edm::InputTag zCands_, muonIsolations_;
  double ptMin_, etaMin_, etaMax_, massMin_, massMax_, isoMax_;
};

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
using namespace edm;
using namespace std;
using namespace reco;

ZToMuMuFilter::ZToMuMuFilter( const ParameterSet & cfg ) :
  zCands_(cfg.getParameter<InputTag>("zCands")),
  muonIsolations_(cfg.getParameter<InputTag>("muonIsolations")),
  ptMin_(cfg.getParameter<double>("ptMin")),
  etaMin_(cfg.getParameter<double>("etaMin")),
  etaMax_(cfg.getParameter<double>("etaMax")),
  massMin_(cfg.getParameter<double>("massMin")),
  massMax_(cfg.getParameter<double>("massMax")),
  isoMax_(cfg.getParameter<double>("isoMax")) {
}

bool ZToMuMuFilter::filter (Event & ev, const EventSetup &) {
  Handle<CandidateCollection> zToMuMu;
  ev.getByLabel(zCands_, zToMuMu);
  Handle<CandDoubleAssociations> muIso;
  ev.getByLabel(muonIsolations_, muIso);
  size_t nZ = zToMuMu->size();
  if ( nZ == 0) return false;
  for( size_t i = 0; i < nZ; ++ i ) {
    const Candidate & zCand = (*zToMuMu)[i];
    double zMass = zCand.mass();
    if ( zMass < massMin_ || zMass > massMax_ ) return false;
    if(zCand.numberOfDaughters()!=2) return false;
    const Candidate * dau0 = zCand.daughter(0);
    const Candidate * dau1 = zCand.daughter(1);
    double pt0 = dau0->pt(), pt1 = dau1->pt();
    if ( pt0 < ptMin_ || pt1 < ptMin_ ) return false;
    double eta0 = dau0->eta(), eta1 = dau1->eta();
    if( eta0 < etaMin_ || eta0 > etaMax_ ) return false;
    if( eta1 < etaMin_ || eta1 > etaMax_ ) return false;
    CandidateRef mu0 = dau0->masterClone().castTo<CandidateRef>();
    CandidateRef mu1 = dau1->masterClone().castTo<CandidateRef>();
    double iso0 = (*muIso)[mu0];
    double iso1 = (*muIso)[mu1];
    if (iso0 > isoMax_) return false;
    if (iso1 > isoMax_) return false;
  }
  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( ZToMuMuFilter );
