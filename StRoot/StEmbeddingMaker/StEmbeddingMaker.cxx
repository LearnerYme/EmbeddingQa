#include "StEmbeddingMaker.h"

#include <TMath.h>

#include <algorithm>
#include <fstream>
#include <vector>
#include <map>

#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoMcTrack.h"
#include "StPicoEvent/StPicoMcVertex.h"
#include "StThreeVectorF.hh"
#include "StLorentzVector.hh"
#include "Stiostream.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "phys_constants.h"

StEmbeddingMaker::StEmbeddingMaker(
	const char* name, 
	StPicoDstMaker* picoMaker,
    const char* outName
) : StMaker(name) {
	mOutputName = outName;
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
}

StEmbeddingMaker::~StEmbeddingMaker() {}

Int_t StEmbeddingMaker::Init() {
  	mFileOut = new TFile(mOutputName, "recreate");

	// event wise
	h1Ev_Vz = new TH1F(
		"h1Ev_Vz", ";v_{z} cm;Counts",
		300, -150, 150
	);
	h2Ev_VxVy = new TH2F(
		"h2Ev_VxVy", ";v_{x} cm;v_{y} cm;Counts",
		500, -2.4, 2.5,
		500, -2.4, 2.5
	);
	h2Ev_VxVz = new TH2F(
		"h2Ev_VxVz", ";v_{x} cm;v_{z} cm;Counts",
		500, -2.4, 2.5,
		600, -150, 150
	);
	h2Ev_VyVz = new TH2F(
		"h2Ev_VyVz", ";v_{y} cm;v_{z} cm;Counts",
		500, -2.4, 2.5,
		600, -150, 150
	);

	// mult distributions
	h1Ev_RefMult = new TH1F(
		"h1Ev_RefMult", ";RefMult;Counts",
		1000, -0.5, 999.5
	);
	h1Ev_RefMult3 = new TH1F(
		"h1Ev_RefMult3", ";RefMult3;Counts",
		1000, -0.5, 999.5
	);
	h1Ev_NumberOfMcTracks = new TH1F(
		"h1Ev_NumberOfMcTracks", ";N;Counts", 
		200, 0, 200
	);
	h1Ev_NumberOfRcTracks = new TH1F(
		"h1Ev_NumberOfRcTracks", ";N;Counts", 
		2400, 0, 2400
	);

	// track wise
	h1MC_GID = new TH1F(
		"h1MC_GID", ";Geant ID;Counts", 
		21, -0.5, 20.5
	);

	// All MC tracks
	h1AM_p = new TH1F(
		"h1AM_p", "All MC tracks;p GeV/c;Counts",
		100, 0, 10
	);
	h1AM_px= new TH1F(
		"h1AM_px", "All MC tracks;p_{x} GeV/c;Counts",
		50, 0, 5
	);
	h1AM_py = new TH1F(
		"h1AM_py", "All MC tracks;p_{y} GeV/c;Counts",
		50, 0, 5
	);
	h1AM_pz = new TH1F(
		"h1AM_pz", "All MC tracks;p_{z} GeV/c;Counts",
		50, 0, 5
	);
	h1AM_pt = new TH1F(
		"h1AM_pt", "All MC tracks;p_{t} GeV/c;Counts",
		50, 0, 5
	);

	h1AM_eta = new TH1F(
		"h1AM_eta", "All MC tracks;#eta;Counts", 
		500, -2.5, 2.5
	);
	h1AM_phi = new TH1F(
		"h1AM_phi", "All MC tracks;#phi;Counts",
		640, -3.2, 3.2
	);
	h1AM_y = new TH1F(
		"h1AM_y", "All MC tracks;#eta;Counts", 
		500, -2.5, 2.5
	);

	h2AM_eta_pt = new TH2F(
		"h2AM_eta_pt", "All MC tracks;#eta;p_{T} GeV/c;Counts",
		500, -2.5, 2.5, 
		500, 0, 5
	);
	h2AM_y_pt = new TH2F(
		"h2AM_y_pt", "All MC tracks;y;p_{T} GeV/c;Counts",
		500, -2.5, 2.5,
		500, 0, 5
	);

	// All RC tracks
	h1AR_p = new TH1F(
		"h1AR_p", "All Rec. tracks;p GeV/c;Counts",
		100, 0, 10
	);
	h1AR_px= new TH1F(
		"h1AR_px", "All Rec. tracks;p_{x} GeV/c;Counts",
		50, 0, 5
	);
	h1AR_py = new TH1F(
		"h1AR_py", "All Rec. tracks;p_{y} GeV/c;Counts",
		50, 0, 5
	);
	h1AR_pz = new TH1F(
		"h1AR_pz", "All Rec. tracks;p_{z} GeV/c;Counts",
		50, 0, 5
	);
	h1AR_pt = new TH1F(
		"h1AR_pt", "All Rec. tracks;p_{t} GeV/c;Counts",
		50, 0, 5
	);

	h1AR_eta = new TH1F(
		"h1AR_eta", "All Rec. tracks;#eta;Counts", 
		500, -2.5, 2.5
	);
	h1AR_phi = new TH1F(
		"h1AR_phi", "All Rec. tracks;#phi;Counts",
		640, -3.2, 3.2
	);

	h1AR_nHitsFit = new TH1F(
		"h1AR_nHitsFit", "All Rec. tracks;nHitsFit;Counts",
		100, -0.5, 99.5
	);
	h1AR_nHitsFitWithCut = new TH1F(
		"h1AR_nHitsFitWithCut", "All Rec. tracks;nHitsFit;Counts",
		100, -0.5, 99.5
	);
	h1AR_nHitsRatio = new TH1F(
		"h1AR_nHitsRatio", "All Rec. tracks;nHitsRatio;Counts",
		100, 0, 1
	);
	h1AR_nHitsDedx = new TH1F(
		"h1AR_nHitsDedx", "All Rec. tracks;nHitsDedx;Counts",
		100, -0.5, 99.5
	);
	h1AR_nHitsDedxWithCut = new TH1F(
		"h1AR_nHitsDedxWithCut", "All Rec. tracks;nHitsDedx;Counts",
		100, -0.5, 99.5
	);

	h1AR_DCA = new TH1F(
		"h1AR_DCA", "All Rec. tracks;DCA [cm];Counts",
		100, 0, 5
	);
	h1AR_DCAz = new TH1F(
		"h1AR_DCAz", "All Rec. tracks;DCAz [cm];Counts",
		100, -5, 5
	);
	h1AR_sDCAxy = new TH1F(
		"h1AR_sDCAxy", "All Rec. tracks;sDCAxy [cm];Counts",
		100, -5, 5
	);

	h2AR_eta_pt = new TH2F(
		"h2AR_eta_pt", "All Rec. tracks;#eta;p_{T} [GeV/c];Counts", 
		500, -2.5, 2.5,
		500, 0, 5
	);

	// Matched MC tracks
	h1MM_p = new TH1F(
		"h1MM_p", "Matched MC tracks;p GeV/c;Counts",
		100, 0, 10
	);
	h1MM_px= new TH1F(
		"h1MM_px", "Matched MC tracks;p_{x} GeV/c;Counts",
		50, 0, 5
	);
	h1MM_py = new TH1F(
		"h1MM_py", "Matched MC tracks;p_{y} GeV/c;Counts",
		50, 0, 5
	);
	h1MM_pz = new TH1F(
		"h1MM_pz", "Matched MC tracks;p_{z} GeV/c;Counts",
		50, 0, 5
	);
	h1MM_pt = new TH1F(
		"h1MM_pt", "Matched MC tracks;p_{t} GeV/c;Counts",
		50, 0, 5
	);

	h1MM_eta = new TH1F(
		"h1MM_eta", "Matched MC tracks;#eta;Counts", 
		500, -2.5, 2.5
	);
	h1MM_phi = new TH1F(
		"h1MM_phi", "Matched MC tracks;#phi;Counts",
		640, -3.2, 3.2
	);
	h1MM_y = new TH1F(
		"h1MM_y", "Matched MC tracks;#eta;Counts", 
		500, -2.5, 2.5
	);

	h2MM_eta_pt = new TH2F(
		"h2MM_eta_pt", "Matched MC tracks;#eta;p_{T} GeV/c;Counts",
		500, -2.5, 2.5, 
		500, 0, 5
	);
	h2MM_y_pt = new TH2F(
		"h2MM_y_pt", "Matched MC tracks;y;p_{T} GeV/c;Counts",
		500, -2.5, 2.5,
		500, 0, 5
	);

	// Matched RC tracks
	h1MR_p = new TH1F(
		"h1MR_p", "Matched Rec. tracks;p GeV/c;Counts",
		100, 0, 10
	);
	h1MR_px= new TH1F(
		"h1MR_px", "Matched Rec. tracks;p_{x} GeV/c;Counts",
		50, 0.0, 5.0
	);
	h1MR_py = new TH1F(
		"h1MR_py", "Matched Rec. tracks;p_{y} GeV/c;Counts",
		50, 0.0, 5.0
	);
	h1MR_pz = new TH1F(
		"h1MR_pz", "Matched Rec. tracks;p_{z} GeV/c;Counts",
		50, 0.0, 5.0
	);
	h1MR_pt = new TH1F(
		"h1MR_pt", "Matched Rec. tracks;p_{t} GeV/c;Counts",
		50, 0.0, 5.0
	);

	h1MR_eta = new TH1F(
		"h1MR_eta", "Matched Rec. tracks;#eta;Counts", 
		500, -2.5, 2.5
	);
	h1MR_phi = new TH1F(
		"h1MR_phi", "Matched Rec. tracks;#phi;Counts",
		640, -3.2, 3.2
	);
	h1MR_y = new TH1F(
		"h1MR_y", "Matched Rec. tracks;#eta;Counts", 
		500, -2.5, 2.5
	);

	h1MR_nHitsFit = new TH1F(
		"h1MR_nHitsFit", "Matched Rec. tracks;nHitsFit;Counts",
		100, -0.5, 99.5
	);
	h1MR_nHitsFitWithCut = new TH1F(
		"h1MR_nHitsFitWithCut", "Matched Rec. tracks;nHitsFit;Counts",
		100, -0.5, 99.5
	);
	h1MR_nHitsRatio = new TH1F(
		"h1MR_nHitsRatio", "Matched Rec. tracks;nHitsRatio;Counts",
		100, 0, 1
	);
	h1MR_nHitsDedx = new TH1F(
		"h1MR_nHitsDedx", "Matched Rec. tracks;nHitsDedx;Counts",
		100, -0.5, 99.5
	);
	h1MR_nHitsDedxWithCut = new TH1F(
		"h1MR_nHitsDedxWithCut", "Matched Rec. tracks;nHitsDedx;Counts",
		100, -0.5, 99.5
	);

	h1MR_DCA = new TH1F(
		"h1MR_DCA", "Matched Rec. tracks;DCA [cm];Counts",
		100, 0, 5
	);
	h1MR_DCAz = new TH1F(
		"h1MR_DCAz", "Matched Rec. tracks;DCAz [cm];Counts",
		100, -5, 5
	);
	h1MR_sDCAxy = new TH1F(
		"h1MR_sDCAxy", "Matched Rec. tracks;sDCAxy [cm];Counts",
		100, -5, 5
	);

	h2MR_eta_pt = new TH2F(
		"h2MR_eta_pt", "Matched Rec. tracks;#eta;p_{T} [GeV/c];Counts", 
		500, -2.5, 2.5,
		500, 0, 5
	);
	h2MR_y_pt = new TH2F(
		"h2MR_y_pt", "Matched Rec. tracks;y;p_{T} [GeV/c];Counts", 
		500, -2.5, 2.5,
		500, 0, 5
	);

	// Matched MC Vs Rc
	h2Mcmp_p = new TH2F(
		"h2Mcmp_p", "Matched protons;p [GeV/c] (MC);p [GeV/c] (RC);Counts",
		100, 0, 5, 
		100, 0, 5
	);
	h2Mcmp_pt = new TH2F(
		"h2Mcmp_pt", "Matched protons;p_{T} [GeV/c] (MC);p_{T} [GeV/c] (RC);Counts", 
		100, 0, 5, 
		100, 0, 5
	);
	h2Mcmp_eta = new TH2F(
		"h2Mcmp_eta", "Matched protons;#eta (MC);#eta (RC);Counts", 
		100, -2.5, 2.5,
		100, -2.5, 2.5
	);
	h2Mcmp_phi = new TH2F(
		"h2Mcmp_phi", "Matched protons;#phi (MC);#phi (RC);Counts", 
		100, -3.2, 3.2,
		100, -3.2, 3.2
	);
	h2Mcmp_y = new TH2F(
		"h2Mcmp_y", "Matched protons;y (MC);y (RC);Counts", 
		100, -2.5, 2.5,
		100, -2.5, 2.5
	);


	nEvents = 0;

	return kStOK;
}

//---------------------------------------------------------
Int_t StEmbeddingMaker::Finish() {
	std::cout << "[LOG] Number of events: " << nEvents << ".\n";
	mFileOut->cd();
	mFileOut->Write();
	mFileOut->Close();
	std::cout << "[LOG] This is the end of this job.\n";
	return kStOK;
}

void StEmbeddingMaker::Clear(Option_t* opt) {}

//---------------------------------------------------------------
Int_t StEmbeddingMaker::Make() {
	if (!mPicoDstMaker) {
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();
	if (!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	if (!mPicoDst) {
		return kStOK;
	}

	// Load event
	event = (StPicoEvent*)mPicoDst->event();
	if (!event) {
		cerr << "Error opening picoDst Event, skip!" << endl;
		return kStOK;
	}
	nEvents += 1;

	TVector3 pVtx = event->primaryVertex();
	Double_t vx = pVtx.X();
	Double_t vy = pVtx.Y();
	Double_t vz = pVtx.Z();
	Double_t vr = sqrt(vx * vx + vy * vy);
	if (vr > 2 || fabs(vz) > 150) {
		return kStOK;
	}

	h1Ev_Vz->Fill(vz);
	h2Ev_VxVy->Fill(vx, vy);
	h2Ev_VxVz->Fill(vx, vz);
	h2Ev_VyVz->Fill(vy, vz);

	h1Ev_RefMult->Fill(event->refMult());
	h1Ev_RefMult3->Fill(event->refMult3());
	Int_t numberOfMcTracks = mPicoDst->numberOfMcTracks();
	Int_t numberOfRcTracks = mPicoDst->numberOfTracks();
	h1Ev_NumberOfMcTracks->Fill(numberOfMcTracks);
	h1Ev_NumberOfRcTracks->Fill(numberOfRcTracks);

	Int_t numberOfMcVertices = mPicoDst->numberOfMcVertices();
	if (!numberOfMcVertices || !numberOfMcTracks) { 
		// this event has no MC information, skip it
		return kStOK;
	}

	// Reconstructed track loop
	// to construct the map of RcTrack ID -> McTrack ID
	TVector3 pMom;
	StThreeVector<Float_t> pMom3;
	StLorentzVector<Float_t> pMom4;
	std::map<Int_t, Int_t> mMc2Rc;
	Float_t MP = 0.938272;
	Float_t EP;
	Int_t nHitsFit;
	Double_t nHitsRatio;
	Int_t nHitsDedx;
	Double_t dca, dcaz, sdcaxy;
	Double_t mBField = event->bField();

	for (Int_t iRcTrk=0; iRcTrk<numberOfRcTracks; iRcTrk++) {
		rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk);
		if (!rcTrack) {
			continue;
		}
		if (!rcTrack->isPrimary()) {
			continue;
		}

		Int_t idTruth = rcTrack->idTruth(); // index of corresponding MC track
		if (idTruth <= 0 || idTruth > 10000) {
			continue;
		}

		// note that, here idTruth - 1 will be the quantity iMcTrk in later codes
		// not just idTruth. VERY IMPORTANT!
		mMc2Rc.insert(std::pair<Int_t, Int_t>(idTruth - 1, iRcTrk)); 

		// Fill All Rc track histograms

		nHitsFit = rcTrack->nHitsFit();
		nHitsRatio = nHitsFit*1.0 / rcTrack->nHitsPoss();
		nHitsDedx = rcTrack->nHitsDedx();
		h1AR_nHitsFit->Fill(nHitsFit);
		h1AR_nHitsDedx->Fill(nHitsDedx);
		h1AR_nHitsRatio->Fill(nHitsRatio);
		if (nHitsRatio > 0.52) {
			h1AR_nHitsFitWithCut->Fill(nHitsFit);
			h1AR_nHitsDedxWithCut->Fill(nHitsDedx);
		}

		dca = rcTrack->gDCA(vx, vy, vz);
		dcaz = rcTrack->gDCAz(vz);
		sdcaxy = rcTrack->helix(mBField).geometricSignedDistance(vx, vy);

		h1AR_DCA->Fill(dca);
		h1AR_DCAz->Fill(dcaz);
		h1AR_sDCAxy->Fill(sdcaxy);

		// do track quality cut before get momentum
		if (dca > 1.0 || nHitsFit < 20 || nHitsDedx < 5 || nHitsRatio < 0.52) {
			continue;
		}

		pMom = rcTrack->pMom();
		// do not need them for rapidity calculation for all rc tracks
		// pMom3 = StThreeVector<Float_t>(pMom.X(), pMom.Y(), pMom.Z());
		// // calculate rapidity
		// EP = sqrt(pMom3.mag2() + MP*MP);
		// rcMom4 = StLorentzVector<Float_t>(pMom3, EP);

		h1AR_p->Fill(pMom.Mag());
		h1AR_px->Fill(pMom.X());
		h1AR_py->Fill(pMom.Y());
		h1AR_pz->Fill(pMom.Z());
		h1AR_pt->Fill(pMom.Perp());
		h1AR_eta->Fill(pMom.PseudoRapidity());
		h1AR_phi->Fill(pMom.Phi());

		h2AR_eta_pt->Fill(pMom.PseudoRapidity(), pMom.Perp());
	}

	// MC track loop
	// to record the MC track and its RC track
	bool flag = false;
	TVector3 rcMom; // this MC block, will use rcMom for RC tarcks
	StThreeVector<Float_t> rcMom3;
	StLorentzVector<Float_t> rcMom4;
	for (Int_t iMcTrk=0; iMcTrk<numberOfMcTracks; iMcTrk++){
		mcTrack = (StPicoMcTrack*)mPicoDst->mcTrack(iMcTrk);
		if (!is_McTrack_from_PV()) {
			continue;
		}

		// find the reconstruct track
		// Indeed, mcTrack->id() is iMcTrk + 1;
		// cannot just use iRcTrk = mMc2Rc.find(iMcTrk)->second to get the RC track ID
		pair<map<Int_t, Int_t>::iterator, map<Int_t, Int_t>::iterator> ret;
		ret = mMc2Rc.equal_range(iMcTrk);
		map<Int_t, Int_t>::iterator iter;
		Int_t count = 0;
		Int_t iRcTrk = -1;
		Int_t iRcTrk_best = -1;
		Int_t qaTruth_best = -1;
		for (iter=ret.first; iter!=ret.second; iter++, count++) { // loop over the possible reconstructed tracks and find the best one with greatest qaTruth
			iRcTrk = iter->second;
			rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk);
			if (!rcTrack) {
				continue;
			}
			if (!rcTrack->isPrimary()){
				continue;
			}
			if (rcTrack->qaTruth() > qaTruth_best){
				qaTruth_best = rcTrack->qaTruth();
				iRcTrk_best = iRcTrk;
			}
		} 
		if (count > 0) { // indicates that rc track was found, and set the best one
			rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk_best);
			if (!rcTrack || !rcTrack->isPrimary()){
				flag = false;
			} else {
				flag = true;
			}
		} else {
			rcTrack = 0;
			flag = false;
		}

		// Fill Geant ID histogram
		h1MC_GID->Fill(mcTrack->geantId());
		if (!is_target_particle()) {
			// only fill Geant ID histogram itself when the track is not proton
			continue;
		}
		if (mcTrack->idVtxStart() != 1) {
			// make sure this track is from PV
			continue;
		}
		
		// Fill all MC histograms
		// in this MC loop, pMom will be used for MC,
		// and rc tracks will have their rcMom
		pMom = mcTrack->p();
		h1AM_p->Fill(pMom.Mag());
		h1AM_px->Fill(pMom.X());
		h1AM_py->Fill(pMom.Y());
		h1AM_pz->Fill(pMom.Z());
		h1AM_pt->Fill(pMom.Perp());
		h1AM_eta->Fill(mcTrack->eta());
		h1AM_phi->Fill(pMom.Phi());
		h1AM_y->Fill(mcTrack->rapidity());
		h2AM_eta_pt->Fill(mcTrack->eta(), pMom.Perp());
		h2AM_y_pt->Fill(mcTrack->rapidity(), pMom.Perp());

		// Get RC information
		if (flag) {
			rcMom = rcTrack->pMom();
			rcMom3 = StThreeVector<Float_t>(rcMom.X(), rcMom.Y(), rcMom.Z());
			// calculate rapidity
			EP = sqrt(rcMom3.mag2() + MP*MP);
			rcMom4 = StLorentzVector<Float_t>(rcMom3, EP);
		}

		if (flag) {
			// fill matched MC tracks

			h1MM_p->Fill(pMom.Mag());
			h1MM_px->Fill(pMom.X());
			h1MM_py->Fill(pMom.Y());
			h1MM_pz->Fill(pMom.Z());
			h1MM_pt->Fill(pMom.Perp());
			h1MM_eta->Fill(mcTrack->eta());
			h1MM_phi->Fill(pMom.Phi());
			h1MM_y->Fill(mcTrack->rapidity());
			h2MM_eta_pt->Fill(mcTrack->eta(), pMom.Perp());
			h2MM_y_pt->Fill(mcTrack->rapidity(), pMom.Perp());

			// fill macthed RC tracks

			h1MR_p->Fill(rcMom.Mag());
			h1MR_px->Fill(rcMom.X());
			h1MR_py->Fill(rcMom.Y());
			h1MR_pz->Fill(rcMom.Z());
			h1MR_pt->Fill(rcMom.Perp());
			h1MR_eta->Fill(rcMom.PseudoRapidity());
			h1MR_phi->Fill(rcMom.Phi());
			h1MR_y->Fill(rcMom4.rapidity());

			nHitsFit = rcTrack->nHitsFit();
			nHitsRatio = nHitsFit*1.0 / rcTrack->nHitsPoss();
			nHitsDedx = rcTrack->nHitsDedx();
			h1MR_nHitsFit->Fill(nHitsFit);
			h1MR_nHitsDedx->Fill(nHitsDedx);
			h1MR_nHitsRatio->Fill(nHitsRatio);
			if (nHitsRatio > 0.52) {
				h1MR_nHitsFitWithCut->Fill(nHitsFit);
				h1MR_nHitsDedxWithCut->Fill(nHitsDedx);
			}

			dca = rcTrack->gDCA(vx, vy, vz);
			dcaz = rcTrack->gDCAz(vz);
			sdcaxy = rcTrack->helix(mBField).geometricSignedDistance(vx, vy);

			h1MR_DCA->Fill(dca);
			h1MR_DCAz->Fill(dcaz);
			h1MR_sDCAxy->Fill(sdcaxy);

			h2MR_eta_pt->Fill(rcMom.PseudoRapidity(), rcMom.Perp());
			h2MR_y_pt->Fill(rcMom4.rapidity(), rcMom.Perp());

			// matched MC vs RC heatmap
			h2Mcmp_p->Fill(pMom.Mag(), rcMom.Mag());
			h2Mcmp_pt->Fill(pMom.Perp(), rcMom.Perp());
			h2Mcmp_eta->Fill(mcTrack->eta(), rcMom.PseudoRapidity());
			h2Mcmp_phi->Fill(pMom.Phi(), rcMom.Phi());
			h2Mcmp_y->Fill(mcTrack->rapidity(), rcMom4.rapidity());
		}
	}

	return kStOK;
}

bool StEmbeddingMaker::is_target_particle() {
	return (mcTrack->geantId() == targetID);
}

bool StEmbeddingMaker::is_McTrack_from_PV() {
	if (!mcTrack) {
		return false;
	}
	Int_t idMcVx = mcTrack->idVtxStart();
	while (idMcVx != 1) {
		StPicoMcVertex* mcVertex = (StPicoMcVertex*)mPicoDst->mcVertex(idMcVx - 1);
		Int_t idMcTrack = mcVertex->idOfParentTrack();
		if (!idMcTrack) {
			break;
		}
		StPicoMcTrack* mcTrackP = (StPicoMcTrack*)mPicoDst->mcTrack(idMcTrack - 1);
		idMcVx = mcTrackP->idVtxStart();
		if (!idMcVx) {
			break;
		}
	}
	if (idMcVx != 1) {
		return false;
	}
	return true;
}
