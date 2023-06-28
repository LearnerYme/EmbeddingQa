#ifndef _StEmbeddingMaker_head
#define _StEmbeddingMaker_head
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "TVector3.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "TNtuple.h"

#include "cent_util.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoMcTrack;
class StPicoDstMaker;
class StPicoMcVertex;
class TH1F;
class TH2F;
class TProfile;
class TTree;
class TNtuple;

class Corr;

class StEmbeddingMaker : public StMaker {
	public:
		StEmbeddingMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName="tofMatchTree.root");
		virtual ~StEmbeddingMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

		void set_target_ID(Int_t val) {
			targetID = val;
		}
		bool is_target_particle();
		bool is_McTrack_from_PV();


	private:
		StPicoDstMaker *mPicoDstMaker;
		StPicoDst      *mPicoDst;
		StPicoEvent	   *event;
		StPicoTrack    *rcTrack;
		StPicoMcTrack  *mcTrack;

		Corr* cent_corr;

		// target particle ID
		Int_t targetID;

		// event wise 
		TH1F* h1Ev_Vz;
		TH2F* h2Ev_VxVy;
		TH2F* h2Ev_VxVz;
		TH2F* h2Ev_VyVz;

		// multiplicity distributions
		TH1F* h1Ev_RefMult;
		TH1F* h1Ev_RefMult3;
		TH1F* h1Ev_NumberOfMcTracks; // MC tracks per event
		TH1F* h1Ev_NumberOfRcTracks; // RC tracks per event

		// track wise 
		TH1F* h1MC_GID; // MC track Geant ID

		// All MC tracks
		TH1F* h1AM_p;
		TH1F* h1AM_px;
		TH1F* h1AM_py;
		TH1F* h1AM_pz;
		TH1F* h1AM_pt;
		TH1F* h1AM_eta;
		TH1F* h1AM_phi;
		TH1F* h1AM_y;
		TH2F* h2AM_eta_pt;
		TH2F* h2AM_y_pt;

		// All RC tracks
		TH1F* h1AR_p;
		TH1F* h1AR_px;
		TH1F* h1AR_py;
		TH1F* h1AR_pz;
		TH1F* h1AR_pt;
		TH1F* h1AR_eta;
		TH1F* h1AR_phi;
		TH1F* h1AR_nHitsFit;
		TH1F* h1AR_nHitsFitWithCut;
		TH1F* h1AR_nHitsRatio;
		TH1F* h1AR_nHitsDedx;
		TH1F* h1AR_nHitsDedxWithCut;
		TH1F* h1AR_DCA;
		TH1F* h1AR_DCAz;
		TH1F* h1AR_sDCAxy;
		TH2F* h2AR_eta_pt;

		// Matched MC tracks (protons only)
		TH1F* h1MM_p;
		TH1F* h1MM_px;
		TH1F* h1MM_py;
		TH1F* h1MM_pz;
		TH1F* h1MM_pt;
		TH1F* h1MM_eta;
		TH1F* h1MM_phi;
		TH1F* h1MM_y;
		TH2F* h2MM_eta_pt;
		TH2F* h2MM_y_pt;

		// Matched RC tracks (protons only)
		TH1F* h1MR_p;
		TH1F* h1MR_px;
		TH1F* h1MR_py;
		TH1F* h1MR_pz;
		TH1F* h1MR_pt;
		TH1F* h1MR_eta;
		TH1F* h1MR_phi;
		TH1F* h1MR_y; // as proton
		TH1F* h1MR_nHitsFit;
		TH1F* h1MR_nHitsFitWithCut;
		TH1F* h1MR_nHitsRatio;
		TH1F* h1MR_nHitsDedx;
		TH1F* h1MR_nHitsDedxWithCut;
		TH1F* h1MR_DCA;
		TH1F* h1MR_DCAz;
		TH1F* h1MR_sDCAxy;
		TH2F* h2MR_eta_pt;
		TH2F* h2MR_y_pt;

		// Matched MC Vs Rc
		TH2F* h2Mcmp_p;
		TH2F* h2Mcmp_pt;
		TH2F* h2Mcmp_eta; // as proton
		TH2F* h2Mcmp_phi;
		TH2F* h2Mcmp_y; // as proton

		TString mOutputName;
		TFile* mFileOut;

		Int_t nEvents;

		ClassDef(StEmbeddingMaker, 1)
};


ClassImp(StEmbeddingMaker)

#endif
