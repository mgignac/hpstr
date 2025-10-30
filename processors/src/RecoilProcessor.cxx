// ============================================================================
// File: src/hpstr/processors/src/RecoilProcessor.cxx
// Author: Emrys Peets, SLAC, Stanford
// Truth matching priority:
//   (A) Truth-Track-ID path (uses Track::getID() set from lc_truth_track->id())
//   (B) Cluster MCPartIDs (ID-based per-hit association)
//   (C) Layer-overlap fallback (no distances)
// Denominator: truth-layer findability stays in MCParticle space.
// Match rule: half-the-hits + truth findable at min threshold; front-layer single-hit exception.
// PDG==13 kept.  Added TruthTracks reading + TT→MC mapping by layer overlap.
// ============================================================================

#include "RecoilProcessor.h"

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

#include <EVENT/SimTrackerHit.h>
#include "TrackerHit.h"
// NOTE: Track class include should already be pulled in by RecoilProcessor.h in this codebase.
// If not, include the appropriate header for your Track type.

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClass.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TString.h" // Form

namespace {
    // Global diagnostics
    TH1F* g_match_completeness     = nullptr;
    TH1F* g_match_matchedLayers    = nullptr;
    TH2F* g_comp_vs_purity         = nullptr;
    TH1F* g_fakes_per_event        = nullptr;
    TH1F* g_duplicates_per_event   = nullptr;
    TH1F* g_missedLayers_per_event = nullptr;
    TH1F* g_mixedHits_per_event    = nullptr;

    static std::vector<std::string> list_branches(TTree* t) {
        std::vector<std::string> out; if (!t) return out;
        auto* lb = t->GetListOfBranches(); if (!lb) return out;
        for (int i = 0; i < lb->GetEntries(); ++i) out.emplace_back(lb->At(i)->GetName());
        return out;
    }
    //static std::string resolve_branch(TTree* t, const std::vector<std::string>& candidates) {
    //    for (const auto& c : candidates) if (t && t->GetBranch(c.c_str())) return c;
    //    return {};
    //}
    //
    //
    //
    //
    
   // static bool bind_branch_loose(TTree* t,
   //                           const char* name,
   //                           void* addr,
   //                           TBranch** out_b,
   //                           int debug,
   //                           bool list_once = true) {
   // static bool printed = false;
   // if (!t || !name || !*name) return false;

   // TBranch* br = t->GetBranch(name);
   // if (!br) {
   //     if (debug > 0) {
   //         std::cout << "[RecoilProcessor] Branch '" << name << "' not found.\n";
   //         if (list_once && !printed) {
   //             printed = true;
   //             std::cout << "[RecoilProcessor] Available branches:\n";
   //             if (auto* lb = t->GetListOfBranches()) {
   //                 for (int i = 0; i < lb->GetEntries(); ++i)
   //                     std::cout << "  - " << lb->At(i)->GetName() << "\n";
   //             }
   //         }
   //     }
   //     if (out_b) *out_b = nullptr;
   //     return false;
   // }

   // // Do NOT probe expected type — it can crash without dictionaries.
   // int rc = t->SetBranchAddress(name, addr, out_b);
   // if (rc != 0) {
   //     if (debug > 0) {
   //         std::cout << "[RecoilProcessor] SetBranchAddress('" << name << "') failed (rc=" << rc << ").\n";
   //     }
   //     if (out_b) *out_b = nullptr;
   //     return false;
   // }

   // if (debug > 0) {
   //     const char* cln = "<unknown>";
   //     if (auto* be = dynamic_cast<TBranchElement*>(br)) {
   //         if (be->GetClassName() && be->GetClassName()[0]) cln = be->GetClassName();
   //     }
   //     std::cout << "[RecoilProcessor] Using branch: " << name << " (class=" << cln << ")\n";
   // }
   // return true;
   // }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // Minimal safe binder (WHY: avoid segfault on missing branch or missing dictionary)
    //static bool bind_branch(TTree* t,
    //                        const char* name,
    //                        void* addr,           // pass &ptrVar
    //                        TBranch** out_b,      // pass &bVar
    //                        int debug,
    //                        bool list_once = true) {
    //    static bool printed = false;
    //    if (!t || !name || !*name) return false;
    //    TBranch* br = t->GetBranch(name);
    //    if (!br) {
    //        if (debug > 0) {
    //            std::cout << "[RecoilProcessor] Branch '" << name << "' not found.\n";
    //            if (list_once && !printed) {
    //                std::cout << "[RecoilProcessor] Available branches:\n";
    //                for (const auto& b : list_branches(t)) std::cout << "  - " << b << "\n";
    //                printed = true;
    //            }
    //        }
    //        return false;
    //    }
    //    TClass* cl = nullptr; EDataType dt = kOther_t;
    //    // NOTE: some split/leaflist branches may not report a useful class; we only *skip* if class==nullptr.
    //    const bool has_type = (br->GetExpectedType(cl, dt) == 0) && (cl != nullptr);
    //    if (!has_type) {
    //        if (debug > 0) {
    //            std::cout << "[RecoilProcessor] Branch '" << name
    //                      << "' has no dictionary/type; skipping SetBranchAddress to avoid crash.\n";
    //        }
    //        return false;
    //    }
    //    const int rc = t->SetBranchAddress(name, addr, out_b);
    //    if (rc != 0 && debug > 0) {
    //        std::cout << "[RecoilProcessor] SetBranchAddress('" << name << "') failed (rc=" << rc << ").\n";
    //    }
    //    if (rc == 0 && debug > 0) {
    //        std::cout << "[RecoilProcessor] Using branch: " << name << " (class=" << cl->GetName() << ")\n";
    //    }
    //    return (rc == 0);
    //}
}

RecoilProcessor::RecoilProcessor(const std::string& name, Process& process)
    : Processor(name,process) {}

RecoilProcessor::~RecoilProcessor() {}

void RecoilProcessor::configure(const ParameterSet& parameters) {
    try {
        debug_     = parameters.getInteger("debug", debug_);
        anaName_   = parameters.getString("anaName", anaName_);
        beamE_     = parameters.getDouble("beamE", (beamE_>0 ? beamE_ : 3.74));
        sampleID_  = parameters.getInteger("sampleID", sampleID_);
        isData     = parameters.getInteger("isData", isData);

        // Explicit cluster collection name
        sclusColl_ = parameters.getString("sclusColl", sclusColl_);
        
	// Name of the ROOT truth-track branch written by TrackingProcessor
        truthTracksCollRoot_ = parameters.getString("truthTrackCollRoot", truthTracksCollRoot_);
        
	// NOTE: thrList_ and kNThr_/kNLayersMax_ are fixed in header; we do not modify them here.
    } catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    }
}




//void RecoilProcessor::configure(const ParameterSet& parameters) {
//    try {
//        debug_     = parameters.getInteger("debug", debug_);
//        anaName_   = parameters.getString("anaName", anaName_);
//        beamE_     = parameters.getDouble("beamE", (beamE_>0 ? beamE_ : 3.74));
//        sampleID_  = parameters.getInteger("sampleID", sampleID_);
//        isData     = parameters.getInteger("isData", isData);
//        sclusColl_ = parameters.getString("sclusColl", sclusColl_);
//        // matches your recotuple_cfg.py key:
//        truthTracksCollRoot_ = parameters.getString("truthTrackCollRoot", truthTracksCollRoot_);
//    } catch (std::runtime_error& e) {
//        std::cout << e.what() << std::endl;
//    }
//}

void RecoilProcessor::initialize(TTree* tree) {
    tree_ = tree;

    // ---- SAFE BIND: required inputs
    evth_ = nullptr;  bevth_ = nullptr;  // bcvth_ in case your header uses bEvth_/bevth_; keep consistent
    mcParts_ = nullptr; bmcParts_ = nullptr;
    trks_ = nullptr; btrks_ = nullptr;
    mcTrkrHits_ = nullptr; bmcTrkrHits_ = nullptr;

    tree_->SetBranchAddress("EventHeader",       &evth_,        &bevth_);
    tree_->SetBranchAddress("MCParticle",        &mcParts_,      &bmcParts_);
    tree_->SetBranchAddress("KalmanFullTracks",  &trks_,         &btrks_);
    tree_->SetBranchAddress( "TrackerSimHits",    &mcTrkrHits_,   &bmcTrkrHits_);


	/* ---- BEGIN SAFE TRUTH BRANCH (super-minimal) ---- */
	
    truthTrks_  = nullptr;    // must be nullptr before binding
    bTruthTrks_ = nullptr;

    // New (optional): TruthTracks with LCIO Truth Track IDs
    if (tree_->GetBranch("Truth_KalmanTracks")) {
        tree_->SetBranchAddress("Truth_KalmanTracks", &truthTrks_, &bTruthTrks_);
        if (debug_>0) std::cout << "[RecoilProcessor] Using truth branch: Truth_KalmanTracks\n";
    } else {
        truthTrks_  = nullptr;
        bTruthTrks_ = nullptr;
        if (debug_>0) std::cout << "[RecoilProcessor] Truth_KalmanTracks branch not found; primary Truth-Track-ID path will be skipped.\n";
    }
    
    
    
    //if (!truthTracksCollRoot_.empty()) {
//    if (tree_->GetBranch(truthTracksCollRoot_.c_str())) {
//        int rc = tree_->SetBranchAddress(truthTracksCollRoot_.c_str(),
//                                         &truthTrks_, &bTruthTrks_);
//        if (rc != 0) { truthTrks_ = nullptr; bTruthTrks_ = nullptr; }
//    }
//}

      //  if (!truthTracksCollRoot_.empty()) {
      //      if (tree_->GetBranch(truthTracksCollRoot_.c_str())) {
      //          tree_->SetBranchAddress(truthTracksCollRoot_.c_str(), &truthTrks_, &bTruthTrks_);
     //       if (debug_>0) {
     //           std::cout << "[RecoilProcessor] Using truth-track branch: "
     //                     << truthTracksCollRoot_ << "\n";
    //        }
    //        } else if (debug_>0) {
    //        std::cout << "[RecoilProcessor] Requested truth-track branch '"
    //                  << truthTracksCollRoot_ << "' not found. Available branches:\n";
    //        for (auto& b : list_branches(tree_)) std::cout << "  - " << b << "\n";
    //    }
    //} else if (debug_>0) {
    //    std::cout << "[RecoilProcessor] truthTrackCollRoot not set; TT-space disabled (MC fallback).\n";
    //}




	    //if (!truthTracksCollRoot_.empty()) {
	    //    TBranch* br = tree_ ? tree_->GetBranch(truthTracksCollRoot_.c_str()) : nullptr;
	    //    if (br) {
	    //        const int rc = tree_->SetBranchAddress(truthTracksCollRoot_.c_str(),
	    //    					   &truthTrks_, &bTruthTrks_);
	    //        if (rc != 0) {
	    //    	std::cout << "[RecoilProcessor] SetBranchAddress('" << truthTracksCollRoot_
	    //    		  << "') failed (rc=" << rc << "); leaving unbound.\n";
	    //    	truthTrks_ = nullptr;  // stay safe
	    //    	bTruthTrks_ = nullptr;
	    //        } else if (debug_ > 0) {
	    //    	std::cout << "[RecoilProcessor] Using truth-track branch: "
	    //    		  << truthTracksCollRoot_ << "\n";
	    //        }
	    //    } else if (debug_ > 0) {
	    //        std::cout << "[RecoilProcessor] Requested truth-track branch '" << truthTracksCollRoot_
	    //    	      << "' not found. Available branches:\n";
	    //        auto* lb = tree_->GetListOfBranches();
	    //        if (lb) for (int i=0;i<lb->GetEntries();++i) std::cout << "  - " << lb->At(i)->GetName() << "\n";
	    //    }
	    //} else if (debug_ > 0) {
	    //    std::cout << "[RecoilProcessor] truthTrackCollRoot not set; TT-space disabled (using MC-space).\n";
	    //}
	
	/* ---- END SAFE TRUTH BRANCH ---- */






    // ---- SAFE BIND: truth-track branch (optional)
    //truthTrks_  = nullptr;
    //bTruthTrks_ = nullptr;
    //if (!truthTracksCollRoot_.empty()) {
    //    bind_branch(tree_, truthTracksCollRoot_.c_str(), &truthTrks_, &bTruthTrks_, debug_);
    //} else if (debug_ > 0) {
    //    std::cout << "[RecoilProcessor] truthTrackCollRoot not set; TT-space disabled (using MC-space).\n";
    //}

    // ---- Resolve a cluster collection (unchanged)
    {
        std::vector<std::string> clusterCandidates;
        if (!sclusColl_.empty()) clusterCandidates.push_back(sclusColl_);
        clusterCandidates.push_back("SiClustersOnTrack");
        clusterCandidates.push_back("SiClustersOnTrack_KF");
        clusterCandidates.push_back("SvtClustersOnTrack");
        clusterCandidates.push_back("SvtClustersOnTrack_KF");
        clusterCandidates.push_back("SiClusters");
        clusterCandidates.push_back("SvtClusters");
        auto resolve_branch = [](TTree* t, const std::vector<std::string>& candidates) -> std::string {
            for (const auto& c : candidates) if (t && t->GetBranch(c.c_str())) return c;
            return {};
        };
        const std::string sclusName = resolve_branch(tree_, clusterCandidates);
        if (sclusName.empty()) {
            std::cerr << "[RecoilProcessor] FATAL: no cluster branch found. Available:\n";
            for (auto& b : list_branches(tree_)) std::cerr << "  - " << b << "\n";
        } else {
            if (debug_>0) std::cout << "[RecoilProcessor] Using cluster branch: " << sclusName << "\n";
            // cluster branch: keep regular SetBranchAddress (has dictionary in this codebase)
            tree_->SetBranchAddress(sclusName.c_str(), &svtClusters_, &bsvtClusters_);
        }
    }

    // ---- Output ROOT file + histos (unchanged below)
    outF_ = TFile::Open("recoil_extra.root", "RECREATE");
    if (!outF_ || outF_->IsZombie()) { std::cerr << "[RecoilProcessor] ERROR opening output file.\n"; outF_ = nullptr; return; }
    outF_->cd();

    const double maxP = (beamE_>0 ? beamE_ : 3.74);
    h_rec_p_        = new TH1F("RecoilElectron_Momentum",    ";p (GeV);Counts", 100, 0., maxP);
    h_rec_theta_    = new TH1F("RecoilElectron_Theta",       ";#theta (deg);Counts", 180, 0., 180.);
    h_rec_vtz_      = new TH1F("RecoilElectron_VertexZ",     ";z_{vtx} (mm);Counts", 200, -20., 180.);
    h_rec_tanlam_p_ = new TH2F("RecoilElectron_TanLambda_vs_P",";tan#lambda;p (GeV);Counts",100, -0.3, 0.3, 100, 0., maxP);
    h_rec_nhits_    = new TH1F("RecoilElectron_nHits",       ";N truth hits;Counts", 20, -0.5, 19.5);
    h_rec_tanlam_vs_nhits_ = new TH2F("RecoilElectron_TanLambda_vs_nhits_",";tan#lambda;N Truth Hits; Counts",100, -0.3, 0.3,20, -0.5, 19.5);

    h_rad_p_        = new TH1F("RadiatedElectron_Momentum",  ";p (GeV);Counts", 100, 0., maxP);
    h_rad_theta_    = new TH1F("RadiatedElectron_Theta",     ";#theta (deg);Counts", 180, 0., 180.);
    h_rad_vtz_      = new TH1F("RadiatedElectron_VertexZ",   ";z_{vtx} (mm);Counts", 200, -20., 180.);
    h_rad_tanlam_p_ = new TH2F("RadiatedElectron_TanLambda_vs_P",";tan#lambda;p (GeV);Counts",100, -0.3, 0.3, 100, 0., maxP);

    const int tanLBins = 100;  const double tanLMin=-0.3, tanLMax=0.3;
    const int pBins    = 100;  const double pMin=0.0, pMax=maxP;

    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        eff_den_tanL_[i] = new TH1F(Form("EffDen_tanL_ge%d",thr), Form("Den: findable (>= %d truth layers);tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        eff_num_tanL_[i] = new TH1F(Form("EffNum_tanL_ge%d",thr), Form("Num: matched & >= %d layers;tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        eff_den_p_[i]    = new TH1F(Form("EffDen_p_ge%d",thr),    Form("Den: findable (>= %d truth layers);p (GeV);count",thr), pBins, pMin, pMax);
        eff_num_p_[i]    = new TH1F(Form("EffNum_p_ge%d",thr),    Form("Num: matched & >= %d layers;p (GeV);count",thr), pBins, pMin, pMax);
    }

    h_S_tanL_ = new TH1F("S_tanL", "Matched reco tracks (any);tan#lambda;S", tanLBins, tanLMin, tanLMax);
    h_B_tanL_ = new TH1F("B_tanL", "Unmatched tracks (overlap or pure fake);tan#lambda;B", tanLBins, tanLMin, tanLMax);
    h_S_p_    = new TH1F("S_p",    "Matched reco tracks (any);p (GeV);S",    pBins, pMin, pMax);
    h_B_p_    = new TH1F("B_p",    "Unmatched tracks (overlap or pure fake);p (GeV);B",    pBins, pMin, pMax);

    h_SoverSqrtB_tanL_ = nullptr;
    h_SoverSqrtB_p_    = nullptr;

    h_match_purity_    = new TH1F("MatchPurity", "Track hit purity to matched recoil;purity;tracks", 50, 0.0, 1.0);
    h_match_holes_     = new TH1F("MatchHoles",  "Missing truth layers on matched track;holes;tracks", kNLayersMax_+1, -0.5, kNLayersMax_+0.5);
    h_vz_vs_purity_    = new TH2F("VZ_vs_Purity", "v_{z} vs purity;v_{z} (mm);purity", 200, -20., 180., 50, 0.0, 1.0);
    h_tanL_vs_purity_  = new TH2F("TanL_vs_Purity", "tan#lambda vs purity;tan#lambda;purity", tanLBins, tanLMin, tanLMax, 50, 0.0, 1.0);

    g_match_completeness     = new TH1F("MatchCompleteness", "Truth-layer completeness on matched track;completeness;tracks", 50, 0.0, 1.0);
    g_match_matchedLayers    = new TH1F("MatchMatchedLayers","Matched layers to recoil on matched track;# matched layers;tracks", kNLayersMax_+1, -0.5, kNLayersMax_+0.5);
    g_comp_vs_purity         = new TH2F("Completeness_vs_Purity","Completeness vs purity;purity;completeness", 50, 0., 1., 50, 0., 1.);
    g_fakes_per_event        = new TH1F("FakesPerEvent","Per-event fake tracks (>= minThr hits, not matched);count;events", 50, -0.5, 49.5);
    g_duplicates_per_event   = new TH1F("DuplicatesPerEvent","Per-event duplicates (extra matches per recoil);count;events", 50, -0.5, 49.5);
    g_missedLayers_per_event = new TH1F("MissedLayersPerEvent","Per-event sum of holes over matched tracks;holes;events", 200, -0.5, 199.5);
    g_mixedHits_per_event    = new TH1F("MixedHitsPerEvent","Per-event sum of mixed hits over matched tracks;non-recoil hits;events", 500, -0.5, 499.5);

    h_B_overlap_tanL_ = new TH1F("Boverlap_tanL","Unmatched with truth overlap;tan#lambda;B", tanLBins, tanLMin, tanLMax);
    h_B_overlap_p_    = new TH1F("Boverlap_p",   "Unmatched with truth overlap;p (GeV);B", pBins, pMin, pMax);
    h_B_zero_tanL_    = new TH1F("Bzero_tanL",   "Unmatched with no truth (pure fake);tan#lambda;B", tanLBins, tanLMin, tanLMax);
    h_B_zero_p_       = new TH1F("Bzero_p",      "Unmatched with no truth (pure fake);p (GeV);B", pBins, pMin, pMax);

    recoilTree_ = new TTree("recoilTree","Truth‐level muon info");
    recoilTree_->Branch("sampleID",   &tree_sampleID_,  "sampleID/I");
    recoilTree_->Branch("category",   &tree_category_,  "category/I");
    recoilTree_->Branch("pdg",        &tree_pdg_,       "pdg/I");
    recoilTree_->Branch("OriginPDG",  &tree_OriginPDG_, "OriginPDG/I");
    recoilTree_->Branch("MomPDG",     &tree_MomPDG_,    "MomPDG/I");
    recoilTree_->Branch("charge",     &tree_charge_,    "charge/I");
    recoilTree_->Branch("ID",         &tree_ID_,        "ID/I");
    recoilTree_->Branch("px", &tree_px_,"px/F"); recoilTree_->Branch("py",&tree_py_,"py/F"); recoilTree_->Branch("pz",&tree_pz_,"pz/F");
    recoilTree_->Branch("p",&tree_p_,"p/F"); recoilTree_->Branch("pt",&tree_pt_,"pt/F");
    recoilTree_->Branch("tanLambda",&tree_tanLambda_,"tanLambda/F"); recoilTree_->Branch("theta",&tree_theta_,"theta/F");
    recoilTree_->Branch("energy",&tree_energy_,"energy/F"); recoilTree_->Branch("mass",&tree_mass_,"mass/F");
    recoilTree_->Branch("time",&tree_time_,"time/F"); recoilTree_->Branch("vz",&tree_vz_,"vz/F");
    recoilTree_->Branch("nHits",&tree_nHits_,"nHits/I");
    recoilTree_->Branch("x",&tree_x_,"x/F"); recoilTree_->Branch("phi",&tree_phi_,"phi/F");
    recoilTree_->Branch("eventID",&tree_eventID_,"eventID/I");
}














//void RecoilProcessor::initialize(TTree* tree) {
//    tree_ = tree;
//
//    // --- Branches we need
//    tree_->SetBranchAddress("EventHeader", &evth_, &bevth_);
//    tree_->SetBranchAddress("MCParticle",  &mcParts_, &bmcParts_);
//    tree_->SetBranchAddress("KalmanFullTracks", &trks_, &btrks_);
//    tree_->SetBranchAddress("TrackerSimHits",   &mcTrkrHits_, &bmcTrkrHits_);
//
//    
//    
//    
//
//
//// Replace ONLY the SAFE TRUTH BRANCH BIND block in RecoilProcessor::initialize with this:
//
///* ---- BEGIN SAFE TRUTH BRANCH BIND (minimal) ---- */
//{
//    // WHY: avoid segfault if ROOT can’t map the type or branch is missing.
//    truthTrks_  = nullptr;
//    bTruthTrks_ = nullptr;
//
//    static bool printed_branch_diag = false;
//
//    if (!truthTracksCollRoot_.empty()) {
//        TBranch* br = tree_ ? tree_->GetBranch(truthTracksCollRoot_.c_str()) : nullptr;
//        if (!br) {
//            if (debug_ > 0) {
//                std::cout << "[RecoilProcessor] Requested truth-track branch '" << truthTracksCollRoot_
//                          << "' not found.\n";
//                if (!printed_branch_diag) {
//                    std::cout << "[RecoilProcessor] Available branches:\n";
//                    for (const auto& bname : list_branches(tree_)) std::cout << "  - " << bname << "\n";
//                    printed_branch_diag = true;
//                }
//            }
//        } else {
//            TClass* cl = nullptr; EDataType dt = kOther_t;
//            const bool hasType = (br->GetExpectedType(cl, dt) == 0) && (cl != nullptr);
//            if (!hasType) {
//                if (debug_ > 0) {
//                    std::cout << "[RecoilProcessor] Branch '" << truthTracksCollRoot_
//                              << "' has no dictionary/type; skipping binding.\n";
//                }
//            } else {
//                const int rc = tree_->SetBranchAddress(truthTracksCollRoot_.c_str(), &truthTrks_, &bTruthTrks_);
//                if (rc != 0) {
//                    std::cout << "[RecoilProcessor] SetBranchAddress('" << truthTracksCollRoot_
//                              << "') failed (rc=" << rc << "); leaving unbound.\n";
//                } else if (debug_ > 0) {
//                    std::cout << "[RecoilProcessor] Using truth-track branch: "
//                              << truthTracksCollRoot_ << " (class=" << cl->GetName() << ")\n";
//                }
//            }
//        }
//    } else if (debug_ > 0) {
//        std::cout << "[RecoilProcessor] truthTrackCollRoot not set; TT-space disabled (using MC-space).\n";
//    }
//}
///* ---- END SAFE TRUTH BRANCH BIND (minimal) ---- */
//
//
//
/////* ---- BEGIN SAFE TRUTH BRANCH BIND ---- */
////{
////    //truthTrks_  = nullptr;
////    //bTruthTrks_ = nullptr;
////
////    if (!truthTracksCollRoot_.empty()) {
////        if (tree_->GetBranch(truthTracksCollRoot_.c_str())) {
////            tree_->SetBranchAddress(truthTracksCollRoot_.c_str(), &truthTrks_, &bTruthTrks_);
////            if (debug_>0) std::cout << "[RecoilProcessor] Using truth-track branch: "
////                                    << truthTracksCollRoot_ << "\n";
////        } else {
////            // One-time branch listing to help diagnose
////            if (debug_>0) {
////                std::cout << "[RecoilProcessor] Requested truth-track branch '" << truthTracksCollRoot_
////                          << "' not found. Available branches:\n";
////                for (const auto& bname : list_branches(tree_)) std::cout << "  - " << bname << "\n";
////            }
////            // Leave truthTrks_ null; code will fall back to MC-space safely.
////        }
////    } else if (debug_>0) {
////        std::cout << "[RecoilProcessor] truthTrackCollRoot not set; TT-space disabled (using MC-space).\n";
////    }
////}
/////* ---- END SAFE TRUTH BRANCH BIND ---- */
//
//
//
//
//
////    // TruthTracks branch written by TrackingProcessor under truthTracksCollRoot_
////    auto resolveTruthBranch = [&](TTree* t)->std::string {
////        if (!truthTracksCollRoot_.empty() && t->GetBranch(truthTracksCollRoot_.c_str()))
////            return truthTracksCollRoot_;
////        // Fallback: first branch containing both "truth" and "track"
////        auto names = list_branches(t);
////        for (const auto& nm : names) {
////            std::string lo = nm; std::transform(lo.begin(), lo.end(), lo.begin(), ::tolower);
////            if (lo.find("truth")!=std::string::npos && lo.find("track")!=std::string::npos) return nm;
////        }
////        return {};
////    };
////    const std::string ttName = resolveTruthBranch(tree_);
////    if (!ttName.empty()) {
////        tree_->SetBranchAddress(ttName.c_str(), &truthTrks_, &bTruthTrks_);
////        if (debug_>0) std::cout << "[RecoilProcessor] Using truth-track branch: " << ttName << "\n";
////    } else {
////        truthTrks_  = nullptr; bTruthTrks_ = nullptr;
////        if (debug_>0) std::cout << "[RecoilProcessor] No truth-track branch; TT-space efficiency disabled (fallback to MC).\n";
////    } 
//    
//    
//    
//    
//    
//   
//    // New TruthTracks with LCIO Truth Track IDs
//    //if (tree_->GetBranch("TruthTracks")) {
//    //    tree_->SetBranchAddress("TruthTracks", &truthTrks_, &bTruthTrks_);
//    //    if (debug_>0) std::cout << "[RecoilProcessor] Using truth branch: TruthTracks\n";
//    //} else {
//    //    truthTrks_  = nullptr;
//    //    bTruthTrks_ = nullptr;
//    //    if (debug_>0) std::cout << "[RecoilProcessor] TruthTracks branch not found; primary Truth-Track-ID path will be skipped.\n";
//    //}
//
//    // Resolve a cluster collection with MCPartIDs (prefer on-track)
//    {
//        std::vector<std::string> clusterCandidates;
//        if (!sclusColl_.empty()) clusterCandidates.push_back(sclusColl_);
//        clusterCandidates.push_back("SiClustersOnTrack");
//        clusterCandidates.push_back("SiClustersOnTrack_KF");
//        clusterCandidates.push_back("SvtClustersOnTrack");
//        clusterCandidates.push_back("SvtClustersOnTrack_KF");
//        clusterCandidates.push_back("SiClusters");
//        clusterCandidates.push_back("SvtClusters");
//        const std::string sclusName = resolve_branch(tree_, clusterCandidates);
//        if (sclusName.empty()) {
//            std::cerr << "[RecoilProcessor] FATAL: no cluster branch found. Available:\n";
//            for (auto& b : list_branches(tree_)) std::cerr << "  - " << b << "\n";
//        } else {
//            if (debug_>0) std::cout << "[RecoilProcessor] Using cluster branch: " << sclusName << "\n";
//            tree_->SetBranchAddress(sclusName.c_str(), &svtClusters_, &bsvtClusters_);
//        }
//    }
//
//    // --- Output ROOT file
//    outF_ = TFile::Open("recoil_extra.root", "RECREATE");
//    if (!outF_ || outF_->IsZombie()) {
//        std::cerr << "[RecoilProcessor] ERROR opening output file.\n";
//        outF_ = nullptr; return;
//    }
//    outF_->cd();
//
//    // --- Histograms (unchanged) ------------------------------------------------
//    const double maxP = (beamE_>0 ? beamE_ : 3.74);
//    h_rec_p_        = new TH1F("RecoilElectron_Momentum",    ";p (GeV);Counts", 100, 0., maxP);
//    h_rec_theta_    = new TH1F("RecoilElectron_Theta",       ";#theta (deg);Counts", 180, 0., 180.);
//    h_rec_vtz_      = new TH1F("RecoilElectron_VertexZ",     ";z_{vtx} (mm);Counts", 200, -20., 180.);
//    h_rec_tanlam_p_ = new TH2F("RecoilElectron_TanLambda_vs_P",";tan#lambda;p (GeV);Counts",100, -0.3, 0.3, 100, 0., maxP);
//    h_rec_nhits_    = new TH1F("RecoilElectron_nHits",       ";N truth hits;Counts", 20, -0.5, 19.5);
//    h_rec_tanlam_vs_nhits_ = new TH2F("RecoilElectron_TanLambda_vs_nhits_",";tan#lambda;N Truth Hits; Counts",100, -0.3, 0.3,20, -0.5, 19.5);
//
//    h_rad_p_        = new TH1F("RadiatedElectron_Momentum",  ";p (GeV);Counts", 100, 0., maxP);
//    h_rad_theta_    = new TH1F("RadiatedElectron_Theta",     ";#theta (deg);Counts", 180, 0., 180.);
//    h_rad_vtz_      = new TH1F("RadiatedElectron_VertexZ",   ";z_{vtx} (mm);Counts", 200, -20., 180.);
//    h_rad_tanlam_p_ = new TH2F("RadiatedElectron_TanLambda_vs_P",";tan#lambda;p (GeV);Counts",100, -0.3, 0.3, 100, 0., maxP);
//
//    const int tanLBins = 100;  const double tanLMin=-0.3, tanLMax=0.3;
//    const int pBins    = 100;  const double pMin=0.0, pMax=maxP;
//
//    for (int i = 0; i < kNThr_; ++i) {
//        const int thr = thrList_[i];
//        eff_den_tanL_[i] = new TH1F(Form("EffDen_tanL_ge%d",thr), Form("Den: findable (>= %d truth layers);tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
//        eff_num_tanL_[i] = new TH1F(Form("EffNum_tanL_ge%d",thr), Form("Num: matched & >= %d layers;tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
//        eff_den_p_[i]    = new TH1F(Form("EffDen_p_ge%d",thr),    Form("Den: findable (>= %d truth layers);p (GeV);count",thr), pBins, pMin, pMax);
//        eff_num_p_[i]    = new TH1F(Form("EffNum_p_ge%d",thr),    Form("Num: matched & >= %d layers;p (GeV);count",thr), pBins, pMin, pMax);
//    }
//
//    h_S_tanL_ = new TH1F("S_tanL", "Matched reco tracks (any);tan#lambda;S", tanLBins, tanLMin, tanLMax);
//    h_B_tanL_ = new TH1F("B_tanL", "Unmatched tracks (overlap or pure fake);tan#lambda;B", tanLBins, tanLMin, tanLMax);
//    h_S_p_    = new TH1F("S_p",    "Matched reco tracks (any);p (GeV);S",    pBins, pMin, pMax);
//    h_B_p_    = new TH1F("B_p",    "Unmatched tracks (overlap or pure fake);p (GeV);B",    pBins, pMin, pMax);
//
//    h_SoverSqrtB_tanL_ = nullptr;
//    h_SoverSqrtB_p_    = nullptr;
//
//    h_match_purity_    = new TH1F("MatchPurity", "Track hit purity to matched recoil;purity;tracks", 50, 0.0, 1.0);
//    h_match_holes_     = new TH1F("MatchHoles",  "Missing truth layers on matched track;holes;tracks", kNLayersMax_+1, -0.5, kNLayersMax_+0.5);
//    h_vz_vs_purity_    = new TH2F("VZ_vs_Purity", "v_{z} vs purity;v_{z} (mm);purity", 200, -20., 180., 50, 0.0, 1.0);
//    h_tanL_vs_purity_  = new TH2F("TanL_vs_Purity", "tan#lambda vs purity;tan#lambda;purity", tanLBins, tanLMin, tanLMax, 50, 0.0, 1.0);
//
//    g_match_completeness     = new TH1F("MatchCompleteness", "Truth-layer completeness on matched track;completeness;tracks", 50, 0.0, 1.0);
//    g_match_matchedLayers    = new TH1F("MatchMatchedLayers","Matched layers to recoil on matched track;# matched layers;tracks", kNLayersMax_+1, -0.5, kNLayersMax_+0.5);
//    g_comp_vs_purity         = new TH2F("Completeness_vs_Purity","Completeness vs purity;purity;completeness", 50, 0., 1., 50, 0., 1.);
//    g_fakes_per_event        = new TH1F("FakesPerEvent","Per-event fake tracks (>= minThr hits, not matched);count;events", 50, -0.5, 49.5);
//    g_duplicates_per_event   = new TH1F("DuplicatesPerEvent","Per-event duplicates (extra matches per recoil);count;events", 50, -0.5, 49.5);
//    g_missedLayers_per_event = new TH1F("MissedLayersPerEvent","Per-event sum of holes over matched tracks;holes;events", 200, -0.5, 199.5);
//    g_mixedHits_per_event    = new TH1F("MixedHitsPerEvent","Per-event sum of mixed hits over matched tracks;non-recoil hits;events", 500, -0.5, 499.5);
//
//    // Optional diagnostics (B split)
//    h_B_overlap_tanL_ = new TH1F("Boverlap_tanL","Unmatched with truth overlap;tan#lambda;B", tanLBins, tanLMin, tanLMax);
//    h_B_overlap_p_    = new TH1F("Boverlap_p",   "Unmatched with truth overlap;p (GeV);B", pBins, pMin, pMax);
//    h_B_zero_tanL_    = new TH1F("Bzero_tanL",   "Unmatched with no truth (pure fake);tan#lambda;B", tanLBins, tanLMin, tanLMax);
//    h_B_zero_p_       = new TH1F("Bzero_p",      "Unmatched with no truth (pure fake);p (GeV);B", pBins, pMin, pMax);
//
//    // --- Output ntuple of truth recoils
//    recoilTree_ = new TTree("recoilTree","Truth‐level muon info");
//    recoilTree_->Branch("sampleID",   &tree_sampleID_,  "sampleID/I");
//    recoilTree_->Branch("category",   &tree_category_,  "category/I");
//    recoilTree_->Branch("pdg",        &tree_pdg_,       "pdg/I");
//    recoilTree_->Branch("OriginPDG",  &tree_OriginPDG_, "OriginPDG/I");
//    recoilTree_->Branch("MomPDG",     &tree_MomPDG_,    "MomPDG/I");
//    recoilTree_->Branch("charge",     &tree_charge_,    "charge/I");
//    recoilTree_->Branch("ID",         &tree_ID_,        "ID/I");
//    recoilTree_->Branch("px", &tree_px_,"px/F"); recoilTree_->Branch("py",&tree_py_,"py/F"); recoilTree_->Branch("pz",&tree_pz_,"pz/F");
//    recoilTree_->Branch("p",&tree_p_,"p/F"); recoilTree_->Branch("pt",&tree_pt_,"pt/F");
//    recoilTree_->Branch("tanLambda",&tree_tanLambda_,"tanLambda/F"); recoilTree_->Branch("theta",&tree_theta_,"theta/F");
//    recoilTree_->Branch("energy",&tree_energy_,"energy/F"); recoilTree_->Branch("mass",&tree_mass_,"mass/F");
//    recoilTree_->Branch("time",&tree_time_,"time/F"); recoilTree_->Branch("vz",&tree_vz_,"vz/F");
//    recoilTree_->Branch("nHits",&tree_nHits_,"nHits/I");
//    recoilTree_->Branch("x",&tree_x_,"x/F"); recoilTree_->Branch("phi",&tree_phi_,"phi/F");
//    recoilTree_->Branch("eventID",&tree_eventID_,"eventID/I");
//}
//
bool RecoilProcessor::process(IEvent* /*ievent*/) {
    if (!evth_ || !mcParts_ || !svtClusters_ || !trks_) {
        std::cerr << "[RecoilProcessor] Missing required collections (MC/Clusters/Tracks)\n";
        return false;
    }

    // --- Compute min threshold from fixed thrList_
    int minThrLocal = thrList_[0];
    for (int i = 1; i < kNThr_; ++i) minThrLocal = std::min(minThrLocal, thrList_[i]);

    // ---- (1) Select recoils: PDG=13, remove beam-energy lepton
    const double ETOL  = 1e-3;
    const double BEAME = (beamE_>0 ? beamE_ : 3.74);

    std::vector<MCParticle*> recoils;
    recoils.reserve(mcParts_->size());
    for (auto* part : *mcParts_) {
        if (!part) continue;
        if (std::abs(part->getPDG()) != 13) continue;                       // keep muons
        if (std::abs(part->getEnergy() - BEAME) < ETOL) continue;           // drop beam lepton
        recoils.push_back(part);
    }
    if (recoils.empty()) return true;

    // ---- (2) Truth masks (prefer clusters; fallback to sim-hits if clusters lack MCPartIDs)
    struct TruthInfo {
        int nTruthLayers{0};
        std::array<bool, 64> truthMask{};    // generous upper bound
        double tl_true{999.};
        double p_true{0.};
    };
    std::unordered_map<int,TruthInfo> truthByRecoil; truthByRecoil.reserve(recoils.size());
    for (auto* recoil : recoils) {
        TruthInfo ti; ti.truthMask.fill(false);
        const auto& p = recoil->getMomentum();
        // HPS convention for tanλ truth: py/pz
        ti.tl_true = (p[2]!=0.0) ? (p[1]/p[2]) : 999.;
        ti.p_true  = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        truthByRecoil[recoil->getID()] = std::move(ti);
    }

    // Detect if cluster collection has MCPartIDs
    int clusters_with_mc = 0;
    for (size_t i=0;i<svtClusters_->size();++i) {
        auto* c = svtClusters_->at(i);
        if (c && !c->getMCPartIDs().empty()) { clusters_with_mc++; break; }
    }

    if (clusters_with_mc > 0) {
        for (auto* clu : *svtClusters_) {
            if (!clu) continue;
            const int L = clu->getLayer();
            if (L < 0 || L >= kNLayersMax_) continue;
            const auto& ids = clu->getMCPartIDs();
            if (ids.empty()) continue;
            for (int rid : ids) {
                auto it = truthByRecoil.find(rid);
                if (it == truthByRecoil.end()) continue;
                it->second.truthMask[L] = true;
            }
        }
    } else if (mcTrkrHits_) {
        // Fallback: use sim-hits to mark truth layers (IDs must match MCParticle::getID()).
        for (auto* sh : *mcTrkrHits_) {
            if (!sh) continue;
            const int rid = 0; //sh->getID(); Disabled
            const int L   = sh->getLayer();
            if (L < 0 || L >= kNLayersMax_) continue;
            auto it = truthByRecoil.find(rid);
            if (it == truthByRecoil.end()) continue;
            it->second.truthMask[L] = true;
        }
    }

    // finalize counts + fill denominators
    //for (auto& kv : truthByRecoil) {
    //    auto& ti = kv.second;
    //    ti.nTruthLayers = 0;
    //    for (int L=0; L<kNLayersMax_; ++L) ti.nTruthLayers += (ti.truthMask[L] ? 1 : 0);
    //    for (int t=0; t<kNThr_; ++t) {
    //        const int thr = thrList_[t];
    //        if (ti.nTruthLayers >= thr) {
    //            eff_den_tanL_[t]->Fill(ti.tl_true);
    //            eff_den_p_[t]->Fill(ti.p_true);
    //        }
    //    }
    //}


    // Recompute counts on MC masks (used later), but don't fill DEN here
    for (auto& kv : truthByRecoil) {
        auto& ti = kv.second;
        ti.nTruthLayers = 0;
        for (int L=0; L<kNLayersMax_; ++L) ti.nTruthLayers += (ti.truthMask[L] ? 1 : 0);
    }

    // Build truth-track masks + TT→MC map by layer overlap
    //struct TruthInfo { int nTruthLayers{0}; std::array<bool,64> truthMask{}; double tl_true{999.}; double p_true{0.}; };
    std::unordered_map<int, TruthInfo> truthByTruthTrackId;   // TT id -> mask/kin
    std::unordered_map<int, int>       tt2mc_best;            // TT id -> MCParticle id
    if (truthTrks_) {
        truthByTruthTrackId.reserve(truthTrks_->size());
        for (auto* tt : *truthTrks_) {
            if (!tt) continue;
            TruthInfo ti; ti.truthMask.fill(false);
            for (int L : tt->getHitLayers()) if (0<=L && L<kNLayersMax_) ti.truthMask[L] = true;
            ti.nTruthLayers = std::count(ti.truthMask.begin(), ti.truthMask.begin()+kNLayersMax_, true);
            ti.tl_true = tt->getTanLambda();
            ti.p_true  = tt->getP();
            truthByTruthTrackId[ tt->getID() ] = std::move(ti);
        }
        // map each TT to best recoil MC by overlap
        for (const auto& kvTT : truthByTruthTrackId) {
            const int ttId = kvTT.first; const auto& tiTT = kvTT.second;
            int bestMc=-1, bestOv=-1;
            for (const auto& kvMC : truthByRecoil) {
                int ov=0; for (int L=0; L<kNLayersMax_; ++L) if (tiTT.truthMask[L] && kvMC.second.truthMask[L]) ov++;
                if (ov>bestOv) { bestOv=ov; bestMc=kvMC.first; }
            }
            if (bestMc>=0 && bestOv>0) tt2mc_best[ttId]=bestMc;
        }
        // DEN in TT-space restricted to TTs that map to one of our selected recoils
        for (const auto& kv : truthByTruthTrackId) {
            const int ttId = kv.first; const auto& tiTT = kv.second;
            auto itMap = tt2mc_best.find(ttId); if (itMap==tt2mc_best.end()) continue;
            if (truthByRecoil.find(itMap->second)==truthByRecoil.end()) continue;
            for (int t=0; t<kNThr_; ++t) { const int thr = thrList_[t];
                if (tiTT.nTruthLayers >= thr) { eff_den_tanL_[t]->Fill(tiTT.tl_true); eff_den_p_[t]->Fill(tiTT.p_true); }
            }
        }
    } else {
        // fallback: keep old DEN in MC space if TT absent
        for (const auto& kv : truthByRecoil) {
            const auto& ti = kv.second;
            for (int t=0; t<kNThr_; ++t) { const int thr = thrList_[t];
                if (ti.nTruthLayers >= thr) { eff_den_tanL_[t]->Fill(ti.tl_true); eff_den_p_[t]->Fill(ti.p_true); }
            }
        }
    }






    // ---- New: build truth-track space + TT->MC mapping by layer overlap ----------
    //std::unordered_map<int, TruthInfo> truthByTruthTrackId;   // key: TruthTrack ID
    //std::unordered_map<int, int>       tt2mc_best;            // key: TruthTrack ID -> MCParticle ID
    if (truthTrks_) {
        truthByTruthTrackId.reserve(truthTrks_->size());
        for (auto* tt : *truthTrks_) {
            if (!tt) continue;
            TruthInfo ti; ti.truthMask.fill(false);
            // Layers from truth track object
            for (int L : tt->getHitLayers()) {
                if (0 <= L && L < kNLayersMax_) ti.truthMask[L] = true;
            }
            ti.nTruthLayers = std::count(ti.truthMask.begin(), ti.truthMask.begin()+kNLayersMax_, true);

            // Truth-track kinematics (kept for potential debugging)
            ti.tl_true = tt->getTanLambda();
            ti.p_true  = tt->getP();
            truthByTruthTrackId[tt->getID()] = std::move(ti);
        }

        // Map each truth-track to the MCParticle with maximum layer overlap.
        for (const auto& kvTT : truthByTruthTrackId) {
            const int ttId = kvTT.first;
            const auto& tiTT = kvTT.second;
            int bestMc = -1, bestOv = -1;
            for (const auto& kvMC : truthByRecoil) {
                const int mcId = kvMC.first;
                const auto& tiMC = kvMC.second;
                int ov = 0;
                for (int L=0; L<kNLayersMax_; ++L) {
                    if (tiTT.truthMask[L] && tiMC.truthMask[L]) ov++;
                }
                if (ov > bestOv) { bestOv = ov; bestMc = mcId; }
            }
            if (bestMc >= 0 && bestOv > 0) tt2mc_best[ttId] = bestMc;
        }
    }

    // ---- (3) Map clusterID -> cluster* (for per-track hit lookup)
    std::unordered_map<ULong64_t, ::TrackerHit*> id2cluster; id2cluster.reserve(svtClusters_->size());
    for (size_t i=0;i<svtClusters_->size();++i) {
        auto* c = svtClusters_->at(i);
        if (!c) continue;
        id2cluster[(ULong64_t)c->getID()] = c;
    }
    const bool useClusterMCIDs = (clusters_with_mc > 0);
    if (!useClusterMCIDs && debug_>0) {
        std::cerr << "[RecoilProcessor] Using layer-overlap fallback (no MCPartIDs on clusters).\n";
    }

    // ---- (4) Per-track association
    struct AssocTrack {
        const Track* trk{nullptr};
        int bestRid{-1};                 // MCParticle ID (by design)
        int totalHits{0};
        int matchedHits{0};
        int matchedLayers{0};
        double purity{0.0};
        double completeness{0.0};
        double tanL{0.};
        double p{0.};
        bool hadAnyMC{false};
        std::array<bool, 64> layerHasBest{};
        int firstMCLayer{-1};
        int totalHitsWithAnyMC{0};
    };
    auto score = [&](const AssocTrack& a)->double{
        return 0.6*a.purity + 0.3*a.completeness + 0.1*(double(a.matchedLayers)/double(kNLayersMax_));
    };
    const double kPurityCut = 0.0;  // OFF

    std::unordered_map<int, std::vector<AssocTrack>> candidatesByRecoil; // mc rid -> candidates
    candidatesByRecoil.reserve(recoils.size());
    int eventFakes = 0, eventDuplicates = 0, eventMissed = 0, eventMixed = 0;

    for (auto* trk : *trks_) {
        if (!trk) continue;
        const int T = trk->getSvtHits().GetEntriesFast();
        if (T <= 0) continue;

        // Build track layer mask and per-hit layer list
        std::array<bool, 64> trackLayerMask{}; trackLayerMask.fill(false);
        std::vector<int> trackHitLayers; trackHitLayers.reserve(T);
        for (int ih=0; ih<T; ++ih) {
            auto* th = static_cast<::TrackerHit*>(trk->getSvtHits().At(ih));
            if (!th) continue;
            const int L = th->getLayer();
            if (L >= 0 && L < kNLayersMax_) {
                trackLayerMask[L] = true;
                trackHitLayers.push_back(L);
            }
        }

        struct Acc { int hits=0; std::array<bool,64> L{}; Acc(){L.fill(false);} };
        std::unordered_map<int, Acc> accByRid; // mc rid -> counts
        bool hadAnyTruth = false;
        int  assocUnits  = 0;
        int  onlyAssocLayer = -1;

        // --- (A) Truth-Track-ID path (now using TT→MC mapping)
        bool usedTruthIdPath = false;
        if (truthTrks_) {
            const int ttId_from_track = trk->getID(); // set by TrackingProcessor to truth track id
            auto itMap = tt2mc_best.find(ttId_from_track);
            if (itMap != tt2mc_best.end()) {
                const int mcRid = itMap->second;
                auto itTruthMC = truthByRecoil.find(mcRid);
                if (itTruthMC != truthByRecoil.end()) {
                    usedTruthIdPath = true;
                    const auto& tiMC = itTruthMC->second;

                    Acc a;
                    for (int L = 0; L < kNLayersMax_; ++L) {
                        if (trackLayerMask[L] && tiMC.truthMask[L]) a.L[L] = true;
                    }
                    int matchedHitsCount = 0;
                    for (int L_hit : trackHitLayers) {
                        if (L_hit >= 0 && L_hit < kNLayersMax_ && a.L[L_hit]) matchedHitsCount++;
                    }
                    a.hits = matchedHitsCount;

                    if (a.hits > 0) {
                        hadAnyTruth = true;
                        assocUnits  = a.hits;
                        if (a.hits == 1) {
                            for (int L=0; L<kNLayersMax_; ++L) { if (a.L[L]) { onlyAssocLayer = L; break; } }
                        }
                    }
                    accByRid[mcRid] = a; // downstream remains in MCParticle id space
                }
            }
        }

        // --- (B) Cluster MCPartIDs path (unchanged)
        if (!usedTruthIdPath && useClusterMCIDs) {
            int total_hits_with_any_mc = 0;
            int only_mc_hit_layer = -1;
            for (int ih=0; ih<T; ++ih) {
                auto* th = static_cast<::TrackerHit*>(trk->getSvtHits().At(ih));
                if (!th) continue;
                auto itc = id2cluster.find((ULong64_t)th->getID());
                if (itc == id2cluster.end()) continue;
                ::TrackerHit* clu = itc->second;
                const int L = clu->getLayer();
                const auto& mcids = clu->getMCPartIDs();
                if (!mcids.empty()) {
                    hadAnyTruth = true;
                    total_hits_with_any_mc++;
                    if (total_hits_with_any_mc == 1) only_mc_hit_layer = L;
                }
                for (int rid : mcids) {
                    auto itTruth = truthByRecoil.find(rid);
                    if (itTruth == truthByRecoil.end()) continue;
                    auto& slot = accByRid[rid];
                    slot.hits += 1;
                    if (L>=0 && L<kNLayersMax_) slot.L[L] = true;
                }
            }
            assocUnits     = total_hits_with_any_mc;
            onlyAssocLayer = only_mc_hit_layer;
        }

        // --- (C) Layer-overlap fallback across all recoils (unchanged)
        if (!usedTruthIdPath && !useClusterMCIDs) {
            int bestRidLoc = -1, bestOverlap = -1;
            for (const auto& kvTruth : truthByRecoil) {
                const int rid = kvTruth.first;
                const auto& ti = kvTruth.second;
                int overlap = 0;
                for (int L=0; L<kNLayersMax_; ++L) {
                    if (trackLayerMask[L] && ti.truthMask[L]) overlap++;
                }
                if (overlap > bestOverlap) { bestOverlap = overlap; bestRidLoc = rid; }
            }
            if (bestRidLoc >= 0 && bestOverlap > 0) {
                hadAnyTruth = true;
                assocUnits  = bestOverlap;
                if (bestOverlap == 1) {
                    for (int L=0; L<kNLayersMax_; ++L) {
                        if (trackLayerMask[L] && truthByRecoil[bestRidLoc].truthMask[L]) { onlyAssocLayer = L; break; }
                    }
                }
                Acc a;
                for (int L=0; L<kNLayersMax_; ++L) {
                    if (trackLayerMask[L] && truthByRecoil[bestRidLoc].truthMask[L]) a.L[L] = true;
                }
                int matchedHitsCount = 0;
                for (int L_hit : trackHitLayers) {
                    if (L_hit >= 0 && L_hit < kNLayersMax_ && a.L[L_hit]) matchedHitsCount++;
                }
                a.hits = matchedHitsCount;
                accByRid[bestRidLoc] = a;
            }
        }

        // --- Choose best rid by hits then layers
        int bestRid = -1, bestHits = -1, bestLayers = -1;
        for (const auto& kv : accByRid) {
            const int rid = kv.first;
            const auto& a = kv.second;
            int layers=0; for (int L=0; L<kNLayersMax_; ++L) layers += (a.L[L]?1:0);
            if (a.hits > bestHits || (a.hits==bestHits && layers>bestLayers)) {
                bestRid = rid; bestHits = a.hits; bestLayers = layers;
            }
        }

        // --- Build AssocTrack and downstream logic
        AssocTrack AT;
        AT.trk            = trk;
        AT.bestRid        = bestRid;                    // MCParticle id
        AT.totalHits      = T;
        AT.matchedHits    = (bestRid>=0 ? bestHits : 0);
        AT.matchedLayers  = (bestRid>=0 ? bestLayers : 0);
        AT.purity         = (T>0 ? double(AT.matchedHits)/double(T) : 0.0);
        AT.hadAnyMC       = hadAnyTruth;
        AT.tanL           = trk->getTanLambda();
        AT.p              = trk->getP();
        AT.layerHasBest.fill(false);
        AT.firstMCLayer   = onlyAssocLayer;
        AT.totalHitsWithAnyMC = assocUnits;

        if (bestRid>=0) {
            const auto& Lmask = accByRid[bestRid].L;
            for (int L=0; L<kNLayersMax_; ++L) AT.layerHasBest[L] = Lmask[L];
            const auto& ti = truthByRecoil[bestRid];
            AT.completeness = (ti.nTruthLayers>0 ? double(AT.matchedLayers)/double(ti.nTruthLayers) : 0.0);
        }

        // ---- Decide per-track S/B and candidate status ----
        bool anyTruth = (AT.matchedLayers > 0);
        int  overlapLayers = AT.matchedLayers;

        bool truthFindableAtMin = (AT.bestRid>=0) && (truthByRecoil[AT.bestRid].nTruthLayers >= minThrLocal);
        bool halfHits           = (AT.bestRid>=0) && (AT.matchedHits >= (AT.totalHits+1)/2);
        bool matched            = (halfHits && truthFindableAtMin);

        bool front_single_mc = (AT.totalHitsWithAnyMC == 1) && (AT.firstMCLayer == 0 || AT.firstMCLayer == 1);

        if (matched) {
            h_S_tanL_->Fill(AT.tanL);
            h_S_p_->Fill(AT.p);
        } else {
            if ( anyTruth || (!anyTruth && AT.totalHits >= minThrLocal) ) {
                h_B_tanL_->Fill(AT.tanL);
                h_B_p_->Fill(AT.p);
            }
            if (anyTruth) { h_B_overlap_tanL_->Fill(AT.tanL); h_B_overlap_p_->Fill(AT.p); }
            else if (AT.totalHits >= minThrLocal) { h_B_zero_tanL_->Fill(AT.tanL); h_B_zero_p_->Fill(AT.p); }
        }

        if (matched) {
            candidatesByRecoil[AT.bestRid].push_back(AT);
        } else {
            if (!front_single_mc && AT.totalHits >= minThrLocal && anyTruth) eventFakes += 1;
        }
    }

    // Duplicates: extra candidates per recoil
    for (const auto& kv : candidatesByRecoil) {
        const int n = (int)kv.second.size();
        if (n>1) eventDuplicates += (n-1);
    }

    // ---- (5) Choose best per recoil; fill NUM and diagnostics
    for (auto* recoil : recoils) {
        const int rid = recoil->getID();
        auto trIt = truthByRecoil.find(rid);
        if (trIt == truthByRecoil.end()) continue;
        const auto& ti = trIt->second;

        auto candIt = candidatesByRecoil.find(rid);
        if (candIt == candidatesByRecoil.end() || candIt->second.empty()) continue;

        const AssocTrack* best = &candIt->second[0];
        for (const auto& a : candIt->second) if (score(a) > score(*best)) best = &a;

        int holes = 0;
        for (int L=0; L<kNLayersMax_; ++L) {
            if (!ti.truthMask[L]) continue;
            if (!best->layerHasBest[L]) holes++;
        }

        h_match_purity_->Fill(best->purity);
        h_match_holes_->Fill(holes);
        h_tanL_vs_purity_->Fill(best->tanL, best->purity);
        h_vz_vs_purity_->Fill(recoil->getVertexPosition()[2], best->purity);

        g_match_completeness->Fill(best->completeness);
        g_match_matchedLayers->Fill(best->matchedLayers);
        g_comp_vs_purity->Fill(best->purity, best->completeness);

        eventMissed += holes;
        eventMixed  += (best->totalHits - best->matchedHits);

        // Fill numerator for ALL thresholds satisfied by BOTH truth and matched layers
       // for (int t=0; t<kNThr_; ++t) {
       //     const int need = thrList_[t];
       //     if (ti.nTruthLayers >= need && best->matchedLayers >= need) {
       //         eff_num_tanL_[t]->Fill(ti.tl_true);
       //         eff_num_p_[t]->Fill(ti.p_true);
       //     }
       // }
       
        for (int t=0; t<kNThr_; ++t) {
            const int need = thrList_[t];
            if (best->matchedLayers < need) continue;
            bool filled=false;
            if (truthTrks_) {
                const int ttId = best->trk ? best->trk->getID() : -1; // Track::getID()==TruthTrack id
                auto itMap = tt2mc_best.find(ttId);
                if (ttId>=0 && itMap!=tt2mc_best.end() && itMap->second==rid) {
                    auto itTT = truthByTruthTrackId.find(ttId);
                    if (itTT!=truthByTruthTrackId.end() && itTT->second.nTruthLayers>=need) {
                        eff_num_tanL_[t]->Fill(itTT->second.tl_true);
                        eff_num_p_[t]->Fill(itTT->second.p_true);
                        filled=true;
                    }
                }
            }
            if (!filled && ti.nTruthLayers>=need) {
                // fallback to MC kinematics (keeps behavior if TT missing)
                eff_num_tanL_[t]->Fill(ti.tl_true);
                eff_num_p_[t]->Fill(ti.p_true);
            }
        }








    }

    // ---- (6) Trees & simple kinematics
    std::unordered_map<int,MCParticle*> id2part; id2part.reserve(mcParts_->size());
    for (auto* p : *mcParts_) if (p) id2part[p->getID()] = p;

    std::unordered_map<MCParticle*,int> mcTruthHits;
    if (mcTrkrHits_) {
        for (auto* hit : *mcTrkrHits_) {
            if (!hit) continue;
            int pid = 0; // hit->getID(); // Disabled
            auto it = id2part.find(pid);
            if (it != id2part.end()) mcTruthHits[it->second]++;
        }
    }

    auto fillTree = [&](MCParticle* e, int category){
        tree_sampleID_  = sampleID_;
        tree_category_  = category;
        tree_pdg_       = e->getPDG();
        tree_OriginPDG_ = e->getOriginPDG();
        tree_MomPDG_    = e->getMomPDG();
        tree_charge_    = e->getCharge();
        tree_ID_        = e->getID();

        auto mom = e->getMomentum();
        tree_px_       = mom[0];
        tree_py_       = mom[1];
        tree_pz_       = mom[2];
        tree_p_        = std::sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
        tree_pt_       = std::sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
        tree_tanLambda_= (tree_pz_!=0 ? tree_py_/tree_pz_ : 999.);
        tree_theta_    = std::atan2(tree_pt_,tree_pz_) * TMath::RadToDeg();

        tree_energy_   = e->getEnergy();
        tree_mass_     = e->getMass();
        tree_time_     = e->getTime();
        tree_vz_       = e->getVertexPosition()[2];

        tree_nHits_    = mcTruthHits[e];
        tree_eventID_  = evth_->getEventNumber();
        tree_x_        = 1.0 - (e->getEnergy() / beamE_);
        tree_phi_      = std::atan2(tree_py_, tree_px_);

        recoilTree_->Fill();
    };

    for (auto* e : recoils)  fillTree(e, 0); // category=0 for recoils

    auto fillKin = [&](MCParticle* ele, bool isRecoil) {
        const auto& mom = ele->getMomentum();
        double px=mom[0], py=mom[1], pz=mom[2];
        double p  = std::sqrt(px*px + py*py + pz*pz);
        double pt = std::sqrt(px*px + py*py);
        double tanL = (pz != 0 ? py/pz : 999.);
        double thetaDeg = std::atan2(pt,pz)*TMath::RadToDeg();
        const auto& vtx = ele->getVertexPosition(); double vz = (vtx.size()>2)? vtx[2] : 0.;
        if (isRecoil) { h_rec_p_->Fill(p); h_rec_theta_->Fill(thetaDeg); h_rec_vtz_->Fill(vz); h_rec_tanlam_p_->Fill(tanL, p); }
        else          { h_rad_p_->Fill(p); h_rad_theta_->Fill(thetaDeg); h_rad_vtz_->Fill(vz); h_rad_tanlam_p_->Fill(tanL, p); }
    };
    for (auto* e : recoils)  fillKin(e, true);

    for (auto* e : recoils) {
        int nh = mcTruthHits[e];
        h_rec_nhits_->Fill(nh);
        auto mom=e->getMomentum();
        double tanL = (mom[2]!=0. ? mom[1]/mom[2] : 999.);
        h_rec_tanlam_vs_nhits_->Fill(tanL, nh);
    }

    // ---- (7) Per-event diagnostics
    g_fakes_per_event->Fill(eventFakes);
    g_duplicates_per_event->Fill(eventDuplicates);
    g_missedLayers_per_event->Fill(eventMissed);
    g_mixedHits_per_event->Fill(eventMixed);

    if (debug_>0) {
        std::cout << "[RecoilProcessor] Event " << evth_->getEventNumber()
                  << "  recoils=" << recoils.size()
                  << "  clusters_with_mc=" << clusters_with_mc
                  << "  fakes=" << eventFakes
                  << "  duplicates=" << eventDuplicates
                  << "  missed(sum)=" << eventMissed
                  << "  mixed(sum)=" << eventMixed << "\n";
    }

    // Extra ID debug (first few only)
    if (debug_>1) {
        std::cerr << "Example trk IDs (TruthTrack IDs): ";
        int shown=0; for (auto* trk : *trks_) { std::cerr << trk->getID() << " "; if (++shown==10) break; }
        std::cerr << "\nExample TruthTrack IDs (from branch): ";
        shown=0; if (truthTrks_) { for (auto* tt : *truthTrks_) { std::cerr << tt->getID() << " "; if (++shown==10) break; } }
        std::cerr << "\nExample MCParticle IDs: ";
        shown=0; for (auto* p : *mcParts_) { std::cerr << p->getID() << " "; if (++shown==10) break; }
        std::cerr << "\n";
    }

    return true;
}

void RecoilProcessor::finalize() {
    if (!outF_ || !outF_->IsOpen()) { std::cerr << "[RecoilProcessor] finalize(): output not open\n"; return; }
    outF_->cd();
    auto writeIf = [](TObject* o){ if (o) o->Write(); };

    // Make S/sqrt(B) overlays
    auto makeSoverSqrtB = [](TH1F* S, TH1F* B, const char* name, const char* title)->TH1F*{
        TH1F* R = (TH1F*)S->Clone(name); R->SetTitle(title);
        for (int i=1;i<=R->GetNbinsX();++i) {
            const double s=S->GetBinContent(i), b=B->GetBinContent(i);
            R->SetBinContent(i, (b>0.0)? s/std::sqrt(b) : 0.0);
            R->SetBinError(i, 0.0);
        }
        return R;
    };
    h_SoverSqrtB_tanL_ = makeSoverSqrtB(h_S_tanL_, h_B_tanL_, "SoverSqrtB_tanL", "S/#sqrt{B} vs tan#lambda");
    h_SoverSqrtB_p_    = makeSoverSqrtB(h_S_p_,    h_B_p_,    "SoverSqrtB_p" ,   "S/#sqrt{B} vs p");

    // Write trees and histos
    writeIf(recoilTree_);

    writeIf(h_rec_p_); writeIf(h_rec_theta_); writeIf(h_rec_vtz_);
    writeIf(h_rec_tanlam_p_); writeIf(h_rec_nhits_); writeIf(h_rec_tanlam_vs_nhits_);

    writeIf(h_rad_p_); writeIf(h_rad_theta_); writeIf(h_rad_vtz_); writeIf(h_rad_tanlam_p_);

    for (int i=0;i<kNThr_;++i) {
        writeIf(eff_den_tanL_[i]); writeIf(eff_num_tanL_[i]);
        writeIf(eff_den_p_[i]);    writeIf(eff_num_p_[i]);
    }

    writeIf(h_S_tanL_);  writeIf(h_B_tanL_);
    writeIf(h_S_p_);     writeIf(h_B_p_);
    writeIf(h_SoverSqrtB_tanL_); writeIf(h_SoverSqrtB_p_);

    writeIf(h_match_purity_); writeIf(h_match_holes_);
    writeIf(h_vz_vs_purity_); writeIf(h_tanL_vs_purity_);

    writeIf(g_match_completeness);
    writeIf(g_match_matchedLayers);
    writeIf(g_comp_vs_purity);
    writeIf(g_fakes_per_event);
    writeIf(g_duplicates_per_event);
    writeIf(g_missedLayers_per_event);
    writeIf(g_mixedHits_per_event);

    writeIf(h_B_overlap_tanL_);
    writeIf(h_B_overlap_p_);
    writeIf(h_B_zero_tanL_);
    writeIf(h_B_zero_p_);

    // TEfficiency outputs
    for (int i=0;i<kNThr_;++i) {
        const int thr = thrList_[i];
        if (eff_den_tanL_[i] && eff_num_tanL_[i] && TEfficiency::CheckConsistency(*eff_num_tanL_[i], *eff_den_tanL_[i])) {
            auto* eTL = new TEfficiency(*eff_num_tanL_[i], *eff_den_tanL_[i]);
            eTL->SetName(Form("Eff_tanL_ge%d",thr));
            eTL->SetTitle(Form("Tracking efficiency vs tan#lambda (>= %d layers)",thr));
            eTL->SetStatisticOption(TEfficiency::kFCP);
            eTL->Write();
        }
        if (eff_den_p_[i] && eff_num_p_[i] && TEfficiency::CheckConsistency(*eff_num_p_[i], *eff_den_p_[i])) {
            auto* eP = new TEfficiency(*eff_num_p_[i], *eff_den_p_[i]);
            eP->SetName(Form("Eff_p_ge%d",thr));
            eP->SetTitle(Form("Tracking efficiency vs p (>= %d layers)",thr));
            eP->SetStatisticOption(TEfficiency::kFCP);
            eP->Write();
        }
    }

    outF_->Close();
}

DECLARE_PROCESSOR(RecoilProcessor);

