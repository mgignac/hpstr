// ============================================================================
// File: src/hpstr/processors/src/RecoilProcessor.cxx
// Cleaned: DEN by truth-hit counts; NUM by ID match + reco distinct-layer thresholds.
// S/B per-threshold uses identical gates. Robust branch resolution & diagnostics.
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

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TString.h"

namespace {
    // Minimal helper to resolve first existing branch from candidates.
    static std::string resolve_branch(TTree* t, const std::vector<std::string>& candidates) {
        if (!t) return {};
        for (const auto& c : candidates) {
            if (t->GetBranch(c.c_str())) return c;
        }
        return {};
    }
}

RecoilProcessor::RecoilProcessor(const std::string& name, Process& process)
    : Processor(name, process) {}

RecoilProcessor::~RecoilProcessor() {}

void RecoilProcessor::configure(const ParameterSet& parameters) {
    try {
        debug_     = parameters.getInteger("debug", debug_);
        anaName_   = parameters.getString("anaName", anaName_);
        beamE_     = parameters.getDouble("beamE", beamE_);
        sampleID_  = parameters.getInteger("sampleID", sampleID_);
        isData     = parameters.getInteger("isData", isData);

        // Optional overrides
        sclusColl_ = parameters.getString("sclusColl", sclusColl_);
        truthTracksCollRoot_ = parameters.getString("truthTrackCollRoot", truthTracksCollRoot_);
        trkColl_   = parameters.getString("trkCollRoot", trkColl_); // allow configure to override default "GBLTracks"
    } catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    }
}

void RecoilProcessor::initialize(TTree* tree) {
    tree_ = tree;

    // Bind must-haves
    evth_     = nullptr;         bevth_     = nullptr;
    mcParts_  = nullptr;         bmcParts_  = nullptr;
    trks_     = nullptr;         btrks_     = nullptr;
    mcTrkrHits_ = nullptr;       bmcTrkrHits_ = nullptr;

    // EventHeader
    tree_->SetBranchAddress("EventHeader", &evth_, &bevth_);

    // MCParticle (required)
    tree_->SetBranchAddress("MCParticle", &mcParts_, &bmcParts_);

    // Tracks (resolve)
    {
        std::vector<std::string> trkCandidates;
        if (!trkColl_.empty()) trkCandidates.push_back(trkColl_);
        trkCandidates.push_back("KalmanFullTracks");
        trkCandidates.push_back("GBLTracks");
        trkCandidates.push_back("Tracks");
        const std::string trkName = resolve_branch(tree_, trkCandidates);
        if (trkName.empty()) {
            std::cerr << "[RecoilProcessor] FATAL: no track branch found.\n";
        } else {
            if (debug_>0) std::cout << "[RecoilProcessor] Using track branch: " << trkName << "\n";
            tree_->SetBranchAddress(trkName.c_str(), &trks_, &btrks_);
        }
    }

    // SimTrackerHits (resolve) — needed for DEN by truth HITS
    {
        const std::vector<std::string> simHitCandidates{
            "TrackerSimHits","MCTrackerHits","SvtSimHits","TrackerHits"
        };
        const std::string simHitName = resolve_branch(tree_, simHitCandidates);
        if (simHitName.empty()) {
            mcTrkrHits_ = nullptr; bmcTrkrHits_ = nullptr;
            std::cerr << "[RecoilProcessor] WARNING: no sim-hit branch found; DEN by truth hits will be empty.\n";
        } else {
            tree_->SetBranchAddress(simHitName.c_str(), &mcTrkrHits_, &bmcTrkrHits_);
            if (debug_>0) std::cout << "[RecoilProcessor] Using sim-hit branch: " << simHitName << "\n";
        }
    }

    // Optional clusters (not required for logic)
    {
        std::vector<std::string> clusterCandidates;
        if (!sclusColl_.empty()) clusterCandidates.push_back(sclusColl_);
        clusterCandidates.push_back("SiClustersOnTrack");
        clusterCandidates.push_back("SiClustersOnTrack_KF");
        clusterCandidates.push_back("SvtClustersOnTrack");
        clusterCandidates.push_back("SvtClustersOnTrack_KF");
        clusterCandidates.push_back("SiClusters");
        clusterCandidates.push_back("SvtClusters");
        const std::string sclusName = resolve_branch(tree_, clusterCandidates);
        if (!sclusName.empty()) {
            tree_->SetBranchAddress(sclusName.c_str(), &svtClusters_, &bsvtClusters_);
            if (debug_>0) std::cout << "[RecoilProcessor] Using cluster branch: " << sclusName << "\n";
        }
    }

    // Output file + histos
    outF_ = TFile::Open("recoil_extra.root", "RECREATE");
    if (!outF_ || outF_->IsZombie()) {
        std::cerr << "[RecoilProcessor] ERROR opening output file.\n";
        outF_ = nullptr;
        return;
    }
    outF_->cd();

    const double maxP = (beamE_>0 ? beamE_ : 3.74);
    const int tanLBins = 100;  const double tanLMin=-0.3, tanLMax=0.3;
    const int pBins    = 100;  const double pMin=0.0, pMax=maxP;

    h_rec_p_        = new TH1F("Recoil_Momentum",    ";p (GeV);Counts", 100, 0., maxP);
    h_rec_theta_    = new TH1F("Recoil_Theta",       ";#theta (deg);Counts", 180, 0., 180.);
    h_rec_vtz_      = new TH1F("Recoil_VertexZ",     ";z_{vtx} (mm);Counts", 200, -20., 180.);
    h_rec_tanlam_p_ = new TH2F("Recoil_TanLambda_vs_P",";tan#lambda;p (GeV);Counts",100, -0.3, 0.3, 100, 0., maxP);
    h_rec_nhits_    = new TH1F("Recoil_nTruthHits",  ";N truth hits;Counts", 40, -0.5, 39.5);
    h_rec_tanlam_vs_nhits_ = new TH2F("Recoil_TanLambda_vs_nTruthHits",";tan#lambda;N Truth Hits; Counts",100, -0.3, 0.3,40, -0.5, 39.5);

    h_rad_p_        = new TH1F("Radiated_Momentum",  ";p (GeV);Counts", 100, 0., maxP);
    h_rad_theta_    = new TH1F("Radiated_Theta",     ";#theta (deg);Counts", 180, 0., 180.);
    h_rad_vtz_      = new TH1F("Radiated_VertexZ",   ";z_{vtx} (mm);Counts", 200, -20., 180.);
    h_rad_tanlam_p_ = new TH2F("Radiated_TanLambda_vs_P",";tan#lambda;p (GeV);Counts",100, -0.3, 0.3, 100, 0., maxP);

    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        eff_den_tanL_[i] = new TH1F(Form("EffDen_tanL_ge%d",thr), Form("Den: truth hits >= %d;tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        eff_num_tanL_[i] = new TH1F(Form("EffNum_tanL_ge%d",thr), Form("Num: matched & reco layers >= %d;tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        eff_den_p_[i]    = new TH1F(Form("EffDen_p_ge%d",thr),    Form("Den: truth hits >= %d;p (GeV);count",thr), pBins, pMin, pMax);
        eff_num_p_[i]    = new TH1F(Form("EffNum_p_ge%d",thr),    Form("Num: matched & reco layers >= %d;p (GeV);count",thr), pBins, pMin, pMax);
    }

    h_S_tanL_ = new TH1F("S_tanL", "Matched (ID in recoils);tan#lambda;S", tanLBins, tanLMin, tanLMax);
    h_B_tanL_ = new TH1F("B_tanL", "Unmatched;tan#lambda;B", tanLBins, tanLMin, tanLMax);
    h_S_p_    = new TH1F("S_p",    "Matched (ID in recoils);p (GeV);S",    pBins, pMin, pMax);
    h_B_p_    = new TH1F("B_p",    "Unmatched;p (GeV);B",                  pBins, pMin, pMax);
    h_SoverSqrtB_tanL_ = nullptr;
    h_SoverSqrtB_p_    = nullptr;

    // Per-threshold S/B
    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        h_S_tanL_thr_[i] = new TH1F(Form("S_tanL_ge%d",thr),  Form("S (truth hits & reco layers >= %d);tan#lambda;S",thr), tanLBins, tanLMin, tanLMax);
        h_B_tanL_thr_[i] = new TH1F(Form("B_tanL_ge%d",thr),  Form("B (else) >= %d;tan#lambda;B",thr), tanLBins, tanLMin, tanLMax);
        h_S_p_thr_[i]    = new TH1F(Form("S_p_ge%d",thr),     Form("S (truth hits & reco layers >= %d);p (GeV);S",thr),    pBins, pMin, pMax);
        h_B_p_thr_[i]    = new TH1F(Form("B_p_ge%d",thr),     Form("B (else) >= %d;p (GeV);B",thr),    pBins, pMin, pMax);
        h_SoverSqrtB_tanL_thr_[i] = nullptr;
        h_SoverSqrtB_p_thr_[i]    = nullptr;
    }

    // Matching-quality (kept for compatibility; not heavily used here)
    h_match_purity_   = new TH1F("MatchPurity", "purity;purity;tracks", 50, 0.0, 1.0);
    h_match_holes_    = new TH1F("MatchHoles",  "holes;holes;tracks", kNLayersMax_+1, -0.5, kNLayersMax_+0.5);
    h_vz_vs_purity_   = new TH2F("VZ_vs_Purity",";v_{z} (mm);purity", 200, -20., 180., 50, 0.0, 1.0);
    h_tanL_vs_purity_ = new TH2F("TanL_vs_Purity",";tan#lambda;purity", tanLBins, tanLMin, tanLMax, 50, 0.0, 1.0);

    h_B_overlap_tanL_ = new TH1F("Boverlap_tanL","Unmatched with overlap (not used here);tan#lambda;B", tanLBins, tanLMin, tanLMax);
    h_B_overlap_p_    = new TH1F("Boverlap_p",   "Unmatched with overlap (not used here);p (GeV);B", pBins, pMin, pMax);
    h_B_zero_tanL_    = new TH1F("Bzero_tanL",   "Unmatched pure fake (not used here);tan#lambda;B", tanLBins, tanLMin, tanLMax);
    h_B_zero_p_       = new TH1F("Bzero_p",      "Unmatched pure fake (not used here);p (GeV);B", pBins, pMin, pMax);

    // Tree
    recoilTree_ = new TTree("recoilTree","Truth-level recoil info");
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

bool RecoilProcessor::process(IEvent* /*ievent*/) {
    if (!evth_ || !mcParts_ || !trks_) {
        std::cerr << "[RecoilProcessor] Missing required branches (EventHeader/MCParticle/Tracks)\n";
        return false;
    }

    // === Select recoils (PDG filter). Change 13->11 if needed for electrons. ===
    const double ETOL  = 1e-3;
    const double BEAME = (beamE_>0 ? beamE_ : 3.74);
    std::vector<MCParticle*> recoils;
    recoils.reserve(mcParts_->size());
    for (auto* part : *mcParts_) {
        if (!part) continue;
        if (std::abs(part->getPDG()) != 13) continue;      // <- swap to 11 for electrons if desired
        if (std::abs(part->getEnergy() - BEAME) < ETOL) continue;
        recoils.push_back(part);
    }
    if (recoils.empty()) return true;


    // -- (1) Map pid -> MCParticle*
std::unordered_map<int, MCParticle*> id2part;
id2part.reserve(mcParts_->size());
for (auto* p : *mcParts_) { if (p) id2part[p->getID()] = p; }

// -- (2) nhit_denom: MCParticle* -> # MC tracker hits (truth hits)
std::unordered_map<MCParticle*, int> nHits;
if (mcTrkrHits_) {
    nHits.reserve(id2part.size());
    for (auto* hit : *mcTrkrHits_) {
        if (!hit) continue;
        const int pid = hit->getID();                // MCParticle ID on sim-hit
        auto it = id2part.find(pid);
        if (it != id2part.end()) nHits[it->second]++; // pointer-keyed (fast recoil lookup)
    }
} else if (debug_ > 0) {
    std::cout << "[RecoilProcessor] WARNING: mcTrkrHits_ is null; nHits[e]=0 for all.\n";
}

// Helper: denominator nhits for a given MCParticle*
auto nhitDenom = [&](MCParticle* e)->int {
    auto it = nHits.find(e);
    return (it != nHits.end()) ? it->second : 0;
};

// -- (3) Best distinct reco-layer count per Track::getID() (unchanged NUM logic)
std::unordered_map<int,int> maxRecoLayersByTrkID;
maxRecoLayersByTrkID.reserve(trks_->size());
for (auto* trk : *trks_) {
    if (!trk) continue;
    std::array<bool,64> mask{}; mask.fill(false);
    const int Nh = trk->getSvtHits().GetEntriesFast();
    for (int ih=0; ih<Nh; ++ih) {
        auto* th = static_cast<::TrackerHit*>(trk->getSvtHits().At(ih));
        if (!th) continue;
        const int L = th->getLayer();
        if (0 <= L && L < kNLayersMax_) mask[L] = true;
    }
    int nL=0; for (int L=0; L<kNLayersMax_; ++L) if (mask[L]) ++nL;
    const int tid = trk->getID();
    auto it = maxRecoLayersByTrkID.find(tid);
    if (it==maxRecoLayersByTrkID.end() || nL>it->second) maxRecoLayersByTrkID[tid] = nL;
}

// -- (4) Recoil ID set (for S/B and quick membership)
std::unordered_set<int> recoilIDs;
recoilIDs.reserve(recoils.size());
for (auto* r : recoils) if (r) recoilIDs.insert(r->getID());

// -- (5) Efficiency DEN/NUM (DEN: nhit_denom only; NUM: reco-layers only)
for (auto* e : recoils) {
    if (!e) continue;

    const auto& mom = e->getMomentum();
    const double tl = (mom[2] != 0.0) ? (mom[1] / mom[2]) : 999.;   // keep efficiency x-axis continuity
    const double p  = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);

    const int pid        = e->getID();
    const int nh         = nhitDenom(e);                            // <-- DEN gate
    const int bestRecoL  = (maxRecoLayersByTrkID.count(pid) ? maxRecoLayersByTrkID[pid] : 0);

    // DEN: nhit_denom thresholds ONLY
    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        if (nh >= thr) {
            eff_den_tanL_[i]->Fill(tl);
            eff_den_p_[i]->Fill(p);
        }
    }

    // NUM: ID + best distinct reco-layers thresholds ONLY (decoupled from DEN)
    if (bestRecoL > 0) {
        for (int i = 0; i < kNThr_; ++i) {
            const int thr = thrList_[i];
            if (bestRecoL >= thr) {
                eff_num_tanL_[i]->Fill(tl);
                eff_num_p_[i]->Fill(p);
            }
        }
    }
}

// -- (6) Global S/B (unchanged: membership only)
for (auto* trk : *trks_) {
    if (!trk) continue;
    const bool isSig = (recoilIDs.count(trk->getID()) > 0);
    const double tl = trk->getTanLambda();
    const double pp = trk->getP();
    if (isSig) { h_S_tanL_->Fill(tl); h_S_p_->Fill(pp); }
    else       { h_B_tanL_->Fill(tl); h_B_p_->Fill(pp); }
}

// -- (7) Per-threshold S/B (by reco-layer thresholds only; mirrors NUM)
for (auto* trk : *trks_) {
    if (!trk) continue;
    const int tid = trk->getID();
    const bool isRecoil = (recoilIDs.count(tid) > 0);
    const int bestRecoL = (maxRecoLayersByTrkID.count(tid) ? maxRecoLayersByTrkID[tid] : 0);
    const double tl = trk->getTanLambda();
    const double pp = trk->getP();

    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        if (isRecoil && bestRecoL >= thr) {
            h_S_tanL_thr_[i]->Fill(tl);  h_S_p_thr_[i]->Fill(pp);
        } else {
            h_B_tanL_thr_[i]->Fill(tl);  h_B_p_thr_[i]->Fill(pp);
        }
    }
}

// -- (8) Diagnostic denominator plots exactly like your snippet
for (auto* e : recoils) {
    if (!e) continue;
    const int nh = nhitDenom(e);
    h_rec_nhits_->Fill(nh);
    const auto& mom = e->getMomentum();
    const double pt   = std::hypot(mom[0], mom[1]);
    const double tanL = (mom[2]!=0. ? pt/mom[2] : 999.);
    h_rec_tanlam_vs_nhits_->Fill(tanL, nh);
}

// -- (9) Extra debug
if (debug_ > 0) {
    int shown = 0;
    std::cout << "[RecoilProcessor] DEN nhit_denom summary (first few recoils): ";
    for (auto* e : recoils) {
        if (!e) continue;
        std::cout << e->getID() << ":" << nhitDenom(e) << " ";
        if (++shown == 8) break;
    }
    std::cout << "\n";
}










// -- (G) Tree: store nhit_denom exactly as requested --------------------------
for (auto* e : recoils) {
    if (!e) continue;

    tree_sampleID_  = sampleID_;
    tree_category_  = 0;
    tree_pdg_       = e->getPDG();
    tree_OriginPDG_ = e->getOriginPDG();
    tree_MomPDG_    = e->getMomPDG();
    tree_charge_    = e->getCharge();
    tree_ID_        = e->getID();

    auto mom = e->getMomentum();
    tree_px_ = mom[0]; tree_py_ = mom[1]; tree_pz_ = mom[2];
    tree_p_  = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
    tree_pt_ = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1]);
    tree_tanLambda_ = (tree_pz_!=0 ? tree_py_/tree_pz_ : 999.);
    tree_theta_     = std::atan2(tree_pt_,tree_pz_) * TMath::RadToDeg();
    tree_energy_    = e->getEnergy();
    tree_mass_      = e->getMass();
    tree_time_      = e->getTime();
    const auto& vtx = e->getVertexPosition();
    tree_vz_        = (vtx.size()>2) ? vtx[2] : 0.0f;

    const int nh = nhitDenom(e);             // <- nhit_denom stored in tree
    tree_nHits_   = nh;
    tree_eventID_ = evth_->getEventNumber();
    tree_x_       = 1.0f - float(e->getEnergy() / beamE_);
    tree_phi_     = std::atan2(tree_py_, tree_px_);
    recoilTree_->Fill();

    const double tanL = (mom[2] != 0.0) ? (mom[1]/mom[2]) : 999.;
    h_rec_p_->Fill(tree_p_);
    h_rec_theta_->Fill(tree_theta_);
    h_rec_vtz_->Fill(tree_vz_);
    h_rec_tanlam_p_->Fill(tanL, tree_p_);
    h_rec_nhits_->Fill(nh);
    if (nh > 0) h_rec_tanlam_vs_nhits_->Fill(tanL, nh);
}

// -- (H) Debug ---------------------------------------------------------------
if (debug_ > 0) {
    std::cout << "[RecoilProcessor][nhit_denom] Event " << evth_->getEventNumber()
              << " recoils=" << recoils.size()
              << " trks="    << (trks_ ? trks_->size() : 0)
              << " mcTrkrHits=" << (mcTrkrHits_ ? mcTrkrHits_->size() : 0)
              << " nHitsMap=" << nHits.size()
              << "\n";
}
// ============================================================================
// End patch






//// === DROP-IN: place inside RecoilProcessor::process(), after you build `recoils` ===
//// Assumes: trks_, mcTrkrHits_, eff_den_*[], eff_num_*[], h_S_*[], h_B_*[], h_*_thr_[],
////          thrList_[kNThr_], kNLayersMax_, and your existing histos are already defined.
//
///* 0) Fast bailouts */
////if (!trks_ || recoils.empty()) return true;
//
///* 1) Fast ID set for simple matching */
//std::unordered_set<int> trkIDs;
//trkIDs.reserve(trks_->size());
//for (auto* trk : *trks_) if (trk) trkIDs.insert(trk->getID());
//
///* 2) Best distinct reco-layer count per Track::getID() (handles duplicates) */
//std::unordered_map<int,int> maxRecoLayersByTrkID;
//maxRecoLayersByTrkID.reserve(trks_->size());
//for (auto* trk : *trks_) {
//    if (!trk) continue;
//    std::array<bool,64> mask{}; mask.fill(false);
//    const int Nh = trk->getSvtHits().GetEntriesFast();
//    for (int ih=0; ih<Nh; ++ih) {
//        auto* th = static_cast<::TrackerHit*>(trk->getSvtHits().At(ih));
//        if (!th) continue;
//        const int L = th->getLayer();
//        if (0 <= L && L < kNLayersMax_) mask[L] = true;
//    }
//    int nL = 0; for (int L=0; L<kNLayersMax_; ++L) if (mask[L]) ++nL;
//    const int tid = trk->getID();
//    auto it = maxRecoLayersByTrkID.find(tid);
//    if (it==maxRecoLayersByTrkID.end() || nL>it->second) maxRecoLayersByTrkID[tid] = nL;
//}
//
///* 3) Truth *hit* counts per MC pid for DEN (safe fallback if missing) */
//std::unordered_map<int,int> truthHitCountByPid;
//bool haveTruthHits = false;
//if (mcTrkrHits_) {
//    truthHitCountByPid.reserve(4096);
//    for (auto* hit : *mcTrkrHits_) {
//        if (!hit) continue;
//        const int pid = hit->getID();      // expected: MCParticle ID
//        truthHitCountByPid[pid] += 1;
//    }
//    haveTruthHits = !truthHitCountByPid.empty();
//}
//// Gate helper: if no sim-hits, DEN gate returns true so histos still fill.
//auto passTruthThr = [&](int pid, int thr)->bool {
//    if (!haveTruthHits) return true;              // fallback: always in DEN
//    auto it = truthHitCountByPid.find(pid);
//    const int n = (it!=truthHitCountByPid.end() ? it->second : 0);
//    return n >= thr;
//};
//auto truthHits = [&](int pid)->int {
//    auto it = truthHitCountByPid.find(pid);
//    return (it!=truthHitCountByPid.end() ? it->second : 0);
//};
//
///* 4) Recoil MC ID set (for S/B) */
//std::unordered_set<int> recoilIDs;
//recoilIDs.reserve(recoils.size());
//for (auto* r : recoils) if (r) recoilIDs.insert(r->getID());
//
///* 5) Efficiency DEN/NUM with thresholds (DEN by truth hits, NUM by reco layers) */
//for (auto* part : recoils) {
//    if (!part) continue;
//
//    const auto& mom = part->getMomentum();
//    const double tl = (mom[2] != 0.0) ? (mom[1] / mom[2]) : 999.;
//    const double p  = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
//
//    const int pid = part->getID();
//    const bool matched = (trkIDs.find(pid) != trkIDs.end());
//    const int nRecoLay = (maxRecoLayersByTrkID.count(pid) ? maxRecoLayersByTrkID[pid] : 0);
//
//    for (int i = 0; i < kNThr_; ++i) {
//        const int thr = thrList_[i];
//
//        if (passTruthThr(pid, thr)) {               // DEN gate (truth hits or fallback)
//            eff_den_tanL_[i]->Fill(tl);
//            eff_den_p_[i]->Fill(p);
//
//            if (matched && nRecoLay >= thr) {       // NUM gate (ID match + reco layers)
//                eff_num_tanL_[i]->Fill(tl);
//                eff_num_p_[i]->Fill(p);
//            }
//        }
//    }
//}
//
///* 6) Global S/B (ID membership only; unchanged semantics) */
//for (auto* trk : *trks_) {
//    if (!trk) continue;
//    const bool isSig = (recoilIDs.count(trk->getID()) > 0);
//    const double tl = trk->getTanLambda();
//    const double pp = trk->getP();
//    if (isSig) { h_S_tanL_->Fill(tl); h_S_p_->Fill(pp); }
//    else       { h_B_tanL_->Fill(tl); h_B_p_->Fill(pp); }
//}
//
///* 7) Per-threshold S/B (same gates as NUM) */
//for (auto* trk : *trks_) {
//    if (!trk) continue;
//    const int tid = trk->getID();
//    const bool isRecoil = (recoilIDs.count(tid) > 0);
//    const int nRecoLay  = (maxRecoLayersByTrkID.count(tid) ? maxRecoLayersByTrkID[tid] : 0);
//    const double tl = trk->getTanLambda();
//    const double pp = trk->getP();
//
//    for (int i = 0; i < kNThr_; ++i) {
//        const int thr = thrList_[i];
//        if (isRecoil && passTruthThr(tid, thr) && nRecoLay >= thr) {
//            h_S_tanL_thr_[i]->Fill(tl);  h_S_p_thr_[i]->Fill(pp);
//        } else {
//            h_B_tanL_thr_[i]->Fill(tl);  h_B_p_thr_[i]->Fill(pp);
//        }
//    }
//}
//
//
//auto bestRecoLayersFor = [&](int pid)->int {
//    auto it = maxRecoLayersByTrkID.find(pid);
//    return (it != maxRecoLayersByTrkID.end()) ? it->second : 0;
//};



    //// === Distinct reco-layer count per Track::getID() (best over dups) ===
    //std::unordered_map<int,int> maxRecoLayersByTrkID;
    //maxRecoLayersByTrkID.reserve(trks_->size());
    //for (auto* trk : *trks_) {
    //    if (!trk) continue;
    //    std::array<bool,64> mask{}; mask.fill(false);
    //    const int Nh = trk->getSvtHits().GetEntriesFast();
    //    for (int ih=0; ih<Nh; ++ih) {
    //        auto* th = static_cast<::TrackerHit*>(trk->getSvtHits().At(ih));
    //        if (!th) continue;
    //        const int L = th->getLayer();
    //        if (0 <= L && L < kNLayersMax_) mask[L] = true;
    //    }
    //    int nL=0; for (int L=0; L<kNLayersMax_; ++L) if (mask[L]) ++nL;
    //    const int tid = trk->getID();                   // set by TrackingProcessor to truth ID
    //    auto it = maxRecoLayersByTrkID.find(tid);
    //    if (it==maxRecoLayersByTrkID.end() || nL>it->second) maxRecoLayersByTrkID[tid]=nL;
    //}

    //// === Truth *hit* counts per pid (DEN) ===
    //std::unordered_map<int,int> truthHitCountByPid;
    //if (mcTrkrHits_) {
    //    truthHitCountByPid.reserve(4096);
    //    for (auto* hit : *mcTrkrHits_) {
    //        if (!hit) continue;
    //        const int pid = hit->getID();               // MCParticle ID carried on hit
    //        truthHitCountByPid[pid] += 1;
    //    }
    //} else if (debug_>0) {
    //    std::cout << "[RecoilProcessor] No sim-hits bound; denominators will be zero.\n";
    //}
    //auto nTruthHits = [&](int pid)->int {
    //    auto it = truthHitCountByPid.find(pid);
    //    return (it != truthHitCountByPid.end()) ? it->second : 0;
    //};

    //// === Set of recoil IDs (for S/B) ===
    //std::unordered_set<int> recoilIDs;
    //recoilIDs.reserve(recoils.size());
    //for (auto* r : recoils) if (r) recoilIDs.insert(r->getID());

    //// === Efficiency DEN/NUM (DEN by truth-hits; NUM by reco-layer thresholds) ===
    //for (auto* part : recoils) {
    //    if (!part) continue;

    //    const auto& mom = part->getMomentum();
    //    const double tl = (mom[2] != 0.0) ? (mom[1] / mom[2]) : 999.;
    //    const double p  = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);

    //    const int pid = part->getID();
    //    const int hitsTruth = nTruthHits(pid);
    //    const int recoLay   = (maxRecoLayersByTrkID.count(pid) ? maxRecoLayersByTrkID[pid] : 0);

    //    for (int i = 0; i < kNThr_; ++i) {
    //        const int thr = thrList_[i];
    //        if (hitsTruth >= thr) {
    //            eff_den_tanL_[i]->Fill(tl);
    //            eff_den_p_[i]->Fill(p);
    //            if (recoLay >= thr) {
    //                eff_num_tanL_[i]->Fill(tl);
    //                eff_num_p_[i]->Fill(p);
    //            }
    //        }
    //    }
    //}

    //// === Global S/B (ID membership only) ===
    //for (auto* trk : *trks_) {
    //    if (!trk) continue;
    //    const bool isSig = (recoilIDs.count(trk->getID()) > 0);
    //    const double tl = trk->getTanLambda();
    //    const double pp = trk->getP();
    //    if (isSig) { h_S_tanL_->Fill(tl); h_S_p_->Fill(pp); }
    //    else       { h_B_tanL_->Fill(tl); h_B_p_->Fill(pp); }
    //}

    //// === Per-threshold S/B (truth hits + reco layers + ID) ===
    //for (auto* trk : *trks_) {
    //    if (!trk) continue;
    //    const int tid = trk->getID();
    //    const bool isRecoil = (recoilIDs.count(tid) > 0);
    //    const int hitsTruth = nTruthHits(tid);
    //    const int recoLay   = (maxRecoLayersByTrkID.count(tid) ? maxRecoLayersByTrkID[tid] : 0);
    //    const double tl = trk->getTanLambda();
    //    const double pp = trk->getP();

    //    for (int i = 0; i < kNThr_; ++i) {
    //        const int thr = thrList_[i];
    //        if (isRecoil && hitsTruth >= thr && recoLay >= thr) {
    //            h_S_tanL_thr_[i]->Fill(tl);  h_S_p_thr_[i]->Fill(pp);
    //        } else {
    //            h_B_tanL_thr_[i]->Fill(tl);  h_B_p_thr_[i]->Fill(pp);
    //        }
    //    }
    //}



//// === Tree & kinematics (use truth hits; fallback to reco-layer count) ===
//for (auto* e : recoils) {
//    if (!e) continue;
//
//    tree_sampleID_  = sampleID_;
//    tree_category_  = 0;
//    tree_pdg_       = e->getPDG();
//    tree_OriginPDG_ = e->getOriginPDG();
//    tree_MomPDG_    = e->getMomPDG();
//    tree_charge_    = e->getCharge();
//    tree_ID_        = e->getID();
//
//    auto mom = e->getMomentum();
//    tree_px_ = mom[0]; tree_py_ = mom[1]; tree_pz_ = mom[2];
//    tree_p_  = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
//    tree_pt_ = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1]);
//    tree_tanLambda_ = (tree_pz_!=0 ? tree_py_/tree_pz_ : 999.);
//    tree_theta_     = std::atan2(tree_pt_,tree_pz_) * TMath::RadToDeg();
//    tree_energy_    = e->getEnergy();
//    tree_mass_      = e->getMass();
//    tree_time_      = e->getTime();
//    const auto& vtx = e->getVertexPosition();
//    tree_vz_        = (vtx.size()>2) ? vtx[2] : 0.0f;
//
//    const int pid = e->getID();
//    const int nTruth = truthHits(pid);
//    const int nRecoL = bestRecoLayersFor(pid);
//    const int nh     = (haveTruthHits ? nTruth : nRecoL);   // fallback so it always fills
//
//    tree_nHits_     = nh;                                   // stored as "nHits" in tree
//    tree_eventID_   = evth_->getEventNumber();
//    tree_x_         = 1.0f - float(e->getEnergy() / beamE_);
//    tree_phi_       = std::atan2(tree_py_, tree_px_);
//    recoilTree_->Fill();
//
//    const double tanL = (mom[2] != 0.0) ? (mom[1]/mom[2]) : 999.;
//    h_rec_p_->Fill(tree_p_);
//    h_rec_theta_->Fill(tree_theta_);
//    h_rec_vtz_->Fill(tree_vz_);
//    h_rec_tanlam_p_->Fill(tanL, tree_p_);
//    h_rec_nhits_->Fill(nh);
//    if (nh > 0) h_rec_tanlam_vs_nhits_->Fill(tanL, nh);
//}
















    return true;
}

void RecoilProcessor::finalize() {
    if (!outF_ || !outF_->IsOpen()) {
        std::cerr << "[RecoilProcessor] finalize(): output not open\n";
        return;
    }
    outF_->cd();
    auto writeIf = [](TObject* o){ if (o) o->Write(); };

    // Global S/√B
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

    // Per-threshold S/√B
    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        if (h_S_tanL_thr_[i] && h_B_tanL_thr_[i]) {
            h_SoverSqrtB_tanL_thr_[i] = (TH1F*)h_S_tanL_thr_[i]->Clone(Form("SoverSqrtB_tanL_ge%d",thr));
            h_SoverSqrtB_tanL_thr_[i]->SetTitle(Form("S/#sqrt{B} vs tan#lambda (>= %d)",thr));
            for (int b=1; b<=h_SoverSqrtB_tanL_thr_[i]->GetNbinsX(); ++b) {
                const double s = h_S_tanL_thr_[i]->GetBinContent(b);
                const double bg= h_B_tanL_thr_[i]->GetBinContent(b);
                h_SoverSqrtB_tanL_thr_[i]->SetBinContent(b, (bg>0.0)? s/std::sqrt(bg) : 0.0);
                h_SoverSqrtB_tanL_thr_[i]->SetBinError(b, 0.0);
            }
            h_SoverSqrtB_tanL_thr_[i]->Write();
        }
        if (h_S_p_thr_[i] && h_B_p_thr_[i]) {
            h_SoverSqrtB_p_thr_[i] = (TH1F*)h_S_p_thr_[i]->Clone(Form("SoverSqrtB_p_ge%d",thr));
            h_SoverSqrtB_p_thr_[i]->SetTitle(Form("S/#sqrt{B} vs p (>= %d)",thr));
            for (int b=1; b<=h_SoverSqrtB_p_thr_[i]->GetNbinsX(); ++b) {
                const double s = h_S_p_thr_[i]->GetBinContent(b);
                const double bg= h_B_p_thr_[i]->GetBinContent(b);
                h_SoverSqrtB_p_thr_[i]->SetBinContent(b, (bg>0.0)? s/std::sqrt(bg) : 0.0);
                h_SoverSqrtB_p_thr_[i]->SetBinError(b, 0.0);
            }
            h_SoverSqrtB_p_thr_[i]->Write();
        }
    }

    // Write everything
    writeIf(recoilTree_);
    writeIf(h_rec_p_); writeIf(h_rec_theta_); writeIf(h_rec_vtz_);
    writeIf(h_rec_tanlam_p_); writeIf(h_rec_nhits_); writeIf(h_rec_tanlam_vs_nhits_);
    writeIf(h_rad_p_); writeIf(h_rad_theta_); writeIf(h_rad_vtz_); writeIf(h_rad_tanlam_p_);

    for (int i=0;i<kNThr_;++i) {
        writeIf(eff_den_tanL_[i]); writeIf(eff_num_tanL_[i]);
        writeIf(eff_den_p_[i]);    writeIf(eff_num_p_[i]);
        writeIf(h_S_tanL_thr_[i]); writeIf(h_B_tanL_thr_[i]);
        writeIf(h_S_p_thr_[i]);    writeIf(h_B_p_thr_[i]);
    }

    writeIf(h_S_tanL_);  writeIf(h_B_tanL_);
    writeIf(h_S_p_);     writeIf(h_B_p_);
    writeIf(h_SoverSqrtB_tanL_); writeIf(h_SoverSqrtB_p_);

    writeIf(h_match_purity_); writeIf(h_match_holes_);
    writeIf(h_vz_vs_purity_); writeIf(h_tanL_vs_purity_);

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
            eTL->SetTitle(Form("Efficiency vs tan#lambda (>= %d)",thr));
            eTL->SetStatisticOption(TEfficiency::kFCP);
            eTL->Write();
        }
        if (eff_den_p_[i] && eff_num_p_[i] && TEfficiency::CheckConsistency(*eff_num_p_[i], *eff_den_p_[i])) {
            auto* eP = new TEfficiency(*eff_num_p_[i], *eff_den_p_[i]);
            eP->SetName(Form("Eff_p_ge%d",thr));
            eP->SetTitle(Form("Efficiency vs p (>= %d)",thr));
            eP->SetStatisticOption(TEfficiency::kFCP);
            eP->Write();
        }
    }

    outF_->Close();
}

DECLARE_PROCESSOR(RecoilProcessor);

