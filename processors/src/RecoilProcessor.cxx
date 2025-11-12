// =========================== RecoilProcessor.cxx  ==========================
#include "RecoilProcessor.h"

#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TMath.h"

namespace {
    static std::string resolve_branch(TTree* t, const std::vector<std::string>& candidates) {
        if (!t) return {};
        for (const auto& c : candidates) if (t->GetBranch(c.c_str())) return c;
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

        trkColl_   = parameters.getString("trkCollRoot", trkColl_);
        sclusColl_ = parameters.getString("sclusColl",   sclusColl_);
        truthTracksCollRoot_ = parameters.getString("truthTrackCollRoot", truthTracksCollRoot_);

    } catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    }
}

void RecoilProcessor::initialize(TTree* tree) {
    tree_ = tree;

    // Branch ptrs
    evth_        = nullptr; bevth_        = nullptr;
    mcParts_     = nullptr; bmcParts_     = nullptr;
    trks_        = nullptr; btrks_        = nullptr;
    svtClusters_ = nullptr; bsvtClusters_ = nullptr;

    // Required
    tree_->SetBranchAddress("EventHeader", &evth_, &bevth_);
    tree_->SetBranchAddress("MCParticle",  &mcParts_, &bmcParts_);

    // Tracks branch
    {
        std::vector<std::string> trkCandidates;
        if (!trkColl_.empty()) trkCandidates.push_back(trkColl_);
        trkCandidates.push_back("KalmanFullTracks");
        trkCandidates.push_back("GBLTracks");
        trkCandidates.push_back("Tracks");
        const auto trkName = resolve_branch(tree_, trkCandidates);
        if (trkName.empty()) {
            std::cerr << "[RecoilProcessor] FATAL: no track branch found.\n";
        } else {
            if (debug_>0) std::cout << "[RecoilProcessor] Using track branch: " << trkName << "\n";
            tree_->SetBranchAddress(trkName.c_str(), &trks_, &btrks_);
        }
        if (!truthTracksCollRoot_.empty() && tree_->GetBranch(truthTracksCollRoot_.c_str())) {
            tree_->SetBranchAddress(truthTracksCollRoot_.c_str(), &truthTrks_, &btruthTrks_);
            if (debug_>0) std::cout << "[RecoilProcessor] Using truth-track branch: " << truthTracksCollRoot_ << "\n";
        }
    }

    // SiClusters branch (prefer unbiased)
    {
        std::vector<std::string> clusterCandidates;

        if (!sclusColl_.empty() && sclusColl_.find("OnTrack") == std::string::npos)
            clusterCandidates.push_back(sclusColl_);

        clusterCandidates.push_back("SiClusters");
        clusterCandidates.push_back("SvtClusters");
        clusterCandidates.push_back("RotatedHelicalTrackHits");
        clusterCandidates.push_back("Tracker3DHits");

        if (!sclusColl_.empty() && sclusColl_.find("OnTrack") != std::string::npos)
            clusterCandidates.push_back(sclusColl_);
        clusterCandidates.push_back("SiClustersOnTrack");
        clusterCandidates.push_back("SvtClustersOnTrack");

        const auto sclusName = resolve_branch(tree_, clusterCandidates);
        if (sclusName.empty()) {
            std::cerr << "[RecoilProcessor] FATAL: no Si/Tracker cluster branch found; cannot build findable set.\n";
        } else {
            if (debug_>0) {
                std::cout << "[RecoilProcessor] Using cluster branch: " << sclusName << "\n";
                if (sclusName.find("OnTrack") != std::string::npos)
                    std::cout << "[RecoilProcessor] WARNING: Using *OnTrack* clusters biases findable set toward reconstructed objects.\n";
            }
            tree_->SetBranchAddress(sclusName.c_str(), &svtClusters_, &bsvtClusters_);
        }

        if (!truthTracksCollRoot_.empty() && !tree_->GetBranch(truthTracksCollRoot_.c_str()) && debug_ > 0) {
            std::cout << "[RecoilProcessor] WARNING: truth-track branch '"
                      << truthTracksCollRoot_ << "' not found. Only ID+layers matching will be used.\n";
        }
    }

    // Output file
    outF_ = TFile::Open("recoil_extra.root", "RECREATE");
    if (!outF_ || outF_->IsZombie()) { std::cerr << "[RecoilProcessor] ERROR opening output file.\n"; outF_=nullptr; return; }
    outF_->cd();

    // Histograms
    const double maxP = (beamE_>0 ? beamE_ : 3.74);
    const int tanLBins = 100;  const double tanLMin=-0.3, tanLMax=0.3;
    const int pBins    = 100;  const double pMin=0.0, pMax=maxP;

    // Efficiency hists
    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        eff_den_tanL_[i] = new TH1F(Form("EffDen_tanL_ge%d",thr), Form("Efficiency DEN: findable (>= %d);tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        eff_num_tanL_[i] = new TH1F(Form("EffNum_tanL_ge%d",thr), Form("Efficiency NUM: matched & findable (>= %d);tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        eff_den_p_[i]    = new TH1F(Form("EffDen_p_ge%d",thr),    Form("Efficiency DEN: findable (>= %d);p (GeV);count",thr), pBins, pMin, pMax);
        eff_num_p_[i]    = new TH1F(Form("EffNum_p_ge%d",thr),    Form("Efficiency NUM: matched & findable (>= %d);p (GeV);count",thr), pBins, pMin, pMax);
    }

    // Acceptance hists
    for (int i = 0; i < kNThr_; ++i) {
        const int thr = thrList_[i];
        acc_den_tanL_[i] = new TH1F(Form("AccDen_tanL_ge%d",thr), Form("Acceptance DEN: all generated recoil;tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        acc_num_tanL_[i] = new TH1F(Form("AccNum_tanL_ge%d",thr), Form("Acceptance NUM: findable (>= %d);tan#lambda;count",thr), tanLBins, tanLMin, tanLMax);
        acc_den_p_[i]    = new TH1F(Form("AccDen_p_ge%d",thr),    Form("Acceptance DEN: all generated recoil;p (GeV);count",thr), pBins, pMin, pMax);
        acc_num_p_[i]    = new TH1F(Form("AccNum_p_ge%d",thr),    Form("Acceptance NUM: findable (>= %d);p (GeV);count",thr), pBins, pMin, pMax);
    }

    // Helpful match-quality histos already declared in header
    h_match_purity_        = new TH1F("match_purity",        "Best match purity;purity;count", 50, 0.0, 1.0);
    h_match_holes_         = new TH1F("match_holes",         "Missing truth layers on best match;holes;count", 20, 0.0, 20.0);
    h_vz_vs_purity_        = new TH2F("vz_vs_purity",        "v_{z} vs purity;v_{z};purity", 100, -50.0, 50.0, 50, 0.0, 1.0);
    h_tanL_vs_purity_      = new TH2F("tanL_vs_purity",      "tan#lambda vs purity;tan#lambda;purity", tanLBins, tanLMin, tanLMax, 50, 0.0, 1.0);
    h_B_overlap_tanL_      = new TH1F("B_overlap_tanL",      "Findable with >1 ID-matched tracks;tan#lambda;count", tanLBins, tanLMin, tanLMax);
    h_B_overlap_p_         = new TH1F("B_overlap_p",         "Findable with >1 ID-matched tracks;p (GeV);count", pBins, pMin, pMax);
    h_B_zero_tanL_         = new TH1F("B_zero_tanL",         "Findable with zero ID-matched tracks;tan#lambda;count", tanLBins, tanLMin, tanLMax);
    h_B_zero_p_            = new TH1F("B_zero_p",            "Findable with zero ID-matched tracks;p (GeV);count", pBins, pMin, pMax);

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
    recoilTree_->Branch("tanLambda",&tree_tanLambda_,"tanLambda/F");
    recoilTree_->Branch("theta",&tree_theta_,"theta/F");
    recoilTree_->Branch("energy",&tree_energy_,"energy/F"); recoilTree_->Branch("mass",&tree_mass_,"mass/F");
    recoilTree_->Branch("time",&tree_time_,"time/F"); recoilTree_->Branch("vz",&tree_vz_,"vz/F");
    recoilTree_->Branch("nHits",&tree_nHits_,"nHits/I");
    recoilTree_->Branch("eventID",&tree_eventID_,"eventID/I");
    recoilTree_->Branch("x",&tree_x_,"x/F"); recoilTree_->Branch("phi",&tree_phi_,"phi/F");
}

bool RecoilProcessor::process(IEvent* /*ievent*/) {
    if (!evth_ || !mcParts_ || !trks_ || !svtClusters_) {
        std::cerr << "[RecoilProcessor] Missing required branches.\n";
        return false;
    }

    constexpr double ETOL = 5e-3;
    const double BEAME    = (beamE_ > 0 ? beamE_ : 3.74);

    // --- Build distinct truth-layer counts per MC ID from clusters (findable) ---
    std::unordered_map<int, std::bitset<64>> maskByMC;
    maskByMC.reserve(8192);
    size_t clustersWithMC = 0;

    for (auto* cl : *svtClusters_) {
        if (!cl) continue;
        const int L = cl->getLayer();
        if (L < 0 || L >= 64) continue;
        const std::vector<int>& mcids = cl->getMCPartIDs();
        if (!mcids.empty()) ++clustersWithMC;
        for (int mcid : mcids) maskByMC[mcid].set(static_cast<size_t>(L));
    }
    if (debug_ > 0 && clustersWithMC == 0) {
        std::cout << "[RecoilProcessor] WARNING: cluster branch has ZERO MCPartIDs; "
                     "findable=0, acceptance→0, efficiency undefined.\n";
    }

    std::unordered_map<int,int> truthLayersByMCID;
    truthLayersByMCID.reserve(maskByMC.size());
    for (auto& kv : maskByMC) truthLayersByMCID[kv.first] = static_cast<int>(kv.second.count());

    // --- Per-event counters (including matched multiplicity sums for numerators) ---
    int nRecoils = 0;
    std::array<int, kNThr_> accDenCnt{}; accDenCnt.fill(0);
    std::array<int, kNThr_> accNumCnt{}; accNumCnt.fill(0);
    std::array<int, kNThr_> effDenCnt{}; effDenCnt.fill(0);
    std::array<int, kNThr_> effNumCnt{}; effNumCnt.fill(0);
    std::array<long long, kNThr_> effNumMatchedTracksSum{}; effNumMatchedTracksSum.fill(0);

    const int thrMin = thrList_.back(); // e.g. 6

    // --- Loop MC recoils (ELECTRONS) ---
    for (auto* part : *mcParts_) {
        if (!part) continue;
        if (std::abs(part->getPDG()) != 13) continue;              // use Muon
        if (std::abs(part->getEnergy() - BEAME) < ETOL) continue;   // skip beam leg
        ++nRecoils;

        const int pid = part->getID();
        const int nTruthLayers = (truthLayersByMCID.count(pid) ? truthLayersByMCID[pid] : 0);

        const auto& mom = part->getMomentum();
        const double px = mom[0], py = mom[1], pz = mom[2];
        const double pt = std::hypot(px, py);
        const double p  = std::sqrt(px*px + py*py + pz*pz);
        const double tanL = (pz != 0.0) ? (py / pz) : 999.; // analysis convention

        // --- Simple ID-equality matching over reco tracks ---
        int matchedCount = 0;
        int bestRecoLayers = 0;
        double bestPurity = 0.0;

        for (size_t t = 0; t < trks_->size(); ++t) {
            auto* trk = trks_->at(static_cast<int>(t));
            if (!trk) continue;
            if (trk->getID() != pid) continue; // strict ID equality

            ++matchedCount;

            // Distinct layer count on this track
            std::array<bool, 64> layerMask{}; layerMask.fill(false);
            int Nh = trk->getSvtHits().GetEntriesFast();
            int nHitsTotal = 0, nHitsFromPID = 0;
            for (int ih = 0; ih < Nh; ++ih) {
                auto* th = static_cast<::TrackerHit*>(trk->getSvtHits().At(ih));
                if (!th) continue;
                const int L = th->getLayer();
                if (0 <= L && L < kNLayersMax_) layerMask[static_cast<size_t>(L)] = true;
                ++nHitsTotal;

                // Simple purity: does this hit carry the same MC ID?
                const auto& hitMC = th->getMCPartIDs();
                if (!hitMC.empty()) {
                    for (int mcid : hitMC) { if (mcid == pid) { ++nHitsFromPID; break; } }
                }
            }
            int nRecoLayersThis = 0;
            for (int L=0; L<kNLayersMax_; ++L) if (layerMask[static_cast<size_t>(L)]) ++nRecoLayersThis;

            const double purityThis = (nHitsTotal > 0 ? double(nHitsFromPID) / double(nHitsTotal) : 0.0);
            if (nRecoLayersThis > bestRecoLayers || (nRecoLayersThis == bestRecoLayers && purityThis > bestPurity)) {
                bestRecoLayers = nRecoLayersThis;
                bestPurity     = purityThis;
            }
        }

        // --- Global diagnostics based on multiplicity for findable@thrMin ---
        const bool isFindableMin = (nTruthLayers >= thrMin);
        if (isFindableMin) {
            if (matchedCount == 0) {
                if (h_B_zero_tanL_) h_B_zero_tanL_->Fill(tanL);
                if (h_B_zero_p_)    h_B_zero_p_->Fill(p);
            } else if (matchedCount > 1) {
                if (h_B_overlap_tanL_) h_B_overlap_tanL_->Fill(tanL);
                if (h_B_overlap_p_)    h_B_overlap_p_->Fill(p);
            }
        }

        // --- Fill acceptance/efficiency per threshold ---
        for (int i = 0; i < kNThr_; ++i) {
            const int thr = thrList_[i];

            // Acceptance
            if (acc_den_tanL_[i]) acc_den_tanL_[i]->Fill(tanL);
            if (acc_den_p_[i])    acc_den_p_[i]->Fill(p);
            ++accDenCnt[i];

            if (nTruthLayers >= thr) {
                if (acc_num_tanL_[i]) acc_num_tanL_[i]->Fill(tanL);
                if (acc_num_p_[i])    acc_num_p_[i]->Fill(p);
                ++accNumCnt[i];

                // Efficiency DEN
                if (eff_den_tanL_[i]) eff_den_tanL_[i]->Fill(tanL);
                if (eff_den_p_[i])    eff_den_p_[i]->Fill(p);
                ++effDenCnt[i];

                // Efficiency NUM: bestRecoLayers >= thr via ID-equality match
                if (bestRecoLayers >= thr) {
                    if (eff_num_tanL_[i]) eff_num_tanL_[i]->Fill(tanL);
                    if (eff_num_p_[i])    eff_num_p_[i]->Fill(p);
                    ++effNumCnt[i];
                    effNumMatchedTracksSum[i] += matchedCount; // how many tracks matched for this numerator
                }
            }
        }

        // --- Match-quality histos for the best match (if any) ---
        if (matchedCount > 0) {
            const int holes = std::max(0, nTruthLayers - bestRecoLayers);
            if (h_match_holes_)    h_match_holes_->Fill(holes);
            if (h_match_purity_)   h_match_purity_->Fill(bestPurity);
            const auto& vtx = part->getVertexPosition();
            const double vz = (vtx.size()>2) ? vtx[2] : 0.0;
            if (h_vz_vs_purity_)   h_vz_vs_purity_->Fill(vz, bestPurity);
            if (h_tanL_vs_purity_) h_tanL_vs_purity_->Fill(tanL, bestPurity);
        }

        // --- Tree payload ---
        tree_sampleID_  = sampleID_;
        tree_category_  = 0;
        tree_pdg_       = part->getPDG();
        tree_OriginPDG_ = part->getOriginPDG();
        tree_MomPDG_    = part->getMomPDG();
        tree_charge_    = part->getCharge();
        tree_ID_        = pid;

        tree_px_ = px; tree_py_ = py; tree_pz_ = pz;
        tree_p_  = p;  tree_pt_ = pt;
        tree_tanLambda_ = static_cast<float>(tanL);
        tree_theta_     = std::atan2(pt, pz) * TMath::RadToDeg();
        tree_energy_    = part->getEnergy();
        tree_mass_      = part->getMass();
        tree_time_      = part->getTime();
        const auto& vtx = part->getVertexPosition();
        tree_vz_        = (vtx.size()>2) ? vtx[2] : 0.0f;
        tree_nHits_     = nTruthLayers; // truth-layer count
        tree_eventID_   = evth_->getEventNumber();
        tree_x_         = 1.0f - float(part->getEnergy() / (beamE_>0?beamE_:1.0));
        tree_phi_       = std::atan2(py, px);
        if (recoilTree_) recoilTree_->Fill();
    }

    if (debug_ > 0) {
        std::cout << "[RecoilProcessor] Event " << evth_->getEventNumber()
                  << " recoils=" << nRecoils
                  << " clustersWithMCIDs=" << clustersWithMC << "\n";
        if (debug_ > 1) {
            for (int i=0;i<kNThr_;++i) {
                std::cout << "  thr>=" << thrList_[i]
                          << "  ACC: DEN=" << accDenCnt[i] << " NUM=" << accNumCnt[i]
                          << "  EFF: DEN=" << effDenCnt[i] << " NUM=" << effNumCnt[i]
                          << "  (sum matched tracks in NUM=" << effNumMatchedTracksSum[i] << ")\n";
            }
        }
    }

    return true;
}

void RecoilProcessor::finalize() {
    if (!outF_ || !outF_->IsOpen()) { std::cerr << "[RecoilProcessor] finalize(): output not open\n"; return; }
    outF_->cd();
    auto writeIf = [](TObject* o){ if (o) o->Write(); };

    // Efficiencies
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

        if (acc_den_tanL_[i] && acc_num_tanL_[i] && TEfficiency::CheckConsistency(*acc_num_tanL_[i], *acc_den_tanL_[i])) {
            auto* aTL = new TEfficiency(*acc_num_tanL_[i], *acc_den_tanL_[i]);
            aTL->SetName(Form("Acc_tanL_ge%d",thr));
            aTL->SetTitle(Form("Acceptance vs tan#lambda (>= %d)",thr));
            aTL->SetStatisticOption(TEfficiency::kFCP);
            aTL->Write();
        }
        if (acc_den_p_[i] && acc_num_p_[i] && TEfficiency::CheckConsistency(*acc_num_p_[i], *acc_den_p_[i])) {
            auto* aP = new TEfficiency(*acc_num_p_[i], *acc_den_p_[i]);
            aP->SetName(Form("Acc_p_ge%d",thr));
            aP->SetTitle(Form("Acceptance vs p (>= %d)",thr));
            aP->SetStatisticOption(TEfficiency::kFCP);
            aP->Write();
        }

        // Write raw
        writeIf(eff_den_tanL_[i]); writeIf(eff_num_tanL_[i]);
        writeIf(eff_den_p_[i]);    writeIf(eff_num_p_[i]);
        writeIf(acc_den_tanL_[i]); writeIf(acc_num_tanL_[i]);
        writeIf(acc_den_p_[i]);    writeIf(acc_num_p_[i]);
    }

    // Match-quality histos
    writeIf(h_match_purity_);
    writeIf(h_match_holes_);
    writeIf(h_vz_vs_purity_);
    writeIf(h_tanL_vs_purity_);
    writeIf(h_B_overlap_tanL_);
    writeIf(h_B_overlap_p_);
    writeIf(h_B_zero_tanL_);
    writeIf(h_B_zero_p_);

    writeIf(recoilTree_);
    outF_->Close();
}

DECLARE_PROCESSOR(RecoilProcessor);

