#ifndef __TRACK_EFFICIENCYPROCESSOR_H__
#define __TRACK_EFFICIENCYPROCESSOR_H__


//HPSTR
#include "HpsEvent.h"
#include "Collections.h"
#include "EventHeader.h"
#include "Vertex.h"
#include "CalCluster.h"
#include "Track.h"
#include "TrackerHit.h"
#include "MCTrackerHit.h"
#include "MCParticle.h"
#include "Particle.h"
#include "Processor.h"
#include "BaseSelector.h"
#include "TrackEfficHistos.h"
#include "ThreeProngHistos.h"
#include "FlatTupleMaker.h"
//#include "AnaHelpers.h"
//#include "Histos.h"


#include "IMPL/MCParticleImpl.h" 
#include "EVENT/MCParticle.h"  





//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"

//C++ 
#include <memory>
#include <array>
#include <unordered_map>


class TrackerHit; // fwd decl
class TBranch;
class Track;


class RecoilProcessor : public Processor {
    
public:
    RecoilProcessor(const std::string& name, Process& process);
    ~RecoilProcessor();
    virtual bool process(IEvent* ievent);

    virtual void initialize(TTree* tree);
    
    virtual void finalize();
    
    virtual void configure(const ParameterSet& parameters);
    
    //sample identifier BH = 0, Ap = 1
    int   sampleID_{-1};

private:
    
    std::shared_ptr<BaseSelector> cluSelector;
    std::shared_ptr<BaseSelector> trkSelector;
    std::vector<std::string> regionSelections_;
    std::vector<std::string> threeProngSelections_;
    
    std::string selectionCfg_;
    std::string trkSelCfg_;
    std::string cluSelCfg_;



    TBranch* bmcParts_{nullptr}; //!< description
    TBranch* bmcTrkrHits_{nullptr}; //!< description

    std::vector<MCTrackerHit*> * mcTrkrHits_{nullptr};
    std::vector<MCParticle*>   * mcParts_{nullptr}; //!< description

    TFile* outF_{nullptr}; 

    TBranch* bfspart_{nullptr};
    TBranch* bclus_{nullptr};
    TBranch* btrks_{nullptr};
    TBranch* bevth_{nullptr};
    
    std::vector<CalCluster*> * clus_{};
    std::vector<Particle*> * fspart_{};
    std::vector<Track*>  * trks_{};
    EventHeader* evth_{nullptr};
    
    std::string anaName_{"trkEffAna"};
    std::string cluColl_{"ECalClusters"};
    std::string fspartColl_{"FinalStateParticles"};
    std::string trkColl_{"GBLTracks"};
    TTree* tree_{nullptr};

    std::shared_ptr<TrackEfficHistos> _trkeff_histos;
    
    std::map<std::string, std::shared_ptr<BaseSelector> > _reg_trkeff_selectors;
    std::map<std::string, std::shared_ptr<TrackEfficHistos> > _reg_trkeff_histos;
    std::vector<std::string> _regions;

    std::map<std::string, std::shared_ptr<BaseSelector> > _reg_three_prong_trkeff_selectors;
    std::map<std::string, std::shared_ptr<ThreeProngHistos> > _reg_three_prong_trkeff_histos;
    std::vector<std::string> _three_prong_regions;
    

    typedef std::map<std::string,std::shared_ptr<TrackEfficHistos> >::iterator reg_it;
    typedef std::map<std::string,std::shared_ptr<ThreeProngHistos> >::iterator three_prong_reg_it;

    std::string histoCfg_{""};
    std::string thrProngCfg_{""};
    std::string cluHistoCfg_{""};
    double timeOffset_{-999};
    //In GeV. Default is 2016 value;
    double beamE_{2.3};
    int isData{0};
    std::shared_ptr<AnaHelpers> _ah;
    /*
    struct TridentCand{
	std::pair<CalCluster*,Track*> ele;
	std::pair<CalCluster*,Track*> pos;
    };
    */
    struct TridentCand{
	Particle* ele;
	Particle* pos;
    };

    struct ThreeProngCand{
	Particle* ele;
	Particle* pos;
	Particle* recoil;
    };

     
     //trees to best later process information
       // A little tree to hold truth‐level electron info
    TTree* recoilTree_{nullptr};

    // Branches
    Int_t   tree_sampleID_;
    Int_t   tree_category_;    // 0=recoil,1=radiated,2=secondary
    Int_t   tree_pdg_, tree_OriginPDG_, tree_MomPDG_, tree_charge_, tree_ID_;
    Float_t tree_px_, tree_py_, tree_pz_;
    Float_t tree_p_,  tree_pt_, tree_tanLambda_, tree_theta_;
    Float_t tree_energy_, tree_mass_, tree_time_, tree_vz_;
    Int_t   tree_nHits_;
    
    Float_t tree_x_;           // energy‐fraction of the mediator 
    Float_t tree_phi_;         // azimuthal angle of the electron
    Int_t   tree_eventID_;     // for merging/overlaps
    
     //histograms
     TH1F* h_rec_p_{nullptr};
     TH1F* h_rec_theta_{nullptr};
     TH1F* h_rec_vtz_{nullptr};
     TH2F* h_rec_tanlam_p_{nullptr};
     TH1F* h_rec_nhits_{nullptr};
     TH2F* h_rec_tanlam_vs_nhits_{nullptr};
     
     TH1F* h_rad_p_{nullptr};
     TH1F* h_rad_theta_{nullptr};
     TH1F* h_rad_vtz_{nullptr};
     TH2F* h_rad_tanlam_p_{nullptr};












    //more truth info
    // --- TruthTracks branch (LCIO Truth Track ID space)
    std::vector<Track*>* truthTrks_{nullptr};
    TBranch*             btruthTrks_{nullptr};
    
   std::string truthTracksCollRoot_{""};
   //track matching and efficiency functionality


	// New input collection name
	std::string sclusColl_ = "SvtClusters";

	// Branch ptrs
	std::vector<TrackerHit*>* svtClusters_{nullptr};
	TBranch* bsvtClusters_{nullptr};

	// Binning/thresholds
	static constexpr int kNLayersMax_ = 14;
	static constexpr int kNThr_ = 5;
	const std::array<int,kNThr_> thrList_{{10,9,8,7,6}};

	// Efficiency histograms (per-threshold, 1D vs tanλ and p)
	std::array<TH1F*, kNThr_> eff_den_tanL_{};
	std::array<TH1F*, kNThr_> eff_num_tanL_{};
	std::array<TH1F*, kNThr_> eff_den_p_{};
	std::array<TH1F*, kNThr_> eff_num_p_{};

	// Signal/Background shapes + derived S/sqrt(B)
	TH1F* h_S_tanL_{nullptr};
	TH1F* h_B_tanL_{nullptr};
	TH1F* h_S_p_{nullptr};
	TH1F* h_B_p_{nullptr};
	TH1F* h_SoverSqrtB_tanL_{nullptr};
	TH1F* h_SoverSqrtB_p_{nullptr};





     // Per-threshold S/B shapes and derived S/sqrt(B) this takes into account denominator having different number of hits
	std::array<TH1F*, kNThr_> h_S_tanL_thr_{};
	std::array<TH1F*, kNThr_> h_B_tanL_thr_{};
	std::array<TH1F*, kNThr_> h_S_p_thr_{};
	std::array<TH1F*, kNThr_> h_B_p_thr_{};
	std::array<TH1F*, kNThr_> h_SoverSqrtB_tanL_thr_{};
	std::array<TH1F*, kNThr_> h_SoverSqrtB_p_thr_{};
                                                                                                                          	// Matching quality
	TH1F* h_match_purity_{nullptr};   // fraction of hits from matched recoil
	TH1F* h_match_holes_{nullptr};    // # of truth layers missing on track
	TH2F* h_vz_vs_purity_{nullptr};
	TH2F* h_tanL_vs_purity_{nullptr};

        TH1F* h_B_overlap_tanL_{nullptr};
	TH1F* h_B_overlap_p_{nullptr};
	TH1F* h_B_zero_tanL_{nullptr};
	TH1F* h_B_zero_p_{nullptr};



      //acceptance histos
          // Acceptance (findable / all generated) histograms per threshold
    std::array<TH1F*, kNThr_> acc_den_tanL_{};  // DEN: all generated recoils (per-thr copy)
    std::array<TH1F*, kNThr_> acc_num_tanL_{};  // NUM: findable recoils (>=thr)
    std::array<TH1F*, kNThr_> acc_den_p_{};     // DEN: all generated recoils (per-thr copy)
    std::array<TH1F*, kNThr_> acc_num_p_{};     // NUM: findable recoils (>=thr)



    //Debug level
    int debug_{0};
};

#endif
