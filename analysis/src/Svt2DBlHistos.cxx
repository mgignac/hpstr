#include "Svt2DBlHistos.h"
#include <math.h>
#include "TCanvas.h"

Svt2DBlHistos::Svt2DBlHistos(const std::string& inputName) {
    m_name = inputName;
    mmapper_ = new ModuleMapper();
}

Svt2DBlHistos::~Svt2DBlHistos() {

    for (std::map<std::string, TGraphErrors*>::iterator it = baselineGraphs.begin(); 
            it!=baselineGraphs.end(); ++it) {
        if (it->second) {
            delete (it->second);
            it->second = nullptr;
        }
    }
    baselineGraphs.clear();
}

void Svt2DBlHistos::get2DHistoOccupancy(std::vector<std::string> histos2dNames) {
}

void Svt2DBlHistos::FillHistograms(std::vector<RawSvtHit*> *rawSvtHits_,float weight) {

    //std::cout << "[Svt2DBlHistos] FillHistograms" << std::endl;
    int nhits = rawSvtHits_->size();
    std::vector<std::string> hybridStrings={};
    std::string histokey;
    if(Event_number%1000 == 0) std::cout << "Event: " << Event_number 
        << " Number of RawSvtHits: " << nhits << std::endl;
    //std::cout << "Event: " << Event_number 

    //Following Block counts the total number of hits each hybrid records per event
    int svtHybMulti[4][15] = {0};
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        int mod = rawSvtHit->getModule();
        int lay = rawSvtHit->getLayer();
        svtHybMulti[mod][lay]++;

    }
    for (int i =0; i < 4; i++)
    {
        for (int j = 1; j < 15; j++)
        {
            if (!(j<9 && i>1))
            {   
                std::string swTag = mmapper_->getStringFromSw("ly"+std::to_string(j)+"_m"+std::to_string(i));
                hybridStrings.push_back(swTag);
                Fill1DHisto("hitN_"+swTag+"_h", svtHybMulti[i][j],weight);
            }
        }
    }

    Fill1DHisto("svtHitN_h", nhits,weight);
    //End of counting block

    //Populates histograms for each hybrid
    for (int i = 0; i < nhits; i++)
    {
        RawSvtHit* rawSvtHit = rawSvtHits_->at(i);
        auto mod = std::to_string(rawSvtHit->getModule());
        auto lay = std::to_string(rawSvtHit->getLayer());
        std::string swTag= mmapper_->getStringFromSw("ly"+lay+"_m"+mod);
        
        //Manually select which baselines (0 - 6) are included. THIS MUST MATCH THE JSON FILE!

        int ss = 0;
        
            //if(debug_ > 0) std::cout << "Filling Histogram with RawSvtHit" << std::endl; 
            histokey = "baseline"+std::to_string(ss)+"_"+swTag+"_hh";
                        Fill2DHisto(histokey, 
                    (float)rawSvtHit->getStrip(),
                    (float)rawSvtHit->getADCs()[ss], 
                    weight);
            //if(debug_ > 0) std::cout << "Histogram Filled" << std::endl; 
        

        ss = 3;
        
            //if(debug_ > 0) std::cout << "Filling Histogram with RawSvtHit" << std::endl; 
            histokey = "baseline"+std::to_string(ss)+"_"+swTag+"_hh";
                        Fill2DHisto(histokey, 
                    (float)rawSvtHit->getStrip(),
                    (float)rawSvtHit->getADCs()[ss], 
                    weight);
            //if(debug_ > 0) std::cout << "Histogram Filled" << std::endl; 
        
    }

            Event_number++;
}      
