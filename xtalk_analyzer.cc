// Studying cross-talk effect in CMS Pixel Detector Barrel Layer1
// Author: Arash Jofrehei
// Date: 26 November, 2018

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerFwd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include <math.h>

using namespace std;
using namespace edm;

class MonitorElement;
class PixelDigi;
class Phase2TrackerDigi;
class TrackerGeometry;

class xtalk_analyzer : public DQMEDAnalyzer{

public:
  
  explicit xtalk_analyzer(const edm::ParameterSet& conf);
  virtual ~xtalk_analyzer() override;
  void bookHistograms(DQMStore::IBooker & ibooker,
		      edm::Run const &  iRun ,
		      edm::EventSetup const &  iSetup ) override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  MonitorElement* ClusterWidth_XY;
  MonitorElement* x_res;
  MonitorElement* y_res;

private:

  void bookLayerHistos(DQMStore::IBooker & ibooker, unsigned int det_id, const TrackerTopology* tTopo); 
  //void fillClusterWidthHisto(const edm::Handle<edm::DetSetVector<PixelDigi>>  handle, const edm::ESHandle<TrackerGeometry> gHandle, edm::Handle<SiPixelRecHitCollection> recHitColl);
  edm::ParameterSet config_;
  TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  std::string geomType_; 
  edm::InputTag itPixelDigiSrc_; 
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit>> pixelRecHits_token;
  const edm::EDGetTokenT< edm::DetSetVector<PixelDigi> > itPixelDigiToken_;
  edm::ESHandle<TrackerTopology> tTopoHandle_;
};

// Constructer
xtalk_analyzer::xtalk_analyzer(const edm::ParameterSet& iConfig) :
  config_(iConfig),
  trackerHitAssociatorConfig_(iConfig, consumesCollector()),
  geomType_(iConfig.getParameter<std::string>("GeometryType")),
  itPixelDigiSrc_(iConfig.getParameter<edm::InputTag>("InnerPixelDigiSource")),
  pixelRecHits_token(consumes< edmNew::DetSetVector<SiPixelRecHit>>(edm::InputTag("siPixelRecHits"))),
  itPixelDigiToken_(consumes< edm::DetSetVector<PixelDigi> >(itPixelDigiSrc_))
{
  edm::LogInfo("xtalk_analyzer") << ">>> Construct xtalk analyzer ";
}

// destructor
xtalk_analyzer::~xtalk_analyzer() {}

void xtalk_analyzer::analyze(const edm::Event& evt, const edm::EventSetup& evtset)
{
  std::vector<PSimHit> matched;
  const PSimHit* closest_simhit=nullptr;
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  evt.getByToken(pixelRecHits_token, recHitColl);
  TrackerHitAssociator associate(evt, trackerHitAssociatorConfig_);
  
  if((recHitColl.product())->dataSize() > 0) {
    for ( auto recHitIdIterator : *(recHitColl.product())) {
      SiPixelRecHitCollection::DetSet detset = recHitIdIterator;
      if( detset.empty() ) continue;
      for ( auto iterRecHit : detset ) {
        matched.clear();
        matched = associate.associateHit(iterRecHit);
        LocalPoint lp = iterRecHit.localPosition();
	      float rechit_x = lp.x();
	      float rechit_y = lp.y();
        //ClusterWidth_XY->Fill(rechit_x,rechit_y,iterRecHit.cluster()->sizeX());
        if (matched.empty()) continue;
        float closest = 9999.9;
        for(auto const& m : matched){
	        float sim_x1 = m.entryPoint().x();
	        float sim_x2 = m.exitPoint().x();
	        float sim_xpos = 0.5*(sim_x1+sim_x2);
	        float sim_y1 = m.entryPoint().y();
	        float sim_y2 = m.exitPoint().y();
	        float sim_ypos = 0.5*(sim_y1+sim_y2);
	        float x_resi = sim_xpos - rechit_x;
	        float y_resi = sim_ypos - rechit_y;
	        float dist = sqrt(x_resi*x_resi + y_resi*y_resi);
          if ( dist < closest ) {
            closest = dist;
            closest_simhit = &m;
	        }
        } // end of loop on simhits
        float simhit_x = 0.5*(closest_simhit->entryPoint().x() + closest_simhit->exitPoint().x());
        float simhit_y = 0.5*(closest_simhit->entryPoint().y() + closest_simhit->exitPoint().y());
        x_res->Fill(10000*(rechit_x - simhit_x));
        y_res->Fill(10000*(rechit_y - simhit_y));
        float pixelSizeX = 0.0025; //cm
        float pixelSizeY = 0.0100; //cm
        float modified_x = fmod(fmod(fmod(simhit_x-0.5*pixelSizeX,2*pixelSizeX)+0.5*pixelSizeX,2*pixelSizeX)-pixelSizeX,2*pixelSizeX)+0.5*pixelSizeX;
        float modified_y = fmod(simhit_y-0.5*pixelSizeY,pixelSizeY);        
        ClusterWidth_XY->Fill(modified_x,modified_y,iterRecHit.cluster()->sizeX());
        //ClusterWidth_XY->Fill(simhit_x,simhit_y,iterRecHit.cluster()->sizeX());
      } // end of loop on reconstructed hits
    } // end of loop on detset
  } // end of if dataSize > 0
  evtset.get<TrackerTopologyRcd>().get(tTopoHandle_);
  edm::Handle< edm::DetSetVector<PixelDigi> > pixDigiHandle;
  evt.getByToken(itPixelDigiToken_, pixDigiHandle);
  edm::ESHandle<TrackerGeometry> geomHandle;
  evtset.get<TrackerDigiGeometryRecord>().get(geomType_, geomHandle);
  //fillClusterWidthHisto(pixDigiHandle, geomHandle , recHitColl);
} // end of xtalk_analyzer::analyze

/*
void xtalk_analyzer::fillClusterWidthHisto(const edm::Handle<edm::DetSetVector<PixelDigi>>  handle, const edm::ESHandle<TrackerGeometry> gHandle, edm::Handle<SiPixelRecHitCollection> recHitColl) {
  
  const edm::DetSetVector<PixelDigi>* digis = handle.product();

  const TrackerTopology* tTopo = tTopoHandle_.product();
  const TrackerGeometry* tGeom = gHandle.product();  

  for (typename edm::DetSetVector<PixelDigi>::const_iterator DSViter = digis->begin(); DSViter != digis->end(); DSViter++) {
    unsigned int rawid = DSViter->id;  
    int layer = tTopo->getITPixelLayerNumber(rawid);
    if (layer != 1) continue;

    const DetId detId(rawid);

    if (DetId(detId).det() != DetId::Detector::Tracker) continue;
  
    const GeomDetUnit* gDetUnit = tGeom->idToDetUnit(detId);

    const Phase2TrackerGeomDetUnit* tkDetUnit = dynamic_cast<const Phase2TrackerGeomDetUnit*>(gDetUnit);
    int nRows     = tkDetUnit->specificTopology().nrows();
    int nColumns  = tkDetUnit->specificTopology().ncolumns();
    if (nRows*nColumns == 0) continue;  

    int nDigi = 0; 
    int row_last = -1;
    int col_last = -1;
    int nclus = 0;
    int width = 1;
    int position = 0; 
    std::vector<int> charges;
    for (typename edm::DetSet< PixelDigi >::const_iterator di = DSViter->begin(); di != DSViter->end(); di++) {
      int col = di->column(); // column
      int row = di->row();    // row
      int adc = di->adc();    // digi charge 
      nDigi++;
      if (row_last == -1 ) {
        position = row+1;
        nclus++; 
        charges.push_back(adc);
      } else {
	if (abs(row - row_last) == 1 && col == col_last) {
	  position += row+1;
	  width++;
	  charges.push_back(adc);
	} else {
          position /= width;  
          //ClusterWidth_XY->Fill(col,position,width);
          
	  charges.clear();
          charges.push_back(adc);
	  width  = 1;
	  position = row+1;
          nclus++;
	}
      }
      row_last = row;
      col_last = col;
    }
  }
}

*/
// -- Book Histograms
//
void xtalk_analyzer::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const &  iRun , edm::EventSetup const &  iSetup ) {
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle_);
  const TrackerTopology* const tTopo = tTopoHandle_.product();
  edm::ESHandle<TrackerGeometry> geom_handle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geomType_, geom_handle);
  const TrackerGeometry* tGeom = geom_handle.product();
  for (auto const & det_u : tGeom->detUnits()) {
    unsigned int detId_raw = det_u->geographicalId().rawId();
    bookLayerHistos(ibooker,detId_raw, tTopo);
  }
}

void xtalk_analyzer::bookLayerHistos(DQMStore::IBooker & ibooker, unsigned int det_id, const TrackerTopology* tTopo) {
  //int layer = tTopo->getITPixelLayerNumber(det_id);
  //if (layer < 0) return;
  ibooker.cd();
  std::string top_folder = config_.getParameter<std::string>("TopFolderName");
  ibooker.setCurrentFolder(top_folder);
  float pixelSizeX = 0.0025; //cm
  float pixelSizeY = 0.0100; //cm
  //ClusterWidth_XY = ibooker.bookProfile2D("ClusterWidth XY", "cluster width wrt simulated XY - 20% xtalk - layer 1", 20, -0.0025, 0.0025, 40, -0.01, 0.01,0,20);
  ClusterWidth_XY = ibooker.bookProfile2D("ClusterWidth XY", "cluster width wrt simulated XY - 20% xtalk - layer 1", 30, -1.5*pixelSizeX, 1.5*pixelSizeX, 1, -1*pixelSizeY, pixelSizeY,0,20);
  x_res = ibooker.book1D("x_res","Residual along X - eta 0 - 20% xtalk;recoHit_x - simHit_x (um)",100,-100,100);
  y_res = ibooker.book1D("y_res","Residual along Y - eta 0 - 20% xtalk;recoHit_y - simHit_y (um)",100,-100,100);
}

// Define as a plugin
DEFINE_FWK_MODULE(xtalk_analyzer);