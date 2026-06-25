// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   HCalTestbeamTask.cxx
/// \author Ian Pascal Møller <ian.moeller@nbi.ku.dk>
/// \brief Standalone QC task for FoCal-H beam tests
///

#include <boost/crc.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <unordered_set>

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TMath.h>

#include "DataFormatsFOCAL/Constants.h"
#include "QualityControl/QcInfoLogger.h"
#include "FOCAL/HCalTestbeamTask.h"
#include <CommonConstants/Triggers.h>
#include <Framework/InputRecord.h>
#include <Framework/InputRecordWalker.h>
#include <Framework/DataRefUtils.h>
#include <DPLUtils/DPLRawParser.h>
#include <DetectorsRaw/RDHUtils.h>
#include <Headers/DataHeader.h>
#include <Headers/RDHAny.h>

namespace o2::quality_control_modules::focal
{

HCalTestbeamTask::~HCalTestbeamTask()
{
  if (mTFerrorCounter)      { delete mTFerrorCounter;     }
  if (mFEENumberHBF)        { delete mFEENumberHBF;       }
  if (mFEENumberTF)         { delete mFEENumberTF;        }
  if (mNumLinksTF)          { delete mNumLinksTF;         }
  if (mNumHBFPerCRU)        { delete mNumHBFPerCRU;       }
  if (mCRUcounter)          { delete mCRUcounter;         }
  if (mHCalGlobalADCSum)    { delete mHCalGlobalADCSum;   }
  if (mPayloadSizeHCalGBT)  { delete mPayloadSizeHCalGBT; }

  for (auto& hist : mHCalASICChannelADC) { delete hist; }
  for (auto& hist : mHCalASICChannelTOA) { delete hist; }
  for (auto& hist : mHCalASICChannelTOT) { delete hist; }
}

void HCalTestbeamTask::initialize(o2::framework::InitContext& /*ctx*/)
{
  mMapper = o2::focal::HCalMapper();

  auto get_bool = [](const std::string_view input) -> bool {
    return input == "true";
  };

  auto hasDebugParam = mCustomParameters.find("Debug");
  if (hasDebugParam != mCustomParameters.end()) {
    mDebugMode = get_bool(hasDebugParam->second);
  }
  
  std::fill(mHCalASICChannelADC.begin(), mHCalASICChannelADC.end(), nullptr);
  std::fill(mHCalASICChannelTOA.begin(), mHCalASICChannelTOA.end(), nullptr);
  std::fill(mHCalASICChannelTOT.begin(), mHCalASICChannelTOT.end(), nullptr);

  mTFerrorCounter = new TH1F("NumberOfTFerror", "Number of TFbuilder errors", 2, 0.5, 2.5);
  mTFerrorCounter->GetYaxis()->SetTitle("Time Frame Builder Error");
  mTFerrorCounter->GetXaxis()->SetBinLabel(1, "empty");
  mTFerrorCounter->GetXaxis()->SetBinLabel(2, "filled");
  getObjectsManager()->startPublishing(mTFerrorCounter);

  mFEENumberHBF = new TH1F("NumberOfHBFPerFEE", "Number of HBFs per FEE", 100001, -0.5, 100000.5);
  mFEENumberHBF->GetXaxis()->SetTitle("FEE ID");
  mFEENumberHBF->GetYaxis()->SetTitle("Number of HBFs");
  getObjectsManager()->startPublishing(mFEENumberHBF);

  mFEENumberTF = new TH1F("NumberOfTFPerFEE", "Number of TFs per FEE", 100001, -0.5, 100000.5);
  mFEENumberTF->GetXaxis()->SetTitle("FEE ID");
  mFEENumberTF->GetYaxis()->SetTitle("Number of TFs");
  getObjectsManager()->startPublishing(mFEENumberTF);

  mNumLinksTF = new TH1F("NumberOfLinksPerTF", "Number of Links / timeframe", 101, -0.5, 100.5);
  mNumLinksTF->GetXaxis()->SetTitle("Number of links");
  mNumLinksTF->GetYaxis()->SetTitle("Number of TFs");
  getObjectsManager()->startPublishing(mNumLinksTF);

  mNumHBFPerCRU = new TH1F("NumberOfHBFPerCRU", "Number of HBFs / CRU", 10001, -0.5, 10000.5);
  mNumHBFPerCRU->GetXaxis()->SetTitle("Endpoint ID");
  mNumHBFPerCRU->GetYaxis()->SetTitle("Link ID");
  getObjectsManager()->startPublishing(mNumHBFPerCRU);

  mCRUcounter = new TH2D("CRUcounter", "Number of HBFs per CRU link", 2, -0.5, 1.5, 21, -0.5, 20.5);
  mCRUcounter->GetXaxis()->SetTitle("Endpoint ID");
  mCRUcounter->GetYaxis()->SetTitle("Link ID");
  getObjectsManager()->startPublishing(mCRUcounter);

  mPayloadSizeTF = new TH1D("PayloadSizeTF", "Payload size TF", 256, 0., 0.);
  getObjectsManager()->startPublishing(mPayloadSizeTF);
 
  mPayloadSizeHCalGBT = new TH1D("PayloadSizeHCalGBT", "Payload Size GBT Words", 256, 0., 0.);
  getObjectsManager()->startPublishing(mPayloadSizeHCalGBT);


  mHCalGlobalADCSumCanvas = new TCanvas("HCalGlobalADCSum", "Global ADC Sum", 1920, 1080);
  mHCalGlobalADCSum = new THStack("HCalGlobalADCSumStack", "Global ADC Sum Stack");
  for (int sample = 0; sample < 16; ++sample) {
    TH1I* graph = new TH1I(Form("HCalGlobalADCSumSample_%02d", sample),
			                     Form("Sample %02d", sample),
			                     256, 0., 0.);

    mHCalGlobalADCSum->Add(graph);
    mHCalGlobalADCSumContainer[sample] = graph;
  }

  mHCalGlobalADCSumCanvas->cd(1);
  mHCalGlobalADCSum->Draw("plc nostack");
  gPad->BuildLegend();
  getObjectsManager()->startPublishing(mHCalGlobalADCSumCanvas);


  mHCalHeatmapCanvas = new TCanvas("HCalHeatmapCanvas", "HCal Heatmap Canvas", 1920, 1080);
  mHCalHeatmapCanvas->DivideSquare(16, 0.001, 0.001);
  for (int sample = 0; sample < 16; ++sample) {
    mHCalHeatmapCanvas->cd(sample+1);
    TH2D* graph = new TH2D(Form("HCalHeatmapSample_%02d", sample),
                           Form("Sample %d", sample),
                           16, -0.5, 16 - 0.5,
                           12, -0.5, 12 - 0.2);

    graph->SetStats(0);
    graph->Draw("COLZ");
    //graph->SetMinimum(0.);
    //graph->SetMaximum(1024.);
    gPad->SetLogz();
    mHCalHeatmapContainer[sample] = graph;
  }
  getObjectsManager()->startPublishing(mHCalHeatmapCanvas);


  for (int i = 0; i < mNUM_ASICS; ++i) {
    mHCalASICChannelADC[i] = new TH2D(Form("HCalADC_ASIC_%d", i),
                                      Form("ADC per channel (ASIC %d)", i),
                                      mCHANNELS_PER_ASIC, -0.5, mCHANNELS_PER_ASIC -0.5,
                                      mRANGE_ADC, 0., mRANGE_ADC);

    mHCalASICChannelTOT[i] = new TH2D(Form("HCalTOT_ASIC_%d", i),
                                      Form("TOT per channel (ASIC %d)", i),
                                      mCHANNELS_PER_ASIC, -0.5, mCHANNELS_PER_ASIC -0.5,
                                      mRANGE_TOT, 0., mRANGE_TOT);

    mHCalASICChannelTOA[i] = new TH2D(Form("HCalTOA_ASIC_%d", i),
                                      Form("TOA per channel (ASIC %d)", i),
                                      mCHANNELS_PER_ASIC, -0.5, mCHANNELS_PER_ASIC -0.5,
                                      mRANGE_TOA, 0., mRANGE_TOA);
  
    getObjectsManager()->startPublishing(mHCalASICChannelADC[i]);
    getObjectsManager()->startPublishing(mHCalASICChannelTOT[i]);
    getObjectsManager()->startPublishing(mHCalASICChannelTOA[i]);
  }

  for (int i = 0; i < o2::focal::constants::HCAL_NUM_GBT_LINKS; ++i) {
    for (int j = 0; j < o2::focal::constants::HCAL_NUM_ROCS_PER_LINK; ++j) {
      for (int k = 0; k < 2; ++k) {
	      TCanvas* currentCanvas = new TCanvas(Form("WaveformsCanvas_Link_%d_ROC_%d.%d", i, j, k), Form("WaveformsCanvas_Link_%d_ROC_%d.%d", i, j, k), 1920, 1080);
	      currentCanvas->DivideSquare(36, 0.001, 0.001);
	      mHCalWaveforms[i][j][k] = currentCanvas;
	      for (int chn = 0; chn < 36; ++chn) {
	        currentCanvas->cd(chn+1);
	        TH2D* graph = new TH2D(Form("HCalADC_Link_%d_ROC_%d_Half_%d_Chn_%d", i, j, k, chn),
	                       				 Form("ADC per sample (Link %d, ROC %d.%d, chn. %d)", i, j, k, chn),
	        				               16, 0, 15,
	        				               128, 0., 1024.);
	        graph->SetStats(0);
	        gPad->SetLogz();
	        graph->Draw();
	        mHCalWaveformsContainer[i][j][k][chn] = graph;
	      }
	    
      getObjectsManager()->startPublishing(currentCanvas);
      }
    }
  }

  // TODO: set up some enum to hold error types to indices or something
  int numChips = 2 * o2::focal::constants::HCAL_NUM_GBT_LINKS * o2::focal::constants::HCAL_NUM_ROCS_PER_LINK;
  mHCalDataErrors = new TH2I("HCalDataErrors", "HCal Data Errors",
                             numChips, -0.5, numChips - 0.5,  // count errors for each chip
                             5, -0.5, 5 - 0.5);               // number of types of errors

  mHCalDataErrors->GetXaxis()->SetTitle("Data Link (GBTLink:ROC.Half)");
  mHCalDataErrors->GetYaxis()->SetBinLabel(1, "CRC_CHECKSUM_MISMATCH");
  mHCalDataErrors->GetYaxis()->SetBinLabel(2, "DAQH_HEADER_MISMATCH");
  mHCalDataErrors->GetYaxis()->SetBinLabel(3, "DAQH_TRAILER_MISMATCH");
  mHCalDataErrors->GetYaxis()->SetBinLabel(4, "HAMMING_BITS_SET");
  mHCalDataErrors->GetYaxis()->SetBinLabel(5, "DECODER_ERROR");
  int c = 0;
  for (int i = 0; i < o2::focal::constants::HCAL_NUM_GBT_LINKS; ++i) {
    for (int j = 0; j < o2::focal::constants::HCAL_NUM_ROCS_PER_LINK; ++j) {
      for (int k = 0; k < 2; ++k) {
        mHCalDataErrors->GetXaxis()->SetBinLabel(++c, Form("%02d:%02d.%d", i, j, k));
      }
    }
  }
  //mHCalDataErrors->Draw("text colz");
  getObjectsManager()->startPublishing(mHCalDataErrors);

}

void HCalTestbeamTask::startOfActivity(const Activity& activity)
{
  ILOG(Debug, Devel) << "startOfActivity " << activity.mId << ENDM;
  reset();
}

void HCalTestbeamTask::startOfCycle()
{
  ILOG(Debug, Devel) << "startOfCycle" << ENDM;
}

void HCalTestbeamTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  constexpr unsigned short FEE_HCal = 0xBEEF;

  std::vector<char> recordBuffer;
  o2::InteractionRecord currentInteractionRecord;
  for (const auto& rawData : framework::InputRecordWalker(ctx.inputs())) {
    recordBuffer.clear();
    if (rawData.header != nullptr && rawData.payload != nullptr) {
      const auto payloadSize = o2::framework::DataRefUtils::getPayloadSize(rawData); 
      gsl::span<const char> databuffer(rawData.payload, payloadSize);

      long int currentpos = 0;
      while (currentpos < databuffer.size()) {
        auto rdh = reinterpret_cast<const o2::header::RDHAny*>(databuffer.data() + currentpos);
        auto pageSize = o2::raw::RDHUtils::getMemorySize(rdh);
        auto pageHeaderSize = o2::raw::RDHUtils::getHeaderSize(rdh);
        auto pagePayloadSize = pageSize - pageHeaderSize;

        auto feeID        = o2::raw::RDHUtils::getFEEID(rdh);
        auto triggerBC    = o2::raw::RDHUtils::getTriggerBC(rdh);
        auto triggerOrbit = o2::raw::RDHUtils::getTriggerOrbit(rdh);
        auto offset       = o2::raw::RDHUtils::getOffsetToNext(rdh);
        auto stop         = o2::raw::RDHUtils::getStop(rdh);
        
        LOGF(debug, "---- New HBF ----");
        currentInteractionRecord.bc = triggerBC;
        currentInteractionRecord.orbit = triggerOrbit;

        LOGF(debug, "FEE: %04x", feeID);
        LOGF(debug, "BC: %d", triggerBC);
        LOGF(debug, "Orbit: %d", triggerOrbit);
        LOGF(debug, "Offset to next: %d", offset);
        LOGF(debug, "Stop bit: %d", stop);

        auto pagePayload = databuffer.subspan(currentpos + pageHeaderSize, pagePayloadSize);
        std::copy(pagePayload.begin(), pagePayload.end(), std::back_inserter(recordBuffer));
        LOGF(debug, "Record buffer now contains %d bytes of payload data", recordBuffer.size());

        auto triggerType = o2::raw::RDHUtils::getTriggerType(rdh);
        if (triggerType & o2::trigger::HB) {
          if (stop) {
            if (recordBuffer.size() == 0) {
            } else {
              if (feeID == FEE_HCal) {
                processHCalEvent(recordBuffer);
              }
            }
          }
        }

        currentpos += offset;
      }
    }
  }
}

void HCalTestbeamTask::processHCalEvent(const gsl::span<const char> hcalpayload)
{
  mDecoder.reset();
  mDecoder.decodeBuffer(hcalpayload);
  if (not mDecoder.hasEventData()) {
    return;
  }

  if (not mDecoder.isDataValid()) {
    for (int i = 0; i < 8; ++i) {
      mHCalDataErrors->Fill(i, 4., 1.);
    }

    return;
  }

  ++mNumEvents;
  LOGF(info, "events: %d", mNumEvents);
  std::array<int, 2> numSamples = mDecoder.getNumSamplesRead();
  LOGF(debug, "Samples read: %02d %02d", numSamples[0], numSamples[1]);

  std::array<std::array<o2::focal::HCalGBTLink, o2::focal::constants::HCAL_NUM_GBT_LINKS>, o2::focal::constants::HCAL_NUM_SAMPLES_PER_EVENT> links = mDecoder.getData();

  int currentChannel = 0;
  for (int sample = 0; sample < o2::focal::constants::HCAL_NUM_SAMPLES_PER_EVENT; ++sample) {
    int globalADCsum = 0;
    for (int link_id = 0; link_id < 2; ++link_id) {
      o2::focal::HCalGBTLink currentLink = links[sample][link_id];
        for (int roc_id = 0; roc_id < 2; ++roc_id) {
          o2::focal::HCalROC currentROC = currentLink.getROC(roc_id);
          for (int half = 0; half < 2; ++half) {
            o2::focal::HCalROCDataLink currentHalf = currentROC.getChipHalf(half);
            bool errored = false;
            if (not checkCRC(currentHalf)) {
              errored = true;
              LOGF(debug, "CRC mismatch in link %02d:%02d.%d", link_id, roc_id, half);
              int binIndex = link_id * o2::focal::constants::HCAL_NUM_GBT_LINKS     *
                                       o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             roc_id  * o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             half;

              mHCalDataErrors->Fill(binIndex, 0., 1.);
            }
            if (not checkDAQHHeader(currentHalf)) {
              errored = true;
              LOGF(debug, "DAQH header mismatch in link %02d:%02d.%d", link_id, roc_id, half);
              int binIndex = link_id * o2::focal::constants::HCAL_NUM_GBT_LINKS     *
                                       o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             roc_id  * o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             half;

              mHCalDataErrors->Fill(binIndex, 1., 1.);
            }
            if (not checkDAQHTrailer(currentHalf)) {
              errored = true;
              LOGF(debug, "DAQH trailer mismatch in link %02d:%02d.%d", link_id, roc_id, half);
              int binIndex = link_id * o2::focal::constants::HCAL_NUM_GBT_LINKS     *
                                       o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             roc_id  * o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             half;

              mHCalDataErrors->Fill(binIndex, 2., 1.);
            }
            if (not checkHammingBits(currentHalf)) {
              errored = true;
              LOGF(debug, "Hamming decode error bits are set in link %02d:%02d.%d", link_id, roc_id, half);
              int binIndex = link_id * o2::focal::constants::HCAL_NUM_GBT_LINKS     *
                                       o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             roc_id  * o2::focal::constants::HCAL_NUM_ROCS_PER_LINK +
                             half;

              mHCalDataErrors->Fill(binIndex, 3., 1.);
            }

            //if (errored) {
            //  auto words = currentHalf.getWords();
            //  for (const auto& w : words) {
            //    LOGF(debug, "%08X", w);
            //  }
            //}

            for (int chn = 0; chn < 36; ++chn) {
              o2::focal::HCalChannel currentChannel = currentHalf.getChannel(chn);
              int adc = currentChannel.adc;
              int tot = currentChannel.tot;
              int toa = currentChannel.toa;

              // TODO: decide on how to fill these plots, i.e. which sample?
      	      if (sample == 3) {
                mHCalASICChannelADC[roc_id]->Fill(chn + 36 * half, adc);
                mHCalASICChannelTOT[roc_id]->Fill(chn + 36 * half, tot);
                mHCalASICChannelTOA[roc_id]->Fill(chn + 36 * half, toa);
	            }	    

              globalADCsum += adc;
              mHCalWaveformsContainer[link_id][roc_id][half][chn]->Fill(sample, adc);

              std::pair<int, int> coords = mMapper.getRowCol(link_id, roc_id, half, chn);
              if ((coords.first == -1) || (coords.second == -1)) {
                continue;
              } else {

                // Instead of simply filling histograms with ADC values,
                // we want to show the "average" ADC value on the heatmaps
                double previous = mHCalHeatmapContainer[sample]->GetBinContent(coords.second+1, coords.first+1) * mNumEvents; // "unaveraged" value from current cell
                mHCalHeatmapContainer[sample]->SetBinContent(coords.second+1, coords.first+1, (previous + adc) / (mNumEvents + 1)); // new average value
              }
            }
          }
        }
    }

  mHCalGlobalADCSumContainer[sample]->Fill(globalADCsum);
  }
}

void HCalTestbeamTask::endOfCycle()
{
  ILOG(Debug, Devel) << "endOfCycle" << ENDM;
}

void HCalTestbeamTask::endOfActivity(const Activity& /*activity*/)
{
  ILOG(Debug, Devel) << "endOfActivity" << ENDM;
}

void HCalTestbeamTask::reset()
{
  ILOG(Debug, Devel) << "Resetting the histograms" << ENDM;
  mTFerrorCounter ->Reset();
  mFEENumberHBF   ->Reset();
  mFEENumberTF    ->Reset();
  mNumLinksTF     ->Reset();
  mNumHBFPerCRU   ->Reset();
  mCRUcounter     ->Reset();
  mPayloadSizeTF  ->Reset();

  mPayloadSizeHCalGBT ->Reset();

  for (auto& hist : mHCalGlobalADCSumContainer) {
    if (hist) { hist->Reset(); }
  }

  for (auto& hist : mHCalASICChannelADC) {
      if (hist) { hist->Reset(); }
  }
  for (auto& hist : mHCalASICChannelTOT) {
      if (hist) { hist->Reset(); }
  }
  for (auto& hist : mHCalASICChannelTOA) {
      if (hist) { hist->Reset(); }
  }
}

bool HCalTestbeamTask::isLostTimeframe(framework::ProcessingContext& ctx) const
{
  // direct data
  constexpr auto originFOC = header::gDataOriginFOC;
  o2::framework::InputSpec dummy{ "dummy",
                                  framework::ConcreteDataMatcher{ originFOC,
                                                                  header::gDataDescriptionRawData,
                                                                  0xDEADBEEF } };
  for (const auto& ref : o2::framework::InputRecordWalker(ctx.inputs(), { dummy })) {
    const auto dh = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
    const auto payloadSize = o2::framework::DataRefUtils::getPayloadSize(ref);
    if (payloadSize == 0) {
      return true;
    }
  }

  // sampled data
  o2::framework::InputSpec dummyDS{ "dummyDS",
                                    framework::ConcreteDataMatcher{ "DS",
                                                                    "focrawdata0",
                                                                    0xDEADBEEF } };
  for (const auto& ref : o2::framework::InputRecordWalker(ctx.inputs(), { dummyDS })) {
    const auto dh = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
    const auto payloadSize = o2::framework::DataRefUtils::getPayloadSize(ref);
    if (payloadSize == 0) {
      return true;
    }
  }

  return false;
}

bool HCalTestbeamTask::checkDAQHHeader(o2::focal::HCalROCDataLink half) {
  return (half.getHeader().hd == 0b1111);
}

bool HCalTestbeamTask::checkDAQHTrailer(o2::focal::HCalROCDataLink half) {
  return (half.getHeader().tr == 0b0101);
}

bool HCalTestbeamTask::checkHammingBits(o2::focal::HCalROCDataLink half) {
  return (half.getHeader().hm == 0);
}

bool HCalTestbeamTask::checkCRC(o2::focal::HCalROCDataLink half) {
  std::array<unsigned int, o2::focal::constants::HCAL_NUM_GBT_LINES_PER_LINK> words = half.getWords();
  unsigned int crcWord = words[o2::focal::constants::HCAL_NUM_GBT_LINES_PER_LINK - 1];

  boost::crc_basic<32> crc_32(0x04C11DB7, 0x00000000, 0x00000000, false, false);
  for (int i = 0; i < o2::focal::constants::HCAL_NUM_GBT_LINES_PER_LINK - 1; ++i) {
    unsigned char bytes[4] = {
      static_cast<unsigned char>(words[i] >> 24),
      static_cast<unsigned char>(words[i] >> 16),
      static_cast<unsigned char>(words[i] >> 8 ),
      static_cast<unsigned char>(words[i]      ),
    };

    crc_32.process_bytes(bytes, 4);
  }
  
  return (crc_32.checksum() == crcWord);
}

} // namespace o2::quality_control_modules::focal

