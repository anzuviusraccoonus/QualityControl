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
#include <TPaveText.h>
#include <TStyle.h>
#include <TLegend.h>

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

using namespace o2::focal::constants;

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

  delete mHCalGlobalADCSumCanvas;
  delete mHCalGlobalADCSum;
  
  delete mHCalSamplesPerEventCanvas;
  delete mHCalSamplesPerEvent;

  delete mHCalHeatmapCanvas;

  delete mHCalDataErrors;

  for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
    delete mHCalSamplesPerEventContainer[i];

    for (int j = 0; j < HCAL_NUM_ROCS_PER_LINK; ++j) {
      delete mHCalROCADC[i][j];
      delete mHCalROCTOT[i][j];
      delete mHCalROCTOA[i][j];

      for (int k = 0; k < 2; ++k) {
        delete mHCalWaveforms[i][j][k];

        for (int chn = 0; chn < HCAL_NUM_CHANNELS_PER_ROC_HALF; ++chn) {
          delete mHCalWaveformsContainer[i][j][k][chn];
        }
      }
    }
  }

  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    delete mHCalGlobalADCSumContainer[s];
    delete mHCalHeatmapContainer[s];
    delete mHCalADCSaturation[s];
    delete mHCalTOTSaturation[s];
    delete mHCalTOTZero[s];
    delete mHCalTOTBelowHalf[s];
  }
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


  /////////////////////////////////////////////////////////////////////////////////////
  /// Global ADC sum plots

  mHCalGlobalADCSumCanvas = new TCanvas("HCalGlobalADCSum", "Global ADC Sum", 1920, 1080);
  mHCalGlobalADCSum       = new THStack("HCalGlobalADCSumStack", "Global ADC Sum");
  gStyle->SetPalette(kCMYK);
  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    TH1I* graph = new TH1I(Form("HCalGlobalADCSumSample_%d", s),
                           Form("Sample %d", s),
                           256, 0., 0.);

    mHCalGlobalADCSum->Add(graph);
    mHCalGlobalADCSumContainer[s] = graph;
  }

  mHCalGlobalADCSumCanvas->cd(1);
  mHCalGlobalADCSum->Draw("plc nostack");
  TLegend* legend = gPad->BuildLegend();
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  getObjectsManager()->startPublishing(mHCalGlobalADCSumCanvas);


  /////////////////////////////////////////////////////////////////////////////////////
  /// Samples per event plots

  mHCalSamplesPerEventCanvas = new TCanvas("HCalSamplesPerEvent", "Number of Samples per Event", 1920, 1080);
  mHCalSamplesPerEvent       = new THStack("HCalSamplesPerEventStack", "Number of Samples per Event");
  gStyle->SetPalette(kCMYK);
  for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
    TH1I* graph = new TH1I(Form("HCalSamplesPerEventLink_%d", i),
                           Form("Link %d", i),
                           HCAL_NUM_SAMPLES_PER_EVENT + 1, -0.5, HCAL_NUM_SAMPLES_PER_EVENT + 0.5);

    mHCalSamplesPerEvent->Add(graph);
    mHCalSamplesPerEventContainer[i] = graph;
  }

  mHCalSamplesPerEventCanvas->cd(1);
  mHCalSamplesPerEvent->Draw("plc nostack");
  gPad->SetLogy();
  legend = gPad->BuildLegend();
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  getObjectsManager()->startPublishing(mHCalSamplesPerEventCanvas);


  /////////////////////////////////////////////////////////////////////////////////////
  /// Heatmap plots

  mHCalHeatmapCanvas = new TCanvas("HCalHeatmapCanvas", "HCal Heatmap Canvas", 1920, 1080);
  mHCalHeatmapCanvas->DivideSquare(HCAL_NUM_SAMPLES_PER_EVENT, 0.001, 0.001);
  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    mHCalHeatmapCanvas->cd(s+1);
    TH2D* graph = new TH2D(Form("HCalHeatmapSample_%d", s),
                           Form("Sample %d", s),
                           16, -0.5, 16 - 0.5,  // columns
                           12, -0.5, 12 - 0.2); // rows

    graph->SetStats(0);
    graph->Draw("COLZ");
    gPad->SetLogz();
    mHCalHeatmapContainer[s] = graph;
  }

  getObjectsManager()->startPublishing(mHCalHeatmapCanvas);


  /////////////////////////////////////////////////////////////////////////////////////
  /// Channel-wise ADC, per ROC, plots

  for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
    for (int j = 0; j < HCAL_NUM_ROCS_PER_LINK; ++j) {
      mHCalROCADC[i][j] = new TH2D(Form("HCalADC_Link_%d:%d", i, j),
                                        Form("ADC Per Channel (Link %d:%d)", i, j),
                                        2*HCAL_NUM_CHANNELS_PER_ROC_HALF, -0.5, 2*HCAL_NUM_CHANNELS_PER_ROC_HALF -0.5,
                                        256, 0., 1024);

      mHCalROCTOT[i][j] = new TH2D(Form("HCalTOT_Link_%d:%d", i, j),
                                        Form("TOT Per Channel (Link %d:%d)", i, j),
                                        2*HCAL_NUM_CHANNELS_PER_ROC_HALF, -0.5, 2*HCAL_NUM_CHANNELS_PER_ROC_HALF -0.5,
                                        256, 0., 1024);

      mHCalROCTOA[i][j] = new TH2D(Form("HCalTOA_Link_%d:%d", i, j),
                                        Form("TOA Per Channel (Link %d:%d)", i, j),
                                        2*HCAL_NUM_CHANNELS_PER_ROC_HALF, -0.5, 2*HCAL_NUM_CHANNELS_PER_ROC_HALF -0.5,
                                        256, 0., 1024);

      mHCalROCADC[i][j]->GetXaxis()->SetTitle("Channel Index");
      mHCalROCTOT[i][j]->GetXaxis()->SetTitle("Channel Index");
      mHCalROCTOA[i][j]->GetXaxis()->SetTitle("Channel Index");
      
      mHCalROCADC[i][j]->GetYaxis()->SetTitle("ADC");
      mHCalROCTOT[i][j]->GetYaxis()->SetTitle("TOT (Coded)");
      mHCalROCTOA[i][j]->GetYaxis()->SetTitle("TOA (Coded)");
  
      getObjectsManager()->startPublishing(mHCalROCADC[i][j]);
      getObjectsManager()->startPublishing(mHCalROCTOT[i][j]);
      getObjectsManager()->startPublishing(mHCalROCTOA[i][j]);
    }
  }


  /////////////////////////////////////////////////////////////////////////////////////
  /// Channel-wise ADC (waveform plots)

  for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
    for (int j = 0; j < HCAL_NUM_ROCS_PER_LINK; ++j) {
      for (int k = 0; k < 2; ++k) {
        TCanvas* currentCanvas = new TCanvas(Form("WaveformsCanvas_Link_%d_ROC_%d.%d", i, j, k), 
                                             Form("Link %d:%d.%d", i, j, k), 
                                             1920, 1080);

        currentCanvas->DivideSquare(HCAL_NUM_CHANNELS_PER_ROC_HALF, 0.001, 0.001);
        mHCalWaveforms[i][j][k] = currentCanvas;
        for (int chn = 0; chn < HCAL_NUM_CHANNELS_PER_ROC_HALF; ++chn) {
          currentCanvas->cd(chn+1);
          TH2D* graph = new TH2D(Form("HCalADC_Link_%d_ROC_%d_Half_%d_Chn_%d", i, j, k, chn),
                                 Form("%d:%d.%d/%d)", i, j, k, chn),
                                 HCAL_NUM_SAMPLES_PER_EVENT, 0, HCAL_NUM_SAMPLES_PER_EVENT - 1,
                                 128, 0., 1024.);
          graph->SetStats(0);
          graph->GetXaxis()->SetTitle("Sample");
          graph->GetYaxis()->SetTitle("ADC");
          gPad->SetLogz();
          graph->Draw();

          TPaveText* info = new TPaveText(0.55, 0.75, 0.9, 0.9, "NDC");
          info->SetTextSize(0.08);
          info->SetFillColor(0);
          info->SetFillStyle(4000);
          info->SetBorderSize(0);
          info->AddText(Form("CHANNEL %02d", chn));
          info->Draw();

          mHCalWaveformsContainer[i][j][k][chn] = graph;
        }

      currentCanvas->Modified();
      currentCanvas->Update();
      getObjectsManager()->startPublishing(currentCanvas);
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////////////
  /// Data error histogram

  // TODO: set up some enum to hold error types to indices or something
  int numChips = 2 * HCAL_NUM_GBT_LINKS * HCAL_NUM_ROCS_PER_LINK;
  mHCalDataErrors = new TH2I("HCalDataErrors", "HCal Data Errors",
                             numChips, -0.5, numChips - 0.5,  // count errors for each chip
                             5, -0.5, 5 - 0.5);               // number of types of errors

  mHCalDataErrors->GetXaxis()->SetTitle("Data Link (GBTLink:ROC.Half)");
  mHCalDataErrors->GetYaxis()->SetBinLabel(1, "CRC_CHECKSUM_MISMATCH");
  mHCalDataErrors->GetYaxis()->SetBinLabel(2, "DAQH_HEADER_MISMATCH");
  mHCalDataErrors->GetYaxis()->SetBinLabel(3, "DAQH_TRAILER_MISMATCH");
  mHCalDataErrors->GetYaxis()->SetBinLabel(4, "HAMMING_BITS_SET");
  mHCalDataErrors->GetYaxis()->SetBinLabel(5, "DECODER_ERROR");
  mHCalDataErrors->SetStats(0);

  int c = 0;
  for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
    for (int j = 0; j < HCAL_NUM_ROCS_PER_LINK; ++j) {
      for (int k = 0; k < 2; ++k) {
        mHCalDataErrors->GetXaxis()->SetBinLabel(++c, Form("%d:%d.%d", i, j, k));
      }
    }
  }

  getObjectsManager()->startPublishing(mHCalDataErrors);


  /////////////////////////////////////////////////////////////////////////////////////
  /// ADC Saturation plots

  mHCalADCSaturationCanvas = new TCanvas("HCalADCSaturation", "# Channels with ADC > 1015", 1920, 1080);
  mHCalADCSaturationStack  = new THStack("HCalADCSaturation", "# Channels with ADC > 1015");
  gStyle->SetPalette(kCMYK);
  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    TH1I* graph = new TH1I(Form("HCalADCSaturationSample%d", s),
                           Form("Sample %d", s),
                           16, 0., 0.);

    mHCalADCSaturationStack->Add(graph);
    mHCalADCSaturation[s] = graph;
  }

  mHCalADCSaturationCanvas->cd(1);
  mHCalADCSaturationStack->Draw("plc nostack");
  legend = gPad->BuildLegend();
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  getObjectsManager()->startPublishing(mHCalADCSaturationCanvas);
  

  /////////////////////////////////////////////////////////////////////////////////////
  /// TOT Saturation plots

  mHCalTOTSaturationCanvas = new TCanvas("HCalTOTSaturation", "# Channels with TOT > 1015", 1920, 1080);
  mHCalTOTSaturationStack  = new THStack("HCalTOTSaturation", "# Channels with TOT > 1015");
  gStyle->SetPalette(kCMYK);
  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    TH1I* graph = new TH1I(Form("HCalTOTSaturationSample%d", s),
                           Form("Sample %d", s),
                           16, 0., 0.);

    mHCalTOTSaturationStack->Add(graph);
    mHCalTOTSaturation[s] = graph;
  }

  mHCalTOTSaturationCanvas->cd(1);
  mHCalTOTSaturationStack->Draw("plc nostack");
  legend = gPad->BuildLegend();
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  getObjectsManager()->startPublishing(mHCalTOTSaturationCanvas);


  /////////////////////////////////////////////////////////////////////////////////////
  /// TOT Zero plots

  mHCalTOTZeroCanvas = new TCanvas("HCalTOTZero", "# Channels with TOT = 0", 1920, 1080);
  mHCalTOTZeroStack  = new THStack("HCalTOTZero", "# Channels with TOT = 0");
  gStyle->SetPalette(kCMYK);
  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    TH1I* graph = new TH1I(Form("HCalTOTZeroSample%d", s),
                           Form("Sample %d", s),
                           16, 0., 0.);

    mHCalTOTZeroStack->Add(graph);
    mHCalTOTZero[s] = graph;
  }

  mHCalTOTZeroCanvas->cd(1);
  mHCalTOTZeroStack->Draw("plc nostack");
  legend = gPad->BuildLegend();
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  getObjectsManager()->startPublishing(mHCalTOTZeroCanvas);


  /////////////////////////////////////////////////////////////////////////////////////
  /// TOT Below Half plots

  mHCalTOTBelowHalfCanvas = new TCanvas("HCalTOTBelowHalf", "# Channels with TOT < 512", 1920, 1080);
  mHCalTOTBelowHalfStack  = new THStack("HCalTOTBelowHalf", "# Channels with TOT < 512");
  gStyle->SetPalette(kCMYK);
  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    TH1I* graph = new TH1I(Form("HCalTOTBelowHalfSample%d", s),
                           Form("Sample %d", s),
                           16, 0., 0.);

    mHCalTOTBelowHalfStack->Add(graph);
    mHCalTOTBelowHalf[s] = graph;
  }

  mHCalTOTBelowHalfCanvas->cd(1);
  mHCalTOTBelowHalfStack->Draw("plc nostack");
  legend = gPad->BuildLegend();
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  getObjectsManager()->startPublishing(mHCalTOTBelowHalfCanvas);

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

  if (isLostTimeframe(ctx)) {
    mTFerrorCounter->Fill(1);
    return;
  } else {
    mTFerrorCounter->Fill(2);
  }
    
  std::vector<char> recordBuffer;
  o2::InteractionRecord currentInteractionRecord;
  for (const auto& rawData : framework::InputRecordWalker(ctx.inputs())) {
    recordBuffer.clear();
    if (rawData.header != nullptr && rawData.payload != nullptr) {
      const auto payloadSize = o2::framework::DataRefUtils::getPayloadSize(rawData); 
      mPayloadSizeTF->Fill(payloadSize);

      gsl::span<const char> databuffer(rawData.payload, payloadSize);

      long int currentpos = 0;
      while (currentpos < databuffer.size()) {
        auto rdh             = reinterpret_cast<const o2::header::RDHAny*>(databuffer.data() + currentpos);

        auto pageSize        = o2::raw::RDHUtils::getMemorySize(rdh);
        auto pageHeaderSize  = o2::raw::RDHUtils::getHeaderSize(rdh);
        auto pagePayloadSize = pageSize - pageHeaderSize;

        auto feeID           = o2::raw::RDHUtils::getFEEID(rdh);
        auto triggerBC       = o2::raw::RDHUtils::getTriggerBC(rdh);
        auto triggerOrbit    = o2::raw::RDHUtils::getTriggerOrbit(rdh);
        auto offset          = o2::raw::RDHUtils::getOffsetToNext(rdh);
        auto stop            = o2::raw::RDHUtils::getStop(rdh);
        
        currentInteractionRecord.bc    = triggerBC;
        currentInteractionRecord.orbit = triggerOrbit;
        
        LOGF(debug, "---- New HBF ----");

        LOGF(debug, "FEE: %04x",          feeID        );
        LOGF(debug, "BC: %d",             triggerBC    );
        LOGF(debug, "Orbit: %d",          triggerOrbit );
        LOGF(debug, "Offset to next: %d", offset       );
        LOGF(debug, "Stop bit: %d",       stop         );

        auto pagePayload = databuffer.subspan(currentpos + pageHeaderSize, pagePayloadSize);
        std::copy(pagePayload.begin(), pagePayload.end(), std::back_inserter(recordBuffer));
        LOGF(debug, "Record buffer now contains %d bytes of payload data", recordBuffer.size());

        auto triggerType = o2::raw::RDHUtils::getTriggerType(rdh);
        if (triggerType & o2::trigger::HB) {
          if (stop) {
            mNumHBFPerCRU->Fill(o2::raw::RDHUtils::getCRUID(rdh));
            mCRUcounter->Fill(o2::raw::RDHUtils::getEndPointID(rdh), o2::raw::RDHUtils::getLinkID(rdh));
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

// TODO: clean up this mess of a function
// (it works well enough though!)
void HCalTestbeamTask::processHCalEvent(const gsl::span<const char> hcalpayload)
{
  mDecoder.reset();
  mDecoder.decodeBuffer(hcalpayload);
  if (not mDecoder.hasEventData()) {
    LOGF(debug, "Event hasEventData = false; skipping");
    return;
  }

  if (not mDecoder.isDataValid()) {
    LOGF(debug, "Event isDataValid = false; skipping");
    for (int i = 0; i < HCAL_NUM_GBT_LINKS * HCAL_NUM_ROCS_PER_LINK * 2; ++i) {
      mHCalDataErrors->Fill(i, 4., 1.);
    }

    return;
  }

  ++mNumEvents;
  LOGF(debug, "events: %d", mNumEvents);
  for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
    int numSamples = mDecoder.getNumSamplesRead(i);
    LOGF(debug, "Samples read from link %d: %d)", i, numSamples);
    mHCalSamplesPerEventContainer[i]->Fill(numSamples, 1);
  }

  std::array<std::array<o2::focal::HCalGBTLink, HCAL_NUM_GBT_LINKS>, HCAL_NUM_SAMPLES_PER_EVENT> links = mDecoder.getData();

  int currentChannel = 0;
  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    int globalADCsum = 0;
    int numChannelsSaturatedADC = 0;
    int numChannelsSaturatedTOT = 0;
    int numChannelsZeroTOT      = 0;
    int numChannelsBelowHalfTOT = 0;

    for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
      o2::focal::HCalGBTLink currentLink = links[s][i];

      for (int j = 0; j < HCAL_NUM_ROCS_PER_LINK; ++j) {
        o2::focal::HCalROC currentROC = currentLink.getROC(j);

        for (int k = 0; k < 2; ++k) {
          o2::focal::HCalROCDataLink currentHalf = currentROC.getChipHalf(k);
          bool errored = false;
          if (not checkCRC(currentHalf)) {
            errored = true;
            LOGF(debug, "CRC mismatch in link %d:%d.%d", i, j, k);
            int binIndex = i * HCAL_NUM_GBT_LINKS * HCAL_NUM_ROCS_PER_LINK 
                         + j * HCAL_NUM_ROCS_PER_LINK 
                         + k;

            mHCalDataErrors->Fill(binIndex, 0., 1.);
          }
          if (not checkDAQHHeader(currentHalf)) {
            errored = true;
            LOGF(debug, "DAQH header mismatch in link %d:%d.%d", i, j, k);
            int binIndex = i * HCAL_NUM_GBT_LINKS * HCAL_NUM_ROCS_PER_LINK 
                         + j * HCAL_NUM_ROCS_PER_LINK 
                         + k;

            mHCalDataErrors->Fill(binIndex, 1., 1.);
          }
          if (not checkDAQHTrailer(currentHalf)) {
            errored = true;
            LOGF(debug, "DAQH trailer mismatch in link %d:%d.%d", i, j, k);
            int binIndex = i * HCAL_NUM_GBT_LINKS * HCAL_NUM_ROCS_PER_LINK 
                         + j * HCAL_NUM_ROCS_PER_LINK 
                         + k;

            mHCalDataErrors->Fill(binIndex, 2., 1.);
          }
          if (not checkHammingBits(currentHalf)) {
            errored = true;
            LOGF(debug, "Hamming decode error bits are set in link %d:%d.%d", i, j, k);
            int binIndex = i * HCAL_NUM_GBT_LINKS * HCAL_NUM_ROCS_PER_LINK 
                         + j * HCAL_NUM_ROCS_PER_LINK 
                         + k;

            mHCalDataErrors->Fill(binIndex, 3., 1.);
          }

          //if (errored) {
          //  auto words = currentHalf.getWords();
          //  for (const auto& w : words) {
          //    LOGF(debug, "%08X", w);
          //  }
          //}

          for (int chn = 0; chn < HCAL_NUM_CHANNELS_PER_ROC_HALF; ++chn) {
            o2::focal::HCalChannel currentChannel = currentHalf.getChannel(chn);

            int adc = currentChannel.adc;
            int tot = currentChannel.tot;
            int toa = currentChannel.toa;
            
            // TODO: decide on how to fill these plots, i.e. which sample?
            //if (s == 3) {
              mHCalROCADC[i][j]->Fill(chn + HCAL_NUM_CHANNELS_PER_ROC_HALF * k, adc);
              mHCalROCTOT[i][j]->Fill(chn + HCAL_NUM_CHANNELS_PER_ROC_HALF * k, tot);
              mHCalROCTOA[i][j]->Fill(chn + HCAL_NUM_CHANNELS_PER_ROC_HALF * k, toa);
            //}     

            globalADCsum += adc;
            mHCalWaveformsContainer[i][j][k][chn]->Fill(s, adc);

            if (adc > 1015) {
              ++numChannelsSaturatedADC;
            }

            if (tot > 1015) {
              ++numChannelsSaturatedTOT;
            }
            else if (tot == 0) {
              ++numChannelsZeroTOT;
            }
            else if (tot < 512) {
              ++numChannelsBelowHalfTOT;
            }

            std::pair<int, int> coords = mMapper.getRowCol(i, j, k, chn);
            if ((coords.first == -1) || (coords.second == -1)) {
              continue;
            } else {

              // Instead of simply filling histograms with ADC values,
              // we want to show the "average" ADC value on the heatmaps
              double previous = mHCalHeatmapContainer[s]
                                ->GetBinContent(coords.second+1, coords.first+1) 
                                * mNumEvents; // "unaveraged" value from current cell

              mHCalHeatmapContainer[s]->SetBinContent(coords.second+1, 
                                                      coords.first+1, 
                                                      (previous + adc) / (mNumEvents + 1)); // new average value
            }
          }
        }
      }
    }

    mHCalGlobalADCSumContainer[s] ->Fill(globalADCsum);    
    mHCalADCSaturation[s]         ->Fill(numChannelsSaturatedADC);
    mHCalTOTSaturation[s]         ->Fill(numChannelsSaturatedTOT);
    mHCalTOTZero[s]               ->Fill(numChannelsZeroTOT);
    mHCalTOTBelowHalf[s]          ->Fill(numChannelsBelowHalfTOT);
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

  mTFerrorCounter     ->Reset();
  mFEENumberHBF       ->Reset();
  mFEENumberTF        ->Reset();
  mNumLinksTF         ->Reset();
  mNumHBFPerCRU       ->Reset();
  mCRUcounter         ->Reset();
  mPayloadSizeTF      ->Reset();
  mPayloadSizeHCalGBT ->Reset();
  
  for (int i = 0; i < HCAL_NUM_GBT_LINKS; ++i) {
    mHCalSamplesPerEventContainer[i]->Reset();

    for (int j = 0; j < HCAL_NUM_ROCS_PER_LINK; ++j) {
      mHCalROCADC[i][j]->Reset();
      mHCalROCTOT[i][j]->Reset();
      mHCalROCTOA[i][j]->Reset();

      for (int k = 0; k < 2; ++k) {
        for (int chn = 0; chn < HCAL_NUM_CHANNELS_PER_ROC_HALF; ++chn) {
          mHCalWaveformsContainer[i][j][k][chn]->Reset();
        }
      }
    }
  }

  for (int s = 0; s < HCAL_NUM_SAMPLES_PER_EVENT; ++s) {
    mHCalGlobalADCSumContainer[s] ->Reset();
    mHCalHeatmapContainer[s]      ->Reset();
    mHCalADCSaturation[s]         ->Reset();
    mHCalTOTSaturation[s]         ->Reset();
    mHCalTOTZero[s]               ->Reset();
    mHCalTOTBelowHalf[s]          ->Reset();
  }
}

bool HCalTestbeamTask::isLostTimeframe(framework::ProcessingContext& ctx) const
{
  // direct data
  constexpr auto originFOC = header::gDataOriginFOC;
  o2::framework::InputSpec dummy{"dummy", 
                                 framework::ConcreteDataMatcher{originFOC,
                                                                header::gDataDescriptionRawData,
                                                                0xDEADBEEF} };

  for (const auto& ref : o2::framework::InputRecordWalker(ctx.inputs(), { dummy })) {
    const auto dh = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
    const auto payloadSize = o2::framework::DataRefUtils::getPayloadSize(ref);
    if (payloadSize == 0) {
      return true;
    }
  }

  // sampled data
  o2::framework::InputSpec dummyDS{"dummyDS",
                                   framework::ConcreteDataMatcher{"DS",
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
  std::array<unsigned int, HCAL_NUM_GBT_LINES_PER_LINK> words = half.getWords();
  unsigned int crcWord = words[HCAL_NUM_GBT_LINES_PER_LINK - 1];

  boost::crc_basic<32> crc_32(0x04C11DB7, 0x00000000, 0x00000000, false, false);
  for (int i = 0; i < HCAL_NUM_GBT_LINES_PER_LINK - 1; ++i) {
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

