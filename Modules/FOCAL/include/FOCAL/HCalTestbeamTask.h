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
/// \file   HCalTestbeamTask.h
/// \author Ian Pascal Møller <ian.moeller@nbi.ku.dk>
/// \brief Standalone QC task for FoCal-H beam tests
///

#ifndef QC_MODULE_FOCAL_FOCALHCALTESTBEAMTASK_H
#define QC_MODULE_FOCAL_FOCALHCALTESTBEAMTASK_H

#include <array>
#include <unordered_map>

#include "QualityControl/TaskInterface.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTReconstruction/GBTWord.h"
#include "FOCALReconstruction/HCalDecoder.h"
#include "FOCALReconstruction/HCalDataWord.h"
#include "FOCALReconstruction/HCalGBTLink.h"
#include "FOCALReconstruction/HCalMapper.h"

#include "THStack.h"

class TH1;
class TH2;
class TProfile2D;

using namespace o2::quality_control::core;
using namespace o2::focal::constants;

namespace o2::quality_control_modules::focal
{

class HCalTestbeamTask final : public TaskInterface 
{
  public:
    HCalTestbeamTask() = default;
    ~HCalTestbeamTask() override;
    
    void initialize(o2::framework::InitContext& ctx) override;
    void startOfActivity(const Activity& activity) override;
    void startOfCycle() override;
    void monitorData(o2::framework::ProcessingContext& ctx) override;
    void endOfCycle() override;
    void endOfActivity(const Activity& activity) override;
    void reset() override;

  private:
    int mNumEvents = 0; // internal counter for re-calculating average values for heatmap

  bool checkDAQHHeader(o2::focal::HCalROCDataLink half);
  bool checkDAQHTrailer(o2::focal::HCalROCDataLink half);
  bool checkHammingBits(o2::focal::HCalROCDataLink half);
  bool checkCRC(o2::focal::HCalROCDataLink half);
  bool isLostTimeframe(framework::ProcessingContext& ctx) const;
  void processHCalEvent(const gsl::span<const char> gbtpayload);

  o2::focal::HCalDecoder mDecoder;
  o2::focal::HCalMapper  mMapper;

  bool mDebugMode = false;

  /////////////////////////////////////////////////////////////////////////////////////
  /// General histograms
  /////////////////////////////////////////////////////////////////////////////////////

  TH1* mTFerrorCounter  = nullptr; ///< Number of TF builder errors
  TH1* mFEENumberHBF    = nullptr; ///< Number of HBFs per FEE
  TH1* mFEENumberTF     = nullptr; ///< Number of TFs per FEE
  TH1* mNumLinksTF      = nullptr; ///< Number of links per timeframe
  TH1* mNumHBFPerCRU    = nullptr; ///< Number of HBFs per CRU
  TH2* mCRUcounter      = nullptr; ///< CRU counter
  TH1* mPayloadSizeTF   = nullptr; ///< Payload size per timeframe
  

  /////////////////////////////////////////////////////////////////////////////////////
  /// HCal-specific histograms
  /////////////////////////////////////////////////////////////////////////////////////
  
  TH1*      mPayloadSizeHCalGBT = nullptr;

  TH2*      mHCalROCADC[HCAL_NUM_GBT_LINKS][HCAL_NUM_ROCS_PER_LINK];
  TH2*      mHCalROCTOT[HCAL_NUM_GBT_LINKS][HCAL_NUM_ROCS_PER_LINK];
  TH2*      mHCalROCTOA[HCAL_NUM_GBT_LINKS][HCAL_NUM_ROCS_PER_LINK];

  TCanvas*  mHCalWaveforms         [HCAL_NUM_GBT_LINKS][HCAL_NUM_ROCS_PER_LINK][2];
  TH2*      mHCalWaveformsContainer[HCAL_NUM_GBT_LINKS][HCAL_NUM_ROCS_PER_LINK][2][HCAL_NUM_CHANNELS_PER_ROC_HALF];

  TCanvas*  mHCalGlobalADCSumCanvas;
  THStack*  mHCalGlobalADCSum;
  TH1*      mHCalGlobalADCSumContainer[HCAL_NUM_SAMPLES_PER_EVENT];

  TCanvas*  mHCalSamplesPerEventCanvas;
  THStack*  mHCalSamplesPerEvent;
  TH1*      mHCalSamplesPerEventContainer[HCAL_NUM_GBT_LINKS];

  TCanvas*  mHCalADCHeatmapCanvas;
  TH2*      mHCalADCHeatmapContainer[HCAL_NUM_SAMPLES_PER_EVENT];
  
  TCanvas*  mHCalTOTHeatmapCanvas;
  TH2*      mHCalTOTHeatmapContainer[HCAL_NUM_SAMPLES_PER_EVENT];
  
  TCanvas*  mHCalTOAHeatmapCanvas;
  TH2*      mHCalTOAHeatmapContainer[HCAL_NUM_SAMPLES_PER_EVENT];

  TH2*      mHCalDataErrors;

  TCanvas*  mHCalADCSaturationCanvas;
  THStack*  mHCalADCSaturationStack;
  TH1*      mHCalADCSaturation[HCAL_NUM_SAMPLES_PER_EVENT];

  TCanvas*  mHCalTOTSaturationCanvas;
  THStack*  mHCalTOTSaturationStack;
  TH1*      mHCalTOTSaturation[HCAL_NUM_SAMPLES_PER_EVENT];

  TCanvas*  mHCalTOTZeroCanvas;
  THStack*  mHCalTOTZeroStack;
  TH1*      mHCalTOTZero[HCAL_NUM_SAMPLES_PER_EVENT];

  TCanvas*  mHCalTOTBelowHalfCanvas;
  THStack*  mHCalTOTBelowHalfStack;
  TH1*      mHCalTOTBelowHalf[HCAL_NUM_SAMPLES_PER_EVENT];

  TH2*      mHCalADCvsTOT;
  TH2*      mHCalADCvsTOA;
  TH2*      mHCalTOTvsTOA;
  
  TH1*      mHCalNumTOTValues;
  TH1*      mHCalNumTOAValues;
  TH1*      mHCalSampleTOTIdx;
  TH1*      mHCalSampleTOAIdx;

};

} // namespace o2::quality_control_modules::focal

#endif // QC_MODULE_FOCAL_FOCALHCALTESTBEAMTASK_H
