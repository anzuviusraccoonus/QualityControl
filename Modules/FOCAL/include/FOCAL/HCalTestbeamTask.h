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
#include "THStack.h"
class TH1;
class TH2;
class TProfile2D;

using namespace o2::quality_control::core;

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
    static constexpr int mNUM_ASICS = 2;
    static constexpr int mCHANNELS_PER_ASIC = 76;
    static constexpr int mRANGE_ADC = 1024;
    static constexpr int mRANGE_TOA = 1024;
    static constexpr int mRANGE_TOT = 1024;

  bool isLostTimeframe(framework::ProcessingContext& ctx) const;
  void processHCalEvent(const gsl::span<const char> gbtpayload);

  o2::focal::HCalDecoder mDecoder;
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
  
  TH1* mPayloadSizeHCalGBT = nullptr;
  //TH1* mHCalGlobalADCSum   = nullptr;

  std::array<TH2*, mNUM_ASICS> mHCalASICChannelADC;
  std::array<TH2*, mNUM_ASICS> mHCalASICChannelTOA;
  std::array<TH2*, mNUM_ASICS> mHCalASICChannelTOT;

  //TH2* mHCalWaveforms[2][2][2][36];
  TCanvas* mHCalWaveforms[2][2][2];
  TH2* mHCalWaveformsContainer[2][2][2][36];

  TCanvas* mHCalGlobalADCSumCanvas;
  THStack* mHCalGlobalADCSum;
  TH1* mHCalGlobalADCSumContainer[16];

};

} // namespace o2::quality_control_modules::focal

#endif // QC_MODULE_FOCAL_FOCALHCALTESTBEAMTASK_H

