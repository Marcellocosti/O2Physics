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

/// \file correlatorFlowCharmHadrons.cxx
/// \brief CharmHadrons-Hadrons correlator tree creator for data and MC-reco analyses
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, CERN

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

template <typename T>
concept IsTrack = requires(T candidate) {
  candidate.originTrackId();
};

template<IsTrack TTrack>
double getPt(const TTrack& track) {
  return track.ptAssocTrack();
}

template<typename TCand>
double getPt(const TCand& candidate) {
  return candidate.ptCand();
}

template<IsTrack TTrack>
double getEta(const TTrack& track) {
  return track.etaAssocTrack();
}

template<typename TCand>
double getEta(const TCand& candidate) {
  return candidate.etaCand();
}

template<typename T1, typename T2>
double getDeltaEta(const T1& cand1, const T2& cand2) {
  return getEta(cand1) - getEta(cand2);
}

template<IsTrack TTrack>
double getPhi(const TTrack& track) {
  return track.phiAssocTrack();
}

template<typename TCand>
double getPhi(const TCand& candidate) {
  return candidate.phiCand();
}

template<typename T1, typename T2>
double getDeltaPhi(const T1& cand1, const T2& cand2) {
  return getPhi(cand1) - getPhi(cand2);
}

enum PoolBinningStrategy {
  Centrality = 0,
  Multiplicity
};

using BinningTypeDerivedCent = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality>;
using BinningTypeDerivedMult = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity>;

struct HfCorrelatorFlowCharmHadrons {
  Produces<aod::HfcRedChHadPairs> entryCharmHadPair;
  Produces<aod::HfcRedCharmRecs> entryCharmRecoInfo;

  Configurable<bool> fillHistoData{"fillHistoData", true, "Flag for filling histograms in data processes"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{1., 3., 5., 8., 16., 36.}, "pT bin limits for candidate mass plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 1., 2., 50.}, "pT bin limits for assoc particle"};

  SliceCache cache;

  Preslice<aod::HfcRedTrkAssocs> tracksPerCol = aod::hf_candidate_reduced::hfcRedFlowCollId;
  Preslice<aod::HfcRedCharmHads> candPerCol = aod::hf_candidate_reduced::hfcRedFlowCollId;

  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "Z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "Event multiplicity pools (FT0M)"};
  ConfigurableAxis centPoolBins{"centPoolBins", {VARIABLE_WIDTH, 0., 10., 20., 30.}, "Event centrality pools"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {100, 0., 10000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsCent{"binsCent", {100, 0., 100.}, "Centrality bins"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "Primary vertex z coordinate"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "Eta bins"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "Phi bins"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    std::array<bool, 8> doprocess{doprocessSameEventCharmHadWCentMix, doprocessSameEventCharmHadWMultMix, doprocessMixedEventCharmHadWCentMix, doprocessMixedEventCharmHadWMultMix,
                                  doprocessSameEventHadHadWCentMix, doprocessSameEventHadHadWMultMix, doprocessMixedEventHadHadWCentMix, doprocessMixedEventHadHadWMultMix};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
      if ( !((doprocessSameEventCharmHadWCentMix && doprocessMixedEventCharmHadWCentMix) || (doprocessSameEventCharmHadWMultMix && doprocessMixedEventCharmHadWMultMix)) ){
        LOG(fatal) << "Different binning policies between Same Event and Mixed Event";
      }
      if ( !((doprocessSameEventHadHadWCentMix && doprocessMixedEventHadHadWCentMix) || (doprocessSameEventHadHadWMultMix && doprocessMixedEventHadHadWMultMix)) ){
        LOG(fatal) << "Different binning policies between Same Event and Mixed Event";
      }
    }

    AxisSpec axisCent = {binsCent, "Centrality"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec axisEta = {binsEta, "#it{#eta}"};
    AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T} Hadron (GeV/#it{c})"};

    // Histograms for data analysis
    if (fillHistoData) {
      if (doprocessSameEventCharmHadWCentMix || doprocessMixedEventCharmHadWCentMix || doprocessSameEventHadHadWCentMix || doprocessMixedEventHadHadWCentMix) {
        registry.add("hCent", "Centrality", {HistType::kTH1F, {axisCent}});
      } else {
        registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH1F, {axisMultFT0M}});
      }
      registry.add("hZVtx", "z vertex", {HistType::kTH1F, {axisPosZ}});
      registry.add("hCollisionPoolBin", "Collision pool bin", {HistType::kTH1F, {axisPoolBin}});
      registry.add("hCharmPoolBin", "Charm candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
      registry.add("hPhiVsPtCand", "Charm candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtD}}});
      registry.add("hPhiVsPtPartAssoc", "Associated particles phiVsPt", {HistType::kTH3F, {{axisPhi}, {axisPtD}, {axisPtHadron}}});
      registry.add("hEtaVsPtCand", "Charm candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtD}}});
      registry.add("hEtaVsPtPartAssoc", "Associated particles etaVsPt", {HistType::kTH3F, {{axisEta}, {axisPtD}, {axisPtHadron}}});
      registry.add("hTracksPoolBin", "Associated particles pool bin", {HistType::kTH1F, {axisPoolBin}});
    }
  }

  template<typename TTrigParts, typename TAssocParts, typename TBinningType>
  void fillSameEvent(aod::HfcRedFlowColls const& collisions,
                     TTrigParts const& trigCands,
                     TAssocParts const& assocTracks,
                     TBinningType corrBinning)
  {
    for (const auto& collision : collisions) {
      int poolBin{0};
      if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedCent>) {
        poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.centrality()));
        registry.fill(HIST("hCent"), collision.centrality());
      } else if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedMult>) {
        poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multiplicity()));
        registry.fill(HIST("hMultiplicityFT0M"), collision.multiplicity());
      }
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      registry.fill(HIST("hZVtx"), collision.posZ());
      
      auto thisCollId = collision.globalIndex();
      if constexpr (std::is_same_v<TTrigParts, soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels>>) {
        TTrigParts candsTrig = trigCands.sliceBy(tracksPerCol, thisCollId);
        auto candsAssoc = assocTracks.sliceBy(tracksPerCol, thisCollId);

        for (const auto& candidate : candsTrig) {
          registry.fill(HIST("hCharmPoolBin"), poolBin);
          registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(getPhi(candidate), -o2::constants::math::PIHalf), getPt(candidate));
          registry.fill(HIST("hEtaVsPtCand"), getEta(candidate), getPt(candidate));
          for (const auto& track : candsAssoc) {
            // Removing Ds daughters by checking track indices
            if constexpr (std::is_same_v<TTrigParts, soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels>>) {
              if (candidate.originTrackId() == track.originTrackId()) {
                continue;
              }
            } else {
              if ((candidate.prong0Id() == track.originTrackId()) || (candidate.prong1Id() == track.originTrackId()) || (candidate.prong2Id() == track.originTrackId())) {
                continue;
              }
            }
            registry.fill(HIST("hTracksPoolBin"), poolBin);
            registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(getPhi(track), -o2::constants::math::PIHalf), getPt(candidate), getPt(track));
            registry.fill(HIST("hEtaVsPtPartAssoc"), getEta(track), getPt(candidate), getPt(track));

            // LOG(info) << "candidate.globalIndex(): " << candidate.globalIndex() << ", track.globalIndex(): " << track.globalIndex();
            entryCharmHadPair(candidate.globalIndex(), track.globalIndex(),
                              RecoDecay::constrainAngle(getPhi(track), getPhi(candidate), -o2::constants::math::PIHalf),
                              getEta(track) - getEta(candidate),
                              poolBin);
          }
        }
      } else {
        TTrigParts candsTrig = trigCands.sliceBy(candPerCol, thisCollId);
        auto candsAssoc = assocTracks.sliceBy(tracksPerCol, thisCollId);
  
        for (const auto& candidate : candsTrig) {
          registry.fill(HIST("hCharmPoolBin"), poolBin);
          registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(getPhi(candidate), -o2::constants::math::PIHalf), getPt(candidate));
          registry.fill(HIST("hEtaVsPtCand"), getEta(candidate), getPt(candidate));
          for (const auto& track : candsAssoc) {
            // Removing Ds daughters by checking track indices
            if constexpr (std::is_same_v<TTrigParts, soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels>>) {
              if (candidate.originTrackId() == track.originTrackId()) {
                continue;
              }
            } else {
              if ((candidate.prong0Id() == track.originTrackId()) || (candidate.prong1Id() == track.originTrackId()) || (candidate.prong2Id() == track.originTrackId())) {
                continue;
              }
            }
            registry.fill(HIST("hTracksPoolBin"), poolBin);
            registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(getPhi(track), -o2::constants::math::PIHalf), getPt(candidate), getPt(track));
            registry.fill(HIST("hEtaVsPtPartAssoc"), getEta(track), getPt(candidate), getPt(track));

            // LOG(info) << "candidate.globalIndex(): " << candidate.globalIndex() << ", track.globalIndex(): " << track.globalIndex();
            entryCharmHadPair(candidate.globalIndex(), track.globalIndex(),
                              RecoDecay::constrainAngle(getPhi(track), getPhi(candidate), -o2::constants::math::PIHalf),
                              getEta(track) - getEta(candidate),
                              poolBin);
          }
        }
      }
    }
  }

  template<typename TTrigParts, typename TAssocParts, typename TBinningType>
  void fillMixedEvent(aod::HfcRedFlowColls const& collisions,
                      TTrigParts const& trigCands,
                      TAssocParts const& assocTracks,
                      TBinningType corrBinning)
  {
    // LOG(info) << "Entered fillMixedEvent ...";
    for (const auto& collision : collisions) {
      int poolBin{0};
      if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedCent>) {
        poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.centrality()));
        registry.fill(HIST("hCent"), collision.centrality());
      } else if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedMult>) {
        poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multiplicity()));
        registry.fill(HIST("hMultiplicityFT0M"), collision.multiplicity());
      }
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      registry.fill(HIST("hZVtx"), collision.posZ());

      auto thisCollId = collision.globalIndex();
      if constexpr (std::is_same_v<TTrigParts, soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels>>) {
        TTrigParts candsTrig = trigCands.sliceBy(tracksPerCol, thisCollId);
        auto candsAssoc = assocTracks.sliceBy(tracksPerCol, thisCollId);
        for (const auto& candidate : candsTrig) {
          registry.fill(HIST("hCharmPoolBin"), poolBin);
          registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(getPhi(candidate), -o2::constants::math::PIHalf), getPt(candidate));
          registry.fill(HIST("hEtaVsPtCand"), getEta(candidate), getPt(candidate));
          for (const auto& track : candsAssoc) {
            registry.fill(HIST("hTracksPoolBin"), poolBin);
            registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(getPhi(track), -o2::constants::math::PIHalf), getPt(candidate), getPt(track));
            registry.fill(HIST("hEtaVsPtPartAssoc"), getEta(track), getPt(candidate), getPt(track));
          }
        }
      } else {
        TTrigParts candsTrig = trigCands.sliceBy(candPerCol, thisCollId);
        auto candsAssoc = assocTracks.sliceBy(tracksPerCol, thisCollId);
        for (const auto& candidate : candsTrig) {
          registry.fill(HIST("hCharmPoolBin"), poolBin);
          registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(getPhi(candidate), -o2::constants::math::PIHalf), getPt(candidate));
          registry.fill(HIST("hEtaVsPtCand"), getEta(candidate), getPt(candidate));
          for (const auto& track : candsAssoc) {
            registry.fill(HIST("hTracksPoolBin"), poolBin);
            registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(getPhi(track), -o2::constants::math::PIHalf), getPt(candidate), getPt(track));
            registry.fill(HIST("hEtaVsPtPartAssoc"), getEta(track), getPt(candidate), getPt(track));
          }
        }
      }
    }

    // LOG(info) << "About to start mixing with ncolls: " << collisions.size() << ", ntrigparts: " << trigCands.size() << ", nassocparts: " << assocTracks.size();
    auto tracksTuple = std::make_tuple(trigCands, assocTracks);
    Pair<aod::HfcRedFlowColls, TTrigParts, TAssocParts, TBinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    // LOG(info) << "Tuples created";

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      // LOG(info) << "tracks1.size(): " << tracks1.size() << ", tracks2.size(): " << tracks2.size();
      if (tracks1.size() == 0 || tracks2.size() == 0) {
        LOG(info) << "Null size, continue";
        continue;
      }
      int poolBinAssoc{0}, poolBinCharm{0};
      if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedCent>) {
        poolBinAssoc = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.centrality()));
        poolBinCharm = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.centrality()));
      } else if constexpr (std::is_same_v<TBinningType, BinningTypeDerivedMult>) {
        poolBinAssoc = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multiplicity()));
        poolBinCharm = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multiplicity()));
      }

      if (poolBinAssoc != poolBinCharm) {
        LOGF(info, "Error, poolBins are diffrent");
      }

      for (const auto& [trigCand, trackAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", trigCand.index(), trackAssoc.index(), c1.index(), c2.index(), trigCand.hfcRedFlowCollId(), trackAssoc.hfcRedFlowCollId());
        entryCharmHadPair(trigCand.globalIndex(), trackAssoc.globalIndex(),
                          RecoDecay::constrainAngle(getPhi(trackAssoc) - getPhi(trigCand), -o2::constants::math::PIHalf),
                          getEta(trackAssoc) - getEta(trigCand),
                          poolBinCharm);
        entryCharmRecoInfo(false, false);
        // entryCharmCandHadGenInfo(false, false, 0);
      }
    }
  }

  void processSameEventCharmHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                        soa::Join<aod::HfcRedCharmHads, aod::HfcRedCharmMls> const& candidates,
                                        soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    BinningTypeDerivedCent corrBinningCent{{zPoolBins, centPoolBins}, true};
    fillSameEvent(collisions, candidates, tracks, corrBinningCent);
    // fillSameEvent<PoolBinningStrategy::Centrality>(collisions, candidates, tracks, corrBinningCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processSameEventCharmHadWCentMix, "Process Same Event for Charm-Had with centrality pools", true);

  void processSameEventCharmHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                        soa::Join<aod::HfcRedCharmHads, aod::HfcRedCharmMls> const& candidates,
                                        soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    BinningTypeDerivedMult corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillSameEvent(collisions, candidates, tracks, corrBinningMult);
    // fillSameEvent<PoolBinningStrategy::Multiplicity>(collisions, candidates, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processSameEventCharmHadWMultMix, "Process Same Event for Charm-Had with multiplicity pools", false);

  void processMixedEventCharmHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                         soa::Join<aod::HfcRedCharmHads, aod::HfcRedCharmMls> const& candidates,
                                         soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    // LOG(info) << "Processing Mixed Event with centrality";
    BinningTypeDerivedCent corrBinningCent{{zPoolBins, centPoolBins}, true};
    // LOG(info) << "Binning policy set!";
    fillMixedEvent(collisions, candidates, tracks, corrBinningCent);
    // fillMixedEvent<PoolBinningStrategy::Centrality>(collisions, candidates, tracks, corrBinningCent);
    // LOG(info) << "Mixed event filled!";
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processMixedEventCharmHadWCentMix, "Process Mixed Event for Charm-Had with centrality pools", false);

  void processMixedEventCharmHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                         soa::Join<aod::HfcRedCharmHads, aod::HfcRedCharmMls> const& candidates,
                                         soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    BinningTypeDerivedMult corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillMixedEvent(collisions, candidates, tracks, corrBinningMult);
    // fillMixedEvent<PoolBinningStrategy::Multiplicity>(collisions, candidates, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processMixedEventCharmHadWMultMix, "Process Mixed Event for Charm-Had with multiplicity pools", false);

  void processSameEventHadHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                      soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    // LOG(info) << "Processing Same Event with centrality";
    BinningTypeDerivedCent corrBinningCent{{zPoolBins, centPoolBins}, true};
    // LOG(info) << "Binning policy set!";
    fillSameEvent(collisions, tracks, tracks, corrBinningCent);
    // fillSameEvent<PoolBinningStrategy::Centrality>(collisions, tracks, tracks, corrBinningCent);
    // LOG(info) << "Same event filled!";
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processSameEventHadHadWCentMix, "Process Same Event for Had-Had with centrality pools", false);

  void processSameEventHadHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                      soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    BinningTypeDerivedMult corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillSameEvent(collisions, tracks, tracks, corrBinningMult);
    // fillSameEvent<PoolBinningStrategy::Multiplicity>(collisions, tracks, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processSameEventHadHadWMultMix, "Process Same Event for Had-Had with multiplicity pools", false);

  void processMixedEventHadHadWCentMix(aod::HfcRedFlowColls const& collisions,
                                       soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    BinningTypeDerivedCent corrBinningCent{{zPoolBins, centPoolBins}, true};
    fillMixedEvent(collisions, tracks, tracks, corrBinningCent);
    // fillMixedEvent<PoolBinningStrategy::Centrality>(collisions, tracks, tracks, corrBinningCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processMixedEventHadHadWCentMix, "Process Mixed Event for Had-Had with centrality pools", false);

  void processMixedEventHadHadWMultMix(aod::HfcRedFlowColls const& collisions,
                                       soa::Join<aod::HfcRedTrkAssocs, aod::HfcRedTrkSels> const& tracks)
  {
    BinningTypeDerivedMult corrBinningMult{{zPoolBins, multPoolBins}, true};
    fillMixedEvent(collisions, tracks, tracks, corrBinningMult);
    // fillMixedEvent<PoolBinningStrategy::Multiplicity>(collisions, tracks, tracks, corrBinningMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadrons, processMixedEventHadHadWMultMix, "Process Mixed Event for Had-Had with multiplicity pools", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorFlowCharmHadrons>(cfgc)};
}
