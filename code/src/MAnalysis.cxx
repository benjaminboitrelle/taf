// @(#)maf/maf:$Name:  $:$Id: MAnalysis.cxx v.1 2005/10/02 18:03:46 sha Exp $
// Author: A. Shabetai
// Last modified: JB 2009/07/20
// Last modified: JB 2009/07/22
// Last modified: JB 2009/08/26 MimosaPro clean-up
// Last modified: JB 2009/09/01
// Last modified: JB 2009/09/08 MimosaPro investigating bad T-H matching
// Last modified: JB 2009/09/15 MimosaPro bugs correction
// Last modified: JB 2009/05/13 Cluster shape study
// Last modified: JB 2010/06/03 Binary sensor study
// Last modified: JB 2010/08/30 MiniVectors study
// Last modified: JB 2010/09/08 AlignMimosa and update the alignment parameters
// Last modified: JB 2011/06/30 GetAlignment
// Last modified: JB 2011/07/01 MimosaFake
// Last modified: JB 2011/07/07 to localize path names
// Last Modified: JB 2011/07/21 introduce fSession member pointer
// Last Modified: SS 2011/10/26 correction of the resolution calculations
// Last Modified: SS 2011/11/04 new methods Efficiency_
// Last Modified: SS 2011/11/23 new methods HotPixel_
// Last Modified: JB 2011/11/25 AlignMimosa
// Last Modified: SS 2011/12/01 Removed efficiency error calculation for the
// high efficiency case Last Modified: JB 2012/08/18 CreateNewEta with
// conditions on histo entries Last Modified: JB 2012/09/05 CreateNewEta Last
// Modified: JB 2012/09/05 all PixelSize changed to pixelSizeU or V Last
// Modified: JB 2012/10/30 new consistency check in SetPlaneNumber Last
// Modified: JB 2012/11/21 OpenInputFile, GetFileNumber Last Modified: JB
// 2012/11/21 ComputePixelPosition Last Modified: CLM 2012/11/24
// ComputePixelPosition for M32ter-P31 Last Modified: CLM 2013/01/23
// TrackHitPosition_fill Last Modified: JB 2013/05/01 ClusterPosition_compute,
// ClusterPosition_fill Last Modified: JB 2013/07/17 ClusterPosition_compute
// Last Modified: JB 2013/08/21 CreateNewEta, GetParameters
// Last Modified: JB 2013/09/12 GetParameters
// Last Modified: JB 2013/09/19 OpenInputFile
// Last Modified: JB 2013/09/29 CreateNewEta
// Last Modified: JB 2013/10/04 GetMiEta
// Last Modified: JB 2013/11/08 ClusterCharges_fill, GetParameters,
// MiniVector_fill Last Modified: JB 2014/01/21 GetParameters Last Modified: JB
// 2014/02/10 OpenInputFile, GetAlignment, IsPlaneInLadder Last Modified: VR
// 2014/06/29 Add GetRax() method Last Modified: VR 2014/06/30 Replace DTDIR by
// fWorkingDirectory in some methods Last Modified: JH 2014/07/21
// ProjectionImaging_init, ProjectionImaging_Fill Last Modified: JB 2015/01/28
// GetParameters

/////////////////////////////////////////////////////////////
//                                                         //
//  Defines some basic analysis methods,
//   which are called from the main analysis methods
//   like MimosaPro, MimosaCalibration, MimosaMiniVector, ...
//
//                                                         //
////////////////////////////////////////////////////////////

#ifndef _MimosaAnalysis_included_
#include "MAnalysis.h"
#endif
#include "DSetup.h"
#include <filesystem>

ClassImp(MimosaAnalysis)

    //______________________________________________________________________________
    //
    MimosaAnalysis *MimosaAnalysis::fgInstance =
        0; // returns pointer to global object
const std::unordered_map<std::tuple<int, int, int>, MAnalysis::ClusterType, TupleHash, TupleKeyEqual> MimosaAnalysis::clusterMap = {
    // 1 to 4 pixels
    {std::make_tuple(1, 0, 0), MAnalysis::ClusterType::one_pix},
    {std::make_tuple(2, 0, 2), MAnalysis::ClusterType::two_pix_row},     // 2 pixels in a row
    {std::make_tuple(2, 2, 0), MAnalysis::ClusterType::two_pix_col},     // 2 pixels in a column
    {std::make_tuple(3, 2, 2), MAnalysis::ClusterType::three_pix_L},     // L-shape (2x2 or 2x3)
    {std::make_tuple(3, 3, 0), MAnalysis::ClusterType::three_pix_row},   // 3 pixels in a row
    {std::make_tuple(3, 0, 3), MAnalysis::ClusterType::three_pix_col},   // 3 pixels in a column
    {std::make_tuple(4, 2, 2), MAnalysis::ClusterType::four_pix_square}, // Square 2x2
    {std::make_tuple(4, 2, 3), MAnalysis::ClusterType::four_pix_L_row},  // L shape, 2x3
    {std::make_tuple(4, 3, 2), MAnalysis::ClusterType::four_pix_L_col},  // L shape, 3x2
    {std::make_tuple(4, 0, 0), MAnalysis::ClusterType::four_others},     // Any other 4-pixel cluster
    // More than 4 pixels
    {std::make_tuple(5, 3, 2), MAnalysis::ClusterType::five_pix_squarerow}, // 5 pixels, 3x2
    {std::make_tuple(5, 2, 3), MAnalysis::ClusterType::five_pix_squarecol}, // 5 pixels, 2x3
    {std::make_tuple(5, 0, 0), MAnalysis::ClusterType::five_others},        // Any other 5-pixel cluster
    {std::make_tuple(6, 2, 3), MAnalysis::ClusterType::six_pix_3col2row},   // 6 pixels, 3x2
    {std::make_tuple(6, 3, 2), MAnalysis::ClusterType::six_pix_2col3row},   // 6 pixels, 2x3
    {std::make_tuple(6, 0, 0), MAnalysis::ClusterType::six_others},         // Any other 6-pixel cluster
    // More than 6 pixels
    {
        std::make_tuple(0, 0, 0), MAnalysis::ClusterType::more_than_six} // Any cluster larger than 6 pixels
};

MRaw *MimosaAnalysis::GetRaw()
{

  // returns a pointer to Raw data analysis methods

  if (!CheckIfDone("init"))
    return 0;
  return MRaw::InstanceRaw(fSession); // JB 2011/07/21 to pass session pointer
}

MRax *MimosaAnalysis::GetRax()
{
  // returns a pointer to Rax data analysis methods
  //  VR 2014/06/29
  if (!CheckIfDone("init"))
    return 0;
  return MRax::InstanceRax(fSession);
}

//______________________________________________________________________________
//
void MimosaAnalysis::Clear(Option_t * /*option = ""*/)
{
  // Clears all histos and canvases

  MHist::Clear();

  fClearDone = kTRUE;
  fMimosaDone = kFALSE;
  fMimosaProDone = kFALSE;
  fMimosaFakerateDone = kFALSE;
  fMimosaCalibrationDone = kFALSE;
  fMimosaMiniVectorsDone = kFALSE;
  fMimosaPro2PlanesDone = kFALSE;
  fMimosaImagingDone = kFALSE;
  fIfMCGeneration = kFALSE;
}

//______________________________________________________________________________
//
MimosaAnalysis::MimosaAnalysis()
{

  // MimosaAnalysis default constructor
  //
  // Modified: JB 2014/01/11 new flags for analysis purpose

  // ROOTDIR = gSystem->Getenv("ROOTDIR"); // suppressed by JB 2011/04/12

  m_rootFile = "";
  theCorFile = 0;

  fSession = nullptr;
  // MimoDebug = fSession->GetDebug();
  MimoDebug = 0;
  cout << endl;
  if (MimoDebug)
    cout << "MAFDebug = " << MimoDebug << endl;
  //  SettHitMapReadOpt(1);
  fWriteHistos = 0;
  fUserFileNumber = 0;
  fInitDone = kFALSE;
  fMimosaDone = kFALSE;
  fMimosaProDone = kFALSE;
  fMimosaFakerateDone = kFALSE;
  fMimosaCalibrationDone = kFALSE;
  fMimosaMiniVectorsDone = kFALSE;
  fMimosaPro2PlanesDone = kFALSE;
  fMimosaImagingDone = kFALSE;
  fIfMCGeneration = kFALSE;
  fClearDone = kTRUE;

  // Flags to control analysis purpose, JB 2014/01/11
  // Later, might be initialized with DSetup::GetAnalysisPar().AnalysisGoal
  fIfReferenceTrack = kFALSE;
  fIfCalibration = kFALSE;
  fIfLaserSpot = kFALSE;
  fIfMiniVector = kFALSE;
  fIfVertex = kFALSE;
  fIfFake = kFALSE;

  // Cuts for hit rate per pixels are initially set to zero,
  // they are re-initialized later with the HotPixel_init() or
  // SetCut_Max/MinHitRatePerPixel() methods
  CUT_MaxHitRatePerPixel = 0.; // cdritsa: set to 5; if the pixel is a seed too
                               // many times in the run, remove the hit.
  CUT_MinHitRatePerPixel =
      0.; // you can also remove pixels with low occupancy for testing

  if (fgInstance)
    Warning("MimosaAnalysis", "object already instantiated");
  else
  {
    fgInstance = this;
    // gTAF = fgInstance;
  }

  // Initializing the pointers declared in MHist.h
  saved = nullptr;
  dir = nullptr;

  hdummy = nullptr;

  //------------------------------------------------------------------------------
  //----- general control
  //------------------------------------------------------------------------------
  selection = nullptr;
  MainCanvas = nullptr;

  //------------------------------------------------------------------------------
  // Charges, noise, S/N, pixel multiplicity plots for SELECTED clusters
  //------------------------------------------------------------------------------

  // seed
  hChargeInSeed = nullptr;
  hRealTrackNoise = nullptr;
  hsn = nullptr;
  hSNReal = nullptr;

  // 2nd higest S/N pixel
  hS2N2ndRH = nullptr;
  hSNNReal = nullptr;
  hind = nullptr;

  // all pixels
  hChargeInCluster = nullptr;
  for (int i = 0; i < MaxNofPixelsInClusterMax; i++)
  {
    hqcn[i] = nullptr;
    hqcngeom[i] = nullptr;
    hindivq[i] = nullptr;
    hsnn[i] = nullptr;
  }

  hChargeOrder1 = nullptr;
  hChargeOrder2 = nullptr;
  hChargeOrder3 = nullptr;
  hChargeOrder4 = nullptr;
  hChargeOrder5 = nullptr;
  hChargeOrder6 = nullptr;
  hChargeOrder7 = nullptr;
  hChargeOrder8 = nullptr;
  hChargeOrder9 = nullptr;

  // TH2F *hChargeCor_1_2;
  // TH2F *hChargeCor_1_3 ;
  // TH2F *hChargeCor_1_4 ;
  // TH2F *hChargeCor_2_3 ;

  hChargeRap_1_over_2 = nullptr;
  hChargeRap_1_over_3 = nullptr;
  hChargeRap_1_over_4 = nullptr;
  hChargeRap_2_over_3 = nullptr;

  hQofPix3x3 = nullptr;
  hQ3x34 = nullptr;
  hQ3x35 = nullptr;
  hQ3x345 = nullptr;
  hQ3x327 = nullptr;
  hQ3x345r = nullptr;
  hQ3x327r = nullptr;

  // cluster wise
  hsnc = nullptr;
  hChargeSum_4 = nullptr;
  hqc_c = nullptr;
  hqc = nullptr;
  hChargeInCross = nullptr;
  hChargeInX = nullptr;
  hChargeInCrown1 = nullptr;
  hChargeInCrown2 = nullptr;
  hqcel = nullptr;
  hsnc1 = nullptr;
  hoptimalsize = nullptr;
  hSNneighbour = nullptr;
  hSNseedvsSNneighbour = nullptr;
  hQseedvsQcluster = nullptr;
  hQseedvsQneighbour = nullptr;
  hSNseedvsQcluster = nullptr;
  hnpix_c = nullptr;
  hnWindows_c = nullptr;
  h2DWindows = nullptr;
  hnpixCumul_c = nullptr;
  hNomEffic = nullptr;

  hDuplicate_2DMult = nullptr;
  hDuplicate_DeltaTS = nullptr;
  hDuplicate_npixc = nullptr;
  hDuplicate_npixc_vs_TrkDistToDiode = nullptr;

  static const int jpixmax = 8;
  for (int i = 0; i < jpixmax; i++)
  {
    hsn_pix_0[i] = nullptr;
    hsn_pix_1[i] = nullptr;
    hsn_seed_vs_pix_0[i] = nullptr;
    hsn_seed_vs_pix_1[i] = nullptr;
  }

  // canvas
  // TCanvas *cClusterCharge = nullptr;

  //------------------------------------------------------------------------------
  //--- Charges for calibration peak
  //------------------------------------------------------------------------------

  hqSeedCalibration = nullptr;
  hqNeighbourCalibration = nullptr;
  hqSeedVsNeighbourCalibration = nullptr;

  cCalibration = nullptr;
  Calib_distr1 = nullptr;
  Calib_distr2 = nullptr;
  Calib_distr3 = nullptr;

  //------------------------------------------------------------------------------
  //----- Comparebadgoodraw()
  //------------------------------------------------------------------------------
  hraw1badone = nullptr;
  hraw1goodone = nullptr;
  hraw1goodone_time = nullptr;
  hraw1badone_time = nullptr;
  hraw1goodoneNoise_time = nullptr;
  hraw1badoneNoise_time = nullptr;
  hRealTrackNoise_time = nullptr;

  hraw1PFrfr1GOOD = nullptr;
  hraw1PFrfr1BAD = nullptr;
  hraw1PFrfr1GOOD_time = nullptr;
  hraw1PFrfr1BAD_time = nullptr;
  hraw1PFrfr2GOOD = nullptr;
  hraw1PFrfr2BAD = nullptr;
  hraw1PFrfr2GOOD_time = nullptr;
  hraw1PFrfr2BAD_time = nullptr;
  hraw1NoiseGOOD_time = nullptr;
  hraw1NoiseBAD_time = nullptr;

  hraw1Noise = nullptr;
  hraw1Pedestal = nullptr;
  hraw1CDS = nullptr;
  hraw1Signal = nullptr;
  hraw1Noise_time = nullptr;
  hraw1Pedestal_time = nullptr;
  hraw1CDS_time = nullptr;
  hraw1Signal_time = nullptr;

  //------------------------------------------------------------------------------
  //--- efficiency
  // MG 2010/06/04
  //------------------------------------------------------------------------------

  effimap = nullptr;
  effi_vs_TrkHitDist = nullptr;
  effiCorr_vs_TrkHitDist = nullptr;
  goodtracks = nullptr;
  TrkInMimo = nullptr;
  effinpixel = nullptr;
  HitInPixel = nullptr;
  TrkInPixel = nullptr;

  //------------------------------------------------------------------------------
  //--- hit position
  //------------------------------------------------------------------------------

  cHitMap = nullptr;
  cHitMap2 = nullptr;
  cPixEvent = nullptr;
  cPosStudy1 = nullptr;
  cPosStudy2 = nullptr;
  cPosStudy_TrackhitPos_vs_Mult = nullptr;
  cPosStudy_TrackhitPos_vs_Mult2 = nullptr;

  //------------------------------------------------------------------------------
  //--- ?
  //------------------------------------------------------------------------------

  hseedQvsS2NAll = nullptr;
  hseedQvsS2NGood = nullptr;
  hQseedAll = nullptr;
  hIndex2x2 = nullptr;
  hnhit = nullptr;
  vec = nullptr;
  tudv = nullptr;
  tvdu = nullptr;
  hudv = nullptr;
  hvdu = nullptr;
  tuhu1 = nullptr;
  tvhv1 = nullptr;
  tuhu = nullptr;
  tvhv = nullptr;
  tuhv = nullptr;
  tvhu = nullptr;

  hAllHvvsAllTv = nullptr;
  hAllHuvsAllTv = nullptr;
  hAllHvvsAllTu = nullptr;
  hAllHuvsAllTu = nullptr;

  huCG5tu = nullptr;
  hvCG5tv = nullptr;
  huCGtu = nullptr;
  hvCGtv = nullptr;
  huCGtuInPix = nullptr;
  hvCGtvInPix = nullptr;
  huCGtuInPix5 = nullptr;
  hvCGtvInPix5 = nullptr;
  htuvInPix = nullptr;
  ProfUCG = nullptr;
  ProfVCG = nullptr;
  ProfACGn = nullptr;
  hUeta3TuInPix = nullptr;
  hVeta3TvInPix = nullptr;
  hUcorTuInPix = nullptr;
  hVcorTvInPix = nullptr;
  huCGtuInPix4 = nullptr;
  huCG2x2tuInPix = nullptr;
  hvCG2x2tvInPix = nullptr;
  hEta2x2tu1L = nullptr;
  hEta2x2tv1L = nullptr;
  hEta2x2tu2L = nullptr;
  hEta2x2tv2L = nullptr;
  hCG2x2tu1L = nullptr;
  hCG2x2tv1L = nullptr;

  hEta2x2 = nullptr;
  hEta2x2U = nullptr;
  hEta2x2V = nullptr;
  hEta2x2m = nullptr;

  hEta2x2UL = nullptr;
  hEta2x2VL = nullptr;

  hdCGEtaU = nullptr;
  hdCGEtaV = nullptr;
  hdCGDigU = nullptr;
  hdCGDigV = nullptr;
  hdCGDigUV = nullptr;
  htmp5 = nullptr;
  hEta2x2tu2 = nullptr;
  hEta2x2tv2 = nullptr;
  hTHCorr2 = nullptr;
  hDifCorU = nullptr;
  hDifCorV = nullptr;
  hAllHuvsAllTu2 = nullptr;
  hAllHvvsAllTv2 = nullptr;
  hAlignHuTu = nullptr;
  hAlignHvTv = nullptr;
  hAlignHuTuvsTv = nullptr;
  hAlignHvTvvsTu = nullptr;
  hAlignHuTuvsTu = nullptr;
  hAlignHvTvvsTv = nullptr;
  hAlignHuTuvsHv = nullptr;
  hAlignHvTvvsHu = nullptr;
  hAllHvvsAllTv3 = nullptr;
  h2dallhits = nullptr;
  h2dgoodhits = nullptr;
  h2dmatchedhits = nullptr;
  h2DpictureMatched = nullptr;
  hEta2x2vsInd = nullptr;
  hChargeVsPosition = nullptr;
  hChargeVsDistance = nullptr;
  hNorm = nullptr;
  hAllHitsInPixel = nullptr;
  h1RmsOnTheta = nullptr;
  h2RmsOnThetaScanVsThetaCut = nullptr;
  h1ProjectionOnX = nullptr;
  h1ProjectionOnY = nullptr;
  h1NumberOfHitsMult = nullptr;
  for (int i = 0; i < 4; i++)
  {
    h1ProjectionOnXMult[i] = nullptr;
    h1ProjectionOnYMult[i] = nullptr;
    h2GoodHitsMult[i] = nullptr;
  }
  h1Sigma = nullptr;
  for (int i = 0; i < 6; i++)
  {
    hdCGDigUVMult[i] = nullptr;
  }
  h1SigmaSlices = nullptr;
  h1SigmaSlices_Pull = nullptr;
  h2DXprimeVsYprime = nullptr;

  huv = nullptr;
  huv_rate = nullptr;       // AP 2014/09/25
  h_pixels_event = nullptr; // AP 2014/10/07
  hxy = nullptr;
  huvBad = nullptr;
  hxyBad = nullptr;
  tuv = nullptr;
  tuv1 = nullptr;

  hhu = nullptr;
  hhv = nullptr;
  hhx = nullptr;
  hhy = nullptr;
  hhuS = nullptr;
  hhvS = nullptr;
  htu = nullptr;
  htv = nullptr;
  //---ab
  htuhtv = nullptr;
  hGOODqcel = nullptr;

  hgoodSeedPixel = nullptr;
  h2DgoodSeedPixel = nullptr;
  hSeedBetweenDist = nullptr;
  hEta3DigU = nullptr;
  hEta3DigV = nullptr;
  hqc_nc = nullptr;

  hnpix = nullptr;
  hnpix_nc = nullptr;

  FalseHitMap = nullptr;
  hClusterChargeProfile = nullptr;
  hClusterChargeNorm = nullptr;
  for (int i = 0; i < 2; i++)
  {
    etal[i] = nullptr;
    etal1[i] = nullptr;
  }
  for (int i = 0; i < 2; i++)
  {
    vecu[i] = nullptr;
    vecv[i] = nullptr;
  }
  for (int i = 0; i < 64; i++)
    hChargeInSeedInGroup[i] = nullptr;
  DuvCG = nullptr;
  duvall = nullptr;

  hCDSvar = nullptr;
  CDSVarvsTime = nullptr;
  dtime = nullptr;

  //------------------------------------------------------------------------------
  //----- Spatial Resolution
  //------------------------------------------------------------------------------
  hAllHuvsAllTu1 = nullptr;
  hAllHvvsAllTv1 = nullptr;
  huCGtu1 = nullptr;
  hvCGtv1 = nullptr;
  huCGtu1_vs_Mult = nullptr;
  hvCGtv1_vs_Mult = nullptr;
  huvCGtuv = nullptr;
  huvCGtuv_NoAsso = nullptr;
  hDist_Trck_Diode_Asso = nullptr;
  hDist_Trck_Diode_NoAsso = nullptr;
  hEffic_vs_Dist_Trck_Diode = nullptr;
  hDist_Trck_Diode_Asso_vs_Mult = nullptr;
  hnpixc_vs_TrkDistToDiode = nullptr;

  huCGwidth_vs_Mult = nullptr;
  hvCGwidth_vs_Mult = nullptr;
  huCGmean_vs_Mult = nullptr;
  hvCGmean_vs_Mult = nullptr;
  huCGtu1_vs_TracksPerEvent = nullptr;
  hvCGtv1_vs_TracksPerEvent = nullptr;
  grvCGVsMult = nullptr;
  huCGtu2 = nullptr;
  hvCGtv2 = nullptr;
  hCG2x2tu1 = nullptr;
  hCG2x2tv1 = nullptr;
  hCG5URes = nullptr;
  hCG5VRes = nullptr;
  hTuHuCorr = nullptr;
  hTvHvCorr = nullptr;
  hEta2x2tu1 = nullptr;
  hEta2x2tv1 = nullptr;
  hEtaURes = nullptr;
  hEtaVRes = nullptr;
  hEta3URes = nullptr;
  hEta3VRes = nullptr;
  hEtaU_2x2Res = nullptr;
  hEtaV_2x2Res = nullptr;
  hEtaU_5x5Res = nullptr;
  hEtaV_5x5Res = nullptr;
  hAHTURes = nullptr;
  hAHTVRes = nullptr;

  cres = nullptr;
  cresfit_CGWidth_vs_Mult = nullptr;
  cresfit_CGMean_vs_Mult = nullptr;
  leg_CGWidth_vs_Mult = nullptr;
  cresfit_CGWidth_vs_TracksPerEvent = nullptr;

  //------------------------------------------------------------------------------
  // Fake rate study
  //------------------------------------------------------------------------------
  hNhitperpixel = nullptr;
  hNhitRateperpixel = nullptr;
  hPixelsPerFakeRate = nullptr;

  cfake = nullptr;

  //------------------------------------------------------------------------------
  // tracks properties and chi2
  //------------------------------------------------------------------------------

  hNTracksPerEvent = nullptr;
  hNTracksPerEventievt = nullptr;
  hNGoodGeomTracksPerEvent = nullptr;

  hAllTvTu = nullptr;
  hAllTu = nullptr;
  hAllTv = nullptr;
  hGoodChi2TvTu = nullptr;
  hGoodChi2Tu = nullptr;
  hGoodChi2Tv = nullptr;
  hGoodChi2Tx = nullptr;
  hGoodChi2Ty = nullptr;
  hGoodChi2AngleXZ = nullptr;
  hGoodChi2AngleYZ = nullptr;

  hchi2_c = nullptr;
  hchi2_nc = nullptr;
  hchi2 = nullptr;

  hTrackToClusterMinDistance = nullptr;
  hTrackTo2ndclosestClusterDistance = nullptr;
  hWrongAssociationProba = nullptr;
  hMinDistance_vs_2ndDistance = nullptr;

  hNhitsPerTrack_all = nullptr;
  hNhitsPerTrack_good = nullptr;
  hNhitsPerTrack_matched = nullptr;

  static const int NRefPlane = 10;
  for (int i = 0; i < NRefPlane; i++)
  {
    hRef_Tud_vs_Tv[i] = nullptr;
    hRef_Tud_vs_Tu[i] = nullptr;
    hRef_Tud_vs_TDOX[i] = nullptr;
    hRef_Tud_vs_TDOY[i] = nullptr;
    hRef_Tud_vs_Chi2[i] = nullptr;
    hRef_Tud_vs_Tu1[i] = nullptr;
    hRef_Tud_vs_Tk1[i] = nullptr;
  }

  // Char_t titreref[100] ;
  // Char_t nomRef[50];

  //------------------------------------------------------------------------------
  // Good/All hits.
  //------------------------------------------------------------------------------
  hnhitievt = nullptr;
  hnahitievt = nullptr;
  hnGOODhit = nullptr;
  hnGOODhitwhentrack = nullptr;
  hnTracksinSensorPerEvent = nullptr;
  hTrackDensity = nullptr;
  hnTracksInGeomatrixVsTrackPerEvent = nullptr;
  hnGOODhitInGeomatrixVsTrackPerEvent = nullptr;
  hEfficiencyInGeomatrixVsTrackPerEvent = nullptr;
  hResidueUInGeomatrixVsTrackPerEvent = nullptr;
  hResidueVInGeomatrixVsTrackPerEvent = nullptr;
  leg_TracksinSensorPerEvent = nullptr;
  l_average_Effic = nullptr;
  l_average_Effic_PErr = nullptr;
  l_average_Effic_MErr = nullptr;

  hAllS2N = nullptr;
  hallhitSN = nullptr;
  hallSNneighbour = nullptr;
  hgoodSNneighbour = nullptr;
  hgoodSN_vs_SNN = nullptr;
  hallSN_vs_SNN = nullptr;

  hallSN_minus_hgoodSN_vs_SNN = nullptr;
  hdistchi2 = nullptr;

  //------------------------------------------------------------------------------
  // REAL S/N NEIGHBOURS ----ab25072005
  //------------------------------------------------------------------------------

  hSN_vs_SNNReal = nullptr;

  //------------------------------------------------------------------------------
  // Charge by ordered pixels Histos.
  //------------------------------------------------------------------------------
  //-------------
  hGOODChargeOrder1 = nullptr;
  hGOODChargeOrder2 = nullptr;
  hGOODChargeOrder3 = nullptr;
  hGOODChargeOrder4 = nullptr;
  hGOODChargeOrder5 = nullptr;
  hGOODChargeOrder6 = nullptr;
  hGOODChargeOrder7 = nullptr;
  hGOODChargeOrder8 = nullptr;
  hGOODChargeOrder9 = nullptr;
  hGOODChargeCor_1_2 = nullptr;
  hGOODChargeCor_1_3 = nullptr;
  hGOODChargeCor_1_4 = nullptr;
  hGOODChargeCor_2_3 = nullptr;
  hGOODChargeSum_4 = nullptr;
  hGOODChargeRap_1_over_2 = nullptr;
  hGOODChargeRap_1_over_3 = nullptr;
  hGOODChargeRap_1_over_4 = nullptr;
  hGOODChargeRap_2_over_3 = nullptr;

  //------------------------------------------------------------------------------
  //----- ?
  //------------------------------------------------------------------------------

  // histograms needed for charge distribution function:
  hS2N2All = nullptr;
  hS2N2nd = nullptr;
  hS2N2RH = nullptr;
  hPedestal = nullptr;

  //----ALL HITS: Signal / noise seed vs Signal / noise in the pixel range [i-j]
  for (int i = 0; i < jpixmax; i++)
  {
    hsnALL_seed_vs_pix_0[i] = nullptr;
    hsnALL_seed_vs_pix_1[i] = nullptr;
  }
  hChargeIntegral1 = nullptr;
  hChargeNorm1 = nullptr;
  hChargeIntegral2 = nullptr;
  hChargeNorm2 = nullptr;
  hChargeIntegral3 = nullptr;
  hChargeNorm3 = nullptr;
  hChargeIntegral4 = nullptr;
  hChargeNorm4 = nullptr;
  hChargeIntegral = nullptr;
  hChargeIntSmoothed = nullptr;

  //------------------------------------------------------------------------------
  // Check pixel homogeneity depending on real track position
  //------------------------------------------------------------------------------
  hHOM_tu_tv_modulo = nullptr;
  hHOM_ResU_tu = nullptr;
  hHOM_ResV_tv = nullptr;
  ProfHOM_ResU_tu = nullptr;
  ProfHOM_ResV_tv = nullptr;

  hHOM_modUCG_modtu = nullptr;
  hHOM_modVCG_modtv = nullptr;
  hHOM_modUeta3_modtu = nullptr;
  hHOM_modVeta3_modtv = nullptr;
  hHOM_modUeta3_realtu = nullptr;
  hHOM_modVeta3_realtv = nullptr;
  hHOM_modUCG_realtu = nullptr;
  hHOM_modVCG_realtv = nullptr;
  hHOM_modUeta3_Eta3U = nullptr;
  hHOM_modVeta3_Eta3V = nullptr;
  hHOM_modUeta3_modVeta3 = nullptr;
  hHOM_modUCG_modVCG = nullptr;
  hHOM_modUeta3_modUCG = nullptr;
  hHOM_modVeta3_modVCG = nullptr;

  hHOM_Charge_diodedist_alg = nullptr;
  hHOM_Charge_diodedist_alg_u = nullptr; // clm
  hHOM_Charge_diodedist_alg_v = nullptr; // clm
  hHOM_Charge_diodedist = nullptr;
  hHOM_Charge2_diodedist = nullptr;
  hHOM_Charge4_diodedist = nullptr;
  hHOM_Charge9_diodedist = nullptr;
  hHOM_Charge25_diodedist = nullptr;
  hHOM_Noise_diodedist = nullptr;

  // start CLM 2013/01/23
  ProfhGOODCharge_Charge_DiodePosition = nullptr;
  ProfhGOODCharge_Charge_DiodePositionSeedQLT300 = nullptr;
  ProfhGOODCharge_Charge_DiodePositionSeedQGT2000 = nullptr;
  // clm
  ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow = nullptr;
  // clm
  ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow_seed = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow_seed = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow_seed = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow_seed = nullptr;
  // clm
  ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow_1stcrown = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow_1stcrown = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow_1stcrown = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow_1stcrown = nullptr;
  // clm
  ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow_2ndcrown = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow_2ndcrown = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow_2ndcrown = nullptr;
  ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow_2ndcrown = nullptr;

  hDistVSeedOtherOldCalc = nullptr;
  hDistVSeedOtherNewCalc = nullptr;

  h2dCharge_Charge_DiodePosition_Track = nullptr;
  h2dCharge_Charge_DiodePosition_CluSize = nullptr;

  hNpixInClu = nullptr;
  hQpixInClu = nullptr;

  ProfhGOODCharge_Charge_DiodePositionSimpDist = nullptr;

  hHOM_Charge_diodedist3D = nullptr;
  hHOM_Charge2_diodedist3D = nullptr;
  hHOM_Charge4_diodedist3D = nullptr;
  hHOM_Charge9_diodedist3D = nullptr;
  hHOM_Charge25_diodedist3D = nullptr;
  // end CLM 2013/01/23

  // TProfile *ProfHOM_Charge_diodedist_alg = nullptr;  // JB 2010/03/11
  // TProfile *ProfHOM_Charge_diodedist_alg_v = nullptr;  //clm 2013/01/23
  // TProfile *ProfHOM_Charge_diodedist_alg_u = nullptr;  //clm 2013/01/23
  // TProfile *ProfHOM_Charge_diodedist = nullptr;
  // TProfile *ProfHOM_Charge2_diodedist = nullptr;
  // TProfile *ProfHOM_Charge4_diodedist = nullptr;
  // TProfile *ProfHOM_Charge9_diodedist = nullptr;
  // TProfile *ProfHOM_Charge25_diodedist = nullptr;

  hHOM_SNseed_diodedist = nullptr;
  // TProfile *ProfHOM_SNseed_diodedist = nullptr;

  hHOM_Charge_diodedist_00_10 = nullptr;
  hHOM_Charge_diodedist_10_20 = nullptr;
  hHOM_Charge_diodedist_20_30 = nullptr;
  hHOM_Charge_diodedist_30_40 = nullptr;
  hHOM_Charge_diodedist_40_50 = nullptr;
  hHOM_Charge_diodedist_50_60 = nullptr;
  hHOM_Charge_diodedist_60_70 = nullptr;
  hHOM_Charge_diodedist_70_80 = nullptr;
  hHOM_Charge_diodedist_80_90 = nullptr;
  hHOM_Charge_diodedist_90_inf = nullptr;

  hHOM_DU_Nevent = nullptr;
  hHOM_DV_Nevent = nullptr;
  hHOM_modtu_Nevent = nullptr;
  hHOM_modtv_Nevent = nullptr;

  hHOM_modUCG_Nevent = nullptr;
  hHOM_modVCG_Nevent = nullptr;

  hHOM_ResEta25_U = nullptr;
  hHOM_ResEta25_V = nullptr;

  hHOM_modUeta2x2_modVeta2x2 = nullptr;
  hHOM_modUCG2x2_modVCG2x2 = nullptr;
  hHOM_modUeta5x5_modVeta5x5 = nullptr;
  hHOM_modUCG5x5_modVCG5x5 = nullptr;

  //--------------------------------------------------------------------------
  // cluster shape study, JB 2010/04/13 based on cdritsa stuff
  //--------------------------------------------------------------------------

  hClusterMeanForm = nullptr;
  hprojLForm = nullptr;
  hprojCForm = nullptr;
  for (int i = 0; i < 10; i++)
  {
    hCountPixels[i] = nullptr;
    hPixelsOverSNR[i] = nullptr;
    hPixelsOverCharge[i] = nullptr;
    hprojL[i] = nullptr;
    hprojC[i] = nullptr;
    Cluster[i] = nullptr;
    hMultVsFormFactor[i] = nullptr;
    hMultVsFormFactor1D[i] = nullptr; // clm added back 1D clu mult. 2013.08.25
  }
  hChargeCoG_Correl = nullptr;
  hChargeCoG_Correl2 = nullptr;
  h_SNRratioL = nullptr;
  hClusterTest1 = nullptr;
  hClusterTest2 = nullptr;
  hClusterSizeInLines = nullptr; // JB 2014/03/31
  hClusterSizeInColumns = nullptr;
  hClusterSize2DInLineAndColumns = nullptr; // AP 2014/09/25
  hChargeDistrInSeed = nullptr;
  hChargeDistrInLine = nullptr;
  hChargeDistrIn3rdLeftNeigh = nullptr;
  hChargeDistrIn2ndLeftNeigh = nullptr;
  hChargeDistrIn1stLeftNeigh = nullptr;
  hChargeDistrIn3rdRightNeigh = nullptr;
  hChargeDistrIn2ndRightNeigh = nullptr;
  hChargeDistrIn1stRightNeigh = nullptr;
  hClusterTypes = nullptr;        // JB 2014/03/31
  hClusterTypesBeyond4 = nullptr; // JB 2015/05/12

  ang1 = nullptr;
  ang2 = nullptr;
  ang3 = nullptr;
  ang4 = nullptr;
  ang5 = nullptr;
  ang5v2 = nullptr;
  ang6 = nullptr; // clm added back 1D clu mult. 2013.08.25
  shape = nullptr;
  solo = nullptr;
  projL = nullptr;
  projC = nullptr;

  //------------------------------------------------------------------------------
  // Binary output histos.
  //------------------------------------------------------------------------------
  hBinary_test = nullptr;
  hBinary_NumberOf_1_ALL = nullptr;
  hBinary_NumberOf_1_goodhit = nullptr;
  hBinary_NumberOf_1_submatrix = nullptr;
  hBinary_Nhitperpixel_submatrix = nullptr;
  hBinary_NhitRateperpixel_submatrix = nullptr;

  //------------------------------------------------------------------------------
  // MiniVectors Histograms  = nullptr; NCS , 2010/01/21
  //------------------------------------------------------------------------------
  hDiffPosX = nullptr;     // NCS 210110
  hDiffPosY = nullptr;     // NCS 210110
  hDiffAngleX = nullptr;   // NCS 210110
  hDiffAngleY = nullptr;   // NCS 210110
  hDiffAngleX11 = nullptr; // JB 2011/11/01
  hDiffAngleY11 = nullptr;
  hDiffAngleX12 = nullptr; // JB 2011/11/01
  hDiffAngleY12 = nullptr;
  hDiffAngleX21 = nullptr; // JB 2011/11/01
  hDiffAngleY21 = nullptr;
  hDiffAngleX22 = nullptr; // JB 2011/11/01
  hDiffAngleY22 = nullptr;
  hDiffAngleXg1g1 = nullptr; // JB 2011/11/01
  hDiffAngleYg1g1 = nullptr;
  hxtxPL3 = nullptr; // NCS 260110 not CG position only Hu Hv
  hytyPL3 = nullptr; // NCS 260110
  hxtxPL4 = nullptr; // NCS 260110
  hytyPL4 = nullptr; // NCS 260110
  hutuPL3 = nullptr; // NCS 260110 not CG position only Hu Hv
  hvtvPL3 = nullptr; // NCS 260110
  hutuPL4 = nullptr; // NCS 260110
  hvtvPL4 = nullptr; // NCS 260110
  // TH2F* hdiffydiffx = nullptr; // JB 2011/11/01
  // TH2F* hdiffydiffx11 = nullptr; // JB 2011/11/01
  // TH2F* hdiffydiffx12 = nullptr; // JB 2011/11/01
  // TH2F* hdiffydiffx21 = nullptr; // JB 2011/11/01
  // TH2F* hdiffydiffx22 = nullptr; // JB 2011/11/01
  // TH2F* hdiffydiffxg1g1 = nullptr; // JB 2011/11/01

  //------------------------------------------------------------------------------
  // MimosaVertexFinder Histograms  = nullptr; LC , 2012/09/06
  //------------------------------------------------------------------------------

  hVertexPosX = nullptr;
  hVertexPosY = nullptr;
  hVertexPosZ = nullptr;
  hVertexTrackDistance = nullptr;
  hVertexPosXY = nullptr;
  hVertexTrackChi2 = nullptr;

  cVertexFinder = nullptr;

  cMiniVec = nullptr;
  cMiniVec1 = nullptr;

  //------------------------------------------------------------------------------
  // canvas
  //------------------------------------------------------------------------------
  c2 = nullptr;
  casym = nullptr;
  ceffi = nullptr;
  ccomp = nullptr;
  ccomp2 = nullptr;
  ccomp3 = nullptr;
  c3 = nullptr;
  c3_2 = nullptr;
  c4 = nullptr;
  cClusterProperties2 = nullptr;
  cClusterProperties3 = nullptr;
  cClusterProperties4 = nullptr;
  c4_2 = nullptr;
  c5 = nullptr;
  c6 = nullptr;
  c7 = nullptr;

  for (int i = 0; i < numcanvasSN; i++)
    cSN[i] = nullptr; //!
  for (int i = 0; i < numcanvas; i++)
    cRef[i] = nullptr; //!
  for (int i = 0; i < numcanvasOptimize; i++)
    cOptimize[i] = nullptr; //!

  cM8 = nullptr;

  PixHom = nullptr;
  PixHom2 = nullptr;
  PixHom3 = nullptr;
  PixHom4 = nullptr;
  PixHom5 = nullptr;
  PixHom6 = nullptr;

  // graphs
  greff = nullptr;
  grnum = nullptr;
  grevt = nullptr;
  ChargeSpread = nullptr;

  bar2 = nullptr;

  //------------------------------------------------------------------------------
  // MC Geneation for telescope resolution evaluation
  //------------------------------------------------------------------------------
  hTrackChi2_MC = nullptr;
  hTrackSlopeX_MC = nullptr;
  hTrackSlopeY_MC = nullptr;
  hTrackNHits_MC = nullptr;
  hTrackPlanesUsed_MC = nullptr;
  for (int i = 0; i < MynPlanes; i++)
  {
    hResidualU_MC[i] = nullptr;
    hResidualV_MC[i] = nullptr;
  }
  for (int i = 0; i < MyBins_GlobalResolution; i++)
  {
    hTrackResidualAtDUT_U[i] = nullptr;
    hTrackResidualAtDUT_V[i] = nullptr;
  }

  h_TelescopeResolution_Scan_U = nullptr;
  h_TelescopeResolution_Scan_V = nullptr;
  h_ResidualAtDUT_Scan_U = nullptr;
  h_ResidualAtDUT_Scan_V = nullptr;

  lSolResolutionU1 = nullptr;
  lSolResolutionU2 = nullptr;
  lSolResolutionU3 = nullptr;
  lSolResolutionV1 = nullptr;
  lSolResolutionV2 = nullptr;
  lSolResolutionV3 = nullptr;

  cgeometry = nullptr;
  ccumultrack3 = nullptr;
  ccumultrack4 = nullptr;
  ccumultrack5 = nullptr;
  for (int i = 0; i < MyBins_GlobalResolution; i++)
    cresAtDUT[i] = nullptr;
  cresAtDUT_vs_SensorSP = nullptr;

  for (int i = 0; i < MyMaxEventDisplay; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      RectracksXZ[i][j] = nullptr;
      RectracksYZ[i][j] = nullptr;
    }
  }
  leg_TelResMC = nullptr;

  //------------------------------------------------------------------------------
  // User histograms
  //------------------------------------------------------------------------------
  cUser = nullptr;
  hUserHitCorrelationLine = nullptr;
  hUserHitCorrelationCol = nullptr;

  NPages = 0;
}

//______________________________________________________________________________
//
void MimosaAnalysis::MergingPDFFiles(void)
{
  std::stringstream ss;
  ss << "run" << RunNumber << "Pl" << ThePlaneNumber << "_ClCharge";
  std::string fOutName = fTool.LocalizeDirName(ss.str().c_str());
  const auto EPSName_final = CreateGlobalResultDir() + fOutName + ".pdf";

  std::string listOfFiles{""};

  for (int i = 0; i < NPages; i++)
  {
    std::stringstream tmpSS;
    tmpSS << CreateGlobalResultDir() << fOutName << "_tmpFile" << (i + 1) << ".pdf" << "   ";
    listOfFiles += tmpSS.str();
  }
  auto command = "gs -dQUIET -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=" +
                 EPSName_final + " -dBATCH ";
  command += listOfFiles;
  gSystem->Exec(command.c_str());
  command = "rm ";
  command += listOfFiles;
  gSystem->Exec(command.c_str());
}
//______________________________________________________________________________
//
int MimosaAnalysis::OpenInputFile()
{
  // Open the input file containing the TTree:
  // returns the number of entries in the tree
  // or 0 if no file was loaded.
  //
  // JB 2010/08/30
  // Modified: JB 2011/07/07 to localize path names
  // Modified: JB 2012/11/21 check the file number exists
  // Modified: JB 2013/09/19 allow for a pre-definition of the file name
  // Modified: JB 2014/02/10 Error on # events <=0

  const auto fileNumber = GetFileNumber();

  // If RootFile not yet define, use a number to set it
  if (m_rootFile.empty() && fileNumber)
  {
    std::stringstream ss;
    ss << fSession->GetSummaryFilePath() << "/run"
       << fSession->GetRunNumber() << "_0"
       << fileNumber << ".root";
    m_rootFile = fTool.LocalizeDirName(ss.str().c_str()); // JB 2011/07/07
  }

  if (m_rootFile.empty())
  {
    Error("MAnalysis", "No input file specified ! --> STOP.");
    return 0;
  }

  // Open the file if existing

  auto fileInputTree = std::make_unique<TFile>(m_rootFile.c_str());
  if (!fileInputTree || !fileInputTree->IsOpen())
  {
    Error("MimosaAnalysis", "Failed to open input file %s", m_rootFile.c_str());
    return 0;
  }

  Info("MimosaAnalysis", "Input file %s opened", m_rootFile.c_str());

  t = static_cast<TTree *>(fileInputTree->Get("T"));
  if (!t)
  {
    Error("MimosaAnalysis", "TTree 'T' not found in file %s", m_rootFile.c_str());
    return 0;
  }

  Evt = std::make_unique<DEvent>();
  branch = t->GetBranch("fEvent");
  if (!branch)
  {
    Error("MimosaAnalysis", "Branch 'fEvent' not found in the TTree.");
    return 0;
  }
  branch->SetAddress(&Evt);

  m_nEvt = static_cast<int>(t->GetEntries());
  if (m_nEvt <= 0)
  {
    Error("MimosaPro",
          " The input file contains an incorrect number of events %d!", m_nEvt);
    return 0;
  }
  else
  {
    Info("MimosaPro", "There is %d events in the input file.", m_nEvt);
  }

  return m_nEvt;
}

std::string MimosaAnalysis::CreateGlobalResultDir()
{
  // Create Results/#### directory to store results
  //
  // Modified: JB 2011/07/07 to localize path names

  if (!CheckIfDone("init"))
    return ""; // correction from 0;, JB 2009/07/17

  if (MimosaType == 0)
  {
    Error("CreateGlobalResultDir",
          "MimosaType not set! Please run MimosaPro first");
    return "";
  }

  const std::string resultDirName =
      fWorkingDirectory + "/results_ana_M" + std::to_string(MimosaType) + "/";
  fTool.LocalizeDirName(resultDirName.data()); // JB 2011/07/07

  if (MimoDebug)
    std::cout << "Global Result Dir : " << resultDirName << std::endl;

  std::filesystem::create_directories(resultDirName);

  return resultDirName;
}

//______________________________________________________________________________
//
void MimosaAnalysis::PrepareOnlineDisplay()
{
  // Create the canvas for plots displayed online (during analysis)
  //
  // JB 2010/09/30

  gStyle->SetPalette(1);
  MainCanvas = new TCanvas("MainCanvas", "Main display", 550, 10, 750, 550);
  MainCanvas->Draw();
  // TCanvas *MissedHit = new TCanvas("MissedHit","Missed Hits: red-cluster,
  // green-track",200,10,300,350); MissedHit->Draw();
  MainCanvas->cd();

  Char_t GlobalTitle[50];
  sprintf(GlobalTitle, "Run %d Plane %d Gain: %f", RunNumber, ThePlaneNumber,
          calibration);
  TPaveLabel *gt = new TPaveLabel(0.01, 0.94, 0.2, 0.9999, GlobalTitle);
  gt->Draw();

  Info("MimosaAnalysis", "Display prepared");
}

//______________________________________________________________________________
//
void MimosaAnalysis::SetPlaneNumber(int aPlaneNumber)
{

  // JB 2012/10/30

  if (0 < aPlaneNumber &&
      aPlaneNumber <= fSession->GetTracker()->GetPlanesN())
  {
    ThePlaneNumber = aPlaneNumber;
  }
  else
  {
    Warning("SetPlaneNumber",
            "The required plane %d does not exist! Nothing done.",
            aPlaneNumber);
  }
}

//______________________________________________________________________________
//
void MimosaAnalysis::SetLadderNumber(int aNumber)
{

  // JB 2014/01/10

  if (0 < aNumber && aNumber <= fSession->GetTracker()->GetNumberOfLadders())
  {
    TheLadderNumber = aNumber;
  }
  else
  {
    Warning("SetLadderNumber",
            "The required ladder %d does not exist! Nothing done.", aNumber);
  }
}

//______________________________________________________________________________
//
void MimosaAnalysis::GetAlignment()
{
  // To get alignment parameters from the CorPar file
  // JB 2011/06/19
  // Modified JB 2011/06/30 swap two rotation angles and use convolution method
  // for alignment pars Modified JB 2011/10/31 Print corrected into PrintAll
  // Modified JB 2012/10/30 Stronger condition on CorStatus
  // Modified JB 2014/02/10 Deal with Ladder case

  bool gotAlignment = kFALSE;

  if (theCorFile->IsOpen())
  {
    if (theCorFile->FindKey("alignement") &&
        CorStatus != 2)
    { // if the CorPar files is already initialized with
      // proper objects
      Info("GetAlignment", "Getting alignment from CorPar file.");
      alignement = (DPrecAlign *)theCorFile->Get("alignement");
      cout << "===========================================" << endl;
      cout << "The parameters are taken from CorPar" << RunNumber << "_"
           << ThePlaneNumber << ".root" << endl;
      cout << "===========================================" << endl;
      gotAlignment = kTRUE;
    }
  }

  if (!gotAlignment)
  { // take the config file information for alignment
    // Corrected for angles swap and convolution with alignmentU,V,Tilt, JB
    // 2011/06/30
    Info("GetAlignment", "Getting alignment from configuration file.");

    int DPrecAlignMethod =
        fSession->GetSetup()->GetTrackerPar().DPrecAlignMethod; // LC 2015/01/31
    alignement = new DPrecAlign(DPrecAlignMethod);

    // Decide if analysis on Plane or Ladder
    if (ThePlaneNumber)
    {
      alignement->SetTranslation(
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Position(0),
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Position(1),
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Position(2));
      alignement->SetRotations(
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Tilt(0),
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Tilt(1),
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Tilt(2));
      alignement->ConvoluteUVWAlignment(
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).AlignmentU,
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).AlignmentV,
          fSession->GetSetup()->GetPlanePar(ThePlaneNumber).AlignmentTilt);
    }

    else if (TheLadderNumber)
    {
      alignement->SetTranslation(
          fSession->GetSetup()->GetLadderPar(TheLadderNumber).Position);
      alignement->SetRotations(
          fSession->GetSetup()->GetLadderPar(TheLadderNumber).Tilt);
    }

    cout << "===========================================" << endl;
    cout << "Initial values for the parameters are taken" << endl;
    cout << "===========================================" << endl;
  }

  alignement->PrintAll();

  Info("GetAlignment", "Alignment read");
}

//______________________________________________________________________________
//
void MimosaAnalysis::GetMiEta()
{
  // To get eta-function parametrisation for all the
  // cases: 2x2, 3x3, and profile
  //
  // Modified: JB 2013/10/04 Protection against no histo in CorPar file

  if (!theCorFile)
    InitCorPar(fSession->GetRunNumber(), fSession->GetPlaneNumber());
  Info("GetMiEta", "Status of InitCorPar = %d", CorStatus);

#if ROOT_VERSION_CODE < ROOT_VERSION(4, 03, 4)
  if (MimoDebug)
    cout << ROOT_VERSION_CODE << " : bad version" << endl;
  MimoDebug <= 1 ? gErrorIgnoreLevel = kSysError : 0;
#endif

  // cout<<"------getmieta test 0 "<<CorStatus<<endl;

  TProfile *prACG = 0;
  TProfile *prUCG = 0;
  TProfile *prVCG = 0;
  TH1F *hE3Int = 0;
  TH1F *hE3UInt = 0;
  TH1F *hE3VInt = 0;
  TH1F *hE2UInt = 0;
  TH1F *hE2VInt = 0;

  if (CorStatus != 2)
  {
    prACG = (TProfile *)theCorFile->Get("ProfACG");
    prUCG = (TProfile *)theCorFile->Get("ProfUCG");
    prVCG = (TProfile *)theCorFile->Get("ProfVCG");
    hE3Int = (TH1F *)theCorFile->Get("hEta3Int");
    hE3UInt = (TH1F *)theCorFile->Get("hEta3UInt");
    hE3VInt = (TH1F *)theCorFile->Get("hEta3VInt");
    hE2UInt = (TH1F *)theCorFile->Get("hEta2x2IntU");
    hE2VInt = (TH1F *)theCorFile->Get("hEta2x2IntV");
  }

  // Test the histo are there, otherwise quit !
  if (!prACG || !prUCG || !prVCG || !hE3Int || !hE3UInt || !hE3VInt ||
      !hE2UInt || !hE2VInt)
  {
    Warning("GetMiEta",
            "WARNING: the expected histograms to compute the Eta correction "
            "were not found! --> No Eta correction available");
    return;
  }

  dir->cd();
  if (MimoDebug)
    cout << "GetMiEta Drawing" << endl;
  MainCanvas->Divide(2, 2);

  MainCanvas->cd(1);
  if (CorStatus != 2)
  {
    hE3UInt->Draw();
  }
  MainCanvas->Update();
  MainCanvas->cd(2);
  if (CorStatus != 2)
  {
    ;
    hE3VInt->Draw();
  }
  MainCanvas->Update();
  MainCanvas->cd(3);
  if (CorStatus != 2)
  {
    hE2UInt->Draw();
  }
  MainCanvas->Update();
  MainCanvas->cd(4);
  if (CorStatus != 2)
  {
    hE2VInt->Draw();
  }

  if (MimoDebug)
    cout << "GetMiEta Parametring via Fit" << endl;
  TF1 *fitfa;
  TF1 *fitfu;
  TF1 *fitfv;

  if (CorStatus != 2)
  {
    fitfa = prACG->GetFunction("pol6"); // Take parameters  from the file
    fitfu = prUCG->GetFunction("pol6");
    fitfv = prVCG->GetFunction("pol6");
  }
#if ROOT_VERSION_CODE < ROOT_VERSION(4, 03, 4)
  MimoDebug <= 1 ? gErrorIgnoreLevel = 0 : 0;
#endif

  // float FitParam[7];
  if (CorStatus != 2 && fitfa != nullptr && fitfu != nullptr &&
      fitfv != nullptr)
  {
    for (int ii = 0; ii < 7; ii++)
    {
      // FitParam[ii] =fitfa->GetParameter(ii);
      FitParamU[ii] = fitfu->GetParameter(ii);
      FitParamV[ii] = fitfv->GetParameter(ii);
    }
  }
  // Get eta2x2 integral:
  else
  {
    for (int ii = 0; ii < 7; ii++)
    {
      // FitParam[ii] = 0.;
      FitParamU[ii] = 0.;
      FitParamV[ii] = 0.;
    }
    // FitParam[1] = 1.;
    FitParamU[1] = 1.;
    FitParamV[1] = 1.;
  }

  if (CorStatus != 2)
  {
    NBins2x2 = hE2UInt->GetNbinsX();
  }
  else
  {
    NBins2x2 = 40;
  }
  Edges2x2 = new float[NBins2x2];
  float *Contents2x2 = new float[NBins2x2];
  Contents2x2U = new float[NBins2x2];
  Contents2x2V = new float[NBins2x2];

  if (CorStatus != 2)
  {
    for (int i = 0; i < NBins2x2; i++)
    {
      Edges2x2[i] = hE2UInt->GetBinLowEdge(i);
      Contents2x2[i] = hE2UInt->GetBinContent(i);
      Contents2x2U[i] = hE2UInt->GetBinContent(i);
      Contents2x2V[i] = hE2VInt->GetBinContent(i);
    }
  }
  cout << "=============================" << endl;
  // Get eta3x3 integral:
  if (CorStatus != 2)
  {
    NBins3x3 = hE3Int->GetNbinsX();
  }
  else
  {
    NBins3x3 = 80;
  }
  Edges3x3 = new float[NBins3x3];
  float *Contents3x3 = new float[NBins3x3];
  Contents3x3U = new float[NBins3x3];
  Contents3x3V = new float[NBins3x3];
  if (CorStatus != 2)
  {
    for (int i = 0; i < NBins3x3; i++)
    {
      Edges3x3[i] = hE3Int->GetBinLowEdge(i);
      Contents3x3[i] = hE3Int->GetBinContent(i);
      Contents3x3U[i] = hE3UInt->GetBinContent(i);
      Contents3x3V[i] = hE3VInt->GetBinContent(i);
    }
  }
  // float FitParEta3[PolDeg];
  if (CorStatus != 2 && hE3Int->Integral() > 10 && hE3UInt->Integral() > 10 &&
      hE3VInt->Integral() > 10)
  {
    MainCanvas->cd(1);
    if (hE3Int->Integral() > 10)
      hE3Int->Fit("pol6");
    // TF1 *fitEta3 = hE3Int->GetFunction("pol6");
    // for(int ii=0; ii<PolDeg; ii++){
    //   FitParEta3[ii]=fitEta3->GetParameter(ii);
    // }

    MainCanvas->cd(1);
    if (hE3UInt->Integral() > 10)
      hE3UInt->Fit("pol6");
    TF1 *fitEta3U = hE3UInt->GetFunction("pol6");
    for (int ii = 0; ii < PolDeg; ii++)
    {
      FitParEta3U[ii] = fitEta3U->GetParameter(ii);
    }
    MainCanvas->cd(2);
    if (hE3VInt->Integral() > 10)
      hE3VInt->Fit("pol6");
    TF1 *fitEta3V = hE3Int->GetFunction("pol6");
    for (int ii = 0; ii < PolDeg; ii++)
    {
      FitParEta3V[ii] = fitEta3V->GetParameter(ii);
    }
  }
  else
  {
    for (int ii = 0; ii < 7; ii++)
    {
      // FitParEta3[ii] = 0.;
      FitParEta3U[ii] = 0.;
      FitParEta3V[ii] = 0.;
    }
    // FitParEta3[1] = 1.;
    FitParEta3U[1] = 1.;
    FitParEta3V[1] = 1.;
  }

  //--------------Get Eta via TTree
  if (MimoDebug)
    cout << "GetMiEta Parametring via TTree" << endl;
  //  cout<<nListe_CoG<<" "<<ilist<<" "<<Liste_CoGU[ilist]<<"
  //  "<<Liste_CoGV[ilist]<<endl;
  if (CorStatus != 2)
  {
    TTree *Etatree;
    Etatree = (TTree *)theCorFile->Get("TreeEta");
    int nEtatreeentries = (int)Etatree->GetEntries();

    float n_EtaU;
    float n_EtaV;
    READnListe_CoG = nEtatreeentries;
    //   cout<<"nEtatreeentries "<<nEtatreeentries<<" "<<READnListe_CoG<<endl;
    Etatree->SetBranchAddress("n_EtaU", &n_EtaU);
    Etatree->SetBranchAddress("n_EtaV", &n_EtaV);

    //---etaab debut  GetMiEta() read tree
    float n_Eta2x2U;
    float n_Eta2x2V;
    float n_Eta5x5U;
    float n_Eta5x5V;
    //  cout<<" TEST GET BRANCH 1"<<endl;
    if (Etatree->GetBranch("n_Eta2x2U"))
    {
      // cout<<" TEST GET BRANCH 2"<<endl;
      Etatree->SetBranchAddress("n_Eta2x2U", &n_Eta2x2U);
      // cout<<" TEST GET BRANCH 3"<<endl;
    }
    if (Etatree->GetBranch("n_Eta2x2V"))
    {
      Etatree->SetBranchAddress("n_Eta2x2V", &n_Eta2x2V);
    }
    if (Etatree->GetBranch("n_Eta5x5U"))
    {
      Etatree->SetBranchAddress("n_Eta5x5U", &n_Eta5x5U);
    }
    if (Etatree->GetBranch("n_Eta5x5V"))
    {
      Etatree->SetBranchAddress("n_Eta5x5V", &n_Eta5x5V);
    }
    for (int i = 0; i < READnListe_CoG; i++)
    {
      Etatree->GetEntry(i);
      READListe_CoGU.Set(i + 1);
      READListe_CoGV.Set(i + 1);
      READListe_CoGU[i] = n_EtaU;
      READListe_CoGV[i] = n_EtaV;
      READListe_CoGU.AddAt(n_EtaU, i);
      READListe_CoGV.AddAt(n_EtaV, i);
    }
    READnListe_CoG_eta2x2 = nEtatreeentries;
    for (int i = 0; i < READnListe_CoG_eta2x2; i++)
    {
      Etatree->GetEntry(i);
      READListe_CoGU_eta2x2.Set(i + 1);
      READListe_CoGV_eta2x2.Set(i + 1);
      READListe_CoGU_eta2x2[i] = n_Eta2x2U;
      READListe_CoGV_eta2x2[i] = n_Eta2x2V;
      READListe_CoGU_eta2x2.AddAt(n_Eta2x2U, i);
      READListe_CoGV_eta2x2.AddAt(n_Eta2x2V, i);
    }
    READnListe_CoG_eta5x5 = nEtatreeentries;
    for (int i = 0; i < READnListe_CoG_eta5x5; i++)
    {
      Etatree->GetEntry(i);
      READListe_CoGU_eta5x5.Set(i + 1);
      READListe_CoGV_eta5x5.Set(i + 1);
      READListe_CoGU_eta5x5[i] = n_Eta5x5U;
      READListe_CoGV_eta5x5[i] = n_Eta5x5V;
      READListe_CoGU_eta5x5.AddAt(n_Eta5x5U, i);
      READListe_CoGV_eta5x5.AddAt(n_Eta5x5V, i);
    }
    //---etaab fin GetMiEta()
    /*   cout<<"READListe"<<endl;
   for(int i=0 ; i< READnListe_CoG ; i++){
   cout<<i<<" "<<READListe_CoGU[i]<<" "<<READListe_CoGV[i]<<endl;
     }*/
  }

  //---etaab debut nListe_CoG
  nListe_CoG = 0;
  nListe_CoG_eta2x2 = 0;
  nListe_CoG_eta5x5 = 0;
  //---etaab

  Info("GetMiEta", "End of eta preparation %d", NBins3x3);
  MainCanvas->Update();
}

//_____________________________________________________________
//
void MimosaAnalysis::CreateNewEta()
{
  // creates new distributions needed for precise measurement
  //
  // Modified: JB 2012/08/18 Check of entries before proceed
  // Modified: JB 2012/09/05 Display modified
  // Modified: JB 2013/08/21 Condition on histo Integral to allow fit
  // Modified: JB 2013/09/19 Keep a single key in the Root ouput file

  if (NofClMatchTrack < 50)
  {
    Info("CreateNewEta",
         "Not enough data (%d) to estimate eta parameters -> nothing done.",
         NofClMatchTrack);
    return;
  }
  Info("CreateNewEta", "New Eta parameters are to be computed with %d cluster.",
       NofClMatchTrack);

  const int PolN = 7; // Polynom degree
  // TF1 *fitf;
  // float FitParam[PolN];
  if (ProfACGn->Integral() > 10)
  {
    ProfACGn->Fit("pol6", "F");
    // fitf = ProfACGn->GetFunction("pol6");
    // for(int ii=0; ii<PolN; ii++){
    //   FitParam[ii]=fitf->GetParameter(ii);
    // }
    if (MimoDebug > 1)
      printf("CreateNewEta: ProfACGn fitted with %f integral\n",
             ProfACGn->Integral());
  }
  else
  {
    // for(int ii=0; ii<PolN; ii++){
    //   FitParam[ii]=0.;
    // }
  }
  TF1 *fitfu;
  if (ProfUCG->Integral() > 10)
  {
    ProfUCG->Fit("pol6", "F");
    fitfu = ProfUCG->GetFunction("pol6");
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParamU[ii] = fitfu->GetParameter(ii);
    }
    if (MimoDebug > 1)
      printf("CreateNewEta: ProfUCG fitted with %f integral\n",
             ProfUCG->Integral());
  }
  else
  {
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParamU[ii] = 0.;
    }
  }
  TF1 *fitfv;
  if (ProfVCG->Integral() > 10)
  {
    ProfVCG->Fit("pol6", "F");
    fitfv = ProfVCG->GetFunction("pol6");
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParamV[ii] = fitfv->GetParameter(ii);
    }
    if (MimoDebug > 1)
      printf("CreateNewEta: ProfVCG fitted with %f integral\n",
             ProfVCG->Integral());
  }
  else
  {
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParamV[ii] = 0.;
    }
  }

  //----------------------
  // 3x3 eta functions. OLD
  //----------------------
  if (MimoDebug > 1)
    printf("CreateNewEta: computing correction function for 3x3(old)\n");

  // Eta:
  float tmpEta(0.), tmpEtaU(0.), tmpEtaV(0.);
  hEta3DigU->Draw();
  TH1F *hEta3Int = new TH1F(
      "hEta3Int", "Eta for 3x3 cluster", hEta3DigU->GetXaxis()->GetNbins(),
      hEta3DigU->GetXaxis()->GetXmin(), hEta3DigU->GetXaxis()->GetXmax());
  TH1F *hEta3UInt = new TH1F(
      "hEta3UInt", "Eta for 3x3 cluster(U)", hEta3DigU->GetXaxis()->GetNbins(),
      hEta3DigU->GetXaxis()->GetXmin(), hEta3DigU->GetXaxis()->GetXmax());
  TH1F *hEta3VInt = new TH1F(
      "hEta3VInt", "Eta for 3x3 cluster(V)", hEta3DigV->GetXaxis()->GetNbins(),
      hEta3DigV->GetXaxis()->GetXmin(), hEta3DigV->GetXaxis()->GetXmax());

  float EtaTotal = hEta3DigU->Integral();  // Integral of Eta-distribution
  float EtaUTotal = hEta3DigU->Integral(); // Integral of Eta-distribution
  float EtaVTotal = hEta3DigV->Integral(); // Integral of Eta-distribution

  for (int i = 1; i <= hEta3DigU->GetXaxis()->GetNbins(); i++)
  {
    tmpEta += hEta3DigU->GetBinContent(i);
    tmpEtaU += hEta3DigU->GetBinContent(i);

    if (EtaTotal != 0.)
    {
      hEta3Int->SetBinContent(i,
                              tmpEta / EtaTotal * PixelSize - PixelSize / 2.0);
    }
    if (EtaUTotal != 0.)
    {
      hEta3UInt->SetBinContent(i, tmpEtaU / EtaUTotal * PixelSizeU -
                                      PixelSizeU / 2.0);
    }
  }

  for (int i = 1; i <= hEta3DigV->GetXaxis()->GetNbins(); i++)
  {
    tmpEtaV += hEta3DigV->GetBinContent(i);

    if (EtaVTotal != 0.)
    {
      hEta3VInt->SetBinContent(i, tmpEtaV / EtaVTotal * PixelSizeV -
                                      PixelSizeV / 2.0);
    }
  }

  // TF1 *fitEta3;
  // float FitParEta3[PolN]; // PolN=7
  if (hEta3Int->Integral() > 10)
  {
    hEta3Int->Fit("pol6", "F");
    // fitEta3 = hEta3Int->GetFunction("pol6");
    // for(int ii=0; ii<PolN; ii++){
    //   FitParEta3[ii]=fitEta3->GetParameter(ii);
    // }
    if (MimoDebug > 1)
      printf("CreateNewEta: hEta3Int fitted with %f integral\n",
             hEta3Int->Integral());
  }
  else
  {
    // for(int ii=0; ii<PolN; ii++){
    //   FitParEta3[ii]=0.;
    // }
  }
  //---------
  TF1 *fitEta3U;
  if (hEta3UInt->Integral() > 10)
  {
    hEta3UInt->Fit("pol6", "F");
    fitEta3U = hEta3UInt->GetFunction("pol6");
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParEta3U[ii] = fitEta3U->GetParameter(ii);
    }
    if (MimoDebug > 1)
      printf("CreateNewEta: hEta3UInt fitted with %f integral\n",
             hEta3UInt->Integral());
  }
  else
  {
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParEta3U[ii] = 0.;
    }
  }
  //----------
  TF1 *fitEta3V;
  if (hEta3VInt->Integral() > 10)
  {
    hEta3VInt->Fit("pol6", "F");
    fitEta3V = hEta3VInt->GetFunction("pol6");
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParEta3V[ii] = fitEta3V->GetParameter(ii);
    }
    if (MimoDebug > 1)
      printf("CreateNewEta: hEta3VInt fitted with %f integral\n",
             hEta3VInt->Integral());
  }
  else
  {
    for (int ii = 0; ii < PolN; ii++)
    {
      FitParEta3V[ii] = 0.;
    }
  }
  //----------------------
  // now for 2x2: OLD
  //----------------------
  if (MimoDebug > 1)
    printf("CreateNewEta: computing correction function for 2x2(old)\n");

  tmpEta = 0.;
  int NbOfBinInEta2x2 = hEta2x2->GetNbinsX();
  TH1F *hEta2x2IntU =
      new TH1F("hEta2x2IntU", "Eta for 2x2 cluster", NbOfBinInEta2x2, 0., 1.);
  TH1F *hEta2x2IntV =
      new TH1F("hEta2x2IntV", "Eta for 2x2 cluster", NbOfBinInEta2x2, 0., 1.);

  float *EtaContent2x2 = new float[NbOfBinInEta2x2];
  float EtaTotal2x2 = hEta2x2U->GetEntries() * 1.;

  // Integral of Eta-distribution
  for (int i = 0; i < NbOfBinInEta2x2; i++)
  {
    EtaContent2x2[i] = hEta2x2U->GetBinContent(i);
    tmpEta += EtaContent2x2[i];
    if (EtaTotal2x2 != 0.)
    {
      hEta2x2IntU->SetBinContent(i + 1, tmpEta / EtaTotal2x2);
    }
  }

  tmpEta = 0.;
  EtaTotal2x2 = hEta2x2V->GetEntries() * 1.; // Integral of Eta-distribution
  for (int i = 0; i < NbOfBinInEta2x2; i++)
  {
    EtaContent2x2[i] = hEta2x2V->GetBinContent(i);
    tmpEta += EtaContent2x2[i];
    if (EtaTotal2x2 != 0.)
    {
      hEta2x2IntV->SetBinContent(i + 1, tmpEta / EtaTotal2x2);
    }
  }

  //----------------------
  // new eta function with no fit bias. 3x3
  //----------------------
  if (MimoDebug > 1)
    printf("CreateNewEta: computing new correction functions\n");

  //  cout<<nListe_CoG<<" "<<ilist<<" "<<Liste_CoGU[ilist]<<"
  //  "<<Liste_CoGV[ilist]<<endl;

  //  TArrayF Liste_CoGU;
  //  TArrayF Liste_CoGV;
  //  int nListe_CoG;
  //---reorder both lists.
  int iscan, ilowestCoGU, ilowestCoGV;
  iscan = 0;
  ilowestCoGU = 0;
  ilowestCoGV = 0;
  float tempCoG, lowestCoGU, lowestCoGV;
  while (iscan < nListe_CoG)
  {
    lowestCoGU = 10.0 * PixelSizeU;
    lowestCoGV = 10.0 * PixelSizeV;

    for (int i = iscan; i < nListe_CoG; i++)
    {
      if (Liste_CoGU[i] < lowestCoGU)
      {
        lowestCoGU = Liste_CoGU[i];
        ilowestCoGU = i;
      }
    }
    tempCoG = Liste_CoGU[iscan];
    Liste_CoGU[iscan] = lowestCoGU;
    Liste_CoGU[ilowestCoGU] = tempCoG;

    for (int i = iscan; i < nListe_CoG; i++)
    {
      if (Liste_CoGV[i] < lowestCoGV)
      {
        lowestCoGV = Liste_CoGV[i];
        ilowestCoGV = i;
      }
    }
    tempCoG = Liste_CoGV[iscan];
    Liste_CoGV[iscan] = lowestCoGV;
    Liste_CoGV[ilowestCoGV] = tempCoG;
    iscan++;
  }

  //---etaab debut CreateNewEta()
  int iscan2x2, ilowestCoGU2x2, ilowestCoGV2x2;
  iscan2x2 = 0;
  ilowestCoGU2x2 = 0;
  ilowestCoGV2x2 = 0;
  float tempCoG2x2, lowestCoGU2x2, lowestCoGV2x2;
  while (iscan2x2 < nListe_CoG_eta2x2)
  {
    lowestCoGU2x2 = 10.0 * PixelSizeU;
    lowestCoGV2x2 = 10.0 * PixelSizeV;

    for (int i = iscan2x2; i < nListe_CoG_eta2x2; i++)
    {
      if (Liste_CoGU_eta2x2[i] < lowestCoGU2x2)
      {
        lowestCoGU2x2 = Liste_CoGU_eta2x2[i];
        ilowestCoGU2x2 = i;
      }
    }
    tempCoG2x2 = Liste_CoGU_eta2x2[iscan2x2];
    Liste_CoGU_eta2x2[iscan2x2] = lowestCoGU2x2;
    Liste_CoGU_eta2x2[ilowestCoGU2x2] = tempCoG2x2;

    for (int i = iscan2x2; i < nListe_CoG_eta2x2; i++)
    {
      if (Liste_CoGV_eta2x2[i] < lowestCoGV2x2)
      {
        lowestCoGV2x2 = Liste_CoGV_eta2x2[i];
        ilowestCoGV2x2 = i;
      }
    }
    tempCoG2x2 = Liste_CoGV_eta2x2[iscan2x2];
    Liste_CoGV_eta2x2[iscan2x2] = lowestCoGV2x2;
    Liste_CoGV_eta2x2[ilowestCoGV2x2] = tempCoG2x2;
    iscan2x2++;
  }

  int iscan5x5, ilowestCoGU5x5, ilowestCoGV5x5;
  iscan5x5 = 0;
  ilowestCoGU5x5 = 0;
  ilowestCoGV5x5 = 0;
  float tempCoG5x5, lowestCoGU5x5, lowestCoGV5x5;
  while (iscan5x5 < nListe_CoG_eta5x5)
  {
    lowestCoGU5x5 = 10.0 * PixelSizeU;
    lowestCoGV5x5 = 10.0 * PixelSizeV;

    for (int i = iscan5x5; i < nListe_CoG_eta5x5; i++)
    {
      if (Liste_CoGU_eta5x5[i] < lowestCoGU5x5)
      {
        lowestCoGU5x5 = Liste_CoGU_eta5x5[i];
        ilowestCoGU5x5 = i;
      }
    }
    tempCoG5x5 = Liste_CoGU_eta5x5[iscan5x5];
    Liste_CoGU_eta5x5[iscan5x5] = lowestCoGU5x5;
    Liste_CoGU_eta5x5[ilowestCoGU5x5] = tempCoG5x5;

    for (int i = iscan5x5; i < nListe_CoG_eta5x5; i++)
    {
      if (Liste_CoGV_eta5x5[i] < lowestCoGV5x5)
      {
        lowestCoGV5x5 = Liste_CoGV_eta5x5[i];
        ilowestCoGV5x5 = i;
      }
    }
    tempCoG5x5 = Liste_CoGV_eta5x5[iscan5x5];
    Liste_CoGV_eta5x5[iscan5x5] = lowestCoGV5x5;
    Liste_CoGV_eta5x5[ilowestCoGV5x5] = tempCoG5x5;
    iscan5x5++;
  }

  //---etaab fin CreateNewEta()

  /*  for(int i=0; i<  nListe_CoG; i++){

  cout<<nListe_CoG<<" "<<i<<" "<<Liste_CoGU[i]<<" "<<Liste_CoGV[i]<<endl;
  }
  */
  //  TArrayF STORE_Liste_CoGU(Liste_CoGU);
  //  TArrayF STORE_Liste_CoGV(Liste_CoGV);
  //  STORE_Liste_CoGU->Write();
  //  STORE_Liste_CoGV->Write();
  float n_EtaU;
  float n_EtaV;
  int n_Entries;
  n_Entries = Liste_CoGU.GetSize();

  TTree *TreeEta = new TTree("TreeEta", "Eta tree");
  TreeEta->Branch("n_EtaU", &n_EtaU, "n_EtaU/F");
  TreeEta->Branch("n_EtaV", &n_EtaV, "n_EtaV/F");
  //---etaab debut CreateNewEta()
  float n_Eta2x2U;
  float n_Eta2x2V;
  float n_Eta5x5U;
  float n_Eta5x5V;
  TreeEta->Branch("n_Eta2x2U", &n_Eta2x2U, "n_Eta2x2U/F");
  TreeEta->Branch("n_Eta2x2V", &n_Eta2x2V, "n_Eta2x2V/F");
  TreeEta->Branch("n_Eta5x5U", &n_Eta5x5U, "n_Eta5x5U/F");
  TreeEta->Branch("n_Eta5x5V", &n_Eta5x5V, "n_Eta5x5V/F");

  for (int ii = 0; ii < n_Entries; ii++)
  {
    /*  cout<<"n_Entries "<<ii<<" "<<n_Entries<<" "
  <<Liste_CoGU[ii]<<" "<<Liste_CoGV[ii]<<" "
  <<Liste_CoGU_eta2x2[ii]<<" "<<Liste_CoGV_eta2x2[ii]<<" "
  <<Liste_CoGU_eta5x5[ii]<<" "<<Liste_CoGV_eta5x5[ii]<<" "
  <<endl;*/
    n_EtaU = Liste_CoGU[ii];
    n_EtaV = Liste_CoGV[ii];
    n_Eta2x2U = Liste_CoGU_eta2x2[ii];
    n_Eta2x2V = Liste_CoGV_eta2x2[ii];
    n_Eta5x5U = Liste_CoGU_eta5x5[ii];
    n_Eta5x5V = Liste_CoGV_eta5x5[ii];
    TreeEta->Fill();
  }
  //---etaab fin CreateNewEta()
  // T->Print();

  //----------------------
  // display
  //----------------------
  // display redirected to specific canvas and not MainCanvas anymore, JB
  // 2012/09/05
  if (MimoDebug > 1)
    printf("CreateNewEta: display\n");

  TCanvas *cEta =
      new TCanvas("cEta", "Eta correction distribution", 70, 70, 750, 550);
  cEta->Draw();
  cEta->Divide(2, 2);

  cEta->cd(1);
  hEta2x2U->Draw();
  cEta->cd(3);
  hEta2x2IntU->Draw();
  cEta->cd(2);
  hEta2x2V->Draw();
  cEta->cd(4);
  hEta2x2IntV->Draw();

  cEta->Update();

  //----------------------
  // write
  // Use option kOverwrite to keep only a single key of each object (JB
  // 2013/09/19)
  //----------------------
  if (MimoDebug > 1)
    printf("CreateNewEta: saving\n");

  if (!theCorFile->IsOpen())
  {
    if (theCorFile->ReOpen("UPDATE") == -1)
    { // test added, JB 2011/04/12
      Fatal("CreateNewEta", "Cannot reopen CorPar file!");
    }
  }
  theCorFile->cd();

  TreeEta->Write("TreeEta", kOverwrite);

  ProfACGn->Write("ProfACG", kOverwrite);
  ProfUCG->Write("ProfUCG", kOverwrite);
  ProfVCG->Write("ProfVCG", kOverwrite);

  hEta3Int->Write("hEta3Int", kOverwrite);
  hEta3UInt->Write("hEta3UInt", kOverwrite);
  hEta3VInt->Write("hEta3VInt", kOverwrite);

  hEta2x2IntU->Write("hEta2x2IntU", kOverwrite);
  hEta2x2IntV->Write("hEta2x2IntV", kOverwrite);

  cEta->Write("cEta", kOverwrite);

  Info("CreateNewEta", "New Eta parameters stored into %s, enjoy!",
       theCorFile->GetName());
}

//______________________________________________________________________________
//
void MimosaAnalysis::HotPixel_init(int useMap)
{
  // Initialize the tools need for the hot pixel map.
  // If useMap = 0, the map is not used,
  //    useMap = 2, the map is generated,
  //    useMap = 1, the map is used (after it has been generated).
  //
  // Note the two cut values are hardcoded here for the moment.
  //
  // JB 2011/11/23 imported from A. Besson june 2004

  TheUsePixelMap = (useMap > 0);
  Option_read_Pixel_map = (useMap == 1);

  // If the rate limits are still 0, use the hardcoded initialization below,
  // otherwise keep the previous value.
  if (CUT_MaxHitRatePerPixel < 1.e-3)
  {
    CUT_MaxHitRatePerPixel = 0.05; // if the pixel is a seed too many times in
                                   // the run, remove the hit.
    CUT_MinHitRatePerPixel =
        0.; // you can also remove pixels with low occupancy for testing
  }

  sprintf(HotPixelFileName, "hotPixelMap_run%d_pl%d.root", RunNumber,
          ThePlaneNumber);
  sprintf(HotPixelFileName, "%s", fTool.LocalizeDirName(HotPixelFileName));
  gSystem->cd(CreateGlobalResultDir().c_str());

  if (TheUsePixelMap && Option_read_Pixel_map == 1)
  {
    HotPixelFile = new TFile(HotPixelFileName, "READ");
    h2HotPixelMap = (TH2F *)HotPixelFile->Get("h2HotPixelMap");
  }

  Info("HotPixel_init", "End hot pixel map preparation, usage = %d and %d.",
       TheUsePixelMap, Option_read_Pixel_map);
}

//______________________________________________________________________________
//
int MimosaAnalysis::HotPixel_test(int aPixelIndex)
{
  // Return if the pixel is considered as "hot" (1) or not (0).
  // Always returns 0 if TheUsePixelMap==0.
  //
  // JB 2011/11/23 imported from A. Besson june 2004

  int iRow = aPixelIndex % NofPixelInRaw + 1;
  int iCol = aPixelIndex / NofPixelInRaw + 1;
  if (TheUsePixelMap == 1 && Option_read_Pixel_map == 1 &&
      (h2HotPixelMap->GetBinContent(iRow, iCol) > CUT_MaxHitRatePerPixel ||
       h2HotPixelMap->GetBinContent(iRow, iCol) < CUT_MinHitRatePerPixel))
  {
    return 1;
  }

  return 0;
}

//______________________________________________________________________________
//
void MimosaAnalysis::HotPixel_end(int eventsRead)
{
  // If the hot pixel map was not used, do nothing.
  // If the hot pixel map was only read, simply close the file.
  // If the hot pixel map was generated, save the histogram and close the file.
  //
  // JB 2011/11/23 imported from A. Besson june 2004

  if (TheUsePixelMap)
  { // if map used

    if (Option_read_Pixel_map == 1)
    { // map only read
      cout << "-------- The HOT PIXEL MAP HAS BEEN READ ! KEEP "
           << CUT_MinHitRatePerPixel << " < rate < " << CUT_MaxHitRatePerPixel
           << "." << endl;
    }

    else
    { // map generated
      HotPixelFile = new TFile(HotPixelFileName, "RECREATE");
      h2HotPixelMap = (TH2F *)h2DgoodSeedPixel->Clone("h2HotPixelMap");
      h2HotPixelMap->Scale(1. / eventsRead);
      h2HotPixelMap->Write();
      cout << "-------- The HOT PIXEL MAP HAS BEEN FILLED !" << endl;
    }

    HotPixelFile->Close();
  } // end if map used

  else
  { // if map not used
    cout << "-------- The HOT HIT MAP WAS NOT USED !" << endl;
  }
}

//_____________________________________________________________________________
//
DPrecAlign *MimosaAnalysis::AlignMimosa(int aDistance)
{
  // Align your MIMOSA chip
  //
  // Modified: JB 2009/09/01, get rid of GetParameters() call
  // Modified: JB 2010/09/08, add the passage of the Z plane coordinate
  // (alignPar(6)) Modified: JB 2011/11/25, check there are data points to fit,
  // otherwise, returns current alignemnt

  if (!CheckIfDone("init,mimosapro"))
  {
    return 0;
  }

  TVectorD alignPar(6);
  alignPar(0) = fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Position(0) +
                fSession->GetSetup()->GetPlanePar(ThePlaneNumber).AlignmentU;
  alignPar(1) = fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Position(1) +
                fSession->GetSetup()->GetPlanePar(ThePlaneNumber).AlignmentV;
  alignPar(2) = fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Tilt(0) +
                fSession->GetSetup()->GetPlanePar(ThePlaneNumber).AlignmentTilt;
  alignPar(3) = fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Tilt(2);
  alignPar(4) = fSession->GetSetup()->GetPlanePar(ThePlaneNumber).Tilt(1);
  alignPar(5) = fSession->GetSetup()
                    ->GetPlanePar(ThePlaneNumber)
                    .Position(2); // JB 2010/09/08

  // return
  // MimosaAlignAnalysis::Instance()->AlignMimosa(alignPar,GetCorParFile(),save_result,CreateGlobalResultDir());

  if (!(alignement->GetDataPoints()->First()))
  {
    cout << endl
         << "** WARNING: there is no data points to fit in the current "
            "alignment object.";
    cout << "    --> aligning MIMOSA is not possible, keeping the current "
            "parameters!"
         << endl
         << endl;
    return nullptr; // added JB 2012/08/18
  }
  else
  {
    return MimosaAlignAnalysis::Instance(fSession)->AlignMimosa(
        alignement, GetCorParFile(), CreateGlobalResultDir(),
        aDistance); // arguments for alignment changed, JB 2011/06/18
  }
}

//_____________________________________________________________________________
//
int MimosaAnalysis::GetFileNumber()
{
  // Returns the last DSF file number created or 0 if no file.
  // Try 100 file numbers before stopping.
  // Modified: JB 2011/07/07 to localize path names
  // Modified: JB 2012/11/21 to find files event if low numbers do no exist

  Char_t New_File_Name[1000];
  int FileNumber = 1;
  int missedFiles = 0;
  int Returned_FileNumber = 0;
  int exist = 1;
  while (exist)
  {
    sprintf(New_File_Name, "%srun%d_0%d.root",
            (const char *)fSession->GetSummaryFilePath(),
            fSession->GetRunNumber(), FileNumber);
    sprintf(New_File_Name, "%s",
            fTool.LocalizeDirName(New_File_Name)); // JB 2011/07/07
    if (!gSystem->AccessPathName(New_File_Name))
    {
      missedFiles = -1;
      FileNumber++;
    }
    else if (0 <= missedFiles && missedFiles < 100)
    {
      missedFiles++;
      FileNumber++;
    }
    else
    {
      exist = 0;
    }
  };

  // Check if the user asked for a specific number
  if (!fUserFileNumber)
  {
    Returned_FileNumber = FileNumber - 1;
  }
  else
  {
    Returned_FileNumber = TMath::Min(FileNumber - 1, fUserFileNumber);
  }

  // Check the file exists
  sprintf(New_File_Name, "%srun%d_0%d.root",
          (const char *)fSession->GetSummaryFilePath(),
          fSession->GetRunNumber(), Returned_FileNumber);
  sprintf(New_File_Name, "%s", fTool.LocalizeDirName(New_File_Name));
  if (gSystem->AccessPathName(New_File_Name))
  {
    // Info( "GetFileNumber", "WARNING: file %s does not exist!",
    // New_File_Name);
    Returned_FileNumber = 0;
  }

  return Returned_FileNumber;
}

//_____________________________________________________________________________
//
bool MimosaAnalysis::CheckIfDone(const Option_t *Option)
{

  // Checks if the requiered methods have been called before processing
  //
  // Modified: JB 2010/07/27 for Fakerate and Calibration
  // Modified: JB 2010/08/30 for MiniVectos
  // Modified: JB 2010/09/10 for general Mimosa flag

  std::string opts{Option};
  std::transform(opts.begin(), opts.end(), opts.begin(), ::tolower);
  auto aReturnCode{true};

  // Unordered map to store the mapping between options and their associated
  // checks
  const std::unordered_map<std::string,
                           std::pair<Bool_t MimosaAnalysis::*, std::string>>
      checks = {
          {"clear", {&MimosaAnalysis::fClearDone, "Please run Clear() First"}},
          {"init",
           {&MimosaAnalysis::fInitDone, "Please run InitSession() First"}},
          {"mimosall",
           {&MimosaAnalysis::fMimosaDone,
            "Please run an analysis method MimosaXXX() First"}},
          {"mimosapro",
           {&MimosaAnalysis::fMimosaProDone, "Please run MimosaPro() First"}},
          {"mimosafakerate",
           {&MimosaAnalysis::fMimosaFakerateDone,
            "Please run MimosaFakerate() First"}},
          {"mimosacalibration",
           {&MimosaAnalysis::fMimosaCalibrationDone,
            "Please run MimosaCalibration() First"}},
          {"mimosaminivectors",
           {&MimosaAnalysis::fMimosaMiniVectorsDone,
            "Please run MimosaMiniVectors() First"}},
          {"mimosapro2planes",
           {&MimosaAnalysis::fMimosaPro2PlanesDone,
            "Please run MimosaPro2Planes() First"}},
          {"mimosaimaging",
           {&MimosaAnalysis::fMimosaImagingDone,
            "Please run MimosaImaging() First"}},
          {"MCGeneration",
           {&MimosaAnalysis::fIfMCGeneration,
            "Please run MimosaGeneration() First"}}};

  // Iterate through each option in the map
  for (const auto &check : checks)
  {
    if (opts.find(check.first) == std::string::npos ||
        (this->*(check.second.first))) // Option is no part of the given string
                                       // or corresponding flag is true
    {
      continue;
    }
    Warning("CheckIfDone()", "%s",
            check.second.second.c_str()); // Print the warning message
    aReturnCode = false;
  }

  return aReturnCode;
}

//______________________________________________________________________________
//

bool MimosaAnalysis::TrackInMimo(int aGeoMatrix, float tuVal, float tvVal,
                                 int aSubMatrix)
{

  // Decide wether or not the track is in the area of interest of the DUT
  // this area is defined in the config file
  //
  // JB 2009/08/26
  // Modified: JB 2012/08/20 to allow any geomatrix
  // Modified: JB 2015/12/15 to allow any submatrix

  if (aGeoMatrix < 0 || aGeoMatrix > 3)
  {
    printf("TrackInMimo: WARNING, you asked for a geomatrix (%d) which does "
           "not exist, changed to 0\n",
           aGeoMatrix);
    aGeoMatrix = 0;
  }

  // The aGeomatrix may be different from the global one (ThesubmatrixNumber)
  // So we get the limits.
  double gumin =
      fSession->GetSetup()->GetAnalysisPar().Umin[aSubMatrix][aGeoMatrix];
  double gumax =
      fSession->GetSetup()->GetAnalysisPar().Umax[aSubMatrix][aGeoMatrix];
  double gvmin =
      fSession->GetSetup()->GetAnalysisPar().Vmin[aSubMatrix][aGeoMatrix];
  double gvmax =
      fSession->GetSetup()->GetAnalysisPar().Vmax[aSubMatrix][aGeoMatrix];

  bool TrkCrossMimo =
      gvmin < tvVal && tvVal < gvmax && gumin < tuVal && tuVal < gumax;

  if (MimoDebug > 2)
  {
    printf(" TrackInMimo: (Sub %d) (Geo %d), Is %.1f < v=%.1f < %.1f & %.1f < "
           "u=%.1f < %.1f ? %d\n",
           aSubMatrix, aGeoMatrix, gvmin, tvVal, gvmax, gumin, tuVal, gumax,
           TrkCrossMimo);
  }

  return TrkCrossMimo;
}

//_____________________________________________________________________________
//
void MimosaAnalysis::RefPlane_fill()
{
  // Fill observables relevant for reference plane study
  //
  // JB 2010/05/27 based on MPro original stuff
}

//_____________________________________________________________________________
//
void MimosaAnalysis::GoodHit_Fill(DAuthenticHit *ahit)
{
  // Fill the histograms for GOOD hits
  // which are not yet selected or matched to a track
  //
  // JB 2010/10/06

  hnpix->Fill(ahit->HNNS);

  hgoodSNneighbour->Fill(ahit->HSNneighbour);

  h2dgoodhits->Fill(ahit->Hsu, ahit->Hsv);
  h2DgoodSeedPixel->Fill(ahit->Hsk % NofPixelInRaw, ahit->Hsk / NofPixelInRaw);
  hgoodSeedPixel->Fill(ahit->Hsk);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterShape_init()
{
  // Init the variables required for cluster shape studies
  //
  // Called once for all
  //
  // JB 2010/04/13 based on cdritsa original stuff
  // Modified clm 2013/01/23 compute threshold according to noiseMPV input
  // Modified clm 2013/08/25 new histos

  hitCounter = 0;

  nThresholds = 8; // shall be strictly below 11

  snrThresholds[0] = 3;
  snrThresholds[1] = 5;
  snrThresholds[2] = 6;
  snrThresholds[3] = 7;
  snrThresholds[4] = 8;
  snrThresholds[5] = 10;
  snrThresholds[6] = 15;
  snrThresholds[7] = 20;
  snrThresholds[8] = 25;
  snrThresholds[9] = 30;

  // computation of threshold, clm 2013/01/23
  double noiseMPV = -1;
  if (noiseMPV == -1)
  {
    thresholds[0] = 15;
    thresholds[1] = 30;
    thresholds[2] = 45;
    thresholds[3] = 60;
    thresholds[4] = 75;
    thresholds[5] = 90;
    thresholds[6] = 105;
    thresholds[7] = 120;
    thresholds[8] = 150;
    thresholds[9] = 220;
  }
  else
  {
    for (int i = 0; i < 9; i++)
      thresholds[i] = noiseMPV * snrThresholds[i]; // clm 2013.03.19
  }

  for (int ithres = 0; ithres < nThresholds; ithres++)
  {
    clusterMultiplicity[ithres] = 0;
    clusterMultiplicityCounter[ithres] = 0;
    multV[ithres] = 0.;
    multU[ithres] = 0.;
    mult[ithres] = 0.;
    formFactor[ithres] = 0.;
  }
  for (int i = 0; i < 50; i++)
  {
    ChargeTimesPositionU[i] = 0.;
    ChargeTimesPositionV[i] = 0.;
    ChargeTimesPositionUmean[i] = 0.;
    ChargeTimesPositionVmean[i] = 0.;
    ChargeInCluster[i] = 0.;
    sumChargeTimesPositionU[i] = 0.;
    sumChargeTimesPositionV[i] = 0.;
    sumChargeTimesPositionUmean[i] = 0.;
    sumChargeTimesPositionVmean[i] = 0.;
  }

  RMSV = 0;
  RMSU = 0;
  RMSVmean = 0;
  RMSUmean = 0;

  Char_t title[200];
  // Change the title of some histo now that we know the thresholds
  for (int i = 0; i < nThresholds; i++)
  {

    sprintf(title, "cluster multiplicity with Charge>%.1f e", thresholds[i]);
    hCountPixels[i]->SetTitle(title);

    sprintf(title, "Percentage of pixels with SNR>%.1f", snrThresholds[i]);
    hPixelsOverSNR[i]->SetTitle(title);

    sprintf(title, "Percentage of pixels with Charge>%.1f", thresholds[i]);
    hPixelsOverCharge[i]->SetTitle(title);

    sprintf(title, "Form factor vs cluster multiplicity with SNR>%.1f)",
            snrThresholds[i]);
    hMultVsFormFactor[i]->SetTitle(title);

    sprintf(title, "cluster multiplicity with SNR>%.1f", snrThresholds[i]);
    hMultVsFormFactor1D[i]->SetTitle(
        title); // clm added back 1D clu mult. 2013.08.25
  }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterShape_fill(DAuthenticHit *thehit)
{
  // Compute the variables required for cluster shape studies
  // and fill the corresponding histograms
  //
  // Called for each cluster
  //
  // JB 2010/04/13 based on cdritsa original stuff
  // Modified: JB 2010/07/28 to use charge ordered pixel

  hitCounter++;

  // Init
  for (int ithres = 0; ithres < nThresholds; ithres++)
  {
    clusterMultiplicityCounter[ithres] = 0;
    multV[ithres] = 0.;
    multU[ithres] = 0.;
    mult[ithres] = 0.;
    formFactor[ithres] = 0.;
  }

  //===========
  for (int i = 0; i < NofPixelsInCluster; i++)
  { // loop on pixels in cluster

    // poisition of pixel relative to seed in cluster
    hClusterTest1->Fill(LineInClusterqordered[i]);
    hClusterTest2->Fill(ColumnInClusterqordered[i]);
    hClusterMeanForm->Fill(ColumnInClusterqordered[i], LineInClusterqordered[i],
                           qqordered[i]);

    if (LineInClusterqordered[i] == 0)
    { // for the line of the seed only
      hChargeDistrInLine->Fill(ColumnInClusterqordered[i], qqordered[i]);
      if (ColumnInClusterqordered[i] == -3)
      {
        hChargeDistrIn3rdLeftNeigh->Fill(qqordered[i]);
      }
      else if (ColumnInClusterqordered[i] == -2)
      {
        hChargeDistrIn2ndLeftNeigh->Fill(qqordered[i]);
      }
      else if (ColumnInClusterqordered[i] == -1)
      {
        hChargeDistrIn1stLeftNeigh->Fill(qqordered[i]);
      }
      else if (ColumnInClusterqordered[i] == 0)
      {
        hChargeDistrInSeed->Fill(qqordered[i]);
      }
      else if (ColumnInClusterqordered[i] == 1)
      {
        hChargeDistrIn1stRightNeigh->Fill(qqordered[i]);
      }
      else if (ColumnInClusterqordered[i] == 2)
      {
        hChargeDistrIn2ndRightNeigh->Fill(qqordered[i]);
      }
      else if (ColumnInClusterqordered[i] == 3)
      {
        hChargeDistrIn3rdRightNeigh->Fill(qqordered[i]);
      }
    }

    // Visualise the 6 first individual clusters
    if (hitCounter <= 6)
    {
      Cluster[hitCounter - 1]->Fill(ColumnInClusterqordered[i],
                                    LineInClusterqordered[i], qqordered[i]);
    }

    // ?
    if (i > 1 && abs(ColumnInClusterqordered[i]) < 6 &&
        abs(LineInClusterqordered[i]) < 6)
    {
      ChargeInCluster[i] = ChargeInCluster[i - 1] + qqordered[i];
      ChargeTimesPositionU[i] = ChargeTimesPositionU[i - 1] +
                                qqordered[i] * ColumnInClusterqordered[i] *
                                    ColumnInClusterqordered[i];
      ChargeTimesPositionV[i] =
          ChargeTimesPositionV[i - 1] +
          qqordered[i] * LineInClusterqordered[i] * LineInClusterqordered[i];
      ChargeTimesPositionUmean[i] = ChargeTimesPositionUmean[i - 1] +
                                    qqordered[i] * ColumnInClusterqordered[i];
      ChargeTimesPositionVmean[i] = ChargeTimesPositionVmean[i - 1] +
                                    qqordered[i] * LineInClusterqordered[i];
    }
    else
    { // i==0
      ChargeInCluster[i] = qqordered[i];
      ChargeTimesPositionU[i] = qqordered[i] * ColumnInClusterqordered[i] *
                                ColumnInClusterqordered[i];
      ChargeTimesPositionV[i] =
          qqordered[i] * LineInClusterqordered[i] * LineInClusterqordered[i];
      ChargeTimesPositionUmean[i] = qqordered[i] * ColumnInClusterqordered[i];
      ChargeTimesPositionVmean[i] = qqordered[i] * LineInClusterqordered[i];
    }

    // cout << "ChargeTimesPositionU[" << i << "] = " << ChargeTimesPositionU[i]
    // << endl; cout << "ChargeTimesPositionV[" << i << "] = " <<
    // ChargeTimesPositionV[i] << endl; cout << "ChargeInCluster[" << i << "] =
    // " << ChargeInCluster[i] << endl;

    // average cluster multiplicity: how many pixels are over threshold in
    // average?
    for (int ithres = 0; ithres < nThresholds; ithres++)
    { // loop on thresholds

      // Threshold = Q
      // +3 is added because histo must have x-axis between (0,7) .
      if (qqordered[i] > thresholds[ithres])
      {
        hPixelsOverCharge[ithres]->Fill(ColumnInClusterqordered[i],
                                        LineInClusterqordered[i], 100);
        clusterMultiplicity[ithres] += 1;
        clusterMultiplicityCounter[ithres] +=
            1; // counter is for absolute count of pixels above thresh; no
               // average.
      }
      else
      {
        hPixelsOverCharge[ithres]->Fill(ColumnInClusterqordered[i],
                                        LineInClusterqordered[i], 0);
      }

      // Threshold = SNR
      if (snqordered[i] > snrThresholds[ithres])
      {
        hPixelsOverSNR[ithres]->Fill(ColumnInClusterqordered[i],
                                     LineInClusterqordered[i], 100);
        mult[ithres] += 1;
        if (ColumnInClusterqordered[i] == 0)
        {
          multV[ithres]++;
        }
        if (LineInClusterqordered[i] == 0)
        {
          multU[ithres]++;
        }
      }
      else
      {
        hPixelsOverSNR[ithres]->Fill(ColumnInClusterqordered[i],
                                     LineInClusterqordered[i], 0);
      }

    } // end loop on thresholds

  } // end loop on pixels in cluster
  //===========

  //===========
  // Fill histograms

  // JB 2014/03/31
  hClusterSizeInLines->Fill(LineSizeOfCluster);
  hClusterSizeInColumns->Fill(ColumnSizeOfCluster);
  hClusterSize2DInLineAndColumns->Fill(LineSizeOfCluster, ColumnSizeOfCluster);
  hClusterTypes->Fill(static_cast<int>(ClusterGeometricalType));
  if (NofPixelsInCluster >= 5)
    hClusterTypesBeyond4->Fill(static_cast<int>(ClusterGeometricalTypeBeyond4) -
                               static_cast<int>(MAnalysis::ClusterType::more_than_four));

  for (int ithres = 0; ithres < nThresholds; ithres++)
  { // loop on thresholds

    hCountPixels[ithres]->Fill(clusterMultiplicityCounter[ithres]);

    if (multV[ithres] != 0)
    {
      formFactor[ithres] = multU[ithres] / multV[ithres];
    }
    else
    {
      formFactor[ithres] = 0;
    }

    hMultVsFormFactor[ithres]->Fill(formFactor[ithres], mult[ithres]);
    hMultVsFormFactor1D[ithres]->Fill(
        mult[ithres]); // clm added back 1D clu mult. 2013.08.25

  } // end loop on thresholds

  if (ChargeInCluster[25] != 0)
  {
    RMSV = (ChargeTimesPositionV[25] / ChargeInCluster[25]);
    RMSU = (ChargeTimesPositionU[25] / ChargeInCluster[25]);
  }
  else
  {
    RMSV = 0.;
    RMSU = 0.;
  }

  if (ChargeInCluster[24] != 0)
  {
    RMSVmean = (ChargeTimesPositionVmean[24] / ChargeInCluster[24]) *
               (ChargeTimesPositionVmean[24] / ChargeInCluster[24]);
    RMSUmean = (ChargeTimesPositionUmean[24] / ChargeInCluster[24]) *
               (ChargeTimesPositionUmean[24] / ChargeInCluster[24]);
  }
  else
  {
    RMSVmean = 0.;
    RMSUmean = 0.;
  }

  if ((qqordered[0] / nqordered[0]) > 10)
  {
    hChargeCoG_Correl->Fill(sqrt(RMSV), sqrt(RMSU));
  }
  hChargeCoG_Correl2->Fill(sqrt(RMSV - RMSVmean), sqrt(RMSU - RMSUmean));

  //===========
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterShape_end()
{
  // Finalize histograms created and filled for cluster shape studies
  //
  // Called once for all
  //
  // JB 2010/04/13 based on cdritsa original stuff

  cout << "******* hitCounter for cluster shape study ******* " << hitCounter
       << endl;

  if (!hitCounter)
    return; // leave if no hits were selected

  double normfactor = (double)1 / hitCounter;

  hClusterMeanForm->SetEntries(hitCounter);
  hClusterMeanForm->Scale(normfactor);
  // projections
  hprojLForm = hClusterMeanForm->ProjectionY("ClShapeOnLines");
  hprojLForm->SetTitle("Cluster Average Shape projected / lines");
  hprojCForm = hClusterMeanForm->ProjectionX("ClShapeOnColumns");
  hprojCForm->SetTitle("Cluster Average Shape projected / columns");

  for (int ithres = 0; ithres < nThresholds; ithres++)
  { // loop on thresholds

    hPixelsOverSNR[ithres]->Scale(normfactor);
    hPixelsOverSNR[ithres]->SetEntries(hitCounter);

    hPixelsOverCharge[ithres]->Scale(normfactor);

    Char_t text[20];
    /** Project a line **/
    sprintf(text, "hprojL%d", ithres);
    hprojL[ithres] = hPixelsOverCharge[ithres]->ProjectionX(text, 4, 4);
    hprojL[ithres]->SetEntries(hitCounter);

    /** Project a column **/
    sprintf(text, "hprojC%d", ithres);
    hprojC[ithres] = hPixelsOverCharge[ithres]->ProjectionY(text, 4, 4);
    hprojC[ithres]->SetEntries(hitCounter);

  } // end loop on thresholds
}

//_____________________________________________________________________________
//
void MimosaAnalysis::BinarySpecific_HitsOnly_fill(DAuthenticHit *ahit)
{
  // Compute the variables required for binary output sensors
  // and fill the corresponding histograms
  //
  // Called for each selected cluster
  //
  // JB 2010/06/03 based on A.Besson original stuff

  float Charge_all_m8 = 0.0;
  float Charge_sub_m8 = 0.0;
  float hUdigital_m8 = 0.0;
  float hVdigital_m8 = 0.0;

  // Count the number of pixels inside the hit
  for (int i = 0; i < NofPixelsInCluster; i++)
  {
    Charge_all_m8 += ahit->HqM[i];
    hVdigital_m8 =
        ((ahit->Hsk) / NofPixelInRaw) * PixelSizeV - FirstPixelShiftV;
    hUdigital_m8 =
        ((ahit->Hsk) % NofPixelInRaw) * PixelSizeU - FirstPixelShiftU;
    // use track in mimo to see if the hit is inside the window.
    if (TrackInMimo(Thegeomatrix, hUdigital_m8, hVdigital_m8,
                    ThesubmatrixNumber))
    {
      Charge_sub_m8 += ahit->HqM[i];
    }
  }

  hBinary_NumberOf_1_ALL->Fill(Charge_all_m8);       // all hits
  hBinary_NumberOf_1_submatrix->Fill(Charge_sub_m8); // hits in the geomatrix
}

//_____________________________________________________________________________
//
void MimosaAnalysis::BinarySpecific_fill(DAuthenticHit *thehit)
{
  // Compute the variables required for binary output sensors
  // and fill the corresponding histograms
  //
  // Called for each cluster associated to a track
  //
  // JB 2010/06/03 based on A.Besson original stuff

  hBinary_NumberOf_1_goodhit->Fill(qonly[NofPixelsInCluster - 1]);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::BinarySpecific_end()
{
  // Finalize histograms created and filled for binary output sensors
  //
  // Called once for all
  //
  // JB 2010/06/03 based on A.Besson original stuff

  if (MimoDebug)
    Info("BinarySpecific_end", "Finalizing histograms");

  if (NtrkInMimo != 0)
  {

    for (int i = 1; i < NofPixelInColumn + 1; i++)
    {
      for (int j = 1; j < NofPixelInRaw + 1; j++)
      {
        hBinary_Nhitperpixel_submatrix->Fill(
            h2DgoodSeedPixel->GetBinContent(j, i));
        hBinary_NhitRateperpixel_submatrix->Fill(
            h2DgoodSeedPixel->GetBinContent(j, i) / float(NtrkInMimo));
      }
    }
  }

  float hM8_sub_integ;
  hM8_sub_integ = hBinary_NhitRateperpixel_submatrix->Integral();
  if (hM8_sub_integ != 0.0)
  {
    hBinary_NhitRateperpixel_submatrix->Scale(1 / hM8_sub_integ);
  }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterCharges_init()
{
  // Initialize usefull cluster properties for further studies
  //
  // Called once for all
  //
  // JB 2010/06/03 based on original MPro stuff

  q = new float[MaxNofPixelsInCluster];  // cumulative cluster Charge in
                                         // [index+1] pixels
  n = new float[MaxNofPixelsInCluster];  // cumulative cluster Noise in [index+1]
                                         // pixels
  sn = new float[MaxNofPixelsInCluster]; // cumulative signal/noise in [index+1]
                                         // pixels

  qsumsn = new float[MaxNofPixelsInCluster];
  nsumsn = new float[MaxNofPixelsInCluster];
  snsumsn = new float[MaxNofPixelsInCluster];

  qsumspiral = new float[MaxNofPixelsInCluster];
  nsumspiral = new float[MaxNofPixelsInCluster];
  snsumspiral = new float[MaxNofPixelsInCluster];

  UofPix = new float[MaxNofPixelsInCluster];        // u of pixels in a cluster
  VofPix = new float[MaxNofPixelsInCluster];        // v -------" --------
  LineInCluster = new int[MaxNofPixelsInCluster];   // line shift wrt seed
  ColumnInCluster = new int[MaxNofPixelsInCluster]; // column shift wrt seed
  IndexofPix = new int[MaxNofPixelsInCluster];      // index of pixel in the matrix
  IndexInCluster =
      new int[MaxNofPixelsInCluster]; // index of pixel relative to seed

  IndexU2x2[0] = 4; // used for Eta 2x2 computation...not clear, JB 2010/07/27
  IndexU2x2[1] = 4;
  IndexU2x2[2] = 5;
  IndexU2x2[3] = 5;
  IndexV2x2[0] = 2;
  IndexV2x2[1] = 7;
  IndexV2x2[2] = 2;
  IndexV2x2[3] = 7;

  qonly =
      new float[MaxNofPixelsInCluster]; // charge in the pixels of the cluster
  noise =
      new float[MaxNofPixelsInCluster]; // Noise in the pixels of the cluster
  snonly =
      new float[MaxNofPixelsInCluster];                // charge in the pixels of the cluster
  qqordered = new float[MaxNofPixelsInCluster];        // charge in the pixels of the
                                                       // cluster ordered by S/N
  nqordered = new float[MaxNofPixelsInCluster];        // Noise in the pixels of the
                                                       // cluster ordered by S/N
  snqordered = new float[MaxNofPixelsInCluster];       // signal/noise ordered by S/N
  IndexofPixqordered = new int[MaxNofPixelsInCluster]; //
  UofPixqordered = new float[MaxNofPixelsInCluster];   // u of pixels in a cluster
                                                       // ordered by S/N
  VofPixqordered =
      new float[MaxNofPixelsInCluster];                   // v -------" -------- ordered by S/N
  LineInClusterqordered = new int[MaxNofPixelsInCluster]; // line shift wrt seed
  ColumnInClusterqordered =
      new int[MaxNofPixelsInCluster]; // column shift wrt seed

  qsnordered = new float[MaxNofPixelsInCluster];        // charge in the pixels of the
                                                        // cluster ordered by S/N
  nsnordered = new float[MaxNofPixelsInCluster];        // Noise in the pixels of the
                                                        // cluster ordered by S/N
  snsnordered = new float[MaxNofPixelsInCluster];       // signal/noise ordered by S/N
  IndexofPixsnordered = new int[MaxNofPixelsInCluster]; //
  UofPixsnordered = new float[MaxNofPixelsInCluster];   // u of pixels in a
                                                        // cluster ordered by S/N
  VofPixsnordered =
      new float[MaxNofPixelsInCluster]; // v -------" -------- ordered by S/N
  LineInClustersnordered =
      new int[MaxNofPixelsInCluster]; // line shift wrt seed
  ColumnInClustersnordered =
      new int[MaxNofPixelsInCluster]; // column shift wrt seed

  qspiral =
      new float[MaxNofPixelsInCluster]; // charge in the pixels of the cluster
                                        // ordered by geometrical position
  nspiral =
      new float[MaxNofPixelsInCluster];              // Noise in the pixels of the cluster
                                                     // ordered by geometrical position
  snspiral = new float[MaxNofPixelsInCluster];       // signal/noise ordered by
                                                     // geometrical position
  IndexofPixspiral = new int[MaxNofPixelsInCluster]; //
  UofPixspiral =
      new float[MaxNofPixelsInCluster]; // u of pixels in a cluster ordered by
                                        // geometrical position
  VofPixspiral =
      new float[MaxNofPixelsInCluster];                 // v -------" -------- ordered by
                                                        // geometrical position
  LineInClusterspiral = new int[MaxNofPixelsInCluster]; // line shift wrt seed
  ColumnInClusterspiral =
      new int[MaxNofPixelsInCluster]; // column shift wrt seed
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterCharges_compute(DAuthenticHit *thehit)
{
  // compute usefull cluster properties for further studies
  // Note that all charges computed here include the calibration factor.
  //
  // Called for each matched clusters
  //
  // JB 2010/06/03 based on original MPro stuff
  // Modified 2010/07/27: JB, for the "spiral" geomtrical ordering
  // Modified 2010/09/06: JB, usage of NoiseScope
  // Modified 2012/11/21: JB, use ComputePixelPosition method

  if (MimoDebug)
    Info("ClusterCharges_compute",
         " Computing charges for hit %d with %d pixels\n", thehit->Hhk,
         NofPixelsInCluster);

  //------------
  //-- Initialize
  int LineInClusterSeed = (thehit->HkM[0]) / NofPixelInRaw;
  int ColumnInClusterSeed = (thehit->HkM[0]) % NofPixelInRaw;
  int highestLineInCruster = -1000; // init, JB 2014/03/31
  int lowestLineInCruster = 1000;
  int highestColumnInCruster = -1000;
  int lowestColumnInCruster = 1000;
  for (int iPix = 0; iPix < MaxNofPixelsInCluster; iPix++)
  {
    qonly[iPix] = 0.;
    noise[iPix] = 0.;
    snonly[iPix] = 0.;
    IndexofPix[iPix] = 0;
    IndexInCluster[iPix] = 0;
    UofPix[iPix] = 0.;
    VofPix[iPix] = 0.;
    LineInCluster[iPix] = 0;
    ColumnInCluster[iPix] = 0;

    qqordered[iPix] = 0.;
    nqordered[iPix] = 0.;
    snqordered[iPix] = 0.;
    IndexofPixqordered[iPix] = 0;
    UofPixqordered[iPix] = 0.;
    VofPixqordered[iPix] = 0.;
    LineInClusterqordered[iPix] = 0;
    ColumnInClusterqordered[iPix] = 0;

    q[iPix] = 0.;
    n[iPix] = 0.;
    sn[iPix] = 0.;

    qsnordered[iPix] = 0.;
    nsnordered[iPix] = 0.;
    snsnordered[iPix] = 0.;
    IndexofPixsnordered[iPix] = 0;
    UofPixsnordered[iPix] = 0.;
    VofPixsnordered[iPix] = 0.;
    LineInClustersnordered[iPix] = 0;
    ColumnInClustersnordered[iPix] = 0;

    qsumsn[iPix] = 0.;
    nsumsn[iPix] = 0.;
    snsumsn[iPix] = 0.;

    qspiral[iPix] = 0.;
    nspiral[iPix] = 0.;
    snspiral[iPix] = 0.;
    IndexofPixspiral[iPix] = 0;
    UofPixspiral[iPix] = 0.;
    VofPixspiral[iPix] = 0.;
    LineInClusterspiral[iPix] = 0;
    ColumnInClusterspiral[iPix] = 0;

    qsumspiral[iPix] = 0.;
    nsumspiral[iPix] = 0.;
    snsumspiral[iPix] = 0.;
  }

  for (short i = 0; i < 4; i++)
  {
    Qof2x2[i] = 0.;
    UofPix2x2[i] = 0.;
    VofPix2x2[i] = 0.;
  }

  for (short i = 0; i < 9; i++)
  {
    Qof3x3[i] = 0.;
    UofPix3x3[i] = 0.;
    VofPix3x3[i] = 0.;
    IndexofPix3x3[i] = 0;
  }

  for (short i = 0; i < 25; i++)
  {
    Qof5x5[i] = 0.;
    UofPix5x5[i] = 0.;
    VofPix5x5[i] = 0.;
    IndexofPix5x5[i] = 0;
  }

  TotalCharge = 0.;
  TotalCharge2x2 = 0.;
  TotalCharge3x3 = 0.;
  TotalCharge5x5 = 0.;
  ChargeAroundSeed = 0.;
  ChargeInCross = 0.;
  ChargeInX = 0.;
  ChargeInCrown1 = 0.;
  ChargeInCrown2 = 0.;
  snOptimal = 0.;
  optimalSize = 0;
  IndexOfMaxPixel = 0;
  float MinUofPix = 1.e8;
  float MinVofPix = 1.e8;

  //------------
  //-- Fill the initial arrays
  float ChargeOfMaxPixel = 0.;
  for (int iq = 0; iq < NofPixelsInCluster;
       iq++)
  { // loop on pixels inside cluster

    qonly[iq] = thehit->HqM[iq] * calibration;
    // cout << " Init charge " <<  thehit->HqM[iq] << " * calib " << calibration
    // << " = " << qonly[iq] << endl;
    noise[iq] = thehit->HnM[iq] * calibration * NoiseScope;
    if (NoiseScope < 1.e-2)
    {                   // noise not used (NoiseScope=0.), JB 2010/09/06
      noise[iq] = 1e-5; // to avoid potential division by zero
      snonly[iq] = 0.;
    }
    else
    { // noise used
      if (noise[iq] < 1.e-2)
      { // unrealistic small value
        if (MimoDebug > 1)
          printf("ClusterCharges_compute: unrealistic noise %.f set to %.f\n",
                 noise[iq], 1.e-3 * calibration * NoiseScope);
        noise[iq] = 1.e-3 * calibration * NoiseScope;
      }
      snonly[iq] = qonly[iq] / noise[iq];
    }
    IndexInCluster[iq] = thehit->HkM[iq] - thehit->HkM[0];
    // the folowing line was in the original code: wrong ? JB 2010/06/03
    // IndexInCluster[iq]=IndexInCluster[iq]- IndexInCluster[0];
    IndexofPix[iq] = thehit->HkM[iq];
    // Specific index correction for MIMOSA 22 THR, JB 2013/08/13
    // due to error in indexation in DPlane
    if (MimosaType == 22)
    {
      // IndexofPix[iq]--;
      if (IndexofPix[iq] < 0)
      {
        Warning("ClusterCharges_compute", "WARNING negative index %d, set to 0",
                IndexofPix[iq]);
        IndexofPix[iq] = 0;
      }
    }
    // UofPix[iq] = IndexofPix[iq]%NofPixelInRaw*PixelSizeU-FirstPixelShiftU;
    // VofPix[iq] = IndexofPix[iq]/NofPixelInRaw*PixelSizeV-FirstPixelShiftV;
    // cout << "Calling ComputePosition for index " << IndexofPix[iq] << endl;
    ComputePixelPosition(IndexofPix[iq] % NofPixelInRaw,
                         IndexofPix[iq] / NofPixelInRaw, UofPix[iq],
                         VofPix[iq]); // JB 2012/11/21
    LineInCluster[iq] = (thehit->HkM[iq]) / NofPixelInRaw - LineInClusterSeed;
    ColumnInCluster[iq] =
        (thehit->HkM[iq]) % NofPixelInRaw - ColumnInClusterSeed;
    if (LineInCluster[iq] < lowestLineInCruster)
      lowestLineInCruster = LineInCluster[iq];
    if (highestLineInCruster < LineInCluster[iq])
      highestLineInCruster = LineInCluster[iq];
    if (ColumnInCluster[iq] < lowestColumnInCruster)
      lowestColumnInCruster = ColumnInCluster[iq];
    if (highestColumnInCruster < ColumnInCluster[iq])
      highestColumnInCruster = ColumnInCluster[iq];

    TotalCharge += qonly[iq];

    // The pixel with the highest charge defines the seed
    // its index is stored in IndexOfMaxPixel
    if (ChargeOfMaxPixel < qonly[iq])
    {
      ChargeOfMaxPixel = qonly[iq];
      IndexOfMaxPixel = iq;
    }

    if (UofPix[iq] < MinUofPix)
      MinUofPix = UofPix[iq];
    if (VofPix[iq] < MinVofPix)
      MinVofPix = VofPix[iq];

    qqordered[iq] = qonly[iq];
    nqordered[iq] = noise[iq];
    snqordered[iq] = snonly[iq];
    IndexofPixqordered[iq] = IndexofPix[iq];
    UofPixqordered[iq] = UofPix[iq];
    VofPixqordered[iq] = VofPix[iq];
    LineInClusterqordered[iq] = LineInCluster[iq];
    ColumnInClusterqordered[iq] = ColumnInCluster[iq];

    qsnordered[iq] = qonly[iq];
    nsnordered[iq] = noise[iq];
    snsnordered[iq] = snonly[iq];
    IndexofPixsnordered[iq] = IndexofPix[iq];
    UofPixsnordered[iq] = UofPix[iq];
    VofPixsnordered[iq] = VofPix[iq];
    LineInClustersnordered[iq] = LineInCluster[iq];
    ColumnInClustersnordered[iq] = ColumnInCluster[iq];

    if (MimoDebug > 1)
      printf(" un-ordered pix %d, q=%6.1f n=%4.1f sn=%6.1f index=%d u=%6.1f "
             "v=%6.1f line/seed=%d col/seed=%d\n",
             iq, qonly[iq], noise[iq], snonly[iq], IndexofPix[iq], UofPix[iq],
             VofPix[iq], LineInCluster[iq], ColumnInCluster[iq]);

  } // end loop on pixels inside cluster

  // Size of cluster
  // JB 2014/03/31
  LineSizeOfCluster = highestLineInCruster - lowestLineInCruster + 1;
  ColumnSizeOfCluster = highestColumnInCruster - lowestColumnInCruster + 1;

  SetClusterGeometricalType(NofPixelsInCluster, LineSizeOfCluster, ColumnSizeOfCluster);

  //------------
  //-- Reorder pixels by their charge:
  // index 0 is then the "true" seed
  bool continuer = kTRUE;
  float tmpq, tmpnoise, tmpsn, tmpU, tmpV;
  int tmpi, tmplinei, tmpcolumni;
  while (continuer)
  {
    continuer = kFALSE;
    for (int iq = 0; iq < NofPixelsInCluster - 1; iq++)
    {
      if (qqordered[iq] < qqordered[iq + 1])
      {
        continuer = kTRUE;
        tmpq = qqordered[iq];
        tmpnoise = nqordered[iq];
        tmpi = IndexofPixqordered[iq];
        tmpsn = snqordered[iq];
        tmpU = UofPixqordered[iq];
        tmpV = VofPixqordered[iq];
        tmplinei = LineInClusterqordered[iq];
        tmpcolumni = ColumnInClusterqordered[iq];
        qqordered[iq] = qqordered[iq + 1];
        qqordered[iq + 1] = tmpq;
        nqordered[iq] = nqordered[iq + 1];
        nqordered[iq + 1] = tmpnoise;
        snqordered[iq] = snqordered[iq + 1];
        snqordered[iq + 1] = tmpsn;
        IndexofPixqordered[iq] = IndexofPix[iq + 1];
        IndexofPixqordered[iq + 1] = tmpi;
        UofPixqordered[iq] = UofPixqordered[iq + 1];
        UofPixqordered[iq + 1] = tmpU;
        VofPixqordered[iq] = VofPixqordered[iq + 1];
        VofPixqordered[iq + 1] = tmpV;
        LineInClusterqordered[iq] = LineInClusterqordered[iq + 1];
        LineInClusterqordered[iq + 1] = tmplinei;
        ColumnInClusterqordered[iq] = ColumnInClusterqordered[iq + 1];
        ColumnInClusterqordered[iq + 1] = tmpcolumni;
      }
    }
  }

  if (MimoDebug > 1)
  {
    for (int iq = 0; iq < NofPixelsInCluster - 1; iq++)
    {
      printf(" q-ordered pix %d, q=%6.1f n=%4.1f sn=%6.1f index=%d u=%6.1f "
             "v=%6.1f line/seed=%d col/seed=%d\n",
             iq, qqordered[iq], noise[iq], snqordered[iq],
             IndexofPixqordered[iq], UofPixqordered[iq], VofPixqordered[iq],
             LineInClusterqordered[iq], ColumnInClusterqordered[iq]);
    }
  }

  //------------
  //-- Reorder pixels by S/N
  continuer = kTRUE;
  while (continuer)
  {
    continuer = kFALSE;
    for (int iq = 0; iq < NofPixelsInCluster - 1; iq++)
    {
      if (snsnordered[iq] < snsnordered[iq + 1])
      {
        continuer = kTRUE;
        tmpq = qsnordered[iq];
        tmpnoise = nsnordered[iq];
        tmpsn = snsnordered[iq];
        tmpi = IndexofPixsnordered[iq];
        tmpU = UofPixsnordered[iq];
        tmpV = VofPixsnordered[iq];
        tmplinei = LineInClustersnordered[iq];
        tmpcolumni = ColumnInClustersnordered[iq];
        IndexofPixsnordered[iq] = IndexofPixsnordered[iq + 1];
        IndexofPixsnordered[iq + 1] = tmpi;
        qsnordered[iq] = qsnordered[iq + 1];
        qsnordered[iq + 1] = tmpq;
        nsnordered[iq] = nsnordered[iq + 1];
        nsnordered[iq + 1] = tmpnoise;
        snsnordered[iq] = snsnordered[iq + 1];
        snsnordered[iq + 1] = tmpsn;
        UofPixsnordered[iq] = UofPixsnordered[iq + 1];
        UofPixsnordered[iq + 1] = tmpU;
        VofPixsnordered[iq] = VofPixsnordered[iq + 1];
        VofPixsnordered[iq + 1] = tmpV;
        LineInClustersnordered[iq] = LineInClustersnordered[iq + 1];
        LineInClustersnordered[iq + 1] = tmplinei;
        ColumnInClustersnordered[iq] = ColumnInClustersnordered[iq + 1];
        ColumnInClustersnordered[iq + 1] = tmpcolumni;
      }
    }
  }

  //------------
  //-- Reorder pixels geometrically
  // numbering as a spiral in the trigonometric way
  // 34 33 32 31 30 29 28
  // 35 15 14 13 12 11 27
  // 36 16  4  3  2 10 26
  // 37 17  5 .0  1  9 25 etc...
  // 38 18  6  7  8 24 47
  // 39 19 20 21 22 23 46
  // 40 41 41 42 43 44 45
  // The seed is the highest charge in the center
  // Algo tested by JB 2010/07/28, see macros/spiralNumerotation.C
  qspiral[0] = qqordered[0];
  nspiral[0] = nqordered[0];
  snspiral[0] = snqordered[0];
  IndexofPixspiral[0] = IndexofPixqordered[0];
  UofPixspiral[0] = UofPixqordered[0];
  VofPixspiral[0] = VofPixqordered[0];
  double as, dist, angle;
  int nCrown, firstIndex, nofPixels, spiralIndex;
  double pi = TMath::Pi();
  for (int iq = 1; iq < NofPixelsInCluster; iq++)
  {

    // center (seed) pixel is a special case
    if (LineInClusterqordered[iq] == 0 && ColumnInClusterqordered[iq] == 0)
    {
      spiralIndex = 0;
    }

    // pixels on crowns
    else
    {
      dist = sqrt(LineInClusterqordered[iq] * LineInClusterqordered[iq] +
                  ColumnInClusterqordered[iq] *
                      ColumnInClusterqordered[iq]); // distance to seed
      nCrown = TMath::Max(abs(LineInClusterqordered[iq]),
                          abs(ColumnInClusterqordered[iq])); // crown number
      firstIndex = 1 + 4 * nCrown * (nCrown - 1);            // always odd
      nofPixels = 8 * nCrown;                                // always even
      angle = acos(LineInClusterqordered[iq] / dist);        // 0->pi
      as = asin(ColumnInClusterqordered[iq] / dist);         // -pi/2->pi/2
      if (as < 0)
      {
        angle = 2. * pi - angle;
      } // build an angle from 0->2pi
      spiralIndex = firstIndex +
                    (int)floor(nofPixels * angle / 2. / pi +
                               0.5); // +.5 is a trick to round to nearest int
      if (MimoDebug > 1)
      {
        printf("   spiral  pix=%d, dist=%6.1f n=%2d first=%3d pixels=%3d "
               "(asin=%5.1f) angle=%6.2f i=%6.2f",
               iq, dist, nCrown, firstIndex, nofPixels, as, angle,
               nofPixels * angle / 2. / pi);
        printf(" [ %2d, %2d] -> %3d\n", LineInClusterqordered[iq],
               ColumnInClusterqordered[iq], spiralIndex);
      }
    }

    if (spiralIndex < MaxNofPixelsInCluster)
    { // check the index make sense,
      // otherwise do nothing
      qspiral[spiralIndex] = qqordered[iq];
      nspiral[spiralIndex] = nqordered[iq];
      snspiral[spiralIndex] = snqordered[iq];
      IndexofPixspiral[spiralIndex] = IndexofPixqordered[iq];
      UofPixspiral[spiralIndex] = UofPixqordered[iq];
      VofPixspiral[spiralIndex] = VofPixqordered[iq];
      LineInClusterspiral[spiralIndex] = LineInClusterqordered[iq];
      ColumnInClusterspiral[spiralIndex] = ColumnInClusterqordered[iq];
    }
    else if (MimoDebug)
    {
      Info("ClusterCharges_compute",
           "For pixel %d, the spiralIndex %d from line %d and column %d / "
           "seed, overwhelms the max index %d",
           iq, spiralIndex, LineInClusterqordered[iq],
           ColumnInClusterqordered[iq], MaxNofPixelsInCluster);
    }
  }

  //------------
  //-- Reorder pixels by their distance from center:
  /*
    bool continuer = kTRUE ;
    while (continuer) {
    continuer = kFALSE  ;
    for (int iq=0 ; iq<NofPixelsInCluster-2 ; iq++){
    float tmp,tmpnoise;
    int itmp;
    //		  cout << "d[" << iq<< "]=" <<
    sqrt((hUdigital-UofPix[iq])*(hUdigital-UofPix[iq]) +
    //
    (hVdigital-VofPix[iq])*(hVdigital-VofPix[iq])) << endl;
    //		  cout << "d["<< iq+1<< "]=" <<
    sqrt((hUdigital-UofPix[iq+1])*(hUdigital-UofPix[iq+1]) +
    //
    (hVdigital-VofPix[iq+1])*(hVdigital-VofPix[iq+1])) << endl; if
    (sqrt((hUdigital-UofPix[iq])*(hUdigital-UofPix[iq]) +
    (hVdigital-VofPix[iq])*(hVdigital-VofPix[iq]))>
    sqrt((hUdigital-UofPix[iq+1])*(hUdigital-UofPix[iq+1]) +
    (hVdigital-VofPix[iq+1])*(hVdigital-VofPix[iq+1]))
    ){
    continuer = kTRUE ;
    tmp = q[iq] ;
    tmpnoise = noise[iq];
    itmp = IndexInCluster[iq];
    //	    cout << "index= " << IndexInCluster[iq] << endl;
    IndexInCluster[iq]=IndexInCluster[iq+1];
    IndexInCluster[iq+1]=itmp;
    q[iq]  = q[iq+1] ;
    q[iq+1]= tmp     ;
    noise[iq] = noise[iq+1];
    noise[iq+1] = tmpnoise;
    tmp=UofPix[iq];
    UofPix[iq]=UofPix[iq+1];
    UofPix[iq+1]=tmp;
    tmp=VofPix[iq];
    VofPix[iq]=VofPix[iq+1];
    VofPix[iq+1]=tmp;
    //		    cout << " ------- ievt= " << ievt << endl;
    //		    cout << "Index=" << IndexInCluster[iq] << "q[" << iq << "]="
    << q[iq] << endl;
    //		    cout << "Index=" << IndexInCluster[iq+1] << "   q[" << iq+1
    << "]=" << q[iq+1] << endl;
    }
    }
    }
  */

  //------------
  // compute cumulate charge from ordered pixels (with charge, s/n and geometry)
  // noise is added in quadrature
  // also determine "optimal cluster" with the largest S/N

  q[0] = qqordered[0]; // charge ordered
  n[0] = nqordered[0];
  sn[0] = snqordered[0];
  qsumsn[0] = qsnordered[0]; // s/n ordered
  nsumsn[0] = nsnordered[0];
  snsumsn[0] = snsnordered[0];
  qsumspiral[0] = qspiral[0]; // geometry (spiral numerotation) ordered
  nsumspiral[0] = nspiral[0];
  snsumspiral[0] = snspiral[0];
  snOptimal = sn[0]; // init largest s/n cluster
  optimalSize = 1;
  for (int iq = 1; iq < NofPixelsInCluster; iq++)
  { // loop on pixels
    q[iq] = q[iq - 1] + qqordered[iq];
    n[iq] = sqrt(n[iq - 1] * n[iq - 1] + nqordered[iq] * nqordered[iq]);
    if (n[iq] > 1.e-2)
    {
      sn[iq] = q[iq] / n[iq];
    }
    else
    {
      sn[iq] = 0.;
    }
    if (sn[iq] > snOptimal)
    {
      snOptimal = sn[iq];
      optimalSize = iq + 1;
    }
    qsumsn[iq] += qsnordered[iq];
    nsumsn[iq] =
        sqrt(nsumsn[iq - 1] * nsumsn[iq - 1] + nsnordered[iq] * nsnordered[iq]);
    if (nsumsn[iq] > 1.e-2)
    {
      snsumsn[iq] = qsumsn[iq] / nsumsn[iq];
    }
    else
    {
      snsumsn[iq] = 0.;
    }
    qsumspiral[iq] = qsumspiral[iq - 1] + qspiral[iq];
    nsumspiral[iq] = sqrt(nsumspiral[iq - 1] * nsumspiral[iq - 1] +
                          nspiral[iq] * nspiral[iq]);
    if (nsumspiral[iq] > 1.e-2)
    {
      snsumspiral[iq] = qsumspiral[iq] / nsumspiral[iq];
    }
    else
    {
      snsumspiral[iq] = 0.;
    }
  } // end loop on pixels

  //------------
  // compute some charges relevant for the cluster

  ChargeAroundSeed = TotalCharge - qqordered[0]; // JB 2010/07/23
  ChargeInCross = qspiral[1] + qspiral[3] + qspiral[5] + qspiral[7];
  ChargeInX = qspiral[2] + qspiral[4] + qspiral[6] + qspiral[8];
  for (short i = 9; i <= 24; i++)
    ChargeInCrown1 += qspiral[i];
  for (short i = 25; i <= 48; i++)
    ChargeInCrown2 += qspiral[i];

  //-- compute charges limited to 3x3 cluster
  // (Note that if NofPixelsInCluster<25, arrays are updated with 0s)
  // pixel with index 0 is the seed, in the center
  for (short i = 0; i < 25; i++)
  {
    UofPix5x5[i] = UofPixspiral[i];
    VofPix5x5[i] = VofPixspiral[i];
    Qof5x5[i] = qspiral[i];
    TotalCharge5x5 += Qof5x5[i];
  }

  //-- compute charges limited to 3x3 cluster
  // (Note that if NofPixelsInCluster<9, arrays are updated with 0s)
  // pixel with index 0 is the seed, in the center
  for (short i = 0; i < 9; i++)
  {
    UofPix3x3[i] = UofPixspiral[i];
    VofPix3x3[i] = VofPixspiral[i];
    Qof3x3[i] = qspiral[i];
    TotalCharge3x3 += Qof3x3[i];
  }

  //-- compute charges limited to 2x2 cluster
  // Find 2x2 cluster, there are 4 of them
  //  we keep the one with the largest charge:
  float sum[4] = {0., 0., 0., 0.};
  sum[0] = qspiral[1] + qspiral[2] + qspiral[3];
  sum[1] = qspiral[3] + qspiral[4] + qspiral[5];
  sum[2] = qspiral[5] + qspiral[6] + qspiral[7];
  sum[3] = qspiral[7] + qspiral[8] + qspiral[1];
  // find the 2x2 cluster with the largest portion of energy
  ExternalCharge2x2 = sum[0];
  IndexOfCluster2x2 = 0;
  for (int i = 1; i < 4; i++)
  {
    if (ExternalCharge2x2 < sum[i])
    {
      ExternalCharge2x2 = sum[i];
      IndexOfCluster2x2 = i;
    }
  }
  TotalCharge2x2 = ExternalCharge2x2 + qspiral[0];
  // Fill the pixel position of the selected 2x2 cluster
  Qof2x2[0] = qspiral[0];
  UofPix2x2[0] = UofPixspiral[0];
  VofPix2x2[0] = VofPixspiral[0];
  for (short i = 1; i < 4; i++)
  {
    Qof2x2[i] = qspiral[2 * IndexOfCluster2x2 + i];
    UofPix2x2[i] = UofPixspiral[2 * IndexOfCluster2x2 + i];
    VofPix2x2[i] = VofPixspiral[2 * IndexOfCluster2x2 + i];
  }
}

void MimosaAnalysis::SetClusterGeometricalType(int nPixelsInCluster, int lineSizeOfCluster, int columnSizeOfCluster)
{

  // Geometrical type of cluster
  //  see the definition of clusterTypeDef in MAnalysis.h
  // JB 2014/03/31
  // From 1 to 4 pixels

  const auto key = std::make_tuple(nPixelsInCluster, lineSizeOfCluster, columnSizeOfCluster);

  auto it = clusterMap.find(key);

  if (it != clusterMap.end())
  {
    // Found the matching ClusterType
    ClusterGeometricalType = it->second;
  }
  else
  {
    // Handle the case where no match is found, you can throw an error or handle default case
    ClusterGeometricalType = MAnalysis::ClusterType::more_than_six; // Default case
  }
  // rGeometricalTypeBeyond4 = it->second;
  // }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterCharges_fill(DAuthenticHit *thehit, int ievt)
{
  // Fill the different cluster charges
  //
  // Called for each matched clusters
  //
  // JB 2010/06/03 based on original MPro stuff
  // Modified: JB 2013/11/08 some new histograms

  if (MimoDebug)
    Info("ClusterCharges_fill", " Filling histos with charges for hit %d",
         thehit->Hhk);

  // -----------------------------
  // charge of individual pixels or cumulated
  // -----------------------------
  for (int i = 0; i < NofPixelsInCluster; i++)
  {
    hChargeInCluster->Fill(qonly[i]);
    hindivq[i]->Fill(qqordered[i]);
    hsnn[i]->Fill(snsnordered[i]);
    // cumulated charges
    hqcn[i]->Fill(q[i]);
    hqcngeom[i]->Fill(qsumspiral[i]);
  }

  // seed pixel
  hChargeInSeed->Fill(qonly[0]); // seed ?
  hRealTrackNoise->Fill(noise[0]);
  hRealTrackNoise_time->Fill(ievt, noise[0]);
  hSNReal->Fill(snsnordered[0]);

  // pixel with the second highest charge
  hS2N2ndRH->Fill(qqordered[1]);
  hSNNReal->Fill(snsnordered[1]);
  hind->Fill(IndexInCluster[1]);

  /*
   // 9 first pixels, order from TTree
     hChargeOrder1->Fill(qonly[0]);
     hChargeOrder2->Fill(qonly[1]);
     hChargeOrder3->Fill(qonly[2]);
     hChargeOrder4->Fill(qonly[3]);
     hChargeOrder5->Fill(qonly[4]);
     hChargeOrder6->Fill(qonly[5]);
     hChargeOrder7->Fill(qonly[6]);
     hChargeOrder8->Fill(qonly[7]);
     hChargeOrder9->Fill(qonly[8]);
  */

  hChargeCor_1_2->Fill(qqordered[1], qqordered[0]);
  hChargeCor_1_3->Fill(qqordered[2], qqordered[0]);
  hChargeCor_1_4->Fill(qqordered[3], qqordered[0]);
  hChargeCor_2_3->Fill(qqordered[2], qqordered[1]);
  if ((qqordered[1] != 0.0) && (qqordered[2] != 0.0) && (qqordered[2] != 0.0))
  {
    hChargeRap_1_over_2->Fill(qqordered[0] / qqordered[1]);
    hChargeRap_1_over_3->Fill(qqordered[0] / qqordered[2]);
    hChargeRap_1_over_4->Fill(qqordered[0] / qqordered[3]);
    hChargeRap_2_over_3->Fill(qqordered[1] / qqordered[2]);
  }

  if (fabs(snqordered[1]) > 0.)
    h_SNRratioL->Fill(snqordered[0] / snqordered[1]);

  // Restricted to 3x3 pixels cluster around seed
  for (short i = 0; i < 9; i++)
  {
    hQofPix3x3->Fill(Qof3x3[i]);
  }

  // Restricted to 2x2 pixels cluster around seed
  hIndex2x2->Fill(IndexOfCluster2x2);

  // to study "shadows":
  if (Qof3x3[4] + Qof3x3[5] > 0. && (Qof3x3[2] + Qof3x3[7]) > 0.)
  {
    hQ3x345->Fill(Qof3x3[4] - Qof3x3[5]);
    hQ3x327->Fill(Qof3x3[2] - Qof3x3[7]);
    hQ3x345r->Fill((Qof3x3[4] - Qof3x3[5]) / (Qof3x3[4] + Qof3x3[5]));
    hQ3x327r->Fill((Qof3x3[2] - Qof3x3[7]) / (Qof3x3[2] + Qof3x3[7]));
    hQ3x34->Fill(Qof3x3[4]);
    hQ3x35->Fill(Qof3x3[5]);
  }

  // -----------------------------
  // Cluster charge
  // -----------------------------

  hqc_c->Fill(thehit->Hqc); // cluster charge from TTree
  hqc->Fill(thehit->Hqc * calibration);
  hsnc1->Fill(thehit->HSNc1); // cluster S/N from TTree

  hChargeSum_4->Fill(qonly[0] + qonly[1] + qonly[2] + qonly[3]);
  hChargeInCross->Fill(ChargeInCross);
  hChargeInX->Fill(ChargeInX);
  hChargeInCrown1->Fill(ChargeInCrown1);
  hChargeInCrown2->Fill(ChargeInCrown2);

  // cluster with the largest S/N
  hoptimalsize->Fill(optimalSize);
  if (optimalSize > 0)
    hqcel->Fill(q[optimalSize - 1]);
  hsnc->Fill(snOptimal);

  hSNneighbour->Fill(
      thehit->HSNneighbour); // cluster charge without seed from TTree
  hSNseedvsSNneighbour->Fill(snsnordered[0],
                             thehit->HSNneighbour); // JB 2013/11/07
  hQseedvsQcluster->Fill(thehit->Hqc * calibration,
                         ChargeAroundSeed);                  // JB 2014/05/22
  hQseedvsQneighbour->Fill(qsnordered[0], ChargeAroundSeed); // JB 2013/11/08
  hSNseedvsQcluster->Fill(snsnordered[0],
                          thehit->Hqc * calibration); // JB 2014/05/22

  // -----------------------------
  // Pixel multiplicity
  // -----------------------------
  // # pixels in the cluster
  hnpix_c->Fill(thehit->HNNS); // Hsk); JB 2009/09/01
  // Counting how many 4x2 (colxlin) withdows are need to contain cluster.
  for (int isupcol = 0; isupcol < h2DWindows->GetXaxis()->GetNbins();
       isupcol++)
  {
    for (int isuplin = 0; isuplin < h2DWindows->GetYaxis()->GetNbins();
         isuplin++)
    {
      h2DWindows->SetBinContent(isupcol + 1, isuplin + 1, 0);
    }
  }
  // if(thehit->HNNS == 4) {
  for (int ipix = 0; ipix < thehit->HNNS; ipix++)
  {
    int col_pix = (thehit->HkM[ipix] % NofPixelInRaw) / 4;
    int lin_pix = (thehit->HkM[ipix] / NofPixelInRaw) / 2;
    h2DWindows->Fill(col_pix, lin_pix);
  }
  int NWindows = 0;
  for (int isupcol = 0; isupcol < h2DWindows->GetXaxis()->GetNbins();
       isupcol++)
  {
    for (int isuplin = 0; isuplin < h2DWindows->GetYaxis()->GetNbins();
         isuplin++)
    {
      if (h2DWindows->GetBinContent(isupcol + 1, isuplin + 1) > 0)
      {
        NWindows++;
      }
    }
  }
  hnWindows_c->Fill(NWindows);
  //}

  // -----------------------------
  // Correlation between different charges
  // -----------------------------
  //-------Signal / noise for different combinations of pixels: (0=seed) ordered

  hS2N2RH->Fill(snqordered[0], qqordered[1]);
  hseedQvsS2NGood->Fill(snqordered[0], qqordered[0]);
  hSN_vs_SNNReal->Fill(snsnordered[1], snsnordered[0]);

  float temp_pulsesum = 0;
  float temp_noisesqsum = 0;
  for (int ip = 0; ip < 2; ip++)
  {
    for (int jp = 0; jp < jpixmax; jp++)
    {
      temp_pulsesum = 0;
      temp_noisesqsum = 0;
      int counterpulse;
      for (counterpulse = ip; counterpulse < jp + 2; counterpulse++)
      {
        // oredered by charge:
        //		      temp_pulsesum += qqordered[counterpulse];
        //		      temp_noisesqsum += nqordered[counterpulse] *
        // nqordered[counterpulse];
        // ordered by S/N:
        temp_pulsesum += qsnordered[counterpulse];
        temp_noisesqsum += nsnordered[counterpulse] * nsnordered[counterpulse];
      }
      temp_noisesqsum = sqrt(temp_noisesqsum);
      if (temp_noisesqsum != 0)
      {
        if (ip == 0)
        {
          hsn_pix_0[jp]->Fill(temp_pulsesum / temp_noisesqsum);
          hsn_seed_vs_pix_0[jp]->Fill(temp_pulsesum / temp_noisesqsum,
                                      snsnordered[0]);
        }
        if (ip == 1)
        {
          hsn_pix_1[jp]->Fill(temp_pulsesum / temp_noisesqsum);
          hsn_seed_vs_pix_1[jp]->Fill(temp_pulsesum / temp_noisesqsum,
                                      snsnordered[0]);
        }
      }
    }
  }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterCharges_end()
{
  // Finalize histograms created and filled for cluster charges studies
  //
  // Called once for all
  //
  // JB 2010/07/22 based on cdritsa original stuff
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterPosition_init()
{

  // Initializer to be called once for all
  //
  // Created: JB 2014/01/10

  hitCounterPos = 0;
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterPosition_compute(DAuthenticHit *thehit,
                                             DPrecAlign *align)
{
  // compute cluster position
  //
  // Called for each clusters
  //
  // JB 2010/06/03 based on original MPro stuff
  // Modified: JB 2010/09/10, correct digital position from Ttree information
  // (thehit) Modified: JB 2013/05/01, convert position in tracker frame
  // (require alignment) Modified: JB 2013/01/10, increment hit counter

  if (MimoDebug)
    Info("ClusterPosition_compute",
         " Computing positions for hit %d in event, %d /total", thehit->Hhk,
         hitCounterPos);

  hitCounterPos++; // JB 2014/01/10

  // digital position
  hUdigital = thehit->Hsu;
  hVdigital = thehit->Hsv;

  // position as computed by clustering algo
  hU = thehit->Hu;
  hV = thehit->Hv;
  // hW = fSession->GetSetup()->GetPlanePar(thehit->Hpk).Position(3);
  hW = 0.; // not sure yet, JB 2010/09/08

  // Convert in tracker frame
  // JB 2013/05/01
  DR3 aPosition(hU, hV, hW);
  if (align != nullptr)
  {
    hX = align->TransformHitToTracker(aPosition)(0);
    hY = align->TransformHitToTracker(aPosition)(1);
    hZ = align->TransformHitToTracker(aPosition)(2);
  }
  else
  {
    hX = 0.;
    hY = 0.;
    hZ = 0.;
  }

  // Sometimes max S/N doesn't match Signal max pixel.
  // The following is just a correction for such events:
  if (hU - hUdigital > PixelSizeU / 2.)
  {
    hUdigital += PixelSizeU;
  }
  else if (hU - hUdigital < -PixelSizeU / 2.)
  {
    hUdigital -= PixelSizeU;
  }
  if (hV - hVdigital > PixelSizeV / 2.)
  {
    hVdigital += PixelSizeV;
  }
  else if (hV - hVdigital < -PixelSizeV / 2.)
  {
    hVdigital -= PixelSizeV;
  }

  // CoG position
  ClusterPosition_cog(thehit);

  // Eta position
  ClusterPosition_eta(thehit);

  // AHT position
  ClusterPosition_aht(thehit);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterPosition_eta(DAuthenticHit *thehit)
{
  // compute cluster position with different variations of eta algorithm
  //
  // Called for each clusters
  //
  // JB 2010/06/03 based on original MPro stuff
  // Modified JB 2010/08/27 to avoid crash when CorPar not ready

  if (MimoDebug)
    Info("ClusterPosition_eta",
         " Computing eta positions for hit %d (CorPar file staus=%d)",
         thehit->Hhk, CorStatus);

  // Quit if CorPar file not ready (status==2)
  // JB 2010/08/27
  if (CorStatus == 2)
    return;

  // ---------------------------
  // Init for new positions
  UofHitEta3 = 0.; // Eta method for 3 "strips"
  VofHitEta3 = 0.;
  UofHitEta2x2 = 0.; // Eta method for 2x2
  VofHitEta2x2 = 0.;
  UofHitEta2x2_new = 0.;
  VofHitEta2x2_new = 0.;
  UofHitEta5x5_new = 0.;
  VofHitEta5x5_new = 0.;
  UofHitEta2x2_newR = 0.;
  VofHitEta2x2_newR = 0.;
  UofHitEta5x5_newR = 0.;
  VofHitEta5x5_newR = 0.;

  // ---------------------------
  //== ----etaab debut UofHitEta2x2_new
  // cluster 4 neighbours
  float dist_eta2x2 = PixelSize * 1.1;
  const int Nofeta2x2_max = 25;
  float UofPixeta2x2[Nofeta2x2_max];
  float VofPixeta2x2[Nofeta2x2_max];
  float Qofeta2x2[Nofeta2x2_max];
  float TotalCharge_eta2x2 = 0.0;
  int Nofeta2x2 = 0;
  for (int i = 0; i < NofPixelsInCluster; i++)
  {
    float distanceeta2x2 =
        sqrt((UofPix[i] - UofPix[0]) * (UofPix[i] - UofPix[0]) +
             (VofPix[i] - VofPix[0]) * (VofPix[i] - VofPix[0]));
    if (distanceeta2x2 < dist_eta2x2)
    {
      UofPixeta2x2[Nofeta2x2] = UofPix[i];
      VofPixeta2x2[Nofeta2x2] = VofPix[i];
      Qofeta2x2[Nofeta2x2] = qonly[i];
      Nofeta2x2++;
    }
  }
  for (int i = 0; i < Nofeta2x2; i++)
  {
    UofHitEta2x2_new += Qofeta2x2[i] * UofPixeta2x2[i];
    VofHitEta2x2_new += Qofeta2x2[i] * VofPixeta2x2[i];
    TotalCharge_eta2x2 += Qofeta2x2[i];
  }
  if (TotalCharge_eta2x2 != 0.0)
  {
    UofHitEta2x2_new = UofHitEta2x2_new / TotalCharge_eta2x2;
    VofHitEta2x2_new = VofHitEta2x2_new / TotalCharge_eta2x2;
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_eta: Eta_2x2_new = %.1f, %.1f\n", UofHitEta2x2_new,
           VofHitEta2x2_new);

  // ---------------------------
  // Eta 5x5 NEW
  // cluster 24 neighbours
  float dist_eta5x5 = PixelSize * 5;
  // const int Nofeta5x5_max = 25;
  const int Nofeta5x5_max = MaxNofPixelsInCluster;
  float UofPixeta5x5[Nofeta5x5_max];
  float VofPixeta5x5[Nofeta5x5_max];
  float Qofeta5x5[Nofeta5x5_max];
  float TotalCharge_eta5x5 = 0.0;
  int Nofeta5x5 = 0;
  for (int i = 0; i < NofPixelsInCluster; i++)
  {
    float distanceeta5x5 =
        sqrt((UofPix[i] - UofPix[0]) * (UofPix[i] - UofPix[0]) +
             (VofPix[i] - VofPix[0]) * (VofPix[i] - VofPix[0]));
    if (distanceeta5x5 < dist_eta5x5)
    {
      UofPixeta5x5[Nofeta5x5] = UofPix[i];
      VofPixeta5x5[Nofeta5x5] = VofPix[i];
      Qofeta5x5[Nofeta5x5] = qonly[i];
      Nofeta5x5++;
    }
  }
  for (int i = 0; i < Nofeta5x5; i++)
  {
    if ((qonly[i] / noise[i]) > 1.0)
    {
      UofHitEta5x5_new += TMath::Power(Qofeta5x5[i], 1.0) * UofPixeta5x5[i];
      VofHitEta5x5_new += TMath::Power(Qofeta5x5[i], 1.0) * VofPixeta5x5[i];
      TotalCharge_eta5x5 += TMath::Power(Qofeta5x5[i], 1.0);
      // cout<<"** "<<q[i]/noise[i];
    }
  }
  // cout<<" --> "<<Nofeta5x5<<endl;
  if (TotalCharge_eta5x5 != 0.0)
  {
    UofHitEta5x5_new = UofHitEta5x5_new / TotalCharge_eta5x5;
    VofHitEta5x5_new = VofHitEta5x5_new / TotalCharge_eta5x5;
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_eta: Eta_5x5_new = %.1f, %.1f\n", UofHitEta5x5_new,
           VofHitEta5x5_new);

  // ---------------------------
  // 2x2 calculations:
  float Eta2x2U = 0.;
  float Eta2x2V = 0.;
  if (TotalCharge2x2 != 0)
  {
    switch (IndexOfCluster2x2)
    {
    case 0:
      Eta2x2U = (TotalCharge2x2 - ExternalCharge2x2 +
                 Qof3x3[IndexU2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      Eta2x2V = (TotalCharge2x2 - ExternalCharge2x2 +
                 Qof3x3[IndexV2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      break;
    case 1:
      Eta2x2U = (ExternalCharge2x2 - Qof3x3[IndexU2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      Eta2x2V = (TotalCharge2x2 - ExternalCharge2x2 +
                 Qof3x3[IndexV2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      break;
    case 2:
      Eta2x2U = (TotalCharge2x2 - ExternalCharge2x2 +
                 Qof3x3[IndexU2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      Eta2x2V = (ExternalCharge2x2 - Qof3x3[IndexV2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      break;
    case 3:
      Eta2x2U = (ExternalCharge2x2 - Qof3x3[IndexU2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      Eta2x2V = (ExternalCharge2x2 - Qof3x3[IndexV2x2[IndexOfCluster2x2]]) /
                TotalCharge2x2;
      break;
    }
  }
  else
  {
    if (thehit->HNNS > 0)
      cout << "segment WARNING TotalCharge2x2 = 0" << endl;
  }

  int iBin2x2(0);

  while (Eta2x2U > Edges2x2[iBin2x2] && iBin2x2 < NBins2x2 - 1)
  {
    iBin2x2++;
  }
  //  cout << "iBin2x2 = " << iBin2x2 << endl;
  UofHitEta2x2 = Contents2x2U[iBin2x2] * PixelSizeU - PixelSizeU / 2. +
                 UofPix2x2[IndexOfCluster2x2] / 2. + hUdigital;
  //  cout << "UEta2x2 = " << UofHitEta2x2 << "tu=" << tu << endl;
  iBin2x2 = 0;
  while (Eta2x2V > Edges2x2[iBin2x2] && iBin2x2 < NBins2x2 - 1)
  {
    iBin2x2++;
  }
  VofHitEta2x2 = Contents2x2V[iBin2x2] * PixelSizeV - PixelSizeV / 2. +
                 VofPix2x2[IndexOfCluster2x2] / 2. + hVdigital;
  if (MimoDebug > 1)
    printf("ClusterPosition_eta: Eta2x2 = %.1f, %.1f\n", UofHitEta2x2,
           VofHitEta2x2);

  // ---------------------------
  // calculate the coordinate by eta and profile:
  float xxx = UofHitCG - hUdigital;
  float yxx = VofHitCG - hVdigital;
  //--etaab debut xxx
  float xxx_eta2x2 = UofHitEta2x2_new - hUdigital;
  float yxx_eta2x2 = VofHitEta2x2_new - hVdigital;
  float xxx_eta5x5 = UofHitEta5x5_new - hUdigital;
  float yxx_eta5x5 = VofHitEta5x5_new - hVdigital;
  //----etaab fin
  int iBin(0);
  // cout << "xxx= " << xxx << endl;
  //		  while(xxx > hEta3Int->GetBinLowEdge(iBin++));
  while (xxx > Edges3x3[iBin] && iBin < NBins3x3 - 1)
  {
    // cout << xxx << hEta3Int->GetBinLowEdge(iBin) << endl;
    //  cout << "iBin = " << iBin << endl;
    iBin++;
  }
  //		  		  cout << "iBin = " << iBin << endl;
  UofHitEta3 = Contents3x3U[iBin] + hUdigital;
  //    cout << "U = " << UofHitEta3 << endl;
  iBin = 0;
  while (yxx > Edges3x3[iBin] && iBin < NBins3x3 - 1)
  {
    iBin++;
  }
  VofHitEta3 = Contents3x3V[iBin] + hVdigital;
  // JB, Don't get the use of the following lines ?!?
  UofHitEta3 = hUdigital;
  VofHitEta3 = hVdigital;

  for (int i = 0; i < 7; i++)
  {
    UofHitEta3 += FitParEta3U[i] * TMath::Power(xxx, i);
    VofHitEta3 += FitParEta3V[i] * TMath::Power(yxx, i);
    /*		//------------------ETASTUDIES
if((ievt==14404)||(ievt==10937)){
cout<<i <<"  -------------- UofHitEta3 VofHitEta3 "<<UofHitEta3<<"
"<<VofHitEta3<<endl; cout<<" hUdigital hVdigital "<<hUdigital<<"
"<<hVdigital<<endl; cout<<" UofHitCG VofHitCG "<<UofHitCG<<" "<<VofHitCG<<endl;
cout<<" xxx yxx "<<xxx <<" "<<yxx<<endl;
cout<<"FitParEta3U[i] FitParEta3V[i] "<<FitParEta3U[i] <<" "<<FitParEta3V[i]
<<endl; cout<<"TMath::Power(xxx,i) TMath::Power(yxx,i)"<< TMath::Power(xxx,i)
<<" "<<TMath::Power(yxx,i)<<endl; cout<<"  -------------- "<<endl;

}
//--------------------- */
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_eta: Eta3 = %.1f, %.1f\n", UofHitEta3, VofHitEta3);

  // ---------------------------
  //------------------------------------------------New ETA METHOD WITHOUT FIT
  // BIAS
  if (CorStatus != 2 && READnListe_CoG > 0)
  { // if corPar OK and EtaList filled
    UofHitEta3 = hUdigital;
    VofHitEta3 = hVdigital;
    //---U eta 3x3
    int etacounter = 0;
    int etastop = 1;
    while (etastop == 1)
    {
      if (READListe_CoGU[etacounter] > xxx)
      {
        etastop = 0;
      }
      else
      {
        // cout<<"etacounter etastop READnListe_CoG  READListe_CoGU[etacounter]
        // xxx "
        //   <<etacounter<<" "<<etastop <<" "<<READnListe_CoG <<"
        //   "<<READListe_CoGU[etacounter]<<" "<<xxx<<endl;
        etacounter++;
        if (etacounter >= READnListe_CoG)
        {
          etastop = 0;
        }
      }
    }

    UofHitEta3 += PixelSizeU * (float(etacounter) / float(READnListe_CoG)) -
                  PixelSizeU / 2.;
    //---V eta 3x3
    etacounter = 0;
    etastop = 1;
    while (etastop == 1)
    {
      if (READListe_CoGV[etacounter] > yxx)
      {
        etastop = 0;
      }
      else
      {
        // cout<<"etacounter etastop READnListe_CoG  READListe_CoGV[etacounter]
        // yxx "
        //   <<etacounter<<" "<<etastop <<" "<<READnListe_CoG <<"
        //   "<<READListe_CoGV[etacounter]<<" "<<yxx<<endl;
        etacounter++;
        if (etacounter >= READnListe_CoG)
        {
          etastop = 0;
        }
      }
    }
    VofHitEta3 += PixelSizeV * (float(etacounter) / float(READnListe_CoG)) -
                  PixelSizeV / 2.;

    if (MimoDebug > 1)
      printf("ClusterPosition_eta: Eta3(bis) = %.1f, %.1f\n", UofHitEta3,
             VofHitEta3);

    //---etaab debut Look eta list
    //---U eta 2x2
    etacounter = 0;
    etastop = 1;
    UofHitEta2x2_newR = hUdigital;
    while (etastop == 1)
    {
      if (READListe_CoGU_eta2x2[etacounter] > xxx_eta2x2)
      {
        etastop = 0;
      }
      else
      {
        etacounter++;
        if (etacounter >= READnListe_CoG_eta2x2)
        {
          etastop = 0;
        }
      }
    }
    UofHitEta2x2_newR +=
        PixelSizeU * (float(etacounter) / float(READnListe_CoG_eta2x2)) -
        PixelSizeU / 2.;
    //---V eta 2x2
    etacounter = 0;
    etastop = 1;
    VofHitEta2x2_newR = hVdigital;
    while (etastop == 1)
    {
      if (READListe_CoGV_eta2x2[etacounter] > yxx_eta2x2)
      {
        etastop = 0;
      }
      else
      {
        etacounter++;
        if (etacounter >= READnListe_CoG_eta2x2)
        {
          etastop = 0;
        }
      }
    }
    VofHitEta2x2_newR +=
        PixelSizeV * (float(etacounter) / float(READnListe_CoG_eta2x2)) -
        PixelSizeV / 2.;

    if (MimoDebug > 1)
      printf("ClusterPosition_eta: Eta2x2_newR = %.1f, %.1f\n",
             UofHitEta2x2_newR, VofHitEta2x2_newR);

    //---U eta 5x5
    etacounter = 0;
    etastop = 1;
    UofHitEta5x5_newR = hUdigital;
    while (etastop == 1)
    {
      if (READListe_CoGU_eta5x5[etacounter] > xxx_eta5x5)
      {
        etastop = 0;
      }
      else
      {
        etacounter++;
        if (etacounter >= READnListe_CoG_eta5x5)
        {
          etastop = 0;
        }
      }
    }
    UofHitEta5x5_newR +=
        PixelSizeU * (float(etacounter) / float(READnListe_CoG_eta5x5)) -
        PixelSizeU / 2.;

    //---V eta 5x5
    etacounter = 0;
    etastop = 1;
    VofHitEta5x5_newR = hVdigital;
    while (etastop == 1)
    {
      if (READListe_CoGV_eta5x5[etacounter] > yxx_eta5x5)
      {
        etastop = 0;
      }
      else
      {
        etacounter++;
        if (etacounter >= READnListe_CoG_eta5x5)
        {
          etastop = 0;
        }
      }
    }
    VofHitEta5x5_newR +=
        PixelSizeV * (float(etacounter) / float(READnListe_CoG_eta5x5)) -
        PixelSizeV / 2.;

    if (MimoDebug > 1)
      printf("ClusterPosition_eta: Eta2x2_newR = %.1f, %.1f\n",
             UofHitEta5x5_newR, VofHitEta5x5_newR);

  } // end if corPar OK and EtaList filled
  else
  {
    if (MimoDebug)
      printf(
          "ClusterPosition_eta: Cannot compute new Eta with list either "
          "because CorPar has wrong status %d or because list has %d entries\n",
          CorStatus, READnListe_CoG);
  }
  //---------------------------------END OF New ETA METHOD WITHOUT FIT BIAS

  // ---------------------------
  // Update of liste to compute a better Eta correction ???
  if (nListe_CoG < 7500)
  {
    Liste_CoGU.Set(nListe_CoG + 1);
    Liste_CoGV.Set(nListe_CoG + 1);
    Liste_CoGU.AddAt(UofHitCG - hUdigital, nListe_CoG);
    Liste_CoGV.AddAt(VofHitCG - hVdigital, nListe_CoG);
    nListe_CoG++; // NListe_Cog = is now the exact number of stored points.
    Liste_CoGU_eta2x2.Set(nListe_CoG_eta2x2 + 1);
    Liste_CoGV_eta2x2.Set(nListe_CoG_eta2x2 + 1);
    Liste_CoGU_eta2x2.AddAt(UofHitEta2x2_new - hUdigital, nListe_CoG_eta2x2);
    Liste_CoGV_eta2x2.AddAt(VofHitEta2x2_new - hVdigital, nListe_CoG_eta2x2);
    nListe_CoG_eta2x2++; // NListe_Cog_eta2x2 = is now the exact number of
                         // stored points.
    Liste_CoGU_eta5x5.Set(nListe_CoG_eta5x5 + 1);
    Liste_CoGV_eta5x5.Set(nListe_CoG_eta5x5 + 1);
    Liste_CoGU_eta5x5.AddAt(UofHitEta5x5_new - hUdigital, nListe_CoG_eta5x5);
    Liste_CoGV_eta5x5.AddAt(VofHitEta5x5_new - hVdigital, nListe_CoG_eta5x5);
    nListe_CoG_eta5x5++; // NListe_Cog_eta5x5 = is now the exact number of
                         // stored points.
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_eta: Liste for Eta corrections updated std %d, 2x2 "
           "%d, 5x5 %d\n",
           nListe_CoG, nListe_CoG_eta2x2, nListe_CoG_eta5x5);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterPosition_cog(DAuthenticHit *thehit)
{
  // compute cluster position with different variations of the center of gravity
  //
  // Called for each clusters
  //
  // JB 2010/06/04 based on original MPro stuff
  // Modified: SS 2011/10/26

  if (MimoDebug)
    Info("ClusterPosition_cog", " Computing cog positions for hit %d",
         thehit->Hhk);

  // Init for new positions
  UofHitCG = 0.;
  VofHitCG = 0.;
  UCG2x2 = 0.;
  VCG2x2 = 0.;
  UofHitCG3 = 0.;
  VofHitCG3 = 0.;
  UofHitCG5 = 0.;
  VofHitCG5 = 0.;
  UofHitCG = 0.;
  VofHitCG = 0.;
  UCGcorr = 0.;
  VCGcorr = 0.;

  // center of gravity with all pixels in the cluster
  if (TotalCharge > 0.)
  {
    for (int i = 0; i < NofPixelsInCluster; i++)
    {
      UofHitCG += qonly[i] * UofPix[i] / TotalCharge;
      VofHitCG += qonly[i] * VofPix[i] / TotalCharge;
      // printf(" pix %d, q=%.1f, pos=(%.1f,%.1f), totq=%.1f,
      // finalPos=(%.1f,%.1f)\n", i, qonly[i], UofPix[i], VofPix[i],
      // TotalCharge, UofHitCG5, VofHitCG5);
    }
    // printf("track pos=(%.1f,%.1f)\n", tu, tv);
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_cog: CG = %.1f, %.1f\n", UofHitCG, VofHitCG);

  // center of gravity with the sub-cluster of 2x2 pixels with the highest
  // charge
  if (TotalCharge2x2 > 0.)
  {
    for (int i = 0; i < 4; i++)
    { // SS 2011/10/26 correction to the i limit
      UCG2x2 += Qof2x2[i] * UofPix2x2[i] / TotalCharge2x2;
      VCG2x2 += Qof2x2[i] * VofPix2x2[i] / TotalCharge2x2;
    }
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_cog: CG2x2 = %.1f, %.1f\n", UCG2x2, VCG2x2);

  // center of gravity with a sub-cluster of 3x3 pixels
  if (TotalCharge3x3 > 0.)
  {
    for (int i = 0; i < 9; i++)
    {
      UofHitCG3 += Qof3x3[i] * UofPix3x3[i] / TotalCharge3x3;
      VofHitCG3 += Qof3x3[i] * VofPix3x3[i] / TotalCharge3x3;
    }
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_cog: CG3 = %.1f, %.1f\n", UofHitCG3, VofHitCG3);

  // center of gravity with a sub-cluster of 5x5 pixels
  if (TotalCharge5x5 > 0.)
  {
    for (int i = 0; i < 25; i++)
    { // SS 2011/10/26 correction to the i limit
      UofHitCG5 += Qof5x5[i] * UofPix5x5[i] / TotalCharge5x5;
      VofHitCG5 += Qof5x5[i] * VofPix5x5[i] / TotalCharge5x5;
    }
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_cog: CG5 = %.1f, %.1f\n", UofHitCG5, VofHitCG5);

  // center of gravity corrected from profile histograms (?)
  // int PolN = 6 ; // Polynom degree
  VCGcorr = hVdigital;
  UCGcorr = hUdigital;
  for (int i = 0; i < 7; i++)
  {
    UCGcorr += FitParamU[i] * TMath::Power(UofHitCG - hUdigital, i);
    VCGcorr += FitParamV[i] * TMath::Power(VofHitCG - hVdigital, i);
  }
  if (MimoDebug > 1)
    printf("ClusterPosition_cog: CGCorr = %.1f, %.1f\n", UCGcorr, VCGcorr);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterPosition_aht(DAuthenticHit *thehit)
{
  // compute cluster position with the analog head-tail algorithm
  //
  // Called for each clusters
  //
  // JB 2010/06/04 based on original MPro stuff

  if (MimoDebug)
    Info("ClusterPosition_aht", " Computing aht positions for hit %d",
         thehit->Hhk);

  Uaht = 0.;
  Vaht = 0.;

  // the following variables are the indexes of the pixels mostly on
  // the right, the left, the top or the bottom of the cluster
  int iLeft, iRight, iTop, iBottom;
  iLeft = iRight = iTop = iBottom = 0; // initialize to seed

  // this is the cut on SNR to consider the pixel or not
  double snrCutAht = 3.;
  for (int i = 1; i < NofPixelsInCluster; i++)
  {
    if (snonly[i] > snrCutAht)
    {
      if (IndexofPix[i] % NofPixelInRaw < iLeft)
        iLeft = i;
      if (IndexofPix[i] % NofPixelInRaw > iRight)
        iRight = i;
      if (IndexofPix[i] / NofPixelInRaw < iBottom)
        iBottom = i; // % replaced by /
      if (IndexofPix[i] / NofPixelInRaw > iTop)
        iTop = i; // % replaced by /
    }
  }

  // the aht position is the average between the two extreme pixels
  // corrected by the difference of charges
  Uaht = (UofPix[iLeft] + UofPix[iRight]) / 2. +
         (q[iLeft] - q[iRight]) / TotalCharge * PixelSizeU;
  Vaht = (VofPix[iBottom] + VofPix[iTop]) / 2. +
         (q[iBottom] - q[iTop]) / TotalCharge * PixelSizeV;

  if (MimoDebug > 2)
    printf(" pixLR %2d-%2d, pixBT  %2d-%2d, (u,v) = (%.1f, %.1f), du=%.1f, "
           "dv=%.1f\n",
           iLeft, iRight, iBottom, iTop, Uaht, Vaht, Uaht - tu, Vaht - tv);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterPosition_fill(DAuthenticHit *thehit)
{
  // fill the different selected cluster positions
  //
  // Called for each cluster
  //
  // JB 2010/07/21 based on original MPro stuff
  // Modified: JB 2013/05/01 add histo in tracker frame

  if (MimoDebug)
    Info("ClusterPosition_fill", " Filling the positions for hit %d",
         thehit->Hhk);

  // Optimal(?) position
  huv->Fill(hU, hV);
  huv_rate->Fill(hU, hV);
  hhu->Fill(hU);
  hhv->Fill(hV);

  // Position in tracker frame, JB 2013/05/01
  hxy->Fill(hX, hY);
  hhx->Fill(hX);
  hhy->Fill(hY);

  // Digital position
  h2dmatchedhits->Fill(hUdigital, hVdigital);
  hhuS->Fill(hUdigital);
  hhvS->Fill(hVdigital);

  // Eta position
  hEta2x2U->Fill(Eta2x2U);
  hEta2x2V->Fill(Eta2x2V);

  // 2D picture obtained from cumulated charges of hits
  // JB 2014/01/10
  for (int i = 0; i < NofPixelsInCluster; i++)
  { // loop on pixels in cluster
    h2DpictureMatched->Fill(IndexofPix[i] % NofPixelInRaw,
                            IndexofPix[i] / NofPixelInRaw, qqordered[i]);
  }

  // Digital compared to CoG
  hdCGDigU->Fill(UofHitCG - hUdigital);
  hdCGDigV->Fill(VofHitCG - hVdigital);
  hdCGDigUV->Fill(UofHitCG - hUdigital, VofHitCG - hVdigital);

  // Digital compared to Eta3
  hNorm->Fill(UofHitEta3 - hUdigital, VofHitEta3 - hVdigital);
  hEta3DigU->Fill(UofHitEta3 - hUdigital);
  hEta3DigV->Fill(VofHitEta3 - hVdigital);

  // CoG compared to Eta3
  hdCGEtaU->Fill(UofHitCG - UofHitEta3);
  hdCGEtaV->Fill(VofHitCG - VofHitEta3);

  // CoGcorrected compared to Eta3
  hDifCorU->Fill(UCGcorr - UofHitEta3);
  hDifCorV->Fill(VCGcorr - VofHitEta3);

  //-- ??

  etal[1]->Fill(ChargeAroundSeed);
  if ((thehit->Hq0) != 0.0)
  {
    etal1[1]->Fill((thehit->Hqc / thehit->Hq0) - 1.0);
  }
  else
  {
    if (thehit->HNNS > 0)
      cout << "segment WARNING thehit->Hq0 = 0" << endl;
  }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::ClusterPosition_end()
{
  // Finalize histograms created and filled
  //
  // Called once for all
  //
  // JB 2014/01/10

  cout << "******* hitCounterPos for cluster position study ******* "
       << hitCounterPos << endl;

  if (!hitCounterPos)
    return; // leave if no hits were selected

  double normfactor = (double)1. / hitCounterPos;

  h2DpictureMatched->SetEntries(hitCounterPos);
  h2DpictureMatched->Scale(normfactor);

  cout << "Done" << endl;
}

//_____________________________________________________________________________
//
void MimosaAnalysis::Inefficent_end()
{
  // A.B. April 2015
  // compute probability to associate the track to a random hit if the real hit
  // is inefficient.
  cout << "*******  MimosaAnalysis::void Inefficent_end ******* " << endl;
  // cout<<" N bins "<<hTrackTo2ndclosestClusterDistance->Integral()<<"
  // "<<hTrackTo2ndclosestClusterDistance->GetXaxis()->GetNbins()<<endl;
  //  Normalisation is the total number of matched hits:
  double Integral_norm = h2dmatchedhits->GetEntries();
  // hTrackTo2ndclosestClusterDistance->Integral();
  double Current_integral = 0.0;
  // Integral of hTrackTo2ndclosestClusterDistance:
  for (int ibin = 1;
       ibin < hTrackTo2ndclosestClusterDistance->GetXaxis()->GetNbins() + 1;
       ibin++)
  {
    Current_integral += hTrackTo2ndclosestClusterDistance->GetBinContent(ibin);
    if (Current_integral > 0.0)
    {
      hWrongAssociationProba->SetBinContent(ibin, Current_integral);
    }
  }
  // normalize probability:
  if (Integral_norm > 0)
  {
    hWrongAssociationProba->Scale(1. / Integral_norm);
  }

  cout << "Done" << endl;
}
//_____________________________________________________________________________
//
void MimosaAnalysis::TrackHitPosition_fill(DAuthenticHit *thehit,
                                           DTransparentPlane *atrack,
                                           int ievt)
{
  // fill histos comparing track and hit positions
  //
  // Called for each cluster matched with a track
  //
  // JB 2010/07/21 based on original MPro stuff
  // Modified CLM 2013/01/23 additional histos for charge PSF
  // Modified JB 2013/07/16 additional alignment check histo
  // Modified JB 2014/04/01, #hits/track added

  if (MimoDebug)
    Info("TrackHitPosition_fill",
         " Filling the positions for track matched to hit %d", thehit->Hhk);

  //------------------------
  //-- resolution study

  // Digital from TTree
  hAllHuvsAllTu1->Fill(tu - thehit->Hsu);
  hAllHvvsAllTv1->Fill(tv - thehit->Hsv);
  hAllHuvsAllTu2->Fill(tu - hUdigital); // never drawn !!!
  hAllHvvsAllTv2->Fill(tv - hVdigital);

  // CG from TTree, JB 2009/09/07
  huCGtu1->Fill(tu - thehit->HuCG);
  hvCGtv1->Fill(tv - thehit->HvCG);

  double u = (tu + 0.5 * NofPixelInRaw * PixelSizeU) / (2.0 * PixelSizeU);
  u = (u - int(u)) * 2.0 * PixelSizeU;
  double v = (tv + 0.5 * NofPixelInColumn * PixelSizeV) / (2.0 * PixelSizeV);
  v = (v - int(v)) * 2.0 * PixelSizeV;

  if (thehit->HNNS < huCGwidth_vs_Mult->GetXaxis()->GetNbins())
  {
    for (int imult = 0; imult < huCGwidth_vs_Mult->GetXaxis()->GetNbins() - 1;
         imult++)
    {
      if (thehit->HNNS == imult + 1)
      {
        huCGtu1_vs_Mult[imult]->Fill(tu - thehit->HuCG);
        hvCGtv1_vs_Mult[imult]->Fill(tv - thehit->HvCG);

        // double u = (tu + 0.5*NofPixelInRaw    * PixelSizeU)/(2.0*PixelSizeU);
        // u = (u - int(u))*2.0*PixelSizeU;
        // double v = (tv + 0.5*NofPixelInColumn * PixelSizeV)/(2.0*PixelSizeV);
        // v = (v - int(v))*2.0*PixelSizeV;
        huvCGtuv_vs_Mult[imult]->Fill(u, v);
        break;
      }
    }
  }
  else
  {
    huCGtu1_vs_Mult[huCGwidth_vs_Mult->GetXaxis()->GetNbins() - 1]->Fill(
        tu - thehit->HuCG);
    hvCGtv1_vs_Mult[huCGwidth_vs_Mult->GetXaxis()->GetNbins() - 1]->Fill(
        tv - thehit->HvCG);

    // double u = (tu + 0.5*NofPixelInRaw    * PixelSizeU)/(2.0*PixelSizeU);
    // u = (u - int(u))*2.0*PixelSizeU;
    // double v = (tv + 0.5*NofPixelInColumn * PixelSizeV)/(2.0*PixelSizeV);
    // v = (v - int(v))*2.0*PixelSizeV;
    huvCGtuv_vs_Mult[huCGwidth_vs_Mult->GetXaxis()->GetNbins() - 1]->Fill(u, v);
  }
  huvCGtuv->Fill(u, v);

  int idx_TrksInSensorPerEvt =
      hnTracksInGeomatrixVsTrackPerEvent->FindBin(NTracksPerEventInSensor) - 1;
  if (idx_TrksInSensorPerEvt >= 0 &&
      idx_TrksInSensorPerEvt <=
          hnTracksInGeomatrixVsTrackPerEvent->GetXaxis()->GetNbins() - 1)
  {
    huCGtu1_vs_TracksPerEvent[idx_TrksInSensorPerEvt]->Fill(tu - thehit->HuCG);
    hvCGtv1_vs_TracksPerEvent[idx_TrksInSensorPerEvt]->Fill(tv - thehit->HvCG);
  }

  // CG full cluster
  huCGtu2->Fill(tu - UofHitCG);
  hvCGtv2->Fill(tv - VofHitCG); // correction tv-UCG, 2011/10/06

  // CG 2x2 sub-cluster
  hCG2x2tu1->Fill(tu - UCG2x2);
  hCG2x2tv1->Fill(tv - VCG2x2);

  // CG 5x5 sub-cluster
  hCG5URes->Fill(tu - UofHitCG5);
  hCG5VRes->Fill(tv - VofHitCG5);

  // CG corrected(?)
  hTuHuCorr->Fill(tu - UCGcorr);
  hTvHvCorr->Fill(tv - VCGcorr);

  // Eta from TTree
  hEtaURes->Fill(tu - thehit->HuEta);
  hEtaVRes->Fill(tv - thehit->HvEta);

  // Eta 2x2
  hEta2x2tu1->Fill(tu - UofHitEta2x2);
  hEta2x2tv1->Fill(tv - VofHitEta2x2);

  // Eta 2x2_newR
  hEtaU_2x2Res->Fill(tu - UofHitEta2x2_newR);
  hEtaV_2x2Res->Fill(tv - VofHitEta2x2_newR);

  // Eta 5x5_newR
  hEtaU_5x5Res->Fill(tu - UofHitEta5x5_newR);
  hEtaV_5x5Res->Fill(tv - VofHitEta5x5_newR);

  // Eta3
  hEta3URes->Fill(tu - UofHitEta3);
  hEta3VRes->Fill(tv - VofHitEta3);

  // AHT
  hAHTURes->Fill(tu - Uaht);
  hAHTVRes->Fill(tv - Vaht);

  //------------------------
  //-- Track to hit distance study

  DuvCG->Fill(TrackToHitDistance);
  hdistchi2->Fill(chi2, TrackToHitDistance);
  float adist2 = sqrt((hUdigital - tu) * (hUdigital - tu) +
                      (hVdigital - tv) * (hVdigital - tv));
  htmp5->Fill(adist2);

  //------------------------
  //-- Track parameters

  hchi2_c->Fill(chi2);
  hNhitsPerTrack_all->Fill(atrack->TtHn); // JB 2014/04/01

  //------------------------
  //-- Alignment

  hAlignHuTuvsTv->Fill(tv, tu - thehit->Hu);
  hAlignHvTvvsTu->Fill(tu, tv - thehit->Hv);
  hAlignHuTuvsTu->Fill(tu, tu - thehit->Hu); // JB 2013/07/16
  hAlignHvTvvsTv->Fill(tv, tv - thehit->Hv);
  hAlignHuTuvsHv->Fill(thehit->Hv, tu - thehit->Hu);
  hAlignHvTvvsHu->Fill(thehit->Hu, tv - thehit->Hv);
  hAlignHuTu->Fill(tu - thehit->Hu);
  hAlignHvTv->Fill(tv - thehit->Hv);

  //------------------------
  //-- Correlation study

  vec->Fill(tu - hUdigital, tv - hVdigital);
  tudv->Fill(tu, tv - hVdigital);
  tvdu->Fill(tv, tu - hUdigital);
  hudv->Fill(hUdigital, tv - hVdigital);
  hvdu->Fill(hVdigital, tu - hUdigital);

  tuhu->Fill(tu, hUdigital);
  tuhu1->Fill(tu - hUdigital);
  tvhv->Fill(tv, hVdigital);
  tvhv1->Fill(tv - hVdigital);
  tuhv->Fill(tu, hVdigital);
  tvhu->Fill(tv, hUdigital);

  hEta2x2tu2->Fill(UofHitEta2x2 - hUdigital, tu - hUdigital);
  hEta2x2tv2->Fill(VofHitEta2x2 - hVdigital, tv - hVdigital);

  huCGtu->Fill(UofHitCG, tu);
  hvCGtv->Fill(VofHitCG, tv);
  huCG5tu->Fill(UofHitCG5, tu);
  hvCG5tv->Fill(VofHitCG5, tv);

  //   3x3 calculations:
  // Fill Profile of correlation plot:
  huCGtuInPix5->Fill(UofHitCG5 - hUdigital, tu - hUdigital);
  hvCGtvInPix5->Fill(VofHitCG5 - hVdigital, tv - hVdigital);
  huCGtuInPix->Fill(UofHitCG - hUdigital, tu - hUdigital);
  hvCGtvInPix->Fill(VofHitCG - hVdigital, tv - hVdigital);
  ProfUCG->Fill(UofHitCG - hUdigital, tu - hUdigital, 1);
  ProfVCG->Fill(VofHitCG - hVdigital, tv - hVdigital, 1);
  ProfACGn->Fill(VofHitCG - hVdigital, tv - hVdigital, 1);
  ProfACGn->Fill(UofHitCG - hUdigital, tu - hUdigital, 1);

  hUeta3TuInPix->Fill(UofHitEta3 - hUdigital, tu - hUdigital);
  hVeta3TvInPix->Fill(VofHitEta3 - hVdigital, tv - hVdigital);

  hChargeVsPosition->Fill(tu - hUdigital, tv - hVdigital,
                          Qof3x3[0] / TotalCharge3x3);
  hChargeVsDistance->Fill(TMath::Sqrt((tu - hUdigital) * (tu - hUdigital) +
                                      (tv - hVdigital) * (tv - hVdigital)),
                          Qof3x3[0]); // clm 2013/07/16
  hAllHitsInPixel->Fill(tu - hUdigital, tv - hVdigital, 1.);
  hUcorTuInPix->Fill(UCGcorr - hUdigital, tu - hUdigital);
  hVcorTvInPix->Fill(VCGcorr - hVdigital, tv - hVdigital);
  hTHCorr2->Fill(tu - UCGcorr, tv - VCGcorr);

  //---------------------------
  //-- Charge distribution:

  if (NofPixelsInCluster == MaxNofPixelsInCluster)
  {
    int NPixelsInLine = (int)sqrt((float)NofPixelsInCluster);   // a verif
    int NPixelsInColumn = (int)sqrt((float)NofPixelsInCluster); // a verif
    for (int iColumn = 0; iColumn < NPixelsInColumn; iColumn++)
    {
      for (int iLine = 0; iLine < NPixelsInLine; iLine++)
      {
        float uuu = MinUofPix + PixelSizeU / 2. + iLine * PixelSizeU;
        float vvv = MinVofPix + PixelSizeV / 2. + iColumn * PixelSizeV;
        float Charge1 = 0.;
        float Charge2 = 0.;
        float Charge3 = 0.;
        float Charge4 = 0.;
        for (int i = 0; i < NofPixelsInCluster; i++)
        {
          if (UofPix[i] < uuu && VofPix[i] < vvv)
          {
            Charge1 += qonly[i];
          }
          if (UofPix[i] > uuu && VofPix[i] < vvv)
          {
            Charge2 += qonly[i];
          }
          if (UofPix[i] < uuu && VofPix[i] > vvv)
          {
            Charge3 += qonly[i];
          }
          if (UofPix[i] > uuu && VofPix[i] > vvv)
          {
            Charge4 += qonly[i];
          }
        }
        if (TotalCharge > 0.)
        {
          hChargeIntegral1->Fill(uuu - tu, vvv - tv, Charge1 / TotalCharge);
          hChargeNorm1->Fill(uuu - tu, vvv - tv);
          hChargeIntegral2->Fill(-uuu + tu, vvv - tv, Charge2 / TotalCharge);
          hChargeNorm2->Fill(-uuu + tu, vvv - tv);
          hChargeIntegral3->Fill(uuu - tu, -vvv + tv, Charge3 / TotalCharge);
          hChargeNorm3->Fill(uuu - tu, -vvv + tv);
          hChargeIntegral4->Fill(-uuu + tu, -vvv + tv, Charge4 / TotalCharge);
          hChargeNorm4->Fill(-uuu + tu, -vvv + tv);
        }
      }
    }
  } // end of if(NofPixelsInCluster == MaxNofPixelsInCluster)

  //------------------------
  //-- Pixel response homogeneity study

  hHOM_tu_tv_modulo->Fill(tu - hUdigital, tv - hVdigital);
  hHOM_ResU_tu->Fill(tu - hUdigital, tu - UofHitEta3);
  hHOM_ResV_tv->Fill(tv - hVdigital, tv - VofHitEta3);
  ProfHOM_ResU_tu->Fill(tu - hUdigital, tu - UofHitEta3, 1);
  ProfHOM_ResV_tv->Fill(tv - hVdigital, tv - VofHitEta3, 1);
  float diodedist = sqrt((tu - hUdigital) * (tu - hUdigital) +
                         (tv - hVdigital) * (tv - hVdigital));
  hHOM_Charge_diodedist->Fill(diodedist, qonly[0]);
  ProfHOM_Charge_diodedist->Fill(diodedist, qonly[0], 1);

  // Added JB 2010/03/11, then clm 2013/01/23
  double GOOD_The_Total_Charge = 0;
  for (int iPix = 0; iPix < NofPixelsInCluster; iPix++)
  {
    GOOD_The_Total_Charge += qonly[iPix];
  }

  hNpixInClu->Fill(NofPixelsInCluster);
  for (int iPix = 0; iPix < NofPixelsInCluster; iPix++)
  { // NofPixelsInCluster
    // printf("  pix %d, q=%.0f, lin %.1f, col %.1f, tu-tv=%.0f-%.0f,
    // hu-hv=%.0f-%.0f, pitch=%.1f, diodedist(u-v)=%.1f-%.1f\n", iPix,
    // qonly[iPix], LineInCluster[iPix], ColumnInCluster[iPix], tu, tv,
    // hUdigital, hVdigital, PixelSizeU,
    // tu-hUdigital-ColumnInCluster[iPix]*PixelSizeV,
    // tv-hVdigital-LineInCluster[iPix]*PixelSizeV);

    hQpixInClu->Fill(qonly[iPix]);

    // ProfHOM_Charge_diodedist_alg->Fill(tu-hUdigital-ColumnInCluster[iPix]*PixelSizeU,qonly[iPix],1);
    ProfHOM_Charge_diodedist_alg->Fill(
        sqrt(pow(tu - hUdigital - ColumnInCluster[iPix] * PixelSizeU, 2) +
             pow(tv - hVdigital - LineInCluster[iPix] * PixelSizeV, 2)),
        qonly[iPix], 1);
    ProfHOM_Charge_diodedist_alg_u->Fill(tu - hUdigital -
                                             ColumnInCluster[iPix] * PixelSizeU,
                                         qonly[iPix], 1); // clm
    ProfHOM_Charge_diodedist_alg_v->Fill(tv - hVdigital -
                                             LineInCluster[iPix] * PixelSizeV,
                                         qonly[iPix], 1); // clm

    hHOM_Charge_diodedist_alg->Fill(
        sqrt(pow(tu - hUdigital - ColumnInCluster[iPix] * PixelSizeU, 2) +
             pow(tv - hVdigital - LineInCluster[iPix] * PixelSizeV, 2)),
        qonly[iPix]);
    hHOM_Charge_diodedist_alg_u->Fill(tu - hUdigital -
                                          ColumnInCluster[iPix] * PixelSizeU,
                                      qonly[iPix]); // clm
    hHOM_Charge_diodedist_alg_v->Fill(
        tv - hVdigital - LineInCluster[iPix] * PixelSizeV, qonly[iPix]); // clm

    // clm good
    //  ProfhGOODCharge_Charge_DiodePosition        ->Fill(tu-hUdigital-
    //  ColumnInCluster[iPix]*PixelSizeU,tv-hVdigital-LineInCluster[iPix]*PixelSizeV,qonly[iPix]/GOOD_The_Total_Charge);
    //  //clm clm test
    ProfhGOODCharge_Charge_DiodePosition->Fill(
        tu - UofPix[iPix], tv - VofPix[iPix],
        qonly[iPix] / GOOD_The_Total_Charge); // clm

    if (qonly[0] < 300)
      ProfhGOODCharge_Charge_DiodePositionSeedQLT300->Fill(
          tu - hUdigital - ColumnInCluster[iPix] * PixelSizeU,
          tv - hVdigital - LineInCluster[iPix] * PixelSizeV,
          qonly[iPix] / GOOD_The_Total_Charge); // clm

    if (qonly[0] > 2000)
      ProfhGOODCharge_Charge_DiodePositionSeedQGT2000->Fill(
          tu - hUdigital - ColumnInCluster[iPix] * PixelSizeU,
          tv - hVdigital - LineInCluster[iPix] * PixelSizeV,
          qonly[iPix] / GOOD_The_Total_Charge); // clm

    // clm calculate if the row/col is odd or even
    int myLineInClusterSeed = (thehit->HkM[0]) / NofPixelInRaw;
    int myColumnInClusterSeed = (thehit->HkM[0]) % NofPixelInRaw;
    int mylin = LineInCluster[iPix] +
                myLineInClusterSeed;                           //  as in ClusterCharge_compute method
    int mycol = ColumnInCluster[iPix] + myColumnInClusterSeed; //

    // if (myColumnInClusterSeed%2==0 &&  mycol%2==0 ) printf("same even col as
    // seed, V dist: PixelSizeV %lf from VofPix[iPix] %lf
    // \n",(hVdigital-LineInCluster[iPix]*PixelSizeV),(hVdigital-VofPix[iPix]));

    // if (myColumnInClusterSeed%2==0 &&  mycol%2!=0 ) printf("diff even col as
    // seed, V dist: PixelSizeV %lf from VofPix[iPix] %lf
    // \n",(hVdigital-LineInCluster[iPix]*PixelSizeV)/LineInCluster[iPix],(hVdigital-VofPix[iPix])/LineInCluster[iPix]);

    hDistVSeedOtherOldCalc->Fill(hVdigital - LineInCluster[iPix] * PixelSizeV,
                                 LineInCluster[iPix]);
    hDistVSeedOtherNewCalc->Fill(hVdigital - VofPix[iPix], LineInCluster[iPix]);

    // Corrected pixel position
    //  all pixels
    if (mycol % 2 == 0 && mylin % 2 == 0)
      ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow->Fill(
          tu - UofPix[iPix], tv - VofPix[iPix],
          qonly[iPix] / GOOD_The_Total_Charge); // clm
    if (mycol % 2 == 0 && mylin % 2 != 0)
      ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow->Fill(
          tu - UofPix[iPix], tv - VofPix[iPix],
          qonly[iPix] / GOOD_The_Total_Charge); // clm
    if (mycol % 2 != 0 && mylin % 2 == 0)
      ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow->Fill(
          tu - UofPix[iPix], tv - VofPix[iPix],
          qonly[iPix] / GOOD_The_Total_Charge); // clm
    if (mycol % 2 != 0 && mylin % 2 != 0)
      ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow->Fill(
          tu - UofPix[iPix], tv - VofPix[iPix],
          qonly[iPix] / GOOD_The_Total_Charge); // clm
    // seed
    if (iPix == 0)
    {
      if (mycol % 2 == 0 && mylin % 2 == 0)
        ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow_seed->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 == 0 && mylin % 2 != 0)
        ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow_seed->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 != 0 && mylin % 2 == 0)
        ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow_seed->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 != 0 && mylin % 2 != 0)
        ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow_seed->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
    }
    // 1st crown
    if (IsPixelIn1stCrown(ColumnInCluster[iPix], LineInCluster[iPix]) ==
        kTRUE)
    {
      if (mycol % 2 == 0 && mylin % 2 == 0)
        ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow_1stcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 == 0 && mylin % 2 != 0)
        ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow_1stcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 != 0 && mylin % 2 == 0)
        ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow_1stcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 != 0 && mylin % 2 != 0)
        ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow_1stcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
    }
    // 2nd crown
    if (IsPixelIn2ndCrown(ColumnInCluster[iPix], LineInCluster[iPix]) ==
        kTRUE)
    {
      if (mycol % 2 == 0 && mylin % 2 == 0)
        ProfhGOODCharge_Charge_DiodePosition_evencol_evenrow_2ndcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 == 0 && mylin % 2 != 0)
        ProfhGOODCharge_Charge_DiodePosition_evencol_oddrow_2ndcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 != 0 && mylin % 2 == 0)
        ProfhGOODCharge_Charge_DiodePosition_oddcol_evenrow_2ndcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
      if (mycol % 2 != 0 && mylin % 2 != 0)
        ProfhGOODCharge_Charge_DiodePosition_oddcol_oddrow_2ndcrown->Fill(
            tu - UofPix[iPix], tv - VofPix[iPix],
            qonly[iPix] / GOOD_The_Total_Charge); // clm
    }

    h2dCharge_Charge_DiodePosition_Track->Fill(
        tu - hUdigital - ColumnInCluster[iPix] * PixelSizeU,
        tv - hVdigital - LineInCluster[iPix] * PixelSizeV);

    h2dCharge_Charge_DiodePosition_CluSize->Fill(
        tu - hUdigital - ColumnInCluster[iPix] * PixelSizeU,
        tv - hVdigital - LineInCluster[iPix] * PixelSizeV, NofPixelsInCluster);

    ProfhGOODCharge_Charge_DiodePositionSimpDist->Fill(
        tu - hUdigital, tv - hVdigital,
        qonly[iPix] / GOOD_The_Total_Charge); // clm

    hHOM_Charge_diodedist3D->Fill(
        tu - hUdigital - ColumnInCluster[iPix] * PixelSizeU,
        tv - hVdigital - LineInCluster[iPix] * PixelSizeV, qonly[iPix]); // clm

    //                if( tv-hVdigital>PixelSizeV ) {
    //                  printf("  far away pix %d, q=%.0f, lin %.0f, col %.0f,
    //                  tu-tv=%.0f-%.0f, hu-hv=%.0f-%.0f, pitch=%.1f,
    //                  diodedist(u-v)=%.1f-%.1f\n", iPix, qonly[iPix],
    //                  LineInCluster[iPix], ColumnInCluster[iPix], tu, tv,
    //                  hUdigital, hVdigital, PixelSizeV,
    //                  tu-hUdigital-ColumnInCluster[iPix]*PixelSizeU,
    //                  tv-hVdigital-LineInCluster[iPix]*PixelSizeV);
    //                }
  }

  if (NofPixelsInCluster > 1)
  {
    hHOM_Charge2_diodedist->Fill(diodedist, q[1]);
    ProfHOM_Charge2_diodedist->Fill(diodedist, q[1], 1);
    hHOM_Charge2_diodedist3D->Fill(tu - hUdigital, tv - hVdigital, q[1]); // clm
  }
  if (NofPixelsInCluster > 3)
  {
    hHOM_Charge4_diodedist->Fill(diodedist, q[3]);
    ProfHOM_Charge4_diodedist->Fill(diodedist, q[3], 1);
    hHOM_Charge4_diodedist3D->Fill(tu - hUdigital, tv - hVdigital, q[3]); // clm
  }
  if (NofPixelsInCluster > 8)
  {
    hHOM_Charge9_diodedist->Fill(diodedist, q[8]);
    ProfHOM_Charge9_diodedist->Fill(diodedist, q[8], 1);
    hHOM_Charge9_diodedist3D->Fill(tu - hUdigital, tv - hVdigital, q[8]); // clm
  }
  if (NofPixelsInCluster > 24)
  {
    hHOM_Charge25_diodedist->Fill(diodedist, q[24]);
    ProfHOM_Charge25_diodedist->Fill(diodedist, q[24], 1);
    hHOM_Charge25_diodedist3D->Fill(tu - hUdigital, tv - hVdigital,
                                    q[24]); // clm
  }
  hHOM_SNseed_diodedist->Fill(diodedist, snsnordered[0]);
  ProfHOM_SNseed_diodedist->Fill(diodedist, snsnordered[0], 1);

  hHOM_Noise_diodedist->Fill(diodedist, thehit->Hsn * calibration);

  if (diodedist < (0.1 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_00_10->Fill(qonly[0]);
  }
  else if (diodedist < (0.2 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_10_20->Fill(qonly[0]);
  }
  else if (diodedist < (0.3 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_20_30->Fill(qonly[0]);
  }
  else if (diodedist < (0.4 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_30_40->Fill(qonly[0]);
  }
  else if (diodedist < (0.5 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_40_50->Fill(qonly[0]);
  }
  else if (diodedist < (0.6 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_50_60->Fill(qonly[0]);
  }
  else if (diodedist < (0.7 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_60_70->Fill(qonly[0]);
  }
  else if (diodedist < (0.8 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_70_80->Fill(qonly[0]);
  }
  else if (diodedist < (0.9 * PixelSize / sqrt(2.0)))
  {
    hHOM_Charge_diodedist_80_90->Fill(qonly[0]);
  }
  else
  {
    hHOM_Charge_diodedist_90_inf->Fill(qonly[0]);
  }

  hHOM_modUCG_modtu->Fill(tu - hUdigital, UofHitCG - hUdigital);
  hHOM_modVCG_modtv->Fill(tv - hVdigital, VofHitCG - hVdigital);
  hHOM_modUeta3_modtu->Fill(tu - hUdigital, UofHitEta3 - hUdigital);
  hHOM_modVeta3_modtv->Fill(tv - hVdigital, VofHitEta3 - hVdigital);
  hHOM_modUeta3_realtu->Fill(tu, UofHitEta3 - hUdigital);
  hHOM_modVeta3_realtv->Fill(tv, VofHitEta3 - hVdigital);
  hHOM_modUCG_realtu->Fill(tu, UofHitCG - hUdigital);
  hHOM_modVCG_realtv->Fill(tv, VofHitCG - hVdigital);
  hHOM_modUeta3_Eta3U->Fill(UofHitEta3, UofHitEta3 - hUdigital);
  hHOM_modVeta3_Eta3V->Fill(VofHitEta3, VofHitEta3 - hVdigital);
  hHOM_modUeta3_modVeta3->Fill(VofHitEta3 - hVdigital, UofHitEta3 - hUdigital);
  hHOM_modUCG_modVCG->Fill(VofHitCG - hVdigital, UofHitCG - hUdigital);
  hHOM_modUeta3_modUCG->Fill(UofHitCG - hUdigital, UofHitEta3 - hUdigital);
  hHOM_modVeta3_modVCG->Fill(VofHitCG - hVdigital, VofHitEta3 - hVdigital);
  hHOM_DU_Nevent->Fill(ievt, tu - UofHitEta3);
  hHOM_DV_Nevent->Fill(ievt, tv - VofHitEta3);
  hHOM_modtu_Nevent->Fill(ievt, tu - hUdigital);
  hHOM_modtv_Nevent->Fill(ievt, tv - hVdigital);
  hHOM_modUCG_Nevent->Fill(ievt, UofHitCG - hUdigital);
  hHOM_modVCG_Nevent->Fill(ievt, VofHitCG - hVdigital);

  //---etaab	 debut	store eta list
  hHOM_modUeta2x2_modVeta2x2->Fill(VofHitEta2x2_newR - hVdigital,
                                   UofHitEta2x2_newR - hUdigital);
  hHOM_modUCG2x2_modVCG2x2->Fill(VofHitEta2x2_new - hVdigital,
                                 UofHitEta2x2_new - hUdigital);
  hHOM_modUeta5x5_modVeta5x5->Fill(VofHitEta5x5_newR - hVdigital,
                                   UofHitEta5x5_newR - hUdigital);
  hHOM_modUCG5x5_modVCG5x5->Fill(VofHitEta5x5_new - hVdigital,
                                 UofHitEta5x5_new - hUdigital);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::FakeRate_end(int rateNormalisation)
{
  // Finalize histograms created and filled for fake rate studies
  //
  // Called once for all
  //
  // JB 2010/07/26
  // Modified JB 2011/07/01 Corrected protection for 0 normalisation

  if (MimoDebug)
    Info("FakeRate_end", " Finalizing fakerate histograms for %d events",
         rateNormalisation);

  // exit if the normalisation is wrong (=0)
  if (rateNormalisation == 0)
  {
    return;
  } // condition corrected, JB 2011/07/01

  // Get the # hits and rate from the 2D hit distribution
  for (int i = 1; i < NofPixelInColumn + 1; i++)
  {
    for (int j = 1; j < NofPixelInRaw + 1; j++)
    {
      hNhitperpixel->Fill(h2DgoodSeedPixel->GetBinContent(j, i));
      hNhitRateperpixel->Fill(h2DgoodSeedPixel->GetBinContent(j, i) /
                              float(rateNormalisation));
    }
  }

  // Fill the cumulative distribution
  for (short ib = 1; ib <= hNhitRateperpixel->GetNbinsX(); ib++)
  {
    int idx = hNhitRateperpixel->GetNbinsX() - ib - 1;
    // for( short jb=1; jb<=ib; jb++) {
    for (short jb = 1; jb <= idx; jb++)
    {
      hPixelsPerFakeRate->AddBinContent(ib,
                                        hNhitRateperpixel->GetBinContent(jb));
    }
  }

  hNhitperpixel->Scale(1.0 / (NofPixelInColumn * NofPixelInRaw));
  hNhitRateperpixel->Scale(1.0 / (NofPixelInColumn * NofPixelInRaw));
  hPixelsPerFakeRate->Scale(1.0 / (NofPixelInColumn * NofPixelInRaw));
  h2DgoodSeedPixel->Scale(1.0 / rateNormalisation);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::TrackParameters_compute(DTransparentPlane *atrack,
                                             DPrecAlign *align)
{
  // compute track parameters
  //
  // Called for each good track
  //
  // JB 2010/07/23 based on original MPro stuff
  // Modified: JB 2010/08/29, use the precision alignement global pointer
  // Modified: JB 2011/10/30, angles computed from slope

  if (MimoDebug)
    Info("TrackParameters_compute", " Computing parameters for track %d",
         atrack->Ttk);

  chi2 = atrack->Tchi2;
  tdu = atrack->Tdu;
  tdv = atrack->Tdv;
  tk1 = atrack->Tk1;
  tx = atrack->Tx;
  ty = atrack->Ty;
  tz = atrack->Tz;
  tdx = atrack->Tdx;
  tdy = atrack->Tdy;
  trackAngleXZ = (TMath::ATan(tdx)) * 180 / TMath::Pi(); // JB 2011/10/30
  trackAngleYZ =
      (TMath::ATan(tdy)) * 180 / TMath::Pi(); // moved here from MiniVector

  // -- correct track position in MIMOSA
  // from new alignment
  DataPoints myhit;
  DR3 TrackPos(0.0, 0.0, 0.0);
  myhit.Set(0, 0, 0, 0, tx, ty, tz, tdx, tdy);
  // myhit.Print();
  align->CalculateIntersection(&myhit);
  // cout << " track in plane (Teles frame) " << align->GetTrackPosition()(0) <<
  // ", " << align->GetTrackPosition()(1) << ", " <<
  // align->GetTrackPosition()(2) << "." << endl;
  TrackPos = align->TransformTrackToPlane();
  // cout << " track in plane (Plane frame) " << TrackPos(0) << ", " <<
  // TrackPos(1) << ", " << TrackPos(2) << "." << endl;
  tu = TrackPos(0);
  tv = TrackPos(1);

  // Temporary fix, bypass CorPar File alignement JB 2009/08/25
  // tu = atrack->Tu;
  // tv = atrack->Tv;

  if (MimoDebug > 1)
    cout << " Old track pos = " << atrack->Tu << ";" << atrack->Tv
         << " new track pos = " << tu << ";" << tv << endl;
}

//_____________________________________________________________________________
//
void MimosaAnalysis::TrackParameters_allFill(DTransparentPlane *atrack,
                                             int ievt)
{
  // fill histos for all tracks
  //
  // JB 2010/07/23 based on original MPro stuff
  // Modified JB 2014/04/01, #hits/track added

  if (MimoDebug)
    Info("TrackParameters_allfill", " Filling parameters for (all) track %d",
         atrack->Ttk);

  hchi2->Fill(atrack->Tchi2);
  hAllTvTu->Fill(atrack->Tv, atrack->Tu);
  hAllTu->Fill(atrack->Tu);
  hAllTv->Fill(atrack->Tv);
  hNhitsPerTrack_all->Fill(atrack->TtHn); // JB 2014/04/01
}

//_____________________________________________________________________________
//
void MimosaAnalysis::TrackParameters_goodFill(DTransparentPlane *atrack,
                                              int ievt)
{
  // fill histos for good tracks
  //
  // JB 2010/07/23 based on original MPro stuff
  // Modified JB 2011/10/30, angles added
  // Modified JB 2014/04/01, #hits/track added

  if (MimoDebug)
    Info("TrackParameters_goodfill", " Filling parameters for good track %d",
         atrack->Ttk);

  hGoodChi2TvTu->Fill(tu, tv);
  hGoodChi2Tu->Fill(tu);
  hGoodChi2Tv->Fill(tv);
  hGoodChi2Tx->Fill(tx);
  hGoodChi2Ty->Fill(ty);

  hGoodChi2AngleXZ->Fill(trackAngleXZ); // JB 2011/10/30
  hGoodChi2AngleYZ->Fill(trackAngleYZ);

  hNhitsPerTrack_good->Fill(atrack->TtHn); // JB 2014/04/01
}

//_____________________________________________________________________________
//
void MimosaAnalysis::TrackParameters_end()
{

  if (MimoDebug)
  {
    cout << "******* Normalization for cluster charges in track-hit study: "
         << hChargeNorm1->Integral() << ": " << hChargeNorm2->Integral() << ": "
         << hChargeNorm3->Integral() << ": " << hChargeNorm4->Integral()
         << endl;
    cout << "******* Integral for cluster charges in track-hit study: "
         << hChargeIntegral1->Integral() << ": " << hChargeIntegral2->Integral()
         << ": " << hChargeIntegral3->Integral() << ": "
         << hChargeIntegral4->Integral() << endl;
  }

  if (hChargeNorm1->Integral() > 0.)
    hChargeIntegral1->Divide(hChargeNorm1);
  if (hChargeNorm2->Integral() > 0.)
    hChargeIntegral2->Divide(hChargeNorm2);
  if (hChargeNorm3->Integral() > 0.)
    hChargeIntegral3->Divide(hChargeNorm3);
  if (hChargeNorm4->Integral() > 0.)
    hChargeIntegral4->Divide(hChargeNorm4);
  if (hChargeIntegral1->Integral() > 0.)
    hChargeIntegral->Add(hChargeIntegral1);
  if (hChargeIntegral2->Integral() > 0.)
    hChargeIntegral->Add(hChargeIntegral2);
  if (hChargeIntegral3->Integral() > 0.)
    hChargeIntegral->Add(hChargeIntegral3);
  if (hChargeIntegral4->Integral() > 0.)
    hChargeIntegral->Add(hChargeIntegral4);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::MiniVector_compute()
{

  // Compute the parameters (position and slope) of the minivector
  // hit and track positions are expected to have been computed already
  //
  // NCS 2010/01/20
  // Modified: JB 2010/08/30

  if (MimoDebug)
    Info("MiniVector_compute", "Computing the minivectors information.");

  // --- compute minivector parameters
  // taken in the middle of the two planes
  vectorSlopeXZ =
      (hitPosXYZ[1](0) - hitPosXYZ[0](0)) / (hitPosXYZ[1](2) - hitPosXYZ[0](2));
  vectorSlopeYZ =
      (hitPosXYZ[1](1) - hitPosXYZ[0](1)) / (hitPosXYZ[1](2) - hitPosXYZ[0](2));

  vectorAngleXZ = (TMath::ATan(vectorSlopeXZ)) * 180 / TMath::Pi();
  vectorAngleYZ = (TMath::ATan(vectorSlopeYZ)) * 180 / TMath::Pi();

  double x =
      vectorSlopeXZ * (hitPosXYZ[1](2) + hitPosXYZ[0](2)) / 2. +
      (hitPosXYZ[0](0) * hitPosXYZ[1](2) - hitPosXYZ[1](0) * hitPosXYZ[0](2)) /
          (hitPosXYZ[1](2) - hitPosXYZ[0](2));
  double y =
      vectorSlopeYZ * (hitPosXYZ[1](2) + hitPosXYZ[0](2)) / 2. +
      (hitPosXYZ[0](1) * hitPosXYZ[1](2) - hitPosXYZ[1](1) * hitPosXYZ[0](2)) /
          (hitPosXYZ[1](2) - hitPosXYZ[0](2));
  double z = (hitPosXYZ[1](2) + hitPosXYZ[0](2)) / 2.;
  vectorPosXYZ.SetValue(x, y, z);

  // --- compute track extrapolation in the middle point
  trackMeanPosXYZ.SetValue(tx + tdx * vectorPosXYZ(2),
                           ty + tdy * vectorPosXYZ(2), vectorPosXYZ(2));

  if (MimoDebug > 1)
  {
    for (short iplane = 0; iplane < 2; iplane++)
    { // loop on planes
      cout << "MiniVector_compute: position for plane: " << iplane << endl;
      // cout << "  plane pos: " << << ", " << << ", " << << " and tilt " << <<
      // ", " << << ", " << << endl;
      cout << "  hit   UVW: " << hitPosUVW[iplane](0) << ", "
           << hitPosUVW[iplane](1) << ", " << hitPosUVW[iplane](2) << endl;
      cout << "  track UVW: " << trackPosUVW[iplane](0) << ", "
           << trackPosUVW[iplane](1) << ", " << trackPosUVW[iplane](2) << endl;
      cout << "  hit   XYZ: " << hitPosXYZ[iplane](0) << ", "
           << hitPosXYZ[iplane](1) << ", " << hitPosXYZ[iplane](2) << endl;
      cout << "  track XYZ: " << trackPosXYZ[iplane](0) << ", "
           << trackPosXYZ[iplane](1) << ", " << trackPosXYZ[iplane](2) << endl;
    } // end loop on planes
    cout << "Minivector Positions at Z=" << vectorPosXYZ(2) << ":" << endl;
    cout << "  minivector XY: " << vectorPosXYZ(0) << ", " << vectorPosXYZ(1)
         << endl;
    cout << "  track      XY: " << trackMeanPosXYZ(0) << ", "
         << trackMeanPosXYZ(1) << endl;
    cout << "Minivector Slope in XZ and YZ:" << endl;
    cout << "  minivector : " << vectorSlopeXZ << ", " << vectorSlopeYZ << endl;
    cout << "  track      : " << tdx << ", " << tdy << endl;
    cout << "Minivector Angle in XZ and YZ (deg):" << endl;
    cout << "  minivector : " << vectorAngleXZ << ", " << vectorAngleYZ << endl;
    cout << "  track      : " << trackAngleXZ << ", " << trackAngleYZ << endl;
  }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::MiniVector_fill(DAuthenticHit *hit1, DAuthenticHit *hit2)
{
  // fill histos comparing track and hit positions
  //
  // Called for each cluster matched with a track
  //
  // JB 2010/07/21 based on original MPro stuff
  // Modified: JB 2013/11/08 HnS was replaced by HNNS due to new def in TTree

  if (MimoDebug)
    Info("MiniVector_fill", " Filling the minivectors information.");

  // ---- Filling histo ---- //
  hDiffPosX->Fill(trackMeanPosXYZ(0) - vectorPosXYZ(0));
  hDiffPosY->Fill(trackMeanPosXYZ(1) - vectorPosXYZ(1));
  hDiffAngleX->Fill(trackAngleXZ - vectorAngleXZ);
  hDiffAngleY->Fill(trackAngleYZ - vectorAngleYZ);

  hxtxPL3->Fill(trackPosXYZ[0](0) - hitPosXYZ[0](0));
  hytyPL3->Fill(trackPosXYZ[0](1) - hitPosXYZ[0](1));
  hxtxPL4->Fill(trackPosXYZ[1](0) - hitPosXYZ[1](0));
  hytyPL4->Fill(trackPosXYZ[1](1) - hitPosXYZ[1](1));

  hutuPL3->Fill(trackPosUVW[0](0) - hitPosUVW[0](0));
  hvtvPL3->Fill(trackPosUVW[0](1) - hitPosUVW[0](1));
  hutuPL4->Fill(trackPosUVW[1](0) - hitPosUVW[1](0));
  hvtvPL4->Fill(trackPosUVW[1](1) - hitPosUVW[1](1));

  hdiffydiffx->Fill(hitPosXYZ[1](1) - hitPosXYZ[0](1),
                    hitPosXYZ[1](0) - hitPosXYZ[0](0)); // JB 2011/11/01
  if (hit1->HNNS == 1)
  {
    if (hit2->HNNS == 1)
    {
      hdiffydiffx11->Fill(hitPosXYZ[1](1) - hitPosXYZ[0](1),
                          hitPosXYZ[1](0) - hitPosXYZ[0](0));
      hDiffAngleX11->Fill(trackAngleXZ - vectorAngleXZ);
      hDiffAngleY11->Fill(trackAngleYZ - vectorAngleYZ);
    }
    else if (hit2->HNNS == 2)
    {
      hdiffydiffx12->Fill(hitPosXYZ[1](1) - hitPosXYZ[0](1),
                          hitPosXYZ[1](0) - hitPosXYZ[0](0));
      hDiffAngleX12->Fill(trackAngleXZ - vectorAngleXZ);
      hDiffAngleY12->Fill(trackAngleYZ - vectorAngleYZ);
    }
  }
  else if (hit1->HNNS == 2)
  {
    if (hit2->HNNS == 1)
    {
      hdiffydiffx21->Fill(hitPosXYZ[1](1) - hitPosXYZ[0](1),
                          hitPosXYZ[1](0) - hitPosXYZ[0](0));
      hDiffAngleX21->Fill(trackAngleXZ - vectorAngleXZ);
      hDiffAngleY21->Fill(trackAngleYZ - vectorAngleYZ);
    }
    else if (hit2->HNNS == 2)
    {
      hdiffydiffx22->Fill(hitPosXYZ[1](1) - hitPosXYZ[0](1),
                          hitPosXYZ[1](0) - hitPosXYZ[0](0));
      hDiffAngleX22->Fill(trackAngleXZ - vectorAngleXZ);
      hDiffAngleY22->Fill(trackAngleYZ - vectorAngleYZ);
    }
  }
  else if (hit2->HNNS > 2)
  {
    hdiffydiffxg1g1->Fill(hitPosXYZ[1](1) - hitPosXYZ[0](1),
                          hitPosXYZ[1](0) - hitPosXYZ[0](0));
    hDiffAngleXg1g1->Fill(trackAngleXZ - vectorAngleXZ);
    hDiffAngleYg1g1->Fill(trackAngleYZ - vectorAngleYZ);
  }

  // ---- Filling TTree ---- // NCS 250110
  /*
  evt     = ievt;
  Hx3     = hitPosXYZ[0](0);
  Hy3     = hitPosXYZ[0](1);
  Hz3     = hitPosXYZ[0](2);
  Hx4     = hitPosXYZ[1](0);
  Hy4     = hitPosXYZ[1](1);
  Hz4     = hitPosXYZ[1](2);
  Hxm     = vectorPosXYZ(0);
  Hym     = vectorPosXYZ(1);
  Hzm     = vectorPosXYZ(2);
  Tx3     = trackPosXYZ[0](0);
  Ty3     = trackPosXYZ[0](1);
  Tz3     = hitPosXYZ[0](2);
  Tx4     = trackPosXYZ[1](0);
  Ty4     = trackPosXYZ[1](1);
  Tz4     = hitPosXYZ[1](2);
  Txm     = trackMeanPosXYZ(0);
  Tym     = trackMeanPosXYZ(1);
  Tzm     = trackMeanPosXYZ(2);
  HSlopeX = vectorSlopeXZ;
  HSlopeY = vectorSlopeYZ;
  TSlopeX = tdx;
  TSlopeY = tdy;
  HAngleX = vectorAngleXZ;
  HAngleY = vectorAngleYZ;
  TAngleX = trackAngleXZ;
  TAngleY = trackAngleYZ;
  tree_plume->Fill();
  if( MimoDebug) Info("MiniVector_compute","TTree filled.");
  */
}

//_____________________________________________________________________________
//
void MimosaAnalysis::Efficiency_init()
{
  // Initialize counters for efficiency computation
  //
  // Called once before the loop on events
  //
  // JB 2011/1/04 based on original MPro(AB) stuff

  MimosaEfficiency = 0.;
  MimosaEfficiency_ERR = 0.;

  // MimosaEfficiency_Corr     = 0.;
  // MimosaEfficiency_Corr_ERR = 0.;

  // overall conuters
  NtrkInMimo = 0;
  NofClMatchTrack = 0;
  // RfakeAndNFrame_counter = 0.0;

  // counters per bunch of events
  temp_maxarray = 1000;
  temp_Efficiency_array = new float[temp_maxarray];
  temp_NofClMatchTrack = new float[temp_maxarray];
  temp_NtrkInMimo = new float[temp_maxarray];
  NeventRangeForEfficiency = 1000; // size of bunch
  ievt_array = 0;
  for (int i = 0; i < temp_maxarray; i++)
  {
    temp_Efficiency_array[i] = 0.;
    temp_NofClMatchTrack[i] = 0.;
    temp_NtrkInMimo[i] = 0.;
  }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::Efficiency_end(int eventsRead)
{
  // Finalize efficiency computation
  //
  // Called once after the loop on events
  //
  // JB 2011/1/04 based on original MPro(AB) stuff
  // Modified JB 2013/09/12 protection on array size temp_Efficiency_array

  if (MimoDebug)
    Info("Efficiency_end", " Finalizing efficiency computation for %d events",
         eventsRead);

  // -- overall efficiency

  if (NtrkInMimo != 0)
  {
    MimosaEfficiency = 1. * NofClMatchTrack / NtrkInMimo;
    MimosaEfficiency_ERR =
        sqrt(MimosaEfficiency * (1.0 - MimosaEfficiency) / NtrkInMimo);
  }

  // Special case when efficiency is close to 100 %
  // the uncertainty goes like 1/(N) and not anymore follows the binomial
  // variance using the estimated efficiency
  //   if(MimosaEfficiency>=0.99){
  //     MimosaEfficiency_ERR = 1. / NtrkInMimo;
  //   }

  //---ab   running efficiency

  int const numpoints = eventsRead / NeventRangeForEfficiency;

  float *xnum = new float[numpoints];
  float *xerr = new float[numpoints];
  float *yerr = new float[numpoints];
  float *ynumerr = new float[numpoints];
  float *ytraerr = new float[numpoints];

  if (MimoDebug)
    cout << " temp_Efficiency_array (step = " << NeventRangeForEfficiency
         << " ) " << endl;
  for (int ie = 0;
       ie < (eventsRead / NeventRangeForEfficiency) && ie < temp_maxarray;
       ie++)
  {
    if (temp_NofClMatchTrack[ie] > 0)
    {
      temp_Efficiency_array[ie] =
          temp_NofClMatchTrack[ie] / temp_NtrkInMimo[ie];
      yerr[ie] = sqrt(temp_Efficiency_array[ie] *
                      (1.0 - temp_Efficiency_array[ie]) / temp_NtrkInMimo[ie]);
      ynumerr[ie] =
          sqrt(temp_Efficiency_array[ie] * (1.0 - temp_Efficiency_array[ie]) *
               temp_NtrkInMimo[ie]);
      ytraerr[ie] = 0.;
    }
    else
    {
      temp_Efficiency_array[ie] = 2.0;
      yerr[ie] = 0;
      ynumerr[ie] = 0;
      ytraerr[ie] = 0.;
    }
    xnum[ie] = ie;
    xerr[ie] = 0.0;
    if (MimoDebug)
      cout << temp_Efficiency_array[ie] << " ";
  }
  if (MimoDebug)
  {
    cout << endl;
    cout << " temp_NofClMatchTrack =" << endl;
    for (int ie = 0; ie < (eventsRead / NeventRangeForEfficiency); ie++)
    {
      cout << temp_NofClMatchTrack[ie] << " ";
    }
    cout << endl;
    cout << " temp_NtrkInMimo =" << endl;
    for (int ie = 0; ie < (eventsRead / NeventRangeForEfficiency); ie++)
    {
      cout << temp_NtrkInMimo[ie] << " ";
    }
  }

  greff = new TGraphErrors(numpoints, xnum, temp_Efficiency_array, xerr, yerr);
  greff->SetName("geffevt");
  grnum =
      new TGraphErrors(numpoints, xnum, temp_NofClMatchTrack, xerr, ynumerr);
  grnum->SetName("gmatchedtrackevt");
  grevt = new TGraphErrors(numpoints, xnum, temp_NtrkInMimo, xerr, ynumerr);
  grevt->SetName("gmimotrackevt");

  Char_t efficiency_title[100];
  sprintf(efficiency_title, "Efficieny. Run %d Plane %d, range = %d evts",
          RunNumber, ThePlaneNumber, NeventRangeForEfficiency);
  greff->SetTitle(efficiency_title);
  greff->SetMarkerColor(4);
  greff->SetMarkerStyle(21);
  // greff->Draw("ALP");
  cout << endl
       << endl;

  effimap->Add(goodtracks); // MG 2011/07/08
  effimap->Divide(TrkInMimo);

  effinpixel->Add(HitInPixel); // JB 2021/06/11
  effinpixel->Add(TrkInPixel);
}

//_____________________________________________________________________________
//
void MimosaAnalysis::GetAnalysisGoal()
{

  // If the setup parameter is not void,
  // use it to switch on only some flags controling the analysis
  //
  // JB 2014/01/11

  // By default only tracks and clusters are analysed
  fIfReferenceTrack = kTRUE;
  fIfCalibration = kFALSE;
  fIfLaserSpot = kFALSE;
  fIfMiniVector = kFALSE;
  fIfVertex = kFALSE;
  fIfFake = kFALSE;
  fIfImaging = kFALSE;

  Char_t *goal = fSession->GetSetup()->GetAnalysisPar().AnalysisGoal;

  // only cluster are analysed
  if (strstr(goal, "cluster") || strstr(goal, "Cluster"))
  {
    fIfReferenceTrack = kFALSE;
    fIfCalibration = kFALSE;
    fIfLaserSpot = kFALSE;
    fIfMiniVector = kFALSE;
    fIfVertex = kFALSE;
    fIfFake = kFALSE;
    fIfImaging = kFALSE;
    cout
        << "Analysis goal is cluster -> switching off all track related histos."
        << endl;
  }
  // laser spot are analysed
  else if (strstr(goal, "calib") || strstr(goal, "Calib"))
  {
    fIfReferenceTrack = kFALSE;
    fIfCalibration = kTRUE;
    fIfLaserSpot = kFALSE;
    fIfMiniVector = kFALSE;
    fIfVertex = kFALSE;
    fIfFake = kFALSE;
    fIfImaging = kFALSE;
    cout << "Analysis goal is calibration -> switching off all track related "
            "histos."
         << endl;
  }
  // laser spot are analysed
  else if (strstr(goal, "laser") || strstr(goal, "Laser"))
  {
    fIfReferenceTrack = kFALSE;
    fIfCalibration = kFALSE;
    fIfLaserSpot = kTRUE;
    fIfMiniVector = kFALSE;
    fIfVertex = kFALSE;
    fIfFake = kFALSE;
    fIfImaging = kFALSE;
    cout << "Analysis goal is laser -> switching off all track related histos."
         << endl;
  }
  // minivectors are analysed, need also tracks
  else if (strstr(goal, "vector") || strstr(goal, "Vector"))
  {
    fIfReferenceTrack = kTRUE;
    fIfCalibration = kFALSE;
    fIfLaserSpot = kFALSE;
    fIfMiniVector = kTRUE;
    fIfVertex = kFALSE;
    fIfFake = kFALSE;
    fIfImaging = kFALSE;
    cout << "Analysis goal is minivector -> switching on related histos."
         << endl;
  }
  // vertices are analysed, need also tracks
  else if (strstr(goal, "vertex") || strstr(goal, "Vertex") ||
           strstr(goal, "vertice") || strstr(goal, "Vertice"))
  {
    fIfReferenceTrack = kTRUE;
    fIfCalibration = kFALSE;
    fIfLaserSpot = kFALSE;
    fIfMiniVector = kFALSE;
    fIfVertex = kTRUE;
    fIfFake = kFALSE;
    fIfImaging = kFALSE;
    cout << "Analysis goal is vertex -> switching on related histos." << endl;
  }
  // fake hit analysis
  else if (strstr(goal, "fake") || strstr(goal, "Fake"))
  {
    fIfReferenceTrack = kTRUE;
    fIfCalibration = kFALSE;
    fIfLaserSpot = kFALSE;
    fIfMiniVector = kFALSE;
    fIfVertex = kFALSE;
    fIfFake = kTRUE;
    fIfImaging = kFALSE;
    cout << "Analysis goal is fake -> switching on related histos." << endl;
  }
  // imaging analysis
  else if (strstr(goal, "imaging") || strstr(goal, "Imaging"))
  {
    fIfReferenceTrack = kFALSE;
    fIfCalibration = kFALSE;
    fIfLaserSpot = kFALSE;
    fIfMiniVector = kFALSE;
    fIfVertex = kFALSE;
    fIfFake = kFALSE;
    fIfImaging = kTRUE;
    cout << "Analysis goal is imaging -> switching on related histos." << endl;
  }
  else if (fIfReferenceTrack == kTRUE)
  {
    cout << "Analysis goal needs reference tracks -> switching on related "
            "histos."
         << endl;
  }
}

//_____________________________________________________________________________
//
void MimosaAnalysis::GetParameters()
{
  // Get the parameters of the MIMOSA sensor from the config file.
  // These parameters are taken in the "Parameter for Analysis" section
  // and depends on each submatrix.
  //
  // JB 2010/07/23 based on original MPro stuff
  // Modified JB 2013/07/17 new Matrixtype, UserFlag parameters introduced
  // Modified JB 2013/08/21,22 now load cut values
  // Modified JB 2013/09/12 setting of MinNofPixelsInCluster
  // Modified JB 2013/11/08 setting of new cuts on Q_seed and Q_neighbour
  // Modified JB 2014/01/21 addition of new cut Q_cluster
  // Modified JB 2015/01/28 check submatrix index is in proper range

  cout << "ThesubmatrixNumber = " << ThesubmatrixNumber << ", "
       << fSession->GetSetup()->GetAnalysisPar().Submatrices << endl;

  // --------------
  // Check the submatrix does exist
  if (ThesubmatrixNumber < 0 ||
      fSession->GetSetup()->GetAnalysisPar().Submatrices <=
          ThesubmatrixNumber)
  {
    Error("MAnalysis:GetParameters",
          "The submatrix specified is not in the existing range [%d, %d] ! --> "
          "STOP.",
          0, fSession->GetSetup()->GetAnalysisPar().Submatrices);
  }

  // --------------
  // parameters of the submatrix

  PixelSizeU =
      fSession->GetSetup()->GetAnalysisPar().PixelSizeU[ThesubmatrixNumber];
  PixelSizeV =
      fSession->GetSetup()->GetAnalysisPar().PixelSizeV[ThesubmatrixNumber];
  PixelSize = PixelSizeU;
  NofPixelInRaw =
      fSession->GetSetup()->GetAnalysisPar().PixelsInRaw[ThesubmatrixNumber];
  NofPixelInColumn =
      fSession->GetSetup()
          ->GetAnalysisPar()
          .PixelsInColumn[ThesubmatrixNumber];  // bug corrected, JB 2009/09/15
  NofPixels = NofPixelInRaw * NofPixelInColumn; // JB 2010/07/23

  FirstPixelShiftV = (NofPixelInColumn - 1) * PixelSizeV / 2;
  FirstPixelShiftU = (NofPixelInRaw - 1) * PixelSizeU / 2;

  MimosaSizeV = NofPixelInColumn / 2. * PixelSizeV;
  MimosaSizeU = NofPixelInRaw / 2. * PixelSizeU;

  calibration =
      fSession->GetSetup()->GetAnalysisPar().Calibration[ThesubmatrixNumber];
  NoiseScope =
      fSession->GetSetup()->GetAnalysisPar().NoiseScope[ThesubmatrixNumber];

  Matrixtype = fSession->GetSetup()
                   ->GetAnalysisPar()
                   .Matrixtype[ThesubmatrixNumber]; // JB 2013/07/17

  // User flag, JB 2013/07/17
  UserFlag = fSession->GetSetup()->GetAnalysisPar().UserFlag;

  // if(MimoDebug) {
  cout << " DUT - submatrix " << ThesubmatrixNumber << " has:" << endl;
  cout << "     # pixels: " << NofPixelInRaw << " in raw " << NofPixelInColumn
       << " in column" << endl
       << ", total " << NofPixels << endl;
  cout << "     pixel pitch " << PixelSizeU << " x " << PixelSizeV << endl;
  cout << "     matrix type (mapping) " << Matrixtype << endl;
  cout << "     first pixel position: U=" << FirstPixelShiftU
       << " V=" << FirstPixelShiftV << endl;
  cout << "     total size/2: U=" << MimosaSizeU << " V=" << MimosaSizeV
       << endl;
  cout << "     calibration " << calibration << endl;
  cout << "     noise scope " << NoiseScope << endl;
  cout << " UserFlag is " << UserFlag << endl;
  //}

  // --------------
  // cut values

  // on hits
  CUT_MaxNbOfHits = fSession->GetSetup()
                        ->GetAnalysisPar()
                        .MaxNbOfHits; // take 50 to clean picture //for Mimosa
                                      // 5, 400 hits crashes.
  CUT_MinNbOfHits = fSession->GetSetup()
                        ->GetAnalysisPar()
                        .MinNbOfHits; // for tests. put 0 by defaults.
  MaxNofPixelsInCluster = fSession->GetSetup()
                              ->GetAnalysisPar()
                              .MaxNofPixelsInCluster[ThesubmatrixNumber];
  MinNofPixelsInCluster = fSession->GetSetup()
                              ->GetAnalysisPar()
                              .MinNofPixelsInCluster[ThesubmatrixNumber];
  CUT_MinSeedIndex = fSession->GetSetup()
                         ->GetAnalysisPar()
                         .MinSeedIndex[ThesubmatrixNumber]; // JB 2013/08/21
  CUT_MaxSeedIndex = fSession->GetSetup()
                         ->GetAnalysisPar()
                         .MaxSeedIndex[ThesubmatrixNumber]; // JB 2013/08/21
  if (CUT_MinSeedIndex > CUT_MaxSeedIndex)
  {
    Warning("MimosaAnalysis::GetParameters",
            "Seed index range inconsistent, min=%d> max=%d --> both set to 0\n",
            CUT_MinSeedIndex, CUT_MaxSeedIndex);
    CUT_MinSeedIndex = CUT_MaxSeedIndex = 0;
  }
  CUT_MinSeedCol = fSession->GetSetup()
                       ->GetAnalysisPar()
                       .MinSeedCol[ThesubmatrixNumber]; // JB 2013/08/21
  CUT_MaxSeedCol = fSession->GetSetup()
                       ->GetAnalysisPar()
                       .MaxSeedCol[ThesubmatrixNumber]; // JB 2013/08/21
  if (CUT_MinSeedCol > CUT_MaxSeedCol)
  {
    Warning("MimosaAnalysis::GetParameters",
            "Seed col range inconsistent, min=%d> max=%d --> both set to 0\n",
            CUT_MinSeedCol, CUT_MaxSeedCol);
    CUT_MinSeedCol = CUT_MaxSeedCol = 0;
  }
  CUT_MinSeedRow = fSession->GetSetup()
                       ->GetAnalysisPar()
                       .MinSeedRow[ThesubmatrixNumber]; // JB 2013/08/21
  CUT_MaxSeedRow = fSession->GetSetup()
                       ->GetAnalysisPar()
                       .MaxSeedRow[ThesubmatrixNumber]; // JB 2013/08/21
  if (CUT_MinSeedRow > CUT_MaxSeedRow)
  {
    Warning("MimosaAnalysis::GetParameters",
            "Seed row range inconsistent, min=%d> max=%d --> both set to 0\n",
            CUT_MinSeedRow, CUT_MaxSeedRow);
    CUT_MinSeedRow = CUT_MaxSeedRow = 0;
  }
  CUT_Q_seed = fSession->GetSetup()
                   ->GetAnalysisPar()
                   .MinSeedCharge[ThesubmatrixNumber]; // JB 2013/11/08
  CUT_Q_cluster = fSession->GetSetup()
                      ->GetAnalysisPar()
                      .MinClusterCharge[ThesubmatrixNumber]; // JB 2014/01/21
  CUT_MinQ_neighbour =
      fSession->GetSetup()
          ->GetAnalysisPar()
          .MinNeighbourCharge[ThesubmatrixNumber]; // JB 2013/11/08

  // on tracks
  MinHitsPerTrack =
      fSession->GetSetup()->GetAnalysisPar().MinHitsPerTrack; // JB 2013/06/22
  TrackChi2Limit = fSession->GetSetup()->GetAnalysisPar().TrackChi2Limit;
  // the two next cut values are now read from configuration file, JB 2013/06/21
  MaxNbOfTracksInGeom = fSession->GetSetup()->GetAnalysisPar().MaxTracksExGeom;
  GeoMatrixForTrackCut = fSession->GetSetup()->GetAnalysisPar().ExGeomatrix;
  exgeomUmin = fSession->GetSetup()
                   ->GetAnalysisPar()
                   .Umin[ThesubmatrixNumber][GeoMatrixForTrackCut];
  exgeomUmax = fSession->GetSetup()
                   ->GetAnalysisPar()
                   .Umax[ThesubmatrixNumber][GeoMatrixForTrackCut];
  exgeomVmin = fSession->GetSetup()
                   ->GetAnalysisPar()
                   .Vmin[ThesubmatrixNumber][GeoMatrixForTrackCut];
  exgeomVmax = fSession->GetSetup()
                   ->GetAnalysisPar()
                   .Vmax[ThesubmatrixNumber][GeoMatrixForTrackCut];

  // on geometry
  geomUmin = fSession->GetSetup()->GetAnalysisPar().Umin[ThesubmatrixNumber]
                                                        [Thegeomatrix];
  geomUmax = fSession->GetSetup()->GetAnalysisPar().Umax[ThesubmatrixNumber]
                                                        [Thegeomatrix];
  geomVmin = fSession->GetSetup()->GetAnalysisPar().Vmin[ThesubmatrixNumber]
                                                        [Thegeomatrix];
  geomVmax = fSession->GetSetup()->GetAnalysisPar().Vmax[ThesubmatrixNumber]
                                                        [Thegeomatrix];

  cout << "----------- CUTS:" << endl;
  cout << "  # hits (min, max): " << CUT_MinNbOfHits << " - " << CUT_MaxNbOfHits
       << endl;
  cout << "  min charge on seed: " << CUT_Q_seed
       << ", on neighbour: " << CUT_MinQ_neighbour
       << ", on cluster: " << CUT_Q_cluster << endl;
  cout << "  min - max # pixels in a hit: " << MinNofPixelsInCluster << " - "
       << MaxNofPixelsInCluster << endl;
  cout << "  max track chi2 " << TrackChi2Limit << endl;
  cout << "  geomatrix is " << Thegeomatrix << ", limits in U: " << geomUmin
       << ", " << geomUmax << ", limits in V: " << geomVmin << ", " << geomVmax
       << endl;

  Info("MimosaParameters", "Parameters read");
}

//______________________________________________________________________________
//

void MimosaAnalysis::ComputePixelPosition(int col, int lin, float &u,
                                          float &v)
{

  double uu, vv, ww;

  fTool.ComputeStripPosition(Matrixtype, col, lin, uu, vv, ww, NofPixelInRaw,
                             NofPixelInColumn, (double)PixelSizeU,
                             (double)PixelSizeV, 1.);

  u = uu;
  v = vv;

  /*
    // Compute the 2D position of the strip at column "col" and line "lin",
    //  set the values in the "u, v" variables.
    // The pixel position depends on the mapping :
    //  1 = orthogonal pixel network (example = any square pixel)
    //  2 = staggered pixel network (example = MIMOSA 32ter, P31 pixel)
    //  3 = ? pixel network (example = MIMOSA 32, L82 pixel)
    //
    // Same as in DPlane.
    //
    // JB 2012/11/21
    // Modified: 2012/11/24 CLM correct version for P31 M32ter

    //int mapping = 1;// for p32 and for p1 //clm
    //int mapping = 2;// for p31 //clm
    //
    int mapping = Matrixtype; // JB 2013/07/17

    switch (mapping) {

      // When pixels are organized on an orthogonal grid
      case 1:
      default:
        u = col * PixelSizeU - FirstPixelShiftU;
        v = lin * PixelSizeV - FirstPixelShiftV;
        break;

        // When pixels are staggered from one column to the other
      case 2:
        u = col * PixelSizeU - FirstPixelShiftU;
        // v = (lin + 1./4. - (lin%2)/2.) * PixelSizeV - FirstPixelShiftV;
        // v = (lin - 1./4. + (lin%2)/2.) * PixelSizeV - FirstPixelShiftV; //clm
    better resolution first try

        // v center line = (lin * PixelSizeV - FirstPixelShiftV )
        // shift down by FirstPixelShiftV move up to center
        if ( col%2 == 0 ) v =  (lin * PixelSizeV - FirstPixelShiftV ) + 0.30 *
    PixelSizeV ;      //clm 2012.11.24 else              v =  (lin * PixelSizeV
    - FirstPixelShiftV ) - 0.19 * PixelSizeV ;      //clm 2012.11.24

        break;

        //clm mapping for M32 L8_2
      case 3:
        u = (col+1) * PixelSizeU - FirstPixelShiftU;
        v = lin * PixelSizeV - FirstPixelShiftV;
        break;

        // MonteCarlo Simulation. // LC 2014/01/10.
      case 4:
        u = ((2*col - NofPixelInRaw    + 2 ) * PixelSizeU)/2;
        v = ((2*lin - NofPixelInColumn + 2 ) * PixelSizeV)/2;
        break;

      case 5:
        {
    u = (2*col + 1 - NofPixelInRaw) * PixelSizeU/2.0;

    double fraction = 0.25;
    //double fraction = 0.75;
    if ( col%2 == 0 ) v = (NofPixelInColumn - 2*lin - 2*(fraction    )) *
    PixelSizeV/2.0; else              v = (NofPixelInColumn - 2*lin - 2*(1 -
    fraction)) * PixelSizeV/2.0; break;
        }
    case 7:  //M22-THRB6 and B7 -> set the matrix type and Mapping accordingly
    in the config file ! u = col * PixelSizeU - FirstPixelShiftU; if ( col%2 ==
    0 ) v =  (lin * PixelSizeV - FirstPixelShiftV ) - 0.50 * PixelSizeV ;
    //correct +0.5 else              v =  (lin * PixelSizeV - FirstPixelShiftV )
    - 0.0 * PixelSizeV ;     // correct -0.5
      //cout<<" MimosaAnalysis::ComputePixelPosition "<<u<<" "<<v<<" "<<col<<"
    "<<lin<<endl; break;

    }

  */

  return;
}

void MimosaAnalysis::ComputePixelPosition_UVToColRow(double u, double v,
                                                     double &col, double &lin)
{

  fTool.ComputeStripPosition_UVToColRow(Matrixtype, u, v, col, lin,
                                        NofPixelInRaw, NofPixelInColumn,
                                        PixelSizeU, PixelSizeV);

  /*
   // Compute the 2D position of a strip in the variables u and v to the set of
   variables column and line
   // The strip position depends on the mapping, look below.
   // AP 2012/06/26

   int mapping = Matrixtype;

   switch (mapping) {
     // When pixels are organized on an orthogonal grid
     case 1:
     default:
       col = (u/PixelSizeU) + ((NofPixelInRaw - 1)/2.);
       lin = (v/PixelSizeV) + ((NofPixelInColumn - 1)/2.);
       break;
     // When pixels are staggered from one column to the other
     case 2:
       col = (u/PixelSizeU) + ((NofPixelInRaw - 1)/2.);
       if(int(col)%2 == 0) lin = (v/PixelSizeV) + ((NofPixelInColumn - 1)/2.) -
   0.30; else                lin = (v/PixelSizeV) + ((NofPixelInColumn - 1)/2.)
   + 0.19; break; case 3:
       //Clm Mapping For M32 L8_2
       col = (u/PixelSizeU) + ((NofPixelInRaw - 1)/2.);
       lin = (v/PixelSizeV) + ((NofPixelInColumn - 1)/2.);
       break;
     // MonteCarlo Simulation. // LC 2014/01/10.
     case 4:
       col = (u/PixelSizeU) + ((NofPixelInRaw - 2)/2.);
       lin = (v/PixelSizeV) + ((NofPixelInColumn - 2)/2.);
       break;
     case 5:
       {
   col = (u/PixelSizeU) + ((NofPixelInRaw - 1)/2.);

   double fraction = 0.25;
   //double fraction = 0.75;
   if(int(col)%2 == 0) lin = ((NofPixelInColumn - 2*(fraction    ))/2.0)  -
   (v/PixelSizeV); else                lin = ((NofPixelInColumn - 2*(1 -
   fraction))/2.0)  -  (v/PixelSizeV); break;
       }
   case 7:  //M22-THRB6 and B7 -> set the matrix type and Mapping accordingly in
   the config file !
     //cout<<"  MimosaAnalysis::ComputePixelPosition_UVToColRow "<<u<<" "<<v<<"
   "<<col<<" "<<lin<<endl; col = (u/PixelSizeU) + ((NofPixelInRaw - 1)/2.);
     if(int(col)%2 == 0) lin = (v/PixelSizeV) + ((NofPixelInColumn - 1)/2.) +
   0.50;  //correct -0.5 else                lin = (v/PixelSizeV) +
   ((NofPixelInColumn - 1)/2.) + 0.0; //correct +0.5 break;
   }

 */

  return;
}

//______________________________________________________________________________
//
bool MimosaAnalysis::IsPixelIn1stCrown(int u, int v)
{
  //
  //  Check if pixel is in the 1st crown ... simple...
  //
  // clm 2013/01/23

  bool isPixOk = kFALSE;

  if ((u == 1 && v == 0) || (u == 1 && v == 1) || (u == 0 && v == 1) ||
      (u == -1 && v == 1) || (u == -1 && v == 0) || (u == -1 && v == -1) ||
      (u == 0 && v == -1) || (u == 1 && v == -1))
    isPixOk = kTRUE;

  return isPixOk;
}

//______________________________________________________________________________
//
bool MimosaAnalysis::IsPixelIn2ndCrown(int u, int v)
{
  //
  //  Check if pixel is in the 1st crown ... simple...
  //
  // clm 2013/01/23

  bool isPixOk = kFALSE;
  if ((u == 2 && v == 0) || (u == 2 && v == 1) || (u == 2 && v == 2) ||
      (u == 1 && v == 2) || (u == 0 && v == 2) || (u == -1 && v == 2) ||
      (u == -2 && v == 2) || (u == -2 && v == -1) || (u == -2 && v == 0) ||
      (u == -2 && v == -1) || (u == -2 && v == -2) || (u == -1 && v == -2) ||
      (u == 0 && v == -2) || (u == 1 && v == -2) || (u == 2 && v == -2) ||
      (u == 2 && v == -1))
    isPixOk = kTRUE;

  return isPixOk;
}

//______________________________________________________________________________
//
bool MimosaAnalysis::IsPlaneInLadder(int aPlaneNumber, int aLadderNumber)
{

  // Check the plane identified by aPlaneNumber
  // is in the list of planes of the ladder identified by aLadderNumber.
  // Return true or false.
  //
  // JB 2014/02/10

  DLadder *aLadder = fSession->GetTracker()->GetLadder(aLadderNumber);
  bool planeFound = kFALSE;

  int iplane = 0;
  while (iplane < aLadder->GetNumberOfPlanes() && !planeFound)
  {
    if (aLadder->GetPlane(iplane + 1)->GetPlaneNumber() == aPlaneNumber)
    {
      planeFound = kTRUE;
    }
    iplane++;
  }

  return planeFound;
}

//__________________________________________________________________________
//
void MimosaAnalysis::ProjectionImaging_init()
{
  // Init the variables required for spatial resolution studies (with bands)
  //
  // Called once for all
  //
  // JH 2014/07/21

  // for (int i=0; i<5; i++) NOfHitsMult[5];
}

//__________________________________________________________________________
//
void MimosaAnalysis::ProjectionImaging_Fill(DAuthenticHit *thehit,
                                            double angle)
{
  // Compute the variables required for spatial resolution studies (with bands)
  // and fill the corresponding histograms
  //
  // Called for each hit
  //
  // JH 2014/07/21

  double ImHu = (thehit->HuCG) * TMath::Cos(angle * TMath::Pi() / 180.0) +
                (thehit->HvCG) * TMath::Sin(angle * TMath::Pi() / 180.0);
  double ImHv = -(thehit->HuCG) * TMath::Sin(angle * TMath::Pi() / 180.0) +
                (thehit->HvCG) * TMath::Cos(angle * TMath::Pi() / 180.0);

  h1ProjectionOnX->Fill(ImHu);
  h1ProjectionOnY->Fill(ImHv);

  int multiplicity = thehit->HNNS;

  if (multiplicity == 1)
  {
    h1ProjectionOnXMult[0]->Fill(ImHu);
    h1ProjectionOnYMult[0]->Fill(ImHv);
    h2GoodHitsMult[0]->Fill(thehit->Hsu, thehit->Hsv);
    hdCGDigUVMult[0]->Fill(UofHitCG - hUdigital, VofHitCG - hVdigital);
    NOfHitsMult[0]++;
  }
  else if (multiplicity == 2)
  {
    h1ProjectionOnXMult[1]->Fill(ImHu);
    h1ProjectionOnYMult[1]->Fill(ImHv);
    h2GoodHitsMult[1]->Fill(thehit->Hsu, thehit->Hsv);
    hdCGDigUVMult[1]->Fill(UofHitCG - hUdigital, VofHitCG - hVdigital);
    NOfHitsMult[1]++;
  }
  else if (multiplicity == 3)
  {
    h1ProjectionOnXMult[2]->Fill(ImHu);
    h1ProjectionOnYMult[2]->Fill(ImHv);
    h2GoodHitsMult[2]->Fill(thehit->Hsu, thehit->Hsv);
    hdCGDigUVMult[2]->Fill(UofHitCG - hUdigital, VofHitCG - hVdigital);
    NOfHitsMult[2]++;
  }
  else if (multiplicity == 4)
  {
    h1ProjectionOnXMult[3]->Fill(ImHu);
    h1ProjectionOnYMult[3]->Fill(ImHv);
    h2GoodHitsMult[3]->Fill(thehit->Hsu, thehit->Hsv);
    hdCGDigUVMult[3]->Fill(UofHitCG - hUdigital, VofHitCG - hVdigital);
    NOfHitsMult[3]++;
  }
  else if (multiplicity == 5)
  {
    hdCGDigUVMult[4]->Fill(UofHitCG - hUdigital, VofHitCG - hVdigital);
    NOfHitsMult[4]++;
  }
  else
  {
    hdCGDigUVMult[5]->Fill(UofHitCG - hUdigital, VofHitCG - hVdigital);
    NOfHitsMult[4]++;
  }
}

//__________________________________________________________________________
//
void MimosaAnalysis::ProjectionImaging_end(int numberOfHits)
{
  // Finalize histograms created and filled for spatial resolution studies
  //
  // Called once for all
  //
  // JH 2014/07/21

  for (int i = 0; i < 5; i++)
  {

    NOfHitsMult[i] /= 1.0 * numberOfHits;
    h1NumberOfHitsMult->SetBinError(
        i + 1, 100.0 * sqrt(NOfHitsMult[i] * (1.0 - NOfHitsMult[i]) /
                            (1.0 * numberOfHits))); // Poisson's distribution

    NOfHitsMult[i] *= 100.0;
    h1NumberOfHitsMult->SetBinContent(i + 1, NOfHitsMult[i]);
  }

  hnGOODhit->SetLineColor(kBlue);
  hnGOODhit->SetMarkerColor(kBlue);
}
//__________________________________________________________________________
//
float MimosaAnalysis::GetTrackDistantToClosestDiode(float tu, float tv)
{

  // Gives the distance of the track to the closest diode at the DUT
  //(tu,tv) are the local coordinates of the track intersect with the DUT

  double distance = 1.0e+20;
  int col, lin;
  double col_p, lin_p;
  float u, v;

  // Find the pixel which the track intersects
  ComputePixelPosition_UVToColRow((double)tu, (double)tv, col_p, lin_p);
  col = int(col_p);
  lin = int(lin_p);

  // cout << endl;
  // Now loop over the pixels around the pixel which the track intersects and
  // look for the one with the closest diode's distance to the track
  for (int iii = -1; iii < 2; iii++)
  {
    // Looping on the colums to the left and to the right of the main pixel
    int col_tmp = col + iii;

    // Cut to ensure that the columns tested are inbetween 0 - Ncolumns-1
    if (col_tmp < 0)
      continue;
    if (col_tmp > NofPixelInColumn - 1)
      continue;
    for (int jjj = -1; jjj < 2; jjj++)
    {
      int lin_tmp = lin + jjj;

      // Cut to ensure that the lines tested are inbetween 0 - Nlines-1
      if (lin_tmp < 0)
        continue;
      if (lin_tmp > NofPixelInRaw - 1)
        continue;

      // cout << iii << "  " << jjj << endl;

      ComputePixelPosition(col_tmp, lin_tmp, u, v);
      float distance_tmp = sqrt(pow(tu - u, 2) + pow(tv - v, 2));

      if (distance > distance_tmp)
        distance = distance_tmp;
    }
  }
  // cout << endl;

  return distance;
}
//__________________________________________________________________________
//
