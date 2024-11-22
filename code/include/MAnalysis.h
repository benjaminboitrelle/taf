#pragma once

// Choose your environment
// #define LINUX
#define MAC
// #define WINDOWS

/////////////////////////////////////////////////////////////
//                                                         //
//  TAF - Main header file                                 //
//  Defines all user accessible analysis methods            //
//                                                         //
////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TControlBar.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "DEvent.h"
#include "DPrecAlign.h"
#include "DSession.h"
#include "TPaveLabel.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include "TError.h"
#include "TEnv.h"
#include "TInterpreter.h"

#include "DGlobalTools.h" // to have fTool has a data member

// pre
#include "DLadder.h" // JB 2014/02/10
#ifdef UseROOFIT
#include "DXRay2DPdf.h" // AP 2014/12/08
#endif
#include "DTracker.h"
#include "DPlane.h"
#include "DAlign.h"
#include "DHit.h"
#include "DR3.h"

// post
#include "TPad.h"
#include "TFrame.h"
#include "TText.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TLatex.h"
#include "TLegend.h"

// Results store
#include "TNtuple.h"

// align Mimosa
#include "MAlign.h"
#include "MHist.h"

// MRaw/x
#include "MRaw.h"
#include "MRax.h"

//---ADC
#include "DSetup.h"
//---ADC

#include <string>

class MRaw;
class MRax;

//---ADC
class DSetup; // forwards
//---ADC

class DR3;

class MimosaAnalysis : public MHist
{

  //====================
private:
  //====================
  static MimosaAnalysis *fgInstance; // global pointer to the TAF object
#ifdef MAC
  MimosaAnalysis *gTAF; // for MAC
#endif

  // ------------------------
  // Files and directories
  // ------------------------
  TFile *fileInputTree;
  ofstream outFileGoodEvt;    //!
  ofstream outFileBadEvt;     //!
  ofstream csvfile;           //!
  TFile *theCorFile;          //!
  TString DTDIR;              // set from the environment variable DTDIR = TAF installation directory
  TString fLaunchedDirectory; // directory where taf has been launched // VR 2014/06/30
  TString fWorkingDirectory;  // directory where taf will work, choose in MimosaAnalysis::MimosaAnalysis() // VR 2014/06/30
  // TString ROOTDIR; // removed by JB 2011/04/12
  TString ResultDirName; // JB 2011/04/12
  DGlobalTools fTool;    // JB 2011/07/18
  char RootFile[300];    // name of the input DSF file

  // ------------------------
  // Input tree
  // ------------------------
  // JB 2010/08/30
  TTree *t;
  DEvent *Evt;
  TBranch *branch;

  // ------------------------
  // Parameters for sensor, alignment, geometry, analysis
  // ------------------------

  DSession *fSession;   // JB 2011/07/21 to get rid of global var tSession
  double trU, trV, trW; // alignment
  double theta[3];
  DPrecAlign *alignement;

  float calibration;
  float NoiseScope;

  int MimosaType;
  int RunNumber;
  int ThePlaneNumber;
  int ThePlaneNumber2; // JB 2011/10/30
  int TheLadderNumber; // JB 2014/01/10
  int Nevt, User_MaxEvt;

  int ThesubmatrixNumber;
  int ThesubmatrixNumber2; // JB 2015/12/15

  float MimosaSizeU, MimosaSizeV;
  float PixelSize;
  float PixelSizeU;
  float PixelSizeV;
  int Matrixtype;
  int NofPixelInColumn;
  int NofPixelInRaw;
  int NofPixels;
  float FirstPixelShiftU, FirstPixelShiftV;
  int Thegeomatrix;
  int Thegeomatrix2; // JB 2015/12/15
  //  float geomUmin, geomUmax;
  //  float geomVmin, geomVmax;
  float geomUmin2, geomUmax2;
  float geomVmin2, geomVmax2;
  float exgeomUmin, exgeomUmax; // JB 2013/06/22
  float exgeomVmin, exgeomVmax;

  int UserFlag; // parameter transmitted through config file, free for user, JB 2013/07/17

  // ------------------------
  // flags
  // ------------------------
  int UsedMacro[50]; // to be modified
  int MimoDebug;
  int fWriteHistos;
  int fUserFileNumber;
  bool fInitDone;
  bool fMimosaDone; // for general purpose histo JB 2010/09/10
  bool fMimosaProDone;
  bool fMimosaFakerateDone;    // JB 2010/07/26
  bool fMimosaCalibrationDone; // JB 2010/07/27
  bool fMimosaMiniVectorsDone; // JB 2010/08/30
  bool fMimosaPro2PlanesDone;  // JB 2013/05/01
  bool fMimosaImagingDone;     // JH 2014/06/16
  bool fClearDone;
  // The following flags indicate the purpose of the analysis
  // JB 2014/01/11
  bool fIfReferenceTrack;
  bool fIfCalibration;
  bool fIfImaging;
  bool fIfLaserSpot;
  bool fIfMiniVector;
  bool fIfVertex;
  bool fIfFake;
  bool fIfBinary;       // for binary output sensor, JB 2014/03/15
  bool fIfMCGeneration; // for MC generation, AP 2015/04/02

  // ------------------------
  // Hot pixel map
  // ------------------------
  int TheUsePixelMap;        // use pixel map or not
  int Option_read_Pixel_map; // Hot hit map read=1/write=0 option
  TH2F *h2HotPixelMap;       // pointer to the histogram containing pixel hit rate
  Char_t HotPixelFileName[40];
  TFile *HotPixelFile; // pointer to file containing the previous histogram

  // ------------------------
  // cuts
  // ------------------------
  float TrackChi2Limit;
  int MinHitsPerTrack;      // JB 2013/06/22
  int MaxNbOfTracksInGeom;  // JB 2012/08/30
  int GeoMatrixForTrackCut; // JB 2012/08/30
  int CUT_MaxNbOfHits;
  int CUT_MinNbOfHits;
  float TrackToHitDistanceLimit;
  float TrackToHitDistanceLimit2; // JB 2015/12/15
  float CUT_S2N_seed;
  float CUT_S2N_neighbour;
  float CUT_Q_seed;         // JB 2013/11/08
  float CUT_Q_cluster;      // JB 2014/01/21
  float CUT_MinQ_neighbour; // JB 2010/07/27
  float CUT_MaxQ_neighbour; // JB 2013/11/08
  double CUT_MaxHitRatePerPixel;
  double CUT_MinHitRatePerPixel;
  int CUT_MinSeedIndex; // JB 2013/08/21
  int CUT_MaxSeedIndex; // JB 2013/08/21
  int CUT_MinSeedCol;   // JB 2013/08/22
  int CUT_MaxSeedCol;   // JB 2013/08/22
  int CUT_MinSeedRow;   // JB 2013/08/22
  int CUT_MaxSeedRow;   // JB 2013/08/22

  // ------------------------
  // counters & efficiency
  // JB 2011/11/04
  // ------------------------
  double MimosaEfficiency;
  double MimosaEfficiency_ERR;

  int NofClMatchTrack; // # good tracks in the DUT WITH a matched hit
  int NtrkInMimo;      // Number of tracks in the MIMOSA aperture

  int NeventRangeForEfficiency; // binning for the efficiency during the run
  int temp_maxarray;
  int ievt_array; // counter of bunches
  float *temp_Efficiency_array;
  float *temp_NofClMatchTrack;
  float *temp_NtrkInMimo;

  // ------------------------
  //
  // ------------------------

  float *ResolutionPoint; // better in MHIst ??

  // ------------------------
  // pixel data inside cluster, JB 2010/06/03
  // ------------------------
  int NofPixelsInCluster;
  float *qonly;
  float *snonly;
  float *noise;
  float *UofPix;
  float *VofPix;
  int *LineInCluster;
  int *ColumnInCluster;
  int *IndexofPix;
  int *IndexInCluster;
  int LineSizeOfCluster; // JB 2014/03/31
  int ColumnSizeOfCluster;

  typedef enum
  {
    one_pix = 1,
    two_pix_row,
    two_pix_col,
    three_pix_L,
    three_pix_row,
    three_pix_col,
    four_pix_square,
    four_pix_L_row,
    four_pix_L_col,
    four_others,
    more_than_four,
    five_pix_squarerow,
    five_pix_squarecol,
    five_others,
    six_pix_3col2row,
    six_pix_2col3row,
    six_others,
    more_than_six
  } clusterTypeDef;
  clusterTypeDef ClusterGeometricalType;
  clusterTypeDef ClusterGeometricalTypeBeyond4;

  int IndexOfMaxPixel; // index of the pixel with the highest charge (in qonly array)

  float snOptimal; // Largest S/N obtained by summing q ordered pixel
  int optimalSize; // number of pixels in the largest S/N

  float MinUofPix; // position of lower left pixel
  float MinVofPix;

  float TotalCharge;      // all the pixels in the cluster
  float ChargeAroundSeed; // all the pixels except the seed
  float ChargeInCross;    // the 4 pixels direct neighbours of the seed
  float ChargeInX;        // the 4 pixels in the corner of the seed
  float ChargeInCrown1;   // the 16 pixels after the first neighbours of the seed
  float ChargeInCrown2;   // the 24 pixels next to next neighbours of the seed

  float TotalCharge2x2;    // Cluster charge limited to 2x2 sub-cluster
  float ExternalCharge2x2; // 2x2 cluster charge without the seed (Qof3x3[0])
  float Qof2x2[4];
  float UofPix2x2[4];
  float VofPix2x2[4];
  int IndexU2x2[4];
  int IndexV2x2[4];
  short IndexOfCluster2x2;

  float TotalCharge3x3; // Cluster charge limited to 3x3 sub-cluster
  float Qof3x3[9];
  float UofPix3x3[9];
  float VofPix3x3[9];
  int IndexofPix3x3[9];

  float TotalCharge5x5; // Cluster charge limited to 5x5 sub-cluster
  float Qof5x5[25];
  float UofPix5x5[25];
  float VofPix5x5[25];
  int IndexofPix5x5[25];

  float *qqordered; // charge ordered pixels
  float *nqordered;
  float *snqordered;
  int *IndexofPixqordered;
  float *UofPixqordered;
  float *VofPixqordered;
  int *LineInClusterqordered;
  int *ColumnInClusterqordered;

  float *q; // cumulated q,n,sn with charge ordered pixels
  float *n;
  float *sn;

  float *qsnordered; // sn ordered pixels
  float *nsnordered;
  float *snsnordered;
  int *IndexofPixsnordered;
  float *UofPixsnordered;
  float *VofPixsnordered;
  int *LineInClustersnordered;
  int *ColumnInClustersnordered;

  float *qsumsn; // cumulated q,n,sn with s/n ordered pixels
  float *nsumsn;
  float *snsumsn;

  float *qspiral; // geometry ordered pixels
  float *nspiral;
  float *snspiral;
  int *IndexofPixspiral;
  float *UofPixspiral;
  float *VofPixspiral;
  int *LineInClusterspiral;
  int *ColumnInClusterspiral;

  float *qsumspiral; // cumulated q,n,sn with s/n ordered pixels
  float *nsumspiral;
  float *snsumspiral;

  // ------------------------
  // Track information, JB 2010/06/03
  // ------------------------
  float TrackToHitDistance;
  float TrackToHit2ndDistance;
  float chi2;
  float tu; // track position in plane
  float tv;
  float tw;
  float tdu; // track slope in plane coord.
  float tdv;
  float tk1; // index of nearest strip to seed
  float tx;  // track position in telescope
  float ty;
  float tz;
  float tdx; // track slope in tracker coord.
  float tdy;

  // ------------------------
  // MiniVector information, JB 2010/08/30
  // ------------------------
  DR3 trackPosUVW[2];
  DR3 trackPosXYZ[2];
  DR3 trackMeanPosXYZ;
  double trackAngleXZ;
  double trackAngleYZ;

  DR3 hitPosUVW[2];
  DR3 hitPosXYZ[2];

  DR3 vectorPosXYZ;
  double vectorSlopeXZ;
  double vectorSlopeYZ;
  double vectorAngleXZ;
  double vectorAngleYZ;

  // ------------------------
  // cluster positions
  // ------------------------
  int hitCounterPos; // JB 2014/01/10

  float hX, hY, hZ; // position in tracker frame, JB 2010/08/30
  float hU;         // position in plane frame
  float hV;
  float hW;
  float hUdigital; // seed or digital position
  float hVdigital;
  float UCG2x2; // CoG 2x2 ??
  float VCG2x2;
  float UofHitCG3; // CoG using only 3x3 pixels
  float VofHitCG3;
  float UofHitCG5; // CoG using only 5x5 pixels
  float VofHitCG5;
  float UofHitCG; // CoG using all pixels
  float VofHitCG;
  float UCGcorr; // CoG corrected by correlation plot
  float VCGcorr;
  float UofHitEta3; // Eta algo.
  float VofHitEta3;
  float Eta2x2U;
  float Eta2x2V;
  float UofHitEta2x2;
  float VofHitEta2x2;
  float UofHitEta2x2_new;
  float VofHitEta2x2_new;
  float UofHitEta5x5_new;
  float VofHitEta5x5_new;
  float UofHitEta2x2_newR;
  float VofHitEta2x2_newR;
  float UofHitEta5x5_newR;
  float VofHitEta5x5_newR;
  double Uaht;
  double Vaht;

  // ------------------------
  // used in getmieta and create new eta
  // ------------------------
  int NBins2x2;
  float *Edges2x2;     //!
  float *Contents2x2U; //!
  float *Contents2x2V; //!
  float *Edges3x3;     //!
  int NBins3x3;
  float *Contents3x3U; //!
  float *Contents3x3V; //!
  static const int PolDeg = 7;
  float FitParEta3U[PolDeg];
  float FitParEta3V[PolDeg];
  float FitParamU[PolDeg];
  float FitParamV[PolDeg];

  // ------------------------
  // For cluster shape studies, cdritsa & JB 2010/04/13
  // ------------------------
  int hitCounter;

  int nThresholds; // <10
  float thresholds[10];
  float snrThresholds[10];
  int clusterMultiplicity[10];
  int clusterMultiplicityCounter[10];
  float multV[10];
  float multU[10];
  float mult[10];
  float formFactor[10];

  double ChargeTimesPositionU[50];
  double ChargeTimesPositionV[50];
  double ChargeTimesPositionUmean[50];
  double ChargeTimesPositionVmean[50];

  double ChargeInCluster[50];
  double sumChargeTimesPositionU[50];
  double sumChargeTimesPositionV[50];
  double sumChargeTimesPositionUmean[50];
  double sumChargeTimesPositionVmean[50];

  double RMSV;
  double RMSU;
  double RMSVmean;
  double RMSUmean;
  // end for cluster shape

  int NPages;

  int NTracksPerEventInSensor;
  double AverageEffic[2];

  int GetFileNumber();
  void GetParameters();
  void GetAnalysisGoal();                                                             // JB 2014/01/17
  void ComputePixelPosition(int col, int lin, float &u, float &v);                    // JB 2012/11/21
  void ComputePixelPosition_UVToColRow(double u, double v, double &col, double &lin); // JH 2014/07/30
  void PrepareOnlineDisplay();                                                        // JB 2010/08/30
  void BookingHistograms() { MHist::BookingHistograms(RunNumber, PixelSizeU, PixelSizeV, NofPixelInRaw, NofPixelInColumn, TrackToHitDistanceLimit, fSession->GetSetup(), fIfReferenceTrack, fIfCalibration, fIfLaserSpot, fIfMiniVector, fIfVertex, fIfFake, fIfImaging); };
  void CreateNewEta();
  bool TrackInMimo(int aGeoMatrix, float tuVal, float tvVal, int aSubMatrix = 0);
  // void Zero(TObject** ptr,int  size=0) {/*for(int ind=0;ind<size;ind++){ if(!ptr) ptr=0;}*/};
  void Zero(TObject *ptr)
  {
    if (ptr)
      delete ptr;
    if (!ptr)
      ptr = 0;
  };
  const char *CreateGlobalResultDir();
  void PreparePost();
  int CorStatus; // corpar file status.
  void SetCalibration(float acalibration = 1) { calibration = acalibration; }
  void CreateConfig(int RunNumberBigin, int NumberOfFilesPerRun = 0, int RunNumberEnd = 0, int RunNumberStep = 1);
  void CreateLinks(int RunNumberBigin, int NumberOfFilesPerRun = 0, int RunNumberEnd = 0, int RunNumberStep = 1);
  int OpenInputFile(); // JB 2010/09/30

  // ************************************
  // Analysis functions

  void GoodHit_Fill(DAuthenticHit *ahit); // JB 2010/10/06

  void ClusterCharges_init(); // JB 2010/06/03
  void ClusterCharges_compute(DAuthenticHit *thehit);
  void ClusterCharges_fill(DAuthenticHit *thehit, int ievt);
  void ClusterCharges_end();

  void ClusterPosition_init();                                            // JB 2014/01/10
  void ClusterPosition_compute(DAuthenticHit *thehit, DPrecAlign *align); // align added, JB 2013/05/01
  void ClusterPosition_cog(DAuthenticHit *thehit);
  void ClusterPosition_eta(DAuthenticHit *thehit);
  void ClusterPosition_aht(DAuthenticHit *thehit);
  void ClusterPosition_fill(DAuthenticHit *thehit);
  void ClusterPosition_end(); // JB 2014/01/10

  void TrackParameters_init() { ; } // JB 2010/07/23
  void TrackParameters_compute(DTransparentPlane *atrack, DPrecAlign *align);
  void TrackParameters_allFill(DTransparentPlane *atrack, int ievt);
  void TrackParameters_goodFill(DTransparentPlane *atrack, int ievt);
  void TrackHitPosition_fill(DAuthenticHit *thehit, DTransparentPlane *atrack, int ievt);
  void TrackParameters_end(); // JB 2014/01/17

  void ClusterShape_init(); // JB 2010/04/13
  void ClusterShape_fill(DAuthenticHit *thehit);
  void ClusterShape_end();

  void RefPlane_init() { ; }
  void RefPlane_fill(); // JB 2010/05/27
  void RefPlane_end() { ; }

  void BinarySpecific_init() { ; } // JB 2010/06/03
  void BinarySpecific_HitsOnly_fill(DAuthenticHit *thehit);
  void BinarySpecific_fill(DAuthenticHit *thehit);
  void BinarySpecific_end();

  // FakeRate_compute and _fill are kept for consistency
  // but there are so far useless because _end uses histos
  // filled by
  // JB 2010/07/22
  void FakeRate_init() { ; }
  void FakeRate_compute() { ; }
  void FakeRate_fill() { ; }
  void FakeRate_end(int rateNormalisation);

  // MiniVector analysis
  // JB 2010/08/30 (from NCS 2009/10)
  void MiniVector_init() { ; }
  void MiniVector_compute();
  void MiniVector_fill(DAuthenticHit *hit1, DAuthenticHit *hit2);

  // Proba of fake association with inefficient events
  void Inefficent_end();

  // Efficiency computation
  // JB 2011/11/04 (from AB)
  void Efficiency_init();
  void Efficiency_end(int eventsRead);

  // Hot pixel removing
  void HotPixel_init(int useMap);
  int HotPixel_test(int aPixelIndex);
  void HotPixel_end(int eventsRead);

  // Check what crown the pixel belongs to //clm 2013/01/23
  bool IsPixelIn1stCrown(int lin, int col);
  bool IsPixelIn2ndCrown(int lin, int col);

  // Imaging analysis
  double NOfHitsMult[5];
  void ProjectionImaging_init();
  void ProjectionImaging_Fill(DAuthenticHit *thehit, double angle);
  void ProjectionImaging_end(int numberOfHits);

  //====================
public:
  //====================
  MimosaAnalysis();
  virtual ~MimosaAnalysis() { ; }
  bool CheckIfDone(const Option_t *Option);

  static MimosaAnalysis *&Instance()
  {
    if (!fgInstance)
    {
      // cout<<"No current TAF session detected! Creatung a new one..."<<endl;
      fgInstance = new MimosaAnalysis();
    }
    return fgInstance;
  }

  void MergingPDFFiles(void);

  // ************************************
  // Internal tool methods

  void GetMiEta();
  // void       StoreEta() ;  // replaced by CreateNewEta, JB 2011/06/19
  void GetAlignment(); // JB 2011/06/19
  int GetHitMapReadOpt() { return Option_read_Pixel_map; }
  // void       SettHitMapReadOpt(int anOpt)      {Option_read_Pixel_map = anOpt;}
  int GetDebug() { return MimoDebug; }
  int GetRunNumber() { return RunNumber; }
  DSession *GefSession() { return fSession; }
  void SetRunNumber(int aRunNumber) { RunNumber = aRunNumber; }
  void InitRunNumber() { RunNumber = fSession->GetRunNumber(); }
  int GetPlaneNumber() { return ThePlaneNumber; }
  void SetPlaneNumber(int aPlaneNumber);
  void SetLadderNumber(int aNumber); // JB 2014/01/10
  void InitPlaneNumber() { ThePlaneNumber = fSession->GetPlaneNumber(); }
  int GetMimosaType() { return MimosaType; }
  void SetMimosaType(int aMimosaType) { MimosaType = aMimosaType; }
  int GetAlignStatus() { return fSession->GetTracker()->GetAlignmentStatus(); }
  void SetDebug(int aMimoDebug);
  void SetAlignStatus(int aStatus) { fSession->GetTracker()->SetAlignmentStatus(aStatus); } // See DTracker

  double GetCUT_MaxHitRatePerPixel() { return CUT_MaxHitRatePerPixel; }
  double GetCUT_MinHitRatePerPixel() { return CUT_MinHitRatePerPixel; }
  void SetCUT_MaxHitRatePerPixel(double aRate) { CUT_MaxHitRatePerPixel = aRate; }
  void SetCUT_MinHitRatePerPixel(double aRate) { CUT_MinHitRatePerPixel = aRate; }

  void InitMimosaType()
  {
    MimosaType = (int)(RunNumber / 1000.);
    if (RunNumber == 2110)
      MimosaType = 4;
    else if (RunNumber < 1000)
      MimosaType = 99;
  } // JB 2012/05/04 deal with RunNumber<1000 case

  void SetFileNumber(int aUserFileNumber = 0) { fUserFileNumber = aUserFileNumber; };
  void SetDSFFile(const Char_t *aFileName) { sprintf(RootFile, "%s", aFileName); }; // JB 2013/09/19
  DPrecAlign *AlignMimosa(int aDistance);
  // DPrecAlign*  AlignPrecMimosa(DPrecAlign* initAlignment, bool modeAuto); // LC 2012/09/06.
  // TFile*     InitCorPar(int aRun, int aPlane, int &aStatus);
  // void       InitCorPar(int aRun, int aPlane); // Jb 2009/05/18
  void InitCorPar(int aRun, int aPlane, const char *objType = "plane"); // JB 2014/02/10
  bool IsPlaneInLadder(int aPlaneNumber, int aLadderNumber);            // JB 2014/02/10

  int ResultsStrore();

  float GetCalibration() { return calibration; };

  // prepare
  void PrepareRun(int RunNumberBigin, Option_t *aType, int NumberOfFilesPerRun = 0, int RunNumberEnd = 0, int RunNumberStep = 1);

  float GetMimosaEfficiency() { return MimosaEfficiency; }
  float GetMimosaEfficiency_ERR() { return MimosaEfficiency_ERR; }

  int GetNtrkInMimo() { return NtrkInMimo; } // Number of tracks in the MIMOSA aperture

  TFile *GetCorParFile()
  {
    if (!theCorFile)
    {
      InitCorPar(GetRunNumber(), GetPlaneNumber());
    }
    return theCorFile;
  }

  //----ADC
  float tVect[100]; // n res
  // int       *fChannel;                  //! pointer to Channels
  // int&       GetChannel(int aSk)              { return  fChannel[aSk-1];     }
  float *GettVect() { return tVect; }; // n res
  //----ADC

  // MC simulation parameters
  bool UsingTrackerResolution;
  int Bins_GlobalResolution;
  double DUT_residualU_CGDSF[2];
  double DUT_residualV_CGDSF[2];
  double DUT_resolutionU_CGDSF[2];
  double DUT_resolutionV_CGDSF[2];
  double Tel_resolutionU_CGDSF[2];
  double Tel_resolutionV_CGDSF[2];
  double From_mu_to_mm;
  int _PrintFreq;
  std::vector<int> _PlaneList; // Plane indexes in orderer by Z coordinate, from lower to higher.
  std::vector<float> _PlaneListZ;
  std::vector<float> _PlaneThickness;
  std::vector<TString> _PlaneMaterial;
  std::vector<float> _RU_min;
  std::vector<float> _RU_max;
  std::vector<float> _RV_min;
  std::vector<float> _RV_max;
  double MyRX[2];
  double MyRY[2];
  double MyRZ[2];
  double DUTGeomMatrixX[2];
  double DUTGeomMatrixY[2];
  double RGlobalResolution[2];
  Color_t color[500];
  double Nsigma;
  int *TracksPerEvt;
  int *NRectracks;

  // ************************************
  // Arrays for eta correction

  TArrayF Liste_CoGU;
  TArrayF Liste_CoGV;
  int nListe_CoG;
  TArrayF READListe_CoGU;
  TArrayF READListe_CoGV;
  int READnListe_CoG;

  TArrayF Liste_CoGU_eta2x2;
  TArrayF Liste_CoGV_eta2x2;
  int nListe_CoG_eta2x2;
  TArrayF READListe_CoGU_eta2x2;
  TArrayF READListe_CoGV_eta2x2;
  int READnListe_CoG_eta2x2;

  TArrayF Liste_CoGU_eta5x5;
  TArrayF Liste_CoGV_eta5x5;
  int nListe_CoG_eta5x5;
  TArrayF READListe_CoGU_eta5x5;
  TArrayF READListe_CoGV_eta5x5;
  int READnListe_CoG_eta5x5;

  //  int GetREADnListe_CoG()                    { return READnListe_CoG;                  }
  //  void  SetREADnListe_CoG(int aREADnListe_CoG)    { READnListe_CoG  = aREADnListe_CoG;           }

  // ************************************
  // Comonly used commands

  void CDHist()
  {
    if (!CheckIfDone("mimosapro"))
    {
      return;
    }
    else
    {
      MHist::Dir();
    }
  }
  void Help();

  void InitSession(const int TheRun, const unsigned int ThePlane = 0, int EventBuildingMode = -1, TString ConfigFile = "", TString ConfigDir = "", TString DataDir = "", TString DataFile = "");
  // EventBuildingMode key by default is at -1 meaning that it will not override the value taken in DSetup::ReadConfiguration(), SS 2011/11/14
  //  VR 2014/06/30, JB 2021/11/15 Add optional args

  // ************************************
  // Raw data reconstruction commands

  void MakeEta(int NEvt = 10000);

  void SetTrackGeoLimitsForAlign(double xmin, double xmax, double ymin, double ymax) { fSession->SetTrackGeoLimitsForAlign(xmin, xmax, ymin, ymax); } // JB 2013/06/11
  void SetTrackChi2LimitForAlign(double aLimit) { fSession->SetTrackChi2LimitForAlign(aLimit); }                                                      // JB 2013/07/14
  void AlignTracker(const double tiniBoundU,
                    const double tiniBoundV,
                    int nAlignEvents = 4000,
                    int nAdditionalEvents = 2000,
                    bool UseAllHits = false); // 2013/06/10.
                                              // AP 2015/06/10: added bool parameter (UseAllHits) to decide if doing alignement with all hits or the closest one to the track.
  void AlignTracker(const double tiniBound = 480.,
                    int nAlignEvents = 4000,
                    int nAdditionalEvents = 2000,
                    bool UseAllHits = false) { AlignTracker(tiniBound, tiniBound, nAlignEvents, nAdditionalEvents, UseAllHits); }
  void AlignTrackerMinuit(bool modeAuto = 1,
                          const double tiniBound = 480.,
                          int nAlignEvents = 5000,
                          int nAlignHitsInPlane = 4000,
                          int nAdditionalEvents = 2000,
                          double chi2Limit = 0.,
                          bool UseAllHits = true); // LC 2012/09/06
                                                   // AP 2015/06/08 : adding bool parameter (UseAllHits) to decide if doing alignment with all hits or the closest one
  void AlignTrackerMinuitLadder(bool modeAuto = 1, const double tiniBound = 480., int nAlignEvents = 5000, int nAlignHitsInPlane = 4000, int nAdditionalEvents = 2000, double chi2Limit = 0.);
  // void       AlignLadder(bool modeAuto=1, const double tiniBound=480., int nAlignEvents=5000, int nAlignHitsInPlane=4000, int nAdditionalEvents=2000, int ladderFace=1, double chi2Limit=0.);  // Commented 2014/12/18
  void AlignLadderMV(bool modeAuto = 1, const double tiniBound = 480., const double boundSlopes = 100., int nAlignEvents = 5000, int nGoodTracks = 1000, int nAdditionalEvents = 2000, int chi2Limit = 40, int mode = 0);   // LC 2013/09/10
  void AlignTrackerGlobal(const int refPlane = 1, const int nEvents = 4000, const int nIterations = 10, const bool alignConstrainsts = true, const bool trackConstrainsts = true, const bool multipleScatteringFit = true); // LC && LiuQ 2015/02/06.
  void AlignTrackerMillepede(int nAlignEvents = 4000);                                                                                                                                                                      // LC 2012/12/24.
                                                                                                                                                                                                                            //  void        Gener(double* xp, double* yp, double& aX, double& bX, double& aY, double& bY, double* sigX, double* disX, double* sigY, double* disY, double* z, double* phi); // LC 2012/01/07

  void DSFProduction(int NEvt = 500000, int fillLevel = 1);
  void StudyDeformation(const float tiniBound = 480., int nEvents = 2000, bool fitAuto = 0); // BB 2014/05/20

  // ************************************
  // Analysis commands

  DSession *GetSession() { return fSession; } // JB 2013/06/23
  MRaw *GetRaw();
  MRax *GetRax();
  void MimosaPro(int MaxEvt, float TrackHitDist, float S2N_seed, float S2N_neighbour, int submatrix, int GeoMatrix, Option_t *SaveAlign, int UseHitMap = 0, int WriteMissedHits = 0, int WriteGoodHits = 0);
  void MimosaFakerate(int MaxEvt, float S2N_seed, float S2N_neighbour, int submatrix = 0, int GeoMatrix = 0);
  void MimosaCalibration(int MaxEvt, float S2N_seed = 0., float S2N_neighbour = 0., float ChargeAroundSeedCutMax = 5000., int submatrix = 0, int GeoMatrix = 0);
  void MimosaMiniVectors(int MaxEvt, int TrackHitDist, short plane0, float S2N_seed0, float S2N_neighbour0, int submatrix0, int GeoMatrix0, short plane1, float S2N_seed1, float S2N_neighbour1, int submatrix1, int GeoMatrix1);
  void MimosaPro2Planes(int MaxEvt, int TrackHitDist, short plane0, float S2N_seed0, float S2N_neighbour0, int submatrix0, short plane1, float S2N_seed1, float S2N_neighbour1, int submatrix1, int GeoMatrix);
  void MimosaVertex(int MaxEvt, short plane1, short plane2);
  void MimosaVertexFinder(int MaxEvt, int submatrix, int GeoMatrix, float chi2Cut);
  void MimosaCluster(int MaxEvt, float S2N_seed = 0., float S2N_neighbour = 0., int submatrix = 0, int GeoMatrix = 0);                                                                               // JB 2014/01/17
  void MimosaProLadder(int MaxEvt, int TrackHitDist, float S2N_seed, float S2N_neighbour, int submatrix, int GeoMatrix, Option_t *SaveAlign, int UseHitMap, int WriteMissedHits, int WriteGoodHits); // JB 2014/02/1

  void MimosaImaging(int MaxEvt, int submatrix = 0, int GeoMatrix = 0, double theta_init = 0.0, int Npeak = 1, double S = 100, double W = 100, bool FirstLoop = kTRUE, bool chooseFit = kTRUE); // JH 2014/06/13

  void MimosaGeneration_ToyMC(int MaxEvt = 10000,
                              int PlaneNumber = 1,
                              int Submatrix = 0,
                              int Geomatrix = 0,
                              bool DoGaussianMS = true,
                              int seed = 8292894,
                              bool DoDisplay = false,
                              double XMeanDiv = 0.0,
                              double XRMSDiv = 1.0e-20,
                              double YMeanDiv = 0.0,
                              double YRMSDiv = 1.0e-20,
                              bool CalledFromMimosaPro = false,
                              bool verbose = false); // AP 2015/03/11

  void MimosaGeneration_LoicG4Simu(int MaxEvt = 10000,
                                   int PlaneNumber = 1,
                                   int Submatrix = 0,
                                   int Geomatrix = 0,
                                   int seed = 8292894,
                                   bool DoDisplay = false,
                                   bool verbose = false); // AP 2015/07/08

  void MimosaGeneration_APG4Simu(int MaxEvt = 10000,
                                 int PlaneNumber = 1,
                                 int Submatrix = 0,
                                 int Geomatrix = 0,
                                 int seed = 8292894,
                                 bool DoDisplay = false,
                                 bool verbose = false); // AP 2015/07/29

  void BookingMC(int PlaneNumber = 1,
                 int Submatrix = 0,
                 int Geomatrix = 0,
                 bool CalledFromMimosaPro = false); // AP 2015/07/29
  void DoMCDisplay(int NDisplay, bool verbose);     // AP 2015/07/29
  void DoMCDisplay3D(int NDisplay, bool verbose,
                     const char *file_gdml); // AP 2015/07/31
  double DistanceOfPointToPlane(const int PlaneNumber,
                                DR3 Point); // AP 2015/08/03
  DR3 ClosestPointInPlaneToPoint(const int PlaneNumber,
                                 DR3 Point); // AP 2015/08/03

  float GetTrackDistantToClosestDiode(float tu,
                                      float tv);

  // ************************************
  // Plotting functions

  void MimosaResolution();
  void CheckMCSimulation(bool StandAlone = false);
  void MimosaClusterCharge(TString fitOption = "landau");
  void MimosaBinaryplots();
  void CheckMimosaAlign();
  void CheckAsymmetry();
  void InspectEfficiency();
  void checkRefPlane();
  void CompareBadGoodRaw();
  void CheckClusters();
  void Checkreferencetracks();
  void pixelcharge();
  void MimosaSN();
  void ClusterShape(); // cdritsa feb 2008, JB 2010
  void HitMap();       // JB 2010/10/06
  void MimosaPixelHomogeneity();
  void MimosaOptimize();
  void FakeRate();
  void Calibration(); // JB 2010/07/27
  TPad *PlotCalibration(TH1F *h, int manualFitRange = 0);
  void MiniVectors(); // Jb 2010/08/21
  void UserAnalysis();
  void SaveGifFile();
  void Clear(Option_t * /*option = ""*/);
  void CheckImaging(double theta_init = 62.0, int Npeak = 1, double S = 100, double W = 100, bool FirstLoop = kTRUE, bool chooseFit = kTRUE); // JH 2014/06/16

  // ************************************

  ClassDef(MimosaAnalysis, 2)
};

// R__EXTERN MimosaAnalysis *gTAF; // JB, for inline compilation

#ifndef MAC
// #define gTAF MimosaAnalysis::Instance() // to comment for MAC
#endif
