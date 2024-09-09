
// Pre-defined constants
#define SELEX_PRESPECIFIED 1
#define SELEX_COEFFICIENTS 2
#define SELEX_LOGISTIC 3
#define SELEX_KNIFE 4
#define SELEX_CONSTANT1 5
#define SELEX_LOGISTIC_OFFSET 6
#define SELEX_KNIFE_OFFSET 7
#define SELEX_CONSTANT_OFFSET 8
#define MOVE_NONE 0
#define MOVE_CONSTANT 1
#define MOVE_KNIFE 2
#define GROWTH_PRESPECIFIED 1

// Structure that contains the fixed constants
template <class Type>
struct dataSet{
  int Nyear;
  int Nstep;
  int Narea;
  int Nage;
  int Nsex;
  int Nfleet;
  int MaxLen;
  int NselPatterns;
  int NretPatterns;
  int NlegalPatterns;
  int BurnIn;
  vector<int> BurnInVec;
  int Num_Iteration;
  int Tune_Years;
  int First_yr;
  int MaxProjYr;
  int Nproj;
  int DoProject;
  vector<int> Nlen;
  matrix <int> SelSpec;
  matrix <int> RetSpec;
  matrix <int> LegalSpec;
  vector <int> Fleet_area;
  vector <int> Narea_fleet;
  matrix <int> Area_fleet;
  array <int> SelPnt;
  array <int> RetPnt;
  array <int> LegalFleetPnt;
  array <int> SelPntFut;
  array <int> RetPntFut;
  array <int> LegalFleetPntFut;
  array <int> IsRed;
  matrix <Type> TimeStepLen;
  array<Type> Catch;
  matrix<Type> MidLenBin;
  matrix<Type> LowLenBin;
  array <Type> Phi;
  int NmovePatterns;
  matrix <int> MoveSpec;
  array <int> MovePnt;
  int NrecruitPatternsA;
  matrix <int> RecruitSpecsA;
  int NrecruitPatternsB;
  matrix <int> RecruitSpecsB;
  matrix <int> RecruitPnt;
  vector <int> RecruitLenPnt;
  matrix <Type> RecruitFrac;
  int CalcRecruitFrac;
  int NgrowthPatterns;
  matrix <int> GrowthSpecs;
  array <int> GrowthPnt;
  array <Type> TransInp;
  int RecYr1;
  int RecYr2;
  int RecSpatYr1;
  int RecSpatYr2;
  int MatTimeStep;
  matrix<Type> MatFem;
  int InitOpt;
  Type Bias_Ramp_Yr1;
  Type Bias_Ramp_Yr2;
  Type Bias_Ramp_Yr3;
  Type Bias_Ramp_Yr4;
  int NefficPar;
};

// Data that appears in the likelihood
template <class Type>
struct TheData{
  int Ncpue;
  matrix<int> IndexI;
  matrix<Type> IndexR;
  int Nnumbers;
  matrix<int> NumbersI;
  matrix<Type> NumbersR;
  int NlenComp;
  matrix<int> LenCompI;
  vector<Type> Stage1W;
  matrix<Type> LenCompR;
  Type LambdaCpue;
  Type LambdaNumbers;
  Type LambdaLength;
  Type LambdaLarval;
  Type LambdaTag1;
  Type LambdaTag2;
  Type WeightInitialN;
  Type WeightInit3;
  int NcpueDataSeries;
  vector<int> IndexType;
  vector<int> FixSigmaCpue;
  vector<int> TreatQcpue;
  vector<int> EnvIndCpue;
  vector<int> EffCrIndCpue;
  vector<int> EffCrLag;
  Type SigmaCpueOffset;
  int NcatchDataSeries;
  vector<int> FixSigmaCatchN;
  Type SigmaCatchNOffset;
  vector<Type> LambdaCpue2;
  vector<Type> LambdaNumbers2;
  array<Type> LambdaLength2;
  int LarvalLikeOpt;
  int Larval_Offset;
  int NLarvalData;
  matrix<int> Lar_dataI;
  matrix<Type> Lar_dataR;
  array<Type> EnvData;
  Type InitialLoss;
  Type TagLossRate;
  int NrepSplit;
  vector<Type> RepRate;
  vector<int> FitTagSizes;
  int NtagLag;
  int NtagGroups;
  vector<int> Year1Tag;
  vector<int> Year2Tag;
  int TagYr1;
  int TagYr2;
  int NyearTags;
  array<Type> TagRel;
  array<Type> TagRec;
  array<Type> RecapObs;
  matrix<Type> NotReportedObs;
  matrix<Type> NrelTotal;
  array<Type>PropRepSplit;


};
