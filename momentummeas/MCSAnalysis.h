#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

// This is what we are here for
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// Standard ROOT stuff
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "Math/RootFinder.h"
#include "TLine.h"

// Read directly from the MAUS data structure.
#include "src/common_cpp/DataStructure/TOFEvent.hh"
#include "src/common_cpp/DataStructure/SciFiEvent.hh"
#include "src/common_cpp/DataStructure/SciFiTrack.hh"
#include "src/common_cpp/DataStructure/SciFiTrackPoint.hh"
#include "src/common_cpp/DataStructure/SciFiSpacePoint.hh"
#include "src/common_cpp/DataStructure/KLEvent.hh"
#include "src/common_cpp/DataStructure/CkovEvent.hh"
#include "src/common_cpp/DataStructure/EMREvent.hh"
#include "src/common_cpp/DataStructure/MCEvent.hh"
#include "src/common_cpp/DataStructure/VirtualHit.hh"
#include "src/common_cpp/DataStructure/Primary.hh"
#include "src/common_cpp/DataStructure/Spill.hh"
#include "src/common_cpp/DataStructure/Data.hh"

// Borrowed from Chris Rogers script
#include "src/legacy/BeamTools/BTFieldConstructor.hh"
#include "src/common_cpp/Utils/JsonWrapper.hh"
#include "src/common_cpp/Utils/Globals.hh"
#include "src/common_cpp/Globals/GlobalsManager.hh"
#include "src/common_cpp/Recon/Global/GlobalTools.hh"
// #include "src/common_cpp/JsonCppStreamer/IRStream.hh"


// And my own home brew object collection.
#include "Collection.h"

  
class MCSAnalysis {
 public:
  MCSAnalysis(std::string tree, std::string mctree, std::string outname, std::map<std::string, double> histlimits);
  ~MCSAnalysis();

  void Write();

  TChain* GetTree(){ return chain; }
  TChain* GetMCTree(){ return mcchain; }
  void Execute(int mode);
  void dataSelection();
  void referenceSelection();
  void generateMCSResponse();
  void ConvolveWithInputDistribution(std::string distname);
  void DoUnfolding();
  void DoDeconvolution(std::string model, int n_sel);
  void DoFFTDeconvolution();
  void UpdateRunInfo();
  void PlotRunInfo();
  void FitGaussian(std::string outfilename);
  void CalculateChi2(std::string outfilename, std::string distname);
  double myfunc(double pz, double s1, double s2, double E, double delta);
  double myfunc_deriv(double pz, double s1, double s2, double E, double delta);
  double myfunc1(double x);
  double myfunc_deriv1(double x);
  double operator()(double x) const { return x*x; } 
  double Eval(double x) const { return x+x; }
  double Derivative(double x) const { return 2*x; }

  void SetTOFUpperLimit(double a){ TOF_upper_limit=a; }
  void SetTOFLowerLimit(double a){ TOF_lower_limit=a; }
  void SetRadialLimit(double a){ meanp=a; }
  void SetGradientLimit(double a){ sigmap=a; }
  void SetModelFileName(std::string a){ modelfile=a; }
  void SetModelFileName2(std::string a){ modelfile2=a; }
  void SetModelName1(std::string a) {modelname1=a; }
  void SetModelName2(std::string a) {modelname2=a; }
  void SetModelName3(std::string a) {modelname3=a; }
  void SetParentGeometryFile(std::string a) {geometryfile=a; }
  void SetFFTBinLimit(int a) { binlimit=a; }
  void SetFileName(std::string a) {outfilename=a; }

  void SetSysOffset(std::string eff, double val) { _sys[eff] = val; }

 private:
  
  int jUS, jDS, kUS, kDS;
  
  double TOF_lower_limit;
  double TOF_upper_limit;

  double meanp;
  double sigmap;
  int binlimit;
  
  std::string modelfile;
  std::string modelfile2;
  std::string modelname1;
  std::string modelname2;
  std::string modelname3;
  std::string geometryfile;

  double USrefplaneZ;
  double DSrefplaneZ;
  int USrefplaneI;
  int DSrefplaneI;

  Collection _USset;
  Collection _DSset;
  Collection _USAbsset;
  Collection _DSAbsset;
  Collection _USMCset;
  Collection _DSMCset;
  Collection _UStmpset;

  RooUnfoldResponse resp_thetaX;
  RooUnfoldResponse resp_thetaY;
  RooUnfoldResponse resp_thetaScatt;
  RooUnfoldResponse resp_theta2Scatt;

  RooUnfoldResponse tresp_thetaX;
  RooUnfoldResponse tresp_thetaY;
  RooUnfoldResponse tresp_thetaScatt;
  RooUnfoldResponse tresp_theta2Scatt;

  RooUnfoldResponse mresp_thetaX;
  RooUnfoldResponse mresp_thetaY;
  RooUnfoldResponse mresp_thetaScatt;
  RooUnfoldResponse mresp_theta2Scatt;
  
  RooUnfoldResponse fresp_thetaX;
  RooUnfoldResponse fresp_thetaY;
  RooUnfoldResponse fresp_thetaScatt;
  RooUnfoldResponse fresp_theta2Scatt;
  
  int USabsPlaneI;
  int DSabsPlaneI;


  
  std::map<std::string, double> _histlimits;
  std::map<std::string, double> _sys;
  // Chain containing only the data of interest (not necessarily MC).
  TChain* chain;
  // Training tree containing the response to the volume of interest
  TChain* mcchain; 
  
  int runnumber, LastRunNumber; 
  int SpillNumber;
  MAUS::TOFEvent* tofevent;
  MAUS::SciFiEvent* scifievent;
  MAUS::KLEvent* klevent;
  MAUS::CkovEvent* ckovevent;
  MAUS::EMREvent* emrevent;
  MAUS::MCEvent* mcevent;
  MAUS::Primary* primary;
  MAUS::SpecialVirtualHitArray* sphitarray;

  TFile* outfile;
  std::string outfilename;
  TH1D* tof10;
  TH1D* tof10_sel;
  TH1D* tof21;
  TH1D* tof21_sel;
  TH1D* calc_mom;
  TH2D* t_cor;
  TH2D* residual_plot;
  TH2D* TOFcom;
  TH1D* cor_mom;
  TH1D* mctof10;
  TH1D* mctof10_sel;
  TH1D* mctof21;
  TH1D* mctof21_sel;
  TH1D* mccalc_mom;
  TH1D* mctrue_mom;
  TH1D* cuts_accept;
  TH1D* mccuts_accept;

  TH1D* scattering_pos_x;
  TH1D* scattering_pos_y;
  TH1D* theta_true_x_graph;
  TH1D* theta_true_x_bin;
  TH1D* theta_true_y_graph;

  TH2D* scattering_proj_x_R;
  TH2D* scattering_proj_y_R;
  
  TH2D* scattering_proj_x_resp;
  TH2D* scattering_proj_y_resp;

  TH1D* theta_meas_y_graph;
  TH1D* theta_meas_x_graph;

  std::vector<int> _run;
  std::vector<double> _emit;
  std::vector<double> _emitError;
  std::vector<double> _meanX;
  std::vector<double> _meanXerror;
  std::vector<double> _meanY;
  std::vector<double> _meanYerror;
  std::vector<double> _meandXdz;
  std::vector<double> _meandXdzError;
  std::vector<double> _meandYdz;
  std::vector<double> _meandYdzError;
  std::vector<double> path_length;
  
  bool MatchUSDS();
  bool PIDSelection(bool isdata);
  bool RadialSelection(double pz);
  std::vector<double> DefineProjectionAngles(Vars US, Vars DS);
  double MomentumFromTOF(bool isdata);
  double BetheBloch(double pz);
  double TimeofFlight();
  std::vector<double> CalculatePathLength(double pz);
  double CorMomFromTOF(double pz);
  bool findVirtualPlanes();
  void FillMuScattResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC);
  void FillMCSResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC);
  void FillVarsVirtual(Vars& tmpvar, int j);
  void FillCollectionSciFi(Collection& Set, int j, int k, double pz, int isDS, bool project=false);
  void FillVarsSciFi(Vars& tmpvar, int j, int k, double pz, int isDS);
  void make_beam_histograms(Collection Set, std::string desc, std::string suffix);
  void make_acceptance_histograms(Collection USset, Collection DSset, 
				  std::string desc, std::string suffix);
  void make_scattering_acceptance_histograms(Collection USset,
					     Collection DSset,
					     Collection DSrefset,
					     std::string desc,
					     std::string suffix);
  Json::Value SetupConfig(int verbose_level);
  Vars PropagateVarsMu(Vars event, double z0);
  /*
  static TVectorD p_vec;
  static TVectorD res;
  static TVectorD pStart_vec;
  static TVectorD pStart_vec_y;
  static TVectorD theta_true_x;
  static TVectorD theta_true_y;
  */
  double pStart[19];
  // double p_vec[19];
  
  int counter;
  int nSize_true;
  int nSize_true_y;
  double scattering_proj_x_entry;
  double scattering_proj_x_new;
  double scattering_proj_x_centre;
  double scattering_proj_y_entry;
  double scattering_proj_y_new;
  double scattering_proj_y_centre;
  
  struct SumDistance2 {
    TVectorD D_x_func;
    TMatrixD e_x_func;
    TMatrixD R_x_func;

  SumDistance2(TVectorD D, TMatrixD e, TMatrixD R) :
    D_x_func(D), e_x_func(e), R_x_func(R) {}
    double operator() (const double * par) {
      double sum = 0;
      TVectorD p_vec(19);
      int z = 1;
      for (int j=0; j<19; j++)
	p_vec[j] = par[j];
      for (int i=10; i<19; i++){
	p_vec[i] = par[i-z];
	z += 2;
      }
      TVectorD res = R_x_func * e_x_func * p_vec;
      for (int i = 0; i<30; i++) {
	if (D_x_func(i) != 0) {
	  double d = pow(D_x_func(i)-res(i),2) / D_x_func(i);
	  sum += d;
	}
	/*
	  if (first) {
	  std::cout << "point " << i << "\t"
	  << D_x_func(i) << "\t"
	  << res(i) << "\t" << std::endl;
	  }
	  }
	  if (first) {
	  std::cout << "Total Initial distance square = " << sum << std::endl;
	  first = false;
	*/
      }
      return sum;
   }
  };




	     
};
    
class MyFunction1D { 
	public:
		double E;
		double delta;
		double s1;
		double s2;
		double operator()(double pz) const { return s1*sqrt(pow(pz+delta,2)+pow(105.65,2))/(pz+delta)+s2*sqrt(pow(pz-delta,2)+pow(105.65,2))/(pz-delta); } 
		double Derivative(double pz) const { return -(s1 - s2) * delta * pow(105.65,2)/ (pow(pz,2)*E); }
};
 
