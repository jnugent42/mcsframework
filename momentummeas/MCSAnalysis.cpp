#include "MCSAnalysis.h"

// #include "TVirtualFFT.h"

MCSAnalysis::MCSAnalysis(std::string tree, std::string mctree, std::string outname, std::map<std::string, double> histlimits)
//: p_vec(TVectorD(19)), res(TVectorD(30)), pStart_vec(TVectorD(19)),
//    pStart_vec_y(TVectorD(19)), theta_true_x(TVectorD(19)), theta_true_y(TVectorD(19))
{
  chain   = new TChain(tree.c_str());
  mcchain = new TChain(mctree.c_str());
  outfile = new TFile(outname.c_str(), "recreate");

  _sys["alXUS"] = 0.0; //-0.0356;
  _sys["alYUS"] = 0.0; //1.13;
  _sys["alXDS"] = 0.0; //-0.126;
  _sys["alYDS"] = 0.0; //0.118;
  _sys["thXUS"] = 0.0; //-0.0031;
  _sys["thYUS"] = 0.0; //-0.005;
  _sys["thXDS"] = 0.0; //-0.003;
  _sys["thYDS"] = 0.0; //-0.0151;
  _sys["resX"] = 1.0;
  _sys["resY"] = 1.0;
  _sys["psel_lcut"] = -1.;
  _sys["psel_ucut"] = -1.;
  _sys["TOF0_z"]  = 5287.24720607;
  _sys["TOF1_z"]  = 12929.2608098;
  _sys["TOF2_z"]  = 21138.306769;
  _sys["niter"]   = 10;
  _sys["abspos"]  = 16952.5;
  _sys["Nevents"]  = -1.;
  _sys["FracEvents"] = 1.;


  modelname2 = "Cobb";
  modelname1 = "GEANT";
  modelname3 = "Moliere";
 
  _histlimits["NbinsXY"] = histlimits.count("NbinsXY") != 0 ? histlimits["NbinsXY"]: 200;
  _histlimits["minXY"] = histlimits.count("minXY") != 0 ? histlimits["minXY"]: -0.2;
  _histlimits["maxXY"] = histlimits.count("maxXY") != 0 ? histlimits["maxXY"]:  0.2;
  _histlimits["NbinsTh"] = histlimits.count("NbinsTh") != 0 ? histlimits["NbinsTh"]: 100;
  _histlimits["minTh"] = histlimits.count("minTh") != 0 ? histlimits["minTh"]:  0.0;
  _histlimits["maxTh"] = histlimits.count("maxTh") != 0 ? histlimits["maxTh"]:  0.2;
  _histlimits["NbinsTh2"] = histlimits.count("NbinsTh2") != 0 ? histlimits["NbinsTh2"]: 400;
  _histlimits["minTh2"] = histlimits.count("minTh2") != 0 ? histlimits["minTh2"]:  0.0;
  _histlimits["maxTh2"] = histlimits.count("maxTh2") != 0 ? histlimits["maxTh2"]:  0.004;
  
  resp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  resp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  resp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  resp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  tresp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tresp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tresp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tresp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  mresp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  mresp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  mresp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  mresp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  
  fresp_thetaX.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  fresp_thetaY.Setup(_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  fresp_thetaScatt.Setup(_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  fresp_theta2Scatt.Setup(_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  
  const int NBINS = 19;
  const int NPBINS = 30;
  double scat_bin_array[NBINS+1] =
    {-0.1151, -0.0938, -0.0754, -0.0597,  -0.0463,
     -0.0347, -0.0248, -0.0162, -0.00895, -0.00269,
     0.00269,  0.00895, 0.0162,  0.0248,   0.0347,
     0.0463,   0.0597,  0.0754,  0.0938,   0.1151};
  double pos_array[NPBINS+1];
  for(int i=0; i<NPBINS+1; i++)  pos_array[i] = -100 + (6.67 * i);

  runnumber  = 0;
  LastRunNumber = 0;
  tofevent   = new MAUS::TOFEvent();
  scifievent = new MAUS::SciFiEvent();
  klevent    = new MAUS::KLEvent();
  ckovevent  = new MAUS::CkovEvent();
  emrevent   = new MAUS::EMREvent();
  mcevent    = new MAUS::MCEvent();
  primary    = new MAUS::Primary();
  sphitarray = new MAUS::SpecialVirtualHitArray();

  // set some default values that may be adjusted later.
  TOF_lower_limit = 27;
  TOF_upper_limit = 28;
  meanp = 300.0;
  sigmap = 0.035;
  binlimit = 22;
  int counter = 0;
  // histograms of selection variables
  
  t_cor = new TH2D("t_cor","t_cor", 1000, 0, 1000,100, 8.05, 9.15);
  residual_plot = new TH2D("MCTruth vs p_corrected","MCTruth vs p_corrected", 100, 200, 300,100, 200, 300);
  TOFcom = new TH2D("TOF01 vs TOF12","TOF01 vs TOF12", 200, 210, 222,200, 210, 222);
  tof10 = new TH1D("tof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  tof10_sel = new TH1D("tof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  tof21 = new TH1D("tof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  tof21_sel = new TH1D("tof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  mctof10 = new TH1D("mctof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  mctof10_sel = new TH1D("mctof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
  mctof21 = new TH1D("mctof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  mctof21_sel = new TH1D("mctof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
  calc_mom = new TH1D("calc_mom","Momentum Calculated from TOF; Momentum (MeV/c)", 400, 0, 400);
  cor_mom = new TH1D("cor_mom","Cor Momentum Calculated from TOF; Momentum (MeV/c)", 400, 0, 400);
  mccalc_mom = new TH1D("mccalc_mom","Momentum Calculated from TOF; Momentum (MeV/c)", 400, 0, 400);
  mctrue_mom = new TH1D("mctrue_mom","Momentum from Virtual Planes; Momentum (MeV/c)", 400, 0, 400);
  cuts_accept   = new TH1D("cuts_accept", 
			   ";Selection Criteria; Surviving Events", 5, 0, 5);
  mccuts_accept = new TH1D("mccuts_accept", 
			   ";Selection Criteria; Surviving Events", 6, 0, 6);

  scattering_pos_x =
    new TH1D("scattering_pos_x",
	     "Change in Position (X);#Delta X; Events per mm", 
	     NBINS, scat_bin_array);
  scattering_pos_y =
    new TH1D("scattering_pos_y",
	     "Change in Position (Y);#Delta Y; Events per mm", 
	     NBINS, scat_bin_array);

  theta_true_x_graph =
    new TH1D("theta_true_x_graph",
	     "Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
	     NBINS, scat_bin_array);
  theta_true_y_graph =
    new TH1D("theta_true_y_graph",
	     "Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
	     NBINS, scat_bin_array);
  
  scattering_proj_x_resp =
    new TH2D("scattering_proj_x_resp",
	     "Change in Projected Angle (X);#Delta#theta^{Measured}_{X}; #Delta#theta^{True}_{X}", 
	     NBINS, scat_bin_array, NBINS, scat_bin_array);
  scattering_proj_y_resp =
    new TH2D("scattering_proj_y_resp",
	     "Change in Projected Angle (Y);#Delta#theta^{Measured}_{Y}; #Delta#theta^{True}_{Y}", 
	     NBINS, scat_bin_array, NBINS, scat_bin_array);
  
  scattering_proj_x_R =
    new TH2D("scattering_proj_x_R",
	     "Change in Projected Angle (X);#Delta#theta^{Measured}_{X}; #Delta#theta^{True}_{Y}", 
	     NPBINS, pos_array, NBINS, scat_bin_array);
  scattering_proj_y_R =
    new TH2D("scattering_proj_y_R",
	     "Change in Projected Angle (Y);#Delta#theta^{Measured}_{Y}; #Delta#theta^{True}_{X}", 
	     NPBINS, pos_array, NBINS, scat_bin_array);
  theta_meas_x_graph =
    new TH1D("theta_meas_x_graph","Accepted Events;#Delta#theta^{True}_X",
	     NBINS, scat_bin_array);
  theta_meas_y_graph =
    new TH1D("theta_meas_y_graph","Accepted Events;#Delta#theta^{True}_Y",
	     NBINS, scat_bin_array);
  jUS=-1;
  jDS=-1;
  kUS=-1;
  kDS=-1;
  
  USrefplaneI = -1;
  DSrefplaneI = -1;
  USrefplaneZ = 0.;
  DSrefplaneZ = 0.;
  USabsPlaneI = -1;
  DSabsPlaneI = -1;
}

MCSAnalysis::~MCSAnalysis(){
  delete chain;
  delete mcchain;
  // delete tof10;
  // delete tof21;
  // delete tof10_sel;
  // delete tof21_sel;
  // delete calc_mom;
}

void MCSAnalysis::Write(){
  
  outfile->cd();
  tof10->Write();
  tof21->Write();
  tof10_sel->Write();
  tof21_sel->Write();
  calc_mom->Write();
  TCanvas *c1 = new TCanvas();
  calc_mom->GetXaxis()->SetRangeUser(120,280);
  calc_mom->Draw();
  cor_mom->SetLineColor(2);
  cor_mom->Draw("SAME");
  c1->SaveAs("calc_mom.pdf");
  cor_mom->Write();
  c1->Clear();
  calc_mom->GetXaxis()->SetRangeUser(120,280);
  cor_mom->Draw();
  c1->SaveAs("cor_mom.pdf");
  c1->Clear();
  residual_plot->GetXaxis()->SetTitle("MCTruth pz");
  residual_plot->GetYaxis()->SetTitle("p_corrected");
  residual_plot->Draw();
  TLine *line = new TLine(200,200,300,300);
  line->SetLineColor(kRed);
  line->Draw();
  c1->SaveAs("residual_plot.pdf");
  c1->Clear();
  TOFcom->GetXaxis()->SetTitle("p TOF01");
  TOFcom->GetYaxis()->SetTitle("p TOF12");
  TOFcom->Draw();
  c1->SaveAs("TOFcom.root");
  c1->SaveAs("TOFcom.pdf");
  mctof10->Write();
  mctof21->Write();
  mctof10_sel->Write();
  mctof21_sel->Write();
  mccalc_mom->Write();
  mctrue_mom->Write();
  cuts_accept->Write();
  mccuts_accept->Write();
  //outfile->Close();
  delete c1;
}

void MCSAnalysis::Execute(int mode=0){
  
  int verbose_level = 1;
  // MAUS::GlobalsManager::InitialiseGlobals(JsonWrapper::JsonToString(SetupConfig(verbose_level)));

  dataSelection();
  
  if (mode == 0){
    generateMCSResponse();
    // DoUnfolding();
    DoDeconvolution(modelname2.c_str(), 1);
  }
  else if (mode >= 1){
    referenceSelection();
    DoUnfolding();
    DoFFTDeconvolution();
    ConvolveWithInputDistribution(modelname1.c_str());
    DoDeconvolution(modelname1.c_str(), 1);
    ConvolveWithInputDistribution(modelname2.c_str());
    DoDeconvolution(modelname2.c_str(), 1);
    //FitGaussian(outfilename.c_str());
    //CalculateChi2(outfilename.c_str(), modelname2.c_str());
    //Write();
  }
  else if (mode == -1){
    referenceSelection();
    DoUnfolding();
    DoFFTDeconvolution();
    ConvolveWithInputDistribution(modelname1.c_str());
    DoDeconvolution(modelname1.c_str(), 1);
    ConvolveWithInputDistribution(modelname2.c_str());
    DoDeconvolution(modelname2.c_str(), 1);
    ConvolveWithInputDistribution(modelname3.c_str());
    DoDeconvolution(modelname3.c_str(), 1);
    //FitGaussian(outfilename.c_str());
    //CalculateChi2(outfilename.c_str(), modelname2.c_str());
    //Write();
  }
  /* else if (mode > 1){
    referenceSelection();
    ConvolveWithInputDistribution(modelname1.c_str());
    DoDeconvolution(modelname1.c_str(), mode);
    } */
  else {
    std::cout<<"Unknown operation mode"<<std::endl;
  } 
}

void MCSAnalysis::dataSelection(){
  // Set addresses for tree selection
  chain->SetBranchAddress("RunNumber", &runnumber);
  chain->SetBranchAddress("TOFBranch", &tofevent);
  chain->SetBranchAddress("SciFiBranch", &scifievent);
  chain->SetBranchAddress("CkovBranch", &ckovevent);
  chain->SetBranchAddress("KLBranch", &klevent);
  chain->SetBranchAddress("EMRBranch", &emrevent);
  chain->SetBranchAddress("MCEvent", &mcevent);

  // Restrict the number of entries to less than or equal to the mcchain entries
  
  int Nentries = chain->GetEntries();
  // < mcchain->GetEntries() ? 
  //  chain->GetEntries() : mcchain->GetEntries();
  
  // Loop over all tree entries.
  Collection USAllTOF, DSAllTOF, USPreRadSel, DSPreRadSel;
  double pz = 0.;
  int Nevents = 0;
  // Nentries *= int(_sys["FracEvents"]) > 0 ? _sys["FracEvents"] : 1.0;
  for (int i=0; i<100000; i++){
    chain->GetEntry(i);
    if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 
    // Set cuts based on the TOFs, ckov, kl, and EMR information
    // Locate the tracker reference planes. To be agnostic locate
    // the downstream most station of the upstream tracker and the
    // upstream most station of the downstream tracker.
    cuts_accept->Fill("All Events",1);
    // if (scifievent->scifitracks().size() != 2) continue;
    // if ( !SelectMomentum() ) continue;
    if ( !MatchUSDS() ) {
      if (jUS == -1 || kUS == -1){
	continue;
      } 
    }
    cuts_accept->Fill("US Track Found",1);    
    FillCollectionSciFi(USAllTOF, jUS, kUS, pz, 0);
    FillCollectionSciFi(DSAllTOF, jDS, kDS, pz, 1);
    if ( !PIDSelection(true) ) continue;
    pz = MomentumFromTOF(true);
    cuts_accept->Fill("TOF Selection",1);
    FillCollectionSciFi(USPreRadSel, jUS, kUS, pz, 0);
    FillCollectionSciFi(DSPreRadSel, jDS, kDS, pz, 1);
    if ( !RadialSelection(pz) ) continue;
    cuts_accept->Fill("Fiducial Selection",1);
    if (jDS != -1 && kDS != -1) { 
        double pz_cor = CorMomFromTOF(pz);
    }
    if ( _sys["psel_lcut"] > 0 && _sys["psel_ucut"] > 0 ){
      if ( pz < _sys["psel_lcut"] || pz > _sys["psel_ucut"] ) 
	// Cut out events to evaluate systematic uncertainties.
	continue;
    }
    Nevents++;
    //
    if(_sys["Nevents"] > 0 && Nevents > int(_sys["Nevents"]))
      break;
    FillCollectionSciFi(_USset, jUS, kUS, pz, 0);
    FillCollectionSciFi(_DSset, jDS, kDS, pz, 1);
    FillCollectionSciFi(_UStmpset, jUS, kUS, pz, 0);
    if (runnumber != LastRunNumber){
      UpdateRunInfo();
    }
  }
  make_beam_histograms(USAllTOF, "Upstream, Data", "dataUS_alltof");
  make_beam_histograms(DSAllTOF, "Downstream, Data", "dataDS_alltof");
  make_acceptance_histograms(USAllTOF, DSAllTOF, 
			     "Data Projection", "dataAcc_alltof");
  make_beam_histograms(USPreRadSel, "Upstream, Data", "dataUS_prerad");
  make_beam_histograms(DSPreRadSel, "Downstream, Data", "dataDS_prerad");
  make_beam_histograms(_USset, "Upstream, Data", "dataUS");
  make_beam_histograms(_DSset, "Downstream, Data", "dataDS");
  make_acceptance_histograms(_USset, _DSset, "Data Projection", "dataAcc");
  PlotRunInfo();
}
void MCSAnalysis::referenceSelection(){
  // Set addresses for tree selection
  mcchain->SetBranchAddress("TOFBranch", &tofevent);
  mcchain->SetBranchAddress("SciFiBranch", &scifievent);
  mcchain->SetBranchAddress("CkovBranch", &ckovevent);
  mcchain->SetBranchAddress("KLBranch", &klevent);
  mcchain->SetBranchAddress("EMRBranch", &emrevent);

  // Restrict the number of entries to less than or equal to the mcchain entries
  
  int Nentries = mcchain->GetEntries();
  // < mcchain->GetEntries() ? 
  //  chain->GetEntries() : mcchain->GetEntries();
  Collection USAllTOF, DSAllTOF, USPreRadSel, DSPreRadSel;
  double pz = 0.;
  // Loop over all tree entries.
  for (int i=0; i<100000; i++){
    mcchain->GetEntry(i);
    if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 
    // Set cuts based on the TOFs, ckov, kl, and EMR information
    
    // Locate the tracker reference planes. To be agnostic locate
    // the downstream most station of the upstream tracker and the
    //// upstream most station of the downstream tracker.
    cuts_accept->Fill("All Events",1);
    if ( !MatchUSDS() ) {
      if (jUS == -1 || kUS == -1){
	continue;
      } 
    }
    cuts_accept->Fill("US Track Found",1);
    FillCollectionSciFi(USAllTOF, jUS, kUS, pz, 0);
    FillCollectionSciFi(DSAllTOF, jDS, kDS, pz, 1);
    if ( !PIDSelection(false) ) continue;
    cuts_accept->Fill("TOF Selection",1);
    pz = MomentumFromTOF(false);
    // if (scifievent->scifitracks().size() != 2) continue;
    // if ( !SelectMomentum() ) continue;
    FillCollectionSciFi(USPreRadSel, jUS, kUS, pz, 0);
    FillCollectionSciFi(DSPreRadSel, jDS, kDS, pz, 1);
    if ( !RadialSelection(pz) ) continue;
    cuts_accept->Fill("Fiducial Selection",1);
    
    FillCollectionSciFi(_USMCset, jUS, kUS, pz, 0);
    FillCollectionSciFi(_DSMCset, jDS, kDS, pz, 1);
  }
  make_beam_histograms(USAllTOF, "Upstream, Data Reference", "dataUSref_alltof");
  make_beam_histograms(DSAllTOF, "Downstream, Data Reference", "dataDSref_alltof");
  make_acceptance_histograms(USAllTOF, DSAllTOF, "Data Reference", "dataProj_alltof");
  make_beam_histograms(USPreRadSel, "Upstream, Data Reference", "dataUSref_prerad");
  make_beam_histograms(DSPreRadSel, "Downstream, Data Reference", "dataDSref_prerad");
  make_acceptance_histograms(USPreRadSel, DSPreRadSel, "Data Reference", "dataProj_prerad");
  make_beam_histograms(_USMCset, "Upstream, Data Reference", "dataUSref");
  make_beam_histograms(_DSMCset, "Downstream, Data Reference", "dataDSref");
  make_acceptance_histograms(_USMCset, _DSMCset, "Data Reference", "dataProj");
}

 
void MCSAnalysis::generateMCSResponse(){

  mcchain->SetBranchAddress("TOFBranch", &tofevent);
  mcchain->SetBranchAddress("SciFiBranch", &scifievent);
  mcchain->SetBranchAddress("CkovBranch", &ckovevent);
  mcchain->SetBranchAddress("KLBranch", &klevent);
  mcchain->SetBranchAddress("EMRBranch", &emrevent);
  mcchain->SetBranchAddress("MCEvent", &mcevent);
  // Loop over the training sample.
  Collection USRecSet, DSRecSet, USVirtSet, DSVirtSet;
  Collection USAllTOF, DSAllTOF, USPreRadSel, DSPreRadSel;
  Collection USVAllTOF, DSVAllTOF, USVPreRadSel, DSVPreRadSel;
  // Collection USeRecSet, DSeRecSet, USmuRecSet, DSmuRecSet;
  int nMuAbsSel=0, nEAbsSel=0, nMuAbsAll=0, nEAbsAll=0;
  TH1D* mom_resUS = new TH1D("mom_resUS", ";p_{rec} - p_{MC} (MeV/c)",
			     2000, -100, 100);
  TH1D* mom_resDS = new TH1D("mom_resDS", ";p_{rec} - p_{MC} (MeV/c)",
			     2000, -100, 100);
  TH2D* mom_responseUS = new TH2D("mom_responseUS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)",
				  300, 0, 300, 300, 0, 300);
  TH2D* mom_responseDS = new TH2D("mom_responseDS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)",
				  300, 100, 300, 300, 0, 300);
  TH2D* mom_responseABS = new TH2D("mom_responseABS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)",
				  300, 100, 300, 300, 0, 300);
  int ngood=0;
  for (int j=0; j<mcchain->GetEntries(); j++){
    mcchain->GetEntry(j);
    if (j%100000==0) std::cout<<"MC Event "<<j<<", selected "<<ngood<<" events.\n"; 
    // if(fabs(mcevent->GetPrimary()->GetParticleId()) != 13) continue;
    // Select events that produce virtual plane hits and pass the data selection cuts.
    Vars USAbsHit, DSAbsHit, USTrackerRefHit, DSTrackerRefHit;
    double pz = 0.;
    mccuts_accept->Fill("All Events",1);
    if ( !PIDSelection(false) ) continue; // event_ok=false;
    if( !findVirtualPlanes() ) continue;
    mccuts_accept->Fill("Found Virtual Planes",1);
    FillVarsVirtual(USAbsHit, USabsPlaneI);
    FillVarsVirtual(DSAbsHit, DSabsPlaneI);
    FillVarsVirtual(USTrackerRefHit, USrefplaneI);
    FillVarsVirtual(DSTrackerRefHit, DSrefplaneI);
    if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13 &&
	mcevent->GetVirtualHits()->at(DSrefplaneI).GetParticleId()==-11)
      nEAbsAll++;
    if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13)
      nMuAbsAll++;
    // Apply selection cuts as with the data
    bool event_ok=true;
    if ( !MatchUSDS() ) {
      if (jUS == -1 || kUS == -1){
	cuts_accept->Fill("US Tracker Found", 1);
	event_ok=false;
      } 
    }
    if (event_ok) mccuts_accept->Fill("US Track Found",1);
    if (event_ok) pz = MomentumFromTOF(false);
    if (event_ok) mccuts_accept->Fill("TOF Selection",1);
    // if (scifievent->scifitracks().size() != 2) continue;
    FillCollectionSciFi(USPreRadSel, jUS, kUS, pz, 0);
    FillCollectionSciFi(DSPreRadSel, jDS, kDS, pz, 1);
    USVPreRadSel.append_instance(USAbsHit);
    DSVPreRadSel.append_instance(DSAbsHit);
    if ( !RadialSelection(pz) ) event_ok=false;
    if (event_ok) mccuts_accept->Fill("Fiducial Selection",1);
    if (event_ok) mctrue_mom->Fill(mcevent->GetVirtualHits()->at(USrefplaneI).GetMomentum().z());
    if (event_ok) mctrue_mom->Fill(mcevent->GetVirtualHits()->at(USrefplaneI).GetMomentum().z());
    Vars USSciFiRec, DSSciFiRec;
    if (event_ok){
      FillVarsSciFi(USSciFiRec, jUS, kUS, pz, 0);
      FillVarsSciFi(DSSciFiRec, jDS, kDS, pz, 1);
      USRecSet.append_instance(USSciFiRec);
      DSRecSet.append_instance(DSSciFiRec);
      if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13 &&
	  mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-11)
	nEAbsSel++;
      if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13)
	nMuAbsSel++;
      
      mom_resUS->Fill(USSciFiRec.pz - USTrackerRefHit.pz);
      mom_resDS->Fill(DSSciFiRec.pz - DSTrackerRefHit.pz);
      mom_responseUS->Fill(USSciFiRec.pz, USTrackerRefHit.pz);
      mom_responseDS->Fill(DSSciFiRec.pz, DSTrackerRefHit.pz);
      mom_responseABS->Fill(DSSciFiRec.pz, USAbsHit.pz);
    }
    USVirtSet.append_instance(USAbsHit);
    DSVirtSet.append_instance(DSAbsHit);
    FillMCSResponse(event_ok, USSciFiRec, DSSciFiRec, USAbsHit, DSAbsHit);
    FillMuScattResponse(event_ok, USSciFiRec, DSSciFiRec, USAbsHit, DSAbsHit);
    ngood++;
  }
  std::cout<<"For all events simulated there are "<<nEAbsAll<<" decay electrons and "<<nMuAbsAll<<" muons.\n";
  std::cout<<"For the selected MC events there are "<<nEAbsSel<<" decay electrons and "<<nMuAbsSel<<" muons.\n";
  make_beam_histograms(USVirtSet, "Upstream, Data", "VirtMCUS");
  make_beam_histograms(DSVirtSet, "Downstream, Data", "VirtMCDS");
  //make_beam_histograms(USAllTOF, "Upstream, Reconstructed Simulation", "recMCUS_alltof");
  //make_beam_histograms(DSAllTOF, "Downstream, Reconstructed Simulation", "recMCDS_alltof");
  make_beam_histograms(USPreRadSel, "Upstream, Reconstructed Simulation", "recMCUS_prerad");
  make_beam_histograms(DSPreRadSel, "Downstream, Reconstructed Simulation", "recMCDS_prerad");
  make_scattering_acceptance_histograms(USVPreRadSel,
					DSVPreRadSel, DSPreRadSel,
					"Virtual Projection","VirtPreRad");

  make_beam_histograms(USRecSet, "Upstream, Reconstructed Simulation", "recMCUS");
  make_beam_histograms(DSRecSet, "Downstream, Reconstructed Simulation", "recMCDS");
  make_scattering_acceptance_histograms(USVirtSet,
					DSVirtSet,DSRecSet,
					"Virtual Projection","VirtProj");
  mom_resUS->Write();
  mom_resDS->Write();
  mom_responseUS->Write();
  mom_responseDS->Write();
  mom_responseABS->Write();
}

bool MCSAnalysis::MatchUSDS(){

  bool trackmatched = true;
  jUS=-1;
  jDS=-1;
  kUS=-1;
  kDS=-1;
  
  if( scifievent->scifitracks().size() == 1 || 
      scifievent->scifitracks().size() == 2){
    for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
      int npoints = scifievent->scifitracks()[j]->scifitrackpoints().size();
      double maxUS=0.0, minDS=44000;
      int tracker = scifievent->scifitracks()[j]->tracker();
      for ( int k=0; k<npoints; k++){
	double zpos = 
	  scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
	if(tracker==0 && zpos > maxUS){
	  maxUS = zpos;
	  kUS = k;
	  jUS = j;
	}
	if(tracker==1 && zpos < minDS){      
	  minDS = zpos;
	  kDS = k;
	  jDS = j;
	}
      }
      if (jUS != -1 && kUS != -1 && jDS != -1 && kDS != 1) {
	USrefplaneZ = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
	DSrefplaneZ = minDS;
	// std::cout<<"Identified US track by ["<<jUS<<", "<<kUS<<"] at z = "<<USrefplaneZ<<" mm and DS track by ["<<jDS<<", "<<kDS<<"] at z = "<<DSrefplaneZ<<" mm."<<std::endl; 
      }
    }
    if (jUS == -1 || kUS == -1) {
      //std::cout<<"Failed US track by ["<<jUS<<", "<<kUS<<"] and DS track by ["<<jDS<<", "<<kDS<<"]"<<std::endl; 
      trackmatched = false;
    }
  }
  
  else trackmatched = false;
  return trackmatched;  
}
bool MCSAnalysis::PIDSelection(bool isdata=true){

  // Consider the TOF hits first 
  double rawTOF0HitTime = -1., rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
  /*
  if ( tofevent->GetTOFEventSlabHit().GetTOF0SlabHitArraySize() == 2 ) 
    rawTOF0HitTime  = ( tofevent->GetTOFEventSlabHit().GetTOF0SlabHitArray()[0].GetRawTime() +
			tofevent->GetTOFEventSlabHit().GetTOF0SlabHitArray()[1].GetRawTime() ) / 2.;
  else
    return false;
  if ( tofevent->GetTOFEventSlabHit().GetTOF1SlabHitArraySize() == 2 ) 
    rawTOF1HitTime  = ( tofevent->GetTOFEventSlabHit().GetTOF1SlabHitArray()[0].GetRawTime() +
			tofevent->GetTOFEventSlabHit().GetTOF1SlabHitArray()[1].GetRawTime() ) / 2.;
  else
    return false;
  if ( tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArraySize() == 2 ) 
    rawTOF2HitTime  = ( tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArray()[0].GetRawTime() +
			tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArray()[1].GetRawTime() ) / 2.;
  else
    return false;
  */
  if ( tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArraySize() >= 1 )
    rawTOF0HitTime = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  else
    return false;
  if ( tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArraySize() >= 1 )
    rawTOF1HitTime = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  else
    return false;
  if (  tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArraySize() >= 1 )
    rawTOF2HitTime = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
  else // Allow for TOF2 to not be hit
    rawTOF2HitTime = rawTOF1HitTime + 100.0 * 8.22475 / 0.299792458; // ns. 
  // return false;
  if(isdata){
    tof10->Fill(rawTOF1HitTime - rawTOF0HitTime);
    tof21->Fill(rawTOF2HitTime - rawTOF1HitTime);
    if ( rawTOF1HitTime - rawTOF0HitTime < TOF_lower_limit ||
	 rawTOF1HitTime - rawTOF0HitTime > TOF_upper_limit) return false;
    // if ( rawTOF2HitTime - rawTOF1HitTime < TOF_lower_limit ||
    //      rawTOF2HitTime - rawTOF1HitTime > TOF_upper_limit) return false;
    tof10_sel->Fill(rawTOF1HitTime - rawTOF0HitTime);
    tof21_sel->Fill(rawTOF2HitTime - rawTOF1HitTime);
  }
  else {
    mctof10->Fill(rawTOF1HitTime - rawTOF0HitTime);
    mctof21->Fill(rawTOF2HitTime - rawTOF1HitTime);
    if ( rawTOF1HitTime - rawTOF0HitTime < TOF_lower_limit ||
	 rawTOF1HitTime - rawTOF0HitTime > TOF_upper_limit) return false;
    // if ( rawTOF2HitTime - rawTOF1HitTime < TOF_lower_limit ||
    //      rawTOF2HitTime - rawTOF1HitTime > TOF_upper_limit) return false;
    mctof10_sel->Fill(rawTOF1HitTime - rawTOF0HitTime);
    mctof21_sel->Fill(rawTOF2HitTime - rawTOF1HitTime);
  }
  return true;
  
}

bool MCSAnalysis::RadialSelection(double pz){
  bool selected = true;
  if (jUS == -1 || kUS == -1) 
    selected = false;
  else {
    Vars USplane;
    USplane.X = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().x() + _sys["alXUS"];
    USplane.Y = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().y() + _sys["alYUS"];
    USplane.Z = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
    // if ( sqrt(xpos*xpos + ypos*ypos) > meanp) selected = false;
    USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0);
    USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
      tan(_sys["thYUS"] * atan(1.)/45.0);
    USplane.pz   = pz; // scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().z();
    USplane.px   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    USplane.py   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    //      scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().z();
    // if ( sqrt(dXdz*dXdz + dYdz*dYdz) > sigmap) selected = false;
    double abspos = _sys["abspos"];
    double phi = atan2(USplane.dYdz, USplane.dXdz);
    double zdiff = 2*fabs(USplane.Z - abspos);
    // double xabs  = dXdz > 0 ? (dXdz + sigmap*cos(phi)) * zdiff + xpos : (dXdz - sigmap*cos(phi)) * zdiff + xpos;
    // double yabs  = dYdz > 0 ? (dYdz + sigmap*sin(phi)) * zdiff + ypos : (dYdz - sigmap*sin(phi)) * zdiff + ypos;
    USplane.dXdz += sigmap*cos(phi);
    USplane.dYdz += sigmap*sin(phi);
    USplane.px   += sigmap*cos(phi)*USplane.pz;
    USplane.py   += sigmap*sin(phi)*USplane.pz;
    // double xabs  = (dXdz + sigmap*cos(phi)) * zdiff + xpos;
    // double yabs  = (dYdz + sigmap*sin(phi)) * zdiff + ypos;
    Vars DSproj = PropagateVarsMu(USplane, _sys["abspos"] + 549.95);
    if ( sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y) > meanp) selected = false;
  }
  return selected;
}

std::vector<double> MCSAnalysis::DefineProjectionAngles(Vars US, Vars DS){

  std::vector<double> projTheta;
  double USnorm = 1./sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz);
  double u[3] = {US.dXdz*USnorm, US.dYdz*USnorm, USnorm};
  double w[3] = {-u[0]*u[1], (u[0]*u[0] + u[2]*u[2]), -u[1]*u[2]};
  double Wnorm  = 1./sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  w[0] *= Wnorm;
  w[1] *= Wnorm;
  w[2] *= Wnorm;
  double DSnorm = 1./sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz);
  double d[3] = {DS.dXdz*DSnorm, DS.dYdz*DSnorm, DSnorm};
  projTheta.push_back( atan( (d[0]*w[0] + d[1]*w[1] + d[2]*w[2])/
			     (d[0]*u[0] + d[1]*u[1] + d[2]*u[2]) ));
  projTheta.push_back( atan( (d[0]*u[2] - u[0]*d[2])\
			     /(d[0]*u[0] + d[1]*u[1] + d[2]*d[2]) * 
			     1./sqrt(u[2]*u[2] + u[0]*u[0])) );
  // projTheta.push_back( atan( (d[0]*u[2] - u[0]*d[2])			\
  ///(d[0]*u[0] + d[1]*u[1] + d[2]*u[2]) * 
  //			     1./sqrt(u[0]*u[0] + u[2]*u[2])) );

  projTheta.push_back( acos( ( (1 + US.dXdz * DS.dXdz + US.dYdz * DS.dYdz )/
			       sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz)/
			       sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz))) );
 
  return projTheta;
}

std::vector<double> MCSAnalysis::CalculatePathLength(double pz){
        
    std::vector<double> path_length;
    std::vector<double> vpath_length;

    // Path length for muon between TOF1 and absorber
    Vars USplane;
    /*
    std::cout << "jUS " << jUS << std::endl;
    std::cout << "kUS " << kUS << std::endl;
    */
    USplane.X = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().x() + _sys["alXUS"];
    USplane.Y = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().y() + _sys["alYUS"];
    USplane.Z = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
    USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0);
    USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
      tan(_sys["thYUS"] * atan(1.)/45.0);
    USplane.pz   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().z();
    USplane.px   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().x() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    USplane.py   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
      tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
    Vars USabsoproj = PropagateVarsMu(USplane, _sys["abspos"]);
    Vars USTOF1proj = PropagateVarsMu(USplane, _sys["TOF1_z"]);
    /*
    std::cout << USabsoproj.X << std::endl;
    std::cout << USplane.X << std::endl;
    std::cout << USTOF1proj.X << std::endl;
    */
    vpath_length.push_back(USabsoproj.X - USTOF1proj.X);
    vpath_length.push_back(USabsoproj.Y - USTOF1proj.Y);
    vpath_length.push_back(USabsoproj.Z - USTOF1proj.Z);
    path_length.push_back(sqrt(vpath_length[0]*vpath_length[0] + vpath_length[1]*vpath_length[1] + vpath_length[2]*vpath_length[2]));

    //Path length for muon between absorber and TOF2
    vpath_length.clear();
    Vars DSplane;
    /*
    std::cout << "jDS " << jDS << std::endl;
    std::cout << "kDS " << kDS << std::endl;
    */
    DSplane.X = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().x() + _sys["alXDS"];
    DSplane.Y = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().y() + _sys["alYDS"];
    DSplane.Z = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().z();
    DSplane.dXdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().x() +
      tan(_sys["thXDS"] * atan(1.)/45.0);
    DSplane.dYdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().y() +
      tan(_sys["thYDS"] * atan(1.)/45.0);
    DSplane.pz   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().z();
    DSplane.px   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().x() +
      tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
    DSplane.py   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().y() +
      tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
    Vars DSabsoproj = PropagateVarsMu(DSplane, _sys["abspos"]);
    Vars DSTOF2proj = PropagateVarsMu(DSplane, _sys["TOF2_z"]);
    vpath_length.push_back(DSTOF2proj.X - DSabsoproj.X);
    vpath_length.push_back(DSTOF2proj.Y - DSabsoproj.Y);
    vpath_length.push_back(DSTOF2proj.Z - DSabsoproj.Z);
    path_length.push_back(sqrt(vpath_length[0]*vpath_length[0] + vpath_length[1]*vpath_length[1] + vpath_length[2]*vpath_length[2]));

    return path_length;

}

double MCSAnalysis::CorMomFromTOF(double pz){

   // Initialise and collect initial TOF and pz
   double t_path = 0;
   double p_corrected = 0;
   double t_initial = TimeofFlight();
   if (t_initial != 100) {
   t_initial = t_initial*0.299792458;
   std::vector<double> path_length = CalculatePathLength(pz);
   if (path_length.size() == 2) {
   double s1 = path_length.at(0)/1000;
   double s2 = path_length.at(1)/1000;
   pz = 105.65*(s1+s2)/sqrt(pow(t_initial,2)-pow((s1+s2),2));
  
   // delta - from Bethe-Bloch (units are cm)
   double z = 3.25;
   double rho = 0.694;
   double dEdx = BetheBloch(pz);
   double Eloss = dEdx * z * rho;
   double E = sqrt(pow(pz,2)+pow(105.65,2));
   E -= Eloss;
   double ploss = sqrt(pow(E,2) - pow(105.65,2));
   E += Eloss;
   double delta = pz - ploss;
   /*
   std::cout << "delta " << delta << std::endl;
   std::cout << "E " << E << std::endl;
   */

   // time0
   double t1 = s1*sqrt(pow(pz+delta,2)+pow(105.65,2))/(pz+delta);
   double t2 = s2*sqrt(pow(pz-delta,2)+pow(105.65,2))/(pz-delta);
   double t0 = t1 + t2;
   /*
   std::cout << "t2 " << t2 << std::endl;
   std::cout << "t1 " << t1 << std::endl;
   std::cout << "t0 " << t0 << std::endl;
   */
   
   // Root finding
   /*
   MyFunction1D myf1;
   myf1.E = E;
   myf1.delta = delta;
   myf1.s1 = s1;
   myf1.s2 = s2;
   std::cout << "myf1.operator()(pz) " << myf1.operator()(pz) << std::endl;
   std::cout << "myf1.Derivative(pz) " << myf1.Derivative(pz) << std::endl;
   ROOT::Math::GradFunctor1D  f1(myf1); 
   ROOT::Math::RootFinder rfn(ROOT::Math::RootFinder::kGSL_STEFFENSON);
   rfn.SetFunction(f1, pz);
   rfn.Solve();
   cout << rfn.Root() << endl;
   */

   // delta_t
   double term3i = 0.5* (s1 + s2)*pow(105.65,2)*pow(delta,2)/(pow(pz,2)*pow(E,2));
   double term3ii = pz/E + 2*E/pz;
   double term2 = -(s1 - s2) * delta * pow(105.65,2)/ (pow(pz,2)*E);
   double dtime = term2 + term3i*term3ii;
   /*
   std::cout << "term3i " << term3i << std::endl;
   std::cout << "term3ii " << term3ii << std::endl;
   std::cout << "term2 " << term2 << std::endl;
   std::cout << "dtime " << dtime << std::endl;
   */

   // Corrected P
   double time = t0 + dtime; 
   p_corrected = pz - ((1/(s1-s2))*E*pow(pz,2)*dtime)/pow(105.65,2); 
   //std::cout << "pz " << pz << std::endl;
   pz = p_corrected;
   t_path = t0+dtime;
   //std::cout << "p_corrected " << p_corrected << std::endl;
   /*
   std::cout << "p_corrected " << p_corrected << std::endl;
   std::cout << "counter " << counter << std::endl;
   std::cout << "iteration " << i << std::endl;
   std::cout << "t_path " << t_path << std::endl;
   */
   
   }

   // Iteration
   /*
   counter += 1;
   int i = 0;
   for (i;i<1000;i++) 
   while (t0_first*1.1-t_path>0.1)
   t_cor->Fill(i,t_path);
   if (counter==3) { 
	   std::cout << "doing graph" << std::endl;
	   TCanvas *c1 = new TCanvas();
   t_cor->Draw();
   c1->SaveAs("t_cor.pdf");
   c1->Clear();
   }
   */
   
   cor_mom->Fill(p_corrected);
   double MCTruth_pz_up = 0;
   double MCTruth_pz_down = 0;
   for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-16803.7<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-16803.7>-50) {
		   MCTruth_pz_up = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
	   }
	   if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-17101.3<50 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-17101.3>-50) {
		   MCTruth_pz_down = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
	   }

	   }
   std::cout << "MCTruth_pz_up " << MCTruth_pz_up << std::endl;
   std::cout << "MCTruth_pz_down " << MCTruth_pz_down << std::endl;
   if (MCTruth_pz_up != 0 && MCTruth_pz_down != 0){
   double true_delta = MCTruth_pz_up - MCTruth_pz_down;	   
   double MCTruth_pz = MCTruth_pz_up+MCTruth_pz_down/2;
   double residual = p_corrected - MCTruth_pz;
   std::cout << "true_delta " << true_delta << std::endl;
   std::cout << "MCTruth_pz " << MCTruth_pz_up+MCTruth_pz_down/2 << std::endl;
   std::cout << "residual " << residual << std::endl;
   residual_plot->Fill(MCTruth_pz,p_corrected);
   }
   }

}

double MCSAnalysis::TimeofFlight(){
  double rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
    rawTOF1HitTime  = 
      tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
    rawTOF2HitTime = 
      tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
  double dt  = 100.0; // something rediculously large as a default number
  if ( rawTOF1HitTime != -1 && rawTOF2HitTime != -1 ){
    dt  = rawTOF2HitTime - rawTOF1HitTime; 
  }
  return dt;
}

double MCSAnalysis::BetheBloch(double pz){

   double beta = pow(pz,2)/(pow(105.65,2)+pow(pz,2));
   double gamma = 1/sqrt(1-pow(beta,2));
   double W = 2*0.511*pow(beta,2)*pow(gamma,2)/(1+2*gamma*0.511/105.65+pow(0.511/105.65,2));
   double I = 36.5e-6;
   double density = log(18.51e-6/I)+log(beta*gamma)-1/2;
   double dEdxpre = 0.307075*2/(pow(beta,2)*7.94894);   
   double dEdxterm1 = 0.5*log(2*0.511*pow(beta,2)*pow(gamma,2)*W/pow(I,2)) - pow(beta,2) - density;
   double dEdx = dEdxpre * dEdxterm1;
   /*
   std::cout << "beta " << beta << std::endl;
   std::cout << "gamma " << gamma << std::endl;
   std::cout << "W " << W << std::endl;
   std::cout << "density " << density << std::endl;
   std::cout << "dEdxpre " << dEdxpre << std::endl;
   std::cout << "dEdxterm1 " << dEdxterm1 << std::endl;
   std::cout << "dEdx " << dEdx << std::endl;
*/

   return dEdx;
}

double MCSAnalysis::MomentumFromTOF(bool isdata=true){
    // Cuts remove events where the following statements do not make sense so we proceed without cuts.
  double rawTOF0HitTime = -1., rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0)
    rawTOF0HitTime  = 
      tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
    rawTOF1HitTime  = 
      tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  if( int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
    rawTOF2HitTime = 
      tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
  // tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArray()[1].GetRawTime() ) / 2.;
  // need to hard code a few things here (unfortunately) pertaining to the geometry and the time of flight.
  double dt0 = (_sys["TOF1_z"] - _sys["TOF0_z"]) / 0.299792458 / 1000.; // ns. 
  double dt  = 10.0 * dt0; // something rediculously large as a default number
  double pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
  if ( rawTOF1HitTime != -1 && rawTOF0HitTime != -1 ){
    dt  = rawTOF1HitTime - rawTOF0HitTime; 
    pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0) - 36.1;
  }
  /*
  if ( rawTOF1HitTime != -1 && rawTOF2HitTime != -1 ){
    // Better estimate of the longitudinal momentum
    dt0 = (_sys["TOF2_z"] - _sys["TOF1_z"]) / 0.299792458 / 1000.; // ns. 
    dt  = rawTOF2HitTime - rawTOF1HitTime; 
    double pz1 = pz;
    pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
    TOFcom->Fill(pz1,pz);
  }
  */
  if(isdata)
    calc_mom->Fill(pz);
  else
    mccalc_mom->Fill(pz);
  return pz;
}

void MCSAnalysis::ConvolveWithInputDistribution(std::string distname){
  int isfirst = 0;
  bool isGEANT;
  bool isCobb;
  bool isMoliere;
  if (distname.find(modelname1.c_str()) != std::string::npos)
    isGEANT = true;
  if (distname.find(modelname2.c_str()) != std::string::npos)
    isCobb = true;
  if (distname.find(modelname3.c_str()) != std::string::npos)
    isMoliere = true;

  TFile* infile = new TFile(modelfile.c_str());

  std::string tmpname = "thetaX_refconv_";
  tmpname += distname;
  TH1D* thetaX_refconv = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmpname = "thetaY_refconv_";
  tmpname += distname;
  TH1D* thetaY_refconv = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmpname = "thetaScatt_refconv_";
  tmpname += distname;
  TH1D* thetaScatt_refconv = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",
	     _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);


  tmpname = "theta2Scatt_refconv_";
  tmpname += distname;
  TH1D* theta2Scatt_refconv = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",
	     _histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  tmpname = "thetaScatt_refconv_vp_";
  tmpname += distname;
  TH2D* thetaScatt_refconv_vp = 
    new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #theta_{Scatt}", 
	     200, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
    
  TH2D* thetaXUS_thetaXDS = 
    new TH2D("thetaXUS_thetaXDS","Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], 
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH2D* thetaYUS_thetaYDS = 
    new TH2D("thetaYUS_thetaYDS","Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmpname = "thetaX_";
  tmpname += distname;

  std::cout<<"Convolution with "<<tmpname<<" from "<<modelfile<<std::endl;
  TH1D* hx = (TH1D*)infile->Get(tmpname.c_str());
  TH1D* hy = (TH1D*)infile->Get(tmpname.c_str());

  // Collection DSConvSet;
  // for (int l=0; l<10; l++){
  for (int i=0; i<_USMCset.N(); i++){
    for (int j=0; j<20; j++){
      double dthetaX = hx->GetRandom() * _sys["resX"];
      double dthetaY = hy->GetRandom() * _sys["resY"];
      // First project the upstream track to the absorber 
      double zabspos = _sys["abspos"] + 0.0;
      Vars projvarAbs = PropagateVarsMu(_USMCset.E(i), zabspos);
      double xabs = projvarAbs.X;  /// _USMCset.E(i).X + _USMCset.E(i).dXdz * dzabsUS;
      double yabs = projvarAbs.Y;  /// _USMCset.E(i).Y + _USMCset.E(i).dYdz * dzabsUS;
      // Now add the angle from the model to the downstream measurement.
      double dXdz_abs = _DSMCset.E(i).dXdz + tan(dthetaY);
      double dYdz_abs = _DSMCset.E(i).dYdz + tan(dthetaX);
      // double d_thetaY    = atan(dXdz_abs) - atan(_DSMCset.E(i).dXdz);
      // double d_thetaX    = atan(dYdz_abs) - atan(_DSMCset.E(i).dYdz);
      
      // Project the track into the downstream reference plane
      // double xref = xabs + dXdz_abs * dzabsDS;
      // double yref = yabs + dYdz_abs * dzabsDS;
      Vars tmpvar = _DSMCset.E(i);
      tmpvar.X = xabs;
      tmpvar.Y = yabs;
      tmpvar.Z = _sys["abspos"];
      tmpvar.dXdz = dXdz_abs;
      tmpvar.dYdz = dYdz_abs;
      tmpvar.px   = dXdz_abs * _USMCset.E(i).pz;
      tmpvar.py   = dYdz_abs * _USMCset.E(i).pz;
      tmpvar.pz   = _USMCset.E(i).pz;
      tmpvar.TOF12= _USMCset.E(i).TOF12;
      tmpvar.TOF01= _USMCset.E(i).TOF01;
      std::vector<double> projDTheta = DefineProjectionAngles(tmpvar, _DSMCset.E(i));
      double d_thetaX = projDTheta[0];
      double d_thetaY = projDTheta[1];

      //double cosdthetaScatt = ( (1 + dXdz_abs * _DSMCset.E(i).dXdz + dYdz_abs*_DSMCset.E(i).dYdz) / 
      //			sqrt( 1 + dXdz_abs*dXdz_abs + dYdz_abs*dYdz_abs) / 
      //			sqrt(1 + _DSMCset.E(i).dXdz*_DSMCset.E(i).dXdz + _DSMCset.E(i).dYdz * _DSMCset.E(i).dYdz) );
      double dthetaScatt = projDTheta[2];  /// acos(cosdthetaScatt);
      Vars projvar = PropagateVarsMu(tmpvar, _sys["abspos"] + 2993.05);
      // Remove events that do not pass through the tracker at its centre
      // double xtracker = xabs + dXdz_abs * (_sys["abspos"] + 549.95);
      // double ytracker = yabs + dYdz_abs * (dzabsDS + 549.95);
      
      if( sqrt(projvar.X*projvar.X + projvar.Y*projvar.Y) > meanp ){
	tmpvar.X = 0.0;
	tmpvar.Y = 0.0;
	tmpvar.Z = 0.0;
	tmpvar.dXdz = 1./2.;
	tmpvar.dYdz = 1./2.;
      }
      std::vector<double> projTheta = DefineProjectionAngles(_USMCset.E(i), tmpvar);
      double thetaY = projTheta[1];   /// atan(tmpvar.dXdz) - atan(_USMCset.E(i).dXdz);
      double thetaX = projTheta[0];   /// atan(tmpvar.dYdz) - atan(_USMCset.E(i).dYdz);
      // double cosScatt = ( (1 + _USMCset.E(i).dXdz * tmpvar.dXdz +
      //		   _USMCset.E(i).dYdz * tmpvar.dYdz )/
      //		  sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
      //		       _USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
      //		  sqrt(1 + tmpvar.dXdz*tmpvar.dXdz +
      //		       tmpvar.dYdz*tmpvar.dYdz));
      double thetaScatt = projTheta[2];  /// acos(cosScatt);
      thetaXUS_thetaXDS->Fill(atan(_USMCset.E(i).dXdz), atan(tmpvar.dXdz));
      thetaYUS_thetaYDS->Fill(atan(_USMCset.E(i).dYdz), atan(tmpvar.dYdz));
      thetaX_refconv->Fill(thetaX);
      thetaY_refconv->Fill(thetaY);
      thetaScatt_refconv->Fill(thetaScatt);
      theta2Scatt_refconv->Fill(thetaScatt*thetaScatt);
      thetaScatt_refconv_vp->Fill(_USMCset.E(i).pz, thetaScatt);
     
      if (isGEANT) {
    	  resp_thetaX.Fill(thetaX, d_thetaX);
    	  resp_thetaY.Fill(thetaY, d_thetaY);
    	  resp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      if (isCobb) {
    	  tresp_thetaX.Fill(thetaX, d_thetaX);
    	  tresp_thetaY.Fill(thetaY, d_thetaY);
    	  tresp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      if (isMoliere) {
    	  mresp_thetaX.Fill(thetaX, d_thetaX);
    	  mresp_thetaY.Fill(thetaY, d_thetaY);
    	  mresp_thetaScatt.Fill(thetaScatt, dthetaScatt);
    	  mresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
      }
      /*
      isfirst == 1 ? tresp_thetaX.Fill(thetaX, d_thetaX) : resp_thetaX.Fill(thetaX, d_thetaX);
      std::cout << "thetaX " << thetaX << std::endl;
      isfirst == 1 ? tresp_thetaY.Fill(thetaY, d_thetaY) : resp_thetaY.Fill(thetaY, d_thetaY);
      std::cout << "thetaY " << thetaY << std::endl;
      isfirst == 1 ? tresp_thetaScatt.Fill(thetaScatt, dthetaScatt) : resp_thetaScatt.Fill(thetaScatt, dthetaScatt); 
      std::cout << "thetaScatt " << thetaScatt << std::endl;
      isfirst == 1 ? tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt) : resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt); 
      */
    }
  }

  outfile->cd();
  thetaXUS_thetaXDS->Write();
  thetaYUS_thetaYDS->Write();
  thetaX_refconv->Write();
  thetaY_refconv->Write();
  thetaScatt_refconv->Write();
  theta2Scatt_refconv->Write();
  thetaScatt_refconv_vp->Write();

}

void MCSAnalysis::DoUnfolding(){

  std::string  tmpname = "thetaX_recoGold";
  TH1D* thetaX_recoGold = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaY_recoGold";
  TH1D* thetaY_recoGold = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaScatt_recoGold";
  TH1D* thetaScatt_recoGold = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "theta2Scatt_recoGold";
  TH1D* theta2Scatt_recoGold = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  TH1D* thetaX_response = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaY_response";
  TH1D* thetaY_response = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaScatt_response";
  TH1D* thetaScatt_response = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "theta2Scatt_response";
  TH1D* theta2Scatt_response = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);


  for(int i=0; i<_DSset.N(); i++){
    
    double thetaY = atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
    double thetaX = atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
    double cosScatt = ( (1 + _USset.E(i).dXdz * _DSset.E(i).dXdz +
			 _USset.E(i).dYdz * _DSset.E(i).dYdz )/
			sqrt(1 + _USset.E(i).dXdz*_USset.E(i).dXdz +
			     _USset.E(i).dYdz*_USset.E(i).dYdz)/
			sqrt(1 + _DSset.E(i).dXdz*_DSset.E(i).dXdz +
			     _DSset.E(i).dYdz*_DSset.E(i).dYdz));
    double thetaScatt = acos(cosScatt);

    thetaX_recoGold->Fill(thetaX);
    thetaY_recoGold->Fill(thetaY);
    thetaScatt_recoGold->Fill(thetaScatt);
    theta2Scatt_recoGold->Fill(thetaScatt*thetaScatt);

  }
  float* source_thetaX = new float[int(_histlimits["NbinsXY"])];
  float* source_thetaY = new float[int(_histlimits["NbinsXY"])];
  float* source_thetaScatt = new float[int(_histlimits["NbinsTh"])];
  float* source_theta2Scatt = new float[int(_histlimits["NbinsTh2"])];

  for(int i=0;i<_DSMCset.N(); i++){
    std::vector<double> projTheta = DefineProjectionAngles(_USMCset.E(i), _DSMCset.E(i));
    double thetaMCY = projTheta[1]; /// atan(_DSMCset.E(i).dXdz) - atan(_USMCset.E(i).dXdz);
    double thetaMCX = projTheta[0]; /// atan(_DSMCset.E(i).dYdz) - atan(_USMCset.E(i).dYdz);
    // double cosScattMC = ( (1 + _USMCset.E(i).dXdz * _DSMCset.E(i).dXdz +
    //_USMCset.E(i).dYdz * _DSMCset.E(i).dYdz )/
    //sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
    //_USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
    // sqrt(1 + _DSMCset.E(i).dXdz*_DSMCset.E(i).dXdz +
    // _DSMCset.E(i).dYdz*_DSMCset.E(i).dYdz));
    double thetaScattMC = projTheta[2];  /// acos(cosScattMC);
    thetaX_response->Fill(thetaMCX);
    thetaY_response->Fill(thetaMCY);
    thetaScatt_response->Fill(thetaScattMC);
    theta2Scatt_response->Fill(thetaScattMC*thetaScattMC);

  }
  
  float* response_thetaX = new float[int(_histlimits["NbinsXY"])];
  float* response_thetaY = new float[int(_histlimits["NbinsXY"])];
  float* response_thetaScatt = new float[int(_histlimits["NbinsTh"])];
  float* response_theta2Scatt = new float[int(_histlimits["NbinsTh2"])];

  for(int i=0; i<_histlimits["NbinsXY"]; i++){
    source_thetaX[i] = thetaX_recoGold->GetBinContent(i+1);
    source_thetaY[i] = thetaY_recoGold->GetBinContent(i+1);
    response_thetaX[i] = thetaX_response->GetBinContent(i+1);
    response_thetaY[i] = thetaY_response->GetBinContent(i+1);
  }
  for(int i=0; i<_histlimits["NbinsTh"]; i++){
    source_thetaScatt[i] = thetaScatt_recoGold->GetBinContent(i+1);
    response_thetaScatt[i] = thetaScatt_response->GetBinContent(i+1);
  }
  for(int i=0; i<_histlimits["NbinsTh2"]; i++){
    source_theta2Scatt[i] = theta2Scatt_recoGold->GetBinContent(i+1);
    response_theta2Scatt[i] = theta2Scatt_response->GetBinContent(i+1);
  }

  /*
  TSpectrum *sX = new TSpectrum();
  sX->Deconvolution(source_thetaX, response_thetaX, _histlimits["NbinsXY"], 20, 10, 0.0);

  TSpectrum *sY = new TSpectrum();
  sY->Deconvolution(source_thetaY, response_thetaY, _histlimits["NbinsXY"], 20, 10, 0.0);
  for(int i=int(_histlimits["NbinsXY"]*3./4.); 
      i<int(_histlimits["NbinsXY"]*(3./4. + 1./5.)); i++){
    thetaX_recoGold->SetBinContent(i-148, source_thetaX[i]);
    thetaY_recoGold->SetBinContent(i-148, source_thetaY[i]);
  }
  double maxX = thetaX_recoGold->GetMaximum();
  double maxY = thetaY_recoGold->GetMaximum();
  for(int i = 1; i <= _histlimits["NbinsXY"]; i++){
    thetaX_recoGold->SetBinContent(i, thetaX_recoGold->GetBinContent(i)/maxX);
    thetaY_recoGold->SetBinContent(i, thetaY_recoGold->GetBinContent(i)/maxY);
  }
  TSpectrum *sScatt = new TSpectrum();
  sScatt->Deconvolution(source_thetaScatt, response_thetaScatt, _histlimits["NbinsTh"], 20, 10, 0.0);
  for(int i=0; i<int(_histlimits["NbinsTh"]); i++){
    thetaScatt_recoGold->SetBinContent(i, source_thetaScatt[i]);
  }
  TSpectrum *sScatt2 = new TSpectrum();
  sScatt2->Deconvolution(source_theta2Scatt, response_theta2Scatt, _histlimits["NbinsTh2"], 20, 10, 0.0);
  
  for(int i=0; i<int(_histlimits["NbinsTh2"]*(3./4.  + 1./5.)); i++){
    theta2Scatt_recoGold->SetBinContent(i, source_theta2Scatt[i]);
  }
  // double maxX = thetaScatt_recoGold->GetMaximum();
  //for(int i = 1; i <= _histlimits["NbinsTh"]; i++){
  //  thetaScatt_recoGold->SetBinContent(i, thetaScatt_recoGold->GetBinContent(i)/maxX);
  // }

  outfile->cd();
  thetaX_recoGold->Write();
  thetaY_recoGold->Write();
  thetaScatt_recoGold->Write();
  theta2Scatt_recoGold->Write();

  delete thetaX_recoGold;
  delete thetaY_recoGold;
  delete thetaScatt_recoGold;
  delete theta2Scatt_recoGold;

  delete thetaX_response;
  delete thetaY_response;
  delete thetaScatt_response;
  delete theta2Scatt_response;
  */
}


void MCSAnalysis::DoDeconvolution(std::string model, int n_sel=1){
  // The basic methods associated with the RooUnfolding package
  // First generate a histogram of the measured data.
  
  int isfirst = 0;
  if (model.find(modelname1.c_str()) != std::string::npos)
    isfirst = 1;

  // int n_base = int(2 * _DSset.N() / (n_sel * (n_sel + 1)));
  // int k=0;
  // for (int j=0; j<n_sel; j++){
  // define the number of events under the current selection
  int curr_sel = int(_DSset.N());  // n_base * j;
  int curr_k = 0;

  
    
  std::string  tmpname = "thetaX_measdata";
  tmpname += model;
  TH1D* thetaX_measdata = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaY_measdata";
  tmpname += model;
  TH1D* thetaY_measdata = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaScatt_measdata";
  tmpname += model;
  TH1D* thetaScatt_measdata = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "theta2Scatt_measdata";
  tmpname += model;
  TH1D* theta2Scatt_measdata = 
    new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  tmpname = "thetaScatt_measdata_vp";
  tmpname += model;
  TH2D* thetaScatt_measdata_vp = 
    new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #theta_{Scatt}", 
	     400, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  tmpname = "thetaXUS_thetaXDS";
  tmpname += model;
  
  TH2D* thetaXUS_thetaXDS = 
    new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmpname = "thetaYUS_thetaYDS";
  tmpname += model;
  TH2D* thetaYUS_thetaYDS = 
    new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#theta_{X}^{US}; #theta_{X}^{DS}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  
  tmpname = "projposUSDSdiff_";
  tmpname += model;
  TH2D* projposUSDSdiff = 
    new TH2D(tmpname.c_str(),
	     "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
	     180, -400, 400, 180, -400, 400);
  
  
  const Int_t NBINS = 19;
  Double_t scat_bin_array[NBINS + 1] = {-0.1151,-0.0938,-0.0754,-0.0597,-0.0463,-0.0347,-0.0248,-0.0162,-0.00895,-0.00269,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,0.0754,0.0938,0.1151};
  TH1D* scattering_proj_x = new TH1D("scattering_proj_x_DC","Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
				     NBINS, scat_bin_array);
  TH1D* scattering_proj_y = new TH1D("scattering_proj_y_DC","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
				     NBINS, scat_bin_array);
  
  
  for(int i=curr_k; i<curr_sel; i++){
    std::vector<double> projTheta = DefineProjectionAngles(_USset.E(i), _DSset.E(i));
    double thetaY = projTheta[1];  /// atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
    double thetaX = projTheta[0];  /// atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
    // double cosScatt 
    double thetaScatt = projTheta[2];  /// acos(cosScatt);
    double xUSabs = _USset.E(i).X - _USset.E(i).dXdz * (_USset.E(i).Z - _sys["abspos"]);
    double yUSabs = _USset.E(i).Y - _USset.E(i).dYdz * (_USset.E(i).Z - _sys["abspos"]);
    double xDSabs = _DSset.E(i).X - _DSset.E(i).dXdz * (_DSset.E(i).Z - _sys["abspos"]);
    double yDSabs = _DSset.E(i).Y - _DSset.E(i).dYdz * (_DSset.E(i).Z - _sys["abspos"]);
    projposUSDSdiff->Fill(xDSabs - xUSabs, yDSabs - yUSabs);
    thetaXUS_thetaXDS->Fill(atan(_USset.E(i).dXdz), atan(_DSset.E(i).dXdz));
    thetaYUS_thetaYDS->Fill(atan(_USset.E(i).dYdz), atan(_DSset.E(i).dYdz));
    thetaX_measdata->Fill(thetaX);
    thetaY_measdata->Fill(thetaY);
    thetaScatt_measdata->Fill(thetaScatt);
    theta2Scatt_measdata->Fill(thetaScatt*thetaScatt);
    // double dt0 = 7.64186 / 0.299792458; // ns. 
    // double dt = _USset.E(i).TOF01;
    // if (_USset.E(i).TOF12 < 99.0 * 8.22475 / 0.299792458){
    //  dt0 = 8.22475 / 0.299792458; // ns.
    //  dt = _USset.E(i).TOF12;
    //}
    double pz  = _USset.E(i).pz;
    // 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
    thetaScatt_measdata_vp->Fill(pz, thetaScatt);
    scattering_proj_x->Fill(thetaX);
    scattering_proj_y->Fill(thetaY);
    // k++;
  }
  /*
    for ( int itr=1; itr<25; itr++){    
    RooUnfoldBayes unfold_thetaX(isfirst==1 ? &tresp_thetaX : &resp_thetaX, thetaX_measdata, itr);
    unfold_thetaX.Print();
    char res[16];
    sprintf( res, "%d", itr);
    TH1D* thetaX_reco = (TH1D*)unfold_thetaX.Hreco();
    tmpname = "thetaX_reco";
    tmpname += model;
    tmpname += "_";
    tmpname += res;
    thetaX_reco->SetName(tmpname.c_str());
    thetaX_reco->SetTitle(";#Delta #theta_{X}");
    outfile->cd();
    thetaX_reco->Write();
    
    delete thetaX_reco;
    }
    */
  if (model.find(modelname1.c_str()) != std::string::npos) {
	  fresp_thetaX = resp_thetaX;
	  fresp_thetaY = resp_thetaY;
	  fresp_thetaScatt = resp_thetaScatt;
	  fresp_theta2Scatt = resp_theta2Scatt;
  }
  if (model.find(modelname2.c_str()) != std::string::npos) {
	  fresp_thetaX = tresp_thetaX;
	  fresp_thetaY = tresp_thetaY;
	  fresp_thetaScatt = tresp_thetaScatt;
	  fresp_theta2Scatt = tresp_theta2Scatt;
  }
  if (model.find(modelname3.c_str()) != std::string::npos) {
	  fresp_thetaX = mresp_thetaX;
	  fresp_thetaY = mresp_thetaY;
	  fresp_thetaScatt = mresp_thetaScatt;
	  fresp_theta2Scatt = mresp_theta2Scatt;
  }

  RooUnfoldBayes unfold_thetaX(&fresp_thetaX, thetaX_measdata, int(_sys["niter"]));
  unfold_thetaX.Print();
  RooUnfoldBayes unfold_thetaY(&fresp_thetaY, thetaY_measdata, int(_sys["niter"]));
  unfold_thetaY.Print();
  RooUnfoldBayes 
    unfold_thetaScatt(&fresp_thetaScatt, thetaScatt_measdata, int(_sys["niter"]));
  unfold_thetaScatt.Print();
  RooUnfoldBayes 
    unfold_theta2Scatt(&fresp_theta2Scatt, theta2Scatt_measdata, int(_sys["niter"]));
  unfold_theta2Scatt.Print();
  
  for (int i=1; i<=19; i++){
    double x_content = scattering_proj_x->GetBinContent(i);
    scattering_proj_x->Fill(i, x_content / scat_bin_array[i-1]);
    double y_content = scattering_proj_y->GetBinContent(i);
    scattering_proj_y->Fill(i, y_content / scat_bin_array[i-1]);
  }
  
  TH1D* thetaX_reco = (TH1D*)unfold_thetaX.Hreco();
  tmpname = "thetaX_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_reco->SetName(tmpname.c_str());
  thetaX_reco->SetTitle(";#Delta #theta_{X}");
  TH1D* thetaY_reco = (TH1D*)unfold_thetaY.Hreco();
  tmpname = "thetaY_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_reco->SetName(tmpname.c_str());
  thetaY_reco->SetTitle(";#Delta #theta_{Y}");
  TH1D* thetaScatt_reco = (TH1D*)unfold_thetaScatt.Hreco();
  tmpname = "thetaScatt_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_reco->SetName(tmpname.c_str());
  thetaScatt_reco->SetTitle(";#theta_{Scatt}");
  TH1D* theta2Scatt_reco = (TH1D*)unfold_theta2Scatt.Hreco();
  tmpname = "theta2Scatt_reco";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_reco->SetName(tmpname.c_str());
  theta2Scatt_reco->SetTitle(";#theta^{2}_{Scatt}");
  
  TH2D* thetaX_response = (TH2D*)unfold_thetaX.response()->Hresponse();
  tmpname = "thetaX_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_response->SetName(tmpname.c_str());
  thetaX_response->SetTitle(";#Delta #theta^{rec}_{X}; #Delta #theta^{true}_{X}");
  TH2D* thetaY_response = (TH2D*)unfold_thetaY.response()->Hresponse();
  tmpname = "thetaY_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_response->SetName(tmpname.c_str());
  thetaY_response->SetTitle(";#Delta #theta^{rec}_{Y}; #Delta #theta^{true}_{Y}");
  TH2D* thetaScatt_response = (TH2D*)unfold_thetaScatt.response()->Hresponse();
  tmpname = "thetaScatt_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_response->SetName(tmpname.c_str());
  thetaScatt_response->SetTitle(";#theta^{tracker}_{Scatt};#theta^{abs}_{Scatt}");
  TH2D* theta2Scatt_response = (TH2D*)unfold_theta2Scatt.response()->Hresponse();
  tmpname = "theta2Scatt_response";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_response->SetName(tmpname.c_str());
  theta2Scatt_response->SetTitle(";(#theta^{tracker}_{Scatt})^{2}; (#theta^{abs}_{Scatt})^{2}");
  
  TH1D* thetaX_truth = (TH1D*)unfold_thetaX.response()->Htruth();
  tmpname = "thetaX_";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_truth->SetName(tmpname.c_str());
  thetaX_truth->SetTitle(";#Delta #theta_{X}");
  TH1D* thetaY_truth = (TH1D*)unfold_thetaY.response()->Htruth();
  tmpname = "thetaY_";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_truth->SetName(tmpname.c_str());
  thetaX_truth->SetTitle(";#Delta #theta_{Y}");
  TH1D* thetaScatt_truth = (TH1D*)unfold_thetaScatt.response()->Htruth();
  tmpname = "thetaScatt_";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_truth->SetName(tmpname.c_str());
  thetaScatt_truth->SetTitle(";#theta_{Scatt}");
  TH1D* theta2Scatt_truth = (TH1D*)unfold_theta2Scatt.response()->Htruth();
  tmpname = "theta2Scatt_";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_truth->SetName(tmpname.c_str());
  theta2Scatt_truth->SetTitle(";#theta^{2}_{Scatt}");
  
  TH1D* thetaX_measured = (TH1D*)unfold_thetaX.response()->Hmeasured();
  tmpname = "thetaX_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaX_measured->SetName(tmpname.c_str());
  thetaX_measured->SetTitle(";#Delta #theta_{X}");
  TH1D* thetaY_measured = (TH1D*)unfold_thetaY.response()->Hmeasured();
  tmpname = "thetaY_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaY_measured->SetName(tmpname.c_str());
  thetaY_measured->SetTitle(";#Delta #theta_{Y}");
  TH1D* thetaScatt_measured = (TH1D*)unfold_thetaScatt.response()->Hmeasured();
  tmpname = "thetaScatt_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  thetaScatt_measured->SetName(tmpname.c_str());
  thetaScatt_measured->SetTitle(";#theta_{Scatt}");
  TH1D* theta2Scatt_measured = (TH1D*)unfold_theta2Scatt.response()->Hmeasured();
  tmpname = "theta2Scatt_measured";
  tmpname += model;
  // if(j>0) tmpname += j;
  theta2Scatt_measured->SetName(tmpname.c_str());
  theta2Scatt_measured->SetTitle(";#theta^{2}_{Scatt}");
  
  outfile->cd();
  thetaXUS_thetaXDS->Write();
  thetaYUS_thetaYDS->Write();
  thetaX_measdata->Write();
  thetaY_measdata->Write();
  thetaScatt_measdata->Write();
  theta2Scatt_measdata->Write();
  thetaScatt_measdata_vp->Write();
  thetaX_reco->Write();
  thetaY_reco->Write();
  thetaScatt_reco->Write();
  theta2Scatt_reco->Write();
  thetaX_response->Write();
  thetaY_response->Write();
  thetaScatt_response->Write();
  theta2Scatt_response->Write();
  thetaX_truth->Write();
  thetaY_truth->Write();
  thetaScatt_truth->Write();
  theta2Scatt_truth->Write();
  thetaX_measured->Write();
  thetaY_measured->Write();
  thetaScatt_measured->Write();
  theta2Scatt_measured->Write();
  scattering_proj_x->Write();
  scattering_proj_y->Write();
  projposUSDSdiff->Write();
  
  TCanvas* c1 = new TCanvas();
  thetaX_measdata->Draw();
  c1->Print("thetaX_measdata.pdf");
  c1->Clear();
  thetaY_measdata->Draw();
  c1->Print("thetaY_measdata.pdf");
  c1->Clear();
  thetaScatt_measdata->Draw();
  c1->Print("thetaScatt_measdata.pdf");
  c1->Clear();
  thetaScatt_measdata_vp->Draw("colz");
  c1->Print("thetaScatt_measdata_vp.pdf");
  c1->Clear();
  thetaX_reco->Draw();
  c1->Print("thetaX_reco.pdf");
  c1->Clear();
  thetaY_reco->Draw();
  c1->Print("thetaY_reco.pdf");
  c1->Clear();
  thetaScatt_reco->Draw();
  c1->Print("thetaScatt_reco.pdf");
  c1->Clear();
  thetaX_response->Draw();
  c1->Print("thetaX_response.pdf");
  c1->Clear();
  thetaY_response->Draw();
  c1->Print("thetaY_response.pdf");
  c1->Clear();
  thetaScatt_response->Draw();
  c1->Print("thetaScatt_response.pdf");
  c1->Clear();
  thetaX_measured->Draw();
  c1->Print("thetaX_measured.pdf");
  c1->Clear();
  thetaY_measured->Draw();
  c1->Print("thetaY_measured.pdf");
  c1->Clear();
  thetaScatt_measured->Draw();
  c1->Print("thetaScatt_measured.pdf");
  c1->Clear();
  thetaX_truth->Draw();
  c1->Print("thetaX_truth.pdf");
  c1->Clear();
  thetaY_truth->Draw();
  c1->Print("thetaY_truth.pdf");
  c1->Clear();
  thetaScatt_truth->Draw();
  c1->Print("thetaScatt_truth.pdf");
  c1->Clear();
  c1->SetLogy();
  scattering_proj_x->Draw();
  c1->Print("scattering_proj_x.pdf");
  c1->Clear();
  scattering_proj_y->Draw();
  c1->Print("scattering_proj_y.pdf");
  c1->Clear();
  
  delete c1;
  delete thetaX_measdata;
  delete thetaY_measdata;
  delete thetaScatt_measdata;
  delete thetaScatt_measdata_vp;
  delete thetaX_reco;
  delete thetaY_reco;
  delete thetaScatt_reco;
  delete thetaX_response;
  delete thetaY_response;
  delete thetaScatt_response;
  delete thetaX_truth;
  delete thetaY_truth;
  delete thetaScatt_truth;
  delete thetaX_measured;
  delete thetaY_measured;
  delete thetaScatt_measured;
  delete scattering_proj_x;
  delete scattering_proj_y;
  delete projposUSDSdiff;

}

void MCSAnalysis::FitGaussian(std::string outfilename){

	TFile *MyFile = new TFile(outfilename.c_str());
	TH1F * h1 = (TH1F*)MyFile->Get("thetaX_measdataCobb");
	h1->SetDirectory(0);
	h1->Draw();
	std::cout << "Fit result for thetaX_measdata " << std::endl;
	h1->Fit("gaus","","",-0.040,0.040);

	TH1F * h2 = (TH1F*)MyFile->Get("thetaY_measdataCobb");
	h2->SetDirectory(0);
	h2->Draw();
	std::cout << "Fit result for thetaY_measdata " << std::endl;
	h2->Fit("gaus","","",-0.040,0.040);

	TH1F * h3 = (TH1F*)MyFile->Get("thetaX_recoCobb");
	h3->SetDirectory(0);
	h3->Draw();
	std::cout << "Fit result for thetaX_reco " << std::endl;
	h3->Fit("gaus","","",-0.040,0.040);
	
	TH1F * h4 = (TH1F*)MyFile->Get("thetaY_recoCobb");
	h4->SetDirectory(0);
	h4->Draw();
	std::cout << "Fit result for thetaY_reco " << std::endl;
	h4->Fit("gaus","","",-0.040,0.040);
	
	TH1F * h5 = (TH1F*)MyFile->Get("thetaX_refconv_Cobb");
	h5->SetDirectory(0);
	h5->Draw();
	std::cout << "Fit result for thetaX_refconv_Cobb " << std::endl;
	h5->Fit("gaus","","",-0.040,0.040);
	
	TH1F * h6 = (TH1F*)MyFile->Get("thetaY_refconv_Cobb");
	h6->SetDirectory(0);
	h6->Draw();
	std::cout << "Fit result for thetaY_refconv_Cobb " << std::endl;
	h6->Fit("gaus","","",-0.040,0.040);
	
	MyFile->Close();
	delete h1;
	delete h2;
	delete h3;
	delete h4;
	delete h5;
	delete h6;
}

void MCSAnalysis::CalculateChi2(std::string outfilename, std::string distname){

	std::string tmpname = "thetaX_refconv_";
        tmpname += distname;
	TFile *MyFile = new TFile(outfilename.c_str());
        TH1F * hx = (TH1F*)MyFile->Get(tmpname.c_str());
	tmpname = "thetaX_measdata";
        tmpname += distname;
	TH1F * h1 = (TH1F*)MyFile->Get(tmpname.c_str());
	h1->SetDirectory(0);
	hx->SetDirectory(0);
	//hx->Rebin(GetNbinsX/40);
	h1->Draw();
	hx->Draw();
	std::cout << "Chi^2 for thetaX_reco and " << distname << std::endl;
	h1->Chi2Test(hx,"P");
	
        tmpname = "thetaY_refconv_";
        tmpname += distname;
	TH1F * hy = (TH1F*)MyFile->Get(tmpname.c_str());
        tmpname = "thetaY_measdata";
        tmpname += distname;
	TH1F * h2 = (TH1F*)MyFile->Get(tmpname.c_str());
	h2->SetDirectory(0);
	hy->SetDirectory(0);
	//hy->Rebin(GetNbinsX/40);
	h2->Draw();
	hy->Draw();
	std::cout << "Chi^2 for thetaY_reco and " << distname << std::endl;
	h2->Chi2Test(hy,"P");
	
	delete hx;
	delete h1;
	delete hy;
	delete h2;
        MyFile->Close();

}

void MCSAnalysis::DoFFTDeconvolution(){
  
  // The basic methods associated with the RooUnfolding package
  // First generate a histogram of the measured data.
  TH1D* thetaX_data = 
    new TH1D("thetaX_data","Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaY_data = 
    new TH1D("thetaY_data","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaScatt_data = 
    new TH1D("thetaScatt_data","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  TH1D* theta2Scatt_data = 
    new TH1D("theta2Scatt_data","Scattering Angle between Momentum Vectors;#theta_{Scatt}^{2}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh2"], _histlimits["maxTh2"]);
  
  TH1D* thetaX_ref = 
    new TH1D("thetaX_ref","Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaY_ref = 
    new TH1D("thetaY_ref","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaScatt_ref = 
    new TH1D("thetaScatt_ref","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  TH1D* theta2Scatt_ref = 
    new TH1D("theta2Scatt_ref","Scattering Angle between Momentum Vectors;#theta_{Scatt}^{2}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

  /*
  TH1D* thetaX_fft = 
    new TH1D("thetaX_fft","Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaY_fft = 
    new TH1D("thetaY_fft","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  TH1D* thetaScatt_fft = 
    new TH1D("thetaScatt_fft","Scattering Angle between Momentum Vectors;#theta_{Scatt}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

  
  TH1D* fftX_data = 
    new TH1D("fftX_data","Change in Projected Angle (X);mode_{X}; Events per 4 mrad",400,0, 400);
  TH1D* fftY_data = 
    new TH1D("fftY_data","Change in Projected Angle (Y);mode_{Y}; Events per 4 mrad",400,0, 400);
  TH1D* fftScatt_dataRe = 
    new TH1D("fftScatt_dataRe","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0,400);
  TH1D* fftScatt_dataIm = 
    new TH1D("fftScatt_dataIm","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0,400);
  
  TH1D* fftX_ref = 
    new TH1D("fftX_ref","Change in Projected Angle (X);mode_{X}; Events per 4 mrad",400,0, 400);
  TH1D* fftY_ref = 
    new TH1D("fftY_ref","Change in Projected Angle (Y);mode_{Y}; Events per 4 mrad",400,0, 400);
  TH1D* fftScatt_refRe = 
    new TH1D("fftScatt_refRe","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0, 400);
  TH1D* fftScatt_refIm = 
    new TH1D("fftScatt_refIm","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400,0, 400);

  TH1D* fftX_div = 
    new TH1D("fftX_div","Change in Projected Angle (X);mode_{X}; Events per 4 mrad",400, 0, 400);
  TH1D* fftY_div = 
    new TH1D("fftY_div","Change in Projected Angle (Y);mode_{Y}; Events per 4 mrad",400, 0, 400);
  TH1D* fftScatt_divRe = 
    new TH1D("fftScatt_divRe","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400, 0, 400);
  TH1D* fftScatt_divIm = 
    new TH1D("fftScatt_divIm","Scattering Angle between Momentum Vectors;mode_{Scatt}; Events per mrad",400, 0, 400);
  */
  TH2D* thetaScatt_measdata_vp = 
    new TH2D("thetaScatt_measdata_vp","Scattering Angle between Momentum Vectors;Momentum (MeV/c); #theta_{Scatt}", 
	     400, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
  
  const Int_t NBINS = 19;
  Double_t scat_bin_array[NBINS + 1] = {-0.1151,-0.0938,-0.0754,-0.0597,-0.0463,-0.0347,-0.0248,-0.0162,-0.00895,-0.00269,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,0.0754,0.0938,0.1151};
  TH1D* scattering_proj_x = new TH1D("scattering_proj_x_DC","Change in Projected Angle (X);#Delta#theta_{X}; Events per radian", 
				     NBINS, scat_bin_array);
  TH1D* scattering_proj_y = new TH1D("scattering_proj_y_DC","Change in Projected Angle (Y);#Delta#theta_{Y}; Events per radian", 
				     NBINS, scat_bin_array);
  
  TH2D* projposUSDSdiff = 
    new TH2D("projposUSDSdiff_ref",
	     "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
	     180, -400, 400, 180, -400, 400);

  for(int i=0; i<_DSset.N(); i++){
    std::vector<double> projTheta = DefineProjectionAngles(_USset.E(i), _DSset.E(i));
    double thetaY = projTheta[1];  // atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
    double thetaX = projTheta[0];  // atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
    // double cosScatt = ( (1 + _USset.E(i).dXdz * _DSset.E(i).dXdz +
    //_USset.E(i).dYdz * _DSset.E(i).dYdz )/
    //			sqrt(1 + _USset.E(i).dXdz*_USset.E(i).dXdz +
    //			      _USset.E(i).dYdz*_USset.E(i).dYdz)/
    //			sqrt(1 + _DSset.E(i).dXdz*_DSset.E(i).dXdz +
    //			     _DSset.E(i).dYdz*_DSset.E(i).dYdz));
    double thetaScatt = projTheta[2];  /// acos(cosScatt);
    
    thetaX_data->Fill(thetaX);
    thetaY_data->Fill(thetaY);
    thetaScatt_data->Fill(thetaScatt);
    theta2Scatt_data->Fill(thetaScatt*thetaScatt);
    /*
      double dt0 = 7.64186 / 0.299792458; // ns. 
    double dt = _USset.E(i).TOF01;
    if (_USset.E(i).TOF12 < 99.0 * 8.22475 / 0.299792458){
      dt0 = 8.22475 / 0.299792458; // ns.
      dt = _USset.E(i).TOF12;
    }
    */
    double pz  = _USset.E(i).pz;
    // 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
    thetaScatt_measdata_vp->Fill(pz, thetaScatt);
    scattering_proj_x->Fill(thetaX);
    scattering_proj_y->Fill(thetaY);
  }

  for(int i=0; i<_DSMCset.N(); i++){
    
    double thetaY = atan(_DSMCset.E(i).dXdz) - atan(_USMCset.E(i).dXdz);
    double thetaX = atan(_DSMCset.E(i).dYdz) - atan(_USMCset.E(i).dYdz);
    double cosScatt = ( (1 + _USMCset.E(i).dXdz * _DSMCset.E(i).dXdz +
			 _USMCset.E(i).dYdz * _DSMCset.E(i).dYdz )/
			sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
			      _USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
			sqrt(1 + _DSMCset.E(i).dXdz*_DSMCset.E(i).dXdz +
			     _DSMCset.E(i).dYdz*_DSMCset.E(i).dYdz));
    double thetaScatt = acos(cosScatt);
    
    double xUSabs = _USMCset.E(i).X - _USMCset.E(i).dXdz * (_USMCset.E(i).Z - _sys["abspos"]);
    double yUSabs = _USMCset.E(i).Y - _USMCset.E(i).dYdz * (_USMCset.E(i).Z - _sys["abspos"]);
    double xDSabs = _DSMCset.E(i).X - _DSMCset.E(i).dXdz * (_DSMCset.E(i).Z - _sys["abspos"]);
    double yDSabs = _DSMCset.E(i).Y - _DSMCset.E(i).dYdz * (_DSMCset.E(i).Z - _sys["abspos"]);
    projposUSDSdiff->Fill(xDSabs - xUSabs, yDSabs - yUSabs);

    thetaX_ref->Fill(thetaX);
    thetaY_ref->Fill(thetaY);
    thetaScatt_ref->Fill(thetaScatt);
    theta2Scatt_ref->Fill(thetaScatt*thetaScatt);
  }
  
  /*
  int NXFFT=thetaX_data->Integral();
  thetaX_data->Sumw2();
  thetaX_data->Scale(1./NXFFT);
  thetaX_data->FFT(fftX_data,"Re R2R_3");
  int NYFFT=thetaY_data->Integral();
  thetaY_data->Sumw2();
  thetaY_data->Scale(1./NYFFT);
  thetaY_data->FFT(fftY_data,"Re R2R_3");
  int NScattFFT=thetaScatt_data->Integral();
  thetaScatt_data->Sumw2();
  thetaScatt_data->Scale(1./NScattFFT);
  thetaScatt_data->FFT(fftScatt_dataRe,"Re R2C P");
  thetaScatt_data->FFT(fftScatt_dataIm,"Im R2C P");
    
  for (int i=1; i<fftX_data->GetNbinsX() + 1; i++){
    double c = fftX_data->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftX_data->SetBinContent(i, fabs(c));
  }
  for (int i=1; i<fftY_data->GetNbinsX() + 1; i++){
    double c = fftY_data->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftY_data->SetBinContent(i, fabs(c));
  }
  // Supress modes higher than 10 if modes go to negative values.
  // bool iszero=false;
  //for (int i=1; i<fftScatt_data->GetNbinsX() + 1; i++){
  //  if(fftScatt_data->GetBinContent(i) < 1e-3)
      // iszero=true;
      // if(iszero)
  //    fftScatt_data->SetBinContent(i, 0.0);
  //}
  // Smooth the data
  
  TH1D* thetaScatt_smooth = (TH1D*)thetaScatt_data->Clone();
  thetaScatt_smooth->SetName("thetaScatt_smooth");
  fftScatt_data->FFT(thetaScatt_smooth,"Re R2R_7");
  thetaScatt_smooth->Scale(1./2./thetaScatt_data->GetNbinsX());
  thetaScatt_smooth->FFT(fftScatt_data,"Re R2R_7");
  
  NXFFT=thetaX_ref->Integral();
  thetaX_ref->Scale(1./NXFFT);
  thetaX_ref->FFT(fftX_ref,"Re R2R_3");
  NYFFT=thetaY_ref->Integral();
  thetaY_ref->Scale(1./NYFFT);
  thetaY_ref->FFT(fftY_ref,"Re R2R_3");
  NScattFFT=thetaScatt_ref->Integral();
  thetaScatt_ref->Scale(1./NScattFFT);
  thetaScatt_ref->FFT(fftScatt_refRe,"Re R2C P");
  thetaScatt_ref->FFT(fftScatt_refIm,"Im R2C P");
  

  // Correct the normalization after the FFT
  // fftX_ref->Scale(1./2./fftX_ref->GetNbinsX());
  // fftY_ref->Scale(1./2./fftY_ref->GetNbinsX());
  // fftScatt_ref->Scale(1./2./fftScatt_ref->GetNbinsX());
  
  for (int i=1; i<fftX_ref->GetNbinsX() + 1; i++){
    double c = fftX_ref->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftX_ref->SetBinContent(i, fabs(c));
  }
  for (int i=1; i<fftY_ref->GetNbinsX() + 1; i++){
    double c = fftY_ref->GetBinContent(i);
    if (i > 40 ) c = 0.;
    fftY_ref->SetBinContent(i, fabs(c));
  }
  
  // Supress modes higher than 10 if modes go to negative values.
  // iszero=false;
  
  double *re_full = new double[fftScatt_refRe->GetNbinsX()];
  double *im_full = new double[fftScatt_refRe->GetNbinsX()];
  for (int i=1; i<fftScatt_refRe->GetNbinsX() + 1; i++){
    // if(fftScatt_ref->GetBinContent(i) < 1e-3)
      // iszero=true;
      // if(iszero)
    //  fftScatt_data->SetBinContent(i, 0.0);
    double f0 = fftScatt_dataRe->GetBinContent(i);
    double g0 = fftScatt_dataIm->GetBinContent(i);
    double f1 = fftScatt_refRe->GetBinContent(i);
    double g1 = fftScatt_refIm->GetBinContent(i);
    double mag2 = f1*f1 + g1*g1;
    double rediv = (f1*f0 + g1*g0)/mag2;
    double imdiv = (g0*f1 - f0*g1)/mag2;
    if (mag2 < 1e-4 || fabs(f0) < 0.001 || fabs(g0) < 0.001 || fabs(f1) < 0.001 || fabs(g1) < 0.001){
      rediv = 0.0;
      imdiv = 0.0;
    }
    if ( i > binlimit && i < fftScatt_refRe->GetNbinsX()-binlimit ){	
      rediv = 0.0;
      imdiv = 0.0;
    }
    fftScatt_divRe->SetBinContent(i,rediv);
    fftScatt_divIm->SetBinContent(i,imdiv);
    re_full[i-1] = rediv;
    im_full[i-1] = imdiv;
  }
  int n = fftScatt_divRe->GetNbinsX();
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = 0;
  hb = TH1::TransformHisto(fft_back, hb, "Re");
  for (int i=1; i<fftScatt_refRe->GetNbinsX() + 1;i++){
    thetaScatt_fft->SetBinContent(i,hb->GetBinContent(i));
  }
  
  thetaScatt_smooth->SetName("thetaScatt_smooth");
  fftScatt_ref->FFT(thetaScatt_smooth,"Re R2R_7");
  thetaScatt_smooth->Scale(1./2./thetaScatt_ref->GetNbinsX());
  thetaScatt_smooth->FFT(fftScatt_ref,"Re R2R_7");
  
  // Divide the fft of the data by the fft of the reference
  fftX_div->Divide(fftX_data, fftX_ref);
  fftY_div->Divide(fftY_data, fftY_ref);

  // Undo the fft
  fftX_div->FFT(thetaX_fft,"Re R2R_3");
  fftY_div->FFT(thetaY_fft,"Re R2R_3");

  // Correct the normalization after the FFT
  thetaX_fft->Scale(1./2./fftX_data->GetNbinsX());
  thetaY_fft->Scale(1./2./fftY_data->GetNbinsX());
  thetaScatt_fft->Scale(1./2./thetaScatt_data->GetNbinsX());
  */
  outfile->cd();

  thetaX_data->Write();
  thetaY_data->Write();
  thetaScatt_data->Write();
  theta2Scatt_data->Write();
  thetaX_ref->Write();
  thetaY_ref->Write();
  thetaScatt_ref->Write();
  theta2Scatt_ref->Write();
  projposUSDSdiff->Write();
  /*
  thetaX_fft->Write();
  thetaY_fft->Write();
  thetaScatt_fft->Write();
 
  fftX_data->Write();
  fftY_data->Write();
  fftScatt_dataRe->Write();
  fftScatt_dataIm->Write();
  fftX_ref->Write();
  fftY_ref->Write();
  fftScatt_refRe->Write();
  fftScatt_refIm->Write();
  fftX_div->Write();
  fftY_div->Write();
  fftScatt_divRe->Write();
  fftScatt_divIm->Write();
  */
  delete projposUSDSdiff;
}

void MCSAnalysis::UpdateRunInfo(){
  _run.push_back(runnumber);
  _emit.push_back(_UStmpset.emittance(4, 1));
  _emitError.push_back(_emit.back() / sqrt(_UStmpset.N()));
  _meanX.push_back(_UStmpset.mean().X);  
  _meanXerror.push_back(sqrt(_UStmpset.cov()[0][0]));  
  _meanY.push_back(_UStmpset.mean().Y);  
  _meanYerror.push_back(sqrt(_UStmpset.cov()[2][2]));  
  _meandXdz.push_back(_UStmpset.mean().dXdz);  
  _meandXdzError.push_back(sqrt(_UStmpset.cov()[1][1]));  
  _meandYdz.push_back(_UStmpset.mean().dYdz);  
  _meandYdzError.push_back(sqrt(_UStmpset.cov()[3][3]));  
  
  _UStmpset.Set().clear();
  LastRunNumber = runnumber;
}

void MCSAnalysis::PlotRunInfo(){
  TGraphErrors* gemitt = new TGraphErrors(_run.size());
  gemitt->SetName("gemitt");
  gemitt->SetTitle(";Run Number;Emittance (mm)");
  TGraphErrors* gmeanX = new TGraphErrors(_run.size());
  gmeanX->SetName("gmeanX");
  gmeanX->SetTitle(";Run Number;Mean X (mm)");
  TGraphErrors* gmeanY = new TGraphErrors(_run.size());
  gmeanY->SetName("gmeanY");
  gmeanY->SetTitle(";Run Number;Mean Y (mm)");
  TGraphErrors* gmeandXdz = new TGraphErrors(_run.size());
  gmeandXdz->SetName("gmeandXdz");
  gmeandXdz->SetTitle(";Run Number;Mean dXdz (mm)");
  TGraphErrors* gmeandYdz = new TGraphErrors(_run.size());
  gmeandYdz->SetName("gmeandYdz");
  gmeandYdz->SetTitle(";Run Number;Mean dYdz (mm)");

  for (size_t q=0; q<_run.size(); q++){
    gemitt->SetPoint(q, _run.at(q), _emit.at(q));
    gemitt->SetPointError(q, 0, _emitError.at(q));
    gmeanX->SetPoint(q, _run.at(q), _meanX.at(q));
    gmeanX->SetPointError(q, 0, _meanXerror.at(q));
    gmeanY->SetPoint(q, _run.at(q), _meanY.at(q));
    gmeanY->SetPointError(q, 0, _meanYerror.at(q));
    gmeandXdz->SetPoint(q, _run.at(q), _meandXdz.at(q));
    gmeandXdz->SetPointError(q, 0, _meandXdzError.at(q));
    gmeandYdz->SetPoint(q, _run.at(q), _meandYdz.at(q));
    gmeandYdz->SetPointError(q, 0, _meandYdzError.at(q));
  }
  outfile->cd();
  gemitt->Write();
  gmeanX->Write();
  gmeanY->Write();
  gmeandXdz->Write();
  gmeandYdz->Write();

  _run.clear();
  _emit.clear();
  _meanX.clear();
  _meanY.clear();
  _meandXdz.clear();
  _meandYdz.clear();
  
}

bool MCSAnalysis::findVirtualPlanes(){
  bool planesfound = false;
  // The virtual planes should be the same event to event so once they
  // are found once they should be consistent
  if (USrefplaneI > 0 && USabsPlaneI >0 &&
      DSrefplaneI > 0 && DSabsPlaneI >0 ){
    if (USrefplaneI < mcevent->GetVirtualHits()->size() &&
	DSrefplaneI < mcevent->GetVirtualHits()->size() &&
	DSabsPlaneI < mcevent->GetVirtualHits()->size() &&
	USabsPlaneI < mcevent->GetVirtualHits()->size()){
      planesfound = true;
      return planesfound;
    }
  }
  if (1) { 
    double USmindiff = 9999, DSmindiff=9999;
    int USminI=-1, DSminI=-1;
    double USrefpos = _sys["abspos"] - 149.;
    double DSrefpos = _sys["abspos"] + 149.;
    for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
      double ztest = mcevent->GetVirtualHits()->at(j).GetPosition().z();
      double USzdiff = fabs(ztest - USrefpos);
      double DSzdiff = fabs(ztest - DSrefpos);
      // std::cout<<j<<"\t"<<ztest<<"\t"<<USzdiff<<"\t"<<DSzdiff<<std::endl;
      if(USzdiff < USmindiff){
	
	USmindiff = USzdiff;
	USminI = j;
      }
      if(DSzdiff < DSmindiff){
	//std::cout<<j<<"\t"<<ztest<<"\t"<<DSzdiff<<std::endl;
	DSmindiff = DSzdiff;
	DSminI = j;
      }
    }
    if (USmindiff < 150.){
      planesfound = true;
      USabsPlaneI = USminI;
      // std::cout<<"US absorber plane.\t"<<USrefpos<<"\t"<<USmindiff<<"\t"<<USabsPlaneI<<"\t"<<mcevent->GetVirtualHits()->at(USabsPlaneI).GetPosition().z()<<std::endl;
    }
    if (DSmindiff < 150. && planesfound==true){
      DSabsPlaneI = DSminI;
      // std::cout<<"DS absorber plane.\t"<<DSrefpos<<"\t"<<DSmindiff<<"\t"<<DSabsPlaneI<<"\t"<<mcevent->GetVirtualHits()->at(DSabsPlaneI).GetPosition().z()<<std::endl;
    }
    else 
      planesfound = false;
    
    // Now I will fill indicies corresponding to the tracker
    // reference planes
    USrefpos = USrefplaneZ;
    DSrefpos = DSrefplaneZ;
    
    USmindiff = 9999, DSmindiff=9999;
    USminI=-1, DSminI=-1;
    
    for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
      double ztest = mcevent->GetVirtualHits()->at(j).GetPosition().z();
      double USzdiff = fabs(ztest - USrefpos);
      if(USzdiff < USmindiff){
	USmindiff = USzdiff;
	USminI = j;
    }
      double DSzdiff = fabs(ztest - DSrefpos);
      if(DSzdiff < DSmindiff){
	DSmindiff = DSzdiff;
	DSminI = j;
      }
    }
    if (planesfound == true){
      if (USmindiff < 50. ){
	USrefplaneI = USminI;
	// std::cout<<"US reference plane.\t"<<USrefplaneZ<<"\t"<<USmindiff<<"\t"<<USrefplaneI<<"\t"<<mcevent->GetVirtualHits()->at(USrefplaneI).GetPosition().z()<<std::endl;
      }
      else  
	planesfound = false;
      if (DSmindiff < 50. ){
	DSrefplaneI = DSminI;
	// std::cout<<"DS reference plane.\t"<<DSrefplaneZ<<"\t"<<DSmindiff<<"\t"<<DSrefplaneI<<"\t"<<mcevent->GetVirtualHits()->at(DSrefplaneI).GetPosition().z()<<std::endl;
      }
      else 
	planesfound = false;
    }
    return planesfound;
  }
}
void MCSAnalysis::FillMuScattResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC){
  double thetaYMC = atan(DSMC.dXdz) - atan(USMC.dXdz);
  double thetaXMC = atan(DSMC.dYdz) - atan(USMC.dYdz);
  double XMC = DSMC.X;
  double YMC = DSMC.Y;
  theta_true_x_graph->Fill(thetaXMC);
  theta_true_y_graph->Fill(thetaYMC);
  if ( event_ok ){
    scattering_proj_y_R->Fill(XMC, thetaYMC);
    scattering_proj_x_R->Fill(YMC, thetaXMC);
    theta_meas_x_graph->Fill(thetaXMC);
    theta_meas_y_graph->Fill(thetaYMC);
  }
}

void MCSAnalysis::FillMCSResponse(bool event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC){
  
  double thetaYMC = atan(DSMC.dXdz) - atan(USMC.dXdz);
  double thetaXMC = atan(DSMC.dYdz) - atan(USMC.dYdz);
  double costhetaScattMC = ( (1 + USMC.dXdz * DSMC.dXdz + USMC.dYdz * DSMC.dYdz )/
			      sqrt(1 + USMC.dXdz*USMC.dXdz + USMC.dYdz*USMC.dYdz)/
			      sqrt(1 + DSMC.dXdz*DSMC.dXdz + DSMC.dYdz*DSMC.dYdz));
  double thetaScattMC = acos(costhetaScattMC);
  if(event_ok){
    double thetaY = atan(DS.dXdz) - atan(US.dXdz);
    double thetaX = atan(DS.dYdz) - atan(US.dYdz);
    double cosScatt = ( (1 + US.dXdz * DS.dXdz + US.dYdz * DS.dYdz )/
			sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz)/
			sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz));
    double thetaScatt = acos(cosScatt);
    resp_thetaX.Fill(thetaX, thetaXMC);
    resp_thetaY.Fill(thetaY, thetaYMC);
    resp_thetaScatt.Fill(thetaScatt, thetaScattMC);
    resp_theta2Scatt.Fill(thetaScatt, thetaScattMC * thetaScattMC);

    scattering_proj_x_R->Fill(DS.X - US.X, thetaYMC);
    scattering_proj_y_R->Fill(DS.Y - US.Y, thetaXMC);
    
    theta_true_x_graph->Fill(thetaXMC);
    theta_true_y_graph->Fill(thetaYMC);
  }
  else{
    resp_thetaX.Miss(thetaXMC);
    resp_thetaY.Miss(thetaYMC);
    resp_thetaScatt.Miss(thetaScattMC);
    resp_theta2Scatt.Miss(thetaScattMC * thetaScattMC);

    theta_true_x_graph->Fill(thetaXMC);
    theta_true_y_graph->Fill(thetaYMC);
  }
  
}
 
void MCSAnalysis::FillVarsVirtual(Vars &tmpvar, int j){
  
  tmpvar.X  = mcevent->GetVirtualHits()->at(j).GetPosition().x();
  tmpvar.Y  = mcevent->GetVirtualHits()->at(j).GetPosition().y();
  tmpvar.Z  = mcevent->GetVirtualHits()->at(j).GetPosition().z();
  tmpvar.dXdz = mcevent->GetVirtualHits()->at(j).GetMomentum().x()/
    mcevent->GetVirtualHits()->at(j).GetMomentum().z();
  tmpvar.dYdz = mcevent->GetVirtualHits()->at(j).GetMomentum().y()/
    mcevent->GetVirtualHits()->at(j).GetMomentum().z();
  double px = mcevent->GetVirtualHits()->at(j).GetMomentum().x();
  double py = mcevent->GetVirtualHits()->at(j).GetMomentum().y();
  double pz = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
  tmpvar.TOF12 = 27.5 * sqrt(1 + 105.65*105.65/pz/pz);
  tmpvar.TOF01 = 26.5 * sqrt(1 + 105.65*105.65/pz/pz);
  tmpvar.pz = pz;
  tmpvar.pid = mcevent->GetVirtualHits()->at(j).GetParticleId();
  tmpvar.isgood = true;
    // mcevent->GetVirtualHits()->at(j).GetMomentum().z();
  /*
  while(_SPsets.size() <= j){
    Ensemble tmpens;`
    _SPsets.push_back(tmpens);
  }
  _SPsets[j].append_instance(tmpvar); 
  */
}

void MCSAnalysis::FillCollectionSciFi(Collection& Set, int j, int k, double pz, int isDS, bool project){      

  if(j < int(scifievent->scifitracks().size()) && j != -1){
    if(k < int(scifievent->scifitracks()[j]->scifitrackpoints().size()) && k != -1){
      Vars tmpvar;
      FillVarsSciFi(tmpvar, j, k, pz, isDS);
      if (project){
	Vars newvar = PropagateVarsMu(tmpvar, _sys["abspos"]);
	tmpvar = newvar;
      }
      Set.append_instance(tmpvar);
    }
  }
  if( j == -1 || k == -1){
    Vars tmpvar;
    tmpvar.X = 0.0;
    tmpvar.Y = 0.0;
    tmpvar.Z = 0.0;
    tmpvar.dXdz = 1./2.;
    tmpvar.dYdz = 1./2.;
    tmpvar.px   = pz/2.;
    tmpvar.py   = pz/2.;
    tmpvar.pz   = pz;
    tmpvar.pid  = -13;
    tmpvar.isgood = false;

    if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0 &&
	int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) > 0)
      tmpvar.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
	- tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
    else
      tmpvar.TOF12 = 100.0 * 8.22475 / 0.299792458;
    if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0 &&
	int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0)
      tmpvar.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
	- tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
    else
      tmpvar.TOF01 = 100.0 * 7.64186 / 0.299792458;
    
    Set.append_instance(tmpvar);
  }
}

void MCSAnalysis::FillVarsSciFi(Vars& tmpvar, int j, int k, double pz, int isDS){
  double alX=0.0, alY=0.0, thX=0.0, thY=0.0;
  if(isDS==0){
    alX = _sys["alXUS"];
    alY = _sys["alYUS"];
    thX = _sys["thYUS"];
    thY = _sys["thXUS"];
  } else {
    alX = _sys["alXDS"];
    alY = _sys["alYDS"];
    thX = _sys["thYDS"];
    thY = _sys["thXDS"];
  }

  if(j < int(scifievent->scifitracks().size()) && j != -1){
    if(k < int(scifievent->scifitracks()[j]->scifitrackpoints().size()) && k != -1){
      tmpvar.X  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().x() + alX;
      tmpvar.Y  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().y() + alY;
      tmpvar.Z  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
      tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x() + tan(thX*atan(1.)/45.0);
      // scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().z();
      tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y() + tan(thY*atan(1.)/45.0);
      tmpvar.pz = pz * sqrt(1 + pow(tmpvar.dXdz, 2) + pow(tmpvar.dYdz, 2));
      tmpvar.px = scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().x() + tan(thX*atan(1.)/45.0)*tmpvar.pz;
      tmpvar.py = scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().y() + tan(thY*atan(1.)/45.0)*tmpvar.pz;
      tmpvar.pid = -13;
      tmpvar.isgood = true;
      // scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().z();
    }
  }
  if( j == -1 || k == -1){
    // Vars tmpvar;
    tmpvar.X = 0.0;
    tmpvar.Y = 0.0;
    tmpvar.Z = 0.0;
    tmpvar.dXdz = 1./2.;
    tmpvar.dYdz = 1./2.;
    tmpvar.px   = pz/2.;
    tmpvar.py   = pz/2.;
    tmpvar.pz   = pz;
    tmpvar.pid = -13;
    tmpvar.isgood = false;
  }
  if( tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size() > 0 &&
      tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size() > 0)
    tmpvar.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
      - tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
  else
    tmpvar.TOF12 = 100.0 * 8.22475 / 0.299792458;
  if( tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size() > 0 &&
      tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size() > 0)
    tmpvar.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
      - tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
  else
    tmpvar.TOF01 = 100.0 * 7.64186 / 0.299792458;
}

void MCSAnalysis::make_beam_histograms(Collection Set, std::string desc, std::string suffix){
  
  std::string tmptitle = desc + ";X (mm); Y (mm)";
  std::string tmpname  = suffix + "_XY";
  TH2D* XY = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -225, 225);
  tmptitle = desc + ";X (mm); dXdz";
  tmpname  = suffix + "_XdXdz";
  TH2D* XdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";X (mm); dYdz";
  tmpname  = suffix + "_XdYdz";
  TH2D* XdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";Y (mm); dXdz";
  tmpname  = suffix + "_YdXdz";
  TH2D* YdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";Y (mm); dYdz";
  tmpname  = suffix + "_YdYdz";
  TH2D* YdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
  tmptitle = desc + ";dXdz; dYdz";
  tmpname  = suffix + "_dXdzdYdz";
  TH2D* dXdzdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";X (mm); TOF12 (ns)";
  tmpname  = suffix + "_XTOF12";
  TH2D* XTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
  tmptitle = desc + ";Y (mm); TOF12 (ns)";
  tmpname  = suffix + "_YTOF12";
  TH2D* YTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
  tmptitle = desc + ";dXdz; TOF12 (ns)";
  tmpname  = suffix + "_dXdzTOF12";
  TH2D* dXdzTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
  tmptitle = desc + ";dYdz; TOF12 (ns)";
  tmpname  = suffix + "_dYdzTOF12";
  TH2D* dYdzTOF12 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
  
  for(int i=0; i<Set.N(); i++){
    if(Set.E(i).X==0 && Set.E(i).Y==0) continue;
    XY->Fill(Set.E(i).X,Set.E(i).Y);
    XdXdz->Fill(Set.E(i).X,Set.E(i).dXdz);
    XdYdz->Fill(Set.E(i).X,Set.E(i).dYdz);
    YdXdz->Fill(Set.E(i).Y,Set.E(i).dXdz);
    YdYdz->Fill(Set.E(i).Y,Set.E(i).dYdz);
    dXdzdYdz->Fill(Set.E(i).dXdz,Set.E(i).dYdz);
    XTOF12->Fill(Set.E(i).X,Set.E(i).TOF12);
    YTOF12->Fill(Set.E(i).Y,Set.E(i).TOF12);
    dXdzTOF12->Fill(Set.E(i).dXdz,Set.E(i).TOF12);
    dYdzTOF12->Fill(Set.E(i).dYdz,Set.E(i).TOF12);
  }
  
  XY->GetXaxis()->SetLabelSize(0.05);
  XY->GetXaxis()->SetTitleSize(0.05);
  XY->GetYaxis()->SetLabelSize(0.05);
  XY->GetYaxis()->SetTitleSize(0.05);
  XdXdz->GetXaxis()->SetLabelSize(0.05);
  XdXdz->GetXaxis()->SetTitleSize(0.05);
  XdXdz->GetYaxis()->SetLabelSize(0.05);
  XdXdz->GetYaxis()->SetTitleSize(0.05);
  XdYdz->GetXaxis()->SetLabelSize(0.05);
  XdYdz->GetXaxis()->SetTitleSize(0.05);
  XdYdz->GetYaxis()->SetLabelSize(0.05);
  XdYdz->GetYaxis()->SetTitleSize(0.05);
  YdXdz->GetXaxis()->SetLabelSize(0.05);
  YdXdz->GetXaxis()->SetTitleSize(0.05);
  YdXdz->GetYaxis()->SetLabelSize(0.05);
  YdXdz->GetYaxis()->SetTitleSize(0.05);
  YdYdz->GetXaxis()->SetLabelSize(0.05);
  YdYdz->GetXaxis()->SetTitleSize(0.05);
  YdYdz->GetYaxis()->SetLabelSize(0.05);
  YdYdz->GetYaxis()->SetTitleSize(0.05);
  dXdzdYdz->GetXaxis()->SetLabelSize(0.05);
  dXdzdYdz->GetXaxis()->SetTitleSize(0.05);
  dXdzdYdz->GetYaxis()->SetLabelSize(0.05);
  dXdzdYdz->GetYaxis()->SetTitleSize(0.05);
  XTOF12->GetXaxis()->SetLabelSize(0.05);
  XTOF12->GetXaxis()->SetTitleSize(0.05);
  XTOF12->GetYaxis()->SetLabelSize(0.05);
  XTOF12->GetYaxis()->SetTitleSize(0.05);
  YTOF12->GetXaxis()->SetLabelSize(0.05);
  YTOF12->GetXaxis()->SetTitleSize(0.05);
  YTOF12->GetYaxis()->SetLabelSize(0.05);
  YTOF12->GetYaxis()->SetTitleSize(0.05);
  dXdzTOF12->GetXaxis()->SetLabelSize(0.05);
  dXdzTOF12->GetXaxis()->SetTitleSize(0.05);
  dXdzTOF12->GetYaxis()->SetLabelSize(0.05);
  dXdzTOF12->GetYaxis()->SetTitleSize(0.05);
  dYdzTOF12->GetXaxis()->SetLabelSize(0.05);
  dYdzTOF12->GetXaxis()->SetTitleSize(0.05);
  dYdzTOF12->GetYaxis()->SetLabelSize(0.05);
  dYdzTOF12->GetYaxis()->SetTitleSize(0.05);

  outfile->cd();
  XY->Write();
  XdXdz->Write();
  XdYdz->Write();
  XTOF12->Write();
  YdXdz->Write();
  YdYdz->Write();
  YTOF12->Write();
  dXdzdYdz->Write();
  dXdzTOF12->Write();
  dYdzTOF12->Write();
  
  TCanvas* c1 = new TCanvas();
  XY->Draw("colz");
  std::string tmpfile = suffix + "_XY.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  XdXdz->Draw("colz");
  tmpfile = suffix + "_XdXdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  XdYdz->Draw("colz");
  tmpfile = suffix + "_XdYdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  YdXdz->Draw("colz");
  tmpfile = suffix + "_YdXdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  YdYdz->Draw("colz");
  tmpfile = suffix + "_YdYdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  dXdzdYdz->Draw("colz");
  tmpfile = suffix + "_dXdzdYdz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  XTOF12->Draw("colz");
  tmpfile = suffix + "_XTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  YTOF12->Draw("colz");
  tmpfile = suffix + "_YTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  dXdzTOF12->Draw("colz");
  tmpfile = suffix + "_dXdzTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  dYdzTOF12->Draw("colz");
  tmpfile = suffix + "_dYdzTOF12.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  
  
  delete c1;
  delete XY;
  delete XdXdz;
  delete XdYdz;
  delete YdXdz;
  delete YdYdz;
  delete XTOF12;
  delete YTOF12;
  delete dXdzdYdz;
  delete dXdzTOF12;
  delete dYdzTOF12;
}

void MCSAnalysis::make_acceptance_histograms(Collection USset, Collection DSset, 
					     std::string desc, std::string suffix){
  std::string tmptitle = desc + ";X (mm); Y (mm)";
  std::string tmpname  = suffix + "_posaccXY";
  TH2D* posaccXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_divaccXY";
  TH2D* divaccXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_posresXY";
  TH2D* posresXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_divresXY";
  TH2D* divresXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_posres2XY";
  TH2D* posres2XY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_divres2XY";
  TH2D* divres2XY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";X (mm); Y (mm)";
  tmpname  = suffix + "_posXY";
  TH2D* posXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posaccDXDY";
  TH2D* posaccDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_divaccDXDY";
  TH2D* divaccDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posresDXDY";
  TH2D* posresDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_divresDXDY";
  TH2D* divresDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posres2DXDY";
  TH2D* posres2DXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_divres2DXDY";
  TH2D* divres2DXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
  tmptitle = desc + ";dX/dz ; dY/dz";
  tmpname  = suffix + "_posDXDY";
  TH2D* posDXDY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);

  tmptitle = desc + ";|#vec{r}^{proj#to DS}_{US} - #vec{r}^{meas}_{DS}";
  tmpname  = suffix + "_appRef";
  TH1D* appRef = new TH1D(tmpname.c_str(), tmptitle.c_str(), 200, 0, 200);
  tmptitle = desc + ";X_{US}^{proj#to DS} - X_{DS}^{meas};Y_{US}^{proj#to DS} - Y_{DS}^{meas}";
  tmpname  = suffix + "_projRefDiff";
  TH2D* projRefDiff = new TH2D(tmpname.c_str(), tmptitle.c_str(), 200, -100, 100, 200, -100, 100);

  double latchZ = 0.0;
  Vars DSproj;
  for (int i=0; i<USset.N(); i++){
    if (DSset.E(i).Z != 0 && latchZ == 0.0){
      latchZ = DSset.E(i).Z;
    } 
    if (latchZ != 0.0 && fabs(DSset.E(i).Z - latchZ) < 5){
      DSproj = PropagateVarsMu(USset.E(i), DSset.E(i).Z);
    } else if(latchZ != 0.0){
      DSproj = PropagateVarsMu(USset.E(i), latchZ);
    } else {
      continue;
    }
    
    Vars predDiff = DSproj - DSset.E(i);
    double approach = sqrt(predDiff.X * predDiff.X + 
			   predDiff.Y * predDiff.Y +
			   predDiff.Z * predDiff.Z);
    appRef->Fill(approach); // can be very large if the ds track is lost
    projRefDiff->Fill(predDiff.X, predDiff.Y);
    posXY->Fill(DSproj.X, DSproj.Y);
    posDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
    
    if ( DSset.E(i).Z > USset.E(i).Z ){ // lost ds tracks are set to z = 0
      posaccXY->Fill(DSproj.X, DSproj.Y);
      posaccDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
      if ( approach < 36.) { // 36 mm is approximately the resolution of the projection
	posresXY->Fill(DSproj.X, DSproj.Y);
	posresDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
      }
      if ( approach < 2.*36.) { // 36 mm is approximately the resolution of the projection
	posres2XY->Fill(DSproj.X, DSproj.Y);
	posres2DXDY->Fill(DSproj.dXdz, DSproj.dYdz);
      }
    }
  }
  posXY->Sumw2();
  posaccXY->Sumw2();
  divaccXY->Divide(posaccXY, posXY);
  posresXY->Sumw2();
  divresXY->Divide(posresXY, posXY);
  posres2XY->Sumw2();
  divres2XY->Divide(posres2XY, posXY);
  posaccDXDY->Sumw2();
  divaccDXDY->Divide(posaccDXDY, posDXDY);
  posresDXDY->Sumw2();
  divresDXDY->Divide(posresDXDY, posDXDY);
  posres2DXDY->Sumw2();
  divres2DXDY->Divide(posres2DXDY, posDXDY);

  posXY->GetXaxis()->SetLabelSize(0.05);
  posXY->GetXaxis()->SetTitleSize(0.05);
  posXY->GetYaxis()->SetLabelSize(0.05);
  posXY->GetYaxis()->SetTitleSize(0.05);
  posDXDY->GetXaxis()->SetLabelSize(0.05);
  posDXDY->GetXaxis()->SetTitleSize(0.05);
  posDXDY->GetYaxis()->SetLabelSize(0.05);
  posDXDY->GetYaxis()->SetTitleSize(0.05);
  posaccXY->GetXaxis()->SetLabelSize(0.05);
  posaccXY->GetXaxis()->SetTitleSize(0.05);
  posaccXY->GetYaxis()->SetLabelSize(0.05);
  posaccXY->GetYaxis()->SetTitleSize(0.05);
  divaccXY->GetXaxis()->SetLabelSize(0.05);
  divaccXY->GetXaxis()->SetTitleSize(0.05);
  divaccXY->GetYaxis()->SetLabelSize(0.05);
  divaccXY->GetYaxis()->SetTitleSize(0.05);
  posresXY->GetXaxis()->SetLabelSize(0.05);
  posresXY->GetXaxis()->SetTitleSize(0.05);
  posresXY->GetYaxis()->SetLabelSize(0.05);
  posresXY->GetYaxis()->SetTitleSize(0.05);
  divresXY->GetXaxis()->SetLabelSize(0.05);
  divresXY->GetXaxis()->SetTitleSize(0.05);
  divresXY->GetYaxis()->SetLabelSize(0.05);
  divresXY->GetYaxis()->SetTitleSize(0.05);
  posres2XY->GetXaxis()->SetLabelSize(0.05);
  posres2XY->GetXaxis()->SetTitleSize(0.05);
  posres2XY->GetYaxis()->SetLabelSize(0.05);
  posres2XY->GetYaxis()->SetTitleSize(0.05);
  divres2XY->GetXaxis()->SetLabelSize(0.05);
  divres2XY->GetXaxis()->SetTitleSize(0.05);
  divres2XY->GetYaxis()->SetLabelSize(0.05);
  divres2XY->GetYaxis()->SetTitleSize(0.05);
  posaccDXDY->GetXaxis()->SetLabelSize(0.05);
  posaccDXDY->GetXaxis()->SetTitleSize(0.05);
  posaccDXDY->GetYaxis()->SetLabelSize(0.05);
  posaccDXDY->GetYaxis()->SetTitleSize(0.05);
  divaccDXDY->GetXaxis()->SetLabelSize(0.05);
  divaccDXDY->GetXaxis()->SetTitleSize(0.05);
  divaccDXDY->GetYaxis()->SetLabelSize(0.05);
  divaccDXDY->GetYaxis()->SetTitleSize(0.05);
  posresDXDY->GetXaxis()->SetLabelSize(0.05);
  posresDXDY->GetXaxis()->SetTitleSize(0.05);
  posresDXDY->GetYaxis()->SetLabelSize(0.05);
  posresDXDY->GetYaxis()->SetTitleSize(0.05);
  divresDXDY->GetXaxis()->SetLabelSize(0.05);
  divresDXDY->GetXaxis()->SetTitleSize(0.05);
  divresDXDY->GetYaxis()->SetLabelSize(0.05);
  divresDXDY->GetYaxis()->SetTitleSize(0.05);
  posres2DXDY->GetXaxis()->SetLabelSize(0.05);
  posres2DXDY->GetXaxis()->SetTitleSize(0.05);
  posres2DXDY->GetYaxis()->SetLabelSize(0.05);
  posres2DXDY->GetYaxis()->SetTitleSize(0.05);
  divres2DXDY->GetXaxis()->SetLabelSize(0.05);
  divres2DXDY->GetXaxis()->SetTitleSize(0.05);
  divres2DXDY->GetYaxis()->SetLabelSize(0.05);
  divres2DXDY->GetYaxis()->SetTitleSize(0.05);

  outfile->cd();
  posXY->Write();
  posaccXY->Write();
  posresXY->Write();
  divres2XY->Write();
  divaccXY->Write();
  divresXY->Write();
  divres2XY->Write();
  posDXDY->Write();
  posaccDXDY->Write();
  posresDXDY->Write();
  divres2DXDY->Write();
  divaccDXDY->Write();
  divresDXDY->Write();
  divres2DXDY->Write();
  appRef->Write();
  projRefDiff->Write();

  delete posXY;
  delete posDXDY;
  delete posaccXY;
  delete posaccDXDY;
  delete posresXY;
  delete posresDXDY;
  delete posres2XY;
  delete posres2DXDY;
  delete divaccXY;
  delete divaccDXDY;
  delete divresXY;
  delete divresDXDY;
  delete divres2XY;
  delete divres2DXDY;
  delete appRef;
  delete projRefDiff;
}

void MCSAnalysis::make_scattering_acceptance_histograms(Collection USset,
							Collection DSset,
							Collection DSrefset,
							std::string desc,
							std::string suffix){
  
  std::string tmptitle = desc + ";#theta_Y (radians)";
  std::string tmpname = "thetaY_all_";
  tmpname += suffix;
  TH1D* thetaY_all = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_all_";
  tmpname += suffix;
  TH1D* thetaX_all = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_acc_";
  tmpname += suffix;
  TH1D* thetaY_acc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_acc_";
  tmpname += suffix;
  TH1D* thetaX_acc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_divacc_";
  tmpname += suffix;
  TH1D* thetaY_divacc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_divacc_";
  tmpname += suffix;
  TH1D* thetaX_divacc = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_res_";
  tmpname += suffix;
  TH1D* thetaY_res = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_res_";
  tmpname += suffix;
  TH1D* thetaX_res = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_divres_";
  tmpname += suffix;
  TH1D* thetaY_divres = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_divres_";
  tmpname += suffix;
  TH1D* thetaX_divres = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_res2_";
  tmpname += suffix;
  TH1D* thetaY_res2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_res2_";
  tmpname += suffix;
  TH1D* thetaX_res2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

  tmptitle = desc + ";#theta_Y (radians)";
  tmpname = "thetaY_divres2_";
  tmpname += suffix;
  TH1D* thetaY_divres2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#Delta#theta_{Y}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  tmptitle = desc + ";#theta_X (radians)";
  tmpname = "thetaX_divres2_";
  tmpname += suffix;
  TH1D* thetaX_divres2 = 
    new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#Delta#theta_{X}; Events per 4 mrad",
	     _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
  
  for (int i=0; i<USset.N(); i++){
    if (i >= DSrefset.N()) 
      break;
    std::vector<double> projTheta = DefineProjectionAngles(USset.E(i), DSset.E(i));
    double thetaX = projTheta[0];
    double thetaY = projTheta[1];
    // double thetaScatt = projTheta[2];
    Vars DSdiff = DSset.E(i) - DSrefset.E(i);
    double approach = sqrt(DSdiff.X * DSdiff.X + 
			   DSdiff.Y * DSdiff.Y +
			   DSdiff.Z * DSdiff.Z);
    thetaX_all->Fill(thetaX);
    thetaY_all->Fill(thetaY);
    if (DSrefset.E(i).Z != 0){
      thetaX_acc->Fill(thetaX);
      thetaY_acc->Fill(thetaY);
    
      if (approach < 36.){
	thetaX_res->Fill(thetaX);
	thetaY_res->Fill(thetaY);
      }
      if (approach < 2* 36.){
	thetaX_res2->Fill(thetaX);
	thetaY_res2->Fill(thetaY);
      }
    }
  }
  thetaX_all->Sumw2();
  thetaY_all->Sumw2();
  thetaX_acc->Sumw2();
  thetaY_acc->Sumw2();
  thetaX_divacc->Divide(thetaX_acc, thetaX_all);
  thetaY_divacc->Divide(thetaY_acc, thetaY_all);
  thetaX_res->Sumw2();
  thetaY_res->Sumw2();
  thetaX_divres->Divide(thetaX_res, thetaX_all);
  thetaY_divres->Divide(thetaY_res, thetaY_all);
  thetaX_res2->Sumw2();
  thetaY_res2->Sumw2();
  thetaX_divres2->Divide(thetaX_res2, thetaX_all);
  thetaY_divres2->Divide(thetaY_res2, thetaY_all);

  thetaX_all->GetXaxis()->SetLabelSize(0.05);
  thetaX_all->GetXaxis()->SetTitleSize(0.05);
  thetaX_all->GetYaxis()->SetLabelSize(0.05);
  thetaX_all->GetYaxis()->SetTitleSize(0.05);
  thetaX_acc->GetXaxis()->SetLabelSize(0.05);
  thetaX_acc->GetXaxis()->SetTitleSize(0.05);
  thetaX_acc->GetYaxis()->SetLabelSize(0.05);
  thetaX_acc->GetYaxis()->SetTitleSize(0.05);
  thetaX_res->GetXaxis()->SetLabelSize(0.05);
  thetaX_res->GetXaxis()->SetTitleSize(0.05);
  thetaX_res->GetYaxis()->SetLabelSize(0.05);
  thetaX_res->GetYaxis()->SetTitleSize(0.05);
  thetaX_res2->GetXaxis()->SetLabelSize(0.05);
  thetaX_res2->GetXaxis()->SetTitleSize(0.05);
  thetaX_res2->GetYaxis()->SetLabelSize(0.05);
  thetaX_res2->GetYaxis()->SetTitleSize(0.05);

  thetaX_divacc->GetXaxis()->SetLabelSize(0.05);
  thetaX_divacc->GetXaxis()->SetTitleSize(0.05);
  thetaX_divacc->GetYaxis()->SetLabelSize(0.05);
  thetaX_divacc->GetYaxis()->SetTitleSize(0.05);
  thetaX_divres->GetXaxis()->SetLabelSize(0.05);
  thetaX_divres->GetXaxis()->SetTitleSize(0.05);
  thetaX_divres->GetYaxis()->SetLabelSize(0.05);
  thetaX_divres->GetYaxis()->SetTitleSize(0.05);
  thetaX_divres2->GetXaxis()->SetLabelSize(0.05);
  thetaX_divres2->GetXaxis()->SetTitleSize(0.05);
  thetaX_divres2->GetYaxis()->SetLabelSize(0.05);
  thetaX_divres2->GetYaxis()->SetTitleSize(0.05);


  thetaY_all->GetXaxis()->SetLabelSize(0.05);
  thetaY_all->GetXaxis()->SetTitleSize(0.05);
  thetaY_all->GetYaxis()->SetLabelSize(0.05);
  thetaY_all->GetYaxis()->SetTitleSize(0.05);
  thetaY_acc->GetXaxis()->SetLabelSize(0.05);
  thetaY_acc->GetXaxis()->SetTitleSize(0.05);
  thetaY_acc->GetYaxis()->SetLabelSize(0.05);
  thetaY_acc->GetYaxis()->SetTitleSize(0.05);
  thetaY_res->GetXaxis()->SetLabelSize(0.05);
  thetaY_res->GetXaxis()->SetTitleSize(0.05);
  thetaY_res->GetYaxis()->SetLabelSize(0.05);
  thetaY_res->GetYaxis()->SetTitleSize(0.05);
  thetaY_res2->GetXaxis()->SetLabelSize(0.05);
  thetaY_res2->GetXaxis()->SetTitleSize(0.05);
  thetaY_res2->GetYaxis()->SetLabelSize(0.05);
  thetaY_res2->GetYaxis()->SetTitleSize(0.05);

  thetaY_divacc->GetXaxis()->SetLabelSize(0.05);
  thetaY_divacc->GetXaxis()->SetTitleSize(0.05);
  thetaY_divacc->GetYaxis()->SetLabelSize(0.05);
  thetaY_divacc->GetYaxis()->SetTitleSize(0.05);
  thetaY_divres->GetXaxis()->SetLabelSize(0.05);
  thetaY_divres->GetXaxis()->SetTitleSize(0.05);
  thetaY_divres->GetYaxis()->SetLabelSize(0.05);
  thetaY_divres->GetYaxis()->SetTitleSize(0.05);
  thetaY_divres2->GetXaxis()->SetLabelSize(0.05);
  thetaY_divres2->GetXaxis()->SetTitleSize(0.05);
  thetaY_divres2->GetYaxis()->SetLabelSize(0.05);
  thetaY_divres2->GetYaxis()->SetTitleSize(0.05);

  outfile->cd();
  thetaX_all->Write();
  thetaY_all->Write();
  thetaX_acc->Write();
  thetaY_acc->Write();
  thetaX_res->Write();
  thetaY_res->Write();
  thetaX_res2->Write();
  thetaY_res2->Write();

  thetaX_divacc->Write();
  thetaY_divacc->Write();
  thetaX_divres->Write();
  thetaY_divres->Write();
  thetaX_divres2->Write();
  thetaY_divres2->Write();

  delete thetaX_all;
  delete thetaY_all;
  delete thetaX_acc;
  delete thetaY_acc;
  delete thetaX_res;
  delete thetaY_res;
  delete thetaX_res2;
  delete thetaY_res2;
  delete thetaX_divacc;
  delete thetaY_divacc;
  delete thetaX_divres;
  delete thetaY_divres;
  delete thetaX_divres2;
  delete thetaY_divres2;
}

Json::Value MCSAnalysis::SetupConfig(int verbose_level) {
  std::cerr << "Running with verbose level " << verbose_level << std::endl;
  Json::Value config(Json::objectValue);
  config["maximum_module_depth"] = 50;
  config["check_volume_overlaps"] = true;
  config["reconstruction_geometry_filename"] = geometryfile;
  config["simulation_geometry_filename"] = geometryfile;
  config["simulation_geometry_debug"] = false;
  config["maximum_number_of_steps"] = 10000;
  config["will_do_stack_trace"] = true;
  config["keep_tracks"] = true;
  config["keep_steps"] = true;
  config["verbose_level"] = verbose_level;
  config["geant4_visualisation"] = false;
  config["physics_model"] = "QGSP_BERT";
  config["reference_physics_processes"] = "mean_energy_loss";
  config["physics_processes"] = "standard";
  config["particle_decay"] = true;
  config["polarised_decay"] = true;
  config["charged_pion_half_life"] = -1.;
  config["muon_half_life"] = -1.;
  config["production_threshold"] = 0.5;
  config["fine_grained_production_threshold"] = Json::Value(Json::objectValue);
  config["default_keep_or_kill"] = true;
  config["spin_tracking"] = true;
  config["keep_or_kill_particles"] = "{\"neutron\":False}";
  config["kinetic_energy_threshold"] = 0.1;
  config["max_step_length"] = 100.;
  config["max_track_time"] = 1.e9;
  config["max_track_length"] = 1.e8;
  config["simulation_reference_particle"] = 
    JsonWrapper::StringToJson(
			      std::string("{\"position\":{\"x\":0.0,\"y\":-0.0,\"z\":-5500.0},")+
			      std::string("\"spin\":{\"x\":0.0,\"y\":-0.0,\"z\":1.0},")+
			      std::string("\"momentum\":{\"x\":0.0,\"y\":0.0,\"z\":1.0},")+
			      std::string("\"particle_id\":-13,\"energy\":226.0,\"time\":0.0,")+
			      std::string("\"random_seed\":10}")
			      );
  config["stepping_algorithm"] = "ClassicalRK4";
  config["delta_one_step"] = -1.;
  config["delta_intersection"] = -1.;
  config["epsilon_min"] = -1.;
  config["epsilon_max"] = -1.;
  config["miss_distance"] = -1.;
  config["everything_special_virtual"] = false;
  config["field_tracker_absolute_error"] = 1.e-4;
  config["field_tracker_relative_error"] = 1.e-4;
  config["data_maximum_reference_count"] = 200;
  return config;
}

Vars MCSAnalysis::PropagateVarsMu(Vars event, double z0){
  
  
  double mass = 105.658; // muon mass [MeV/c^2]
  double time = 0.; // time in lab frame [ns]
  double x = event.X; // horizontal position [mm]
  double y = event.Y; // vertical position [mm]
  double z = event.Z; // longitudinal position [mm]
  double px = event.px; // horizontal momentum component [MeV/c]
  double py = event.py; // vertical momentum component [MeV/c]
  double pz = event.pz;
  
  double energy = std::sqrt(pz*pz + px*px + py*py + mass*mass); // Total energy [MeV]
  
  double event_vector[8] = { time, x, y, z, energy, px, py, pz };
  /*
  BTField* field = dynamic_cast<BTField*>(MAUS::Globals::GetMCFieldConstructor());
  try {
    MAUS::GlobalTools::propagate(event_vector, z0, field, 10., 
				 MAUS::DataStructure::Global::kMuPlus, true);
  } catch (...){
  */
  // Assume a straight track
  event_vector[1] = px/pz * (z0 - z) + x;
  event_vector[2] = py/pz * (z0 - z) + y;
  event_vector[3] = z0;
  event_vector[5] = px;
  event_vector[6] = py;
  event_vector[7] = pz;
  // }
  Vars prop;
  prop.X = event_vector[1];
  prop.Y = event_vector[2];
  prop.Z = event_vector[3];
  prop.px = event_vector[5];
  prop.py = event_vector[6];
  prop.pz = event_vector[7];
  prop.dXdz = event_vector[5]/event_vector[7];
  prop.dYdz = event_vector[6]/event_vector[7];
  prop.TOF12 = event.TOF12;
  prop.TOF01 = event.TOF01;

  return prop;
}

