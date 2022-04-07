#include "MCSAnalysis.h"

#include "TRint.h"
#include "TROOT.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using boost::property_tree::ptree;
using boost::property_tree::write_xml;
using boost::property_tree::xml_writer_settings;

MCSAnalysis::MCSAnalysis(std::string tree, std::string mctree, std::string reftree, std::string outname, std::map<std::string, double> histlimits)
    //: p_vec(TVectorD(19)), res(TVectorD(30)), pStart_vec(TVectorD(19)),
    //    pStart_vec_y(TVectorD(19)), theta_true_x(TVectorD(19)), theta_true_y(TVectorD(19))
{
    chain   = new TChain(tree.c_str());
    refchain   = new TChain(reftree.c_str());
    mcchain = new TChain(mctree.c_str());
    mcemptychain = new TChain(mctree.c_str());
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
    _sys["USref"]  = 15062;
    _sys["DSref"]  = 18849;
    _sys["niter"]   = 10;
    _sys["abspos"]  = 16952.5;
    _sys["diffpos"]  = 13620;
    _sys["DStrkr5"]  = 19948.8;
    _sys["Nevents"]  = -1.;
    _sys["FracEvents"] = 1.;

    dte01data = 25.4042;
    dte12data = 27.505;
    dte01MC = 25.554;
    /* dte12MC = 27.4543; */
    dte12MC = dte12data;
    weight172 = 0;
    weight200 = 0;
    weight240 = 0;
    isEmpty=false;
    modelname2 = "Cobb";
    modelname1 = "GEANT";
    modelname3 = "Moliere";
    tracker0 = 0;
    tracker1 = 1;

    double stops[9] = {0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
    double red[9]   = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
    double green[9] = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
    double blue[9]  = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
    double ncontours = 255;
    TColor::CreateGradientColorTable(9, stops, red, green, blue, ncontours);
    gStyle->SetNumberContours(ncontours);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetNdivisions(505, "x");
    gStyle->SetLabelSize(0.05, "x");
    gStyle->SetTitleSize(0.06, "x");
    gStyle->SetTitleOffset(1.10, "x");
    gStyle->SetNdivisions(505, "y");
    gStyle->SetLabelSize(0.05, "y");
    gStyle->SetTitleSize(0.06, "y");
    gStyle->SetTitleOffset(1.10, "y");

    _histlimits["NbinsXY"] = histlimits.count("NbinsXY") != 0 ? histlimits["NbinsXY"]: 200.0;
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
    TOF_lower_limit_ref = 27;
    TOF_upper_limit_ref = 28;

    DSentries = 0;
    DSrefentries = 0;
    DScor = 0;
    DScor = 0;
    difradius = 90;
    meanp = 300.0;
    sigmap = 0.035;
    binlimit = 22;
    semom = 0;
    meanmom = 0;
    rmsmom = 0;
    errrmsmom = 0;
    bwmom = 0;
    tX_rG_int = 0;

    // histograms of selection variables

    chi2pern = new TH1D("chi2n","chi2n", 40, 0, 10);
    diffradius = new TH1D("diffradius","diffradius", 500, 0, 500);
    projradius = new TH1D("projradius","projradius", 500, 0, 500);
    nocutprojradius = new TH1D("nocutprojradius","nocutprojradius", 500, 0, 500);
    trackno = new TH2D("trackno","trackno",5,0,4,5,0,4);
    eventnovspz = new TH2D("eventnovspz","eventnovspz",50,150,250,60000,0,59999);
    tofhitno = new TH2D("tofhitno","tofhitno",5,0,4,5,0,4);
    pathlengthabs = new TH1D("pathlengthabs","pathlengthabs", 500, 0, 500);
    t_cor = new TH2D("t_cor","t_cor", 1000, 0, 1000,100, 8.05, 9.15);
    TOF0Energy = new TH1D("TOF0Energy","TOF0Energy", 200, 150, 350);
    TOF1Energy = new TH1D("TOF1Energy","TOF1Energy", 200, 150, 350);
    TOF2Energy = new TH1D("TOF2Energy","TOF2Energy", 200, 150, 350);
    difEloss = new TH1D("difEloss","difEloss", 10, 0, 10);
    energylossproj = new TH1D("energylossproj","energylossproj", 40, 5000, 24000);
    //energylosspro = new TH1D("energylosspro","energylosspro", 40, 5000, 24000);
    residual = new TH1D("residual","residual", 100, -25, 25);
    residual12 = new TH1D("residual12","residual12", 100, -50, 50);
    residual01 = new TH1D("residual01","residual01", 100, -50, 50);
    residualUD = new TH1D("residualUD","residualUD", 100, -50, 50);
    refresidual = new TH1D("refresidual","refresidual", 100, -50, 50);
    refresidual01 = new TH1D("refresidual01","refresidual01", 100, -50, 50);
    refresidual12 = new TH1D("refresidual12","refresidual12", 100, -50, 50);
    refresidualUD = new TH1D("refresidualUD","refresidualUD", 100, -50, 50);
    pzdEdx = new TH2D("pzdEdx","pzdEdz", 100, 1, 100,250, -500, 500);
    TOFcom = new TH2D("TOF01 vs TOF12","TOF01 vs TOF12", 33, 100, 350,33, 100, 350);
    energyloss = new TH2D("energyloss","energyloss", 40, 5000, 24000,350, 100, 450);
    TOFvsMCTruth = new TH2D("TOFvsMCTruth","", 200, 150, 250,200, 150, 250);
    refTOFvsMCTruth = new TH2D("refTOFvsMCTruth","", 14, 150, 250,14, 150, 250);
    MCTruth = new TH1D("MCTruth","MCTruth", 33, 100, 350);
    refMCTruth = new TH1D("refMCTruth","refMCTruth",33, 100, 350);
    TOF01vsTOF12 = new TH2D("TOF01vsTOF12","TOF01vsTOF12", 40, 100, 300,40, 100, 300);
    refTOF01vsTOF12 = new TH2D("refTOF01vsTOF12","refTOF01vsTOF12", 40, 100, 300,40, 100, 300);
    TOF01vsMCTruth = new TH2D("TOF01vsMCTruth","TOF01vsMCTruth", 40, 100, 300, 40, 100, 300);
    TOF12longPaul6thforvsMCTruth = new TH2D("TOF12longforupvsMCTruth","TOF12longforupvsMCTruth", 33, 100, 350, 33, 100, 350);
    TOF01Paul6thforvsMCTruth = new TH2D("TOF01forupvsMCTruth","TOF01forupvsMCTruth", 33, 100, 350, 33, 100, 350);
    TOF01shortPaul6thforvsMCTruth = new TH2D("TOF01shortforupvsMCTruth","TOF01shortforupvsMCTruth", 33, 100, 350, 33, 100, 350);
    TOF01fordownvsMCTruth = new TH2D("TOF01fordownvsMCTruth","TOF01fordownvsMCTruth", 33, 100, 350,33, 100, 350);
    TOF12fordownvsMCTruth = new TH2D("TOF12fordownvsMCTruth","TOF12fordownvsMCTruth", 33, 100, 350,33, 100, 350);
    TOF12cobbvsMCTruth = new TH2D("TOF12cobbvsMCTruth","TOF12cobbvsMCTruth", 33, 100, 350,33, 100, 350);
    TOF01forupvsTOF01fordown = new TH2D("TOF01forupvsTOF01fordown","TOF01forupvsTOF01fordown", 33, 100, 350,33, 100, 350);
    tof10 = new TH1D("tof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
    tof10_sel = new TH1D("tof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
    tof21 = new TH1D("tof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
    tof21_sel = new TH1D("tof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
    mctof10 = new TH1D("mctof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
    mctof10_sel = new TH1D("mctof10_sel","TOF Selection Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
    mctof21 = new TH1D("mctof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
    mctof21_sel = new TH1D("mctof21_sel","TOF Selection Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
    rawtime = new TH1D("rawtime","TOF Between Stations 0 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 40);
    rawtime12 = new TH1D("rawtime21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
    noweighttof10 = new TH1D("noweighttof10","TOF Between Stations 1 and 0; t_{TOF1} - t_{TOF0} (ns)", 150, 10, 40);
    noweighttof21 = new TH1D("noweighttof21","TOF Between Stations 2 and 1; t_{TOF2} - t_{TOF1} (ns)", 150, 0, 50);
    calc_mom = new TH1D("calc_mom","Momentum Calculated from TOF; Momentum (MeV/c)", 100, 0, 400);
    cor_mom = new TH1D("LiH_mom","; #it{p} (MeV/c); Normalized Entries", 180, 120, 300);
    uncor_mom = new TH1D("uncor_mom","; #it{p} (MeV/c); Normalized Entries", 180, 120, 300);
    uncor_momsel = new TH1D("uncor_momsel","Cor Momentum Calculated from TOF; Momentum (MeV/c)", 180, 120, 300);
    mccalc_mom = new TH1D("empty_mom","; #it{p} (MeV/c); Normalized Entries", 100, 0, 400);
    mctrue_mom = new TH1D("mctrue_mom","Momentum from Virtual Planes; Momentum (MeV/c)", 400, 0, 400);
    cuts_accept   = new TH1D("cuts_accept", 
            ";Selection Criteria; Surviving Events", 6, 0, 6);
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

    sym_theta_true_x_graph =
        new TH1D("thetaX_graphsym","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    sym_theta_true_y_graph =
        new TH1D("thetaY_graphsym","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    theta_true_x_graph =
        new TH1D("thetaX_graph","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    theta_true_y_graph =
        new TH1D("thetaY_graph","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    minustheta_true_x_graph =
        new TH1D("minusthetaX_graph","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    minustheta_true_y_graph =
        new TH1D("minusthetaY_graph","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    theta_true_scat_graph = 
        new TH1D("thetaScatt_graph","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",
                _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
    theta_true_scat2_graph = 
        new TH1D("theta2Scatt_graph","Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",
                _histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

    mctheta_true_x_graph =
        new TH1D("thetaX_graphmc","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    mctheta_true_y_graph =
        new TH1D("thetaY_graphmc","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    mctheta_true_scat_graph = 
        new TH1D("theta_true_scat_graphmc","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",
                _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
    scattering_proj_x_resp =
        new TH2D("scattering_proj_x_resp",
                "Change in Projected Angle (X);#theta^{Measured}_{X}; #theta^{True}_{X}", 
                NBINS, scat_bin_array, NBINS, scat_bin_array);
    scattering_proj_y_resp =
        new TH2D("scattering_proj_y_resp",
                "Change in Projected Angle (Y);#theta^{Measured}_{Y}; #theta^{True}_{Y}", 
                NBINS, scat_bin_array, NBINS, scat_bin_array);

    scattering_proj_x_R =
        new TH2D("scattering_proj_x_R",
                "Change in Projected Angle (X);#theta^{Measured}_{X}; #theta^{True}_{Y}", 
                NPBINS, pos_array, NBINS, scat_bin_array);
    scattering_proj_y_R =
        new TH2D("scattering_proj_y_R",
                "Change in Projected Angle (Y);#theta^{Measured}_{Y}; #theta^{True}_{X}", 
                NPBINS, pos_array, NBINS, scat_bin_array);
    theta_meas_x_graph =
        new TH1D("theta_meas_x_graph","Accepted Events;#theta^{True}_X",
                NBINS, scat_bin_array);
    theta_meas_y_graph =
        new TH1D("theta_meas_y_graph","Accepted Events;#theta^{True}_Y",
                NBINS, scat_bin_array);

    thetaX_all =
        new TH1D("thetaX_all","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaX_tof =
        new TH1D("thetaX_tof","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaX_proj =
        new TH1D("thetaX_proj","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaX_diff =
        new TH1D("thetaX_diff","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaX_chi =
        new TH1D("thetaX_chi","Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

    thetaY_all =
        new TH1D("thetaY_all","Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaY_tof =
        new TH1D("thetaY_tof","Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaY_proj =
        new TH1D("thetaY_proj","Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaY_diff =
        new TH1D("thetaY_diff","Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
    thetaY_chi =
        new TH1D("thetaY_chi","Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

    jUS=-1;
    jDS=-1;
    kUS=-1;
    kDS=-1;
    lUS=-1;
    lDS=-1;
    dte01;
    dte12;


    USrefplaneI = -1;
    DSrefplaneI = -1;
    USrefplaneZ = 0.;
    DSrefplaneZ = 0.;
    USabsPlaneI = -1;
    DSabsPlaneI = -1;
    centre = 0;
}

MCSAnalysis::~MCSAnalysis(){
    delete chain;
    delete refchain;
    delete mcchain;
    // delete tof10;
    // delete tof21;
    // delete tof10_sel;
    // delete tof21_sel;
    // delete calc_mom;
}

void MCSAnalysis::Write(){

    TText t1 = TText(0.18,0.8,"MICE Preliminary");
    TText t3 = TText(0.18,0.765,"ISIS cycle 2015/04");
    TText t2 = TText(0.18,0.73,"LiH, MAUS v3.3.2");
    t1.SetNDC(1);
    t1.SetTextSize(0.04);
    t1.SetTextFont(42);
    t2.SetNDC(1);
    t2.SetTextSize(0.03);
    t2.SetTextFont(42);
    t3.SetNDC(1);
    t3.SetTextSize(0.04);
    t3.SetTextFont(42);
    outfile->cd();
    theta_true_scat_graph->Write();
    theta_true_scat2_graph->Write();
    tof10->Write();
    tof21->Write();
    noweighttof10->Write();
    noweighttof21->Write();
    rawtime->Write();
    rawtime12->Write();
    tof10_sel->Write();
    tof21_sel->Write();
    calc_mom->Write();
    mctof10->Write();
    mctof21->Write();
    mctof10_sel->Write();
    mctof21_sel->Write();
    mccalc_mom->Write();
    mctrue_mom->Write();
    cuts_accept->Write();
    mccuts_accept->Write();
    trackno->Write();
    tofhitno->Write();
    diffradius->Write();
    projradius->Write();
    nocutprojradius->Write();
    cor_mom->Write();
    uncor_mom->Write();
    uncor_momsel->Write();
    thetaX_all->Write();
    thetaX_tof->Write();
    thetaX_diff->Write();
    thetaX_proj->Write();
    thetaX_chi->Write();
    thetaY_all->Write();
    thetaY_tof->Write();
    thetaY_diff->Write();
    thetaY_proj->Write();
    thetaY_chi->Write();
    residual->Write();
    residual12->Write();
    residual01->Write();
    residualUD->Write();
    refresidual->Write();
    refresidual12->Write();
    refresidual01->Write();
    refresidualUD->Write();
    TVectorD NoEntriestruth(1);
    NoEntriestruth[0] = theta_true_x_graph->GetEntries();
    NoEntriestruth.Write("graph");
    theta_true_x_graph->Write();
    theta_true_y_graph->Write();
    chi2pern->Write();
    mctheta_true_x_graph->Write();
    mctheta_true_y_graph->Write();
    mctheta_true_scat_graph->Write();
    MCTruth->Write();
    refMCTruth->Write();
    /* outfile->Close(); */
    
    TCanvas *c1 = new TCanvas();
    //calc_mom->GetXaxis()->SetRangeUser(120,280);
    //calc_mom->Draw();
    cor_mom->SetLineColor(2);
    cor_mom->Draw();
    meanmom = cor_mom->GetMean();
    semom = cor_mom->GetMeanError();
    rmsmom = cor_mom->GetRMS();
    errrmsmom = cor_mom->GetRMSError();
    empty_meanmom = mccalc_mom->GetMean();
    empty_semom = mccalc_mom->GetMeanError();
    empty_rmsmom = mccalc_mom->GetRMS();
    empty_errrmsmom = mccalc_mom->GetRMSError();
    cor_mom->Fit("gaus", "N");
    c1->SaveAs("calc_mom.pdf");
    c1->SaveAs("calc_mom.root");
    c1->Clear();
    calc_mom->GetXaxis()->SetRangeUser(120,280);
    cor_mom->Draw();
    c1->SaveAs("cor_mom.pdf");
    c1->Clear();
    MCTruth->GetXaxis()->SetTitle("P_{z} TOF MeV/c");
    MCTruth->Draw();
    c1->SaveAs("MCTruth.pdf");
    c1->Clear();
    refMCTruth->GetXaxis()->SetTitle("P_{z} TOF MeV/c");
    refMCTruth->Draw();
    c1->SaveAs("refMCTruth.pdf");
    c1->Clear();
    TOF01vsTOF12->GetXaxis()->SetTitle("P_{z} US (MeV/c)");
    TOF01vsTOF12->GetYaxis()->SetTitle("P_{z} DS (MeV/c)");
    TOF01vsTOF12->Draw("colz");
    float col1 = 0.2082*0.9;
    float col2 = 0.1664*0.9;
    float col3 = 0.5293*0.9;
    double fillcolour = TColor::GetColor(col1, col2, col3);
    c1->GetFrame()->SetFillColor(fillcolour);
    c1->SetFrameFillColor(fillcolour);
    TLine *line2 = new TLine(100,100,350,350);
    line2->SetLineColor(kRed);
    line2->Draw();
    c1->SaveAs("TOF01vsTOF12.pdf");
    c1->Clear();
    gStyle->SetOptStat(0);
    TOFvsMCTruth->GetXaxis()->SetTitle("Reconstructed #it{p} (MeV/c)");
    TOFvsMCTruth->GetXaxis()->SetLabelSize(0.05);
    TOFvsMCTruth->GetYaxis()->SetTitle("True #it{p} at abs. centre (MeV/c)");
    TOFvsMCTruth->GetYaxis()->SetLabelSize(0.05);
    TOFvsMCTruth->Draw("colz");
    fillcolour = TColor::GetColor(col1, col2, col3);
    c1->GetFrame()->SetFillColor(fillcolour);
    c1->SetFrameFillColor(fillcolour);
    TLine *line4 = new TLine(150,150,250,250);
    line4->SetLineColor(kRed);
    line4->Draw();
    t1.SetTextColor(0);
    t2.SetTextColor(0);
    t3.SetTextColor(0);
    t1.Draw();
    t3.Draw();
    t2.Draw();
    c1->SaveAs("TOFvsMCTruth.pdf");
    c1->Clear();
    t1.SetTextColor(1);
    t2.SetTextColor(1);
    t3.SetTextColor(1);
    refTOFvsMCTruth->GetXaxis()->SetTitle("Reconstructed P_{z} (MeV/c)");
    refTOFvsMCTruth->GetXaxis()->SetLabelSize(0.05);
    refTOFvsMCTruth->GetYaxis()->SetTitle("True P_{z} at abs. centre (MeV/c)");
    refTOFvsMCTruth->GetYaxis()->SetLabelSize(0.05);
    refTOFvsMCTruth->Draw("colz");
    c1->GetFrame()->SetFillColor(fillcolour);
    c1->SetFrameFillColor(fillcolour);
    line4->SetLineColor(kRed);
    line4->Draw();
    c1->SaveAs("refTOFvsMCTruth.pdf");

    c1->Clear();
    c1->SetFrameFillColor(0);
    energyloss->GetXaxis()->SetTitle("position (mm)");
    energyloss->GetYaxis()->SetTitle("MCTruth energy (MeV)");
    energyloss->Draw("colz");
    c1->SaveAs("energyloss.pdf");
    c1->Clear();
    energylossproj=energyloss->ProjectionX();
    energylossproj->SetTitle("position (mm)");
    energylossproj->Draw();
    c1->SaveAs("energylossproj.pdf");
    c1->Clear();
    energylosspro=energyloss->ProfileX();
    //energylosspro->SetTitle("position (mm)");
    energylosspro->Draw();
    c1->SaveAs("energylosspro.pdf");
    c1->Clear();
    gStyle->SetOptStat(0);
    TOF01vsMCTruth->SetTitle("TOF01 vs MCTruth");
    TOF01vsMCTruth->GetXaxis()->SetTitle("pz TOF01");
    TOF01vsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
    TOF01vsMCTruth->Draw("colz");
    fillcolour = TColor::GetColor(col1, col2, col3);
    c1->GetFrame()->SetFillColor(fillcolour);
    c1->SetFrameFillColor(fillcolour);
    t1.Draw();
    t3.Draw();
    t2.Draw();
    TLine *line5 = new TLine(100,100,350,350);
    line5->SetLineColor(kRed);
    line5->Draw();
    c1->SaveAs("TOF01vsMCTruth.pdf");
    c1->Clear();
    TOF12longPaul6thforvsMCTruth->GetXaxis()->SetTitle("pz TOF12 long Paul 6th for");
    TOF12longPaul6thforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
    TOF12longPaul6thforvsMCTruth->Draw("colz");
    TLine *line25 = new TLine(100,100,350,350);
    line25->SetLineColor(kRed);
    line25->Draw();
    c1->SaveAs("TOF12longPaul6thforvsMCTruth.pdf");
    c1->Clear();
    TOF01Paul6thforvsMCTruth->GetXaxis()->SetTitle("pz TOF01 Paul 6th for");
    TOF01Paul6thforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
    TOF01Paul6thforvsMCTruth->Draw("colz");
    TLine *line16 = new TLine(100,100,350,350);
    line16->SetLineColor(kRed);
    line16->Draw();
    c1->SaveAs("TOF01Paul6thforvsMCTruth.pdf");
    c1->Clear();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    //TOF01shortPaul6thforvsMCTruth->Fit("pol1","","",170,270);
    //TF1 *f3 = TOF01shortPaul6thforvsMCTruth->GetFunction("pol1");
    TOF01shortPaul6thforvsMCTruth->SetTitle("Corrected P upstream vs MC Truth");
    TOF01shortPaul6thforvsMCTruth->GetXaxis()->SetTitle("pz TOF01");
    TOF01shortPaul6thforvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
    TOF01shortPaul6thforvsMCTruth->SetFillColor(kViolet+1);

    //double fillcolour =GetColor(0.2082*0.9, 0.1664*0.9, 0.5293*0.9)
    c1->GetFrame()->SetFillColor(fillcolour);
    c1->SetFrameFillColor(fillcolour);
    TOF01shortPaul6thforvsMCTruth->Draw("colz");
    t1.Draw();
    t3.Draw();
    t2.Draw();
    TLine *line46 = new TLine(100,100,350,350);
    line46->SetLineColor(kRed);
    line46->Draw();
    c1->Update();
    c1->SaveAs("TOF01shortPaul6thforvsMCTruth.pdf");
    c1->Clear();
    c1->SetFrameFillColor(0);
    TOF01fordownvsMCTruth->GetXaxis()->SetTitle("pz John's formula TOF01");
    TOF01fordownvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
    TOF01fordownvsMCTruth->Draw("colz");
    TLine *line6 = new TLine(100,100,350,350);
    line6->SetLineColor(kRed);
    line6->Draw();
    c1->SaveAs("TOF01fordownvsMCTruth.pdf");
    c1->Clear();
    TOF12fordownvsMCTruth->GetXaxis()->SetTitle("pz John's formula TOF12");
    TOF12fordownvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
    TOF12fordownvsMCTruth->Draw("colz");
    TLine *line26 = new TLine(100,100,350,350);
    line26->SetLineColor(kRed);
    line26->Draw();
    c1->SaveAs("TOF12fordownvsMCTruth.pdf");
    c1->Clear();
    TOF12cobbvsMCTruth->GetXaxis()->SetTitle("pz cobb's formula TOF12");
    TOF12cobbvsMCTruth->GetYaxis()->SetTitle("pz MCTruth");
    TOF12cobbvsMCTruth->Draw("colz");
    TLine *line36 = new TLine(100,100,350,350);
    line36->SetLineColor(kRed);
    line36->Draw();
    c1->SaveAs("TOF12cobbvsMCTruth.pdf");
    c1->Clear();
    TOF01forupvsTOF01fordown->GetXaxis()->SetTitle("pz Ryan's formula downstream TOF12");
    TOF01forupvsTOF01fordown->GetYaxis()->SetTitle("pz John's formula TOF01");
    TOF01forupvsTOF01fordown->Draw("colz");
    TLine *line7 = new TLine(100,100,350,350);
    line7->SetLineColor(kRed);
    line7->Draw();
    c1->SaveAs("TOF01forupvsTOF01fordown.pdf");
    c1->Clear();
    gStyle->SetOptStat(0);
    residual->SetTitle("");
    residual->GetXaxis()->SetTitle("Reconstructed #it{p} #minus True #it{p} (MeV/c)");
    residual->GetXaxis()->SetLabelSize(0.05);
    residual->GetYaxis()->SetTitle("No. of events");
    residual->GetYaxis()->SetLabelSize(0.05);
    //residual->GetYaxis()->SetTitleOffset(1.6);
    residual->Draw();
    t1.Draw();
    t3.Draw();
    t2.Draw();
    c1->SaveAs("residual.pdf");
    c1->Clear();
    TText t11 = TText(0.63,0.785,"MICE Preliminary");
    TText t13 = TText(0.63,0.75,"ISIS cycle 2015/04");
    TText t12 = TText(0.63,0.715,"LiH, MAUS v3.3.2");
    c1->Clear();

    residual12->GetXaxis()->SetTitle("residual12");
    residual12->GetYaxis()->SetTitle("No. of events");
    residual12->Draw();
    /* t11.Draw(); */
    /* t13.Draw(); */
    /* t12.Draw(); */
    residual12->Fit("gaus", "N");
    c1->SaveAs("residual12.pdf");
    c1->Clear();
    refresidual12->GetXaxis()->SetTitle("P_{z} recon - MCTruth (MeV/c)");
    refresidual12->GetYaxis()->SetTitle("No. of events");
    refresidual12->Draw();
    t11.Draw();
    t13.Draw();
    t12.Draw();
    refresidual12->Fit("gaus", "N");
    c1->SaveAs("refresidual12.pdf");
    c1->Clear();
    t11.SetNDC(1);
    t11.SetTextSize(0.04);
    t11.SetTextFont(42);
    t12.SetNDC(1);
    t12.SetTextSize(0.03);
    t12.SetTextFont(42);
    t13.SetNDC(1);
    t13.SetTextSize(0.04);
    t13.SetTextFont(42);
    c1->Clear();
    gStyle->SetOptStat(1);
    residual01->GetXaxis()->SetTitle("residual01");
    residual01->GetYaxis()->SetTitle("No. of events");
    residual01->Draw();
    residual01->Fit("gaus", "N");
    c1->SaveAs("residual01.pdf");
    c1->Clear();
    gStyle->SetOptStat(1);
    refresidual01->GetXaxis()->SetTitle("refresidual01");
    refresidual01->GetYaxis()->SetTitle("No. of events");
    refresidual01->Draw();
    refresidual01->Fit("gaus", "N");
    c1->SaveAs("refresidual01.pdf");
    c1->Clear();
    gStyle->SetOptStat(1);
    refresidual->GetXaxis()->SetTitle("refresidual");
    refresidual->GetYaxis()->SetTitle("No. of events");
    refresidual->Draw();
    refresidual->Fit("gaus", "N");
    c1->SaveAs("refresidual.pdf");
    c1->Clear();
    residualUD->SetTitle("");
    residualUD->GetXaxis()->SetTitle("Momentum residual TOF01 vs TOF12 (MeV/c)");
    residualUD->GetYaxis()->SetTitle("No. of events");
    residualUD->GetYaxis()->SetTitleOffset(1.5);
    residualUD->Draw();
    t11.Draw();
    t13.Draw();
    t12.Draw();
    c1->SaveAs("residualUD.pdf");
    c1->Clear();
    refresidualUD->SetTitle("");
    refresidualUD->GetXaxis()->SetTitle("Momentum residual TOF01 vs TOF12");
    refresidualUD->GetYaxis()->SetTitle("No. of events");
    refresidualUD->GetYaxis()->SetTitleOffset(1.5);
    refresidualUD->Draw();
    t11.Draw();
    t13.Draw();
    t12.Draw();
    c1->SaveAs("refresidualUD.pdf");
    c1->Clear();
    TOF0Energy->GetXaxis()->SetTitle("TOF0Energy");
    TOF0Energy->GetYaxis()->SetTitle("No. of events");
    TOF0Energy->Draw();
    c1->SaveAs("TOF0Energy.pdf");
    c1->Clear();
    TOF1Energy->GetXaxis()->SetTitle("TOF1Energy");
    TOF1Energy->GetYaxis()->SetTitle("No. of events");
    TOF1Energy->Draw();
    c1->SaveAs("TOF1Energy.pdf");
    c1->Clear();
    TOF2Energy->GetXaxis()->SetTitle("TOF2Energy");
    TOF2Energy->GetYaxis()->SetTitle("No. of events");
    TOF2Energy->Draw();
    c1->SaveAs("TOF2Energy.pdf");
    c1->Clear();
    difEloss->GetXaxis()->SetTitle("difEloss");
    difEloss->GetYaxis()->SetTitle("No. of events");
    difEloss->Draw();
    c1->SaveAs("difEloss.pdf");
    c1->Clear();
    gStyle->SetOptStat(1);
    pathlengthabs->GetXaxis()->SetTitle("Path Length in Absorber (mm)");
    pathlengthabs->GetYaxis()->SetTitle("No. of events");
    pathlengthabs->Draw();
    c1->SaveAs("pathlengthabs.pdf");
    c1->SaveAs("pathlengthabs.root");
    c1->Clear();
    gStyle->SetOptStat(1);
    gStyle->SetOptStat("nemR");
    theta_true_x_graph->GetXaxis()->SetTitle("#it{#theta_{x}} (mrad)");
    theta_true_x_graph->GetYaxis()->SetTitle("No. of events");
    //theta_true_x_graph->Sumw2();
    theta_true_x_graph->SetFillColor(kOrange-2);
    theta_true_x_graph->DrawNormalized();

    //TH1D* thetaX_reco = (TH1D*)outfile->Get("thetaX_recoGEANT");
    //thetaX_reco->SetMarkerStyle(20);
    //thetaX_reco->Sumw2();
    //thetaX_reco->DrawNormalized("e1 p SAMES");
    TLegend * legend20 = new TLegend(0.1,0.7,0.34,0.9);
    legend20->AddEntry(theta_true_x_graph,"MC","f");
    //legend20->AddEntry(thetaX_reco,"Reco sim","p");
    legend20->Draw();
    gPad->Update();
    //TPaveStats *st = (TPaveStats*)thetaX_reco->FindObject("stats");
    //st->SetY1NDC(0.4); 
    //st->SetY2NDC(0.6); 
    //double chi2 = theta_true_x_graph->Chi2Test(thetaX_reco,"WUP");
    //TText t10 = TText(0.12,0.585,Form("Chi^2 = %f",chi2));
    /* t10.SetNDC(1); */
    /* t10.SetTextSize(0.03); */
    /* t10.SetTextFont(42); */
    /* t10.Draw(); */
    //double Kol = theta_true_x_graph->KolmogorovTest(thetaX_reco);
    //std::cout << "Kolmogorov " << Kol << std::endl;
    //TText t20 = TText(0.12,0.55,Form("Kolmogorov = %f",Kol));
    /* t20.SetNDC(1); */
    /* t20.SetTextSize(0.03); */
    /* t20.SetTextFont(42); */
    /* t20.Draw(); */
    TText t100 = TText(0.652,0.585,"MICE Preliminary");
    TText t300 = TText(0.652,0.55,"ISIS cycle 2015/04");
    TText t200 = TText(0.652,0.515,"LiH, MAUS v3.3.2");
    t100.SetNDC(1);
    t100.SetTextSize(0.04);
    t100.SetTextFont(42);
    t200.SetNDC(1);
    t200.SetTextSize(0.03);
    t200.SetTextFont(42);
    t300.SetNDC(1);
    t300.SetTextSize(0.04);
    t300.SetTextFont(42);
    t100.Draw();
    t300.Draw();
    t200.Draw();
    c1->SaveAs("theta_true_x_graph.pdf");
    c1->SaveAs("theta_true_x_graph.root");
    c1->Clear();
    gStyle->SetOptStat(1);
    gStyle->SetOptStat("nemR");
    theta_true_y_graph->GetXaxis()->SetTitle("#it{#theta_y} (mrad)");
    theta_true_y_graph->GetYaxis()->SetTitle("No. of events");
    theta_true_y_graph->DrawNormalized();

    theta_true_y_graph->SetFillColor(kOrange-2);
    //theta_true_y_graph->Sumw2();
    //TH1D* thetaY_reco = (TH1D*)outfile->Get("thetaY_recoGEANT");
    //thetaY_reco->SetMarkerStyle(20);
    //thetaY_reco->Sumw2();
    //thetaY_reco->DrawNormalized("e1 p SAMES");
    TLegend * legend21 = new TLegend(0.1,0.7,0.34,0.9);
    legend21->AddEntry(theta_true_y_graph,"MC","f");
    //legend21->AddEntry(thetaY_reco,"Reco sim","p");
    legend21->Draw();
    gPad->Update();
    //TPaveStats *st1 = (TPaveStats*)thetaY_reco->FindObject("stats");
    //st1->SetY1NDC(0.4); 
    //st1->SetY2NDC(0.6); 
    //chi2 = theta_true_y_graph->Chi2Test(thetaY_reco,"WUP");
    //TText t101 = TText(0.12,0.585,Form("Chi^2 = %.17g",chi2));
    /* t101.SetNDC(1); */
    /* t101.SetTextSize(0.03); */
    /* t101.SetTextFont(42); */
    /* t101.Draw(); */
    //Kol = theta_true_y_graph->KolmogorovTest(thetaY_reco);
    //std::cout << "Kolmogorov " << Kol << std::endl;
    //TText t21 = TText(0.12,0.55,Form("Kolmogorov = %.17g",Kol));
    /* t21.SetNDC(1); */
    /* t21.SetTextSize(0.03); */
    /* t21.SetTextFont(42); */
    /* t21.Draw(); */
    c1->SaveAs("theta_true_y_graph.pdf");
    c1->SaveAs("theta_true_y_graph.root");
    c1->Clear();
    eventnovspz->Draw("colz");
    c1->SaveAs("eventnovspz.pdf");
    c1->SaveAs("eventnovspz.root");


    delete c1;
}

void MCSAnalysis::print(){
    std::cout << "mean momentum " << meanmom << std::endl;
    std::cout << "standard error on mean momentum " << semom << std::endl;
    std::cout << "rms of momentum " << rmsmom << std::endl;
    std::cout << "Error on rms of momentum " << errrmsmom << std::endl;
    std::cout << "empty mean momentum " << empty_meanmom << std::endl;
    std::cout << "empty standard error on mean momentum " << empty_semom << std::endl;
    std::cout << "empty rms of momentum " << empty_rmsmom << std::endl;
    std::cout << "empty Error on rms of momentum " << empty_errrmsmom << std::endl;
    bwmom = (1.357e3/pow(TOF_upper_limit - 24.28,2)) - (1.357e3/pow(TOF_lower_limit - 24.28,2));
    std::cout << "width of bin mom " << bwmom << std::endl;
    std::cout << "bin width sqrt12 " << bwmom/sqrt(12) << std::endl;
    std::cout << "tX_rG_int " << tX_rG_int << std::endl;
    std::cout << "% events reaching TOF2 " << percentTOF2 << std::endl;
    std::cout << "% events reaching TOF2 " << (TOF2Hit/TOF1Hit)*100 << std::endl;
    std::cout << "No. of selected events US & " << int(_USset.N()) << " & " << int(DSentries) << " \\" << std::endl;
    std::cout << "No. of selected events DS " << int(_DSset.N()) << std::endl;
    std::cout << "No. of selected events US & " << int(_USMCset.N()) << " & " << int(DSrefentries) << " \\" << std::endl;
    std::cout << "No. of selected events DS " << int(_DSMCset.N()) << std::endl;
}

void MCSAnalysis::Execute(int mode){
    
        ReadOffset(0);
    if(offset_needed==0){
        ReadOffset(1);
    }
    /* UScor = 0; */
    /* refUScor = 0; */
    /* DScor = 0; */
    /* refDScor = 0; */

    int verbose_level = 1;
    // MAUS::GlobalsManager::InitialiseGlobals(JsonWrapper::JsonToString(SetupConfig(verbose_level)));
    /* toffitting(); */

    std::cout << "mode " << mode << std::endl;

    if (mode == 0){
    dataevents(mode);
    dataSelection(mode);
        //generateMCSResponse();
        // DoUnfolding();
        DoDeconvolution(modelname1.c_str(), mode, 1);
        TruthData(mode);
    }
    else if (mode == 1){
    dataevents(mode);
    dataSelection(mode);
        referenceSelection();
        DoUnfolding(_USset, _DSset, _USMCset, _DSMCset);
        DoFFTDeconvolution();
        ConvolveWithInputDistribution(modelname1.c_str());
        DoDeconvolution(modelname1.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname2.c_str());
        DoDeconvolution(modelname2.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname3.c_str());
        DoDeconvolution(modelname3.c_str(), mode, 1);
        //FitGaussian(outfilename.c_str());
        //CalculateChi2(outfilename.c_str(), modelname2.c_str());
        //Write();
    }
    else if (mode == -1){

    dataevents(mode);
    dataSelection(mode);
        referenceSelection();
        DoUnfolding(_USset, _DSset, _USMCset, _DSMCset);
        DoFFTDeconvolution();
        ConvolveWithInputDistribution(modelname1.c_str());
        DoDeconvolution(modelname1.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname2.c_str());
        DoDeconvolution(modelname2.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname3.c_str());
        DoDeconvolution(modelname3.c_str(), mode, 1);
        TruthData(mode);
        //FitGaussian(outfilename.c_str());
        //CalculateChi2(outfilename.c_str(), modelname2.c_str());
        complots(mode);
        Write();
        print();
    }
    else if (mode == -7){
    dataevents(mode);
    dataSelection(mode);
        referenceSelection();
        DoUnfolding(_USset, _DSset, _USMCset, _DSMCset);
        DoFFTDeconvolution();
        ConvolveWithInputDistribution(modelname1.c_str());
        DoDeconvolution(modelname1.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname2.c_str());
        DoDeconvolution(modelname2.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname3.c_str());
        DoDeconvolution(modelname3.c_str(), mode, 1);
        TruthData(mode);
        //FitGaussian(outfilename.c_str());
        //CalculateChi2(outfilename.c_str(), modelname2.c_str());
        Write();
        print();
    }
    else if (mode == -2){
    dataevents(mode);
    dataSelection(mode);
        referenceSelection();
        DoUnfolding(_USset, _DSset, _USMCset, _DSMCset);
        DoFFTDeconvolution();
        ConvolveWithInputDistribution(modelname1.c_str());
        DoDeconvolution(modelname1.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname2.c_str());
        DoDeconvolution(modelname2.c_str(), mode, 1);
        ConvolveWithInputDistribution(modelname3.c_str());
        DoDeconvolution(modelname3.c_str(), mode, 1);
        TruthData(mode);
        //FitGaussian(outfilename.c_str());
        //CalculateChi2(outfilename.c_str(), modelname2.c_str());
        complots(mode);
        Write();
        print();
    }
    else if (mode == 2){
    dataevents(mode);
    dataSelection(mode);
        TruthData(mode);
    } 
    else if (mode == 10){
    dataevents(mode);
    dataSelection(mode);
        TruthData(mode);
    } 
    else if (mode == 5){
    dataevents(mode);
    dataSelection(mode);
        referenceSelection();
        DoUnfolding(_USset, _DSset, _USMCset, _DSMCset);
        DoFFTDeconvolution();
        TruthData(mode);
        ConvolveWithVirtualInputDistribution(modelname1.c_str());
        DoDeconvolution(modelname1.c_str(), mode, 1);
        ConvolveWithVirtualInputDistribution(modelname2.c_str());
        DoDeconvolution(modelname2.c_str(), mode, 1);
        ConvolveWithVirtualInputDistribution(modelname3.c_str());
        DoDeconvolution(modelname3.c_str(), mode, 1);
    }
    else if (mode == 4){
    dataevents(mode);
    dataSelection(mode);
        TruthData(mode);
        referenceSelection();
        DoFFTDeconvolution();
        TruthData(mode-1);
        DoUnfolding(USTruthSetLiH, DSTruthSetLiH, USTruthSet, DSTruthSet);
        ConvolveWithVirtualInputDistribution(modelname1.c_str());
        DoVirtualDeconvolution(modelname1.c_str(), 1);
        ConvolveWithVirtualInputDistribution(modelname2.c_str());
        DoVirtualDeconvolution(modelname2.c_str(), 1);
        ConvolveWithVirtualInputDistribution(modelname3.c_str());
        DoVirtualDeconvolution(modelname3.c_str(), 1);
    }
    else if (mode == 6){
        toffitting();
        Write();
    }
    else {
        std::cout<<"Unknown operation mode"<<std::endl;
    } 
    if(offset_needed!=0){
        WriteOffset(isMC);
    }
}

void MCSAnalysis::toffitting(){
    chain->SetBranchAddress("RunNumber", &runnumber);
    chain->SetBranchAddress("TOFBranch", &tofevent);
    chain->SetBranchAddress("SciFiBranch", &scifievent);
    chain->SetBranchAddress("CkovBranch", &ckovevent);
    chain->SetBranchAddress("KLBranch", &klevent);
    chain->SetBranchAddress("EMRBranch", &emrevent);
    chain->SetBranchAddress("MCEvent", &mcevent);

    int Nentries = chain->GetEntries();
    std::cout << Nentries << std::endl;
    for (int i=0; i<Nentries; i++){
        chain->GetEntry(i);
        if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 
        MomentumFromTOF(true);
        TimeofFlight(true);
        TimeofFlight12(true);
    }
    TF1* fittof01 = new TF1("fit01","gaus",24.9,25.6);
    if (isMC==1) fittof01 = new TF1("fit01","gaus",25,26.2);
    rawtime->Fit("fit01", "R");
    dte01 = fittof01->GetParameter(1);
    TF1* fittof12 = new TF1("fit12","gaus",27,28.2);
    if (isMC==1) fittof12 = new TF1("fit12","gaus",27,28.2);
    rawtime12->Fit("fit12", "R");
    dte12 = fittof12->GetParameter(1);
    float depeak = dte01 - dte01data;
    /* dte01 -= dte01data; */
    if (isMC==1) TOF_lower_limit = TOF_lower_limit*dte01/dte01data; 
    if (isMC==1) TOF_upper_limit = TOF_upper_limit*dte01/dte01data; 
}

void MCSAnalysis::dataevents(int mode){
    if (isMC==1) dte01 = dte01MC; 
    else dte01 = dte01data; 
    if (isMC==1) dte12 = dte12MC; 
    else dte12 = dte12data; 
    if (isMC==1) TOF_lower_limit = TOF_lower_limit*dte01MC/dte01data; 
    if (isMC==1) TOF_upper_limit = TOF_upper_limit*dte01MC/dte01data; 
    if (isMC==1) TOF_lower_limit_ref = TOF_lower_limit_ref*dte01MC/dte01data; 
    if (isMC==1) TOF_upper_limit_ref = TOF_upper_limit_ref*dte01MC/dte01data; 
    if (isMC==1) TOF_lower_limit = TOF_lower_limit+0.05; 
    if (isMC==1) TOF_upper_limit = TOF_upper_limit+0.05; 
    if (isMC==1) TOF_lower_limit_ref = TOF_lower_limit_ref+0.05; 
    if (isMC==1) TOF_upper_limit_ref = TOF_upper_limit_ref+0.05; 
    chain->SetBranchAddress("RunNumber", &runnumber);
    chain->SetBranchAddress("TOFBranch", &tofevent);
    chain->SetBranchAddress("SciFiBranch", &scifievent);
    chain->SetBranchAddress("CkovBranch", &ckovevent);
    chain->SetBranchAddress("KLBranch", &klevent);
    chain->SetBranchAddress("EMRBranch", &emrevent);
    chain->SetBranchAddress("MCEvent", &mcevent);

    int Nentries = chain->GetEntries();
    Collection comUSeventset, comDSeventset;
    Collection outsideUSeventset, outsideDSeventset;
    Collection insideUSeventset, insideDSeventset;
    double pz = 0.;
    double ptruth = 0.;
    Vars projdif;
    Vars proj;
    std::vector<float> chi2n;
    std::cout << Nentries << std::endl;
    int Nevents = 0;
    int a172eventcount = 0;
    int a200eventcount = 0;
    int a240eventcount = 0;

    for (int i=0; i<Nentries; i++){
        /* int beamtype = 0; */
    /* std::cout << i << std::endl; */
        int pid = -13;
        bool difcut = false;
        bool chicut = false;
        bool TOFcut = false;
        bool fidcut = false;
        bool a172event = false;
        bool a200event = false;
        bool a240event = false;
        isEmpty=false;
        chain->GetEntry(i);
        if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 
        /* MomentumFromTOF(true); */
        /* TimeofFlight(true); */
        /* TimeofFlight12(true); */

        if (isMC==1) {

            findVirtualPlanes();
        }


        if (mode==-2) { 
            if (mcevent->GetVirtualHits()->size() >= 49) {
                if (mcevent->GetVirtualHits()->at(48).GetParticleId()!=-13) continue;
            }
        }



        FillCollectionSciFi(comUSeventset, jUS, kUS, lUS, pid, pz, tracker0, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        FillCollectionSciFi(comDSeventset, jDS, kDS, lDS, pid, pz, tracker1, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        if ( !MatchUSDS() ) {
            continue;
            /* if (jUS == -1 || kUS == -1){ */
            /*     continue; */
            /* } */ 
        }


        bool trackmatched = MatchUSDS();

        double diff = 0;
        pz = 0;
        projdif = RadialSelection(pz,_sys["diffpos"],difradius);
        if ( sqrt(projdif.X*projdif.X+projdif.Y*projdif.Y)<90 ) difcut = true;
        proj = RadialSelection(pz,_sys["DStrkr5"],meanp);
        if ( sqrt(proj.X*proj.X+proj.Y*proj.Y)<meanp ) fidcut = true;

        Nevents++;
        chi2n = FillCrossCheck();
        for (int kk=0;kk<chi2n.size();kk++) if ( chi2n[kk]<4 ) chicut = true;
        if ( chi2n[0]<4 ) chicut = true;
        //pz = CorMomFromTOF(USeventset.E(i), 0, diff, true);
        if ( PIDSelection(true) ) TOFcut = true ;

        for (int kk=0;kk<26;kk++) {
            if (runnumber==list172[kk]) a172event=true;
        }
        for (int kk=0;kk<28;kk++) {
            if (runnumber==list200[kk]) a200event=true;
        }
        for (int kk=0;kk<35;kk++) {
            if (runnumber==list240[kk]) a240event=true;
        }

        TOFtotal += 1;
        if (a172event) {
            beamtype = 172;
            a172eventcount += 1 ;
        }
        if (a200event) {
            beamtype = 200;
            a200eventcount += 1;
        }
        if (a240event) {
            beamtype = 240;
            a240eventcount += 1;
        }

        /* if (isMC==1) { */
        /* std::cout << "beamtype " << beamtype << std::endl; */
        /* std::cout << "runnumber " << runnumber << std::endl; */
        /* } */

        Vars mceventUSTruthAbsHit, mceventDSTruthAbsHit;
        centre = USabsPlaneI + 1;
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            FillVarsVirtual(mceventUSTruthAbsHit, USabsPlaneI, centre);
            FillVarsVirtual(mceventDSTruthAbsHit, DSabsPlaneI, centre);
            eventMCUSAllTOF.append_instance(mceventUSTruthAbsHit);
            eventMCDSAllTOF.append_instance(mceventDSTruthAbsHit);
        }

        ptruth = mceventUSTruthAbsHit.pz;
        FillCollectionSciFi(USeventset, jUS, kUS, lUS, pid, pz, tracker0, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        FillCollectionSciFi(DSeventset, jDS, kDS, lDS, pid, pz, tracker1, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        jUS = -1;
        kUS = -1;
        lUS = -1;
        jDS = -1;
        kDS = -1;
        lDS = -1;
    }
    data_make_beam_histograms(comUSeventset, "Upstream, Data", "datacomUSeventset", comUSeventset);
    data_make_beam_histograms(comDSeventset, "Downstream, Data", "datacomDSeventset", comDSeventset);
    data_make_beam_histograms(USeventset, "Upstream, Data", "dataUSeventset", USeventset);
    data_make_beam_histograms(DSeventset, "Downstream, Data", "dataDSeventset", DSeventset);
    std::cout << "a172eventcount " << a172eventcount << std::endl;
    std::cout << "a200eventcount " << a200eventcount << std::endl;
    std::cout << "a240eventcount " << a240eventcount << std::endl;
    std::cout << "TOFtotal " << TOFtotal << std::endl;
    weight172 = TOFtotal/a172eventcount;
    weight200 = TOFtotal/a200eventcount;
    weight240 = TOFtotal/a240eventcount;
    if (isMC==1) {
        weight172 = 1;
        weight200 = 1;
        weight240 = 1;
    }
}

void MCSAnalysis::dataSelection(int mode){

    int Nentries = USeventset.N();
    Collection USAllTOF, DSAllTOF, USDifSel, DSDifSel, USFidSel, DSFidSel, USChiSel, DSChiSel;
    Collection MCUSAllTOF, MCDSAllTOF,  MCUSDifSel, MCDSDifSel,  MCUSFidSel, MCDSFidSel, MCUSChiSel, MCDSChiSel;
    double pz = 0.;
    int Nevents = 0;
    Vars projdif;
    Vars proj;
    Vars USholdvar;
    Vars DSholdvar;
    std::vector<float> chi2n;
    std::cout << Nentries << std::endl;

    for (int i=0; i<Nentries; i++){
        isEmpty=false;

        if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 
        if (isMC==1) {
            findVirtualPlanes();
        }

        /* if (mode==-7) { */
        /*     if (USeventset.E(i).jDS == -1 || USeventset.E(i).kDS == -1){ */
        /*         continue; */
        /*     } */
        /* } */

        cuts_accept->Fill("US Track Found",1);    
        USholdvar = USeventset.E(i);
        DSholdvar = DSeventset.E(i);
        FillEventCollection(USAllTOF, USholdvar);
        FillEventCollection(DSAllTOF, DSholdvar);
        Vars mcUSTruthAbsHit, mcDSTruthAbsHit;
        centre = USabsPlaneI + 1;
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            FillVarsVirtual(mcUSTruthAbsHit, USabsPlaneI, centre);
            FillVarsVirtual(mcDSTruthAbsHit, DSabsPlaneI, centre);
            MCUSAllTOF.append_instance(USeventset.E(i));
            MCDSAllTOF.append_instance(DSeventset.E(i));
        }

        cuts_accept->Fill("All Events",1);
        if ( USeventset.E(i).difcut==0 ) continue;
        cuts_accept->Fill("Diffuser Cut",1);

        FillEventCollection(USDifSel, USholdvar);
        FillEventCollection(DSDifSel, DSholdvar);
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            MCUSDifSel.append_instance(USeventset.E(i));
            MCDSDifSel.append_instance(DSeventset.E(i));
        }
        if ( USeventset.E(i).fidcut==0 ) continue;
        FillEventCollection(USFidSel, USholdvar);
        FillEventCollection(DSFidSel, DSholdvar);
        cuts_accept->Fill("Fiducial Selection",1);
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            MCUSFidSel.append_instance(USeventset.E(i));
            MCDSFidSel.append_instance(DSeventset.E(i));
        }

        if ( USeventset.E(i).chicut==0 ) continue;
        FillEventCollection(USChiSel, USholdvar);
        FillEventCollection(DSChiSel, DSholdvar);
        cuts_accept->Fill("chi2ndf",1);
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            MCUSChiSel.append_instance(USeventset.E(i));
            MCDSChiSel.append_instance(DSeventset.E(i));
        }


        if ( USeventset.E(i).TOFcut==0 ) continue;
        cuts_accept->Fill("TOF Selection",1);
        double diff = 0;
        USholdvar = USeventset.E(i);
        pz = CorMomFromTOF(USholdvar, diff, true);

        FillEventCollection(_USset, USholdvar);
        FillEventCollection(_DSset, DSholdvar);
        FillEventCollection(_UStmpset, USholdvar);
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            mcreconUSTruthSet.append_instance(USeventset.E(i));
            mcreconDSTruthSet.append_instance(DSeventset.E(i));
        }

        TOF1Hit += 1;
        if (USeventset.E(i).TOF12 > 0 && USeventset.E(i).TOF12 < 100){
           TOF2Hit += 1;  
        }
        Nevents++;
        cor_mom->Fill(pz);
        tof10_sel->Fill(USeventset.E(i).TOF01);
        tof21_sel->Fill(USeventset.E(i).TOF12);
    }
    TruthGraph(mcreconUSTruthSet, mcreconDSTruthSet, true);
    data_make_beam_histograms(USAllTOF, "Upstream, Data", "dataUS_alltof", MCUSAllTOF);
    data_make_beam_histograms(DSAllTOF, "Downstream, Data", "dataDS_alltof", MCDSAllTOF);
    data_make_beam_histograms(USDifSel, "Upstream, Data", "dataUS_dif", MCUSDifSel);
    data_make_beam_histograms(DSDifSel, "Downstream, Data", "dataDS_dif", MCDSDifSel);
    data_make_beam_histograms(USFidSel, "Upstream, Data", "dataUS_fid", MCUSFidSel);
    data_make_beam_histograms(DSFidSel, "Downstream, Data", "dataDS_fid", MCDSFidSel);
    data_make_beam_histograms(USChiSel, "Upstream, Data", "dataUS_chi", MCUSChiSel);
    data_make_beam_histograms(DSChiSel, "Downstream, Data", "dataDS_chi", MCDSChiSel);
    data_make_beam_histograms(_USset, "Upstream, Data", "dataUS", mcreconUSTruthSet);
    data_make_beam_histograms(_DSset, "Downstream, Data", "dataDS", mcreconDSTruthSet);
    make_acceptance_histograms(_USset, _DSset, "", "dataAcc");
    PlotRunInfo();
    make_beam_histograms(MCUSAllTOF, "Upstream, Data", "MCUS_alltof");
    make_beam_histograms(MCDSAllTOF, "Downstream, Data", "MCDS_alltof");
    make_beam_histograms(MCUSDifSel, "Upstream, Data", "MCUS_dif");
    make_beam_histograms(MCDSDifSel, "Downstream, Data", "MCDS_dif");
    make_beam_histograms(MCUSFidSel, "Upstream, Data", "MCUS_fid");
    make_beam_histograms(MCDSFidSel, "Downstream, Data", "MCDS_fid");
    make_beam_histograms(MCUSChiSel, "Upstream, Data", "MCUS_chi");
    make_beam_histograms(MCDSChiSel, "Downstream, Data", "MCDS_chi");
    make_beam_histograms(mcreconUSTruthSet, "Upstream, Data", "MCUS");
    make_beam_histograms(mcreconDSTruthSet, "Downstream, Data", "MCDS");
}

void MCSAnalysis::complots(int mode){

    int Nentries = USeventset.N();
    double pz = 0.;
    Collection outsideUSeventset, outsideDSeventset, insideUSeventset, insideDSeventset;
    int Nevents = 0;
    Vars USholdvar;
    Vars DSholdvar;
    for (int i=0; i<Nentries; i++){
        nocutprojradius->Fill(sqrt(USeventset.E(i).projX*USeventset.E(i).projX + USeventset.E(i).projY*USeventset.E(i).projY));
        std::vector<double> projTheta = RotDefineProjectionAngles(USeventset.E(i), DSeventset.E(i), angdef);
        thetaX_all->Fill(projTheta[1]);
        thetaY_all->Fill(projTheta[0]);
        isEmpty=false;
        if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 


        if ( USeventset.E(i).difcut && USeventset.E(i).fidcut && USeventset.E(i).chicut) {
            noweighttof10->Fill(USeventset.E(i).TOF01);
            noweighttof21->Fill(USeventset.E(i).TOF12);
            if (USeventset.E(i).beamtype==172) {
                /* std::cout << "beamtype " << beamtype << std::endl; */
                /* std::cout << "USeventset.E(i).TOF01 " << USeventset.E(i).TOF01 << std::endl; */
                tof10->Fill(USeventset.E(i).TOF01,weight172);
                tof21->Fill(USeventset.E(i).TOF12,weight172);
            }
            if (USeventset.E(i).beamtype==200) {
                tof10->Fill(USeventset.E(i).TOF01,weight200);
                tof21->Fill(USeventset.E(i).TOF12,weight200);
            }
            if (USeventset.E(i).beamtype==240) {
                tof10->Fill(USeventset.E(i).TOF01,weight240);
                tof21->Fill(USeventset.E(i).TOF12,weight240);
            }
                /* compz->Fill(USeventset.E(i).pz); */
            std::vector<double> projTheta = RotDefineProjectionAngles(USeventset.E(i), DSeventset.E(i), angdef);
            thetaX_tof->Fill(projTheta[1]);
            thetaY_tof->Fill(projTheta[0]);
        }

        if ( USeventset.E(i).difcut && USeventset.E(i).chicut && USeventset.E(i).TOFcut) {
            projradius->Fill(sqrt(USeventset.E(i).projX*USeventset.E(i).projX + USeventset.E(i).projY*USeventset.E(i).projY));
            std::vector<double> projTheta = RotDefineProjectionAngles(USeventset.E(i), DSeventset.E(i), angdef);
            thetaX_proj->Fill(projTheta[1]);
            thetaY_proj->Fill(projTheta[0]);
        }

        if ( USeventset.E(i).difcut && USeventset.E(i).TOFcut && USeventset.E(i).fidcut) {
            chi2pern->Fill(USeventset.E(i).chi2nUS);
            std::vector<double> projTheta = RotDefineProjectionAngles(USeventset.E(i), DSeventset.E(i), angdef);
            thetaX_chi->Fill(projTheta[1]);
            thetaY_chi->Fill(projTheta[0]);
        }

        if ( USeventset.E(i).TOFcut && USeventset.E(i).fidcut && USeventset.E(i).chicut) {
            diffradius->Fill(sqrt(USeventset.E(i).projdifX*USeventset.E(i).projdifX + USeventset.E(i).projdifY*USeventset.E(i).projdifY));
            std::vector<double> projTheta = RotDefineProjectionAngles(USeventset.E(i), DSeventset.E(i), angdef);
            thetaX_diff->Fill(projTheta[1]);
            thetaY_diff->Fill(projTheta[0]);

        }
        if ( USeventset.E(i).difcut && USeventset.E(i).fidcut && USeventset.E(i).chicut && USeventset.E(i).TOFcut) {
            std::vector<double> projTheta = RotDefineProjectionAngles(USeventset.E(i), DSeventset.E(i), angdef);
            if ((projTheta[1]>0.02 || projTheta[1]<-0.02) && (projTheta[0]>0.02 || projTheta[0]<-0.02)){
                USholdvar = USeventset.E(i);
                DSholdvar = DSeventset.E(i);
                FillEventCollection(outsideUSeventset, USholdvar);
                FillEventCollection(outsideDSeventset, DSholdvar);
            }
            if (projTheta[1]<0.02 && projTheta[1]>-0.02 && projTheta[0]<0.02 && projTheta[0]>-0.02){
                USholdvar = USeventset.E(i);
                DSholdvar = DSeventset.E(i);
                FillEventCollection(insideUSeventset, USholdvar);
                FillEventCollection(insideDSeventset, DSholdvar);
            }
        }
    }
    data_make_beam_histograms(outsideUSeventset, "Upstream, Data", "dataoutsideUSeventset", outsideUSeventset);
    data_make_beam_histograms(outsideDSeventset, "Downstream, Data", "dataoutsideDSeventset", outsideDSeventset);
    data_make_beam_histograms(insideUSeventset, "Upstream, Data", "datainsideUSeventset", insideUSeventset);
    data_make_beam_histograms(insideDSeventset, "Downstream, Data", "datainsideDSeventset", insideDSeventset);
}

void MCSAnalysis::referenceSelection(){
    // Set addresses for tree selection
    refchain->SetBranchAddress("TOFBranch", &tofevent);
    refchain->SetBranchAddress("SciFiBranch", &scifievent);
    refchain->SetBranchAddress("CkovBranch", &ckovevent);
    refchain->SetBranchAddress("KLBranch", &klevent);
    refchain->SetBranchAddress("EMRBranch", &emrevent);
    refchain->SetBranchAddress("MCEvent", &mcevent);

    // Restrict the number of entries to less than or equal to the refchain entries

    int Nentries = refchain->GetEntries();
    Collection USAllTOFEmpty, DSAllTOFEmpty, USPreRadSelEmpty, DSPreRadSelEmpty;
    Collection MCUSAllTOFEmpty, MCDSAllTOFEmpty, MCUSPreRadSelEmpty, MCDSPreRadSelEmpty, MCUSEmpty, MCDSEmpty;
    double pz = 0.;
    double ptruth = 0.;
    Vars projdif;
    Vars proj;
    Vars USholdvar;
    Vars DSholdvar;
    std::vector<float> chi2n;
    // Loop over all tree entries.
    std::cout << Nentries << std::endl;
    std::cout<<"reference"<<std::endl; 
    for (int i=0; i<Nentries; i++){
        int pid = -13;
        bool difcut = false;
        bool chicut = false;
        bool TOFcut = false;
        bool fidcut = false;

        isEmpty=true;

        refchain->GetEntry(i);
        if (i%100000==0) std::cout<<"Event "<<i<<"\n"; 
        if (isMC==1) {
            findVirtualPlanes();
        }
        // Set cuts based on the TOFs, ckov, kl, and EMR information

        // Locate the tracker reference planes. To be agnostic locate
        // the downstream most station of the upstream tracker and the
        //// upstream most station of the downstream tracker.
        /* if (mode==-2) { */ 
        /* 	if (mcevent->GetVirtualHits()->size() >= 49) { */
        /* 		if (mcevent->GetVirtualHits()->at(48).GetParticleId()!=-13) continue; */
        /* 	} */
        /* } */
        if (mode==-2) {  
            if (mcevent->GetVirtualHits()->size() >= 49) { 
                if (pid = mcevent->GetVirtualHits()->at(48).GetParticleId()); 
            }
        }

        mccuts_accept->Fill("All Events",1);
        if ( !MatchUSDS() ) {
            continue;
            /* if (jUS == -1 || kUS == -1){ */
            /* 	continue; */
            /* } */ 
        }
        mccuts_accept->Fill("US Track Found",1);
        FillCollectionSciFi(USAllTOFEmpty, jUS, kUS, lUS, pid, pz, tracker0, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        FillCollectionSciFi(DSAllTOFEmpty, jDS, kDS, lDS, pid, pz, tracker1, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        Vars mcUSTruthAbsHitEmpty, mcDSTruthAbsHitEmpty;
        int centre = USabsPlaneI + 1;
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            FillVarsVirtual(mcUSTruthAbsHitEmpty, USabsPlaneI, centre);
            FillVarsVirtual(mcDSTruthAbsHitEmpty, DSabsPlaneI, centre);
            MCUSAllTOFEmpty.append_instance(mcUSTruthAbsHitEmpty);
            MCDSAllTOFEmpty.append_instance(mcDSTruthAbsHitEmpty);
        }
        ptruth = mcUSTruthAbsHitEmpty.pz;

        projdif = RadialSelection(pz,_sys["diffpos"],difradius);
        if ( sqrt(projdif.X*projdif.X+projdif.Y*projdif.Y)>90 ) continue;
        mccuts_accept->Fill("Diffuser Cut",1);

        FillCollectionSciFi(USPreRadSelEmpty, jUS, kUS, lUS, pid, pz, tracker0, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        FillCollectionSciFi(DSPreRadSelEmpty, jDS, kDS, lDS, pid, pz, tracker1, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            MCUSPreRadSelEmpty.append_instance(mcUSTruthAbsHitEmpty);
            MCDSPreRadSelEmpty.append_instance(mcDSTruthAbsHitEmpty);
        }
        proj = RadialSelection(pz,_sys["DStrkr5"],meanp);
        if ( sqrt(proj.X*proj.X+proj.Y*proj.Y)>meanp ) continue;
        mccuts_accept->Fill("Fiducial Selection",1);

        chi2n = FillCrossCheck();
        if ( chi2n[0]>4 && chi2n[1]>4 ) continue;
        mccuts_accept->Fill("chi2ndf",1);

        if ( !PIDSelection(false) ) continue;
        mccuts_accept->Fill("TOF Selection",1);
        double diff = 1;
        USholdvar = USPreRadSelEmpty.E(0);
        pz = CorMomFromTOF(USholdvar, diff, false);

        FillCollectionSciFi(_USMCset, jUS, kUS, lUS, pid, pz, tracker0, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        FillCollectionSciFi(_DSMCset, jDS, kDS, lDS, pid, pz, tracker1, projdif, proj, chi2n, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            MCUSEmpty.append_instance(mcUSTruthAbsHitEmpty);
            MCDSEmpty.append_instance(mcDSTruthAbsHitEmpty);
        }
        mccalc_mom->Fill(pz);

    }
    data_make_beam_histograms(USAllTOFEmpty, "Upstream, Data Reference", "dataUSref_alltof", MCUSAllTOFEmpty);
    data_make_beam_histograms(DSAllTOFEmpty, "Downstream, Data Reference", "dataDSref_alltof", MCDSAllTOFEmpty);
    data_make_beam_histograms(USPreRadSelEmpty, "Upstream, Data Reference", "dataUSref_prerad", MCUSPreRadSelEmpty);
    data_make_beam_histograms(DSPreRadSelEmpty, "Downstream, Data Reference", "dataDSref_prerad", MCDSPreRadSelEmpty);
    data_make_beam_histograms(_USMCset, "Upstream, Data Reference", "dataUSref", MCUSEmpty);
    data_make_beam_histograms(_DSMCset, "Downstream, Data Reference", "dataDSref", MCDSEmpty);
    make_acceptance_histograms(_USMCset, _DSMCset, "Data Reference", "dataProj");
}


/* void MCSAnalysis::generateMCSResponse(){ */

/*   chain->SetBranchAddress("TOFBranch", &tofevent); */
/*   chain->SetBranchAddress("SciFiBranch", &scifievent); */
/*   chain->SetBranchAddress("CkovBranch", &ckovevent); */
/*   chain->SetBranchAddress("KLBranch", &klevent); */
/*   chain->SetBranchAddress("EMRBranch", &emrevent); */
/*   chain->SetBranchAddress("MCEvent", &mcevent); */
/*   // Loop over the training sample. */
/*   Collection USRecSet, DSRecSet, USVirtSet, DSVirtSet; */
/*   Collection USAllTOF, DSAllTOF, USPreRadSel, DSPreRadSel; */
/*   Collection USVAllTOF, DSVAllTOF, USVPreRadSel, DSVPreRadSel; */
/*   // Collection USeRecSet, DSeRecSet, USmuRecSet, DSmuRecSet; */
/*   int nMuAbsSel=0, nEAbsSel=0, nMuAbsAll=0, nEAbsAll=0; */
/*   TH1D* mom_resUS = new TH1D("mom_resUS", ";p_{rec} - p_{MC} (MeV/c)", */
/* 			     2000, -100, 100); */
/*   TH1D* mom_resDS = new TH1D("mom_resDS", ";p_{rec} - p_{MC} (MeV/c)", */
/* 			     2000, -100, 100); */
/*   TH2D* mom_responseUS = new TH2D("mom_responseUS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)", */
/* 				  300, 0, 300, 300, 0, 300); */
/*   TH2D* mom_responseDS = new TH2D("mom_responseDS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)", */
/* 				  300, 100, 300, 300, 0, 300); */
/*   TH2D* mom_responseABS = new TH2D("mom_responseABS", ";p_{rec} (MeV/c); p_{MC} (MeV/c)", */
/* 				  300, 100, 300, 300, 0, 300); */
/*   int ngood=0; */
/*   int 100000 = chain->GetEntries(); */
/*   for (int j=0; j<10000; j++){ */
/*     chain->GetEntry(j); */
/*     //if (j%100000==0) std::cout<<"MC Event "<<j<<", selected "<<ngood<<" events.\n"; */ 
/*     // if(fabs(mcevent->GetPrimary()->GetParticleId()) != 13) continue; */
/*     // Select events that produce virtual plane hits and pass the data selection cuts. */
/*     Vars USAbsHit, DSAbsHit, USTrackerRefHit, DSTrackerRefHit; */
/*     double pz = 0.; */
/*     mccuts_accept->Fill("All Events",1); */
/*     if ( !PIDSelection(false) ) continue; // event_ok=false; */
/*     if( !findVirtualPlanes() ) continue; */
/*     mccuts_accept->Fill("Found Virtual Planes",1); */
/*     int centre = USabsPlaneI + 1; */
/*     FillVarsVirtual(USAbsHit, USabsPlaneI, centre); */
/*     FillVarsVirtual(DSAbsHit, DSabsPlaneI, centre); */
/*     FillVarsVirtual(USTrackerRefHit, USrefplaneI, centre); */
/*     FillVarsVirtual(DSTrackerRefHit, DSrefplaneI, centre); */
/*     if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13 && */
/* 	mcevent->GetVirtualHits()->at(DSrefplaneI).GetParticleId()==-11) */
/*       nEAbsAll++; */
/*     if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13) */
/*       nMuAbsAll++; */
/*     // Apply selection cuts as with the data */
/*     bool event_ok=true; */
/*     if ( !MatchUSDS() ) { */
/*       if (jUS == -1 || kUS == -1){ */
/* 	cuts_accept->Fill("US Tracker Found", 1); */
/* 	event_ok=false; */
/*       } */ 
/*     } */
/*     if (event_ok) mccuts_accept->Fill("US Track Found",1); */

/*     if (event_ok) pz = MomentumFromTOF(false); */
/*     if (event_ok) mccuts_accept->Fill("TOF Selection",1); */
/*     // if (scifievent->scifitracks().size() != 2) continue; */
/*     FillCollectionSciFi(USPreRadSel, jUS, kUS, lUS, pz, 0); */
/*     FillCollectionSciFi(DSPreRadSel, jDS, kDS, lDS, pz, 1); */
/*     USVPreRadSel.append_instance(USAbsHit); */
/*     DSVPreRadSel.append_instance(DSAbsHit); */
/*     if ( !RadialSelection(pz, 19948.8,meanp) ) event_ok=false; */
/*     if (event_ok) mccuts_accept->Fill("Fiducial Selection",1); */
/*     if ( !RadialSelection(pz,_sys["diffpos"],90) ) event_ok=false; */
/*     if (event_ok) mctrue_mom->Fill(mcevent->GetVirtualHits()->at(USrefplaneI).GetMomentum().z()); */
/*     if (event_ok) mctrue_mom->Fill(mcevent->GetVirtualHits()->at(USrefplaneI).GetMomentum().z()); */
/*     Vars USSciFiRec, DSSciFiRec; */
/*     if (event_ok){ */
/*       FillVarsSciFi(USSciFiRec, jUS, kUS, lUS, pz, 0); */
/*       FillVarsSciFi(DSSciFiRec, jDS, kDS, lDS, pz, 1); */
/*       USRecSet.append_instance(USSciFiRec); */
/*       DSRecSet.append_instance(DSSciFiRec); */
/*       if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13 && */
/* 	  mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-11) */
/* 	nEAbsSel++; */
/*       if (mcevent->GetVirtualHits()->at(USrefplaneI).GetParticleId()==-13) */
/* 	nMuAbsSel++; */

/*       mom_resUS->Fill(USSciFiRec.pz - USTrackerRefHit.pz); */
/*       mom_resDS->Fill(DSSciFiRec.pz - DSTrackerRefHit.pz); */
/*       mom_responseUS->Fill(USSciFiRec.pz, USTrackerRefHit.pz); */
/*       mom_responseDS->Fill(DSSciFiRec.pz, DSTrackerRefHit.pz); */
/*       mom_responseABS->Fill(DSSciFiRec.pz, USAbsHit.pz); */
/*     } */
/*     USVirtSet.append_instance(USAbsHit); */
/*     DSVirtSet.append_instance(DSAbsHit); */
/*     FillMCSResponse(event_ok, USSciFiRec, DSSciFiRec, USAbsHit, DSAbsHit); */
/*     FillMuScattResponse(event_ok, USSciFiRec, DSSciFiRec, USAbsHit, DSAbsHit); */
/*     ngood++; */
/*   } */
/*   std::cout<<"For all events simulated there are "<<nEAbsAll<<" decay electrons and "<<nMuAbsAll<<" muons.\n"; */
/*   std::cout<<"For the selected MC events there are "<<nEAbsSel<<" decay electrons and "<<nMuAbsSel<<" muons.\n"; */
/*   make_beam_histograms(USVirtSet, "Upstream, Data", "VirtMCUS"); */
/*   make_beam_histograms(DSVirtSet, "Downstream, Data", "VirtMCDS"); */
/*   //make_beam_histograms(USAllTOF, "Upstream, Reconstructed Simulation", "recMCUS_alltof"); */
/*   //make_beam_histograms(DSAllTOF, "Downstream, Reconstructed Simulation", "recMCDS_alltof"); */
/*   make_beam_histograms(USPreRadSel, "Upstream, Reconstructed Simulation", "recMCUS_prerad"); */
/*   make_beam_histograms(DSPreRadSel, "Downstream, Reconstructed Simulation", "recMCDS_prerad"); */
/*   make_scattering_acceptance_histograms(USVPreRadSel, */
/* 					DSVPreRadSel, DSPreRadSel, */
/* 					"Virtual Projection","VirtPreRad"); */

/*   make_beam_histograms(USRecSet, "Upstream, Reconstructed Simulation", "recMCUS"); */
/*   make_beam_histograms(DSRecSet, "Downstream, Reconstructed Simulation", "recMCDS"); */
/*   make_scattering_acceptance_histograms(USVirtSet, */
/* 					DSVirtSet,DSRecSet, */
/* 					"Virtual Projection","VirtProj"); */
/*   mom_resUS->Write(); */
/*   mom_resDS->Write(); */
/*   mom_responseUS->Write(); */
/*   mom_responseDS->Write(); */
/*   mom_responseABS->Write(); */
/* } */

void MCSAnalysis::TruthData(int mode){

    mcchain->SetBranchAddress("MCEvent", &mcevent);
    mcemptychain->SetBranchAddress("MCEvent", &mcevent);
    Collection USTruthzero, DSTruthzero;
    Vars USTruthzeroHit, DSTruthzeroHit;
    int ngood=0;
    bool event_ok=true;
    int Nentries;
    std::cout << mode << std::endl;
    if (mode == 5 || mode == 4 || mode == 10) {
        Nentries = mcemptychain->GetEntries();
    }
    else {
        Nentries = mcchain->GetEntries();
    }

    std::cout << "start of truth" << std::endl;
    std::cout << Nentries << std::endl;
    int zz=0;
    for (int j=0; j<Nentries; j++){
        if (mode == 5 || mode == 4 || mode == 10) {
            mcemptychain->GetEntry(j);
        }
        else {
            mcchain->GetEntry(j);
        }

        if (j%100000==0) std::cout<<"MC Event "<<j<<", selected "<<ngood<<" events.\n"; 

        if( !findVirtualPlanes() ) continue;
        if ( !TruthMatchUSDS() ) {
            if (jUS == -1 || kUS == -1){
                continue;
            }
        }

        double pz = 0.;

        if ( !TruthRadialSelection(pz, _sys["diffpos"],difradius,USrefplaneI) ) continue;
        if ( !TruthRadialSelection(pz, _sys["DStrkr5"],meanp,USrefplaneI) ) continue;
        if ( !TruthTime(true) ) continue; // event_ok=false;

        Vars USTruthAbsHit, DSTruthAbsHit;
        int centre = USabsPlaneI + 1;
        FillVarsVirtual(USTruthAbsHit, USabsPlaneI, centre);
        FillVarsVirtual(DSTruthAbsHit, DSabsPlaneI, centre);
        if (mode == 5 || mode == 4 || mode == 10) {
            USTruthSet.append_instance(USTruthAbsHit);
            DSTruthSet.append_instance(DSTruthAbsHit);
        }
        else {
            USTruthSetLiH.append_instance(USTruthAbsHit);
            DSTruthSetLiH.append_instance(DSTruthAbsHit);
        }
        eventnovspz->Fill(USTruthAbsHit.pz,USTruthSet.N());

        zz++;

    }
    if (mode == 5 || mode == 4 || mode == 10) {
        TruthGraph(USTruthSet, DSTruthSet, false);
    }
    else {
        TruthGraph(USTruthSetLiH, DSTruthSetLiH, false);
    }
    /* trkr_eff(USTruthSet,DSTruthSet,_USset,_DSset,"Sel_acc","Sel_acc"); */
    make_beam_histograms(USTruthSet, "Upstream, Truth", "TruthUS");
    make_beam_histograms(DSTruthSet, "Downstream, Truth", "TruthDS");
    make_beam_histograms(USTruthSetLiH, "Upstream, Truth LiH", "LiHTruthUS");
    make_beam_histograms(DSTruthSetLiH, "Downstream, Truth LiH", "LiHTruthDS");
}

bool MCSAnalysis::MatchUSDS(){

    bool trackmatched = true;
    jUS=-1;
    jDS=-1;
    kUS=-1;
    kDS=-1;
    lUS=-1;
    lDS=-1;

    /* int notrackUS = 0; */
    /* int notrackDS = 0; */
    /* for ( size_t j=0; j<->scifitracks()scifievent->scifitracks().size(); j++){ */
    /*     double maxUS=0.0, minDS=44000; */
    /*     int tracker = scifievent->scifitracks()[j]->tracker(); */
    /*     if(tracker==0){ */
    /*         notrackUS += 1; */
    /*     } */
    /*     if(tracker==1){ */
    /*         notrackDS += 1; */
    /*     } */
    /* } */
    /* trackno->Fill(notrackUS,notrackDS); */

    if( scifievent->scifitracks().size() == 1 || 
            scifievent->scifitracks().size() == 2){
        for ( size_t j=0; j<scifievent->scifitracks().size(); j++){
            int npoints = scifievent->scifitracks()[j]->scifitrackpoints().size();
            double maxUS=0.0, minDS=44000;
            int tracker = scifievent->scifitracks()[j]->tracker();
            /* std::cout << tracker << std::endl; */
            /* std::cout << npoints << std::endl; */
            for ( int k=0; k<npoints; k++){
                double zpos = 
                    scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
                /* std::cout << zpos << std::endl; */
                if(tracker==0 && zpos > maxUS){
                    maxUS = zpos;
                    kUS = k;
                    jUS = j;
                }
                if(tracker==0 && zpos < minDS){
                    lUS = k;
                }
                if(tracker==1 && zpos < minDS){      
                    minDS = zpos;
                    kDS = k;
                    jDS = j;
                }
                if(tracker==1 && zpos > maxUS){      
                    lDS = k;
                }
            }
            /* if (jUS != -1 && kUS != -1 && jDS != -1 && kDS != 1) { */
            /*     USrefplaneZ = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z(); */
            /*     DSrefplaneZ = minDS; */
            /*     // std::cout<<"Identified US track by ["<<jUS<<", "<<kUS<<"] at z = "<<USrefplaneZ<<" mm and DS track by ["<<jDS<<", "<<kDS<<"] at z = "<<DSrefplaneZ<<" mm."<<std::endl; */ 
            /* } */
        }
        if (jUS == -1 || kUS == -1) {
            //std::cout<<"Failed US track by ["<<jUS<<", "<<kUS<<"] and DS track by ["<<jDS<<", "<<kDS<<"]"<<std::endl; 
            /* std::cout << scifievent->scifitracks().size() << std::endl; */
            trackmatched = false;
        }
    }


    else trackmatched = false;
    return trackmatched;  
}

bool MCSAnalysis::TruthMatchUSDS(){

    bool trackmatched = true;
    jUS=-1;
    jDS=-1;
    kUS=-1;
    kDS=-1;

    for (int jj=0; jj<mcevent->GetVirtualHitsSize();jj++) {
        int station_id = mcevent->GetVirtualHits()->at(jj).GetStationId();
        if (station_id == 47){
            kUS = 1;
            jUS = 1;
        }

        if (station_id == 49){
            kDS = 1;
            jDS = 1;
        }
    }
    if (jUS == -1 || kUS == -1) {
        trackmatched = false;
    }
    return trackmatched;  
}

std::vector<float> MCSAnalysis::FillCrossCheck(){

    bool stop = true;
    float chi2n;
    std::vector<float> chi2nv;
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
                    chi2n = scifievent->scifitracks()[j]->chi2()/scifievent->scifitracks()[j]->ndf();
                    chi2nv.push_back(chi2n);
                }
                /* if(tracker==1 && zpos < minDS){ */      
                /*     chi2n = scifievent->scifitracks()[j]->chi2()/scifievent->scifitracks()[j]->ndf(); */
                /*     chi2nv.push_back(chi2n); */
                /* } */
            }
        }
        //nocluster->Fill();a
    }  
    return chi2nv;
}

bool MCSAnalysis::PIDSelection(const bool& isdata=true){

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
    if ( tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size() == 1){
        rawTOF0HitTime = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
    }
    else{
        return false;
    }
    if ( tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size() == 1){
        rawTOF1HitTime = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
    }
    else {
        return false;
    }

    if (  tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArraySize() == 1 )
        rawTOF2HitTime = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
    else // Allow for TOF2 to not be hit
        rawTOF2HitTime = rawTOF1HitTime + 100.0 * 8.22475 / 0.299792458; // ns. 
    if(isdata){

        if ( rawTOF1HitTime - rawTOF0HitTime < TOF_lower_limit ||
                rawTOF1HitTime - rawTOF0HitTime > TOF_upper_limit) return false;
 
        /* if (tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArraySize() == 1 ) TOF1Hit += 1; */
        /* if (tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArraySize() == 1 )  TOF2Hit += 1; */
        /* percentTOF2 = TOF2Hit/TOF1Hit; */
        /* double pcor = MomCalc(rawTOF1HitTime - rawTOF0HitTime); */
        /* // if ( pcor < P_lower_limit || */
        //      pcor > P_upper_limit) return false;
    }
    else {
        mctof10->Fill(rawTOF1HitTime - rawTOF0HitTime);
        mctof21->Fill(rawTOF2HitTime - rawTOF1HitTime);
        /* if ( rawTOF1HitTime - rawTOF0HitTime > TOF_lower_limit+0.1 && */
        /* 	rawTOF1HitTime - rawTOF0HitTime < TOF_upper_limit+0.1) { */
        if ( rawTOF1HitTime - rawTOF0HitTime > TOF_lower_limit_ref &&
                rawTOF1HitTime - rawTOF0HitTime < TOF_upper_limit_ref) {
            mctof10_sel->Fill(rawTOF1HitTime - rawTOF0HitTime);
            mctof21_sel->Fill(rawTOF2HitTime - rawTOF1HitTime);
        }
        if ( rawTOF1HitTime - rawTOF0HitTime < TOF_lower_limit_ref ||
                rawTOF1HitTime - rawTOF0HitTime > TOF_upper_limit_ref) {
            return false;
        }
        /* double pcor = MomCalc(rawTOF1HitTime - rawTOF0HitTime); */
        // if ( pcor < P_lower_limit ||
        //      pcor > P_upper_limit) return false;
    }
    return true;

    }

    /* bool MCSAnalysis::PSelection(double p){ */

    /* 	if ( p < p_lower_limit || p > p_upper_limit) return false; */
    /* 	return true; */

    /* } */

    /* bool MCSAnalysis::PSelection(double p){ */

    /* 	if ( p < p_lower_limit || p > p_upper_limit) return false; */
    /* 	return true; */

    /* } */

    bool MCSAnalysis::TruthTime(const bool& isdata=true){

        MAUS::TOFHitArray* tofhit = mcevent->GetTOFHits();
        float tof0;
        float tof1;
        for (int i=0;i<tofhit->size();i++) {
            if (tofhit->at(i).GetPosition().Z() - 5287 < 20 and tofhit->at(i).GetPosition().Z() - 5287 > -20){  
                tof0 = tofhit->at(i).GetTime();   
            }
            if (tofhit->at(i).GetPosition().Z() - 12929 < 20 and tofhit->at(i).GetPosition().Z() - 12929 > -20){        tof1 = tofhit->at(i).GetTime();
            }
        }
        if ((tof1-tof0) > TOF_upper_limit || (tof1-tof0) < TOF_lower_limit) return false; 
        /* double pz = mcevent->GetVirtualHits()->at(centre).GetMomentum().z(); */
        /* if (pz > P_ul || pz < P_ll) return false; */ 
        return true;
    }

    bool MCSAnalysis::TruthTime12(const bool& isdata){

        MAUS::TOFHitArray* tofhit = mcevent->GetTOFHits();
        float tof2;
        float tof1;
        for (int i=0;i<tofhit->size();i++) {
            if (tofhit->at(i).GetPosition().Z() - 21183 < 20 and tofhit->at(i).GetPosition().Z() - 21138 > -20){  
                tof2 = tofhit->at(i).GetTime();   
            }
            if (tofhit->at(i).GetPosition().Z() - 12929 < 20 and tofhit->at(i).GetPosition().Z() - 12929 > -20){        tof1 = tofhit->at(i).GetTime();
            }
        }
        if(isdata){
            if ((tof2-tof1) > TOF_upper_limit || (tof2-tof1) < TOF_lower_limit) return false; 
        }
        else {
            if ((tof2-tof1) > TOF_upper_limit || (tof2-tof1) < TOF_lower_limit) return false; 
        }
        return true;
    }

    Vars MCSAnalysis::RadialSelection(double& pz, double& pos, double& radius){
        Vars USplane;
        USplane.X = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().x() + _sys["alXUS"];
        USplane.Y = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().y() + _sys["alYUS"];
        /* std::cout << "before" << std::endl; */
        /* std::cout << "X " << USplane.X << std::endl; */
        /* std::cout << "Y " << USplane.Y << std::endl; */
        USplane.Z = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
        USplane.px = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().x() +
            tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
        USplane.py = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
            tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
        USplane.pz = pz;// * sqrt(1 + pow(USplane.py, 2) + pow(USplane.px, 2));
        // if ( sqrt(xpos*xpos + ypos*ypos) > meanp) selected = false;
        USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
            tan(_sys["thXUS"] * atan(1.)/45.0);
        USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
            tan(_sys["thYUS"] * atan(1.)/45.0);

        double abspos = _sys["abspos"];
        double phi = atan2(USplane.dYdz, USplane.dXdz);
        double zdiff = 2*fabs(USplane.Z - pos);
        Vars DSproj;
        if (pos == _sys["diffpos"]) {
            DSproj = PropagateVarsMu(USplane, pos);
        }
        else {
            DSproj = PropagateVarsMu(USplane, _sys["abspos"]);
            //USplane.dXdz += sigmap*cos(phi);
            //USplane.dYdz += sigmap*sin(phi);
            //USplane.px   += sigmap*cos(phi)*USplane.pz;
            //USplane.py   += sigmap*sin(phi)*USplane.pz;
            //double xabs  = (dXdz + sigmap*cos(phi)) * zdiff + xpos;
            //double yabs  = (dYdz + sigmap*sin(phi)) * zdiff + ypos;
            //std::cout << "USplane.px " << USplane.px << std::endl;
            //USplane.px -= 10.2384;
            //std::cout << "USplane.px " << USplane.px << std::endl;
            DSproj = PropagateVarsMu(USplane, pos);
            //DSproj.X += 5.402;
            //DSproj.dXdz -= 0.02384;
            //DSproj.dXdz = DSproj.dXdz*cos(-0.002384)-sin(-0.002384);
        }
        /* std::cout << "after" << std::endl; */
        /* std::cout << "X " << DSproj.X << std::endl; */
        /* std::cout << "Y " << DSproj.Y << std::endl; */
        return DSproj;
    }

    bool MCSAnalysis::TruthRadialSelection(double& pz, double& pos, double& radius, int& j){
        bool selected = true;
        if (jUS == -1 || kUS == -1) 
            selected = false;
        else {
            Vars USplane;
            USplane.X = mcevent->GetVirtualHits()->at(j).GetPosition().x();
            USplane.Y = mcevent->GetVirtualHits()->at(j).GetPosition().y();
            USplane.Z = mcevent->GetVirtualHits()->at(j).GetPosition().z();
            USplane.dXdz = mcevent->GetVirtualHits()->at(j).GetMomentum().x()/mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            USplane.dYdz = mcevent->GetVirtualHits()->at(j).GetMomentum().y()/mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            USplane.pz   = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            USplane.px   = mcevent->GetVirtualHits()->at(j).GetMomentum().x();
            USplane.py   = mcevent->GetVirtualHits()->at(j).GetMomentum().y();
            double phi = atan2(USplane.dYdz, USplane.dXdz);
            double zdiff = 2*fabs(USplane.Z - pos);
            Vars DSproj;
            if (pos == _sys["diffpos"]) {
                DSproj = PropagateVarsMu(USplane, pos);
            }
            else {
                DSproj = PropagateVarsMu(USplane, _sys["abspos"]);
                //USplane.dXdz += sigmap*cos(phi);
                //USplane.dYdz += sigmap*sin(phi);
                //USplane.px   += sigmap*cos(phi)*USplane.pz;
                //USplane.py   += sigmap*sin(phi)*USplane.pz;
                DSproj = PropagateVarsMu(USplane, pos);
            }

            if ( sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y) > radius) selected = false;
            /* if (pos == 19948.8) projradius->Fill(sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y)); */
            /* if (pos == _sys["diffpos"]) diffradius->Fill(sqrt(DSproj.X*DSproj.X + DSproj.Y*DSproj.Y)); */
        }
        return selected;
    }

    std::vector<double> MCSAnalysis::DefineProjectionAngles(Vars& US, Vars& DS){

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

    std::vector<double> MCSAnalysis::RotDefineProjectionAngles(const Vars& US, const Vars& DS,int i){

        std::vector<double> projTheta;
        double USnorm = 1./sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz);
        TVector3 u(US.dXdz*USnorm, US.dYdz*USnorm, USnorm);
        //double y = -cos(i*3.14/180);
        double x = sin(i*3.14/180);
        TVector3 s(0, -1, 0);
        s.RotateZ(i*TMath::Pi()/180);
        double snorm = 1./sqrt(s[2]*s[2] + s[1]*s[1] + s[0]*s[0]);
        TVector3 s_(s[0]*snorm,s[1]*snorm,s[2]*snorm);
        TVector3 v = s.Cross(u);
        double vnorm = 1./sqrt(v[2]*v[2] + v[1]*v[1] + v[0]*v[0]); 
        TVector3 v_(v[0]*vnorm,v[1]*vnorm,v[2]*vnorm);
        TVector3 w = v_.Cross(u);
        double Wnorm  = 1./sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
        w[0] *= Wnorm;
        w[1] *= Wnorm;
        w[2] *= Wnorm;
        double DSnorm = 1./sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz);
        TVector3 d(DS.dXdz*DSnorm, DS.dYdz*DSnorm, DSnorm);
        //projTheta.push_back(atan(d.Dot(w)/d.Mag()*w.Mag()));
        //projTheta.push_back(atan(d.Dot(v_)/v_.Mag()*d.Mag()));
        double hold = (atan(d.Dot(w)/d.Dot(u)));
        projTheta.push_back(hold);
        //std::cout << "hold " << hold << std::endl;
        projTheta.push_back(atan(d.Dot(v_)/d.Dot(u)));

        /*
           if (i==45 || i==42 || i ==0){
           std::cout << "i " << i << std::endl;
           std::cout << "w[0] " << w[0] << std::endl;
           std::cout << "w[1] " << w[1] << std::endl;
           std::cout << "w[2] " << w[2] << std::endl;
           std::cout << "d[0] " << d[0] << std::endl;
           std::cout << "d[1] " << d[1] << std::endl;
           std::cout << "d[2] " << d[2] << std::endl;
           std::cout << "d.Dot(w) " << d.Dot(w) << std::endl;
        //std::cout << "atan(d.Dot(w)/d.Dot(u)) " << atan(d.Dot(w)/d.Dot(u)) << std::endl;
        //std::cout << "d.Dot(v_) " << d.Dot(v_) << std::endl;
        }
        */

        projTheta.push_back( acos( ( (1 + US.dXdz * DS.dXdz + US.dYdz * DS.dYdz )/
                        sqrt(1 + US.dXdz*US.dXdz + US.dYdz*US.dYdz)/
                        sqrt(1 + DS.dXdz*DS.dXdz + DS.dYdz*DS.dYdz))) );
        projTheta.push_back(d.Dot(w)); 
        projTheta.push_back(d.Dot(v_)); 
        TVector3 yvec(0,-1,0);
        //std::cout << "yvec.Angle(d) " << yvec.Angle(d) << std::endl;
        //std::cout << "w.Angle(d) " << w.Angle(d) << std::endl;
        projTheta.push_back(yvec.Angle(d)); 
        projTheta.push_back(w.Angle(d)); 
        projTheta.push_back(v_.Angle(d)); 
        projTheta.push_back(w[0]*d[0]);
        projTheta.push_back(w[1]*d[1]);
        //std::cout << "w[0] " << w[0] << std::endl;
        //std::cout << "w[1] " << w[1] << std::endl;
        projTheta.push_back(w[2]*d[2]);
        projTheta.push_back(w[0]*d[0]+w[1]*d[1]+w[2]*d[2]);
        projTheta.push_back(u[2]);
        projTheta.push_back(w.Mag());
        //std::cout << "hold " << hold << std::endl;
        projTheta.push_back(hold); 
        //std::cout << "d[0] " << d[0] << std::endl;
        //std::cout << "d[1] " << d[1] << std::endl;
        return projTheta;
    }

    std::vector<double> MCSAnalysis::mcRotDefineProjectionAngles(const Vars& US, const Vars& DS,int i){

        std::vector<double> projTheta;
        double USnorm = 1./sqrt(1 + US.mcdXdz*US.mcdXdz + US.mcdYdz*US.mcdYdz);
        TVector3 u(US.mcdXdz*USnorm, US.mcdYdz*USnorm, USnorm);
        //double y = -cos(i*3.14/180);
        double x = sin(i*3.14/180);
        TVector3 s(0, -1, 0);
        s.RotateZ(i*TMath::Pi()/180);
        double snorm = 1./sqrt(s[2]*s[2] + s[1]*s[1] + s[0]*s[0]);
        TVector3 s_(s[0]*snorm,s[1]*snorm,s[2]*snorm);
        TVector3 v = s.Cross(u);
        double vnorm = 1./sqrt(v[2]*v[2] + v[1]*v[1] + v[0]*v[0]); 
        TVector3 v_(v[0]*vnorm,v[1]*vnorm,v[2]*vnorm);
        TVector3 w = v_.Cross(u);
        double Wnorm  = 1./sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
        w[0] *= Wnorm;
        w[1] *= Wnorm;
        w[2] *= Wnorm;
        double DSnorm = 1./sqrt(1 + DS.mcdXdz*DS.mcdXdz + DS.mcdYdz*DS.mcdYdz);
        TVector3 d(DS.mcdXdz*DSnorm, DS.mcdYdz*DSnorm, DSnorm);
        //projTheta.push_back(atan(d.Dot(w)/d.Mag()*w.Mag()));
        //projTheta.push_back(atan(d.Dot(v_)/v_.Mag()*d.Mag()));
        double hold = (atan(d.Dot(w)/d.Dot(u)));
        projTheta.push_back(hold);
        //std::cout << "hold " << hold << std::endl;
        projTheta.push_back(atan(d.Dot(v_)/d.Dot(u)));

        /*
           if (i==45 || i==42 || i ==0){
           std::cout << "i " << i << std::endl;
           std::cout << "w[0] " << w[0] << std::endl;
           std::cout << "w[1] " << w[1] << std::endl;
           std::cout << "w[2] " << w[2] << std::endl;
           std::cout << "d[0] " << d[0] << std::endl;
           std::cout << "d[1] " << d[1] << std::endl;
           std::cout << "d[2] " << d[2] << std::endl;
           std::cout << "d.Dot(w) " << d.Dot(w) << std::endl;
        //std::cout << "atan(d.Dot(w)/d.Dot(u)) " << atan(d.Dot(w)/d.Dot(u)) << std::endl;
        //std::cout << "d.Dot(v_) " << d.Dot(v_) << std::endl;
        }
        */

        projTheta.push_back( acos( ( (1 + US.mcdXdz * DS.mcdXdz + US.mcdYdz * DS.mcdYdz )/
                        sqrt(1 + US.mcdXdz*US.mcdXdz + US.mcdYdz*US.mcdYdz)/
                        sqrt(1 + DS.mcdXdz*DS.mcdXdz + DS.mcdYdz*DS.mcdYdz))) );
        projTheta.push_back(d.Dot(w)); 
        projTheta.push_back(d.Dot(v_)); 
        TVector3 yvec(0,-1,0);
        //std::cout << "yvec.Angle(d) " << yvec.Angle(d) << std::endl;
        //std::cout << "w.Angle(d) " << w.Angle(d) << std::endl;
        projTheta.push_back(yvec.Angle(d)); 
        projTheta.push_back(w.Angle(d)); 
        projTheta.push_back(v_.Angle(d)); 
        projTheta.push_back(w[0]*d[0]);
        projTheta.push_back(w[1]*d[1]);
        //std::cout << "w[0] " << w[0] << std::endl;
        //std::cout << "w[1] " << w[1] << std::endl;
        projTheta.push_back(w[2]*d[2]);
        projTheta.push_back(w[0]*d[0]+w[1]*d[1]+w[2]*d[2]);
        projTheta.push_back(u[2]);
        projTheta.push_back(w.Mag());
        //std::cout << "hold " << hold << std::endl;
        projTheta.push_back(hold); 
        //std::cout << "d[0] " << d[0] << std::endl;
        //std::cout << "d[1] " << d[1] << std::endl;
        return projTheta;
    }

    double MCSAnalysis::PathLengthInLH2(double& pz){

        if (jUS != -1 && kUS != -1 && jDS != -1 && kDS != -1){
            float x_up;
            float y_up;
            float z_up;
            float x_down;
            float y_down;
            float z_down;
            float x_up_it;
            float y_up_it;
            float z_up_it;
            float x_down_it;
            float y_down_it;
            float z_down_it;
            float rad_up = 0;
            float curvey = 0;
            float curvey_down = 0;
            float rad_down = 0;
            int i = 0;

            double poly[5];
            poly[0] = 5.00929921969443e-6;
            poly[1] = -0.0006657044;
            poly[2] = 0.0338057657;
            poly[3] = -0.8148399317;
            poly[4] = 10.8671889437;

            std::vector<double> path_length;
            std::vector<double> vpath_length;

            // Path length for muon between TOF1 and absorber
            Vars USplane;
            USplane.X = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().x() + _sys["alXUS"];
            USplane.Y = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().y() + _sys["alYUS"];
            USplane.Z = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->pos().z();
            tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
            USplane.py   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
                tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
            USplane.pz   = pz; //* ;sqrt(1 + pow(USplane.py, 2) + pow(USplane.px, 2));
            USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
                tan(_sys["thXUS"] * atan(1.)/45.0);
            USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
                tan(_sys["thYUS"] * atan(1.)/45.0);
            Vars USabsoproj = PropagateVarsMu(USplane, _sys["abspos"]);
            Vars USabsfront = PropagateVarsMu(USplane, 16729.03);

            //Path length for muon between absorber and TOF2
            //vpath_length.clear();
            Vars DSplane;
            DSplane.X = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().x() + _sys["alXDS"];
            DSplane.Y = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().y() + _sys["alYDS"];
            DSplane.Z = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->pos().z();
            DSplane.px   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().x() +
                tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
            DSplane.py   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().y() +
                tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
            DSplane.pz   = pz* sqrt(1 + pow(DSplane.py, 2) + pow(DSplane.px, 2));
            DSplane.dXdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().x() +
                tan(_sys["thXDS"] * atan(1.)/45.0);
            DSplane.dYdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().y() +
                tan(_sys["thYDS"] * atan(1.)/45.0);
            Vars DSabsoproj = PropagateVarsMu(DSplane, _sys["abspos"]);
            Vars DSabsback = PropagateVarsMu(DSplane, 17179.33);

            x_up = USabsfront.X;
            y_up = USabsfront.Y;
            z_up = 16779.18;
            rad_up = sqrt(pow(USabsfront.X,2)+pow(USabsfront.Y,2));
            while (rad_up >= curvey && i < 50) {
                curvey = (poly[0] * pow(i,5)) + (poly[1] * pow(i,4)) + (poly[2] * pow(i,3)) + (poly[3] * pow(i,2)) + (poly[4] * i);
                x_up_it = USabsfront.X + (USabsfront.dXdz * (i));
                y_up_it = USabsfront.Y + (USabsfront.dYdz * (i));
                rad_up = sqrt(pow(x_up_it,2) + pow(y_up_it,2));
                if (rad_up < curvey) {
                    x_up = x_up_it;
                    y_up = y_up_it;
                    z_up = i+16779.18;
                }
                i++;
            }

            x_down = DSabsback.X;
            y_down = DSabsback.Y;
            z_down = 17129.18;
            i = 0;
            rad_down = sqrt(pow(DSabsback.X,2)+pow(DSabsback.Y,2));
            while (rad_down >= curvey_down && i < 50) {
                curvey_down = sqrt(pow(-(poly[0] * pow(i,5)) - (poly[1] * pow(i,4)) - (poly[2] * pow(i,3)) - (poly[3] * pow(i,2)) - (poly[4] * i),2));
                x_down_it = DSabsback.X + (DSabsback.dXdz * (-i));
                y_down_it = DSabsback.Y + (DSabsback.dYdz * (-i));
                rad_down = sqrt(pow(x_down_it,2) + pow(y_down_it,2));
                if (rad_down < curvey_down) {
                    x_down = x_down_it;
                    y_down = y_down_it;
                    z_down = 17129.18-i;
                }
                i++;
            }

            float dist_absorb = sqrt(pow(x_down-x_up,2)+pow(y_down-y_up,2)+pow(z_down-z_up,2)); 
            /*
               std::cout << "x_down " << x_down << std::endl;
               std::cout << "x_up " << x_up << std::endl;
               std::cout << "y_down " << y_down << std::endl;
               std::cout << "y_up " << y_up << std::endl;
               std::cout << "z_down " << z_down << std::endl;
               std::cout << "z_up " << z_up << std::endl;
               std::cout << "dist_absorb " << dist_absorb << std::endl;
               */
            pathlengthabs->Fill(dist_absorb);
        }
    }

    std::vector<double> MCSAnalysis::CalculatePathLength(double& pz){

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
        USplane.px   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().x() +
            tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
        USplane.py   = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->mom().y() +
            tan(_sys["thXUS"] * atan(1.)/45.0) * USplane.pz;
        USplane.pz   = pz;//* sqrt(1 + pow(USplane.py, 2) + pow(USplane.px, 2));
        USplane.dXdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().x() +
            tan(_sys["thXUS"] * atan(1.)/45.0);
        USplane.dYdz = scifievent->scifitracks()[jUS]->scifitrackpoints()[kUS]->gradient().y() +
            tan(_sys["thYUS"] * atan(1.)/45.0);
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
        DSplane.px   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().x() +
            tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
        DSplane.py   = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->mom().y() +
            tan(_sys["thXDS"] * atan(1.)/45.0) * DSplane.pz;
        DSplane.pz   = pz* sqrt(1 + pow(DSplane.py, 2) + pow(DSplane.px, 2));
        DSplane.dXdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().x() +
            tan(_sys["thXDS"] * atan(1.)/45.0);
        DSplane.dYdz = scifievent->scifitracks()[jDS]->scifitrackpoints()[kDS]->gradient().y() +
            tan(_sys["thYDS"] * atan(1.)/45.0);
        Vars DSabsoproj = PropagateVarsMu(DSplane, _sys["abspos"]);
        Vars DSTOF2proj = PropagateVarsMu(DSplane, _sys["TOF2_z"]);
        vpath_length.push_back(DSTOF2proj.X - DSabsoproj.X);
        vpath_length.push_back(DSTOF2proj.Y - DSabsoproj.Y);
        vpath_length.push_back(DSTOF2proj.Z - DSabsoproj.Z);
        path_length.push_back(sqrt(vpath_length[0]*vpath_length[0] + vpath_length[1]*vpath_length[1] + vpath_length[2]*vpath_length[2]));

        return path_length;

    }

    std::vector<double> MCSAnalysis::rCalculatePathLength(double& pz){

        std::vector<double> path_length;
        std::vector<double> vpath_length;

        // Path length for muon between TOF1 and absorber
        Vars TOF0;
        Vars TOF1;
        /*
           std::cout << "jUS " << jUS << std::endl;
           std::cout << "kUS " << kUS << std::endl;
           */    
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) > 0) {
            TOF0.X = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetGlobalPosX();
            TOF0.Y = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetGlobalPosY();
            TOF0.Z = tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetGlobalPosZ();
        }
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) > 0) {
            TOF1.X = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetGlobalPosX();
            TOF1.Y = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetGlobalPosY();
            TOF1.Z = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetGlobalPosZ();
        }
        /*
           std::cout << USabsoproj.X << std::endl;
           std::cout << USplane.X << std::endl;
           std::cout << USTOF1proj.X << std::endl;
           */
        vpath_length.push_back(TOF1.X - TOF0.X);
        vpath_length.push_back(TOF1.Y - TOF0.Y);
        vpath_length.push_back(TOF1.Z - TOF0.Z);
        path_length.push_back(sqrt(vpath_length[0]*vpath_length[0] + vpath_length[1]*vpath_length[1] + vpath_length[2]*vpath_length[2]));

        //Path length for muon between absorber and TOF2
        vpath_length.clear();

        return path_length;

    }

    double MCSAnalysis::MomCalc(double& TOF, const bool& isdata){

        double t_BB;
        t_BB = TOF;
        double dt0 = (_sys["TOF1_z"] - _sys["TOF0_z"]) / 0.299792458 / 1000.;
        double pzBB = 105.65/sqrt(pow(t_BB,2)/pow(dte01,2)-1);
        double Eini = sqrt(pzBB*pzBB + 105*105);

        double Iair = 85.7e-6;
        double Iscin = 64.7e-6;
        double IAl = 166e-6;
        double IHe = 41.6e-6;
        double ILH2 = 21.8e-6;
        double ILiH = 36.5e-6;

        double zair = 897.16;
        double zscin = 5.09;
        double zAl = 16e-4;
        double zHe = 11.3;
        double zLH2 = 17.5;
        double zLiH = 3.25;

        double Zair = 0.49919;
        double Zscin = 0.54141;
        double ZAl = 13;
        double ZHe = 2;
        double ZLH2 = 1;
        double ZLiH = 0.50321;

        double Aair = 1;
        double Ascin = 1;
        double AAl = 26.9815385;
        double AHe = 4.002602;
        double ALH2 = 1.008;
        double ALiH = 1;

        double Rair = 1.205e-3;
        double Rscin = 1.032;
        double RAl = 2.699;
        double RHe = 1.663e-4;
        double RLH2 = 0.0708;
        double RLiH = 0.696;

        double hwair = 0.71e-6;
        double hwscin = 21.54e-6;
        double hwAl = 32.86e-6;
        double hwHe = 0.26e-6;
        double hwLH2 = 7.64e-6;
        double hwLiH = 18.51e-6;

        double mostprobBBair = MostProbBB(pzBB,Iair,Zair,Aair,hwair,Rair,zair);
        double mostprobBBscin = MostProbBB(pzBB,Iscin,Zscin,Ascin,hwscin,Rscin,zscin);
        double mostprobBBAl = MostProbBB(pzBB,IAl,ZAl,AAl,hwAl,RAl,zAl);
        double mostprobBBHe = MostProbBB(pzBB,IHe,ZHe,AHe,hwHe,RHe,zHe);
        double mostprobBBLH2 = MostProbBB(pzBB,ILH2,ZLH2,ALH2,hwLH2,RLH2,zLH2);
        double mostprobBBLiH = MostProbBB(pzBB,ILiH,ZLiH,ALiH,hwLiH,RLiH,zLiH);

        double BBair = BetheBloch(pzBB,Iair,Zair,Aair,hwair);
        double BBscin = BetheBloch(pzBB,Iscin,Zscin,Ascin,hwscin);
        double BBAl = BetheBloch(pzBB,IAl,ZAl,AAl,hwAl);
        double BBHe = BetheBloch(pzBB,IHe,ZHe,AHe,hwHe);
        double BBLH2 = BetheBloch(pzBB,ILH2,ZLH2,ALH2,hwLH2);

        double ELair = BBair*zair*Rair;
        double ELscin = BBscin*zscin*Rscin;
        double ELAl = BBAl*zAl*RAl;
        double ELHe = BBHe*zHe*RHe;
        double ELLH2 = BBLH2*zLH2*RLH2;

        double pcor;
        double pUS;
        double pDS;
        if (material=="LiH") {
            double totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe+mostprobBBLiH;
            double Efin = Eini - totEL;
            if (isdata) {
                pcor = sqrt(pow(Efin,2)-pow(105,2)) - UScor;
            }
            if (!isdata) {
                pcor = sqrt(pow(Efin,2)-pow(105,2)) - refUScor;
            }
        }

        double t_ds;
        return pcor;
    }
    
    double MCSAnalysis::CorMomFromTOF(Vars& set, double& diff, const bool& isdata){

        double t_BB;
        if (isdata) t_BB = set.TOF01;
        else t_BB = TimeofFlight(false);
        double dt0 = (_sys["TOF1_z"] - _sys["TOF0_z"]) / 0.299792458 / 1000.;
        double pzBB = 105.65/sqrt(pow(t_BB,2)/pow(dte01,2)-1);
        double Eini = sqrt(pzBB*pzBB + 105*105);
        uncor_momsel->Fill(pzBB);

        double Iair = 85.7e-6;
        double Iscin = 64.7e-6;
        double IAl = 166e-6;
        double IHe = 41.6e-6;
        double ILH2 = 21.8e-6;
        double ILiH = 36.5e-6;
        double ICu = 322e-6;
        double IW = 727e-6;

        double zair = 897.16;
        double zscin = 5.09;
        double zAl = 16e-4;
        double zHe = 11.3;
        double zLH2 = 17.5;
        double zLiH = 3.25;
        double zCu = 0.891;
        double zW = 0.84;

        double Zair = 0.49919;
        double Zscin = 0.54141;
        double ZAl = 13;
        double ZHe = 2;
        double ZLH2 = 1;
        double ZLiH = 0.50321;
        double ZCu = 29;
        double ZW = 74;

        double Aair = 1;
        double Ascin = 1;
        double AAl = 26.9815385;
        double AHe = 4.002602;
        double ALH2 = 1.008;
        double ALiH = 1;
        double ACu = 64;
        double AW = 184;

        double Rair = 1.205e-3;
        double Rscin = 1.032;
        double RAl = 2.699;
        double RHe = 1.663e-4;
        double RLH2 = 0.0708;
        double RLiH = 0.696;
        /* double RLiH = 0.82; */
        double RCu = 8.96;
        double RW = 19.3;

        double hwair = 0.71e-6;
        double hwscin = 21.54e-6;
        double hwAl = 32.86e-6;
        double hwHe = 0.26e-6;
        double hwLH2 = 7.64e-6;
        double hwLiH = 18.51e-6;
        double hwCu = 58.27e-6;
        double hwW = 80.32e-6;

        double mostprobBBair = MostProbBB(pzBB,Iair,Zair,Aair,hwair,Rair,zair);
        /* double totEL = mostprobBBair; */
        /* double Efin = Eini - totEL; */
        /* pzBB = sqrt(pow(Efin,2)-pow(105,2)); */
        double mostprobBBscin = MostProbBB(pzBB,Iscin,Zscin,Ascin,hwscin,Rscin,zscin);
        double mostprobBBAl = MostProbBB(pzBB,IAl,ZAl,AAl,hwAl,RAl,zAl);
        double mostprobBBHe = MostProbBB(pzBB,IHe,ZHe,AHe,hwHe,RHe,zHe);
        double mostprobBBLH2 = MostProbBB(pzBB,ILH2,ZLH2,ALH2,hwLH2,RLH2,zLH2);
        /* totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe; */
        /* Efin = Eini - totEL; */
        /* pzBB = sqrt(pow(Efin,2)-pow(105,2)); */
        double mostprobBBLiH = MostProbBB(pzBB,ILiH,ZLiH,ALiH,hwLiH,RLiH,zLiH);
        double mostprobBBCu = MostProbBB(pzBB,ICu,ZCu,ACu,hwCu,RCu,zCu);
        double mostprobBBW = MostProbBB(pzBB,IW,ZW,AW,hwW,RW,zW);

        std::cout << "mostprobBBair " << mostprobBBair << std::endl;
        std::cout << "mostprobBBscin " << mostprobBBscin << std::endl;
        std::cout << "mostprobBBAl " << mostprobBBAl << std::endl;
        std::cout << "mostprobBBHe " << mostprobBBHe << std::endl;
        std::cout << "mostprobBBLiH " << mostprobBBLiH << std::endl;
        /* std::cout << "mostprobBBLH2 " << mostprobBBLH2 << std::endl; */
        /* std::cout << "mostprobBBCu " << mostprobBBCu << std::endl; */
        /* std::cout << "mostprobBBW " << mostprobBBW << std::endl; */

        double BBair = BetheBloch(pzBB,Iair,Zair,Aair,hwair);
        double BBscin = BetheBloch(pzBB,Iscin,Zscin,Ascin,hwscin);
        double BBAl = BetheBloch(pzBB,IAl,ZAl,AAl,hwAl);
        double BBHe = BetheBloch(pzBB,IHe,ZHe,AHe,hwHe);
        double BBLH2 = BetheBloch(pzBB,ILH2,ZLH2,ALH2,hwLH2);
        /*
           std::cout << "BBair " << BBair << std::endl;
           std::cout << "BBscin " << BBscin << std::endl;
           std::cout << "BBAl " << BBAl << std::endl;
           std::cout << "BBHe " << BBHe << std::endl;
           std::cout << "BBLH2 " << BBLH2 << std::endl;
           */

        double ELair = BBair*zair*Rair;
        double ELscin = BBscin*zscin*Rscin;
        double ELAl = BBAl*zAl*RAl;
        double ELHe = BBHe*zHe*RHe;
        double ELLH2 = BBLH2*zLH2*RLH2;
        /*
           std::cout << "ELair " << ELair << std::endl;
           std::cout << "ELscin " << ELscin << std::endl;
           std::cout << "ELAl " << ELAl << std::endl;
           std::cout << "ELHe " << ELHe << std::endl;
           std::cout << "ELLH2 " << ELLH2 << std::endl;
           */

        //std::cout << diff << std::endl;
        double pcor;
        double pUS;
        double pDS;
        /* if (material=="LH2" && mat==0) { */
        /* double totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe+mostprobBBLH2; */
        /* if (diff==1) totEL += mostprobBBCu + mostprobBBW; */
        /* double Efin = Eini - totEL; */
        /* pcor = sqrt(pow(Efin,2)-pow(105,2)); */
        /* } */

        if (material=="LiH") {
            double totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe+mostprobBBLiH;
            /* if (diff==1) totEL += mostprobBBCu + mostprobBBW; */
            if (diff==1) totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe;
            double Efin = Eini - totEL;
            if (isdata) {
                pcor = sqrt(pow(Efin,2)-pow(105,2)) - UScor;
                pUS = sqrt(pow(Efin,2)-pow(105,2)) - UScor;
            }
            if (!isdata) {
                pcor = sqrt(pow(Efin,2)-pow(105,2)) - refUScor;
                pUS = sqrt(pow(Efin,2)-pow(105,2)) - refUScor;
                /* pcor = sqrt(pow(Efin,2)-pow(105,2)); */
                /* pUS = sqrt(pow(Efin,2)-pow(105,2)); */
            }
        }
        /* if (mat==1) { */
        /* double totEL = mostprobBBair+mostprobBBscin+mostprobBBAl+mostprobBBHe; */
        /* if (diff==1) totEL += mostprobBBCu + mostprobBBW; */
        /* double Efin = Eini - totEL; */
        /* pcor = sqrt(pow(Efin,2)-pow(105,2)); */
        /* } */

        /* std::cout << "pzBB " << pzBB << std::endl; */
        /* std::cout << "pcor " << pcor << std::endl; */

        double MCTruth_pz_mid = 0;
        double MCTruth_pz_up = 0;
        double MCTruth_pz_down = 0;
        double TOF0E;
        double TOF1E;
        double TOF2E;
        for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
            energyloss->Fill(mcevent->GetVirtualHits()->at(j).GetPosition().z(),sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2)));
            if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-_sys["abspos"]<10 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-_sys["absopos"]>-10) {
                MCTruth_pz_mid = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            }
            if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-16803.7<10 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-16803.7>-10) {
                MCTruth_pz_up = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            }
            if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-17101.3<10 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-17101.3>-10) {
                MCTruth_pz_down = mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            }
            if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-5000<10 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-5000>-10) {
                TOF0E = sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2));
            }
            if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-12500<10 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-12500>-10) {
                TOF1E = sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2));
            }
            if (mcevent->GetVirtualHits()->at(j).GetPosition().z()-21039<10 && mcevent->GetVirtualHits()->at(j).GetPosition().z()-21039>-10) {
                TOF2E = sqrt(pow(mcevent->GetVirtualHits()->at(j).GetMomentum().z(),2)+pow(105,2));
            }

        }
        if (isdata) MCTruth_pz_mid = set.ptruth;
        TOF0Energy->Fill(TOF0E);
        TOF1Energy->Fill(TOF1E);
        TOF2Energy->Fill(TOF2E);
        /* if (MCTruth_pz_up != 0 && MCTruth_pz_down != 0){ */
        /*     double true_delta = MCTruth_pz_up - MCTruth_pz_down; */	   
        /*     double MCTruth_pz = (MCTruth_pz_up+MCTruth_pz_down)/2; */
        /* } */
        double res01 = pcor - MCTruth_pz_mid;

        double t_ds;
        double res12 = -999;
        if (isdata) t_ds = set.TOF12;
        else t_ds = TimeofFlight12(false);
        if (t_ds != 10000) {
            dt0 = (_sys["TOF2_z"] - _sys["TOF1_z"]) / 0.299792458 / 1000.;
            if (isdata){
            pcor = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1) - DScor;
            pDS = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1) - DScor;
            /* pcor = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1); */
            /* pDS = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1); */
            }
            if (!isdata){
                /* std::cout << "refDScor " << refDScor << std::endl; */
                pcor = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1) - refDScor;
            pDS = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1) - refDScor;
                /* std::cout << "pcor " << pcor << std::endl; */
                /* pcor = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1); */
            /* pDS = 105.65/sqrt(pow(t_ds,2)/pow(dte12,2)-1); */
                /* std::cout << "pcor " << pcor << std::endl; */
            }
            res12 = pcor - MCTruth_pz_mid;
            //t_ds = t_ds*0.299792458;
            //pcor = 105.65*(8.209)/sqrt(pow(t_ds,2)-pow((8.209),2));
            /* std::cout << "p12 " << pcor << std::endl; */
        }
        /* std::cout << "pTruth " << MCTruth_pz_mid << std::endl; */

        /* TOFvsMCTruth->Fill(pcor, MCTruth_pz_mid); */
        double resUD = pUS - pDS;
        double res = pcor - MCTruth_pz_mid;

        if(isdata){
            residual->Fill(res);
            residual12->Fill(res12);
            residual01->Fill(res01);
            residualUD->Fill(resUD);
            TOFvsMCTruth->Fill(pcor,MCTruth_pz_mid);
            TOF01vsTOF12->Fill(pUS,pDS);
            TOF01vsMCTruth->Fill(pUS,MCTruth_pz_mid);
            MCTruth->Fill(MCTruth_pz_mid);
        }
        if(!isdata){
            refresidual->Fill(res);
            refresidual12->Fill(res12);
            refresidual01->Fill(res01);
            refresidualUD->Fill(resUD);
            refTOFvsMCTruth->Fill(pcor,MCTruth_pz_mid);
            refTOF01vsTOF12->Fill(pUS,pDS);
            refMCTruth->Fill(MCTruth_pz_mid);
        }

        PathLengthInLH2(pcor);
        return pcor;
    }


    double MCSAnalysis::WriteOffset(bool file){
    if (isMC==1) TOF_lower_limit = TOF_lower_limit-0.05; 
    if (isMC==1) TOF_upper_limit = TOF_upper_limit-0.05; 
    if (isMC==1) TOF_lower_limit_ref = TOF_lower_limit_ref-0.05; 
    if (isMC==1) TOF_upper_limit_ref = TOF_upper_limit_ref-0.05; 
    if (isMC==1) TOF_lower_limit = TOF_lower_limit*dte01data/dte01MC; 
    if (isMC==1) TOF_upper_limit = TOF_upper_limit*dte01data/dte01MC; 
    if (isMC==1) TOF_lower_limit_ref = TOF_lower_limit_ref*dte01data/dte01MC; 
    if (isMC==1) TOF_upper_limit_ref = TOF_upper_limit_ref*dte01data/dte01MC; 

        stringstream ss;
        if (file==1) ss << residual12->GetMean();
        if (file==0) ss << residualUD->GetMean();
        const char* str = ss.str().c_str();
        const xmlChar* mean = xmlCharStrdup(str);
        stringstream ssref;
        if (file==1) ssref << refresidual12->GetMean();
        if (file==0) ssref << refresidualUD->GetMean();
        const char* refstr = ssref.str().c_str();
        const xmlChar* refmean = xmlCharStrdup(refstr);
        stringstream sstof;
        sstof << TOF_lower_limit;
        const char* tofstr = sstof.str().c_str();
        const xmlChar* TOF_ll_xml = xmlCharStrdup(tofstr);
        stringstream sstofref;
        sstofref << TOF_lower_limit_ref;
        const char* tofrefstr = sstofref.str().c_str();
        const xmlChar* TOF_ll_ref_xml = xmlCharStrdup(tofrefstr);

        bool countLiH = 1;
        bool countref = 1;


        xmlInitParser();
        xmlDocPtr doc=NULL;
        if (file==1) doc = xmlReadFile("/data/neutrino03/jnugent/Unfolding/mom_offset.xml", NULL, 0);
        if (file==0) doc = xmlReadFile("/data/neutrino03/jnugent/Unfolding/USDS_mom_offset.xml", NULL, 0);
        if (doc != NULL) { 
            xmlNode * a_node = xmlDocGetRootElement(doc);

            char * pEnd;
            xmlNode *cur_node = a_node->children;
            for (cur_node= a_node->children; cur_node; cur_node = cur_node->next) { 
                if (cur_node->type == XML_ELEMENT_NODE) {
                    const xmlChar* id = xmlCharStrdup("id");
                    const xmlChar* value = xmlCharStrdup("value");
                    char* val1 = (char*)xmlGetProp(cur_node, id);
                    char* val2 = (char*)xmlGetProp(cur_node, value);
                    printf("node type: Element, name: %s, with prop id: %s, prop value: %s\n", 
                            cur_node->name, val1, val2);
                    std::cout << TOF_ll_xml << std::endl;
                    if (xmlStrEqual(xmlCharStrdup(val1), TOF_ll_xml)) std::cout << "match" << std::endl;
                    if (xmlStrEqual(cur_node->name, xmlCharStrdup("LiH")) && xmlStrEqual(xmlCharStrdup(val1), TOF_ll_xml)){
                        countLiH = 0;
                    std::cout << "Writing offset" << std::endl;
                        (char*)xmlSetProp(cur_node, xmlCharStrdup("value"), mean); 
                    }
                    if (xmlStrEqual(cur_node->name, xmlCharStrdup("ref")) && xmlStrEqual(xmlCharStrdup(val1), TOF_ll_ref_xml)){
                        countref = 0;
                        (char*)xmlSetProp(cur_node, xmlCharStrdup("value"), refmean); 
                    }
                }
            }
            if (countLiH==1){
                xmlNode* node = xmlNewChild(a_node, NULL, BAD_CAST "LiH", BAD_CAST "");
                xmlNewProp(node, BAD_CAST "id", BAD_CAST TOF_ll_xml);
                xmlNewProp(node, BAD_CAST "value", BAD_CAST mean);
                xmlAddChild(a_node, node);
                xmlDocSetRootElement(doc, a_node);
            }
            if (countref==1){
                xmlNode* node = xmlNewChild(a_node, NULL, BAD_CAST "ref", BAD_CAST "");
                xmlNewProp(node, BAD_CAST "id", BAD_CAST TOF_ll_ref_xml);
                xmlNewProp(node, BAD_CAST "value", BAD_CAST refmean);
                xmlAddChild(a_node, node);
                xmlDocSetRootElement(doc, a_node);
            }
            if (file==1) xmlSaveFormatFile("/data/neutrino03/jnugent/Unfolding/mom_offset.xml",doc,1);
            if (file==0) xmlSaveFormatFile("/data/neutrino03/jnugent/Unfolding/USDS_mom_offset.xml",doc,1);
            xmlFreeDoc(doc);
        }
    }

    double MCSAnalysis::getoffset(xmlNode * a_node, bool file){

        stringstream sstof;
        sstof << TOF_lower_limit;
        const char* tofstr = sstof.str().c_str();
        const xmlChar* TOF_ll_xml = xmlCharStrdup(tofstr);
        stringstream sstofref;
        sstofref << TOF_lower_limit_ref;
        const char* tofrefstr = sstofref.str().c_str();
        const xmlChar* TOF_ll_ref_xml = xmlCharStrdup(tofrefstr);
        char *fDScor = "0";
        char *frefDScor= "0";

        char * pEnd;
        char * pEnd2;
        char * pEnd3;
        char * pEnd4;
        xmlNode *cur_node = NULL;
        for (cur_node= a_node->children; cur_node; cur_node = cur_node->next) { 
            if (cur_node->type == XML_ELEMENT_NODE) {
                const xmlChar* id = xmlCharStrdup("id");
                const xmlChar* value = xmlCharStrdup("value");
                char* val1 = (char*)xmlGetProp(cur_node, id);
                char* val2 = (char*)xmlGetProp(cur_node, value);
                printf("node type: Element, name: %s, with prop id: %s, prop name: %s\n", 
                        cur_node->name, val1, val2);
                if (xmlStrEqual(cur_node->name, xmlCharStrdup("LiH")) && xmlStrEqual(xmlCharStrdup(val1), TOF_ll_xml)){
                    fDScor = (char*)xmlGetProp(cur_node, value);
                }
                if (xmlStrEqual(cur_node->name, xmlCharStrdup("ref")) && xmlStrEqual(xmlCharStrdup(val1), TOF_ll_ref_xml)){
                    frefDScor = (char*)xmlGetProp(cur_node, value);
                }
            }
        }
        /* print_element_names(cur_node->children); */
        if (file==0){
        DScor = strtod(fDScor, &pEnd);
        std::cout << "DScor " << DScor << std::endl;
        refDScor = strtod(frefDScor, &pEnd2);
        std::cout << "refDScor " << refDScor << std::endl;
        }
        if (file==1){
        UScor = strtod(fDScor, &pEnd3);
        refUScor = strtod(frefDScor, &pEnd4);
        std::cout << "UScor " << UScor << std::endl;
        std::cout << "refUScor " << refUScor << std::endl;
        }
    }

    double MCSAnalysis::ReadOffset(bool file){
        xmlInitParser();
        xmlDocPtr doc=NULL;
        if (file==0) doc = xmlReadFile("/data/neutrino03/jnugent/Unfolding/mom_offset.xml", NULL, 0);
        if (file==0) std::cout << "Reading from /data/neutrino03/jnugent/Unfolding/mom_offset.xml" << std::endl;
        if (file==1) doc = xmlReadFile("/data/neutrino03/jnugent/Unfolding/USDS_mom_offset.xml", NULL, 0);
        if (file==1) std::cout << "Reading from /data/neutrino03/jnugent/Unfolding/USDS_mom_offset.xml" << std::endl;
        xmlNode * root_element = xmlDocGetRootElement(doc);
        getoffset(root_element, file);
        xmlFreeDoc(doc);
    }

    double MCSAnalysis::TimeofFlight(const bool& raw){
        double rawTOF1HitTime = -1., rawTOF0HitTime = -1.;
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) == 1)
            rawTOF1HitTime  = 
                tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) == 1)
            rawTOF0HitTime = 
                tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
        double dt  = 10000.0; // something rediculously large as a default number
        if ( rawTOF1HitTime != -1 && rawTOF0HitTime != -1 ){
            dt  = rawTOF1HitTime - rawTOF0HitTime; 
        }
        if (raw==true) rawtime->Fill(dt);
        return dt;
    }

    double MCSAnalysis::TimeofFlight12(const bool& raw){
        double rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) == 1)
            rawTOF1HitTime  = 
                tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) == 1)
            rawTOF2HitTime = 
                tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
        double dt  = 10000.0; // something rediculously large as a default number
        if ( rawTOF1HitTime != -1. && rawTOF2HitTime != -1. ){
            dt  = rawTOF2HitTime - rawTOF1HitTime; 
        }
        if (raw==true) rawtime12->Fill(dt);
        return dt;
    }

    double MCSAnalysis::BetheBloch(double& pz, double& Imat, double& Z, double& A, double& hw){

        double beta = pow(pz,2)/(pow(105.65,2)+pow(pz,2));
        double gamma = 1/sqrt(1-pow(beta,2));
        double W = 2*0.511*pow(beta,2)*pow(gamma,2)/(1+2*gamma*0.511/105.65+pow(0.511/105.65,2));
        double I = Imat;
        double density = log(hw/I)+log(beta*gamma)-1/2;
        double dEdxpre = 0.307075*Z/(pow(beta,2)*A);   
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

    double MCSAnalysis::MostProbBB(double& pz, double& Imat, double& Z, double& A, double& hw, double& R, double& z){

        double beta = pow(pz,2)/(pow(105.65,2)+pow(pz,2));
        double gamma = 1/sqrt(1-pow(beta,2));
        double I = Imat;
        double density = log(hw/I)+log(beta*gamma)-1/2;
        double E = 0.307075*Z*R*z/(2*pow(beta,2)*A);   
        double mostprobBB = E*(log(2*105.65*pow(beta,2)*pow(gamma,2)/I)+log(E/I)+0.2-pow(beta,2)-density);

        return mostprobBB;
    }


    TH1D* MCSAnalysis::trkreffix(TH1D* h1, bool muldiv){

        TFile *f=new TFile(trkreffiname.c_str());
        TGraphAsymmErrors* efficiency_scat_x = (TGraphAsymmErrors*)f->Get("Effx_graph");
        TGraphAsymmErrors* mirror_efficiency_scat_x = (TGraphAsymmErrors*)f->Get("Effx_graph");

        TCanvas* c1 = new TCanvas();

//c1->SetBottomMargin(0.15);
        /* TF1*  tophatl = new TF1("mtophatl", "(1/2)*(1+TMath::Erf([0]*((x-(-0.08+0.025))/0.025)/sqrt(1-((x-(-0.08+0.025))/0.025)*((x-(-0.08+0.025))/0.025))))",-0.08,-0.03); */
        /* TF1*  tophatr = new TF1("mtophatr", "(1/2)*(1+TMath::Erf([0]*((x-(0.08-0.025))/0.025)/sqrt(1-((x-(0.08-0.025))/0.025)*((x-(0.08-0.025))/0.025))))",0.03,0.08); */
        /* TF1*  tophat = new TF1("mtophat", "[0]",-0.03,0.03); */
        /* efficiency_scat_x->Fit("mtophatl", "R"); */
        /* efficiency_scat_x->Fit("mtophat", "R+"); */
        /* efficiency_scat_x->Fit("mtophatr", "R+"); */
        
        /* TF1*  tophatl = new TF1("mtophatl", "[0]*x*x+[1]*x+[2]",-0.07,-0.03); */
        /* TF1*  tophatr = new TF1("mtophatr", "[0]*x*x+[1]*x+[2]",0.03,0.07); */
        /* efficiency_scat_x->Fit("mtophatl", "R","",-0.07,-0.03); */
        /* efficiency_scat_x->Fit("mtophatr", "R+","",0.03,0.07); */


        outfile->cd();
        efficiency_scat_x->Write();
        double x;
        double y;
        double t = 2;

        h1->Sumw2();
        Double_t* before = h1->GetIntegral();
        /* for (int i=1; i<25; i++) { */
        /*     efficiency_scat_x->GetPoint(i,x,y); */
        /*     std::cout << i << std::endl; */
        /*     mirror_efficiency_scat_x->SetPoint(i,1,1); */
        /*     std::cout << i << std::endl; */
        /* } */
        for (int i=25; i<h1->GetNbinsX(); i++) {
            efficiency_scat_x->GetPoint(i,x,y);
            mirror_efficiency_scat_x->SetPoint(i,-x,y);
            t += 1;
        }
        TF1* mirror_tophatl = new TF1("mmirrortophatl", "[1]*pow(x,6)+[2]*pow(x,4)+[3]*pow(x,2)+[4]",-0.07,0.07);
        mirror_tophatl->SetParameters(-2.70626e+05,-1.25699e+02,-4.02662e+04,2.63579e+00,9.94541e-01);
        /* TF1* mirror_tophatl = new TF1("mmirrortophatl", "[1]*pow(x,4)+[2]*pow(x,2)+[3]",-0.07,0.07); */
        efficiency_scat_x->Draw();
        efficiency_scat_x->Fit("mmirrortophatl", "R","",-0.07,0.07);
        /* efficiency_scat_x->Draw(); */
        /* TF1*  mirror_tophatl = new TF1("mmirrortophatl", "[0]*x*x+[1]*x+[2]",-0.07,-0.03); */
        /* mirror_efficiency_scat_x->Fit("mmirrortophatl", "R","",-0.07,-0.03); */
        /* mirror_efficiency_scat_x->Draw(); */
        /* c1->SaveAs("mirror_Effx_graph.pdf"); */
        /* c1->Clear(); */
        efficiency_scat_x->GetYaxis()->SetTitle("Acceptance");
        efficiency_scat_x->GetYaxis()->SetTitleSize(0.06);
        efficiency_scat_x->GetXaxis()->SetTitle("#it{#theta_{X}} (mrad)");
        efficiency_scat_x->GetXaxis()->SetTitleSize(0.06);
        efficiency_scat_x->GetXaxis()->SetRangeUser(-0.0465,0.0465);
        efficiency_scat_x->Draw("AP");
        mirror_tophatl->Draw("SAME");
    TText t1 = TText(0.18,0.305,"MICE Preliminary");
    TText t3 = TText(0.18,0.27,"ISIS cycle 2015/04");
    TText t2 = TText(0.18,0.235,"LiH, MAUS v3.3.2");
    t1.SetNDC(1);
    t1.SetTextSize(0.04);
    t1.SetTextFont(42);
    t2.SetNDC(1);
    t2.SetTextSize(0.03);
    t2.SetTextFont(42);
    t3.SetNDC(1);
    t3.SetTextSize(0.04);
    t3.SetTextFont(42);
            t1.Draw("SAME");
            t2.Draw("SAME");
            t3.Draw("SAME");
        /* double zero = mirror_tophatl->GetParameter(0); */
        /* double one = mirror_tophatl->GetParameter(1); */
        /* double two = mirror_tophatl->GetParameter(2); */
        /* TF1*  mirror_tophatr = new TF1("mmirrortophatr", "[0]*x*x+[1]*x+[2]",0.03,0.07); */
        /* mirror_tophatr->SetParameter(0,zero); */
        /* mirror_tophatr->SetParameter(1,-one); */
        /* mirror_tophatr->SetParameter(2,two); */
        /* mirror_tophatr->Draw("SAME"); */
        c1->SaveAs("Effx_graph.pdf");
        for (int i=1; i<h1->GetNbinsX(); i++) {
            efficiency_scat_x->GetPoint(i,x,y);
            double ey = efficiency_scat_x->GetErrorYlow(i);
            /* std::cout << "x " << x << std::endl; */
            /* std::cout << "y " << y << std::endl; */
            /* std::cout << "ey " << ey << std::endl; */
            if  (y!=0 && h1->GetBinContent(i)!=0){
                /* std::cout << "sqrt(pow(h1->GetBinContent(i)/(y-ey)-h1->GetBinContent(i)/y,2)+pow(h1->GetBinError(i),2)) " << sqrt(pow(h1->GetBinContent(i)/(y-ey)-h1->GetBinContent(i)/y,2)+pow(h1->GetBinError(i),2)) << std::endl; */
                /* std::cout << "h1->GetBinContent(i)/(y-ey)-h1->GetBinContent(i)/y " << h1->GetBinContent(i)/(y-ey)-h1->GetBinContent(i)/y << std::endl; */
                /* std::cout << "h1->GetBinContent(i)/(y-ey) " << h1->GetBinContent(i)/(y-ey) << std::endl; */
                /* std::cout << "h1->GetBinContent(i)/y " << h1->GetBinContent(i)/y << std::endl; */
                /* std::cout << "h1->GetBinError(i) " << h1->GetBinError(i) << std::endl; */
                /* std::cout << "h1->GetBinContent(i) " << h1->GetBinContent(i) << std::endl; */
                h1->SetBinError(i,sqrt(pow(h1->GetBinContent(i)/(y-ey)-h1->GetBinContent(i)/y,2)+pow(h1->GetBinError(i),2)));
                if (muldiv==0){
                    h1->SetBinContent(i,h1->GetBinContent(i)/mirror_tophatl->Eval(x));
                }
                else {
                    h1->SetBinContent(i,h1->GetBinContent(i)*mirror_tophatl->Eval(x));
                }
                    /* std::cout << "mirror_tophatl->Eval(x) " << mirror_tophatl->Eval(x) << std::endl; */
                /* } */
                /* else if (x>0.03){ */
                    /* std::cout << "mirror_tophatr->Eval(x) " << mirror_tophatr->Eval(x) << std::endl; */
                    /* h1->SetBinContent(i,h1->GetBinContent(i)/mirror_tophatr->Eval(x)); */
                /* } */
                /* else { */
                    /* h1->SetBinContent(i,h1->GetBinContent(i)/y); */
                /* } */
                /* h1->SetBinError(i,sqrt(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)*h1->GetBinContent(i)+ey*ey/y*y)); */
            }
        }
        /* outfile->Close(); */
        Double_t* after = h1->GetIntegral();
        float change = (*before-*after)/2;
        h1->SetBinContent(0,h1->GetBinContent(0)-change);
        h1->SetBinContent(h1->GetNbinsX()+1,h1->GetBinContent(h1->GetNbinsX()+1)-change);
        f->Close();
        return h1;
    }

    TH1D* MCSAnalysis::trkreffiy(TH1D* h1, bool muldiv){

        TFile *f=new TFile(trkreffiname.c_str());
        TGraphAsymmErrors* efficiency_scat_y = (TGraphAsymmErrors*)f->Get("Effy_graph");
        TGraphAsymmErrors* mirror_efficiency_scat_y = (TGraphAsymmErrors*)f->Get("Effy_graph");
        outfile->cd();
        efficiency_scat_y->Write();
        double x;
        double y;
        double t = 2;
        for (int i=25; i<h1->GetNbinsX(); i++) {
            efficiency_scat_y->GetPoint(i,x,y);
            mirror_efficiency_scat_y->SetPoint(i,-x,y);
            t += 1;
        }
        TCanvas* c1 = new TCanvas();
        c1->SetBottomMargin(0.15);
        TF1* mirror_tophatl = new TF1("mmirrortophatl", "[1]*pow(x,6)+[2]*pow(x,4)+[3]*pow(x,2)+[4]",-0.07,0.07);
        mirror_tophatl->SetParameters(-2.70626e+05,-1.25699e+02,-4.02662e+04,2.63579e+00,9.94541e-01);
        efficiency_scat_y->Fit("mmirrortophatl", "R","",-0.07,0.07);
        /* TF1*  mirror_tophatl = new TF1("mmirrortophatl", "[0]*x*x+[1]*x+[2]",-0.07,-0.03); */
        /* mirror_efficiency_scat_y->Fit("mmirrortophatl", "R","",-0.07,-0.03); */
        /* mirror_efficiency_scat_y->Draw(); */
        /* c1->SaveAs("mirror_Effy_graph.pdf"); */
        /* c1->Clear(); */
        efficiency_scat_y->GetYaxis()->SetTitle("Acceptance");
        efficiency_scat_y->GetYaxis()->SetTitleSize(0.06);
        efficiency_scat_y->GetXaxis()->SetTitle("#it{#theta_{Y}} (mrad)");
        efficiency_scat_y->GetXaxis()->SetTitleSize(0.06);
        efficiency_scat_y->GetXaxis()->SetRangeUser(-0.0465,0.0465);
        efficiency_scat_y->Draw("AP");
        efficiency_scat_y->Draw("AP");
        mirror_tophatl->Draw("SAME");
    TText t1 = TText(0.18,0.305,"MICE Preliminary");
    TText t3 = TText(0.18,0.27,"ISIS cycle 2015/04");
    TText t2 = TText(0.18,0.235,"LiH, MAUS v3.3.2");
    t1.SetNDC(1);
    t1.SetTextSize(0.04);
    t1.SetTextFont(42);
    t2.SetNDC(1);
    t2.SetTextSize(0.03);
    t2.SetTextFont(42);
    t3.SetNDC(1);
    t3.SetTextSize(0.04);
    t3.SetTextFont(42);
            t1.Draw("SAME");
            t2.Draw("SAME");
            t3.Draw("SAME");
        /* double zero = mirror_tophatl->GetParameter(0); */
        /* double one = mirror_tophatl->GetParameter(1); */
        /* double two = mirror_tophatl->GetParameter(2); */
        /* TF1*  mirror_tophatr = new TF1("mmirrortophatr", "[0]*x*x+[1]*x+[2]",0.03,0.07); */
        /* mirror_tophatr->SetParameter(0,zero); */
        /* mirror_tophatr->SetParameter(1,-one); */
        /* mirror_tophatr->SetParameter(2,two); */
        /* mirror_tophatr->Draw("SAME"); */
        c1->SaveAs("Effy_graph.pdf");
        h1->Sumw2();
        Double_t* before = h1->GetIntegral();
        for (int i=1; i<h1->GetNbinsX(); i++) {
            efficiency_scat_y->GetPoint(i,x,y);
            double ey = efficiency_scat_y->GetErrorYlow(i);
            if  (y!=0 && h1->GetBinContent(i)!=0){
                h1->SetBinError(i,sqrt(pow(h1->GetBinContent(i)/(y-ey)-h1->GetBinContent(i)/y,2)+pow(h1->GetBinError(i),2)));
                /* if (x<-0.03){ */
                if (muldiv==0){
                    h1->SetBinContent(i,h1->GetBinContent(i)/mirror_tophatl->Eval(x));
                }
                else {
                    h1->SetBinContent(i,h1->GetBinContent(i)*mirror_tophatl->Eval(x));
                }
                /* } */
                /* else if (x>0.03){ */
                /*     h1->SetBinContent(i,h1->GetBinContent(i)/mirror_tophatr->Eval(x)); */
                /* } */
                /* else { */
                /*     h1->SetBinContent(i,h1->GetBinContent(i)/y); */
                /* } */
            }
        }
        Double_t* after = h1->GetIntegral();
        float change = (before-after)/2;
        h1->SetBinContent(0,h1->GetBinContent(0)-change);
        h1->SetBinContent(h1->GetNbinsX()+1,h1->GetBinContent(h1->GetNbinsX()+1)-change);
        f->Close();

        return h1;
    }

    TH1D* MCSAnalysis::trkreffiscatt(TH1D* h1){

        TFile *f=new TFile(trkreffiname.c_str());
        TGraphAsymmErrors* efficiency_scat_scatt = (TGraphAsymmErrors*)f->Get("Effscatt_graph");
        double x;
        double y;
        TCanvas* c1 = new TCanvas();
        h1->Sumw2();
        std::cout << "9" << std::endl; 
        for (int i=1; i<h1->GetNbinsX(); i++) {
            efficiency_scat_scatt->GetPoint(i,x,y);
            double ey = efficiency_scat_scatt->GetErrorYlow(i);
            if (y!=0 && h1->GetBinContent(i)!=0){
                h1->SetBinContent(i,h1->GetBinContent(i)/y);
                h1->SetBinError(i,sqrt(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)*h1->GetBinContent(i)+ey*ey/y*y));
            }
            /*
               else {
               h1->SetBinContent(i,h1->GetBinContent(i));
               }
               */
        }
        f->Close();

        return h1;
    }

    TH1D* MCSAnalysis::trkreffi2scatt(TH1D* h1){

        TFile *f=new TFile(trkreffiname.c_str());
        TGraphAsymmErrors* efficiency_scat_2scatt = (TGraphAsymmErrors*)f->Get("Eff2scatt_graph");
        double x;
        double y;
        TCanvas* c1 = new TCanvas();
        h1->Sumw2();
        for (int i=1; i<h1->GetNbinsX(); i++) {
            efficiency_scat_2scatt->GetPoint(i,x,y);
            double ey = efficiency_scat_2scatt->GetErrorYlow(i);
            if (y!=0 && h1->GetBinContent(i)!=0){
                h1->SetBinContent(i,h1->GetBinContent(i)/y);
                h1->SetBinError(i,sqrt(h1->GetBinError(i)*h1->GetBinError(i)/h1->GetBinContent(i)*h1->GetBinContent(i)+ey*ey/y*y));
            }
            /*
               else {
               h1->SetBinContent(i,h1->GetBinContent(i));
               }
               */
        }
        f->Close();

        return h1;
    }

    double MCSAnalysis::MomentumFromTOF(const bool& isdata=true){
        // Cuts remove events where the following statements do not make sense so we proceed without cuts.
        double rawTOF0HitTime = -1., rawTOF1HitTime = -1., rawTOF2HitTime = -1.;
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) == 1)
            rawTOF0HitTime  = 
                tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) == 1)
            rawTOF1HitTime  = 
                tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) == 1)
            rawTOF2HitTime = 
                tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime();
        // tofevent->GetTOFEventSlabHit().GetTOF2SlabHitArray()[1].GetRawTime() ) / 2.;
        // need to hard code a few things here (unfortunately) pertaining to the geometry and the time of flight.
        double dt0 = (_sys["TOF1_z"] - _sys["TOF0_z"]) / 0.299792458 / 1000.; // ns. 
        double dt  = 10.0 * dt0; // something rediculously large as a default number
        double pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
        if ( rawTOF1HitTime != -1 && rawTOF0HitTime != -1 ){
            dt  = rawTOF1HitTime - rawTOF0HitTime; 
            double dt0 = (_sys["TOF1_z"] - _sys["TOF0_z"]) / 0.299792458 / 1000.;
            double pz = 105.65/sqrt(pow(dt,2)/pow(dt0,2)-1);
            /* pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0) - 36.1; */
        }
        if ( rawTOF1HitTime != -1 && rawTOF2HitTime != -1 ){
            // Better estimate of the longitudinal momentum
            dt  = rawTOF2HitTime - rawTOF1HitTime; 
            double pz1 = pz;
            /* pz  = 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0); */
            dt0 = (_sys["TOF2_z"] - _sys["TOF1_z"]) / 0.299792458 / 1000.;
            pz = 105.65/sqrt(pow(dt,2)/pow(dt0,2)-1);
            //TOFcom->Fill(pz1,pz);
        }
        if(isdata) uncor_mom->Fill(pz);
        //else
        //   mccalc_mom->Fill(pz);
        return pz;
    }

    void MCSAnalysis::ConvolveWithInputDistribution(std::string distname){
        int isfirst = 0;
        bool isGEANT = false;
        bool isCobb = false;
        bool isMoliere = false;
        if (distname.find(modelname1.c_str()) != std::string::npos)
            isGEANT = true;
        if (distname.find(modelname2.c_str()) != std::string::npos)
            isCobb = true;
        if (distname.find(modelname3.c_str()) != std::string::npos)
            isMoliere = true;

        TFile* infile = new TFile(modelfile.c_str());

        TH1D* hiswX = new TH1D("hiswX","", 1000, -5, 5);
        TH1D* hiswY = new TH1D("hiswY","", 1000, -5, 5);

        // Efficiency plots
        TH1D* scatx = 
            new TH1D("scatx","Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        TH1D* scaty = 
            new TH1D("scaty","Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        TH1D* scatscat = 
            new TH1D("scatscat","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",
                    _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

        TH1D* scat2scatt = 
            new TH1D("scat2scatt","Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

        std::string tmpname = "thetaX_refconv_";
        tmpname += distname;
        TH1D* thetaX_refconv = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_holdrefconv_";
        tmpname += distname;
        TH1D* thetaX_holdrefconv = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_empty_";
        tmpname += distname;
        TH1D* thetaX_empty = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_emptyminus_";
        tmpname += distname;
        TH1D* thetaX_emptyminus = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_emptyasymm_";
        tmpname += distname;
        TH1* thetaX_emptyasymm = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);

        tmpname = "thetaY_refconv_";
        tmpname += distname;
        TH1D* thetaY_refconv = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_holdrefconv_";
        tmpname += distname;
        TH1D* thetaY_holdrefconv = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_empty_";
        tmpname += distname;
        TH1D* thetaY_empty = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname += distname;
        TH1D* thetaY_emptyfold = 
            new TH1D(tmpname.c_str(),";#it{#theta_{Y}};#it{#theta_{Yi}}-#it{#theta_{Y46-i}}",
                    _histlimits["NbinsXY"]/2, _histlimits["minXY"], 0);
        tmpname = "thetaY_emptyminus_";
        tmpname += distname;
        TH1D* thetaY_emptyminus = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_emptyasymm_";
        tmpname += distname;
        TH1* thetaY_emptyasymm = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);

        tmpname = "thetaScatt_refconv_";
        tmpname += distname;
        TH1D* thetaScatt_refconv = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",
                    _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

        tmpname = "thetaScatt_empty_";
        tmpname += distname;
        TH1D* thetaScatt_empty = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",
                    _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        tmpname = "thetaScatt_holdrefconv_";
        tmpname += distname;
        TH1D* thetaScatt_holdrefconv = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",
                    _histlimits["NbinsTh"]*4, _histlimits["minTh"]*2, _histlimits["maxTh"]*4);


        tmpname = "theta2Scatt_refconv_";
        tmpname += distname;
        TH1D* theta2Scatt_refconv = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",
                    _histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);


        tmpname = "thetaScatt_refconv_vp_";
        tmpname += distname;
        TH2D* thetaScatt_refconv_vp = 
            new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #it{#theta_{Scatt}}", 
                    200, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

        TH2D* thetaXUS_thetaXDS = 
            new TH2D("thetaXUS_thetaXDS","Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], 
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        TH2D* thetaYUS_thetaYDS = 
            new TH2D("thetaYUS_thetaYDS","Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "thetaX_";
        tmpname += distname;

        std::cout<<"Convolution with "<<tmpname<<" from "<<modelfile<<std::endl;
        TH1D* hx = (TH1D*)infile->Get(tmpname.c_str());
        tmpname = "thetaY_";
        tmpname += distname;
        TH1D* hy = (TH1D*)infile->Get(tmpname.c_str());
        tmpname = "thetaScatt_";
        tmpname += distname;
        TH1D* hScatt = (TH1D*)infile->Get(tmpname.c_str());
        tmpname = "theta2Scatt_";
        tmpname += distname;
        TH1D* h2Scatt = new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad", _histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
        for (int j=0; j<1e6; j++){
            double angle = hScatt->GetRandom();
            double anglesqu = angle*angle;
            h2Scatt->Fill(anglesqu);
        }



        // Collection DSConvSet;
        // for (int l=0; l<10; l++){
        for (int i=0; i<_USMCset.N(); i++){
            std::vector<double> nominalTheta = RotDefineProjectionAngles(_USMCset.E(i), _DSMCset.E(i), angdef);
            double nomthetaX = nominalTheta[0];
            double nomthetaY = nominalTheta[1];
            double nomthetascat = nominalTheta[2];
            thetaY_empty->Fill(nomthetaY);
            thetaY_emptyminus->Fill(-nomthetaY);
            thetaX_empty->Fill(nomthetaX);
            thetaX_emptyminus->Fill(-nomthetaX);
            thetaScatt_empty->Fill(nomthetascat);

            for (int j=0; j<20; j++){
                double dthetaX = -1;
                double dthetaY = -1;
                dthetaX = hx->GetRandom() * _sys["resX"];
                dthetaY = hy->GetRandom() * _sys["resY"];
                /*
                   while (dthetaX<accpt || dthetaX>-accpt || dthetaY<accpt || dthetaY>-accpt){
                //std::cout << "accpt " << accpt << std::endl;
                //std::cout << "dthetaX " << dthetaX << std::endl;
                //std::cout << "dthetaY " << dthetaY << std::endl;
                dthetaX = hx->GetRandom() * _sys["resX"];
                dthetaY = hy->GetRandom() * _sys["resY"];
                }
                */
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
                tmpvar.isgood = _USMCset.E(i).isgood;
                tmpvar.pid = _USMCset.E(i).pid;
                tmpvar.projX = _USMCset.E(i).projX;
                tmpvar.projY = _USMCset.E(i).projY;
                tmpvar.projdifX = _USMCset.E(i).projdifX;
                tmpvar.projdifY = _USMCset.E(i).projdifY;
                tmpvar.chi2nUS = _USMCset.E(i).chi2nUS;
                tmpvar.chi2nDS = _USMCset.E(i).chi2nDS;
                tmpvar.beamtype = _USMCset.E(i).beamtype;
                tmpvar.jUS = _USMCset.E(i).jUS;
                tmpvar.kUS = _USMCset.E(i).kUS;
                tmpvar.lUS = _USMCset.E(i).lUS;
                tmpvar.jDS = _USMCset.E(i).jDS;
                tmpvar.kDS = _USMCset.E(i).kDS;
                tmpvar.lDS = _USMCset.E(i).lDS;
                tmpvar.difcut = _USMCset.E(i).difcut;
                tmpvar.fidcut = _USMCset.E(i).fidcut;
                tmpvar.chicut = _USMCset.E(i).chicut;
                tmpvar.TOFcut = _USMCset.E(i).TOFcut;
                tmpvar.ptruth = _USMCset.E(i).ptruth;
                tmpvar.mcX = _USMCset.E(i).mcX;
                tmpvar.mcY = _USMCset.E(i).mcY;
                tmpvar.mcZ = _USMCset.E(i).mcZ;
                tmpvar.mcdXdz = _USMCset.E(i).mcdXdz;
                tmpvar.mcdYdz = _USMCset.E(i).mcdYdz;
                tmpvar.mcpx = _USMCset.E(i).mcpx;
                tmpvar.mcpy = _USMCset.E(i).mcpy;
                std::vector<double> projDTheta = RotDefineProjectionAngles(tmpvar, _DSMCset.E(i), angdef);
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
                std::vector<double> projTheta = RotDefineProjectionAngles(_USMCset.E(i), tmpvar, angdef);
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
                /* thetaX_refconv->Fill(thetaX); */
                /* thetaY_refconv->Fill(thetaY); */
                //thetaScatt_refconv->Fill(thetaScatt);
                //theta2Scatt_refconv->Fill(thetaScatt*thetaScatt);
                thetaScatt_refconv_vp->Fill(_USMCset.E(i).pz, thetaScatt);

                // Apply efficiency correction via weighting 
                /*
                   double binx = scatx->Fill(thetaX);
                //std::cout << "binx " << binx << std::endl;
                TFile *f=new TFile(trkreffiname.c_str());
                TGraphAsymmErrors* efficiency_scat_x = (TGraphAsymmErrors*)f->Get("Effx_graph");
                double x;
                double y;
                efficiency_scat_x->GetPoint(binx,x,y);
                double effipred = y; 
                binx = scatx->Fill(nomthetaX);
                TFile *fempty=new TFile(trkreffiemptyname.c_str());
                TGraphAsymmErrors* efficiency_empty_scat_x = (TGraphAsymmErrors*)fempty->Get("Effx_graph");
                double a;
                double b;
                efficiency_empty_scat_x->GetPoint(binx,a,b);
                double effiempdata = b;
                if (effiempdata==0) effiempdata = 0.1;
                //std::cout << "effipred " << effipred << std::endl;
                //std::cout << "binx " << binx << std::endl;
                //std::cout << "effiempdata " << effiempdata << std::endl;
                double weightX = effipred/effiempdata; 
                //std::cout << "weightX " << weightX << std::endl;
                double biny = scaty->Fill(thetaY);
                //std::cout << "biny " << biny << std::endl;
                TGraphAsymmErrors* efficiency_scat_y = (TGraphAsymmErrors*)f->Get("Effy_graph");
                double xy;
                double yy;
                efficiency_scat_y->GetPoint(biny,xy,yy);
                double effipredy = yy; 
                biny = scaty->Fill(nomthetaY);
                TGraphAsymmErrors* efficiency_empty_scat_y = (TGraphAsymmErrors*)fempty->Get("Effy_graph");
                double ay;
                double by;
                efficiency_empty_scat_y->GetPoint(biny,ay,by);
                double effiempdatay = by;
                if (effiempdatay==0) effiempdatay = 0.1;
                double weightY = effipredy/effiempdatay; 
                //std::cout << "effipredy " << effipredy << std::endl;
                //std::cout << "biny " << biny << std::endl;
                //std::cout << "effiempdatay " << effiempdatay << std::endl;
                //std::cout << "weightY " << weightY << std::endl;
                double binscatt = scatscat->Fill(thetaScatt);
                TGraphAsymmErrors* efficiency_scat_scatt = (TGraphAsymmErrors*)f->Get("Effscatt_graph");
                double xscatt;
                double yscatt;
                efficiency_scat_scatt->GetPoint(binscatt,xscatt,yscatt);
                double effipredscatt = yscatt; 
                binscatt = scatscat->Fill(nomthetascat);
                TGraphAsymmErrors* efficiency_empty_scat_scatt = (TGraphAsymmErrors*)fempty->Get("Effscatt_graph");
                double ascatt;
                double bscatt;
                efficiency_empty_scat_scatt->GetPoint(binscatt,ascatt,bscatt);
                double effiempdatascatt = bscatt;
                if (effiempdatascatt==0) effiempdatascatt = 0.1;
                double weightscatt = effipredscatt/effiempdatascatt; 
                //std::cout << "binscatt " << binscatt << std::endl;
                //std::cout << "effipredscatt " << effipredscatt << std::endl;
                //std::cout << "effiempdatascatt " << effiempdatascatt << std::endl;
                //std::cout << "weightscatt " << weightscatt << std::endl;
                double bin2scatt = scat2scatt->Fill(thetaScatt*thetaScatt);
                TGraphAsymmErrors* efficiency_scat_2scatt = (TGraphAsymmErrors*)f->Get("Eff2scatt_graph");
                double x2scatt;
                double y2scatt;
                efficiency_scat_2scatt->GetPoint(bin2scatt,x2scatt,y2scatt);
                double effipred2scatt = y2scatt; 
                bin2scatt = scat2scatt->Fill(nomthetascat*nomthetascat);
                TGraphAsymmErrors* efficiency_empty_scat_2scatt = (TGraphAsymmErrors*)fempty->Get("Eff2scatt_graph");
                double a2scatt;
                double b2scatt;
                efficiency_empty_scat_2scatt->GetPoint(bin2scatt,a2scatt,b2scatt);
                double effiempdata2scatt = b2scatt;
                if (effiempdata2scatt==0) effiempdata2scatt = 0.1;
                if (effipred2scatt==0) effipred2scatt = 1;
                double weight2scatt = effipred2scatt/effiempdata2scatt; 
                f->Close();
                fempty->Close();
                //std::cout << "bin2scatt " << bin2scatt << std::endl;
                //std::cout << "effipred2scatt " << effipred2scatt << std::endl;
                //std::cout << "effiempdata2scatt " << effiempdata2scatt << std::endl;
                //std::cout << "weight2scatt " << weight2scatt << std::endl;
                //std::cout << "weightX " << weightX << std::endl;
                //std::cout << "weightY " << weightY << std::endl;
                //std::cout << "weightscatt " << weightscatt << std::endl;
                //std::cout << "weight2scatt " << weight2scatt << std::endl;

                hiswX->Fill(weightX);
                hiswY->Fill(weightY);
                if (isGEANT) {
                    resp_thetaX.Fill(thetaX, d_thetaX, weightX);
                    resp_thetaY.Fill(thetaY, d_thetaY, weightY);
                    resp_thetaScatt.Fill(thetaScatt, dthetaScatt, weightscatt);
                    resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt, weight2scatt);
                }
                if (isCobb) {
                    tresp_thetaX.Fill(thetaX, d_thetaX, weightX);
                    tresp_thetaY.Fill(thetaY, d_thetaY, weightY);
                    tresp_thetaScatt.Fill(thetaScatt, dthetaScatt, weightscatt);
                    tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt, weight2scatt);
                }
                if (isMoliere) {
                    mresp_thetaX.Fill(thetaX, d_thetaX, weightX);
                    mresp_thetaY.Fill(thetaY, d_thetaY, weightY);
                    mresp_thetaScatt.Fill(thetaScatt, dthetaScatt, weightscatt);
                    mresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt, weight2scatt);
                }
                */
                    if (isGEANT) {
                        resp_thetaX.Fill(thetaX, -d_thetaX);
                        resp_thetaY.Fill(thetaY, -d_thetaY);
                        resp_thetaScatt.Fill(thetaScatt, dthetaScatt);
                        resp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
                    }
                if (isCobb) {
                    tresp_thetaX.Fill(thetaX, -d_thetaX);
                    tresp_thetaY.Fill(thetaY, -d_thetaY);
                    tresp_thetaScatt.Fill(thetaScatt, dthetaScatt);
                    tresp_theta2Scatt.Fill(thetaScatt*thetaScatt, dthetaScatt*dthetaScatt);
                }
                if (isMoliere) {
                    mresp_thetaX.Fill(thetaX, -d_thetaX);
                    mresp_thetaY.Fill(thetaY, -d_thetaY);
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

        thetaX_empty = trkreffix(thetaX_empty, 0);
        thetaY_empty = trkreffiy(thetaY_empty, 0);
        /* thetaScatt_empty = trkreffiscatt(thetaScatt_empty); */
        for (int i=0;i<1e6;i++) {
            Double_t x2 = thetaX_empty->GetRandom();
            x2 += hx->GetRandom();
            Double_t y2 = thetaY_empty->GetRandom();
            y2 += hy->GetRandom();
            Double_t scat2 = thetaScatt_empty->GetRandom();
            scat2 += hScatt->GetRandom();
            thetaX_holdrefconv->Fill(x2);
            thetaY_holdrefconv->Fill(y2);
            thetaScatt_holdrefconv->Fill(scat2);
            x2=0;
            y2=0;
            scat2=0;
        }

        int meanpx = 0;
        int i=0;
        while (meanpx==0) {
            if (thetaX_holdrefconv->GetBinLowEdge(i)>thetaX_holdrefconv->GetMean()) meanpx=i;
            i++;
        }
        int meanpy = 0;
        i=0;
        while (meanpy==0) {
            if (thetaY_holdrefconv->GetBinLowEdge(i)>thetaY_holdrefconv->GetMean()) meanpy=i;
            i++;
        }
        for (i = 1; i < 48; i++) thetaX_refconv->SetBinContent(i,thetaX_holdrefconv->GetBinContent(i+meanpx-25));
        for (i = 1; i < 48; i++) thetaY_refconv->SetBinContent(i,thetaY_holdrefconv->GetBinContent(i+meanpy-25));
        thetaX_refconv = trkreffix(thetaX_refconv,1);
        thetaY_refconv = trkreffiy(thetaY_refconv,1);
        int startscatt = 0;
        while (startscatt==0) {
            if (thetaScatt_holdrefconv->GetBinContent(i)>0) startscatt=i;
            i++;
        }
        for (i = 1; i < 46; i++) thetaScatt_refconv->SetBinContent(i,thetaScatt_holdrefconv->GetBinContent(i+startscatt));
        for (int j=0; j<1e6; j++){
            double angle = thetaScatt_refconv->GetRandom();
            double anglesqu = angle*angle;
            theta2Scatt_refconv->Fill(anglesqu);
        }



        TVectorD NoEntriesrefconv(1);
        NoEntriesrefconv[0] = thetaX_holdrefconv->GetEntries();
        TCanvas *c1 = new TCanvas();
        hiswX->Draw();
        c1->SaveAs("hiswX.pdf");
        c1->Clear();
        hiswY->Draw();
        c1->SaveAs("hiswY.pdf");
        outfile->cd();
        NoEntriesrefconv.Write("conv");
        thetaXUS_thetaXDS->Write();
        thetaYUS_thetaYDS->Write();
        thetaX_holdrefconv->Write();
        thetaX_refconv->Write();
        thetaY_refconv->Write();
        thetaY_empty->Write();
        for (int i=1;i<24;i++) {
            thetaY_emptyfold->SetBinContent(i,thetaY_empty->GetBinContent(i)-thetaY_empty->GetBinContent(46-i));
        }
        thetaY_emptyfold->Write();
        c1->Clear();
        thetaY_emptyasymm=thetaY_empty->GetAsymmetry(thetaY_emptyminus);
        thetaY_emptyasymm->Draw();
        TF1*  fM2D = new TF1("fM2D", "[0]+x*[1]", 0, 0.02);
        //fM2D->SetParameter(0,0);
        //fM2D->SetParameter(1,1);
        //fM2D->SetParNames("A","B");
        //thetaY_emptyasymm->Fit("pol1","RES","",0.,0.02);
        thetaY_emptyasymm->Fit(fM2D,"RES");
        //fM2D->Draw("SAME");
        thetaY_emptyasymm->SetName("thetaY_emptyasymm");
        thetaY_emptyasymm->Write();
        thetaY_emptyasymm->GetXaxis()->SetRange(0,0.06);
        c1->SaveAs("thetaY_emptyasymm.pdf");
        c1->Clear();
        thetaX_emptyasymm=thetaX_empty->GetAsymmetry(thetaX_emptyminus);
        thetaX_emptyasymm->Draw();
        //thetaX_emptyasymm->Fit("pol1","RES","",0,0.02);
        thetaX_emptyasymm->Fit(fM2D,"RES");
        //fM2D->Draw("SAME");
        thetaX_emptyasymm->SetName("thetaX_emptyasymm");
        thetaX_emptyasymm->Write();
        thetaX_emptyasymm->GetXaxis()->SetRange(0,0.06);
        c1->SaveAs("thetaX_emptyasymm.pdf");
        c1->Clear();
        thetaY_emptyfold->Draw();
        c1->Print("thetaY_emptyfold.pdf");
        c1->Print("thetaY_emptyfold.root");
        thetaX_empty->Write();
        thetaScatt_refconv->Write();
        theta2Scatt_refconv->Write();
        thetaScatt_refconv_vp->Write();
        hx->Write();
        hy->Write();
                    if (isGEANT) {
    TVectorD NoEntriesGEANT(1);
    NoEntriesGEANT[0] = hx->GetEntries();
    NoEntriesGEANT.Write("GEANT");
                    }
        hScatt->Write();
        h2Scatt->Write();
        delete c1;

    }

    void MCSAnalysis::ConvolveWithVirtualInputDistribution(std::string distname){
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
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "thetaY_refconv_";
        tmpname += distname;
        TH1D* thetaY_refconv = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "thetaScatt_refconv_";
        tmpname += distname;
        TH1D* thetaScatt_refconv = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",
                    _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);


        tmpname = "theta2Scatt_refconv_";
        tmpname += distname;
        TH1D* theta2Scatt_refconv = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",
                    _histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

        tmpname = "thetaScatt_refconv_vp_";
        tmpname += distname;
        TH2D* thetaScatt_refconv_vp = 
            new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #it{#theta_{Scatt}}", 
                    200, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

        TH2D* thetaXUS_thetaXDS = 
            new TH2D("thetaXUS_thetaXDS","Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], 
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        TH2D* thetaYUS_thetaYDS = 
            new TH2D("thetaYUS_thetaYDS","Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "thetaX_";
        tmpname += distname;

        std::cout<<"Convolution with "<<tmpname<<" from "<<modelfile<<std::endl;
        TH1D* hx = (TH1D*)infile->Get(tmpname.c_str());
        TH1D* hy = (TH1D*)infile->Get(tmpname.c_str());

        // Collection DSConvSet;
        // for (int l=0; l<10; l++){
        for (int i=0; i<USTruthSet.N(); i++){
            std::vector<double> nominalTheta = RotDefineProjectionAngles(USTruthSet.E(i), DSTruthSet.E(i), angdef);
            double nomthetaX = nominalTheta[0];
            double nomthetaY = nominalTheta[1];
            for (int j=0; j<20; j++){
                double dthetaX = hx->GetRandom() * _sys["resX"];
                double dthetaY = hy->GetRandom() * _sys["resY"];
                // First project the upstream track to the absorber 
                double zabspos = _sys["abspos"] + 0.0;
                Vars projvarAbs = PropagateVarsMu(USTruthSet.E(i), zabspos);
                double xabs = projvarAbs.X;  /// _USMCset.E(i).X + _USMCset.E(i).dXdz * dzabsUS;
                double yabs = projvarAbs.Y;  /// _USMCset.E(i).Y + _USMCset.E(i).dYdz * dzabsUS;
                // Now add the angle from the model to the downstream measurement.
                double dXdz_abs = DSTruthSet.E(i).dXdz + tan(dthetaY);
                double dYdz_abs = DSTruthSet.E(i).dYdz + tan(dthetaX);
                // double d_thetaY    = atan(dXdz_abs) - atan(_DSMCset.E(i).dXdz);
                // double d_thetaX    = atan(dYdz_abs) - atan(_DSMCset.E(i).dYdz);

                // Project the track into the downstream reference plane
                // double xref = xabs + dXdz_abs * dzabsDS;
                // double yref = yabs + dYdz_abs * dzabsDS;
                Vars tmpvar = DSTruthSet.E(i);
                tmpvar.X = xabs;
                tmpvar.Y = yabs;
                tmpvar.Z = _sys["abspos"];
                tmpvar.dXdz = dXdz_abs;
                tmpvar.dYdz = dYdz_abs;
                tmpvar.px   = dXdz_abs * USTruthSet.E(i).pz;
                tmpvar.py   = dYdz_abs * USTruthSet.E(i).pz;
                tmpvar.pz   = USTruthSet.E(i).pz;
                tmpvar.TOF12= USTruthSet.E(i).TOF12;
                tmpvar.TOF01= USTruthSet.E(i).TOF01;
                std::vector<double> projDTheta = RotDefineProjectionAngles(tmpvar, DSTruthSet.E(i), angdef);
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
                std::vector<double> projTheta = RotDefineProjectionAngles(USTruthSet.E(i), tmpvar, angdef);
                double thetaY = projTheta[1];   /// atan(tmpvar.dXdz) - atan(_USMCset.E(i).dXdz);
                double thetaX = projTheta[0];   /// atan(tmpvar.dYdz) - atan(_USMCset.E(i).dYdz);
                // double cosScatt = ( (1 + _USMCset.E(i).dXdz * tmpvar.dXdz +
                //		   _USMCset.E(i).dYdz * tmpvar.dYdz )/
                //		  sqrt(1 + _USMCset.E(i).dXdz*_USMCset.E(i).dXdz +
                //		       _USMCset.E(i).dYdz*_USMCset.E(i).dYdz)/
                //		  sqrt(1 + tmpvar.dXdz*tmpvar.dXdz +
                //		       tmpvar.dYdz*tmpvar.dYdz));
                double thetaScatt = projTheta[2];  /// acos(cosScatt);
                thetaXUS_thetaXDS->Fill(atan(USTruthSet.E(i).dXdz), atan(tmpvar.dXdz));
                thetaYUS_thetaYDS->Fill(atan(USTruthSet.E(i).dYdz), atan(tmpvar.dYdz));
                //thetaX_refconv->Fill(thetaX);
                //thetaY_refconv->Fill(thetaY);
                //thetaScatt_refconv->Fill(thetaScatt);
                //theta2Scatt_refconv->Fill(thetaScatt*thetaScatt);
                thetaScatt_refconv_vp->Fill(USTruthSet.E(i).pz, thetaScatt);

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
        //thetaX_refconv->Write();
        //thetaY_refconv->Write();
        //thetaScatt_refconv->Write();
        //theta2Scatt_refconv->Write();
        thetaScatt_refconv_vp->Write();

    }

    void MCSAnalysis::DoUnfolding(Collection& _USsethold, Collection& _DSsethold, Collection& _USMCsethold, Collection& _DSMCsethold){
        // sp, gap = gap after and before distribution for deconvolution
        int sp = 100;
        int sp2 = 100;
        int gap = 20;
        int gap2 = 20;

        std::string  tmpname = "thetaX_recoGold";
        TH1D* thetaX_recoGold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGold";
        TH1D* thetaY_recoGold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_recoGoldhold";
        TH1D* thetaX_recoGoldhold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGoldhold";
        TH1D* thetaY_recoGoldhold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_recoGoldhold2";
        TH1D* thetaX_recoGoldhold2 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGoldhold2";
        TH1D* thetaY_recoGoldhold2 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_recoGoldhold3";
        TH1D* thetaX_recoGoldhold3 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGoldhold3";
        TH1D* thetaY_recoGoldhold3 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaScatt_recoGold";
        TH1D* thetaScatt_recoGold = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        tmpname = "theta2Scatt_recoGold";
        TH1D* theta2Scatt_recoGold = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
        tmpname = "thetaX_response";
        TH1D* thetaX_response = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_response";
        TH1D* thetaY_response = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaScatt_response";
        TH1D* thetaScatt_response = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        tmpname = "theta2Scatt_response";
        TH1D* theta2Scatt_response = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
        tmpname = "theta2Scatt_recoGoldhold";
        TH1D* theta2Scatt_recoGoldhold = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",
                    sp2, _histlimits["minTh2"], 0.02);

        tmpname = "thetaX_recoGold_noeffi";
        TH1D* thetaX_recoGold_noeffi = 
            new TH1D("thetaX_recoGold_noeffi","Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGold_noeffi";
        TH1D* thetaY_recoGold_noeffi = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_recoGold_effi_only_gold";
        TH1D* thetaX_recoGold_effi_only_gold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGold_effi_only_gold";
        TH1D* thetaY_recoGold_effi_only_gold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_recoGold_effi_only_response";
        TH1D* thetaX_recoGold_effi_only_response = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGold_effi_only_response";
        TH1D* thetaY_recoGold_effi_only_response = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "thetaX_recoGoldhold_noeffi";
        TH1D* thetaX_recoGoldhold_noeffi = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGoldhold_noeffi";
        TH1D* thetaY_recoGoldhold_noeffi = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_recoGoldhold_effi_only_gold";
        TH1D* thetaX_recoGoldhold_effi_only_gold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGoldhold_effi_only_gold";
        TH1D* thetaY_recoGoldhold_effi_only_gold = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_recoGoldhold_effi_only_response";
        TH1D* thetaX_recoGoldhold_effi_only_response = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_recoGoldhold_effi_only_response";
        TH1D* thetaY_recoGoldhold_effi_only_response = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "thetaX_response_noeffi";
        TH1D* thetaX_response_noeffi = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_response_noeffi";
        TH1D* thetaY_response_noeffi = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        for(int i=0; i<_DSsethold.N(); i++){
            std::vector<double> projTheta = RotDefineProjectionAngles(_USsethold.E(i), _DSsethold.E(i), angdef);
            //double thetaY = atan(_DSsetfold.E(i).dXdz) - atan(_USsetfold.E(i).dXdz);
            //double thetaX = atan(_DSsetfold.E(i).dYdz) - atan(_USsetfold.E(i).dYdz);
            double thetaY = projTheta[1];
            double thetaX = projTheta[0];
            double thetaScatt = projTheta[2];

            //if (thetaX>0) thetaX_recoGold->Fill(thetaX);
            //if (thetaY>0) thetaY_recoGold->Fill(thetaY);
            thetaX_recoGold->Fill(thetaX);
            thetaY_recoGold->Fill(thetaY);
            thetaX_recoGold_noeffi->Fill(thetaX);
            thetaY_recoGold_noeffi->Fill(thetaY);
            thetaScatt_recoGold->Fill(thetaScatt);
            theta2Scatt_recoGold->Fill(thetaScatt*thetaScatt);
        }

        TVectorD NoEntriesGold(1);
        NoEntriesGold[0] = thetaX_recoGold->GetEntries();
        DSentries = thetaX_recoGold->Integral();

        /* this one */
        thetaX_recoGold = trkreffix(thetaX_recoGold,0);
        thetaY_recoGold = trkreffiy(thetaY_recoGold,0);
        thetaScatt_recoGold = trkreffiscatt(thetaScatt_recoGold);
        theta2Scatt_recoGold = trkreffi2scatt(theta2Scatt_recoGold);
        float* source_thetaX = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* source_thetaY = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* source_thetaScatt = new float[2*(sp+gap)+int(_histlimits["NbinsTh"])];
        float* source_theta2Scatt = new float[2*(sp2+gap2)+int(_histlimits["NbinsTh2"])];
        float* source_thetaX_noeffi = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* source_thetaY_noeffi = new float[sp+gap+int(_histlimits["NbinsXY"])];

        float* dataset_ori_lll = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* dataset_mod_lll = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* dataset_ori_jjj = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* dataset_mod_jjj = new float[sp+gap+int(_histlimits["NbinsXY"])];
        TMatrixDSym cov_calc_matrix(29);
        TMatrixDSym cov_calc_matrix_nonorm(29);
        TMatrixDSym cov_calc_matrix_binnorm(29);
        TMatrixDSym cov_calc_matrix_jjj(29);
        TMatrixD A(29,29);
        TMatrixD modA(29,29);
        TMatrixD myA(29,29);
        TMatrixD epsilon(29,29);

        for(int i=0;i<_DSMCsethold.N(); i++){
            std::vector<double> projTheta = RotDefineProjectionAngles(_USMCsethold.E(i), _DSMCsethold.E(i), angdef);
            //double thetaMCY = atan(_DSMCsetfold.E(i).dXdz) - atan(_USMCsetfold.E(i).dXdz);
            //double thetaMCX = atan(_DSMCsetfold.E(i).dYdz) - atan(_USMCsetfold.E(i).dYdz);
            double thetaMCY = projTheta[1]; 
            double thetaMCX = projTheta[0]; 
            double thetaScattMC = projTheta[2]; 
            //if (thetaMCX>0) thetaX_response->Fill(thetaMCX);
            //if (thetaMCY>0) thetaY_response->Fill(thetaMCY);
            thetaX_response->Fill(thetaMCX);
            thetaY_response->Fill(thetaMCY);
            thetaX_response_noeffi->Fill(thetaMCX);
            thetaY_response_noeffi->Fill(thetaMCY);
            thetaScatt_response->Fill(thetaScattMC);
            theta2Scatt_response->Fill(thetaScattMC*thetaScattMC);

        }

        TVectorD NoEntriesresp(1);
        NoEntriesresp[0] = thetaX_response->GetEntries();
        DSrefentries = thetaX_response->Integral();
        
        thetaX_response = trkreffix(thetaX_response,0);
        thetaY_response = trkreffiy(thetaY_response,0);
        thetaScatt_response = trkreffiscatt(thetaScatt_response);
        theta2Scatt_response = trkreffi2scatt(theta2Scatt_response);
        float* response_thetaX = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* response_thetaY = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* response_thetaScatt = new float[2*(sp+gap)+int(_histlimits["NbinsTh"])];
        float* response_theta2Scatt = new float[2*(sp2+gap2)+int(_histlimits["NbinsTh2"])];
        float* response_thetaX_noeffi = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* response_thetaY_noeffi = new float[sp+gap+int(_histlimits["NbinsXY"])];

        for(int i=0; i<sp+gap+int(_histlimits["NbinsXY"]); i++){
            source_thetaX[i] = 0;
            source_thetaY[i] = 0;
            response_thetaX[i] = 0;
            response_thetaY[i] = 0;
            source_thetaX_noeffi[i] = 0;
            source_thetaY_noeffi[i] = 0;
            response_thetaX_noeffi[i] = 0;
            response_thetaY_noeffi[i] = 0;
        }
        for(int i=0; i<2*(sp+gap)+_histlimits["NbinsTh"]; i++){
            source_thetaScatt[i] = 0;
            response_thetaScatt[i] = 0;
        }
        for(int i=0; i<2*(sp2+gap2)+_histlimits["NbinsTh2"]; i++){
            source_theta2Scatt[i] = 0;
            response_theta2Scatt[i] = 0;
        }
        for(int i=0; i<_histlimits["NbinsXY"]; i++){
            if (thetaX_recoGold->GetBinContent(i+1)<0){
                std::cout << "plot is neg" << std::endl;
            }
            source_thetaX[sp+i] = thetaX_recoGold->GetBinContent(i+1);
            source_thetaY[sp+i] = thetaY_recoGold->GetBinContent(i+1);

            response_thetaX[i] = thetaX_response->GetBinContent(i+1);
            response_thetaY[i] = thetaY_response->GetBinContent(i+1);
            if (thetaX_response->GetBinContent(i+1)<6){
                 response_thetaX[i] = 0;
            }
            if (thetaY_response->GetBinContent(i+1)<6){
                 response_thetaY[i] = 0;
            }
            source_thetaX_noeffi[sp+i] = thetaX_recoGold_noeffi->GetBinContent(i+1);
            source_thetaY_noeffi[sp+i] = thetaY_recoGold_noeffi->GetBinContent(i+1);
            response_thetaX_noeffi[i] = thetaX_response_noeffi->GetBinContent(i+1);
            response_thetaY_noeffi[i] = thetaY_response_noeffi->GetBinContent(i+1);
        }
        for(int i=0; i<_histlimits["NbinsTh"]; i++){
            source_thetaScatt[2*sp+i] = thetaScatt_recoGold->GetBinContent(i+1);
            response_thetaScatt[i] = thetaScatt_response->GetBinContent(i+1);
        }
        for(int i=0; i<_histlimits["NbinsTh2"]; i++){
            source_theta2Scatt[sp2+i+1] = theta2Scatt_recoGold->GetBinContent(i+1);
            source_theta2Scatt[sp2-i] = theta2Scatt_recoGold->GetBinContent(i+1);
            response_theta2Scatt[i+46+1] = theta2Scatt_response->GetBinContent(i+1);
            response_theta2Scatt[46-i] = theta2Scatt_response->GetBinContent(i+1);
        }

        // itMax = number of repetitions tested over
        int itMax = 100; 
        float *kolmogorovX = new float[itMax-1];
        float *kolmogorovY = new float[itMax-1];
        float *kolmogorovScatt = new float[itMax-1];
        float *kolmogorov2Scatt = new float[itMax-1];
        float *ChiSqX = new float[itMax-1];
        float *ChiSqY = new float[itMax-1]; 
        float *kurtX = new float[itMax-1];
        float *kurtY = new float[itMax-1]; 
        float *ChiSqScatt = new float[itMax-1];
        float *ChiSq2Scatt = new float[itMax-1]; 
        // hn* = n iterations distribution - hp* = n-1 iterations distributions 
        TH1D* hpX = new TH1D("hpX","hpX",_histlimits["NbinsXY"],_histlimits["minXY"],_histlimits["maxXY"]);
        TH1D* hnX = new TH1D("hnX","hnX",_histlimits["NbinsXY"],_histlimits["minXY"],_histlimits["maxXY"]);
        TH1D* hpY = new TH1D("hpY","hpY",_histlimits["NbinsXY"],_histlimits["minXY"],_histlimits["maxXY"]);
        TH1D* hnY = new TH1D("hnY","hnY",_histlimits["NbinsXY"],_histlimits["minXY"],_histlimits["maxXY"]);
        TH1D* hpScatt = new TH1D("hpScatt","hpScatt",_histlimits["NbinsTh"],_histlimits["minTh"],_histlimits["maxTh"]);
        TH1D* hnScatt = new TH1D("hnScatt","hnScatt",_histlimits["NbinsTh"],_histlimits["minTh"],_histlimits["maxTh"]);
        TH1D* hp2Scatt = new TH1D("hp2Scatt","hp2Scatt",_histlimits["NbinsTh2"],_histlimits["minTh2"],_histlimits["maxTh2"]);
        TH1D* hn2Scatt = new TH1D("hn2Scatt","hn2Scatt",_histlimits["NbinsTh2"],_histlimits["minTh2"],_histlimits["maxTh2"]);
        // source_* changed each iteration so Source_* is left unchanged
        float* diffbinX = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* diffbinY = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* Source_thetaX = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* Source_thetaY = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* Source_thetaX_noeffi = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* Source_thetaY_noeffi = new float[sp+gap+int(_histlimits["NbinsXY"])];
        float* Source_thetaScatt = new float[2*(sp+gap)+int(_histlimits["NbinsTh"])];
        float* Source_theta2Scatt = new float[2*(sp2+gap2)+int(_histlimits["NbinsTh2"])];

        for( int i=0; i<sp+gap+_histlimits["NbinsXY"]; i++){
            Source_thetaX[i] = source_thetaX[i];
            Source_thetaY[i] = source_thetaY[i];
            Source_thetaX_noeffi[i] = source_thetaX_noeffi[i];
            Source_thetaY_noeffi[i] = source_thetaY_noeffi[i];
        }

        TSpectrum *sX = new TSpectrum();
        TSpectrum *sY = new TSpectrum();
        TSpectrum *sX_noeffi = new TSpectrum();
        TSpectrum *sY_noeffi = new TSpectrum();
        TSpectrum *sX_effi_only_gold = new TSpectrum();
        TSpectrum *sY_effi_only_gold = new TSpectrum();
        TSpectrum *sX_effi_only_response = new TSpectrum();
        TSpectrum *sY_effi_only_response = new TSpectrum();

        for(int i=1; i<=itMax; i++){
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                source_thetaX[n] = Source_thetaX[n];
                source_thetaY[n] = Source_thetaY[n];
            }
            sX->DeconvolutionRL(source_thetaX, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), i, 20, 0);
            sY->DeconvolutionRL(source_thetaY, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), i, 20, 0);

            for(int n=0; n<_histlimits["NbinsXY"]; n++){
                hnX->SetBinContent(n+1,source_thetaX[sp+n]);
                hnY->SetBinContent(n+1,source_thetaY[sp+n]);
            }
            if (i!=1) {
                // test results added to arrays
                ChiSqX[i-2] = float(hnX->Chi2Test(hpX,"CHI2"));
                kolmogorovX[i-2] = float(hnX->KolmogorovTest(hpX,"M"));
                ChiSqY[i-2] = float(hnY->Chi2Test(hpY,"CHI2"));
                kolmogorovY[i-2] = float(hnY->KolmogorovTest(hpY,"M"));
                kurtX[i-2] = float(hnX->GetKurtosis());
                kurtY[i-2] = float(hnY->GetKurtosis());
            }
            hnX->Clear();
            hpX->Clear();
            hnY->Clear();
            hpY->Clear();
            for(int n=0; n<_histlimits["NbinsXY"]; n++){
                hpX->SetBinContent(n+1,source_thetaX[sp+n]);
                hpY->SetBinContent(n+1,source_thetaY[sp+n]);
            }
            sX->Clear();
            sY->Clear();
        }
        
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                source_thetaX[n] = Source_thetaX[n];
                source_thetaY[n] = Source_thetaY[n];
            }

            sX->DeconvolutionRL(source_thetaX, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            sY->DeconvolutionRL(source_thetaY, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                dataset_ori_lll[n] = source_thetaX[n];
            }
        for(int rrr=0; rrr<1000; rrr++){

            TRandom r(0);
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                source_thetaX[n] = Source_thetaX[n]+r.Uniform(-int(sqrt(Source_thetaX[n])),int(sqrt(Source_thetaX[n])));
                source_thetaY[n] = Source_thetaY[n];
            }
            sX->DeconvolutionRL(source_thetaX, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            sY->DeconvolutionRL(source_thetaY, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                dataset_mod_lll[n] = source_thetaX[n];
            }
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                diffbinX[n] += source_thetaX[n];
                diffbinY[n] += source_thetaY[n];
            }

        }
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                diffbinX[n] = dataset_ori_lll[n]-diffbinX[n]/1000;
                diffbinY[n] = dataset_ori_lll[n]-diffbinY[n]/1000;
            }
        
        for(int rrr=10; rrr<39; rrr++){
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                source_thetaX[n] = Source_thetaX[n];
                source_thetaY[n] = Source_thetaY[n];
                /* std::cout << "source_thetaX[rrr] " << source_thetaX[n] << std::endl; */
            }
            /* std::cout << "source_thetaX[sp+rrr] " << source_thetaX[sp+rrr] << std::endl; */

            /* std::cout << "source_thetaX[sp+rrr] " << source_thetaX[sp+rrr] << std::endl; */
            sX->DeconvolutionRL(source_thetaX, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            sY->DeconvolutionRL(source_thetaY, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            /* std::cout << "source_thetaX[sp+rrr] " << source_thetaX[sp+rrr] << std::endl; */
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                dataset_ori_lll[n] = source_thetaX[n];
                /* std::cout << "dataset_ori_lll[rrr] " << dataset_ori_lll[n] << std::endl; */
                /* std::cout << "source_thetaX[rrr] " << source_thetaX[n] << std::endl; */
                /* source_thetaX[n]=source_thetaX[n]*2; */
                /* std::cout << "dataset_ori_lll[rrr] " << dataset_ori_lll[n] << std::endl; */

            }
            /* std::cout << "dataset_ori_lll[sp+rrr] " << dataset_ori_lll[sp+rrr] << std::endl; */

            float dataset_ori_lll_total = 0;
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++) dataset_ori_lll_total += dataset_ori_lll[n];
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                source_thetaX[n] = Source_thetaX[n];
                source_thetaY[n] = Source_thetaY[n];
            }

            source_thetaX[sp+rrr] = Source_thetaX[sp+rrr]+round(sqrt(Source_thetaX[sp+rrr]));
            source_thetaY[sp+rrr] = Source_thetaY[sp+rrr]+round(sqrt(Source_thetaY[sp+rrr]));
            sX->DeconvolutionRL(source_thetaX, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            sY->DeconvolutionRL(source_thetaY, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                dataset_mod_lll[n] = source_thetaX[n];
            }
            /* std::cout << "dataset_mod_lll[sp+rrr] " << dataset_mod_lll[sp+rrr] << std::endl; */
            float dataset_mod_lll_total = 0;
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++) dataset_mod_lll_total += dataset_mod_lll[n];
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                source_thetaX[n] = Source_thetaX[n];
                source_thetaY[n] = Source_thetaY[n];
            }

            /* for(int i=0; i<int(_histlimits["NbinsXY"]+gap); i++){ */
            /*     if (::isnan(dataset_ori_lll[sp+i])==0 && ::isinf(dataset_ori_lll[sp+i])==0) { */
            /*         thetaX_covorihold->SetBinContent(i+1, dataset_ori_lll[sp+i]); */
            /*     } */
            /*     else { */ 
            /*         thetaX_covorihold->SetBinContent(i+1, 0); */
            /*     } */
            /*     if (::isnan(dataset_mod_lll[sp+i])==0 && ::isinf(dataset_mod_lll[sp+i])==0) { */
            /*         thetaX_covmodhold->SetBinContent(i+1, dataset_mod_lll[sp+i]); */
            /*     } */
            /*    else { */ 
            /*         thetaX_covmodhold->SetBinContent(i+1, 0); */

            /*     } */
            /* } */
            /* int meanpx = 0; */
            /* int i=0; */
            /* double meanx = double(thetaX_covorihold->GetMean()); */
            /* while (meanpx==0 && ::isnan(meanx)==0) { */
            /*     if (thetaX_recoGoldhold->GetBinLowEdge(i)>thetaX_covorihold->GetMean()) meanpx=i; */
            /*     i++; */
            /* } */
            /* int meanpy = 0; */
            /* i=0; */
            /* double meany = double(thetaY_covorihold->GetMean()); */
            /* while (meanpy==0 && ::isnan(meany)==0) { */
            /*     if (thetaY_recoGoldhold->GetBinLowEdge(i)>thetaY_covorihold->GetMean()) meanpy=i; */
            /*     i++; */
            /* } */
            // diagonal calc.
            /* std::cout << "dataset_ori_lll[sp+rrr] " << dataset_ori_lll[sp+rrr] << std::endl; */
            /* std::cout << "dataset_mod_lll[sp+rrr] " << dataset_mod_lll[sp+rrr] << std::endl; */
            float cov_calc_ii = 0;
            float cov_calc_jj = 0;
            float cov_calc_ii_ori = 0;
            float cov_calc_ii_bin = 0;
            float A_ij = 0;
            float modA_ij = 0;
            float myA_ij = 0;
            float epsilon_ij = 0;
            myA_ij = (dataset_mod_lll[sp+rrr] - dataset_ori_lll[sp+rrr])/(Source_thetaX[sp+rrr]*1.1);
            myA(rrr-10,rrr-10) = myA_ij;
            modA_ij = pow((dataset_mod_lll[sp+rrr] - dataset_ori_lll[sp+rrr]),2)/(Source_thetaX[sp+rrr]*1.1);
            std::cout << "modA_ij " << modA_ij << std::endl;
            modA(rrr-10,rrr-10) = modA_ij;
            A_ij = (dataset_mod_lll[sp+rrr] - dataset_ori_lll[sp+rrr])/(Source_thetaX[sp+rrr]*1.1-Source_thetaX[sp+rrr]);
            A(rrr-10,rrr-10) = A_ij;
            epsilon_ij = pow((Source_thetaX[sp+rrr]*1.1 - Source_thetaX[sp+rrr]),2)/(Source_thetaX[sp+rrr]*1.1);
            epsilon(rrr-10,rrr-10) = epsilon_ij;
            for( int n=sp+gap+10; n<sp+gap+39; n++ ){ 
                cov_calc_jj += pow((dataset_ori_lll[n] - dataset_mod_lll[n]),2);

            }
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++ ){ 
                cov_calc_ii += pow((dataset_ori_lll[n] - dataset_mod_lll[n]),2);
            }
            /* cov_calc_ii=cov_calc_ii/(dataset_ori_lll_total+dataset_mod_lll_total); */
            cov_calc_ii_ori=cov_calc_ii/(dataset_mod_lll_total);
            cov_calc_matrix(rrr-10,rrr-10) = cov_calc_ii_ori;
            cov_calc_ii=cov_calc_ii;
            cov_calc_matrix_nonorm(rrr-10,rrr-10) = cov_calc_ii;
            cov_calc_ii_bin=cov_calc_ii/29;
            cov_calc_matrix_binnorm(rrr-10,rrr-10) = cov_calc_ii_bin;
            
            cov_calc_jj=cov_calc_jj;
            cov_calc_matrix_jjj(rrr-10,rrr-10) = cov_calc_jj;
            /* std::cout < "rrr-10 " << rrr-10 << std::endl; */
            /* std::cout << "cov_calc_ii " << cov_calc_ii << std::endl; */
            /* std::cout << "rrr-10 " << rrr-10 << std::endl; */

            /* cov_calc_matrix(rrr-10,rrr-10) = dataset_mod_lll_total; */

            /* cov_calc_matrix.Print(); */

            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                dataset_ori_jjj[n] = dataset_ori_lll[n];
            }
            float dataset_ori_jjj_total = 0;
            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++) dataset_ori_jjj_total += dataset_ori_jjj[n];


            for(int j=rrr+1; j<39; j++){
                source_thetaX[sp+j] = Source_thetaX[sp+j]+round(sqrt(Source_thetaX[sp+j]));
                source_thetaY[sp+j] = Source_thetaY[sp+j]+round(sqrt(Source_thetaY[sp+j]));
                source_thetaX[sp+j] = Source_thetaX[sp+j]*1.1;
                source_thetaY[sp+j] = Source_thetaY[sp+j]*1.1;
                sX->DeconvolutionRL(source_thetaX, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
                sY->DeconvolutionRL(source_thetaY, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
                for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                    dataset_mod_jjj[n] = source_thetaX[n];
                }
                float dataset_mod_jjj_total = 0;
                for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++) dataset_mod_jjj_total += dataset_mod_jjj[n];

                float cov_calc_ij = 0;
                float cov_calc_jjj = 0;
                float cov_calc_ij_ori = 0;
                float cov_calc_ij_bin = 0;
                for( int n=sp+gap+10; n<sp+gap+39; n++ ){ 
                    cov_calc_jjj += (dataset_ori_lll[n] - dataset_mod_lll[n])*(dataset_ori_jjj[n] - dataset_mod_jjj[n]);
                }
                for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++ ){
                    cov_calc_ij += (dataset_ori_lll[n] - dataset_mod_lll[n])*(dataset_ori_jjj[n] - dataset_mod_jjj[n]);
                }
                /* cov_calc_ij=cov_calc_ij/(dataset_ori_jjj_total+dataset_mod_jjj_total); */
                cov_calc_ij_ori=cov_calc_ij/(dataset_mod_jjj_total);
                cov_calc_matrix(rrr-10,j-10) = cov_calc_ij_ori;
                cov_calc_matrix(j-10,rrr-10) = cov_calc_ij_ori;
                
                cov_calc_ij=cov_calc_ij;
                cov_calc_matrix_nonorm(rrr-10,j-10) = cov_calc_ij;
                cov_calc_matrix_nonorm(j-10,rrr-10) = cov_calc_ij;
                
                cov_calc_ij_bin=cov_calc_ij/29;
                cov_calc_matrix_binnorm(rrr-10,j-10) = cov_calc_ij_bin;
                cov_calc_matrix_binnorm(j-10,rrr-10) = cov_calc_ij_bin;
                
                cov_calc_jjj=cov_calc_jjj;
                cov_calc_matrix_jjj(rrr-10,j-10) = cov_calc_jjj;
                cov_calc_matrix_jjj(j-10,rrr-10) = cov_calc_jjj;
            
                /* modA_ij = (dataset_mod_jjj[sp+rrr] - dataset_ori_jjj[sp+rrr])/(Source_thetaX[sp+j]*1.1 - Source_thetaX[sp+j]); */
                /* modA(rrr-10,j-10) = modA_ij; */
                /* modA_ij = (dataset_mod_lll[sp+j] - dataset_ori_lll[sp+j])/(Source_thetaX[sp+rrr]*1.1 - Source_thetaX[sp+rrr]); */
                /* modA(j-10,rrr-10) = modA_ij; */
                
                
                myA_ij = (dataset_mod_lll[sp+j] - dataset_ori_lll[sp+j])/(Source_thetaX[sp+j]*1.1);
                myA(j-10,rrr-10) = myA_ij;
                myA(rrr-10,j-10) = myA_ij;
                
                modA_ij = (dataset_mod_lll[sp+rrr] - dataset_ori_lll[sp+rrr])*(dataset_mod_lll[sp+j] - dataset_ori_lll[sp+j])/(Source_thetaX[sp+j]*1.1+Source_thetaX[sp+rrr]*1.1);
                modA(j-10,rrr-10) = modA_ij;
                modA(rrr-10,j-10) = modA_ij;
                
                A_ij = (dataset_mod_lll[sp+j] - dataset_ori_lll[sp+j])/(Source_thetaX[sp+j]*1.1-Source_thetaX[sp+j]);
                A(j-10,rrr-10) = A_ij;
                A(rrr-10,j-10) = A_ij;
                
                /* modA_ij = (dataset_mod_lll[sp+rrr] - dataset_ori_lll[sp+rrr])/(Source_thetaX[sp+j]*1.1 - Source_thetaX[sp+j]); */
                /* modA(rrr-10,j-10) = modA_ij; */
                /* modA_ij = (dataset_mod_jjj[sp+j] - dataset_ori_jjj[sp+j])/(Source_thetaX[sp+rrr]*1.1 - Source_thetaX[sp+rrr]); */
                /* modA(j-10,rrr-10) = modA_ij; */
                /* std::cout << "rrr-10 " << rrr-10 << std::endl; */
                /* std::cout << "j-10 " << j-10 << std::endl; */
                /* cov_calc_matrix(rrr-10,j-10) = cov_calc_ij; */

            for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
                source_thetaX[n] = Source_thetaX[n];
                source_thetaY[n] = Source_thetaY[n];
            }
            }
            sX->Clear();
            sY->Clear();
        }
                cov_calc_matrix.Print();
                std::cout << "modA" << std::endl;
                modA.Print();

        for( int n=0; n<sp+gap+_histlimits["NbinsXY"]; n++){
            source_thetaX[n] = Source_thetaX[n];
            source_thetaY[n] = Source_thetaY[n];
        }
        for(int i=0; i<int(_histlimits["NbinsXY"]+gap); i++){
            thetaX_recoGoldhold2->SetBinContent(i+1, source_thetaX[sp+i]);
            thetaY_recoGoldhold2->SetBinContent(i+1, source_thetaY[sp+i]);
        }
        for(int i=0; i<int(_histlimits["NbinsXY"]); i++){
            thetaX_recoGoldhold3->SetBinContent(i+1, source_thetaX[sp+i]);
            thetaY_recoGoldhold3->SetBinContent(i+1, source_thetaY[sp+i]);
        }
        sX->DeconvolutionRL(source_thetaX, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
        sY->DeconvolutionRL(source_thetaY, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
        sX_noeffi->DeconvolutionRL(source_thetaX_noeffi, response_thetaX_noeffi, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
        sY_noeffi->DeconvolutionRL(source_thetaY_noeffi, response_thetaY_noeffi, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);


        float *XforChi = new float[itMax-1];
        for (int i=0; i<itMax-1; i++){
            XforChi[i] = i+2;
        }

        for(int i=0; i<int(_histlimits["NbinsXY"]+gap); i++){
            if (::isnan(source_thetaX[sp+i])==0 && ::isinf(source_thetaX[sp+i])==0) {
                thetaX_recoGoldhold->SetBinContent(i+1, source_thetaX[sp+i]);
                thetaX_recoGoldhold->SetBinError(i+1, sqrt(pow(thetaX_recoGoldhold->GetBinError(i+1),2)+pow(diffbinX[i],2)));
                thetaX_recoGoldhold_noeffi->SetBinContent(i+1, source_thetaX_noeffi[sp+i]);
                /* std::cout << "has content" << std::endl; */
            }
            else { 
                thetaX_recoGoldhold->SetBinContent(i+1, 0);
                thetaX_recoGoldhold_noeffi->SetBinContent(i+1, 0);
            }
            /* std::cout << "after source_thetaX[sp+i] " << source_thetaX[sp+i] << std::endl; */
            if (::isnan(source_thetaY[sp+i])==0 && ::isinf(source_thetaY[sp+i])==0) {
                thetaY_recoGoldhold->SetBinContent(i+1, source_thetaY[sp+i]);
                thetaY_recoGoldhold->SetBinError(i+1, sqrt(pow(thetaX_recoGoldhold->GetBinError(i+1),2)+pow(diffbinY[i],2)));
                thetaY_recoGoldhold_noeffi->SetBinContent(i+1, source_thetaY_noeffi[sp+i]);
            }
           else { 
                thetaY_recoGoldhold->SetBinContent(i+1, 0);
                thetaY_recoGoldhold_noeffi->SetBinContent(i+1, 0);

            }
        }

        for( int i=0; i<sp+gap+_histlimits["NbinsXY"]; i++){
            source_thetaX[i] = Source_thetaX[i];
            source_thetaY[i] = Source_thetaY[i];
            source_thetaX_noeffi[i] = Source_thetaX_noeffi[i];
            source_thetaY_noeffi[i] = Source_thetaY_noeffi[i];
        }

        sX_effi_only_gold->DeconvolutionRL(source_thetaX, response_thetaX_noeffi, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
        sY_effi_only_gold->DeconvolutionRL(source_thetaY, response_thetaY_noeffi, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
        sX_effi_only_response->DeconvolutionRL(source_thetaX_noeffi, response_thetaX, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);
        sY_effi_only_response->DeconvolutionRL(source_thetaY_noeffi, response_thetaY, sp+gap+int(_histlimits["NbinsXY"]), 10, 20, 0);

        for(int i=0; i<int(_histlimits["NbinsXY"]+gap); i++){
            if (::isnan(source_thetaX[sp+i])==0 && ::isinf(source_thetaX[sp+i])==0) {
                thetaX_recoGoldhold_effi_only_gold->SetBinContent(i+1, source_thetaX[sp+i]);
                thetaX_recoGoldhold_effi_only_response->SetBinContent(i+1, source_thetaX_noeffi[sp+i]);
                /* std::cout << "has content" << std::endl; */
            }
            else { 
                thetaX_recoGoldhold_effi_only_gold->SetBinContent(i+1, 0);
                thetaX_recoGoldhold_effi_only_response->SetBinContent(i+1, 0);
            }
            /* std::cout << "after source_thetaX[sp+i] " << source_thetaX[sp+i] << std::endl; */
            if (::isnan(source_thetaY[sp+i])==0 && ::isinf(source_thetaY[sp+i])==0) {
                thetaY_recoGoldhold_effi_only_gold->SetBinContent(i+1, source_thetaY[sp+i]);
                thetaY_recoGoldhold_effi_only_response->SetBinContent(i+1, source_thetaY_noeffi[sp+i]);
            }
            else { 
                thetaY_recoGoldhold_effi_only_gold->SetBinContent(i+1, 0);
                thetaY_recoGoldhold_effi_only_response->SetBinContent(i+1, 0);

            }
        }

        int meanpx = 0;
        int i=0;
        double meanx = double(thetaX_recoGoldhold->GetMean());
        while (meanpx==0 && ::isnan(meanx)==0) {
            /* std::cout << "in1" << std::endl; */
            /* std::cout << "thetaX_recoGoldhold->GetBinLowEdge(i) " << thetaX_recoGoldhold->GetBinLowEdge(i) << std::endl; */
            /* std::cout << "thetaX_recoGoldhold->GetMean() " << thetaX_recoGoldhold->GetMean() << std::endl; */
            if (thetaX_recoGoldhold->GetBinLowEdge(i)>thetaX_recoGoldhold->GetMean()) meanpx=i;
            i++;
        }
        int meanpy = 0;
        i=0;
        double meany = double(thetaY_recoGoldhold->GetMean());
        while (meanpy==0 && ::isnan(meany)==0) {
            if (thetaY_recoGoldhold->GetBinLowEdge(i)>thetaY_recoGoldhold->GetMean()) meanpy=i;
            i++;
        }
        int meanpx_noeffi = 0;
        i=0;
        double meanx_noeffi = double(thetaX_recoGoldhold_noeffi->GetMean());
        while (meanpx_noeffi==0 && ::isnan(meanx_noeffi)==0) {
            if (thetaX_recoGoldhold_noeffi->GetBinLowEdge(i)>thetaX_recoGoldhold_noeffi->GetMean()) meanpx_noeffi=i;
            i++;
        }
        int meanpy_noeffi = 0;
        i=0;
        double meany_noeffi = double(thetaY_recoGoldhold_noeffi->GetMean());
        while (meanpy_noeffi==0 && ::isnan(meany_noeffi)==0) {
            if (thetaY_recoGoldhold_noeffi->GetBinLowEdge(i)>thetaY_recoGoldhold_noeffi->GetMean()) meanpy_noeffi=i;
            i++;
        }
        int meanpx_only_gold = 0;
        i=0;
        double meanx_only_gold = double(thetaX_recoGoldhold_effi_only_gold->GetMean());
        while (meanpx_only_gold==0 && ::isnan(meanx_only_gold)==0) {
            if (thetaX_recoGoldhold_effi_only_gold->GetBinLowEdge(i)>thetaX_recoGoldhold_effi_only_gold->GetMean()) meanpx_only_gold=i;
            i++;
        }
        int meanpy_only_gold = 0;
        i=0;
        double meany_only_gold = double(thetaY_recoGoldhold_effi_only_gold->GetMean());
        while (meanpy_only_gold==0 && ::isnan(meany_only_gold)==0) {
            if (thetaY_recoGoldhold_effi_only_gold->GetBinLowEdge(i)>thetaY_recoGoldhold_effi_only_gold->GetMean()) meanpy_only_gold=i;
            i++;
        }
        int meanpx_only_response = 0;
        i=0;
        double meanx_only_response = double(thetaX_recoGoldhold_effi_only_response->GetMean());
        while (meanpx_only_response==0 && ::isnan(meanx_only_response)==0) {
            if (thetaX_recoGoldhold_effi_only_response->GetBinLowEdge(i)>thetaX_recoGoldhold_effi_only_response->GetMean()) meanpx_only_response=i;
            i++;
        }
        int meanpy_only_response = 0;
        i=0;
        double meany_only_response = double(thetaY_recoGoldhold_effi_only_response->GetMean());
        while (meanpy_only_response==0 && ::isnan(meany_only_response)==0) {
            if (thetaY_recoGoldhold_effi_only_response->GetBinLowEdge(i)>thetaY_recoGoldhold_effi_only_response->GetMean()) meanpy_only_response=i;
            i++;
        }

        for (i = 1; i < 48; i++) thetaX_recoGold->SetBinContent(i,thetaX_recoGoldhold->GetBinContent(i-1-(48)/2+meanpx-offset));
        for (i = 1; i < 48; i++) std::cout << "entry " << thetaX_recoGoldhold->GetBinContent(i-1-(48)/2+meanpx-offset) << std::endl;
        for (i = 1; i < 48; i++) std::cout << "offset " << offset << std::endl;
        for (i = 1; i < 48; i++) std::cout << "meanpx " << meanpx << std::endl;
        for (i = 1; i < 48; i++) thetaY_recoGold->SetBinContent(i,thetaY_recoGoldhold->GetBinContent(i-1-(48)/2+meanpy-offset));
        for (i = 1; i < 48; i++) thetaX_recoGold_noeffi->SetBinContent(i,thetaX_recoGoldhold_noeffi->GetBinContent(i-1-(48)/2+meanpx_noeffi-offset));
        for (i = 1; i < 48; i++) thetaX_recoGoldhold_noeffi->SetBinContent(i,thetaX_recoGold_noeffi->GetBinContent(i));
        for (i = 1; i < 48; i++) thetaY_recoGold_noeffi->SetBinContent(i,thetaY_recoGoldhold_noeffi->GetBinContent(i-1-(48)/2+meanpy_noeffi-offset));
        for (i = 1; i < 48; i++) thetaY_recoGoldhold_noeffi->SetBinContent(i,thetaY_recoGold_noeffi->GetBinContent(i));
        for (i = 1; i < 48; i++) thetaX_recoGold_effi_only_gold->SetBinContent(i,thetaX_recoGoldhold_effi_only_gold->GetBinContent(i-1-(48)/2+meanpx_only_gold-offset));
        for (i = 1; i < 48; i++) thetaY_recoGold_effi_only_gold->SetBinContent(i,thetaY_recoGoldhold_effi_only_gold->GetBinContent(i-1-(48)/2+meanpy_only_gold-offset));
        for (i = 1; i < 48; i++) thetaX_recoGold_effi_only_response->SetBinContent(i,thetaX_recoGoldhold_effi_only_response->GetBinContent(i-1-(48)/2+meanpx_only_response-offset));
        for (i = 1; i < 48; i++) thetaY_recoGold_effi_only_response->SetBinContent(i,thetaY_recoGoldhold_effi_only_response->GetBinContent(i-1-(48)/2+meanpy_only_response-offset));
        /*
           for (int i = 1; i < int(_histlimits["NbinsXY"])+1; i++) {
           thetaX_recoGold->SetBinContent(i,(thetaX_recoGoldhold->GetBinContent(i)+thetaX_recoGoldhold->GetBinContent(int(_histlimits["NbinsXY"])+1-i))/2);
           thetaY_recoGold->SetBinContent(i,(thetaY_recoGoldhold->GetBinContent(i)+thetaY_recoGoldhold->GetBinContent(int(_histlimits["NbinsXY"])+1-i))/2);
           }
           */

        TSpectrum *sScatt = new TSpectrum();
        //sScatt->Deconvolution(source_thetaScatt, response_thetaScatt, 2*(gap+sp)+int(_histlimits["NbinsTh"]), 20, 10, 0.0);

        for( int i=0; i<2*(sp+gap)+_histlimits["NbinsTh"]; i++){
            Source_thetaScatt[i] = source_thetaScatt[i];
        }
        for(int i=1; i<=itMax; i++){
            for(int n=0; n<2*(sp+gap)+_histlimits["NbinsTh"]; n++){
                source_thetaScatt[n] = Source_thetaScatt[n];
            }
            sScatt->DeconvolutionRL(source_thetaScatt, response_thetaScatt, 2*(gap+sp)+int(_histlimits["NbinsTh"]), i, 20, 0);
            for(int n=0; n<_histlimits["NbinsTh"]; n++){
                hnScatt->SetBinContent(n+1,source_thetaScatt[2*sp+n]);
            }
            if(i!=1){
                ChiSqScatt[i-2] = float(hnScatt->Chi2Test(hpScatt,"CHI2"));
                kolmogorovScatt[i-2] = float(hnScatt->KolmogorovTest(hpScatt,"M"));
            }
            hnScatt->Clear();
            hpScatt->Clear();
            for(int n=0; n<_histlimits["NbinsTh"]; n++){
                hpScatt->SetBinContent(n+1,source_thetaScatt[2*sp+n]);
            }
            sScatt->Clear();
        }
        for(int n=0; n<2*(sp+gap)+_histlimits["NbinsTh"]; n++){
            source_thetaScatt[n] = Source_thetaScatt[n];
        }
        sScatt->DeconvolutionRL(source_thetaScatt, response_thetaScatt, 2*(gap+sp)+int(_histlimits["NbinsTh"]), 10, 20, 0);
        for(int i=0; i< int(_histlimits["NbinsTh"]); i++){
            thetaScatt_recoGold->SetBinContent(i+1, source_thetaScatt[ 2*sp+i]);
        }
        TSpectrum *sScatt2 = new TSpectrum();
        //sScatt2->Deconvolution(source_theta2Scatt, response_theta2Scatt, sp+gap+int(_histlimits["NbinsTh2"]), 20, 10, 0.0);
        //New verion
        for( int i=0; i<2*(sp2+gap2)+_histlimits["NbinsTh2"]; i++){
            Source_theta2Scatt[i] = source_theta2Scatt[i];
        }
        for(int i=1; i<=itMax; i++){
            for(int n=0; n<2*(sp2+gap2)+_histlimits["NbinsTh2"]; n++){
                source_theta2Scatt[n] = Source_theta2Scatt[n];
            }
            sScatt2->DeconvolutionRL(source_theta2Scatt, response_theta2Scatt, 2*(sp2+gap2)+int(_histlimits["NbinsTh2"]), i, 20, 0);
            for(int n=0; n<_histlimits["NbinsTh2"]; n++){
                hn2Scatt->SetBinContent(n+1,source_theta2Scatt[2*sp2+n]);
            }
            if(i!=1){
                ChiSq2Scatt[i-2] = float(hn2Scatt->Chi2Test(hp2Scatt,"CHI2"));
                kolmogorov2Scatt[i-2] = float(hn2Scatt->KolmogorovTest(hp2Scatt,"M"));
            }
            hn2Scatt->Clear();
            hp2Scatt->Clear();
            for(int n=0; n<_histlimits["NbinsTh2"]; n++){
                hp2Scatt->SetBinContent(n+1,source_theta2Scatt[2*sp2+n]);
            }
            sScatt2->Clear();
        }
        for(int n=0; n<2*(sp+gap)+_histlimits["NbinsTh2"]; n++){
            source_theta2Scatt[n] = Source_theta2Scatt[n];
        }
        sScatt2->DeconvolutionRL(source_theta2Scatt, response_theta2Scatt, 2*(sp2+gap2)+int(_histlimits["NbinsTh2"]), 10, 20, 0);
        for(int i=0; i<int(_histlimits["NbinsTh2"]*2); i++){
            theta2Scatt_recoGoldhold->SetBinContent(i+1, source_theta2Scatt[sp2+i-46]);
        }
        int meanp2scatt = 0;
        i=0;
        while (meanp2scatt==0) {
            if (theta2Scatt_recoGoldhold->GetBinLowEdge(i)>theta2Scatt_recoGoldhold->GetMean()) meanp2scatt=i;
            i++;
        }
        for (i = 1; i < 47; i++) theta2Scatt_recoGold->SetBinContent(i,theta2Scatt_recoGoldhold->GetBinContent(i+meanp2scatt));
        /* 
           for( int i=0; i<sp+gap+_histlimits["NbinsTh2"]; i++){
           Source_theta2Scatt[i] = source_theta2Scatt[i];
           }
           for(int i=1; i<=itMax; i++){
           for(int n=0; n<sp+gap+_histlimits["NbinsTh2"]; n++){
           source_theta2Scatt[n] = Source_theta2Scatt[n];
           }
           sScatt2->DeconvolutionRL(source_theta2Scatt, response_theta2Scatt, sp+gap+int(_histlimits["NbinsTh2"]), i, 20, 0);
           for(int n=0; n<_histlimits["NbinsTh2"]; n++){
           hn2Scatt->SetBinContent(n+1,source_theta2Scatt[sp+n]);
           }
           if(i!=1){
           ChiSq2Scatt[i-2] = float(hn2Scatt->Chi2Test(hp2Scatt,"CHI2"));
           kolmogorov2Scatt[i-2] = float(hn2Scatt->KolmogorovTest(hp2Scatt,"M"));
           }
           hn2Scatt->Clear();
           hp2Scatt->Clear();
           for(int n=0; n<_histlimits["NbinsTh2"]; n++){
           hp2Scatt->SetBinContent(n+1,source_theta2Scatt[sp+n]);
           }
           sScatt2->Clear();
           }
           for(int n=0; n<sp+gap+_histlimits["NbinsTh2"]; n++){
           source_theta2Scatt[n] = Source_theta2Scatt[n];
           }
           sScatt2->DeconvolutionRL(source_theta2Scatt, response_theta2Scatt, sp+gap+int(_histlimits["NbinsTh2"]), 10, 20, 0);
           */
        TGraph *ChiSquareX = new TGraph(itMax-1, XforChi, ChiSqX);
        TGraph *KolmogorovX = new TGraph(itMax-1, XforChi, kolmogorovX);
        TGraph *ChiSquareY = new TGraph(itMax-1, XforChi, ChiSqY);
        TGraph *KolmogorovY = new TGraph(itMax-1, XforChi, kolmogorovY);
        TGraph *ChiSquareScatt = new TGraph(itMax-1, XforChi, ChiSqScatt);
        TGraph *KolmogorovScatt = new TGraph(itMax-1, XforChi, kolmogorovScatt);
        TGraph *ChiSquare2Scatt = new TGraph(itMax-1, XforChi, ChiSq2Scatt);
        TGraph *Kolmogorov2Scatt = new TGraph(itMax-1, XforChi, kolmogorov2Scatt);
        TGraph *kurtosisX = new TGraph(itMax-1, XforChi, kurtX);
        TGraph *kurtosisY = new TGraph(itMax-1, XforChi, kurtY);
        /*
           for(int i=0; i<int(_histlimits["NbinsTh2"]); i++){
           theta2Scatt_recoGold->SetBinContent(i+1, source_theta2Scatt[sp+i]);
           }
           */


        //thetaX_recoGold = trkreffix(thetaX_recoGold);
        //thetaY_recoGold = trkreffiy(thetaY_recoGold);
        //thetaScatt_recoGold = trkreffiscatt(thetaScatt_recoGold);
        //theta2Scatt_recoGold = trkreffi2scatt(theta2Scatt_recoGold);
        outfile->cd();
        cov_calc_matrix.Write("cov_matrix");
        cov_calc_matrix_nonorm.Write("cov_matrix_nonorm");
        cov_calc_matrix_binnorm.Write("cov_matrix_binnorm");
        cov_calc_matrix_jjj.Write("cov_matrix_jjj");
        modA.Write("A");
        myA.Write("myA");
        A.Write("rogersA");
        epsilon.Write("epsilon");
        NoEntriesGold.Write("recoGold");
        NoEntriesGold.Write("recoGold");
        NoEntriesGold.Write("measdata");
        NoEntriesGold.Write("Moliere");
        NoEntriesresp.Write("ref");
        tX_rG_int=thetaX_recoGold->Integral();
        
        thetaX_recoGold->Write();
        thetaY_recoGold->Write();
        thetaX_recoGoldhold->Write();
        thetaY_recoGoldhold->Write();
        thetaX_recoGoldhold2->Write();
        thetaY_recoGoldhold2->Write();
        
        thetaX_recoGoldhold_noeffi->Write();
        thetaY_recoGoldhold_noeffi->Write();

        thetaX_recoGold_noeffi->Write();
        thetaY_recoGold_noeffi->Write();
        
        thetaX_recoGold_noeffi->SetName("thetaX_recoGoldhold3");
        thetaY_recoGold_noeffi->SetName("thetaY_recoGoldhold3");
        thetaX_recoGoldhold3 = trkreffix(thetaX_recoGold_noeffi, 0);
        thetaY_recoGoldhold3 = trkreffiy(thetaY_recoGold_noeffi, 0);
        thetaX_recoGoldhold3->Write();
        thetaY_recoGoldhold3->Write();
        
        thetaX_recoGoldhold3->SetName("thetaX_recoGoldhold3_sym");
        thetaY_recoGoldhold3->SetName("thetaY_recoGoldhold3_sym");
        for (int i = 1; i < _histlimits["NbinsXY"]/2+1; i++) thetaX_recoGoldhold3->SetBinContent(i,(thetaX_recoGoldhold3->GetBinContent(i)+thetaX_recoGoldhold3->GetBinContent(_histlimits["NbinsXY"]+1-i))/2);
        for (int i = _histlimits["NbinsXY"]/2+1; i < _histlimits["NbinsXY"]+1; i++) thetaY_recoGoldhold3->SetBinContent(i,(thetaY_recoGoldhold3->GetBinContent(i)+thetaY_recoGoldhold3->GetBinContent(_histlimits["NbinsXY"]+1-i))/2);
        thetaX_recoGoldhold3->Write();
        thetaY_recoGoldhold3->Write();
        
        thetaX_recoGold_effi_only_gold->Write();
        thetaY_recoGold_effi_only_gold->Write();
        thetaX_recoGold_effi_only_response->Write();
        thetaY_recoGold_effi_only_response->Write();
        thetaScatt_recoGold->Write();
        theta2Scatt_recoGold->Write();
        theta2Scatt_recoGoldhold->Write();
        thetaX_response->Write();
        thetaY_response->Write();
        thetaScatt_response->Write();
        theta2Scatt_response->Write();

        ChiSquareX->SetName("ChiSquareX");
        ChiSquareX->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 repetitions of deconvolution of #it{#theta_{X}} (boost = 0.75)");
        ChiSquareX->GetXaxis()->SetTitle("n");
        ChiSquareX->Write();
        ChiSquareY->SetName("ChiSquareY");
        ChiSquareY->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 repetitions of deconvolution of #it{#theta_{Y}} (boost = 0.75)");
        ChiSquareY->GetXaxis()->SetTitle("n");
        ChiSquareY->Write();
        ChiSquareScatt->SetName("ChiSquareScatt");
        ChiSquareScatt->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 repetitions of deconvolution of #it{#theta_{Scatt}} (boost = 0.75)");
        ChiSquareScatt->GetXaxis()->SetTitle("n");
        ChiSquareScatt->Write();
        ChiSquare2Scatt->SetName("ChiSquare2Scatt");
        ChiSquare2Scatt->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 repetitions of deconvolution of #theta^{2}_{Scatt} (boost = 0.75)");
        ChiSquare2Scatt->GetXaxis()->SetTitle("n");
        ChiSquare2Scatt->Write();

        KolmogorovX->SetName("KolmogorovX");
        KolmogorovX->Write();
        KolmogorovY->SetName("KolmogorovY");
        KolmogorovY->Write();
        KolmogorovScatt->SetName("KolmogorovScatt");
        KolmogorovScatt->Write();
        Kolmogorov2Scatt->SetName("Kolmogorov2Scatt");
        Kolmogorov2Scatt->Write();

        TCanvas* c2 = new TCanvas();
        c2->SetLogy();
        TLegend *leg = new TLegend(.65,.65,.85,.85);
        TMultiGraph *Chi2KolmX = new TMultiGraph();
        Chi2KolmX->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 iterations of deconvolution of #it{#theta_{X}};n;");
        ChiSquareX->SetLineColor(kRed);
        Chi2KolmX->Add(ChiSquareX, "c");
        Chi2KolmX->Add(KolmogorovX, "c");
        /* Chi2KolmX->Add(kurtosisX, "c"); */
        leg->AddEntry(ChiSquareX, "#Chi^{2}","l");
        leg->AddEntry(KolmogorovX, "Kol.","l");
        /* leg->AddEntry(kurtosisX, "Kurt.","l"); */
        Chi2KolmX->Draw("a");
        leg->Draw();
        c2->SaveAs("Chi2KolmX.pdf");
        c2->Clear();
        leg->Clear();
        TMultiGraph *Chi2KolmY = new TMultiGraph();
        Chi2KolmY->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 iterations of deconvolution of #it{#theta_{Y}};n;");
        ChiSquareY->SetLineColor(kRed);
        Chi2KolmY->Add(ChiSquareY, "c");
        Chi2KolmY->Add(KolmogorovY, "c");
        Chi2KolmY->Add(kurtosisY, "c");
        leg->AddEntry(ChiSquareY, "#Chi^{2}","l");
        leg->AddEntry(KolmogorovY, "Kol.","l");
        leg->AddEntry(kurtosisY, "kurt.","l");
        Chi2KolmY->Draw("a");
        leg->Draw();
        c2->SaveAs("Chi2KolmY.pdf");
        c2->Clear();
        leg->Clear();
        TMultiGraph *Chi2KolmScatt = new TMultiGraph();
        Chi2KolmScatt->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 iterations of deconvolution of #it{#theta_{Scatt}};n;");
        ChiSquareScatt->SetLineColor(kRed);
        Chi2KolmScatt->Add(ChiSquareScatt, "c");
        Chi2KolmScatt->Add(KolmogorovScatt, "c");
        leg->AddEntry(ChiSquareX, "#Chi^{2}","l");
        leg->AddEntry(KolmogorovX, "Kol.","l");
        Chi2KolmScatt->Draw("a");
        leg->Draw();
        c2->SaveAs("Chi2KolmScatt.pdf");
        c2->Clear();
        leg->Clear();
        TMultiGraph *Chi2Kolm2Scatt = new TMultiGraph();
        Chi2Kolm2Scatt->SetTitle("#chi^{2} and Kolmogorov tests for n and n-1 iterations of deconvolution of #theta^{2}_{Scatt};n;");
        ChiSquare2Scatt->SetLineColor(kRed);
        Chi2Kolm2Scatt->Add(ChiSquare2Scatt, "c");
        Chi2Kolm2Scatt->Add(Kolmogorov2Scatt, "c");
        leg->AddEntry(ChiSquareX, "#Chi^{2}","l");
        leg->AddEntry(KolmogorovX, "Kol.","l");
        Chi2Kolm2Scatt->Draw("a");
        leg->Draw();
        c2->SaveAs("Chi2Kolm2Scatt.pdf");
        c2->Clear();
        thetaX_recoGold_noeffi->Draw();
        c2->SaveAs("thetaX_recoGold_noeffi.pdf");
        c2->Clear();

        TCanvas* c1 = new TCanvas();
        TH1F *hgoldmirrorl = new TH1F("hgoldmirrorl","Mirror",_histlimits["NbinsXY"],-_histlimits["maxXY"],_histlimits["maxXY"]);
        TH1F *hgoldmirrorr = new TH1F("hgoldmirrorr","Mirror",_histlimits["NbinsXY"],-_histlimits["maxXY"],_histlimits["maxXY"]);
        for (int i = 1; i < _histlimits["NbinsXY"]/2+1; i++) hgoldmirrorl->SetBinContent(i,thetaX_recoGold->GetBinContent(i));
        for (int i = 1; i < _histlimits["NbinsXY"]/2+1; i++) hgoldmirrorr->SetBinContent(i,thetaX_recoGold->GetBinContent(_histlimits["NbinsXY"]+1-i));
        hgoldmirrorl->Draw();
        hgoldmirrorr->Draw("SAME");
        hgoldmirrorl->Write();
        hgoldmirrorr->Write();
        TH1F *emptygoldmirrorl = new TH1F("emptygoldmirrorl","Mirror",_histlimits["NbinsXY"],-_histlimits["maxXY"],_histlimits["maxXY"]);
        TH1F *emptygoldmirrorr = new TH1F("emptygoldmirrorr","Mirror",_histlimits["NbinsXY"],-_histlimits["maxXY"],_histlimits["maxXY"]);
        for (int i = 1; i < _histlimits["NbinsXY"]/2+1; i++) emptygoldmirrorl->SetBinContent(i,thetaX_recoGold_noeffi->GetBinContent(i));
        for (int i = 1; i < _histlimits["NbinsXY"]/2+1; i++) emptygoldmirrorr->SetBinContent(i,thetaX_recoGold_noeffi->GetBinContent(_histlimits["NbinsXY"]+1-i));
        emptygoldmirrorl->Write();
        emptygoldmirrorr->Write();


        delete Chi2KolmX;
        delete Chi2KolmY;
        delete Chi2KolmScatt;
        delete Chi2Kolm2Scatt;
        delete thetaX_recoGold;
        delete thetaY_recoGold;
        delete thetaX_recoGoldhold;
        delete thetaY_recoGoldhold;
        delete thetaScatt_recoGold;
        delete theta2Scatt_recoGold;

        delete thetaX_response;
        delete thetaY_response;
        delete thetaScatt_response;
        delete theta2Scatt_response;

    }

    void MCSAnalysis::DoVirtualDeconvolution(std::string model, int n_sel=1){
        // The basic methods associated with the RooUnfolding package
        // First generate a histogram of the measured data.

        int isfirst = 0;
        if (model.find(modelname1.c_str()) != std::string::npos)
            isfirst = 1;

        // int n_base = int(2 * DSTruthSetLiH.N() / (n_sel * (n_sel + 1)));
        // int k=0;
        // for (int j=0; j<n_sel; j++)
        // define the number of events under the current selection
        int curr_sel = int(DSTruthSetLiH.N());  // n_base * j;
        int curr_k = 0;



        std::string  tmpname = "thetaX_measdata";
        tmpname += model;
        TH1D* thetaX_measdata = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_asymm";
        tmpname += model;
        TH1* thetaX_asymm = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{XY}}; Events",_histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);
        tmpname = "thetaX_measdataminus";
        tmpname += model;
        TH1D* thetaX_measdataminus = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_measdata";
        tmpname += model;
        TH1D* thetaY_measdata = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_measdatafold";
        tmpname += model;
        TH1D* thetaY_measdatafold = 
            new TH1D(tmpname.c_str(),";#it{#theta_{Y}};#it{#theta_{Yi}}-#it{#theta_{Y46-i}}",_histlimits["NbinsXY"]/2, _histlimits["minXY"], 0);
        tmpname = "thetaY_asymm";
        tmpname += model;
        TH1* thetaY_asymm = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);
        tmpname = "thetaY_measdataminus";
        tmpname += model;
        TH1D* thetaY_measdataminus = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaScatt_measdata";
        tmpname += model;
        TH1D* thetaScatt_measdata = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        tmpname = "theta2Scatt_measdata";
        tmpname += model;
        TH1D* theta2Scatt_measdata = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
        tmpname = "thetaScatt_measdata_vp";
        tmpname += model;
        TH2D* thetaScatt_measdata_vp = 
            new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #it{#theta_{Scatt}}", 
                    400, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        tmpname = "thetaXUS_thetaXDS";
        tmpname += model;

        TH2D* thetaXUS_thetaXDS = 
            new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaYUS_thetaYDS";
        tmpname += model;
        TH2D* thetaYUS_thetaYDS = 
            new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "projposUSDSdiff_";
        tmpname += model;
        TH2D* projposUSDSdiff = 
            new TH2D(tmpname.c_str(),
                    "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
                    180, -400, 400, 180, -400, 400);
        tmpname = "asymmXYUS_";
        tmpname += model;
        TH2D* asymmXYUS = 
            new TH2D(tmpname.c_str(),";#it{x}(mm);#it{y}(mm)",100, -200, 200, 100, -200, 200);
        tmpname = "asymmXYDS_";
        tmpname += model;
        TH2D* asymmXYDS = 
            new TH2D(tmpname.c_str(),";#it{x}(mm);#it{y}(mm)",100, -200, 200, 100, -200, 200);

        //1D Histogram Settings
        TH1D *defineHist2(const char* name, const char* title, Double_t nbinsx, Double_t xlow, Double_t xup);
        TH1D *Histo[180];
        TH1D *Histoy[180];
        TH1D *Histodotw[180];
        TH1D *Histodotu[180];
        TH1D *Histoyd[180];
        TH1D *Histowd[180];
        TH1D *Histowx[180];
        TH1D *Histowy[180];
        TH1D *Histou[180];
        TH1D *Histos[180];
        TH1D *Histov[180];
        TH1D *Histow[180];
        TH1D *Histodotwdotu[180];

        std::cout << "_histlimits[NbinsXY] " << _histlimits["NbinsXY"] << std::endl;
        std::cout << "_histlimits[minXY] " << _histlimits["minXY"] << std::endl;
        std::cout << "_histlimits[maxXY] " << _histlimits["maxXY"] << std::endl;
        for (int l=0;l<180;l++){
            std::string  tmpname = "thetaX_measdata_";
            tmpname += l;
            Histo[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            //Histo[l] = defineHist2(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "thetaY_measdata_";
            tmpname += l;
            Histoy[l] = defineHist2(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "dotw_";
            tmpname += l;
            Histodotw[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "dotu_";
            tmpname += l;
            Histodotu[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "yd_";
            tmpname += l;
            Histoyd[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "wd_";
            tmpname += l;
            Histowd[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "wx_";
            tmpname += l;
            Histowx[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "wy_";
            tmpname += l;
            Histowy[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "u";
            tmpname += l;
            Histou[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "s";
            tmpname += l;
            Histos[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "dotv";
            tmpname += l;
            Histov[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "w";
            tmpname += l;
            Histow[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "dotwdotu_";
            tmpname += l;
            Histodotwdotu[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        }


        const Int_t NBINS = 19;
        Double_t scat_bin_array[NBINS + 1] = {-0.1151,-0.0938,-0.0754,-0.0597,-0.0463,-0.0347,-0.0248,-0.0162,-0.00895,-0.00269,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,0.0754,0.0938,0.1151};
        TH1D* scattering_proj_x = new TH1D("scattering_proj_x_DC","Change in Projected Angle (X);#it{#theta_{X}}; Events per radian", 
                NBINS, scat_bin_array);
        TH1D* scattering_proj_y = new TH1D("scattering_proj_y_DC","Change in Projected Angle (Y);#it{#theta_{Y}}; Events per radian", 
                NBINS, scat_bin_array);


        for(int i=curr_k; i<curr_sel; i++){
            //double hold = _USset.E(i).dXdz*cos(10)+sin(10);
            //_USset.E(i).dXdz = hold;
            std::vector<double> projTheta = RotDefineProjectionAngles(USTruthSetLiH.E(i), DSTruthSetLiH.E(i), angdef);
            double thetaY = projTheta[1];  /// atan(DSTruthSetLiH.E(i).dXdz) - atan(USTruthSetLiH.E(i).dXdz);
            double thetaX = projTheta[0];  /// atan(DSTruthSetLiH.E(i).dYdz) - atan(USTruthSetLiH.E(i).dYdz);
            // double cosScatt 
            double thetaScatt = projTheta[2];  /// acos(cosScatt);
            double xUSabs = USTruthSetLiH.E(i).X - USTruthSetLiH.E(i).dXdz * (USTruthSetLiH.E(i).Z - _sys["abspos"]);
            double yUSabs = USTruthSetLiH.E(i).Y - USTruthSetLiH.E(i).dYdz * (USTruthSetLiH.E(i).Z - _sys["abspos"]);
            double xDSabs = DSTruthSetLiH.E(i).X - DSTruthSetLiH.E(i).dXdz * (DSTruthSetLiH.E(i).Z - _sys["abspos"]);
            double yDSabs = DSTruthSetLiH.E(i).Y - DSTruthSetLiH.E(i).dYdz * (DSTruthSetLiH.E(i).Z - _sys["abspos"]);
            projposUSDSdiff->Fill(xDSabs - xUSabs, yDSabs - yUSabs);
            thetaXUS_thetaXDS->Fill(atan(USTruthSetLiH.E(i).dXdz), atan(DSTruthSetLiH.E(i).dXdz));
            thetaYUS_thetaYDS->Fill(atan(USTruthSetLiH.E(i).dYdz), atan(DSTruthSetLiH.E(i).dYdz));
            thetaX_measdata->Fill(thetaX);
            thetaY_measdata->Fill(thetaY);
            thetaY_measdataminus->Fill(-thetaY);
            thetaX_measdataminus->Fill(-thetaX);
            if (thetaY>0.02) {
                asymmXYUS->Fill(USTruthSetLiH.E(i).X, USTruthSetLiH.E(i).Y);
                asymmXYDS->Fill(DSTruthSetLiH.E(i).X, DSTruthSetLiH.E(i).Y);
            }
            thetaScatt_measdata->Fill(thetaScatt);
            theta2Scatt_measdata->Fill(thetaScatt*thetaScatt);
            // double dt0 = 7.64186 / 0.299792458; // ns. 
            // double dt = USTruthSetLiH.E(i).TOF01;
            // if (USTruthSetLiH.E(i).TOF12 < 99.0 * 8.22475 / 0.299792458){
            //  dt0 = 8.22475 / 0.299792458; // ns.
            //  dt = USTruthSetLiH.E(i).TOF12;
            //}
            double pz  = USTruthSetLiH.E(i).pz;
            // 105.65 / sqrt(dt*dt/dt0/dt0 - 1.0);
            thetaScatt_measdata_vp->Fill(pz, thetaScatt);
            scattering_proj_x->Fill(thetaX);
            scattering_proj_y->Fill(thetaY);
            // k++;
            for (int l=0;l<180;l++){
                //if (DSTruthSetLiH.E(i).dXdz<0.2 && DSTruthSetLiH.E(i).dYdz<0.2) {
                std::vector<double> RotprojTheta = RotDefineProjectionAngles(USTruthSetLiH.E(i), DSTruthSetLiH.E(i),l);
                double RotthetaY = RotprojTheta[1];
                double RotthetaX = RotprojTheta[0];
                //std::cout << "RotprojTheta[1] " << RotprojTheta[1] << std::endl;
                Histo[l]->Fill(RotthetaX);
                Histoy[l]->Fill(RotthetaY);
                //if (RotprojTheta[3]>0.2) Histodotw[l]->Fill(RotprojTheta[3]); 
                //if (RotprojTheta[3]<-0.2) Histodotw[l]->Fill(RotprojTheta[3]); 
                //if (RotprojTheta[4]>0.2) Histodotv[l]->Fill(RotprojTheta[4]); 
                //if (RotprojTheta[4]<-0.2) Histodotv[l]->Fill(RotprojTheta[4]); 
                Histodotw[l]->Fill(RotprojTheta[3]); 
                //Histodotw[l]->Fill(RotprojTheta[3]); 
                Histodotu[l]->Fill(RotprojTheta[4]); 
                Histoyd[l]->Fill(RotprojTheta[5]); 
                Histowd[l]->Fill(RotprojTheta[6]); 
                //if (i==0) { 
                Histowx[l]->Fill(RotprojTheta[8]); 
                Histowy[l]->Fill(RotprojTheta[9]); 
                Histou[l]->Fill(RotprojTheta[10]);
                Histos[l]->Fill(RotprojTheta[11]);
                //}
                Histov[l]->Fill(RotprojTheta[12]);
                Histow[l]->Fill(RotprojTheta[13]);
                Histodotwdotu[l]->Fill(RotprojTheta[14]); 
                //Histov_x[l]->Fill(RotprojTheta[8[0]]); 
                //Histov_y[l]->Fill(RotprojTheta[8[0]]); 
                //	    }
            }
        }

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

        /*
           TH1D* thetaX_reco_n; 
           TH1D* thetaX_measdata_n; 
           if (model.find(modelname3.c_str()) != std::string::npos) {
           double chi2[100];
           double chi2effi[100];
           double xaxis[100];
           for ( int itr=1; itr<100; itr++){    
           RooUnfoldBayes unfold_thetaX_n(&fresp_thetaX, thetaX_measdata, itr);
           thetaX_reco_n = (TH1D*)unfold_thetaX_n.Hreco();
           chi2[itr]=thetaX_reco_n->Chi2Test(theta_true_x_graph,"WUP");
           xaxis[itr] = itr;
           thetaX_measdata_n = trkreffix(thetaX_measdata);
           RooUnfoldBayes unfold_thetaX_eff(&fresp_thetaX, thetaX_measdata_n, itr);
           thetaX_reco_n = (TH1D*)unfold_thetaX_eff.Hreco();
           chi2effi[itr]=thetaX_reco_n->Chi2Test(theta_true_x_graph,"WUP");
           }
           TCanvas* c2 = new TCanvas("c2","A Simple Graph Example",200,10,1400,1000);
           TGraph* convergence = new TGraph(100,xaxis,chi2);
           convergence->GetXaxis()->SetTitle("iterations");
           convergence->GetYaxis()->SetTitle("#chi^{2} MC reco & Truth");
           convergence->Draw();
           c2->SaveAs("convergence.pdf");
           c2->SaveAs("convergence.root");
           c2->Clear();
           TGraph* convergenceeffi = new TGraph(100,xaxis,chi2effi);
           convergenceeffi->GetXaxis()->SetTitle("iterations");
           convergenceeffi->GetYaxis()->SetTitle("#chi^{2} MC reco & Truth");
           convergenceeffi->Draw();
           c2->SaveAs("convergenceeffi.pdf");
           c2->SaveAs("convergenceeffi.root");
           delete c2;
           delete convergence;
           }
           */
        /*
           thetaX_measdata = trkreffix(thetaX_measdata);
           thetaY_measdata = trkreffiy(thetaY_measdata);
           thetaScatt_measdata = trkreffiscatt(thetaScatt_measdata);
           theta2Scatt_measdata = trkreffi2scatt(theta2Scatt_measdata);
           */

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
        thetaX_reco->SetTitle(";#Delta #it{#theta_{X}}");
        //thetaX_reco = trkreffix(thetaX_reco);
        TH1D* thetaX_reco_noeffi = (TH1D*)unfold_thetaX.Hreco();
        tmpname = "thetaX_reco_noeffi";
        tmpname += model;
        thetaX_reco_noeffi->SetName(tmpname.c_str());
        TH1D* thetaY_reco = (TH1D*)unfold_thetaY.Hreco();
        tmpname = "thetaY_reco";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaY_reco->SetName(tmpname.c_str());
        thetaY_reco->SetTitle(";#Delta #it{#theta_{Y}}");
        tmpname = "thetaY_reco_noeffi";
        TH1D* thetaY_reco_noeffi = (TH1D*)unfold_thetaY.Hreco();
        tmpname += model;
        thetaY_reco_noeffi->SetName(tmpname.c_str());
        //thetaY_reco = , "N"trkreffiy(thetaY_reco);
        TH1D* thetaScatt_reco = (TH1D*)unfold_thetaScatt.Hreco();
        tmpname = "thetaScatt_reco";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaScatt_reco->SetName(tmpname.c_str());
        thetaScatt_reco->SetTitle(";#it{#theta_{Scatt}}");
        //thetaScatt_reco = trkreffiscatt(thetaScatt_reco);
        TH1D* theta2Scatt_reco = (TH1D*)unfold_theta2Scatt.Hreco();
        tmpname = "theta2Scatt_reco";
        tmpname += model;
        // if(j>0) tmpname += j;
        theta2Scatt_reco->SetName(tmpname.c_str());
        theta2Scatt_reco->SetTitle(";#theta^{2}_{Scatt}");
        //theta2Scatt_reco = trkreffi2scatt(theta2Scatt_reco);

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
        thetaX_truth->SetTitle(";#Delta #it{#theta_{X}}");
        TH1D* thetaY_truth = (TH1D*)unfold_thetaY.response()->Htruth();
        tmpname = "thetaY_";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaY_truth->SetName(tmpname.c_str());
        thetaX_truth->SetTitle(";#Delta #it{#theta_{Y}}");
        TH1D* thetaScatt_truth = (TH1D*)unfold_thetaScatt.response()->Htruth();
        tmpname = "thetaScatt_";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaScatt_truth->SetName(tmpname.c_str());
        thetaScatt_truth->SetTitle(";#it{#theta_{Scatt}}");
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
        thetaX_measured->SetTitle(";#Delta #it{#theta_{X}}");
        TH1D* thetaY_measured = (TH1D*)unfold_thetaY.response()->Hmeasured();
        tmpname = "thetaY_measured";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaY_measured->SetName(tmpname.c_str());
        thetaY_measured->SetTitle(";#Delta #it{#theta_{Y}}");
        TH1D* thetaScatt_measured = (TH1D*)unfold_thetaScatt.response()->Hmeasured();
        tmpname = "thetaScatt_measured";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaScatt_measured->SetName(tmpname.c_str());
        thetaScatt_measured->SetTitle(";#it{#theta_{Scatt}}");
        TH1D* theta2Scatt_measured = (TH1D*)unfold_theta2Scatt.response()->Hmeasured();
        tmpname = "theta2Scatt_measured";
        tmpname += model;
        // if(j>0) tmpname += j;
        theta2Scatt_measured->SetName(tmpname.c_str());
        theta2Scatt_measured->SetTitle(";#theta^{2}_{Scatt}");


        TH2D* dotw = new TH2D("dotw","dotw", 180, 0, 180, 200, -1, 1);
        TH2D* dotv = new TH2D("dotv","dotv", 180, 0, 180, 200, -1, 1);
        TH2D* histoyd = new TH2D("histoyd","RMS with plane definition", 180, 0, 180, 100, 1, 3);
        TH2D* histowd = new TH2D("histowd","RMS with plane definition", 180, 0, 180, 100, 1, 3);
        TH2D* histowx = new TH2D("histowx","RMS with plane definition", 180, 0, 180, 100, -2, 2);
        TH2D* histowy = new TH2D("histowy","RMS with plane definition", 180, 0, 180, 100, -2, 2);

        outfile->cd();
        double x[180];
        double xer[180];
        double y[180];
        double yer[180];
        double xaxis[180];
        double xaxiser[180];
        double u[180];
        double s[180];
        double zerox[180];
        double onex[180];
        double v[180];
        double w[180];
        double ddotw[180];
        double ddotu[180];
        double ddotwdotu[180];
        for (int l=0;l<180;l++){
            //std::cout << "Histo[l]->GetRMS() " << Histo[l]->GetRMS() << std::endl;
            //std::cout << "Histoy[l]->GetRMS() " << Histoy[l]->GetRMS() << std::endl;
            xaxis[l]=l;
            xaxiser[l]=0;
            //Histo[l]->Fit("gaus", "","",-0.03,0.03);
            TF1 *myfunc = Histo[l]->GetFunction("gaus");
            x[l]=myfunc->GetParameter(2);
            //x[l]=Histo[l]->GetRMS();
            xer[l]=0.00054;
            y[l]=Histoy[l]->GetRMS();
            yer[l]=0.0006;
            u[l]=Histou[l]->GetRMS();
            s[l]=Histos[l]->GetRMS();
            v[l]=Histov[l]->GetRMS();
            w[l]=Histow[l]->GetRMS();
            ddotw[l]=Histodotw[l]->GetRMS();
            ddotu[l]=Histodotu[l]->GetRMS();
            zerox[l]=Histowx[l]->GetRMS();
            onex[l]=Histowy[l]->GetRMS();
            histoyd->Fill(l,Histoyd[l]->GetMean());
            histowd->Fill(l,Histowd[l]->GetMean());
            histowx->Fill(l,Histowx[l]->GetRMS());
            histowy->Fill(l,Histowy[l]->GetRMS());
            ddotwdotu[l]=Histodotwdotu[l]->GetRMS();
            //Histowx[l]->Write();
            //Histowy[l]->Write();
            //Histoyd[l]->Write(); 
            //Histowd[l]->Write(); 
            if (l>38 && l<45) {
                Histo[l]->Write();
                //Histoy[l]->Write();
            }
            //Histodotw[l]->Write(); 
            //Histodotu[l]->Write();
            //Histodotwdotu[l]->Write(); 
            //Histodotv[l]->Write(); 
        }
        //histowx->Write();
        //histowy->Write();
        //histoyd->Write();
        //histowd->Write();
        //dotw->Write();
        //dotv->Write();
        thetaXUS_thetaXDS->Write();
        thetaYUS_thetaYDS->Write();
        //thetaX_measdata = trkreffix(thetaX_measdata);
        outfile->cd();
        thetaX_measdata->Write();
        //thetaY_measdata = trkreffiy(thetaY_measdata);
        outfile->cd();
        thetaY_measdata->Write();
        for (int i=1;i<24;i++) {
            thetaY_measdatafold->SetBinContent(i,thetaY_measdata->GetBinContent(i)-thetaY_measdata->GetBinContent(46-i));
        }
        thetaY_measdatafold->Write();
        //thetaScatt_measdata = trkreffiscatt(thetaScatt_measdata);
        outfile->cd();
        thetaScatt_measdata->Write();
        //theta2Scatt_measdata = trkreffi2scatt(theta2Scatt_measdata);
        outfile->cd();
        theta2Scatt_measdata->Write();
        thetaScatt_measdata_vp->Write();
        thetaX_reco->Write();
        thetaY_reco->Write();
        thetaX_reco_noeffi->Write();
        thetaY_reco_noeffi->Write();
        thetaScatt_reco->Write();
        theta2Scatt_reco->Write();
        thetaX_response->Write();
        thetaY_response->Write();
        thetaScatt_response->Write();
        theta2Scatt_response->Write();
        //thetaX_truth->Write();
        //thetaY_truth->Write();
        //thetaScatt_truth->Write();
        //theta2Scatt_truth->Write();
        thetaX_measured->Write();
        thetaY_measured->Write();
        thetaScatt_measured->Write();
        theta2Scatt_measured->Write();
        scattering_proj_x->Write();
        scattering_proj_y->Write();
        projposUSDSdiff->Write();
        //convergence->Write();

        TCanvas* c1 = new TCanvas();
        thetaY_asymm=thetaY_measdata->GetAsymmetry(thetaY_measdataminus);
        thetaY_asymm->Draw();
        TF1*  fM2D = new TF1("fM2D", "[0]+x*[1]", 0, 0.02);
        //fM2D->SetParameter(0,0);
        //fM2D->SetParameter(1,1);
        //fM2D->SetParNames("A","B");
        thetaY_asymm->Fit(fM2D,"RES");
        fM2D->DrawCopy("SAME");
        //thetaY_asymm->Fit("pol1","RES","",0,0.02);
        thetaY_asymm->SetName("thetaY_asymm");
        thetaY_asymm->SetTitle("#it{#theta_Y} Asymmetry");
        thetaY_asymm->GetYaxis()->SetTitle("Asymmetry");
        thetaY_asymm->GetXaxis()->SetRange(0,0.06);
        thetaY_asymm->Write();
        c1->SaveAs("thetaY_asymm.pdf");
        c1->Clear();
        thetaX_asymm->SetStats(0);
        thetaX_asymm=thetaX_measdata->GetAsymmetry(thetaX_measdataminus);
        thetaX_asymm->Draw();
        //thetaX_asymm->Fit("pol1","RES","",0,0.02);
        thetaX_asymm->Fit(fM2D,"RES");
        fM2D->DrawCopy("SAME");
        thetaX_asymm->SetName("thetaX_asymm");
        thetaX_asymm->SetTitle("#it{#theta_X} Asymmetry");
        thetaX_asymm->GetYaxis()->SetTitle("Asymmetry");
        thetaX_asymm->GetXaxis()->SetRange(0,0.06);
        thetaX_asymm->Write();
        c1->SaveAs("thetaX_asymm.pdf");
        c1->Clear();
        asymmXYUS->Draw("colz");
        c1->Print("asymmXYUS.pdf");
        c1->Clear();
        asymmXYDS->Draw("colz");
        c1->Print("asymmXYDS.pdf");
        c1->Clear();
        thetaX_measdata->Draw();
        c1->Print("thetaX_measdata.pdf");
        c1->Clear();
        thetaY_measdata->Draw();
        c1->Print("thetaY_measdata.pdf");
        c1->Print("thetaY_measdata.root");
        c1->Clear();
        thetaY_measdatafold->Draw();
        c1->Print("thetaY_measdatafold.pdf");
        c1->Print("thetaY_measdatafold.root");
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
        thetaX_reco_noeffi->Draw();
        c1->Print("thetaX_reco_noeffi.pdf");
        c1->Clear();
        thetaY_reco_noeffi->Draw();
        c1->Print("thetaY_reco_noeffi.pdf");
        c1->Clear();
        thetaScatt_reco->Draw();
        c1->Print("thetaScatt_reco.pdf");
        c1->Clear();
        thetaX_response->Draw("colz");
        c1->Print("thetaX_response.pdf");
        c1->Clear();
        thetaY_response->Draw("colz");
        c1->Print("thetaY_response.pdf");
        c1->Clear();
        thetaScatt_response->Draw("colz");
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
        TGraphErrors* RotDefHis = new TGraphErrors(180, xaxis,x,xaxiser,xer);
        TGraphErrors* RotDefHisy = new TGraphErrors(180, xaxis,y,xaxiser,yer);
        TGraph* ugraph = new TGraph(180, xaxis,u);
        TGraph* sgraph = new TGraph(180, xaxis,s);
        TGraph* vgraph = new TGraph(180, xaxis,zerox);
        TGraph* wgraph = new TGraph(180, xaxis,onex);
        TGraph* ddotwgraph = new TGraph(180, xaxis,ddotw);
        TGraph* ddotugraph = new TGraph(180, xaxis,ddotu);
        TGraph* ddotwdotugraph = new TGraph(180, xaxis,ddotwdotu);
        RotDefHis->SetFillColor(6);
        RotDefHis->SetFillStyle(3005);
        RotDefHis->SetMarkerStyle(34);
        RotDefHis->SetTitle("");
        RotDefHis->Draw("AC*");
        RotDefHis->Write();
        RotDefHis->GetXaxis()->SetTitle("Angle around Z axis (#circ)");
        RotDefHis->GetYaxis()->SetTitle("RMS of scattering distribution");
        c1->SaveAs("RotDefHis.pdf");
        c1->Clear();
        RotDefHisy->SetFillColor(6);
        RotDefHisy->SetFillStyle(3005);
        RotDefHisy->SetMarkerStyle(34);
        RotDefHis->SetTitle("");
        RotDefHisy->Draw("Pa3");
        RotDefHisy->Write();
        RotDefHisy->GetXaxis()->SetTitle("Angle around Z axis (#circ)");
        RotDefHisy->GetYaxis()->SetTitle("RMS of scattering distribution");
        RotDefHisy->GetYaxis()->SetTitleOffset(1.5);
        c1->SaveAs("RotDefHisy.pdf");
        c1->Clear();
        ugraph->Draw("AC*");
        c1->SaveAs("ugraph.pdf");
        c1->Clear();
        sgraph->Draw("AC*");
        c1->SaveAs("sgraph.pdf");
        c1->Clear();
        vgraph->Draw("AC*");
        c1->SaveAs("vgraph.pdf");
        c1->SaveAs("vgraph.root");
        c1->Clear();
        wgraph->Draw("AC*");
        c1->SaveAs("wgraph.pdf");
        c1->SaveAs("wgraph.root");
        c1->Clear();
        ddotwgraph->Draw("AC*");
        c1->SaveAs("ddotwgraph.pdf");
        c1->Clear();
        ddotugraph->Draw("AC*");
        c1->SaveAs("ddotugraph.pdf");
        c1->Clear();
        ddotwdotugraph->Draw("AC*");
        c1->SaveAs("ddotwdotugraph.pdf");
        c1->Clear();
        c1->SetLogy();
        scattering_proj_x->Draw();
        c1->Print("scattering_proj_x.pdf");
        c1->Clear();
        scattering_proj_y->Draw();
        c1->Print("scattering_proj_y.pdf");
        /*
           c1->Clear();
           c1->SetLogy();
           gStyle->SetOptStat(0);
           convergence->Draw();
        //convergence->GetXaxis()->SetTitle("No. of iterations");
        //convergence->GetYaxis()->SetTitle("#chi^{2}_{n}-#chi^{2}_{n-1}");
        c1->SaveAs("convergence.pdf");
        c1->SaveAs("convergence.root");
        */

        delete c1;
        delete thetaX_measdata;
        delete thetaY_measdata;
        delete thetaScatt_measdata;
        delete thetaScatt_measdata_vp;
        delete thetaX_reco;
        delete thetaY_reco;
        //delete thetaX_reco_noeffi;
        //delete thetaY_reco_noeffi;
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
        //delete convergence;

    }

    void MCSAnalysis::DoDeconvolution(std::string model, int mode, int n_sel=1){
        // The basic methods associated with the RooUnfolding package
        // First generate a histogram of the measured data.

        int isfirst = 0;
        if (model.find(modelname1.c_str()) != std::string::npos)
            isfirst = 1;

        // int n_base = int(2 * _DSset.N() / (n_sel * (n_sel + 1)));
        // int k=0;
        // for (int j=0; j<n_sel; j++)
        // define the number of events under the current selection
        int curr_sel = int(_DSset.N());  // n_base * j;
        int curr_k = 0;



        std::string  tmpname = "thetaX_measdata";
        tmpname += model;
        TH1D* thetaX_measdata = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaX_asymm";
        tmpname += model;
        TH1* thetaX_asymm = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{XY}}; Events",_histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);
        tmpname = "thetaX_measdataminus";
        tmpname += model;
        TH1D* thetaX_measdataminus = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_measdata";
        tmpname += model;
        TH1D* thetaY_measdata = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaY_measdatafold";
        tmpname += model;
        TH1D* thetaY_measdatafold = 
            new TH1D(tmpname.c_str(),";#it{#theta_{Y}};#it{#theta_{Yi}}-#it{#theta_{Y46-i}}",_histlimits["NbinsXY"]/2, _histlimits["minXY"], 0);
        tmpname = "thetaY_asymm";
        tmpname += model;
        TH1* thetaY_asymm = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);
        tmpname = "thetaY_measdataminus";
        tmpname += model;
        TH1D* thetaY_measdataminus = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaScatt_measdata";
        tmpname += model;
        TH1D* thetaScatt_measdata = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        tmpname = "theta2Scatt_measdata";
        tmpname += model;
        TH1D* theta2Scatt_measdata = 
            new TH1D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;#theta^{2}_{Scatt}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);
        tmpname = "thetaScatt_measdata_vp";
        tmpname += model;
        TH2D* thetaScatt_measdata_vp = 
            new TH2D(tmpname.c_str(),"Scattering Angle between Momentum Vectors;Momentum (MeV/c); #it{#theta_{Scatt}}", 
                    400, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        tmpname = "thetaXUS_thetaXDS";
        tmpname += model;

        TH2D* thetaXUS_thetaXDS = 
            new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"], _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmpname = "thetaYUS_thetaYDS";
        tmpname += model;
        TH2D* thetaYUS_thetaYDS = 
            new TH2D(tmpname.c_str(),"Upstream vs. Downstream Angle;#it{#theta_{X}^{US}}; #it{#theta_{X}^{DS}}",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"],_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmpname = "projposUSDSdiff_";
        tmpname += model;
        TH2D* projposUSDSdiff = 
            new TH2D(tmpname.c_str(),
                    "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
                    180, -400, 400, 180, -400, 400);
        tmpname = "asymmXYUS_";
        tmpname += model;
        TH2D* asymmXYUS = 
            new TH2D(tmpname.c_str(),";#it{x}(mm);#it{y}(mm)",100, -200, 200, 100, -200, 200);
        tmpname = "asymmXYDS_";
        tmpname += model;
        TH2D* asymmXYDS = 
            new TH2D(tmpname.c_str(),";#it{x}(mm);#it{y}(mm)",100, -200, 200, 100, -200, 200);

        //1D Histogram Settings
        TH1D *defineHist2(const char* name, const char* title, Double_t nbinsx, Double_t xlow, Double_t xup);
        TH1D *Histo[180];
        TH1D *Histoy[180];
        TH1D *Histodotw[180];
        TH1D *Histodotu[180];
        TH1D *Histoyd[180];
        TH1D *Histowd[180];
        TH1D *Histowx[180];
        TH1D *Histowy[180];
        TH1D *Histou[180];
        TH1D *Histos[180];
        TH1D *Histov[180];
        TH1D *Histow[180];
        TH1D *Histodotwdotu[180];

        std::cout << "_histlimits[NbinsXY] " << _histlimits["NbinsXY"] << std::endl;
        std::cout << "_histlimits[minXY] " << _histlimits["minXY"] << std::endl;
        std::cout << "_histlimits[maxXY] " << _histlimits["maxXY"] << std::endl;
        for (int l=0;l<180;l++){
            std::string  tmpname = "thetaX_measdata_";
            tmpname += l;
            Histo[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            /* Histo[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"],_histlimits["minXY"], _histlimits["maxXY"]); */
            //Histo[l] = defineHist2(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "thetaY_measdata_";
            tmpname += l;
            Histoy[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),_histlimits["NbinsXY"],_histlimits["minXY"], _histlimits["maxXY"]);
            /*     tmpname = "dotw_"; */
            /*     tmpname += l; */
            /*     Histodotw[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -0.04, 0.04); */
            /*     tmpname = "dotu_"; */
            /*     tmpname += l; */
            /*     Histodotu[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -0.04, 0.04); */
            /*     tmpname = "yd_"; */
            /*     tmpname += l; */
            /*     Histoyd[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, 1, 3); */
            /*     tmpname = "wd_"; */
            /*     tmpname += l; */
            /*     Histowd[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, 1, 3); */
            /*     tmpname = "wx_"; */
            /*     tmpname += l; */
            /*     Histowx[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, -0.04, 0.04); */
            /*     tmpname = "wy_"; */
            /*     tmpname += l; */
            /*     Histowy[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),100, -0.04, 0.04); */
            /*     tmpname = "u"; */
            /*     tmpname += l; */
            /*     Histou[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -0.04, 0.04); */
            /*     tmpname = "s"; */
            /*     tmpname += l; */
            /*     Histos[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -0.04, 0.04); */
            /*     tmpname = "dotv"; */
            /*     tmpname += l; */
            /*     Histov[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -0.04, 0.04); */
            /*     tmpname = "w"; */
            /*     tmpname += l; */
            /*     Histow[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),200, -0.04, 0.04); */
            /*     tmpname = "dotwdotu_"; */
            /*     tmpname += l; */
            /*     Histodotwdotu[l] = defineHist2(tmpname.c_str(),tmpname.c_str(),46, -0.04, 0.04); */

        }


        const Int_t NBINS = 19;
        Double_t scat_bin_array[NBINS + 1] = {-0.1151,-0.0938,-0.0754,-0.0597,-0.0463,-0.0347,-0.0248,-0.0162,-0.00895,-0.00269,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,0.0754,0.0938,0.1151};
        TH1D* scattering_proj_x = new TH1D("scattering_proj_x_DC","Change in Projected Angle (X);#it{#theta_{X}}; Events per radian", 
                NBINS, scat_bin_array);
        TH1D* scattering_proj_y = new TH1D("scattering_proj_y_DC","Change in Projected Angle (Y);#it{#theta_{Y}}; Events per radian", 
                NBINS, scat_bin_array);


        for(int i=curr_k; i<curr_sel; i++){
            //double hold = _USset.E(i).dXdz*cos(10)+sin(10);
            //_USset.E(i).dXdz = hold;
            std::vector<double> projTheta = RotDefineProjectionAngles(_USset.E(i), _DSset.E(i), angdef);
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
            thetaY_measdataminus->Fill(-thetaY);
            thetaX_measdataminus->Fill(-thetaX);
            if (thetaY>0.02) {
                asymmXYUS->Fill(_USset.E(i).X, _USset.E(i).Y);
                asymmXYDS->Fill(_DSset.E(i).X, _DSset.E(i).Y);
            }
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
            for (int l=0;l<180;l++){
                //if (_DSset.E(i).dXdz<0.2 && _DSset.E(i).dYdz<0.2) {
                std::vector<double> RotprojTheta = RotDefineProjectionAngles(_USset.E(i), _DSset.E(i),l);
                if (mode==-7) {
                    if (_DSset.E(i).jDS == -1 || _DSset.E(i).kDS == -1){
                        continue;
                    }
                }
                double RotthetaY = RotprojTheta[1];
                double RotthetaX = RotprojTheta[0];
                /* //std::cout << "RotprojTheta[1] " << RotprojTheta[1] << std::endl; */
                Histo[l]->Fill(RotthetaX);
                Histoy[l]->Fill(RotthetaY);
                /*         //if (RotprojTheta[3]>0.2) Histodotw[l]->Fill(RotprojTheta[3]); */ 
                /*         //if (RotprojTheta[3]<-0.2) Histodotw[l]->Fill(RotprojTheta[3]); */ 
                /*         //if (RotprojTheta[4]>0.2) Histodotv[l]->Fill(RotprojTheta[4]); */ 
                /*         //if (RotprojTheta[4]<-0.2) Histodotv[l]->Fill(RotprojTheta[4]); */ 
                /*         Histodotw[l]->Fill(RotprojTheta[3]); */ 
                /*         //Histodotw[l]->Fill(RotprojTheta[3]); */ 
                /*         Histodotu[l]->Fill(RotprojTheta[4]); */ 
                /*         Histoyd[l]->Fill(RotprojTheta[5]); */ 
                /*         Histowd[l]->Fill(RotprojTheta[6]); */ 
                /*         //if (i==0) { */ 
                /* Histowx[l]->Fill(RotprojTheta[8]); */ 
                /*         Histowy[l]->Fill(RotprojTheta[9]); */ 
                /* Histou[l]->Fill(RotprojTheta[10]); */
                /* Histos[l]->Fill(RotprojTheta[11]); */
                /* //} */
                /* Histov[l]->Fill(RotprojTheta[12]); */
                /* Histow[l]->Fill(RotprojTheta[13]); */
                /*         Histodotwdotu[l]->Fill(RotprojTheta[14]); */ 
                /*         //Histov_x[l]->Fill(RotprojTheta[8[0]]); */ 
                /*         //Histov_y[l]->Fill(RotprojTheta[8[0]]); */ 
                /* //	    } */
            }
        }

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

        /*
           TH1D* thetaX_reco_n; 
           TH1D* thetaX_measdata_n; 
           if (model.find(modelname3.c_str()) != std::string::npos) {
           double chi2[100];
           double chi2effi[100];
           double xaxis[100];
           for ( int itr=1; itr<100; itr++){    
           RooUnfoldBayes unfold_thetaX_n(&fresp_thetaX, thetaX_measdata, itr);
           thetaX_reco_n = (TH1D*)unfold_thetaX_n.Hreco();
           chi2[itr]=thetaX_reco_n->Chi2Test(theta_true_x_graph,"WUP");
           xaxis[itr] = itr;
           thetaX_measdata_n = trkreffix(thetaX_measdata);
           RooUnfoldBayes unfold_thetaX_eff(&fresp_thetaX, thetaX_measdata_n, itr);
           thetaX_reco_n = (TH1D*)unfold_thetaX_eff.Hreco();
           chi2effi[itr]=thetaX_reco_n->Chi2Test(theta_true_x_graph,"WUP");
           }
           TCanvas* c2 = new TCanvas("c2","A Simple Graph Example",200,10,1400,1000);
           TGraph* convergence = new TGraph(100,xaxis,chi2);
           convergence->GetXaxis()->SetTitle("iterations");
           convergence->GetYaxis()->SetTitle("#chi^{2} MC reco & Truth");
           convergence->Draw();
           c2->SaveAs("convergence.pdf");
           c2->SaveAs("convergence.root");
           c2->Clear();
           TGraph* convergenceeffi = new TGraph(100,xaxis,chi2effi);
           convergenceeffi->GetXaxis()->SetTitle("iterations");
           convergenceeffi->GetYaxis()->SetTitle("#chi^{2} MC reco & Truth");
           convergenceeffi->Draw();
           c2->SaveAs("convergenceeffi.pdf");
           c2->SaveAs("convergenceeffi.root");
           delete c2;
           delete convergence;
           }
           */
        /*
           thetaX_measdata = trkreffix(thetaX_measdata);
           thetaY_measdata = trkreffiy(thetaY_measdata);
           thetaScatt_measdata = trkreffiscatt(thetaScatt_measdata);
           theta2Scatt_measdata = trkreffi2scatt(theta2Scatt_measdata);
           */

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
        thetaX_reco->SetTitle(";#Delta #it{#theta_{X}}");
        //thetaX_reco = trkreffix(thetaX_reco);
        TH1D* thetaX_reco_noeffi = (TH1D*)unfold_thetaX.Hreco();
        tmpname = "thetaX_reco_noeffi";
        tmpname += model;
        thetaX_reco_noeffi->SetName(tmpname.c_str());
        TH1D* thetaY_reco = (TH1D*)unfold_thetaY.Hreco();
        tmpname = "thetaY_reco";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaY_reco->SetName(tmpname.c_str());
        thetaY_reco->SetTitle(";#Delta #it{#theta_{Y}}");
        tmpname = "thetaY_reco_noeffi";
        TH1D* thetaY_reco_noeffi = (TH1D*)unfold_thetaY.Hreco();
        tmpname += model;
        thetaY_reco_noeffi->SetName(tmpname.c_str());
        //thetaY_reco = trkreffiy(thetaY_reco);
        TH1D* thetaScatt_reco = (TH1D*)unfold_thetaScatt.Hreco();
        tmpname = "thetaScatt_reco";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaScatt_reco->SetName(tmpname.c_str());
        thetaScatt_reco->SetTitle(";#it{#theta_{Scatt}}");
        //thetaScatt_reco = trkreffiscatt(thetaScatt_reco);
        TH1D* theta2Scatt_reco = (TH1D*)unfold_theta2Scatt.Hreco();
        tmpname = "theta2Scatt_reco";
        tmpname += model;
        // if(j>0) tmpname += j;
        theta2Scatt_reco->SetName(tmpname.c_str());
        theta2Scatt_reco->SetTitle(";#theta^{2}_{Scatt}");
        //theta2Scatt_reco = trkreffi2scatt(theta2Scatt_reco);

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
        thetaX_truth->SetTitle(";#Delta #it{#theta_{X}}");
        TH1D* thetaY_truth = (TH1D*)unfold_thetaY.response()->Htruth();
        tmpname = "thetaY_";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaY_truth->SetName(tmpname.c_str());
        thetaX_truth->SetTitle(";#Delta #it{#theta_{Y}}");
        TH1D* thetaScatt_truth = (TH1D*)unfold_thetaScatt.response()->Htruth();
        tmpname = "thetaScatt_";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaScatt_truth->SetName(tmpname.c_str());
        thetaScatt_truth->SetTitle(";#it{#theta_{Scatt}}");
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
        thetaX_measured->SetTitle(";#Delta #it{#theta_{X}}");
        TH1D* thetaY_measured = (TH1D*)unfold_thetaY.response()->Hmeasured();
        tmpname = "thetaY_measured";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaY_measured->SetName(tmpname.c_str());
        thetaY_measured->SetTitle(";#Delta #it{#theta_{Y}}");
        TH1D* thetaScatt_measured = (TH1D*)unfold_thetaScatt.response()->Hmeasured();
        tmpname = "thetaScatt_measured";
        tmpname += model;
        // if(j>0) tmpname += j;
        thetaScatt_measured->SetName(tmpname.c_str());
        thetaScatt_measured->SetTitle(";#it{#theta_{Scatt}}");
        TH1D* theta2Scatt_measured = (TH1D*)unfold_theta2Scatt.response()->Hmeasured();
        tmpname = "theta2Scatt_measured";
        tmpname += model;
        // if(j>0) tmpname += j;
        theta2Scatt_measured->SetName(tmpname.c_str());
        theta2Scatt_measured->SetTitle(";#theta^{2}_{Scatt}");


        TH2D* dotw = new TH2D("dotw","dotw", 180, 0, 180, 200, -1, 1);
        TH2D* dotv = new TH2D("dotv","dotv", 180, 0, 180, 200, -1, 1);
        TH2D* histoyd = new TH2D("histoyd","RMS with plane definition", 180, 0, 180, 100, 1, 3);
        TH2D* histowd = new TH2D("histowd","RMS with plane definition", 180, 0, 180, 100, 1, 3);
        TH2D* histowx = new TH2D("histowx","RMS with plane definition", 180, 0, 180, 100, -2, 2);
        TH2D* histowy = new TH2D("histowy","RMS with plane definition", 180, 0, 180, 100, -2, 2);

        outfile->cd();
        double x[180];
        double xer[180];
        double y[180];
        double yer[180];
        double xaxis[180];
        double xaxiser[180];
        /* double u[180]; */
        /* double s[180]; */
        /* double zerox[180]; */
        /* double onex[180]; */
        /* double v[180]; */
        /* double w[180]; */
        /* double ddotw[180]; */
        /* double ddotu[180]; */
        /* double ddotwdotu[180]; */
        for (int l=0;l<180;l++){
            //std::cout << "Histo[l]->GetRMS() " << Histo[l]->GetRMS() << std::endl;
            //std::cout << "Histoy[l]->GetRMS() " << Histoy[l]->GetRMS() << std::endl;
            xaxis[l]=l;
            xaxiser[l]=0;
            Histo[l]->Fit("gaus", "","",-0.03,0.03);
            TF1 *myfunc = Histo[l]->GetFunction("gaus");
            x[l]=myfunc->GetParameter(2);
            x[l]=Histo[l]->GetRMS();
            xer[l]=0.00054;
            y[l]=Histoy[l]->GetRMS();
            yer[l]=0.0006;
            /*    u[l]=Histou[l]->GetRMS(); */
            /*    s[l]=Histos[l]->GetRMS(); */
            /*    v[l]=Histov[l]->GetRMS(); */
            /*    w[l]=Histow[l]->GetRMS(); */
            /*    ddotw[l]=Histodotw[l]->GetRMS(); */
            /*    ddotu[l]=Histodotu[l]->GetRMS(); */
            /*    zerox[l]=Histowx[l]->GetRMS(); */
            /*    onex[l]=Histowy[l]->GetRMS(); */
            /*    histoyd->Fill(l,Histoyd[l]->GetMean()); */
            /*    histowd->Fill(l,Histowd[l]->GetMean()); */
            /*    histowx->Fill(l,Histowx[l]->GetRMS()); */
            /*    histowy->Fill(l,Histowy[l]->GetRMS()); */
            /*    ddotwdotu[l]=Histodotwdotu[l]->GetRMS(); */
            /*    //Histowx[l]->Write(); */
            /*    //Histowy[l]->Write(); */
            /*    //Histoyd[l]->Write(); */ 
            /*    //Histowd[l]->Write(); */ 
            /*    if (l>38 && l<45) { */
            /*    Histo[l]->Write(); */
            /*    //Histoy[l]->Write(); */
            /*    } */
            /*    //Histodotw[l]->Write(); */ 
            /*    //Histodotu[l]->Write(); */
            /*    //Histodotwdotu[l]->Write(); */ 
            /*    //Histodotv[l]->Write(); */ 
        }
        //histowx->Write();
        //histowy->Write();
        //histoyd->Write();
        //histowd->Write();
        //dotw->Write();
        //dotv->Write();
        thetaXUS_thetaXDS->Write();
        thetaYUS_thetaYDS->Write();
        //thetaX_measdata = trkreffix(thetaX_measdata);
        outfile->cd();
        thetaX_measdata->Write();
        //thetaY_measdata = trkreffiy(thetaY_measdata);
        outfile->cd();
        thetaY_measdata->Write();
        for (int i=1;i<24;i++) {
            thetaY_measdatafold->SetBinContent(i,thetaY_measdata->GetBinContent(i)-thetaY_measdata->GetBinContent(46-i));
        }
        thetaY_measdatafold->Write();
        //thetaScatt_measdata = trkreffiscatt(thetaScatt_measdata);
        outfile->cd();
        thetaScatt_measdata->Write();
        //theta2Scatt_measdata = trkreffi2scatt(theta2Scatt_measdata);
        outfile->cd();
        theta2Scatt_measdata->Write();
        thetaScatt_measdata_vp->Write();
        thetaX_reco->Write();
        thetaY_reco->Write();
        thetaX_reco_noeffi->Write();
        thetaY_reco_noeffi->Write();
        thetaScatt_reco->Write();
        theta2Scatt_reco->Write();
        thetaX_response->Write();
        thetaY_response->Write();
        thetaScatt_response->Write();
        theta2Scatt_response->Write();
        //thetaX_truth->Write();
        //thetaY_truth->Write();
        //thetaScatt_truth->Write();
        //theta2Scatt_truth->Write();
        thetaX_measured->Write();
        thetaY_measured->Write();
        thetaScatt_measured->Write();
        theta2Scatt_measured->Write();
        scattering_proj_x->Write();
        scattering_proj_y->Write();
        projposUSDSdiff->Write();
        //convergence->Write();

        TCanvas* c1 = new TCanvas();
        thetaY_asymm=thetaY_measdata->GetAsymmetry(thetaY_measdataminus);
        thetaY_asymm->Draw();
        TF1*  fM2D = new TF1("fM2D", "[0]+x*[1]", 0, 0.02);
        //fM2D->SetParameter(0,0);
        //fM2D->SetParameter(1,1);
        //fM2D->SetParNames("A","B");
        //thetaY_asymm->Fit("pol1","RES","",0,0.02);
        thetaY_asymm->Fit(fM2D,"RES","",0,0.02);
        //fM2D->Draw("SAME");
        thetaY_asymm->SetName("thetaY_asymm");
        thetaY_asymm->SetTitle("#it{#theta_Y} Asymmetry");
        thetaY_asymm->GetYaxis()->SetTitle("Asymmetry");
        thetaY_asymm->GetXaxis()->SetRange(0,0.06);
        thetaY_asymm->Write();
        c1->SaveAs("thetaY_asymm.pdf");
        c1->Clear();
        thetaX_asymm=thetaX_measdata->GetAsymmetry(thetaX_measdataminus);
        thetaX_asymm->Draw();
        //thetaX_asymm->Fit("pol1","RES","",0,0.02);
        thetaX_asymm->Fit(fM2D,"RES","",0,0.02);
        //fM2D->Draw("SAME");
        thetaX_asymm->SetName("thetaX_asymm");
        thetaX_asymm->SetTitle("#it{#theta_X} Asymmetry");
        thetaX_asymm->GetYaxis()->SetTitle("Asymmetry");
        thetaX_asymm->GetXaxis()->SetRange(0,0.06);
        thetaX_asymm->Write();
        c1->SaveAs("thetaX_asymm.pdf");
        c1->Clear();
        asymmXYUS->Draw("colz");
        c1->Print("asymmXYUS.pdf");
        c1->Clear();
        asymmXYDS->Draw("colz");
        c1->Print("asymmXYDS.pdf");
        c1->Clear();
        thetaX_measdata->Draw();
        c1->Print("thetaX_measdata.pdf");
        c1->Clear();
        thetaY_measdata->Draw();
        c1->Print("thetaY_measdata.pdf");
        c1->Print("thetaY_measdata.root");
        c1->Clear();
        thetaY_measdatafold->Draw();
        c1->Print("thetaY_measdatafold.pdf");
        c1->Print("thetaY_measdatafold.root");
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
        thetaX_reco_noeffi->Draw();
        c1->Print("thetaX_reco_noeffi.pdf");
        c1->Clear();
        thetaY_reco_noeffi->Draw();
        c1->Print("thetaY_reco_noeffi.pdf");
        c1->Clear();
        thetaScatt_reco->Draw();
        c1->Print("thetaScatt_reco.pdf");
        c1->Clear();
        thetaX_response->Draw("colz");
        c1->Print("thetaX_response.pdf");
        c1->Clear();
        thetaY_response->Draw("colz");
        c1->Print("thetaY_response.pdf");
        c1->Clear();
        thetaScatt_response->Draw("colz");
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
        TGraphErrors* RotDefHis = new TGraphErrors(180, xaxis,x,xaxiser,xer);
        TGraphErrors* RotDefHisy = new TGraphErrors(180, xaxis,y,xaxiser,yer);
        /* TGraph* ugraph = new TGraph(180, xaxis,u); */
        /* TGraph* sgraph = new TGraph(180, xaxis,s); */
        /* TGraph* vgraph = new TGraph(180, xaxis,zerox); */
        /* TGraph* wgraph = new TGraph(180, xaxis,onex); */
        /* TGraph* ddotwgraph = new TGraph(180, xaxis,ddotw); */
        /* TGraph* ddotugraph = new TGraph(180, xaxis,ddotu); */
        /* TGraph* ddotwdotugraph = new TGraph(180, xaxis,ddotwdotu); */
        RotDefHis->SetFillColor(6);
        RotDefHis->SetFillStyle(3005);
        RotDefHis->SetMarkerStyle(34);
        RotDefHis->SetTitle("");
        RotDefHis->Draw("AC*");
        RotDefHis->Write();
        RotDefHis->GetXaxis()->SetTitle("Angle around Z axis (#circ)");
        RotDefHis->GetYaxis()->SetTitle("RMS of scattering distribution");
        c1->SaveAs("RotDefHis.pdf");
        c1->Clear();
        RotDefHisy->SetFillColor(6);
        RotDefHisy->SetFillStyle(3005);
        RotDefHisy->SetMarkerStyle(34);
        RotDefHis->SetTitle("");
        RotDefHisy->Draw("Pa3");
        RotDefHisy->Write();
        RotDefHisy->GetXaxis()->SetTitle("Angle around Z axis (#circ)");
        RotDefHisy->GetYaxis()->SetTitle("RMS of scattering distribution");
        RotDefHisy->GetYaxis()->SetTitleOffset(1.5);
        c1->SaveAs("RotDefHisy.pdf");
        /* c1->Clear(); */
        /* ugraph->Draw("AC*"); */
        /* c1->SaveAs("ugraph.pdf"); */
        /* c1->Clear(); */
        /* sgraph->Draw("AC*"); */
        /* c1->SaveAs("sgraph.pdf"); */
        /* c1->Clear(); */
        /* vgraph->Draw("AC*"); */
        /* c1->SaveAs("vgraph.pdf"); */
        /* c1->SaveAs("vgraph.root"); */
        /* c1->Clear(); */
        /* wgraph->Draw("AC*"); */
        /* c1->SaveAs("wgraph.pdf"); */
        /* c1->SaveAs("wgraph.root"); */
        /* c1->Clear(); */
        /* ddotwgraph->Draw("AC*"); */
        /* c1->SaveAs("ddotwgraph.pdf"); */
        /* c1->Clear(); */
        /* ddotugraph->Draw("AC*"); */
        /* c1->SaveAs("ddotugraph.pdf"); */
        /* c1->Clear(); */
        /* ddotwdotugraph->Draw("AC*"); */
        /* c1->SaveAs("ddotwdotugraph.pdf"); */
        c1->Clear();
        c1->SetLogy();
        scattering_proj_x->Draw();
        c1->Print("scattering_proj_x.pdf");
        c1->Clear();
        scattering_proj_y->Draw();
        c1->Print("scattering_proj_y.pdf");
        /*
           c1->Clear();
           c1->SetLogy();
           gStyle->SetOptStat(0);
           convergence->Draw();
        //convergence->GetXaxis()->SetTitle("No. of iterations");
        //convergence->GetYaxis()->SetTitle("#chi^{2}_{n}-#chi^{2}_{n-1}");
        c1->SaveAs("convergence.pdf");
        c1->SaveAs("convergence.root");
        */

        delete c1;
        delete thetaX_measdata;
        delete thetaY_measdata;
        delete thetaScatt_measdata;
        delete thetaScatt_measdata_vp;
        delete thetaX_reco;
        delete thetaY_reco;
        //delete thetaX_reco_noeffi;
        //delete thetaY_reco_noeffi;
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
        //delete convergence;

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
        //hx->Rebin(GetNbinsX()/40);
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
            new TH1D("thetaX_data","Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        TH1D* thetaY_data = 
            new TH1D("thetaY_data","Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        TH1D* thetaScatt_data = 
            new TH1D("thetaScatt_data","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        TH1D* theta2Scatt_data = 
            new TH1D("theta2Scatt_data","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}^{2}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh2"], _histlimits["maxTh2"]);

        TH1D* thetaX_ref = 
            new TH1D("thetaX_ref","Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        TH1D* thetaY_ref = 
            new TH1D("thetaY_ref","Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        TH1D* thetaScatt_ref = 
            new TH1D("thetaScatt_ref","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);
        TH1D* theta2Scatt_ref = 
            new TH1D("theta2Scatt_ref","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}^{2}}; Events per mrad",_histlimits["NbinsTh2"], _histlimits["minTh2"], _histlimits["maxTh2"]);

        /*
           TH1D* thetaX_fft = 
           new TH1D("thetaX_fft","Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
           TH1D* thetaY_fft = 
           new TH1D("thetaY_fft","Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
           TH1D* thetaScatt_fft = 
           new TH1D("thetaScatt_fft","Scattering Angle between Momentum Vectors;#it{#theta_{Scatt}}; Events per mrad",_histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);


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
            new TH2D("thetaScatt_measdata_vp","Scattering Angle between Momentum Vectors;Momentum (MeV/c); #it{#theta_{Scatt}}", 
                    400, 100, 300, _histlimits["NbinsTh"], _histlimits["minTh"], _histlimits["maxTh"]);

        const Int_t NBINS = 19;
        Double_t scat_bin_array[NBINS + 1] = {-0.1151,-0.0938,-0.0754,-0.0597,-0.0463,-0.0347,-0.0248,-0.0162,-0.00895,-0.00269,0.00269,0.00895,0.0162,0.0248,0.0347,0.0463,0.0597,0.0754,0.0938,0.1151};
        TH1D* scattering_proj_x = new TH1D("scattering_proj_x_DC","Change in Projected Angle (X);#it{#theta_{X}}; Events per radian", 
                NBINS, scat_bin_array);
        TH1D* scattering_proj_y = new TH1D("scattering_proj_y_DC","Change in Projected Angle (Y);#it{#theta_{Y}}; Events per radian", 
                NBINS, scat_bin_array);

        TH2D* projposUSDSdiff = 
            new TH2D("projposUSDSdiff_ref",
                    "Difference of US and DS projections at absorber;#Delta x_{DS-US} (mm); #Delta y_{DS-US}",
                    180, -400, 400, 180, -400, 400);

        for(int i=0; i<_DSset.N(); i++){
            std::vector<double> projTheta = RotDefineProjectionAngles(_USset.E(i), _DSset.E(i), angdef);
            double thetaY = projTheta[1];  // atan(_DSset.E(i).dXdz) - atan(_USset.E(i).dXdz);
            double thetaX = projTheta[0];  // atan(_DSset.E(i).dYdz) - atan(_USset.E(i).dYdz);
            // double cosScatt = ( (1 + _USset.E(i).dXdz * _DSset.E(i).dXdz +
            //_USset.E(i).dYdz * _DSset.E(i).dYdz )/
            //			sqrt(1 + _USset.E(i).dXdz*_USset.E(i).dXdz +
            //			      _USset.E(i).dYdz*_USset.E(i).dYdz)/
            //			sqrt(1 + _DSset.E(i).dXdz*_DSset.E(i).dXdz +
            //			     _DSset.E(i).dYdz*_DSset.E(i).dYdz));
            double thetaScatt = projTheta[2];  /// acos(cosScatt);

            if (mode==-2) { 
                if (_USset.E(i).pid!=-13) continue;
            }
            thetaX_data->Fill(thetaX);
            thetaY_data->Fill(thetaY);
            /* if (thetaY<0.00127 && thetaY>-0.00127) std::cout << "_USset.E(i).dXdz " << _USset.E(i).dXdz << std::endl; */
            /* if (thetaY<0.00127 && thetaY>-0.00127) std::cout << "_DSset.E(i).dXdz " << _DSset.E(i).dXdz << std::endl; */
            /* if (thetaY<0.00127 && thetaY>-0.00127) std::cout << "thetaY " << thetaY << std::endl; */
            /* if (thetaX<0.00127 && thetaX>-0.00127) std::cout << "_USset.E(i).dYdz " << _USset.E(i).dYdz << std::endl; */
            /* if (thetaX<0.00127 && thetaX>-0.00127) std::cout << "_DSset.E(i).dYdz " << _DSset.E(i).dYdz << std::endl; */
            /* if (thetaX<0.00127 && thetaX>-0.00127) std::cout << "thetaX " << thetaX << std::endl; */
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
        gmeanX->SetTitle(";Run Number;Mean #it{x} (mm)");
        TGraphErrors* gmeanY = new TGraphErrors(_run.size());
        gmeanY->SetName("gmeanY");
        gmeanY->SetTitle(";Run Number;Mean #it{y} (mm)");
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
            double USrefpos = _sys["abspos"] - 33.;
            double DSrefpos = _sys["abspos"] + 33.;
            for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
                double ztest = mcevent->GetVirtualHits()->at(j).GetPosition().z();
                double USzdiff = fabs(ztest - USrefpos);
                double DSzdiff = fabs(ztest - DSrefpos);
                /* std::cout<<j<<"\t"<<ztest<<"\t"<<USzdiff<<"\t"<<DSzdiff<<std::endl; */
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
            if (USmindiff < 5.){
                planesfound = true;
                USabsPlaneI = USminI;
                /* std::cout<<"US absorber plane.\t"<<USrefpos<<"\t"<<USmindiff<<"\t"<<USabsPlaneI<<"\t"<<mcevent->GetVirtualHits()->at(USabsPlaneI).GetPosition().z()<<std::endl; */
            }
            if (DSmindiff < 5. && planesfound==true){
                DSabsPlaneI = DSminI;
                /* std::cout<<"DS absorber plane.\t"<<DSrefpos<<"\t"<<DSmindiff<<"\t"<<DSabsPlaneI<<"\t"<<mcevent->GetVirtualHits()->at(DSabsPlaneI).GetPosition().z()<<std::endl; */
            }
            else 
                planesfound = false;

            // Now I will fill indicies corresponding to the tracker
            // reference planes
            USrefpos = _sys["USref"];
            DSrefpos = _sys["DSref"];

            USmindiff = 9999, DSmindiff=9999;
            USminI=-1, DSminI=-1;

            for ( size_t j=0; j < mcevent->GetVirtualHits()->size(); j++){
                double ztest = mcevent->GetVirtualHits()->at(j).GetPosition().z();
                double USzdiff = fabs(ztest - USrefpos);
                double DSzdiff = fabs(ztest - DSrefpos);
                /* std::cout<<j<<"\t"<<ztest<<"\t"<<USzdiff<<"\t"<<DSzdiff<<std::endl; */
                if(USzdiff < USmindiff){
                    USmindiff = USzdiff;
                    USminI = j;
                }

                if(DSzdiff < DSmindiff){
                    DSmindiff = DSzdiff;
                    DSminI = j;
                }
            }
            if (planesfound == true){
                if (USmindiff < 50. ){
                    USrefplaneI = USminI;
                    /* std::cout<<"US reference plane.\t"<<USrefplaneZ<<"\t"<<USmindiff<<"\t"<<USrefplaneI<<"\t"<<mcevent->GetVirtualHits()->at(USrefplaneI).GetPosition().z()<<std::endl; */
                }
                else  
                    planesfound = false;
                if (DSmindiff < 50. ){
                    DSrefplaneI = DSminI;
                    /* std::cout<<"DS reference plane.\t"<<DSrefplaneZ<<"\t"<<DSmindiff<<"\t"<<DSrefplaneI<<"\t"<<mcevent->GetVirtualHits()->at(DSrefplaneI).GetPosition().z()<<std::endl; */
                }
                else 
                    planesfound = false;
            }
            //std::cout << "planesfound " << planesfound << std::endl;
            return planesfound;
        }
    }

    void MCSAnalysis::FillMuScattResponse(bool& event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC){
        double thetaYMC = atan(DSMC.dXdz) - atan(USMC.dXdz);
        double thetaXMC = atan(DSMC.dYdz) - atan(USMC.dYdz);
        double XMC = DSMC.X;
        double YMC = DSMC.Y;
        //theta_true_x_graph->Fill(thetaXMC);
        //theta_true_y_graph->Fill(thetaYMC);
        if ( event_ok ){
            scattering_proj_y_R->Fill(XMC, thetaYMC);
            scattering_proj_x_R->Fill(YMC, thetaXMC);
            theta_meas_x_graph->Fill(thetaXMC);
            theta_meas_y_graph->Fill(thetaYMC);
        }
    }

    void MCSAnalysis::FillMCSResponse(bool& event_ok, Vars& US, Vars& DS, Vars& USMC, Vars& DSMC){

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

            //theta_true_x_graph->Fill(thetaXMC);
            //theta_true_y_graph->Fill(thetaYMC);
        }
        else{
            resp_thetaX.Miss(thetaXMC);
            resp_thetaY.Miss(thetaYMC);
            resp_thetaScatt.Miss(thetaScattMC);
            resp_theta2Scatt.Miss(thetaScattMC * thetaScattMC);

            //theta_true_x_graph->Fill(thetaXMC);
            //theta_true_y_graph->Fill(thetaYMC);
        }

    }

    void MCSAnalysis::TruthGraph(Collection& USMCTruthset, Collection& DSMCTruthset, bool recon){


        std::string tmpname;
        int curr_sel = int(DSMCTruthset.N());
        int curr_k = 0;
        if (!recon) {
            for(int i=curr_k; i<curr_sel; i++){
                std::vector<double> projDTheta =RotDefineProjectionAngles(USMCTruthset.E(i), DSMCTruthset.E(i), angdef); 
                //double thetaYMC = atan(DSMCTruthset.E(i).dXdz) - atan(USMCTruthset.E(i).dXdz);
                //double thetaXMC = atan(DSMCTruthset.E(i).dYdz) - atan(USMCTruthset.E(i).dYdz);
                double thetaYMC = projDTheta[1];
                double thetaXMC = projDTheta[0];
                double thetascattMC = projDTheta[2];

                //std::cout << "i " << i << std::endl;
                //std::cout << "projDTheta[0] " << projDTheta[0] << std::endl;
                theta_true_x_graph->Fill(thetaXMC);
                theta_true_y_graph->Fill(thetaYMC);
                minustheta_true_x_graph->Fill(-thetaXMC);
                minustheta_true_y_graph->Fill(-thetaYMC);
                theta_true_scat_graph->Fill(thetascattMC);
                theta_true_scat2_graph->Fill(thetascattMC*thetascattMC);
            }
            tmpname = "MCthetaX_asymm";
            TH1* MCthetaX_asymm = 
                new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{XY}}; Events",_histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);
            tmpname = "MCthetaY_asymm";
            TH1* MCthetaY_asymm = 
                new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events",_histlimits["NbinsXY"]/2, 0, _histlimits["maxXY"]);
            TCanvas* c1 = new TCanvas();
            MCthetaY_asymm=theta_true_y_graph->GetAsymmetry(minustheta_true_y_graph);
            MCthetaY_asymm->Draw();
            TF1*  fM2D = new TF1("fM2D", "[0]+x*[1]", 0, 0.02);
            //fM2D->SetParameter(0,0);
            //fM2D->SetParameter(1,1);
            //fM2D->SetParNames("A","B");
            //MCthetaY_asymm->Fit("pol1","RES","",0,0.02);
            MCthetaY_asymm->Fit(fM2D,"RES");
            //fM2D->Draw("SAME");
            MCthetaY_asymm->SetName("MCthetaY_asymm");
            MCthetaY_asymm->SetTitle("#it{#theta_Y} Asymmetry");
            MCthetaY_asymm->GetYaxis()->SetTitle("Asymmetry");
            MCthetaY_asymm->GetXaxis()->SetRange(0,0.06);
            outfile->cd();
            MCthetaY_asymm->Write();
            c1->SaveAs("MCthetaY_asymm.pdf");
            c1->Clear();
            MCthetaX_asymm=theta_true_x_graph->GetAsymmetry(minustheta_true_x_graph);
            MCthetaX_asymm->Draw();
            //MCthetaX_asymm->Fit("pol1","RES","",0,0.02);
            MCthetaX_asymm->Fit(fM2D,"RES");
            //fM2D->Draw("SAME");
            MCthetaX_asymm->SetName("MCthetaX_asymm");
            MCthetaX_asymm->SetTitle("#it{#theta_X} Asymmetry");
            MCthetaX_asymm->GetYaxis()->SetTitle("Asymmetry");
            MCthetaX_asymm->GetXaxis()->SetRange(0,0.06);
            MCthetaX_asymm->Write();
            c1->SaveAs("MCthetaX_asymm.pdf");
            TH1F *htruthmirrorl = new TH1F("htruthmirrorl","Mirror",_histlimits["NbinsXY"],-_histlimits["maxXY"],_histlimits["maxXY"]);
            TH1F *htruthmirrorr = new TH1F("htruthmirrorr","Mirror",_histlimits["NbinsXY"],-_histlimits["maxXY"],_histlimits["maxXY"]);
            for (int i = 1; i < _histlimits["NbinsXY"]/2+1; i++) htruthmirrorl->SetBinContent(i,theta_true_x_graph->GetBinContent(i));
            for (int i = 1; i < _histlimits["NbinsXY"]/2+1; i++) htruthmirrorr->SetBinContent(i,theta_true_x_graph->GetBinContent(_histlimits["NbinsXY"]+1-i));
            htruthmirrorl->Draw();
            htruthmirrorr->Draw();
            htruthmirrorl->Write();
            htruthmirrorr->Write();
            tmpname = "thetaX_graphsym";
            TH1* thetaX_graphsym = new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{XY}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            tmpname = "thetaY_graphsym";
            TH1* thetaY_graphsym = new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{XY}}; Events",_histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
            for (int i = 1; i < _histlimits["NbinsXY"]+1; i++) thetaX_graphsym->SetBinContent(i,(theta_true_x_graph->GetBinContent(i)+theta_true_x_graph->GetBinContent(_histlimits["NbinsXY"]+1-i))/2);
            for (int i = 1; i < _histlimits["NbinsXY"]+1; i++) thetaY_graphsym->SetBinContent(i,(theta_true_y_graph->GetBinContent(i)+theta_true_y_graph->GetBinContent(_histlimits["NbinsXY"]+1-i))/2);
            thetaX_graphsym->Draw();
            thetaY_graphsym->Draw();
            thetaX_graphsym->Write();
            thetaY_graphsym->Write();
                theta_true_x_graph->Write();
                theta_true_y_graph->Write();
        }
        else { 
            for(int i=curr_k; i<curr_sel; i++){
                std::vector<double> projDTheta =mcRotDefineProjectionAngles(USMCTruthset.E(i), DSMCTruthset.E(i), angdef); 
                double thetaYMC = projDTheta[1];
                double thetaXMC = projDTheta[0];
                double thetascattMC = projDTheta[2];

                if (thetaYMC!=0 || thetaXMC!=0 || thetascattMC!=0){
                    mctheta_true_x_graph->Fill(thetaXMC);
                    mctheta_true_y_graph->Fill(thetaYMC);
                    mctheta_true_scat_graph->Fill(thetascattMC);

                }
            }
        }
    }

    void MCSAnalysis::FillVarsVirtual(Vars &tmpvar, int& j, int& centre){

        if (mcevent->GetVirtualHits()->size()>centre && mcevent->GetVirtualHits()->size()>(j+1)) {
            tmpvar.X  = mcevent->GetVirtualHits()->at(j).GetPosition().x();
            tmpvar.Y  = mcevent->GetVirtualHits()->at(j).GetPosition().y();
            tmpvar.Z  = mcevent->GetVirtualHits()->at(j).GetPosition().z();
            tmpvar.dXdz = mcevent->GetVirtualHits()->at(j).GetMomentum().x()/
                mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            tmpvar.dYdz = mcevent->GetVirtualHits()->at(j).GetMomentum().y()/
                mcevent->GetVirtualHits()->at(j).GetMomentum().z();
            tmpvar.px = mcevent->GetVirtualHits()->at(j).GetMomentum().x();
            tmpvar.py = mcevent->GetVirtualHits()->at(j).GetMomentum().y();
            double pz = mcevent->GetVirtualHits()->at(centre).GetMomentum().z();
            tmpvar.TOF12 = 27.5 * sqrt(1 + 105.65*105.65/pz/pz);
            tmpvar.TOF01 = 26.5 * sqrt(1 + 105.65*105.65/pz/pz);
            tmpvar.pz = pz;
            tmpvar.pid = mcevent->GetVirtualHits()->at(j).GetParticleId();
            tmpvar.isgood = true;
        }
        else {
            tmpvar.X  = -9999;
            tmpvar.Y  = -9999;
            tmpvar.Z  = -9999;
            tmpvar.pz = -9999;
            float pz = -9999;
            double v1 = rand() % 100;
            if (v1/100 >0.5) {
                tmpvar.dXdz = -1./2.;
                tmpvar.dYdz = -1./2.;
                tmpvar.px   = -pz/2.;
                tmpvar.py   = -pz/2.;
                tmpvar.pz   = pz;
            }
            else {
                tmpvar.dXdz = 1./2.;
                tmpvar.dYdz = 1./2.;
                tmpvar.mcdXdz = 1./2.;
                tmpvar.mcdYdz = 1./2.;
                tmpvar.px   = pz/2.;
                tmpvar.py   = pz/2.;
                tmpvar.pz   = pz;
            }
            tmpvar.TOF12 = -9999;
            tmpvar.TOF01 = -9999;
            tmpvar.pid = -9999;
            tmpvar.isgood = false;
        }
    }

    void MCSAnalysis::FillEventCollection(Collection& Set, Vars& tmpvar){      
        Set.append_instance(tmpvar);
    }

    void MCSAnalysis::FillCollectionSciFi(Collection& Set, int& j, int& k, int& l, int& pid, double& pz, int& isDS, Vars& projdif, Vars& proj, std::vector<float> chi2nUS, int& beamtype, double& ptruth, const bool& difcut, const bool& fidcut, const bool& chicut, const bool& TOFcut){      

        if(j < int(scifievent->scifitracks().size()) && j != -1){
            if(k < int(scifievent->scifitracks()[j]->scifitrackpoints().size()) && k != -1){
                Vars tmpvar;
                FillVarsSciFi(tmpvar, j, k, l, pid, pz, isDS, proj, projdif, chi2nUS, beamtype, ptruth, difcut, fidcut, chicut, TOFcut);
                /*
                   if (project){
                   std::cout << "2" << std::endl;
                   Vars newvar = PropagateVarsMu(tmpvar, _sys["abspos"]);
                   tmpvar = newvar;
                   }
                   */
                Set.append_instance(tmpvar);

            }
        }
        if( j == -1 || k == -1){
            /* std::cout << " j & k are -1" << std::endl; */
            Vars tmpvar;
            tmpvar.X = -9999;
            tmpvar.Y = -9999;
            tmpvar.Z = -9999;
            double v1 = rand() % 100;
            if (v1/100 >0.5) {
                tmpvar.dXdz = -1./2.;
                tmpvar.dYdz = -1./2.;
                tmpvar.mcdXdz = -1./2.;
                tmpvar.mcdYdz = -1./2.;
                tmpvar.px   = -pz/2.;
                tmpvar.py   = -pz/2.;
                tmpvar.pz   = pz;
            }
            else {
                tmpvar.dXdz = 1./2.;
                tmpvar.dYdz = 1./2.;
                tmpvar.mcdXdz = 1./2.;
                tmpvar.mcdYdz = 1./2.;
                tmpvar.px   = pz/2.;
                tmpvar.py   = pz/2.;
                tmpvar.pz   = pz;
            }
            tmpvar.pid  = -13;
            if (pid!=-13) tmpvar.pid = pid;
            tmpvar.isgood = false;

            if( int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) == 1 &&
                    int(tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size()) == 1)
                tmpvar.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
                    - tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
            else
                tmpvar.TOF12 = 100.0 * 8.22475 / 0.299792458;
            if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) == 1 &&
                    int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) == 1 ) {
                tmpvar.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
                    - tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
            }
            else
                tmpvar.TOF01 = 100.0 * 7.64186 / 0.299792458;

            if(isDS==0){
                _sys["trkrlength"] = 1099.9;
                tmpvar.jUS = j;
                tmpvar.kUS = k;
                tmpvar.lUS = l;
            } else {
                _sys["trkrlength"] = 1099.76;
                tmpvar.jDS = j;
                tmpvar.kDS = k;
                tmpvar.lDS = l;
            }
            Set.append_instance(tmpvar);

        }
    }

    void MCSAnalysis::FillVarsSciFi(Vars& tmpvar, int& j, int& k, int& l, int& pid, double& pz, int& isDS, Vars& proj, Vars& projdif, std::vector<float> chi2nUS, int& beamtype, double& ptruth, const bool& difcut, const bool& fidcut, const bool& chicut, const bool& TOFcut){
        double alX=0.0, alY=0.0, thX=0.0, thY=0.0;
        if(isDS==0){
            _sys["trkrlength"] = 1099.9;
            alX = _sys["alXUS"];
            alY = _sys["alYUS"];
            thX = _sys["thYUS"];
            thY = _sys["thXUS"];
            tmpvar.jUS = j;
            tmpvar.kUS = k;
            tmpvar.lUS = l;
        } else {
            _sys["trkrlength"] = 1099.76;
            alX = _sys["alXDS"];
            alY = _sys["alYDS"];
            thX = _sys["thYDS"];
            thY = _sys["thXDS"];
            tmpvar.jDS = j;
            tmpvar.kDS = k;
            tmpvar.lDS = l;
        }

        if(j < int(scifievent->scifitracks().size()) && j != -1){
            if(k < int(scifievent->scifitracks()[j]->scifitrackpoints().size()) && k != -1){
                tmpvar.X  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().x() + alX;
                tmpvar.Y  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().y() + alY;
                tmpvar.Z  = scifievent->scifitracks()[j]->scifitrackpoints()[k]->pos().z();
                tmpvar.px = scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().x() + tan(thX*atan(1.)/45.0)*tmpvar.pz;
                tmpvar.py = scifievent->scifitracks()[j]->scifitrackpoints()[k]->mom().y() + tan(thY*atan(1.)/45.0)*tmpvar.pz;
                //tmpvar.pz = mcevent->GetVirtualHits()->at(centre).GetMomentum().z();
                tmpvar.pz = pz; //* sqrt(1 + pow(tmpvar.py, 2) + pow(tmpvar.px, 2));
                if (rot_ang==0){
                    tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x() + tan(thX*atan(1.)/45.0);
                    tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y() + tan(thY*atan(1.)/45.0);
                }
                else {

                    if (isMC==1){
                        tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x() + tan(thX*atan(1.)/45.0);
                        tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y() + tan(thY*atan(1.)/45.0);
                    }
                    else {
                        if(isEmpty){
                            if(isDS==0){
                                tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x()*cos(rot_ang_emptyX)-sin(rot_ang_emptyX) + tan(thX*atan(1.)/45.0);
                                tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y()*cos(rot_ang_empty)+sin(rot_ang_empty) + tan(thY*atan(1.)/45.0);
                            } else { 
                                tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x() + tan(thX*atan(1.)/45.0);
                                tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y() + tan(thY*atan(1.)/45.0);
                            }
                        }
                        else {
                            if(isDS==0){
                                tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x()*cos(rot_angX)-sin(rot_angX)+tan(thX*atan(1.)/45.0);
                                tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y()*cos(rot_ang)+sin(rot_ang) + tan(thY*atan(1.)/45.0);
                            } else { 
                                tmpvar.dXdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().x() + tan(thX*atan(1.)/45.0);
                                tmpvar.dYdz = scifievent->scifitracks()[j]->scifitrackpoints()[k]->gradient().y() + tan(thY*atan(1.)/45.0);
                            }
                        }
                    }
                }
                tmpvar.pid = -13;
                if (pid!=-13) tmpvar.pid = pid;
                tmpvar.isgood = true;
            }
        }
        if( j == -1 || k == -1 ){
            // Vars tmpvar;
            tmpvar.X = 0.0;
            tmpvar.Y = 0.0;
            tmpvar.Z = 0.0;
            double v1 = rand() % 100;
            if (v1/100 >0.5) {
                tmpvar.dXdz = -1./2.;
                tmpvar.dYdz = -1./2.;
                tmpvar.px   = -pz/2.;
                tmpvar.py   = -pz/2.;
                tmpvar.pz   = pz;
            }
            else {
                tmpvar.dXdz = 1./2.;
                tmpvar.dYdz = 1./2.;
                tmpvar.px   = pz/2.;
                tmpvar.py   = pz/2.;
                tmpvar.pz   = pz;
            }
            tmpvar.pid = -13;
            tmpvar.isgood = false;
        }
        if( tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size() == 1 &&
                tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray().size() == 1)
            tmpvar.TOF12 = tofevent->GetTOFEventSpacePoint().GetTOF2SpacePointArray()[0].GetTime() 
                - tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime();
        else
            tmpvar.TOF12 = 100.0 * 8.22475 / 0.299792458;
        if( int(tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray().size()) == 1 &&
                int(tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray().size()) == 1 ) {
            tmpvar.TOF01 = tofevent->GetTOFEventSpacePoint().GetTOF1SpacePointArray()[0].GetTime() 
                - tofevent->GetTOFEventSpacePoint().GetTOF0SpacePointArray()[0].GetTime();
        }
        else
            tmpvar.TOF01 = 100.0 * 7.64186 / 0.299792458;
        tmpvar.ptruth = ptruth; 
        tmpvar.difcut = difcut; 
        tmpvar.fidcut = fidcut;
        tmpvar.chicut = chicut;
        tmpvar.TOFcut = TOFcut;
        tmpvar.projX = proj.X;
        tmpvar.projY = proj.Y;
        tmpvar.projdifX = projdif.X;
        tmpvar.projdifY = projdif.Y;
        if (chi2nUS.size()>0) {
            tmpvar.chi2nUS = chi2nUS[0];
            tmpvar.chi2nDS = chi2nUS[1];
        }
        tmpvar.beamtype = beamtype;
        if (USrefplaneI > 0 && USabsPlaneI >0 &&
                DSrefplaneI > 0 && DSabsPlaneI >0 ){
            if (isDS){
                if (mcevent->GetVirtualHits()->size()>DSabsPlaneI) {
                    tmpvar.mcX  = mcevent->GetVirtualHits()->at(DSabsPlaneI).GetPosition().x();
                    tmpvar.mcY  = mcevent->GetVirtualHits()->at(DSabsPlaneI).GetPosition().y();
                    tmpvar.mcZ  = mcevent->GetVirtualHits()->at(DSabsPlaneI).GetPosition().z();
                    tmpvar.mcdXdz = mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().x()/
                        mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().z();
                    tmpvar.mcdYdz = mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().y()/
                        mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().z();
                    tmpvar.mcpx = mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().x();
                    tmpvar.mcpy = mcevent->GetVirtualHits()->at(DSabsPlaneI).GetMomentum().y();

                }
            }
            else{
                if (mcevent->GetVirtualHits()->size()>USabsPlaneI) {
                    tmpvar.mcX  = mcevent->GetVirtualHits()->at(USabsPlaneI).GetPosition().x();
                    tmpvar.mcY  = mcevent->GetVirtualHits()->at(USabsPlaneI).GetPosition().y();
                    tmpvar.mcZ  = mcevent->GetVirtualHits()->at(USabsPlaneI).GetPosition().z();
                    tmpvar.mcdXdz = mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x()/
                        mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().z();
                    tmpvar.mcdYdz = mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().y()/
                        mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().z();
                    tmpvar.mcpx = mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().x();
                    tmpvar.mcpy = mcevent->GetVirtualHits()->at(USabsPlaneI).GetMomentum().y();

                }
            }
        }
    }

    void MCSAnalysis::data_make_beam_histograms(Collection& Set, std::string desc, std::string suffix, Collection& SetMC){

        /* std::cout << desc << std::endl; */
        /* std::cout << Set.N() << std::endl; */
        std::string tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        std::string tmpname  = suffix + "_XY";
        TH2D* XY = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -225, 225);
        tmptitle = desc + ";#it{x} (mm); dXdz";
        tmpname  = suffix + "_XdXdz";
        TH2D* XdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{x} (mm); dYdz";
        tmpname  = suffix + "_XdYdz";
        TH2D* XdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{y} (mm); dXdz";
        tmpname  = suffix + "_YdXdz";
        TH2D* YdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{y} (mm); dYdz";
        tmpname  = suffix + "_YdYdz";
        TH2D* YdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";dXdz; dYdz";
        tmpname  = suffix + "_dXdzdYdz";
        TH2D* dXdzdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{x} (mm); TOF01 (ns)";
        tmpname  = suffix + "_XTOF01";
        TH2D* XTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
        tmptitle = desc + ";#it{y} (mm); TOF01 (ns)";
        tmpname  = suffix + "_YTOF01";
        TH2D* YTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
        tmptitle = desc + ";dXdz; TOF01 (ns)";
        tmpname  = suffix + "_dXdzTOF01";
        TH2D* dXdzTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
        tmptitle = desc + ";dYdz; TOF01 (ns)";
        tmpname  = suffix + "_dYdzTOF01";
        TH2D* dYdzTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
        tmptitle = desc + ";TOF01 (ns)";
        tmpname  = suffix + "_TOF01";
        TH1D* TOF01 = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, 25.0, 35.0);
        tmptitle = desc + ";TOF12 (ns)";
        tmpname  = suffix + "_TOF12";
        TH1D* TOF12 = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, 25.0, 35.0);
        tmptitle = desc + ";dXdz residual";
        tmpname  = suffix + "_dXdzres";
        TH1D* dXdzres = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, -0.125, 0.125);
        tmptitle = desc + ";dYdz residual";
        tmpname  = suffix + "_dYdzres";
        TH1D* dYdzres = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, -0.125, 0.125);
        tmptitle = desc + ";pz residual";
        tmpname  = suffix + "_pzres";
        TH1D* pzres = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, -10.125, 10.125);

        for(int i=0; i<Set.N(); i++){
            /* if(Set.E(i).X==0 && Set.E(i).Y==0) continue; */
            XY->Fill(Set.E(i).X,Set.E(i).Y);
            XdXdz->Fill(Set.E(i).X,Set.E(i).dXdz);
            XdYdz->Fill(Set.E(i).X,Set.E(i).dYdz);
            YdXdz->Fill(Set.E(i).Y,Set.E(i).dXdz);
            YdYdz->Fill(Set.E(i).Y,Set.E(i).dYdz);
            dXdzdYdz->Fill(Set.E(i).dXdz,Set.E(i).dYdz);
            XTOF01->Fill(Set.E(i).X,Set.E(i).TOF01);
            YTOF01->Fill(Set.E(i).Y,Set.E(i).TOF01);
            dXdzTOF01->Fill(Set.E(i).dXdz,Set.E(i).TOF01);
            dYdzTOF01->Fill(Set.E(i).dYdz,Set.E(i).TOF01);
            /* std::cout << Set.E(i).TOF01 << std::endl; */
            TOF01->Fill(Set.E(i).TOF01);
            TOF12->Fill(Set.E(i).TOF12);
            if (isMC==1){
                dXdzres->Fill(Set.E(i).dXdz-SetMC.E(i).px/SetMC.E(i).pz);
                dYdzres->Fill(Set.E(i).dYdz-SetMC.E(i).py/SetMC.E(i).pz);
                pzres->Fill(Set.E(i).pz-SetMC.E(i).pz);
            }
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
        XTOF01->GetXaxis()->SetLabelSize(0.05);
        XTOF01->GetXaxis()->SetTitleSize(0.05);
        XTOF01->GetYaxis()->SetLabelSize(0.05);
        XTOF01->GetYaxis()->SetTitleSize(0.05);
        YTOF01->GetXaxis()->SetLabelSize(0.05);
        YTOF01->GetXaxis()->SetTitleSize(0.05);
        YTOF01->GetYaxis()->SetLabelSize(0.05);
        YTOF01->GetYaxis()->SetTitleSize(0.05);
        dXdzTOF01->GetXaxis()->SetLabelSize(0.05);
        dXdzTOF01->GetXaxis()->SetTitleSize(0.05);
        dXdzTOF01->GetYaxis()->SetLabelSize(0.05);
        dXdzTOF01->GetYaxis()->SetTitleSize(0.05);
        dYdzTOF01->GetXaxis()->SetLabelSize(0.05);
        dYdzTOF01->GetXaxis()->SetTitleSize(0.05);
        dYdzTOF01->GetYaxis()->SetLabelSize(0.05);
        dYdzTOF01->GetYaxis()->SetTitleSize(0.05);
        TOF01->GetXaxis()->SetLabelSize(0.05);
        TOF01->GetXaxis()->SetTitleSize(0.05);
        TOF01->GetYaxis()->SetLabelSize(0.05);
        TOF01->GetYaxis()->SetTitleSize(0.05);
        TOF12->GetXaxis()->SetLabelSize(0.05);
        TOF12->GetXaxis()->SetTitleSize(0.05);
        TOF12->GetYaxis()->SetLabelSize(0.05);
        TOF12->GetYaxis()->SetTitleSize(0.05);


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
        XTOF01->Draw("colz");
        tmpfile = suffix + "_XTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        YTOF01->Draw("colz");
        tmpfile = suffix + "_YTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        dXdzTOF01->Draw("colz");
        tmpfile = suffix + "_dXdzTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        dYdzTOF01->Draw("colz");
        tmpfile = suffix + "_dYdzTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        TOF01->Draw("colz");
        tmpfile = suffix + "_TOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        TOF12->Draw("colz");
        tmpfile = suffix + "_TOF12.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        dXdzres->Draw("colz");
        tmpfile = suffix + "_dXdzres.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        dYdzres->Draw("colz");
        tmpfile = suffix + "_dYdzres.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        pzres->Draw("colz");
        tmpfile = suffix + "_pzres.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();

        outfile->cd();
        XY->Write();
        XdXdz->Write();
        XdYdz->Write();
        XTOF01->Write();
        YdXdz->Write();
        YdYdz->Write();
        YTOF01->Write();
        dXdzdYdz->Write();
        dXdzTOF01->Write();
        dYdzTOF01->Write();
        TOF01->Write();
        TOF12->Write();
        dXdzres->Write();
        dYdzres->Write();
        pzres->Write();

        delete c1;
        delete XY;
        delete XdXdz;
        delete XdYdz;
        delete YdXdz;
        delete YdYdz;
        delete XTOF01;
        delete YTOF01;
        delete dXdzdYdz;
        delete dXdzTOF01;
        delete dYdzTOF01;
        delete TOF01;
        delete TOF12;
        delete dXdzres;
        delete dYdzres;
        delete pzres;
    }

    void MCSAnalysis::make_beam_histograms(Collection& Set, std::string desc, std::string suffix){

        std::string tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        std::string tmpname  = suffix + "_XY";
        TH2D* XY = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -225, 225);
        tmptitle = desc + ";#it{x} (mm); dXdz";
        tmpname  = suffix + "_XdXdz";
        TH2D* XdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{x} (mm); dYdz";
        tmpname  = suffix + "_XdYdz";
        TH2D* XdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{y} (mm); dXdz";
        tmpname  = suffix + "_YdXdz";
        TH2D* YdXdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{y} (mm); dYdz";
        tmpname  = suffix + "_YdYdz";
        TH2D* YdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, -0.125, 0.125);
        tmptitle = desc + ";dXdz; dYdz";
        tmpname  = suffix + "_dXdzdYdz";
        TH2D* dXdzdYdz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
        tmptitle = desc + ";#it{x} (mm); TOF01 (ns)";
        tmpname  = suffix + "_XTOF01";
        TH2D* XTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
        tmptitle = desc + ";#it{y} (mm); TOF01 (ns)";
        tmpname  = suffix + "_YTOF01";
        TH2D* YTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -225, 225, 90, 20.0, 45.0);
        tmptitle = desc + ";dXdz; TOF01 (ns)";
        tmpname  = suffix + "_dXdzTOF01";
        TH2D* dXdzTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
        tmptitle = desc + ";dYdz; TOF01 (ns)";
        tmpname  = suffix + "_dYdzTOF01";
        TH2D* dYdzTOF01 = new TH2D(tmpname.c_str(),tmptitle.c_str(), 90, -0.125, 0.125, 90, 20.0, 45.0);
        tmptitle = desc + ";TOF01 (ns)";
        tmpname  = suffix + "_TOF01";
        TH1D* TOF01 = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, 25.0, 35.0);
        tmptitle = desc + ";TOF12 (ns)";
        tmpname  = suffix + "_TOF12";
        TH1D* TOF12 = new TH1D(tmpname.c_str(),tmptitle.c_str(), 100, 25.0, 35.0);

        for(int i=0; i<Set.N(); i++){
            /* if(Set.E(i).X==0 && Set.E(i).Y==0) continue; */
            XY->Fill(Set.E(i).X,Set.E(i).Y);
            XdXdz->Fill(Set.E(i).X,Set.E(i).dXdz);
            XdYdz->Fill(Set.E(i).X,Set.E(i).dYdz);
            YdXdz->Fill(Set.E(i).Y,Set.E(i).dXdz);
            YdYdz->Fill(Set.E(i).Y,Set.E(i).dYdz);
            dXdzdYdz->Fill(Set.E(i).dXdz,Set.E(i).dYdz);
            XTOF01->Fill(Set.E(i).X,Set.E(i).TOF01);
            YTOF01->Fill(Set.E(i).Y,Set.E(i).TOF01);
            dXdzTOF01->Fill(Set.E(i).dXdz,Set.E(i).TOF01);
            dYdzTOF01->Fill(Set.E(i).dYdz,Set.E(i).TOF01);
            TOF01->Fill(Set.E(i).TOF01);
            TOF12->Fill(Set.E(i).TOF12);
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
        XTOF01->GetXaxis()->SetLabelSize(0.05);
        XTOF01->GetXaxis()->SetTitleSize(0.05);
        XTOF01->GetYaxis()->SetLabelSize(0.05);
        XTOF01->GetYaxis()->SetTitleSize(0.05);
        YTOF01->GetXaxis()->SetLabelSize(0.05);
        YTOF01->GetXaxis()->SetTitleSize(0.05);
        YTOF01->GetYaxis()->SetLabelSize(0.05);
        YTOF01->GetYaxis()->SetTitleSize(0.05);
        dXdzTOF01->GetXaxis()->SetLabelSize(0.05);
        dXdzTOF01->GetXaxis()->SetTitleSize(0.05);
        dXdzTOF01->GetYaxis()->SetLabelSize(0.05);
        dXdzTOF01->GetYaxis()->SetTitleSize(0.05);
        dYdzTOF01->GetXaxis()->SetLabelSize(0.05);
        dYdzTOF01->GetXaxis()->SetTitleSize(0.05);
        dYdzTOF01->GetYaxis()->SetLabelSize(0.05);
        dYdzTOF01->GetYaxis()->SetTitleSize(0.05);
        TOF01->GetXaxis()->SetLabelSize(0.05);
        TOF01->GetXaxis()->SetTitleSize(0.05);
        TOF01->GetYaxis()->SetLabelSize(0.05);
        TOF01->GetYaxis()->SetTitleSize(0.05);
        TOF12->GetXaxis()->SetLabelSize(0.05);
        TOF12->GetXaxis()->SetTitleSize(0.05);
        TOF12->GetYaxis()->SetLabelSize(0.05);
        TOF12->GetYaxis()->SetTitleSize(0.05);


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
        XTOF01->Draw("colz");
        tmpfile = suffix + "_XTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        YTOF01->Draw("colz");
        tmpfile = suffix + "_YTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        dXdzTOF01->Draw("colz");
        tmpfile = suffix + "_dXdzTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        dYdzTOF01->Draw("colz");
        tmpfile = suffix + "_dYdzTOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        TOF01->Draw("colz");
        tmpfile = suffix + "_TOF01.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();
        TOF12->Draw("colz");
        tmpfile = suffix + "_TOF12.pdf";
        c1->Print(tmpfile.c_str());
        c1->Clear();

        outfile->cd();
        XY->Write();
        XdXdz->Write();
        XdYdz->Write();
        XTOF01->Write();
        YdXdz->Write();
        YdYdz->Write();
        YTOF01->Write();
        dXdzdYdz->Write();
        dXdzTOF01->Write();
        dYdzTOF01->Write();
        TOF01->Write();
        TOF12->Write();

        delete c1;
        delete XY;
        delete XdXdz;
        delete XdYdz;
        delete YdXdz;
        delete YdYdz;
        delete XTOF01;
        delete YTOF01;
        delete dXdzdYdz;
        delete dXdzTOF01;
        delete dYdzTOF01;
        delete TOF01;
        delete TOF12;
    }

    void MCSAnalysis::make_acceptance_histograms(Collection& USset, Collection& DSset, 
            std::string desc, std::string suffix){
        std::string tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        std::string tmpname  = suffix + "_posaccXY";
        TH2D* posaccXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
        tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        tmpname  = suffix + "_divaccXY";
        TH2D* divaccXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
        tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        tmpname  = suffix + "_posresXY";
        TH2D* posresXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
        tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        tmpname  = suffix + "_divresXY";
        TH2D* divresXY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
        tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        tmpname  = suffix + "_posres2XY";
        TH2D* posres2XY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
        tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
        tmpname  = suffix + "_divres2XY";
        TH2D* divres2XY = new TH2D(tmpname.c_str(), tmptitle.c_str(), 75, -300, 300, 75, -300, 300);
        tmptitle = desc + ";#it{x} (mm); #it{y} (mm)";
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

        tmptitle = desc + ";|#hat{r}^{proj to LiH}_{US} - #hat{r}^{proj to LiH}_{DS} | (mm); No of events";
        tmpname  = suffix + "_appRef";
        TH1D* appRef = new TH1D(tmpname.c_str(), tmptitle.c_str(), 40, 0, 200);
        tmptitle = desc + ";X_{US}^{proj to LiH} - X_{DS}^{proj to LiH} (mm);Y_{US}^{proj to LiH} - Y_{DS}^{proj to LiH} (mm)";
        tmpname  = suffix + "_projRefDiff";
        TH2D* projRefDiff = new TH2D(tmpname.c_str(), tmptitle.c_str(), 200, -200, 200, 200, -200, 200);
        tmptitle = desc + ";dXdz_{US}^{proj to LiH} - dXdz_{DS}^{proj to LiH};dYdz_{US}^{proj to LiH} - dYdz_{DS}^{proj to LiH}";
        tmpname  = suffix + "_projRefDiffddz";
        TH2D* projRefDiffddz = new TH2D(tmpname.c_str(), tmptitle.c_str(), 90, -0.125, 0.125, 90, -0.125, 0.125);
        tmptitle = desc + ";Position of closest approach (mm); No of events";
        tmpname  = suffix + "_cloesestapp";
        TH1D* closestapp = new TH1D(tmpname.c_str(), tmptitle.c_str(), 40, _sys["USref"], _sys["DSref"]);

        double latchZ = _sys["abspos"];
        //19398.9;
        Vars USproj;
        Vars DSproj;
        Vars USprojC;
        Vars DSprojC;
            double approachclosestappvalue = 99999;
        for (int i=0; i<USset.N(); i++){
            /*
               if (DSset.E(i).Z != 0 && latchZ == 0.0){
               latchZ = DSset.E(i).Z;
               } 
               if (latchZ != 0.0 && fabs(DSset.E(i).Z - latchZ) < 5){
               DSproj = PropagateVarsMu(USset.E(i), DSset.E(i).Z);
               */
            //} else if(latchZ != 0.0){
            USproj = PropagateVarsMu(USset.E(i), latchZ);
            DSproj = PropagateVarsMu(DSset.E(i), latchZ);
            //} else {
            //  continue;
            //}

            Vars predDiff = USproj - DSproj;
            /* std::cout << "USproj.Z " << USproj.Z << std::endl; */
            /* std::cout << "DSproj.Z " << DSproj.Z << std::endl; */
            /* std::cout << "predDiff.dXdz " << predDiff.dXdz << std::endl; */
            /* std::cout << "predDiff.dYdz " << predDiff.dYdz << std::endl; */
            /* std::cout << "predDiff.Z " << predDiff.Z << std::endl; */
            double approach = sqrt(predDiff.X * predDiff.X + 
                    predDiff.Y * predDiff.Y +
                    predDiff.Z * predDiff.Z);
            appRef->Fill(approach); // can be very large if the ds track is lost

            projRefDiff->Fill(predDiff.X, predDiff.Y);
            posXY->Fill(DSproj.X, DSproj.Y);
            posDXDY->Fill(DSproj.dXdz, DSproj.dYdz);
            projRefDiffddz->Fill(predDiff.dXdz, predDiff.dYdz);
            /* projRefDiffddz->Fill(USproj.dXdz-DSproj.dXdz, USproj.dYdz-DSproj.dXdz); */

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
            double approachC = -9999;
            double closestpredz = -9999;

            if (DSset.E(i).X==-9999) closestpredz = -9999;
            else{
            for (int jjj=int(_sys["USref"]); jjj<int(_sys["DSref"]); jjj++){

            USprojC = PropagateVarsMu(USset.E(i), jjj);
            DSprojC = PropagateVarsMu(DSset.E(i), jjj);
            Vars predDiffC = USprojC - DSprojC;
            approachC = sqrt(predDiffC.X * predDiffC.X + 
                    predDiffC.Y * predDiffC.Y +
                    predDiffC.Z * predDiffC.Z);

            /* std::cout << " " << std::endl; */
            /* std::cout << "USprojC.X " << USprojC.X << std::endl; */
            /* std::cout << "USprojC.Y " << USprojC.Y << std::endl; */
            /* std::cout << "USprojC.Z " << USprojC.Z << std::endl; */
            /* std::cout << "DSprojC.X " << DSprojC.X << std::endl; */
            /* std::cout << "DSprojC.Y " << DSprojC.Y << std::endl; */
            /* std::cout << "DSprojC.Z " << DSprojC.Z << std::endl; */
            /* std::cout << "predDiffC.X " << predDiffC.X << std::endl; */
            /* std::cout << "predDiffC.Y " << predDiffC.Y << std::endl; */
            /* std::cout << "predDiffC.Z " << predDiffC.Z << std::endl; */
            /* std::cout << " " << std::endl; */
            
            /* std::cout << "approachC " << approachC << std::endl; */
            /* std::cout << "approachclosestappvalue " << approachclosestappvalue << std::endl; */
            if (approachC<approachclosestappvalue) {
            /* std::cout << "Updating closestpredz to " << USprojC.Z << std::endl; */
                closestpredz = USprojC.Z;
                approachclosestappvalue=approachC;
            }
            }
            }
            /* std::cout << "USset.E(i).X " << USset.E(i).X << std::endl; */
            /* std::cout << "USset.E(i).Y " << USset.E(i).Y << std::endl; */
            /* std::cout << "USset.E(i).Z " << USset.E(i).Z << std::endl; */
            /* std::cout << "DSset.E(i).X " << DSset.E(i).X << std::endl; */
            /* std::cout << "DSset.E(i).Y " << DSset.E(i).Y << std::endl; */
            /* std::cout << "DSset.E(i).Z " << DSset.E(i).Z << std::endl; */
            /* std::cout << "closestpredz " << closestpredz << std::endl; */
            if (closestpredz<18700 && closestpredz>15200) closestapp->Fill(closestpredz);
            approachclosestappvalue = 99999;
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

        TCanvas *c1 = new TCanvas();
        appRef->SetMarkerStyle(20);
        appRef->Draw();
        tmpname = suffix.c_str();
        tmpname += "appRef.pdf";
        c1->SaveAs(tmpname.c_str());
        c1->Clear();
        projRefDiff->Draw("P");
        tmpname = suffix.c_str();
        tmpname += "projRefDiff.pdf";
        c1->SaveAs(tmpname.c_str());
        c1->Clear();
        projRefDiffddz->Draw("P");
        tmpname = suffix.c_str();
        tmpname += "projRefDiffddz.pdf";
        c1->SaveAs(tmpname.c_str());
        c1->Clear();
        closestapp->Draw();
        tmpname = suffix.c_str();
        tmpname += "cloesestapp.pdf";
        c1->SaveAs(tmpname.c_str());
        delete c1;

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
        projRefDiffddz->Write();
        closestapp->Write();

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

    void MCSAnalysis::make_scattering_acceptance_histograms(Collection& USset,
            Collection& DSset,
            Collection& DSrefset,
            std::string desc,
            std::string suffix){

        std::string tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        std::string tmpname = "thetaY_all_";
        tmpname += suffix;
        TH1D* thetaY_all = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}a}} (radians)";
        tmpname = "thetaX_all_";
        tmpname += suffix;
        TH1D* thetaX_all = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        tmpname = "thetaY_acc_";
        tmpname += suffix;
        TH1D* thetaY_acc = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)";
        tmpname = "thetaX_acc_";
        tmpname += suffix;
        TH1D* thetaX_acc = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        tmpname = "thetaY_divacc_";
        tmpname += suffix;
        TH1D* thetaY_divacc = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)";
        tmpname = "thetaX_divacc_";
        tmpname += suffix;
        TH1D* thetaX_divacc = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        tmpname = "thetaY_res_";
        tmpname += suffix;
        TH1D* thetaY_res = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)";
        tmpname = "thetaX_res_";
        tmpname += suffix;
        TH1D* thetaX_res = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        tmpname = "thetaY_divres_";
        tmpname += suffix;
        TH1D* thetaY_divres = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)";
        tmpname = "thetaX_divres_";
        tmpname += suffix;
        TH1D* thetaX_divres = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        tmpname = "thetaY_res2_";
        tmpname += suffix;
        TH1D* thetaY_res2 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)";
        tmpname = "thetaX_res2_";
        tmpname += suffix;
        TH1D* thetaX_res2 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        tmpname = "thetaY_divres2_";
        tmpname += suffix;
        TH1D* thetaY_divres2 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)";
        tmpname = "thetaX_divres2_";
        tmpname += suffix;
        TH1D* thetaX_divres2 = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        for (int i=0; i<USset.N(); i++){
            if (i >= DSrefset.N()) 
                break;
            std::vector<double> projTheta = RotDefineProjectionAngles(USset.E(i), DSset.E(i),angdef);
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

    void MCSAnalysis::trkr_eff(Collection& USTruthSethold,
            Collection& DSTruthSethold,
            Collection& _USsethold,
            Collection& _DSsethold,
            std::string desc,
            std::string suffix){

        std::string tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        std::string tmpname = "thetaY_effi_";
        tmpname += suffix;
        /* TEfficiency* thetaY_effi = */ 
        /*     new TEfficiency(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad", */
        /*             _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]); */
        /* tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)"; */
        /* tmpname = "thetaX_effi_"; */
        /* tmpname += suffix; */
        /* TEfficiency* thetaX_effi = */ 
        /*     new TEfficiency(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad", */
        /*             _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]); */

        tmptitle = desc + ";#it{#theta_{#it{Y}}} (radians)";
        tmpname = "thetaY_reco_";
        tmpname += suffix;
        TH1D* thetaY_reco = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (Y);#it{#theta_{Y}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);
        tmptitle = desc + ";#it{#theta_{#it{X}}} (radians)";
        tmpname = "thetaX_reco_";
        tmpname += suffix;
        TH1D* thetaX_reco = 
            new TH1D(tmpname.c_str(),"Change in Projected Angle (X);#it{#theta_{X}}; Events per 4 mrad",
                    _histlimits["NbinsXY"], _histlimits["minXY"], _histlimits["maxXY"]);

        for (int i=0; i<_DSsethold.N(); i++){
            std::vector<double> projTheta = RotDefineProjectionAngles(_USsethold.E(i), _DSsethold.E(i),angdef);
            double thetaX = projTheta[0];
            double thetaY = projTheta[1];
            if (sqrt(projTheta[0]*projTheta[0]+projTheta[1]*projTheta[1])<eff_cut){
                thetaX_reco->Fill(thetaX);
                thetaY_reco->Fill(thetaY);
            }
        }
        /* thetaX_reco->Fill(0); */
        /*     thetaY_reco->Fill(0); */
    TText t1 = TText(0.12,0.785,"MICE Preliminary");
    TText t3 = TText(0.12,0.75,"ISIS cycle 2015/04");
    TText t2 = TText(0.12,0.715,"LiH, MAUS v3.3.2");
            t1.SetNDC(1);
            t1.SetTextSize(0.04);
            t1.SetTextFont(42);
            t2.SetNDC(1);
            t2.SetTextSize(0.03);
            t2.SetTextFont(42);
            string str;
            if (to_string(eff_cut)!=0.180) {
                str = "tracker_resolution_plots_"+to_string(TOF_lower_limit)+"_"+to_string(eff_cut)+".root";
            }
            else {
                str = "tracker_resolution_plots_"+to_string(TOF_lower_limit)+".root";
            }
            TFile *f = new TFile(str.c_str(),"RECREATE");
            TEfficiency* thetaX_effi = new TEfficiency(*thetaX_reco,*theta_true_x_graph);
            TCanvas *c5 = new TCanvas();
            TGraphAsymmErrors *thetaX_effi_graph = thetaX_effi->CreateGraph();
            thetaX_effi_graph->SetName("Effx_graph");
            thetaX_effi_graph->GetXaxis()->SetTitle("#it{#theta_{X}} (mrad);");
            thetaX_effi_graph->GetYaxis()->SetTitle("Efficiency");
            thetaX_effi_graph->Draw("ap");
            t1.Draw();
            t2.Draw();
            t3.Draw();
            /* t1.Paint(); */
            /* t2.Paint(); */
            c5->SaveAs("pEff_x.pdf");
            f->cd();
            thetaX_effi_graph->Write();

            TEfficiency *thetaY_effi = new TEfficiency(*thetaY_reco,*theta_true_y_graph);
            TCanvas *c7 = new TCanvas();
            TGraphAsymmErrors *thetaY_effi_graph = thetaY_effi->CreateGraph();
            thetaY_effi_graph->SetName("Effy_graph");
            thetaY_effi_graph->GetXaxis()->SetTitle("#it{#theta_{Y}} (mrad);");
            thetaY_effi_graph->GetYaxis()->SetTitle("Efficiency");
            thetaY_effi_graph->Draw("ap");
            t1.Draw();
            t2.Draw();
            t3.Draw();
            /* t1.Paint(); */
            /* t2.Paint(); */
            c7->SaveAs("pEff_y.pdf");
            thetaY_effi_graph->Write();
            f->Close();
    outfile->cd();
    thetaX_reco->Write();
    thetaY_reco->Write();

        /* delete thetaX_effi_graph; */
        /* delete thetaY_effi_graph; */
        /* delete thetaX_reco; */
        /* delete thetaY_reco; */
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

    Vars MCSAnalysis::PropagateVarsMu(const Vars& event, const double & z0){


        double mass = 105.658; // muon mass [MeV/c^2]
        double time = 0.; // time in lab frame [ns]
        double x = event.X; // horizontal position [mm]
        double y = event.Y; // vertical position [mm]
        double z = event.Z; // longitudinal position [mm]
        double px = event.px; // horizontal momentum component [MeV/c]
        double py = event.py; // vertical momentum component [MeV/c]
        double pz = event.pz;
        double dXdz = event.dXdz; 
        double dYdz = event.dYdz; 

        double energy = std::sqrt(pz*pz + px*px + py*py + mass*mass); // Total energy [MeV]

        double event_vector[10] = { time, x, y, z, energy, px, py, pz, dXdz, dYdz};
        /*
           BTField* field = dynamic_cast<BTField*>(MAUS::Globals::GetMCFieldConstructor());
           try {
           MAUS::GlobalTools::propagate(event_vector, z0, field, 10., 
           MAUS::DataStructure::Global::kMuPlus, true);
           } catch (...){
           */
        // Assume a straight track
        /*
           std::cout << "x " << x << std::endl;
           std::cout << "z0 " << z0 << std::endl;
           std::cout << "z " << z << std::endl;
           std::cout << "px " << px << std::endl;
           std::cout << "pz " << pz << std::endl;
           */
        //event_vector[1] = px/pz * (z0 - z) + x;
        //    std::cout << "event_vector[1] " << event_vector[1] << std::endl;
        //event_vector[2] = py/pz * (z0 - z) + y;
        event_vector[1] = dXdz * (z0 - z) + x;
        event_vector[2] = dYdz * (z0 - z) + y;
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
        //prop.dXdz = event_vector[5]/event_vector[7];
        //prop.dYdz = event_vector[6]/event_vector[7];
        prop.dXdz = event_vector[8];
        prop.dYdz = event_vector[9];
        prop.TOF12 = event.TOF12;
        prop.TOF01 = event.TOF01;

        return prop;
    }

    TH1D *defineHist2(const char* fname, const char* ftitle, Double_t fnbinsx, Double_t fxlow, Double_t fxup)
    {
        TH1D *fhis1D = new TH1D(fname, ftitle, fnbinsx, fxlow, fxup);
        fhis1D->SetMinimum(0.001);
        fhis1D->GetXaxis()->SetTitle(ftitle);
        Double_t binning = ( (fhis1D->GetXaxis()->GetXmax()) - (fhis1D->GetXaxis()->GetXmin()) ) / (fhis1D->GetNbinsX());
        fhis1D->GetYaxis()->SetTitle( Form("Events / (%.2f [Mev/c^{2}])", binning) );
        fhis1D->GetYaxis()->SetTitleOffset(1.6);
        fhis1D->GetYaxis()->SetLabelSize(0.035);
        return fhis1D;
    }

