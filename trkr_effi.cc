/* This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
 *
 * MAUS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MAUS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MAUS.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdlib.h>

#include <string>
#include <cstddef> 
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"

// Reconstruction specific headers
#include "src/common_cpp/DataStructure/TOFEvent.hh"
#include "src/common_cpp/DataStructure/SciFiEvent.hh"
#include "src/common_cpp/DataStructure/KLEvent.hh"
#include "src/common_cpp/DataStructure/CkovEvent.hh"
#include "src/common_cpp/DataStructure/EMREvent.hh"
// MC specific headers
#include "src/common_cpp/DataStructure/Primary.hh"
#include "src/common_cpp/DataStructure/VirtualHit.hh"
#include "src/common_cpp/DataStructure/Track.hh"
#include "src/common_cpp/DataStructure/Hit.hh"
#include "src/common_cpp/DataStructure/MCEvent.hh"
// General headers
#include "src/common_cpp/DataStructure/Spill.hh"
#include "src/common_cpp/DataStructure/Data.hh"
#include "src/common_cpp/JsonCppStreamer/IRStream.hh"

#include "Collection.h"

double TOF_upper_limit = 28;
double TOF_lower_limit = 30;
TVector3 x;
TVector3 v_;
TVector3 u;
TVector3 w;
double USXproj;
double USYproj;
double USYdiff;
double USXdiff;
    TH1D MC_theta_x("MC_theta_x", "scat_ang;No. of entries", 47, -0.0705, 0.0705);
    TH1D MC_theta_y("MC_theta_x", "scat_ang;No. of entries", 47, -0.0705, 0.0705);
    TH1D recon_theta_x("recon_theta_x", "scat_ang;No. of entries", 47, -0.0705, 0.0705);
    TH1D recon_theta_y("recon_theta_x", "scat_ang;No. of entries", 47, -0.0705, 0.0705);

void make_beam_histograms(Collection Set, std::string desc, std::string suffix, TFile* outfile){
  
  std::string tmptitle = desc + ";X (mm); Y (mm)";
  std::string tmpname  = suffix + "_XY";
  TH2D* XY = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -1000, 1000, 50, -1000, 1000);
  tmptitle = desc + ";X (mm); Px";
  tmpname  = suffix + "_XPx";
  TH2D* XPx = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -1000, 1000, 50, -100, 100);
  tmptitle = desc + ";X (mm); Py";
  tmpname  = suffix + "_XPy";
  TH2D* XPy = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -1000, 1000, 50, -100, 100);
  tmptitle = desc + ";Y (mm); Px";
  tmpname  = suffix + "_YPx";
  TH2D* YPx = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -1000, 1000, 50, -100, 100);
  tmptitle = desc + ";Y (mm); Py";
  tmpname  = suffix + "_YPy";
  TH2D* YPy = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -1000, 1000, 50, -100, 100);
  tmptitle = desc + ";Px; Py";
  tmpname  = suffix + "_PxPy";
  TH2D* PxPy = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -100, 100, 50, -100, 100);
  tmptitle = desc + ";X (mm); Pz (MeV/c)";
  tmpname  = suffix + "_XPz";
  TH2D* XPz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -1000, 1000, 50, 20.0, 450.0);
  tmptitle = desc + ";Y (mm); Pz (MeV/c)";
  tmpname  = suffix + "_YPz";
  TH2D* YPz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -1000, 1000, 50, 20.0, 450.0);
  tmptitle = desc + ";Px; Pz (MeV/c)";
  tmpname  = suffix + "_PxPz";
  TH2D* PxPz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -100, 100, 50, 20.0, 450.0);
  tmptitle = desc + ";Py; Pz (MeV/c)";
  tmpname  = suffix + "_PyPz";
  TH2D* PyPz = new TH2D(tmpname.c_str(),tmptitle.c_str(), 50, -100, 100, 50, 20.0, 450.0);
  
  for(int i=0; i<Set.N(); i++){
    XY->Fill(Set.E(i).X,Set.E(i).Y);
    XPx->Fill(Set.E(i).X,Set.E(i).dXdz);
    XPy->Fill(Set.E(i).X,Set.E(i).dYdz);
    YPx->Fill(Set.E(i).Y,Set.E(i).dXdz);
    YPy->Fill(Set.E(i).Y,Set.E(i).dYdz);
    PxPy->Fill(Set.E(i).dXdz,Set.E(i).dYdz);
    XPz->Fill(Set.E(i).X,Set.E(i).TOF12);
    YPz->Fill(Set.E(i).Y,Set.E(i).TOF12);
    PxPz->Fill(Set.E(i).dXdz,Set.E(i).TOF12);
    PyPz->Fill(Set.E(i).dYdz,Set.E(i).TOF12);
  }
  
  XY->GetXaxis()->SetLabelSize(0.05);
  XY->GetXaxis()->SetTitleSize(0.05);
  XY->GetYaxis()->SetLabelSize(0.05);
  XY->GetYaxis()->SetTitleSize(0.05);
  XPx->GetXaxis()->SetLabelSize(0.05);
  XPx->GetXaxis()->SetTitleSize(0.05);
  XPx->GetYaxis()->SetLabelSize(0.05);
  XPx->GetYaxis()->SetTitleSize(0.05);
  XPy->GetXaxis()->SetLabelSize(0.05);
  XPy->GetXaxis()->SetTitleSize(0.05);
  XPy->GetYaxis()->SetLabelSize(0.05);
  XPy->GetYaxis()->SetTitleSize(0.05);
  YPx->GetXaxis()->SetLabelSize(0.05);
  YPx->GetXaxis()->SetTitleSize(0.05);
  YPx->GetYaxis()->SetLabelSize(0.05);
  YPx->GetYaxis()->SetTitleSize(0.05);
  YPy->GetXaxis()->SetLabelSize(0.05);
  YPy->GetXaxis()->SetTitleSize(0.05);
  YPy->GetYaxis()->SetLabelSize(0.05);
  YPy->GetYaxis()->SetTitleSize(0.05);
  PxPy->GetXaxis()->SetLabelSize(0.05);
  PxPy->GetXaxis()->SetTitleSize(0.05);
  PxPy->GetYaxis()->SetLabelSize(0.05);
  PxPy->GetYaxis()->SetTitleSize(0.05);
  XPz->GetXaxis()->SetLabelSize(0.05);
  XPz->GetXaxis()->SetTitleSize(0.05);
  XPz->GetYaxis()->SetLabelSize(0.05);
  XPz->GetYaxis()->SetTitleSize(0.05);
  YPz->GetXaxis()->SetLabelSize(0.05);
  YPz->GetXaxis()->SetTitleSize(0.05);
  YPz->GetYaxis()->SetLabelSize(0.05);
  YPz->GetYaxis()->SetTitleSize(0.05);
  PxPz->GetXaxis()->SetLabelSize(0.05);
  PxPz->GetXaxis()->SetTitleSize(0.05);
  PxPz->GetYaxis()->SetLabelSize(0.05);
  PxPz->GetYaxis()->SetTitleSize(0.05);
  PyPz->GetXaxis()->SetLabelSize(0.05);
  PyPz->GetXaxis()->SetTitleSize(0.05);
  PyPz->GetYaxis()->SetLabelSize(0.05);
  PyPz->GetYaxis()->SetTitleSize(0.05);

  outfile->cd();
  XY->Write();
  XPx->Write();
  XPy->Write();
  XPz->Write();
  YPx->Write();
  YPy->Write();
  YPz->Write();
  PxPy->Write();
  PxPz->Write();
  PyPz->Write();
  
  TCanvas* c1 = new TCanvas();
  XY->Draw("colz");
  std::string tmpfile = suffix + "_XY.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  XPx->Draw("colz");
  tmpfile = suffix + "_XPx.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  XPy->Draw("colz");
  tmpfile = suffix + "_XPy.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  YPx->Draw("colz");
  tmpfile = suffix + "_YPx.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  YPy->Draw("colz");
  tmpfile = suffix + "_YPy.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  PxPy->Draw("colz");
  tmpfile = suffix + "_PxPy.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  XPz->Draw("colz");
  tmpfile = suffix + "_XPz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  YPz->Draw("colz");
  tmpfile = suffix + "_YPz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  PxPz->Draw("colz");
  tmpfile = suffix + "_PxPz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  PyPz->Draw("colz");
  tmpfile = suffix + "_PyPz.pdf";
  c1->Print(tmpfile.c_str());
  c1->Clear();
  
  
  delete c1;
  delete XY;
  delete XPx;
  delete XPy;
  delete YPx;
  delete YPy;
  delete XPz;
  delete YPz;
  delete PxPy;
  delete PxPz;
  delete PyPz;
}

/** Perform a toy analysis using ROOT
 *
 *  Plot digits for TOF1 plane 0 and plane 1. See documentation in
 *  maus_user_guide "Accessing Data"
 */
int main(int argc, char* argv[]) {

    std::vector<std::string> filelist;
    if(argc == 1){
        std::cout<<"This is a MAUS tree skim that selects events that \n";
        std::cout<<"produce triggers in TOF2 and writes the reconstructed event \n";
        std::cout<<"objects and the MCevent objects to a tree. The resulting \n";
        std::cout<<"tree may be used with the \"tree->Draw()\" command on the \n";
        std::cout<<"root command line or with the \"run_sp_ensemble_analysis\"\n";
        std::cout<<"executable.\n"<<std::endl;
        std::cout<<"To run the command, provide one or more input files as arguements.\n";
    }
    std::string type  = argv[1];
    std::string type1 = argv[1];
    for (int ifile = 1; ifile < argc; ifile++){
        TFile* file = new TFile(argv[ifile]);
        if(file->IsZombie()) continue;
        TTree* Spill = (TTree*)file->Get("Spill");
        if(!Spill) continue;
        file->Close();
        std::cout<<"Adding "<<argv[ifile]<<" to file list."<<std::endl;
        filelist.push_back(argv[ifile]);
    }

    MAUS::Data *data = new MAUS::Data();
    double tof1pos = 12929.2608098;
    double tof2pos = 20500;
    // Make histograms; to fill later on
    TH1D tof2_digits_0_hist("tof2 digits_0",
            "tof2 digits for plane 0;Slab number",
            10, -0.5, 9.5);
    TH1D tof2_digits_1_hist("tof2 digits_1",
            "tof2 digits for plane 1;Slab number",
            10, -0.5, 9.5);
    TH1D hrate("hrate", 
            ";Position along Z axis (m);Particles at Virtual Plane", \
            240,0.5,24.5);
    TH2D hr("hr", 
            ";Position along Z axis (m); Radial Position (m)",	
            240,0.5,24.5,100,0.0, 0.5);
    TH2D hx("hx", 
            ";Position along Z axis (m); X Position (m)",	
            240,0.5,24.5,200,-1.0, 1.0);
    TH3D hxy("hxy", 
            ";Position along Z axis (m); X Position (m); Y Position (m)",	
            240,0.5,24.5,100,-1.0, 1.0, 100, -1.0, 1.0);
    TH2D hp("hp", 
            ";Position along Z axis (m); Momentum (MeV/c)",	
            240,0.5,24.5,100,0,400);
    TH2D hdp("hdp", 
            ";Position along Z axis (m); Change in Momentum (MeV/c)",	
            240,0.5,24.5,500,0,200);
    // Create a new file to store the output
    std::string fname = "reduced_tree";
    //int found = type.find_last_of("/");
    //fname += type.substr(found-9,9);
    fname += ".root";
    TFile* outfile = new TFile(fname.c_str(),"RECREATE");

    std::cout<<"Make tree; to fill later on"<<std::endl;
    TTree* tree = new TTree("recon_reduced_tree","Reduced MAUS Spill");
    TTree* mctree = new TTree("Truth_reduced_tree","Reduced MAUS MCSpill");
    int RunNumber, SpillNumber;
    MAUS::TOFEvent* tofevent = new MAUS::TOFEvent();
    MAUS::SciFiEvent* scifievent = new MAUS::SciFiEvent();
    MAUS::CkovEvent* ckovevent = new MAUS::CkovEvent();
    MAUS::KLEvent* klevent = new MAUS::KLEvent();
    MAUS::EMREvent* emrevent = new MAUS::EMREvent();
    MAUS::MCEvent* mcevent = new MAUS::MCEvent();
    //MAUS::Primary* primary = new MAUS::Primary();
    //MAUS::SpecialVirtualHitArray* sphitarray = new MAUS::SpecialVirtualHitArray();
    // MAUS::TOFHitArray* tofhitarray = new MAUS::TOFHitArray();
    tree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    tree->Branch("SpillNumber", &SpillNumber, "SpillNumber/I");
    tree->Branch("TOFBranch", "MAUS::TOFEvent", &tofevent, 64000, 10);
    tree->Branch("SciFiBranch", "MAUS::SciFiEvent", &scifievent, 64000, 10);
    tree->Branch("CkovBranch", "MAUS::CkovEvent",&ckovevent, 64000, 10);
    tree->Branch("KLBranch", "MAUS::KLEvent", &klevent, 64000, 10);
    tree->Branch("EMRBranch", "MAUS::EMREvent", &emrevent, 64000, 10);
    tree->Branch("MCEvent", "MAUS::MCEvent", &mcevent, 64000, 10); 
    // mctree->Branch("Primary", "MAUS::Primary", &primary, 6400, 10);
    // mctree->Branch("SpecialVirtuals", "MAUS::SpecialVirtualHitArray", &sphitarray, 64000, 10); 
    mctree->Branch("MCEvent", "MAUS::MCEvent", &mcevent, 64000, 10); 
    // Use MAUS internal routines to generate a ROOT streamer. We are here
    // accessing the Spill tree which contains DAQ output data
    // Other trees are e.g. JobHeader (contains Job information), RunHeader, etc
    std::vector<std::string>::iterator filename;
    Collection muonSet;
    Collection pionSet;
    Collection electronSet;
    Collection USrefSet;
    int TOF0count=0, TOF1count=0, TOF2count=0;
    for (filename = filelist.begin(); filename < filelist.end(); filename++){
        // irstream infile((*filename).c_str(), "Spill");
        TFile infile((*filename).c_str());
        TTree* Spill = (TTree*)infile.Get("Spill");
        Spill->SetBranchAddress("data", &data);
        // Iterate over all events
        std::cout<<"Considering file "<<(*filename).c_str()<<std::endl;
        std::string string = (*filename).c_str();
        std::size_t founds = string.find_last_of("/");
        std::size_t founde = string.find_last_of("/", founds-1);
        std::string runnumber = string.substr(founds, founde);
        for ( int i=0; i<Spill->GetEntries(); i++){
            Spill->GetEntry(i);
            MAUS::Spill* spill = data->GetSpill();
            // Iterate over all events; top level event in the "Spill" tree
            // corresponds to a daq_event; e.g. "start_of_burst" (spill start
            //  signal) or "physics_event" (spill data)
            if (spill != NULL && spill->GetDaqEventType() == "physics_event") {
                // Each recon event corresponds to a particle trigger; data in the
                // recon event should all have the same trigger
                if(type.find("recon")!=std::string::npos){
                    /* RunNumber = spill->GetRunNumber(); */
                    RunNumber = stoi(runnumber.substr(0,4));
                    SpillNumber = spill->GetSpillNumber();
                }
                for (size_t ii = 0;  ii < spill->GetReconEvents()->size(); ++ii) {
                    // TOF event holds TOF information for the particle trigger
                    MAUS::TOFEvent* tof_event =
                        (*spill->GetReconEvents())[ii]->GetTOFEvent();
                    MAUS::TOFEventDigit digits;
                    if(type.find("mc")!=std::string::npos){
                        mcevent    = (*spill->GetMCEvents())[ii];
                        // if(fabs(mcevent->GetPrimary()->GetParticleId()) != 13) continue;
                    }
                    if (tof_event != NULL)
                        digits = tof_event->GetTOFEventDigit();
                    else
                        continue;
                    // std::cout<<"Found TOF digits at event "<<i<<std::endl;
                    // Iterate over TOF digits
                    bool plane0digit = false;
                    bool plane1digit = false;
                    for (size_t j = 0; j < digits.GetTOF0DigitArray().size(); ++j) {
                        MAUS::TOFDigit tof0_digit = digits.GetTOF0DigitArray()[j];
                        // get the slab where the digit was registered
                        if (tof0_digit.GetPlane() == 0) {
                            plane0digit = true;
                        } else {
                            plane1digit = true;
                        }
                    }
                    if( plane0digit && plane1digit ) TOF0count++;
                    plane0digit = false;
                    plane1digit = false;
                    for (size_t j = 0; j < digits.GetTOF1DigitArray().size(); ++j) {
                        MAUS::TOFDigit tof1_digit = digits.GetTOF1DigitArray()[j];
                        // get the slab where the digit was registered
                        if (tof1_digit.GetPlane() == 0) {
                            plane0digit = true;
                        } else {
                            plane1digit = true;
                        }
                    }
                    if( plane0digit && plane1digit ){
                        TOF1count++;
                    }

                    plane0digit = false;
                    plane1digit = false;
                    for (size_t j = 0; j < digits.GetTOF2DigitArray().size(); ++j) {
                        MAUS::TOFDigit tof2_digit = digits.GetTOF2DigitArray()[j];
                        // get the slab where the digit was registered
                        if (tof2_digit.GetPlane() == 0) {
                            tof2_digits_0_hist.Fill(tof2_digit.GetSlab());
                            plane0digit = true;
                        } else {
                            tof2_digits_1_hist.Fill(tof2_digit.GetSlab());
                            plane1digit = true;
                        }
                    }
                    if( plane0digit && plane1digit ) TOF2count++;
                }

                /* if(type.find("sim")!=std::string::npos){ */
                    for (size_t ij=0; ij<spill->GetMCEvents()->size(); ij++){
                        MAUS::TOFHitArray* tofhit = (*spill->GetMCEvents())[ij]->GetTOFHits();
                        for (int i=0;i<tofhit->size();i++) {
                            if (tofhit->at(i).GetPosition().Z() - 12929 < 20 and tofhit->at(i).GetPosition().Z() - 12929 > -20){        
                                mcevent = (*spill->GetMCEvents())[ij];
                                mctree->Fill();
                            }
                        }
                        // find particles that pass through TOF2 regardless of PID
                        int nverthits = (*spill->GetMCEvents())[ij]->GetVirtualHits()->size();
                        bool tof2trigger=false;
                        for(int k=0; k<nverthits; k++){
                            float tof0;
                            float tof1;
                            for (int i=0;i<tofhit->size();i++) {
                                if (tofhit->at(i).GetPosition().Z() - 5287 < 20 and tofhit->at(i).GetPosition().Z() - 5287 > -20){  
                                    tof0 = tofhit->at(i).GetTime();   
                                }
                                if (tofhit->at(i).GetPosition().Z() - 12929 < 20 and tofhit->at(i).GetPosition().Z() - 12929 > -20){        tof1 = tofhit->at(i).GetTime();
                                }
                            }
                            if ((tof1-tof0) > TOF_upper_limit || (tof1-tof0) < TOF_lower_limit) { 
                                    std::vector<double> projTheta;
                                    double z = ((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k).GetPosition().z();
                                    double x = ((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k).GetPosition().x();
                                    double y = ((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k).GetPosition().y();
                                    double pz = ((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k).GetMomentum().z();
                                    double px = ((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k).GetMomentum().x();
                                    double py = ((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k).GetMomentum().y();
                                    double dXdz = px/pz;
                                    double dYdz = py/pz;
                                if (k == 67) {

                                    double dp = 0;
                                    double pmag = sqrt(px*px + py*py + pz*pz);
                                    if(k>0){
                                        double dpx = -((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k-1).GetMomentum().x() + px;
                                        double dpy = -((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k-1).GetMomentum().y() + py;
                                        double dpz = -((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k-1).GetMomentum().z() + pz;
                                        dp = sqrt(dpx*dpx + dpy*dpy + dpz*dpz) *  pz / pmag;
                                    }

                                    double USnorm = 1./sqrt(1 + dXdz*dXdz + dYdz*dYdz);
                                    u[0] = dXdz*USnorm;
                                    u[1] = dYdz*USnorm;
                                    u[2] = USnorm;
                                    /* if (u[0]!=0) std::cout << u[0] << std::endl; */
                                    TVector3 s(0, -1, 0);
                                    double snorm = 1./sqrt(s[2]*s[2] + s[1]*s[1] + s[0]*s[0]);
                                    TVector3 s_(s[0]*snorm,s[1]*snorm,s[2]*snorm);
                                    TVector3 v = s.Cross(u);
                                    double vnorm = 1./sqrt(v[2]*v[2] + v[1]*v[1] + v[0]*v[0]); 
                                    v_[0] = v[0]*vnorm;
                                    v_[1] = v[1]*vnorm;
                                    v_[2] = v[2]*vnorm;
                                    w = v_.Cross(u);
                                    double Wnorm  = 1./sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
                                    w[0] *= Wnorm;
                                    w[1] *= Wnorm;
                                    w[2] *= Wnorm;
                                    double abspos = 16952.5;
                                    double zdiff = (16952.5 - z);
                                    double USXproj = (dXdz) * zdiff + x;
                                    double USYproj = (dYdz) * zdiff + y;
                                    double z0 = 19948.8;
                                    double phi = atan2(dYdz, dXdz);
                                    zdiff = (z0 - 16952.5);
                                    USXproj = dXdz * zdiff + USXproj;
                                    USYproj = dYdz * zdiff + USYproj;
                                    double zdi = 13620-z;
                                    double USXdiff = (dXdz) * zdi + x;
                                    double USYdiff = (dYdz) * zdi + y;
                                }
                                    if (u[0]!=0) std::cout << u[0] << std::endl;
                                    if ((k == 69) && w[0]!=0 && sqrt(USXproj*USXproj + USYproj*USYproj)<90 && sqrt(USXdiff*USXdiff + USYdiff*USYdiff)<90) {
                                    /* std::cout << u[0] << std::endl; */


                                        double DSnorm = 1./sqrt(1 + dXdz*dXdz + dYdz*dYdz);
                                        TVector3 d(dXdz*DSnorm, dYdz*DSnorm, DSnorm);
                                        double hold = (atan(d.Dot(w)/d.Dot(u)));
                                        projTheta.push_back(hold);
                                        projTheta.push_back(atan(d.Dot(v_)/d.Dot(u)));
                                        /* std::cout << projTheta[0] << std::endl; */
                                        if (sqrt(projTheta[0]*projTheta[0]+projTheta[1]*projTheta[1])<0.85){
                                            MC_theta_y.Fill(projTheta[1]);
                                                std::cout << "bndgnbdgfndfn" << std::endl;
                                            MC_theta_x.Fill(projTheta[0]);
                                            if ((*spill->GetReconEvents())[ij]->GetSciFiEvent()->scifitracks().size() == 2) {
                                                recon_theta_y.Fill(projTheta[1]);
                                                recon_theta_x.Fill(projTheta[0]);
                                                std::cout << "bndgnbdgfndfn" << std::endl;
                                            }
                                        }
                                        u.Clear();
                                        w.Clear();
                                        d.Clear();
                                    }
                            }


                                    
                                        int pid = ((*spill->GetMCEvents())[ij]->GetVirtualHits())->at(k).GetParticleId();
                                        if ( k == 0 ){
                                            Vars tmp;
                                            if (abs(pid) == 13)
                                                muonSet.append_instance(tmp);
                                            else if(abs(pid) == 211)
                                                pionSet.append_instance(tmp);
                                            else if(abs(pid) == 11)
                                                electronSet.append_instance(tmp);
                                        }

                                        if ( k == 48 ){
                                            Vars tmp;
                                            USrefSet.append_instance(tmp);
                                        }
                                    }
                                }
                                /*
                                   if(false){
                                   mcevent    = (*spill->GetMCEvents())[ij];
                                // primary    = (*spill->GetMCEvents())[ij]->GetPrimary();
                                // sphitarray = (*spill->GetMCEvents())[ij]->GetSpecialVirtualHits();
                                mctree->Fill();
                                }*/
                                /* } */
                        }
        }
    }
    outfile->cd();
    if(type.find("recon")!=std::string::npos){
        tree->Write();
        tree->Print();
    }
    if(type.find("sim")!=std::string::npos){
        MC_theta_x.Write();
        MC_theta_y.Write();
        recon_theta_x.Write();
        recon_theta_y.Write();
    }
    std::cout<<"TOF0 Count:"<<TOF0count<<std::endl;
    std::cout<<"TOF1 Count:"<<TOF1count<<std::endl;
    std::cout<<"TOF2 Count:"<<TOF2count<<std::endl;
    outfile->Close();
}

