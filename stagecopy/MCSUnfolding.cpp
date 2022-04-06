#include "MCSAnalysis.h"

// xml parser
#include "libxml/tree.h"
#include "libxml/parser.h"
#include "libxml/xpath.h"
#include "libxml/xpathInternals.h"

#if defined(LIBXML_XPATH_ENABLED) && defined(LIBXML_SAX1_ENABLED)

struct specvals {
  std::string outfile;
  double isMC;
  std::string outfilename;
  std::string dataname;
  std::string emptydataname;
  std::string trainname;
  std::string modelname;
  std::string modelname2;
  std::string geometryname;
  std::string model1;			       
  std::string model2;
  std::string model3;
  std::string material;
  std::string trkreffiname;
  std::string trkreffiemptyname;
  std::map<std::string, double> sys;
  std::map<std::string, double> histlimits;
  double angdef;
  double rot_ang;
  double rot_ang_empty;
  double rot_angX;
  double rot_ang_emptyX;
  double accpt;
  double TOF_ll;
  double TOF_ul;
  double fid_grad;
  double fid_rad;
  int    mode;
}; 

static void print_element_names(xmlNode * a_node, specvals& spec)
{
  xmlNode *cur_node = NULL;
  for (cur_node= a_node; cur_node; cur_node = cur_node->next) { 
    if (cur_node->type == XML_ELEMENT_NODE) {
      const xmlChar* id = xmlCharStrdup("id");
      const xmlChar* nm = xmlCharStrdup("name");
      const xmlChar* vl = xmlCharStrdup("value");
      char* val1 = (char*)xmlGetProp(cur_node, id);
      char* val2 = (char*)xmlGetProp(cur_node, nm);
      char* val3 = (char*)xmlGetProp(cur_node, vl);
      printf("node type: Element, name: %s, with prop id: %s, prop name: %s, prop value: %s\n", 
	     cur_node->name, val1, val2, val3);
      if ( xmlStrEqual(cur_node->name, xmlCharStrdup("file")) ) {
	if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("outfile")) ){
	  spec.outfile = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("datafile")) ){
	  spec.dataname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("trainfile")) ){
	  spec.trainname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("isMC")) ){
	  spec.isMC = std::atof((char*)xmlGetProp(cur_node, nm));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("modelfile")) ){
	  spec.modelname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("modelfile2")) ){
	  spec.modelname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("model1")) ) {
	  spec.model1 = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("model2")) ) {
	  spec.model2 = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("model3")) ) {
	  spec.model3 = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("material")) ) {
	  spec.material = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("geofile")) ) {
	  spec.geometryname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("trkreffiname")) ) {
	  spec.trkreffiname = (char*)xmlGetProp(cur_node, nm);
	} else if ( xmlStrEqual(xmlGetProp(cur_node, id), xmlCharStrdup("trkreffiemptyname")) ) {
	  spec.trkreffiemptyname = (char*)xmlGetProp(cur_node, nm);
	}else {
	  std::cout<<"File designation not recognized. Will ignore "
		   << xmlGetProp(cur_node, id) 
		   << std::endl;
	}
      }

      if (xmlStrEqual(cur_node->name, xmlCharStrdup("cuts")) ) {
        if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("TOF_ll")) ){
	  spec.TOF_ll = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("TOF_ul")) ){
	  spec.TOF_ul = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("fid_grad")) ){
	  spec.fid_grad = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("fid_rad")) ){
	  spec.fid_rad = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("mode")) ){
	  spec.mode = std::atof((char*)xmlGetProp(cur_node, vl));
	} else {
	  std::cout<<"Cut designation not recognized. Will ignore "
		   << xmlGetProp(cur_node, nm) 
		   << std::endl;
	}
      }
      if (xmlStrEqual(cur_node->name, xmlCharStrdup("sys")) ){
	if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("angdef")) ){
	  spec.angdef = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("rot_ang")) ){
	  spec.rot_ang = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("rot_ang_empty")) ){
	  spec.rot_ang_empty = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("rot_angX")) ){
	  spec.rot_angX = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("rot_ang_emptyX")) ){
	  spec.rot_ang_emptyX = std::atof((char*)xmlGetProp(cur_node, vl));
	} else if ( xmlStrEqual(xmlGetProp(cur_node, nm), xmlCharStrdup("accpt")) ){
	  spec.accpt = std::atof((char*)xmlGetProp(cur_node, vl));
	} else spec.sys[(char*)xmlGetProp(cur_node, nm)] = 
	  std::atof((char*)xmlGetProp(cur_node, vl));
      }
      if (xmlStrEqual(cur_node->name, xmlCharStrdup("histlimits")) ){
	spec.histlimits[(char*)xmlGetProp(cur_node, nm)] = 
	  std::atof((char*)xmlGetProp(cur_node, vl));
      }
    }

    print_element_names(cur_node->children, spec);
  }
}

int main(int argc, char* argv[]) { 
  // char* mrd = getenv("MAUS_ROOT_DIR");
  
  std::cout<<"Algorithm to analyse reduced MAUS root\n";
  std::cout<<"trees for multiple scattering measurements\n";
  std::cout<<"and to unfold the distribution\n";
  std::cout<<"\nReceived "<<argc<<" arguements.\n";
  if( argc != 9 && argc != 2){
    std::cout<<"2 arguements are required:\n\n";
    std::cout<<"${exec} ${spec file}\n";
    return -1;
  }

  std::string mode_tree = "";

  specvals spec;

  spec.outfile = "";
  spec.isMC = 0;
  spec.dataname  = "";
  spec.trainname = "";
  spec.modelname = "";
  spec.modelname2 = "";
  spec.geometryname = "";
  spec.trkreffiname = "";
  spec.trkreffiemptyname = "";
  spec.angdef = 0;
  spec.rot_ang = 0;
  spec.rot_ang_empty = 0;
  spec.rot_angX = 0;
  spec.rot_ang_emptyX = 0;
  spec.accpt = -0.2;
  spec.TOF_ll = 27.0;
  spec.TOF_ul = 42.0;
  spec.fid_grad = 0.01;
  spec.fid_rad  = 150.0;
  spec.mode     = 0;

  if ( argc == 2 ){
    xmlInitParser();
    xmlDocPtr doc=NULL;


    assert(argv[1]);

    doc = xmlReadFile(argv[1], NULL, 0);
    if (doc == NULL) { 
      printf("error: could not parse file %s\n", argv[1]);
    }
    xmlNode * root_element = xmlDocGetRootElement(doc);

    print_element_names(root_element, spec);

    xmlFreeDoc(doc);

  }
  else if ( argc == 10 ){
    spec.outfile = argv[1];
    spec.dataname  = argv[2];
    spec.trainname = argv[3];
    spec.TOF_ll    = std::atof(argv[4]);
    spec.TOF_ul    = std::atof(argv[5]);
    spec.fid_rad   = std::atof(argv[6]);
    spec.modelname = argv[7];
    spec.modelname2 = argv[7];
    spec.mode      = std::atoi(argv[8]);
    spec.angdef    = std::atof(argv[9]);
    spec.rot_ang   = std::atof(argv[10]);
  }
    
  if (spec.mode > 0){
    mode_tree = "recon_reduced_tree";
  } else {
    mode_tree = "recon_reduced_tree";
  }
  
  std::cout<<"Writing outfile to "<<spec.outfile<<std::endl; 
  std::cout<<"Reading outfilename at "<<spec.outfilename<<std::endl; 
  std::cout<<"Reading data from "<<spec.dataname<<std::endl; 
  std::cout<<"Reading reference data from "<<spec.trainname<<std::endl; 
  std::cout<<"Reading models from "<<spec.modelname<<std::endl; 
  std::cout<<"Reading models from "<<spec.modelname2<<std::endl; 
  std::cout<<"Use "<<spec.geometryname<<" for propagation\n";
  std::cout<<"\n";
  std::cout<<"angdef " <<spec.angdef<<std::endl;
  std::cout<<"rot_ang " <<spec.rot_ang<<std::endl;
  std::cout<<"rot_ang_empty " <<spec.rot_ang_empty<<std::endl;
  std::cout<<"rot_angX " <<spec.rot_angX<<std::endl;
  std::cout<<"rot_ang_emptyX " <<spec.rot_ang_emptyX<<std::endl;
  std::cout<<"accpt " <<spec.accpt<<std::endl;
  std::cout<<"TOF selection between "
	   <<spec.TOF_ll<<" and "<<spec.TOF_ul<<std::endl;
  std::cout<<"Fiducial selection contained by radius "
	   <<spec.fid_rad<<" with scattering "<<spec.fid_grad<<std::endl;
  std::cout<<"isMC "<<spec.isMC<<std::endl; 
  std::cout<<"mode "<<spec.mode<<std::endl; 
  std::cout<<"\n";

  MCSAnalysis anal("recon_reduced_tree", "Truth_reduced_tree", mode_tree, spec.outfile, spec.histlimits);

  /*
  TFile* data = new TFile(spec.dataname.c_str());
  if(data->IsZombie()){
    std::cout<<"Data file is a zombie. Aborting."<<std::endl;
    return -1;
  }
  TFile* training = new TFile(spec.trainname.c_str());
  if(training->IsZombie()){
    std::cout<<"Training file is a zombie. Aborting."<<std::endl;
    return -1;
  }
  */
  anal.Setangdef(spec.angdef);
  anal.Setrot_ang(spec.rot_ang);
  anal.Setrot_ang_empty(spec.rot_ang_empty);
  anal.Setrot_angX(spec.rot_angX);
  anal.Setrot_ang_emptyX(spec.rot_ang_emptyX);
  anal.Setaccpt(spec.accpt);
  anal.SetTOFLowerLimit(spec.TOF_ll);
  anal.SetTOFUpperLimit(spec.TOF_ul);
  anal.SetRadialLimit(spec.fid_rad);
  anal.SetGradientLimit(spec.fid_grad);
  anal.SetModelFileName(spec.modelname.c_str());
  anal.SetModelName1(spec.model1);
  anal.SetModelName2(spec.model2);
  anal.SetModelName3(spec.model3);
  anal.SetMaterial(spec.material);
  anal.SetisMC(spec.isMC);
  anal.Setmode(spec.mode);
  std::ostringstream strs;
  /*
  strs << spec.TOF_ll;
  std::string str = strs.str();
  string TEN = spec.trkreffiname.c_str()+str+".root";
  */
  anal.SetTrkrEffiName(spec.trkreffiname.c_str());
  anal.SetTrkrEffiEmptyName(spec.trkreffiemptyname.c_str());
  anal.SetParentGeometryFile(spec.geometryname.c_str());
  anal.SetFileName(spec.outfilename.c_str());

  void* dir = gSystem->OpenDirectory(spec.trainname.c_str());
  const char *ent;
  TString fn = spec.trainname.c_str();
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
          if (fn.Contains("ZeroAbs")) {
          //if (fn.Contains("output_7652ZeroAbs200_7600_7629.root")) {
          if (fn.Contains(".root")) {
		  std::cout << fn << std::endl;
		  anal.StageRefData("recon_reduced_tree", fn);
	  }
      }
      fn = spec.trainname.c_str();
  }
  gSystem->FreeDirectory(dir);

  dir = gSystem->OpenDirectory(spec.dataname.c_str());
  fn = spec.dataname.c_str();
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
          if (fn.Contains("LiH2") || fn.Contains("LiH1")) {
	  //if (fn.Contains("output_7790LiH240_7600_7629.root")) {
          //if (fn.Contains("7000")) {
          if (fn.Contains(".root")) {
		  std::cout << fn << std::endl;
		  anal.StageData("recon_reduced_tree", fn);
	  }
      }
      fn = spec.dataname.c_str();
  }
  gSystem->FreeDirectory(dir);

  std::cout << "MC selection" << std::endl;
  dir = gSystem->OpenDirectory(spec.dataname.c_str());
  fn = spec.dataname.c_str();
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
          if (fn.Contains("LiH2") || fn.Contains("LiH1")) {
      //if (fn.Contains("output_7790LiH240_7600_7629.root")) {        
      //if (fn.Contains("700000")) {
          if (fn.Contains(".root")) {
		  std::cout << fn << std::endl;
		  anal.StageMCData("Truth_reduced_tree", fn);
	  }
      }
      fn = spec.dataname.c_str();
  }
  gSystem->FreeDirectory(dir);
/*  
  dir = gSystem->OpenDirectory(spec.trainname.c_str());
  fn = spec.trainname.c_str();
  while (ent = gSystem->GetDirEntry(dir)) {
      fn = fn.Append(ent);
          if (fn.Contains("ZeroAbs")) {
          //if (fn.Contains("700")) {
          if (fn.Contains(".root")) {
		  anal.StageMCData("Truth_reduced_tree", fn);
	  }
      }
      fn = spec.trainname.c_str();
  }
  gSystem->FreeDirectory(dir);
  */
  //anal.GetTree()->AddFileInfoList(anal.GetFileInfo()->CreateListMatching(spec.dataname.c_str()));
  //anal.GetTree()->Add(spec.dataname.c_str());
  
  for(std::map<std::string, double>::iterator it=spec.sys.begin();
      it != spec.sys.end(); ++it){
    anal.SetSysOffset(it->first, it->second);
    std::cout<<"Adding systematic effect "<<it->first
	     <<" with value "<<it->second<<std::endl;
  }
  anal.Execute(spec.mode);

  anal.Write();

  return 0;

}

#else
int main(void) {
  fprintf(stderr, "XPath support not compiled in\n");
  exit(1);
}
#endif
