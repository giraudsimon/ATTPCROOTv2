//Unpacks tpc files from /mnt/daqtesting/e18008_attpc_transfer to /mnt/analysis/e18008/rootMerg/

#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m\e[1m"
#define WHITE   "\033[37m"

// Requires the TPC run number
void unpack_new(int runNumberS800, int runNumberATTPC)
{
  //Load the library for unpacking and reconstruction
  gSystem->Load("libATTPCReco.so");
  gSystem->Load("libS800.so");
  gSystem->Load("libXMLParser.so");

  TStopwatch timer;
  timer.Start();

  //Set the input file
  //TString inputFile = TString::Format("hdf5Files/run_%04d.h5", runNumber);
  TString inputFile = TString::Format("/mnt/daqtesting/e18008_attpc_transfer/h5/run_%04d.h5", runNumberATTPC);
  // TString inputFile = TString::Format("/mnt/simulations/ceclub/giraud/attpc/ATTPCROOTv2/macro/Unpack_HDF5/hdf5Files/run_0002.h5", runNumberATTPC);

  //Set the output file
  // TString outputFile = TString::Format("/mnt/analysis/e18008/rootMerg/run_2%03d_%04d.root", runNumberS800, runNumberATTPC);
  TString outputFile = TString::Format("/mnt/analysis/e18008/rootMerg/run_%04d_%04d.root", runNumberS800, runNumberATTPC);
  // TString outputFile = TString::Format("/projects/ceclub/giraud/attpc/ATTPCROOTv2/macro/Simulation/d2He/Analyis_d2He/merged_run0002_48Ca.root", runNumberS800, runNumberATTPC);

  std::cout << "Unpacking AT-TPC run " << runNumberATTPC << " from: " << inputFile << std::endl;
  std::cout << "Saving in: " << outputFile << std::endl;

  //Set the mapping for the TPC
  TString scriptfile = "e12014_pad_mapping.xml";//"Lookup20150611.xml";
 // TString scriptfile = "Lookup20150611.xml";//"Lookup20150611.xml";
  //TString parameterFile = "ATTPC.testnew.par";
  TString parameterFile = "ATTPC.d2He.par";

  //Set directories
  TString dir = gSystem->Getenv("VMCWORKDIR");
  TString scriptdir = dir + "/scripts/"+ scriptfile;
  TString geomDir   = dir + "/geometry/";
  gSystem -> Setenv("GEOMPATH", geomDir.Data());
  TString digiParFile = dir + "/parameters/" + parameterFile;
  //TString geoManFile  = dir + "/geometry/ATTPC_v1.1.root";
  //TString geoManFile  = dir + "/geometry/ATTPC_He1bar.root";
  TString geoManFile  = dir + "/geometry/ATTPC_d2He_07atm.root";


  //Create a run
  FairRunAna* run = new FairRunAna();
  run -> SetOutputFile(outputFile);
  run -> SetGeomFile(geoManFile);

  //Set the parameter file
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1 -> open(digiParFile.Data(), "in");
  rtdb -> setSecondInput(parIo1);

  //Create the unpacker task
  ATHDFParserTask* HDFParserTask = new ATHDFParserTask();
  HDFParserTask->SetPersistence(kFALSE);//kTRUE // do no write Rawevents
  HDFParserTask->SetATTPCMap(scriptdir.Data());
  HDFParserTask->SetFileName(inputFile.Data());
  HDFParserTask->SetOldFormat(false);
  HDFParserTask->SetTimestampIndex(2);
  HDFParserTask->SetBaseLineSubtraction(kTRUE);

 // auto hash = HDFParserTask->CalculateHash(10,0,2,32);
 // HDFParserTask->SetAuxChannel(hash, "IonCb_32");
 // hash = HDFParserTask->CalculateHash(10,0,2,34);
 // HDFParserTask->SetAuxChannel(hash, "IonCb_34");



 // Double_t x0_corr_tof = 0.101259;
 // Double_t afp_corr_tof = 1177.02;
 // Double_t afp_corr_dE = 61.7607;
 // Double_t x0_corr_dE = -0.0403;
 // Double_t rf_offset = 0.0;
 // Double_t corrGainE1up = 0.6754;
 // Double_t corrGainE1down = 1.;

  std::vector<Double_t> S800par;
  S800par.push_back(0.);//x0_corr_tof
  S800par.push_back(0.);//afp_corr_tof
  S800par.push_back(0.);//afp_corr_dE
  S800par.push_back(0.);//x0_corr_dE
  S800par.push_back(0.);//rf_offset
  S800par.push_back(1.);//corrGainE1up
  S800par.push_back(1.);//corrGainE1down
  // TString S800File = "/projects/ceclub/giraud/s800root_jorge/s800root/files/rootFiles/cal/test-runs800-48Ca-CAL-0001-new.root";
  // TString S800File = TString::Format("/mnt/analysis/e18008/rootS800/cal/run-2%03d-00.root", runNumberS800);
  TString S800File = TString::Format("/mnt/analysis/e18008/rootS800/cal/run-%04d-00.root", runNumberS800);
  ATMergeTask *MergeEvt = new ATMergeTask();
  MergeEvt->SetS800File(S800File);
  MergeEvt->SetPersistence(kTRUE);
  // MergeEvt->SetOptiEvtDelta(150);//old way with fit
  MergeEvt->SetOptiEvtDelta(5);//missed few events during the tests with 100
  MergeEvt->SetGlom(2);
  MergeEvt->SetTsDelta(1272);
  //MergeEvt->SetPIDcut("/projects/ceclub/giraud/attpc/ATTPCROOTv2/macro/Unpack_HDF5/e18008_S800/CUT1.root");
  //MergeEvt->SetPIDcut("/projects/ceclub/giraud/attpc/ATTPCROOTv2/macro/Unpack_HDF5/e18008_S800/CUT2.root");
  MergeEvt->SetPIDcut("CUT1.root");
  MergeEvt->SetParameters(S800par);

  //Create PSA task
  ATPSATask *psaTask = new ATPSATask();
  psaTask -> SetPersistence(kTRUE);//kFALSE
  psaTask -> SetThreshold(50);//50
  psaTask -> SetPSAMode(1); //NB: 1 is ATTPC - 2 is pATTPC - 3 Filter for ATTPC - 4: Full Time Buckets
  psaTask -> SetMaxFinder();

/*
  ATRansacTask *RandTask = new ATRansacTask();
  RandTask -> SetPersistence(kTRUE);
  RandTask -> SetVerbose(kFALSE);
  RandTask -> SetDistanceThreshold(10.0);
  RandTask -> SetTiltAngle(0);
  RandTask -> SetMinHitsLine(30);
  //RandTask -> SetFullMode();
*/
  ATRansacTask *RandTask = new ATRansacTask();
  RandTask ->SetPersistence(kTRUE);
  //RandTask ->SetModelType(1);
  //RandTask ->SetFullMode();
  RandTask->SetTiltAngle(0.0);
  RandTask->SetDistanceThreshold(15.0);
  RandTask->SetMinHitsLine(7);
  RandTask->SetAlgorithm(3); // 0=PCL ransac; 1=Homemade Ransac; 2=Homemade Mlesac; 3=Homemade Lmeds;
  RandTask->SetRanSamMode(3);// 0=Uniform; 1=Gaussian; 2=Weighted; 3=Gaussian+Weighted



  //Add unpacker to the run
  run -> AddTask(HDFParserTask);
  run -> AddTask(MergeEvt);
  run -> AddTask(psaTask);
  run -> AddTask(RandTask);

  run -> Init();

  //Get the number of events and unpack the whole run
  auto numEvents = HDFParserTask->GetNumEvents()/2;

  //numEvents = 1700;//217;
  std::cout << "Unpacking " << numEvents << " events. " << std::endl;

  //return;
   run->Run(0,numEvents);
  //run->Run(0,40);


  std::cout << std::endl << std::endl;
  std::cout<<" Check MergeTask :: S800 events "<<YELLOW<<MergeEvt->GetS800TsSize()<<RESET<<" AT-TPC events "<<YELLOW<<numEvents<<RESET<<"  Merged TS size "<<YELLOW<<
    MergeEvt->GetMergedTsSize()<<RESET<<std::endl;
  if((double)MergeEvt->GetMergedTsSize()/(double)MergeEvt->GetS800TsSize()<1)
    std::cout<<" ! WARNING Ratio (merged/S800) "<<RED<<(double)MergeEvt->GetMergedTsSize()/(double)MergeEvt->GetS800TsSize()<<RESET<<std::endl;
  if((double)MergeEvt->GetMergedTsSize()/(double)MergeEvt->GetS800TsSize()==1)
  std::cout<<" Ratio (merged/S800) "<<GREEN<<(double)MergeEvt->GetMergedTsSize()/(double)MergeEvt->GetS800TsSize()<<RESET<<std::endl;
  std::cout << std::endl;
  std::cout << "Done unpacking events"  << std::endl << std::endl;
  std::cout << "- Output file : " << outputFile << std::endl << std::endl;


  //std::cout << "Done unpacking events"  << std::endl << std::endl;
  //std::cout << "- Output file : " << outputFile << std::endl << std::endl;
  // -----   Finish   -------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------

  return 0;
}
