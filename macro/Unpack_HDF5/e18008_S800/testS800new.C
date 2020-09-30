
void testS800new(){

  std::string digiFileName = "run_2016_0026.root";
  TFile* file = new TFile(digiFileName.c_str(),"READ");

  // TTree* tree = (TTree*) file -> Get("caltree");
  TTree* tree = (TTree*) file -> Get("cbmsim");

  Long64_t nEvents = tree -> GetEntries();



  //S800Calc s800cal;

  TH2D *dEup_dEdown = new TH2D ("dEup_dEdown", "dEup_dEdown", 100, 0,250, 100, 0, 250);//test with dE


  S800Calc *s800cal = new S800Calc();

  TBranch *bS800cal = tree->GetBranch("s800cal");
  bS800cal->SetAddress(&s800cal);

  Int_t nb = 0;
  for (Long64_t i=0;i<nEvents;i++) {
  //  for (Long64_t i=0;i<25;i++) {
    bS800cal->GetEntry(i);

    Double_t S800_E1up = s800cal->GetSCINT(0)->GetDEup();
    Double_t S800_E1down = s800cal->GetSCINT(0)->GetDEdown();
    std::cout<<S800_E1up<<"  "<<S800_E1down<<std::endl;


    dEup_dEdown->Fill(S800_E1up,S800_E1down);
    s800cal->Clear();

  }
  dEup_dEdown->Draw("colz");
}
