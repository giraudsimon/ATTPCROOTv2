/*

*/

#include "ATTPCIonDecay.h"

#include "FairPrimaryGenerator.h"
#include "FairRootManager.h"
#include "FairLogger.h"
#include "FairMCEventHeader.h"

#include "FairIon.h"
#include "FairRunSim.h"
#include "FairRunAna.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TObjArray.h"

#include "TRandom.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"
#include "TVirtualMC.h"
#include "TParticle.h"
#include "TClonesArray.h"


#include "FairRunSim.h"
#include "FairIon.h"
#include <iostream>
#include <string>
#include "TParticle.h"

#include "AtStack.h"
#include "AtTpcPoint.h"
#include "ATVertexPropagator.h"

#define amu 931.494

Int_t ATTPCIonDecay::fgNIon = 0;



ATTPCIonDecay::ATTPCIonDecay()
: fMult(0),
fPx(0.), fPy(0.), fPz(0.),
fVx(0.), fVy(0.), fVz(0.),
fIon(0),  fQ(0)
{
  //  cout << "-W- ATTPCIonGenerator: "
  //      << " Please do not use the default constructor! " << endl;
}

// -----   Default constructor   ------------------------------------------
ATTPCIonDecay::ATTPCIonDecay(std::vector<Int_t> *z, std::vector<Int_t> *a, std::vector<Int_t> *q, std::vector<Double_t> *mass,
  Int_t ZB, Int_t AB, Double_t BMass, Double_t SepEne)
  : fMult(0),
  fPx(0.), fPy(0.), fPz(0.),
  fVx(0.), fVy(0.), fVz(0.),
  fIon(0),  fQ(0)
  {

    TParticle* kProton = new TParticle();
    kProton->SetPdgCode(2212);
    TParticle* kNeutron = new TParticle();
    kNeutron->SetPdgCode(2112);
    char buffer[20];

    fgNIon++;
    fMult = a->size();
    fIon.reserve(fMult);

    fPxBeam = 0;
    fPyBeam = 0;
    fPzBeam = 0;
    fZBeam = ZB;
    fABeam = AB;
    fBeamMass = BMass*amu/1000.0;
    fSepEne = SepEne;

    FairRunSim* run = FairRunSim::Instance();
    if ( ! run ) {
      std::cout << "-E- FairIonGenerator: No FairRun instantised!" << std::endl;
      Fatal("FairIonGenerator", "No FairRun instantised!");
    }

    for(Int_t i=0;i<fMult;i++){
      Masses.push_back(mass->at(i));
      FairIon *IonBuff;
      FairParticle *ParticleBuff;
      sprintf(buffer, "Product_Ion_dec%d", i);
      if( a->at(i)!=1  ){
        IonBuff = new FairIon(buffer, z->at(i), a->at(i), q->at(i),0.0,mass->at(i)*amu/1000.0);
        ParticleBuff = new FairParticle("dummyPart",1,1,1.0,0,0.0,0.0);
        fPType.push_back("Ion");
        run->AddNewIon(IonBuff);
        //          std::cout<<" Adding : "<<buffer<<std::endl;

      }else if( a->at(i)==1 && z->at(i)==1  ){
        IonBuff = new FairIon(buffer, z->at(i), a->at(i), q->at(i),0.0,mass->at(i)*amu/1000.0);
        ParticleBuff = new FairParticle(2212,kProton);
        fPType.push_back("Proton");

      }else if( a->at(i)==1 && z->at(i)==0  ){
        IonBuff = new FairIon(buffer, z->at(i), a->at(i), q->at(i),0.0,mass->at(i)*amu/1000.0);
        ParticleBuff = new FairParticle(2112,kNeutron);
        fPType.push_back("Neutron");
      }

      fIon.push_back(IonBuff);
      fParticle.push_back(ParticleBuff);
    }

  }

  // -----   Destructor   ---------------------------------------------------
  ATTPCIonDecay::~ATTPCIonDecay()
  {
    // if (fIon) delete fIon;
  }

  // -----   Public method ReadEvent   --------------------------------------
  Bool_t ATTPCIonDecay::ReadEvent(FairPrimaryGenerator* primGen) {

    Double_t beta;
    Double_t s=0.0;
    Double_t mass_1[10]={0.0};
    Double_t* pMass;
    Double_t M_tot=0;
    TLorentzVector fEnergyImpulsionLab_beam;
    TLorentzVector fEnergyImpulsionLab_Total;
    TLorentzVector fEnergyImpulsionFinal;
    TVector3 fImpulsionLab_beam;
    std::vector<TLorentzVector*> p_vector;
    TGenPhaseSpace event1;

    fPx.clear();
    fPy.clear();
    fPz.clear();

    fPx.resize(fMult);
    fPy.resize(fMult);
    fPz.resize(fMult);

    fIsDecay = kFALSE;

    //AtStack* stack = (AtStack*) gMC->GetStack();

    // === Get ejectile info from previous reaction generator (d2He)

    //std::cout<<"check sizes "<<fMult<<" "<<fPx.size()<<" "<<fIon.size()<<" "<<fParticle.size()<<std::endl;

    //fBeamEnergy = gATVP->GetEnergy()/1000.0;
    fBeamEnergy = gATVP->GetScatterE()/1000.0;
    Double_t ExEject=gATVP->GetScatterEx()/1000.0;//in GeV
    //fPxBeam = gATVP->GetPx();
    //fPyBeam = gATVP->GetPy();
    //fPzBeam = gATVP->GetPz();

    TParticlePDG* thisPart0;
    thisPart0 = TDatabasePDG::Instance()->GetParticle("Product_Ion2");//14N from d2He generator
    int pdgType0 = thisPart0->PdgCode();
    TVector3 ScatP = gATVP->GetScatterP();
    fPxBeam = ScatP.X();
    fPyBeam = ScatP.Y();
    fPzBeam = ScatP.Z();
    std::cout<<" Ejectile info : "<<pdgType0<<" "<<fPxBeam<<" "<<fPyBeam<<" "<<fPzBeam<<" "<<fBeamEnergy<<" "<<ExEject<<std::endl;


    // === Phase Space Calculation

    fImpulsionLab_beam = TVector3(fPxBeam,fPyBeam,fPzBeam);
    fEnergyImpulsionLab_beam = TLorentzVector(fImpulsionLab_beam,fBeamMass+fBeamEnergy+ExEject);

    fEnergyImpulsionLab_Total = fEnergyImpulsionLab_beam;

    s = fEnergyImpulsionLab_Total.M2();
    beta = fEnergyImpulsionLab_Total.Beta();

    for(Int_t i=0;i<fMult;i++){
      M_tot+=Masses.at(i)*amu/1000.0;
      mass_1[i] = Masses.at(i)*amu/1000.0;
      //std::cout<<Masses.at(i)*amu/1000.0<<" "<<M_tot<<" "<<fBeamMass<<" "<<ExEject<<" "<<sqrt(s)<<" "<<(sqrt(s)-M_tot)*1000.0<<std::endl;
    }

    std::cout<<" S : "<<s<<" Pow(M) "<<pow(M_tot,2)<<" "<<gATVP->GetScatterEx()<<" "<<fSepEne<<std::endl;

    if(s>pow(M_tot,2)){
      //if(ExEject*1000.0>fSepEne){
      fIsDecay=kTRUE;
      event1.SetDecay(fEnergyImpulsionLab_Total,fMult, mass_1);
      Double_t weight1 = event1.Generate();

      std::vector<Double_t> KineticEnergy;
      std::vector<Double_t> ThetaLab;

      std::cout<<"  ==== Phase Space Information ==== "<<std::endl;
      for(Int_t i=0;i<fMult;i++){
        p_vector.push_back(event1.GetDecay(i));
        fPx.at(i) = p_vector.at(i)->Px();
        fPy.at(i) = p_vector.at(i)->Py();
        fPz.at(i) = p_vector.at(i)->Pz();
        KineticEnergy.push_back((p_vector.at(i)->E() - mass_1[i]));
        ThetaLab.push_back(p_vector.at(i)->Theta()*180./TMath::Pi());
        std::cout<<" Particle "<<i<<" - TKE (MeV) : "<<KineticEnergy.at(i)*1000<<" - Lab Angle (deg) : "<<ThetaLab.at(i)<<std::endl;
      }

    }// if kinematics condition

    for(Int_t i=0; i<fMult; i++){
      TParticlePDG* thisPart;

      if(fPType.at(i)=="Ion")
      thisPart = TDatabasePDG::Instance()->GetParticle(fIon.at(i)->GetName());
      else if(fPType.at(i)=="Proton")
      thisPart = TDatabasePDG::Instance()->GetParticle(fParticle.at(i)->GetName());
      else if(fPType.at(i)=="Neutron")
      thisPart = TDatabasePDG::Instance()->GetParticle(fParticle.at(i)->GetName());

      if ( ! thisPart ) {
        if(fPType.at(i)=="Ion")
        std::cout << "-W- FairIonGenerator: Ion " << fIon.at(i)->GetName()<< " not found in database!" << std::endl;
        else if(fPType.at(i)=="Proton")
        std::cout << "-W- FairIonGenerator: Particle " << fParticle.at(i)->GetName()<< " not found in database!" << std::endl;
        else if(fPType.at(i)=="Neutron")
        std::cout << "-W- FairIonGenerator: Particle " << fParticle.at(i)->GetName()<< " not found in database!" << std::endl;
        return kFALSE;
      }

      int pdgType = thisPart->PdgCode();

      // Propagate from the vertex of the reaction
      TVector3 d2HeVtx = gATVP->Getd2HeVtx();
      //fVx = gATVP->GetVx();
      //fVy = gATVP->GetVy();
      //fVz = gATVP->GetVz();
      fVx = d2HeVtx.X();
      fVy = d2HeVtx.Y();
      fVz = d2HeVtx.Z();

      //std::cout << "-I- FairIonGenerator: Generating " <<" with mass "<<thisPart->Mass()<<" ions of type "<< fIon.at(i)->GetName() << " (PDG code " << pdgType << ")" << std::endl;
      //std::cout << "    Momentum (" << fPx.at(i) << ", " << fPy.at(i) << ", " << fPz.at(i)
      //<< ") Gev from vertex (" << fVx << ", " << fVy
      //<< ", " << fVz << ") cm" << std::endl;

      if(fIsDecay){
        primGen->AddTrack(pdgType, fPx.at(i), fPy.at(i), fPz.at(i), fVx, fVy, fVz);
      }

    }

    gATVP->IncDecayEvtCnt();

    return kTRUE;

  }


  ClassImp(ATTPCIonDecay)
