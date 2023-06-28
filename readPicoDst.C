#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;
class StTofMatchMaker;

StChain *chain;
void readPicoDst(const Char_t *inputFile = "test.list", TString JobIdName = "tpc")
{
    Int_t nEvents = 15000000;

    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();

    gSystem->Load("StUtilities");
    gSystem->Load("StPicoEvent");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StEmbeddingMaker");

    chain = new StChain();

    TString Name = JobIdName;
    Name.Append(".root");
    StPicoDstMaker *picoMaker = new StPicoDstMaker(2, inputFile, "picoDst");
    StEmbeddingMaker *anaMaker = new StEmbeddingMaker("ana", picoMaker, Name);
    anaMaker->set_target_ID(14); // 15 for antiproton and 14 for proton

    chain->Init();
    cout << "chain->Init();" << endl;
    int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if (nEvents > total)
        nEvents = total;
    for (Int_t i = 0; i < nEvents; i++)
    {
        if (i % 1000 == 0)
            cout << "Working on eventNumber " << i << endl;

        chain->Clear();
        int iret = chain->Make(i);

        if (iret)
        {
            cout << "Bad return code!" << iret << endl;
            break;
        }

        total++;
    }

    cout << "****************************************** " << endl;
    cout << "Work done... now its time to close up shop!" << endl;
    cout << "****************************************** " << endl;
    chain->Finish();
    cout << "****************************************** " << endl;
    cout << "total number of events  " << nEvents << endl;
    cout << "****************************************** " << endl;

    delete chain;
}
