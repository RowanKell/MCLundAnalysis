/******************************************************************************
* macro to test MyProgress                                                      *
******************************************************************************/
//
///////////////////////////
//
// in new root session(s):
//     > .L macroMyProgress.C
//     > Init() 
//     > Test("MyTest", 300)
//
///////////////////////////

//______________________________________________________________________________
void Init()
{
// load libraries
   gSystem->Load("MyProgress_cxx.so");
}//Init

//______________________________________________________________________________
void Test(const char *name, Int_t n)
{
   TFile *file = 0;
   file = new TFile("MyProgress.root", "RECREATE", "Test MyProgress");

   MyProgress *progress = new MyProgress();

   progress->Test(name, n);
   progress->Write("", TObject::kWriteDelete);

   delete progress;
   delete file;
}//Test
