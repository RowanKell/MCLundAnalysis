#include <TSystem.h>

#include "MyProgress.h"

ClassImp(MyProgress);

Int_t MyProgress::fgProgress = 0;

//______________________________________________________________________________
MyProgress::MyProgress():TObject()
{
   fNLoop = 0;
}

//______________________________________________________________________________
MyProgress::~MyProgress()
{
}

//______________________________________________________________________________
void MyProgress::Test(const char *name, Int_t n)
{
   fString = name;
   cout << "  fString = " << name << endl;

   cout << "  fgProgress = " << fgProgress << endl;

   fNLoop = n;
   for (Int_t i=0;i<=n;i++) {
      gSystem->Sleep(2);
      gSystem->ProcessEvents();

      // here is a very long calculation!

      Int_t pc = (Int_t)(100*i/n);
//      if (pc%5 == 0) fgProgress = pc;
      if (pc%5 == 0) {
         fgProgress = pc;
         cerr << "#";
      }//if
   }//for
   cerr << '\n';
   
   cout << "  fgProgress = " << fgProgress << endl;
}
