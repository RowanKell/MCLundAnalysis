#ifndef MYPROGRESS
#define MYPROGRESS

//////////////////////////////////////////////////////////////////////////
// MyProgress                                                           //
//////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "TObject.h"
#include "TString.h"

class MyProgress : public TObject
{
private:
   Int_t        fNLoop;
   TString      fString;
   static Int_t fgProgress; //! indicator for progress of 

public:
   MyProgress();
   virtual ~MyProgress();
   virtual void Test(const char *name, Int_t n);
   
  ClassDef(MyProgress,1)  //MyProgress
};

#endif
