#include <TH1F.h>
#include <TClonesArray.h>

namespace histo {

class histograms {
 public:
  histograms();
  ~histograms();

  //deletes histograms
  void Clear();
  //creates histograms
  int Init();
  //writes histograms to root file
  void Write();

  //getter
  TH1D* GetTest() {return test;}

  bool FillKinematics(double pt);

 private:
  TH1D* test;

};

}
