#include <TH1F.h>
#include <TClonesArray.h>

namespace histo {

histograms::histograms() {
  initialized = false;
  test = 0;
}

histograms::~histograms() {
  Clear();
}

void histograms::Clear() {
  if (test) {
    delete test;
  }
}

int histograms::Init() {
  if (initialized) {
    return 0;
  }
  test = new TH1D("test1", "test1", 1000, -0.5, 999.5);
  initialized = true;
  return 0;
}

void histograms::Write() {
  if (test) {
    test->Write();
  }
}

bool histograms::FillTest(double pt) {
  if (!initialized) {
    cerr << "histograms instance not initialized\n";
  }
  test->Fill(pt);
  return true;
}

}
