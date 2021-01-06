#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <future>

#include "ROOT/RDataFrame.hxx"

using namespace std;
using namespace ROOT;
using namespace ROOT::RDF;

const double LUMI = 21071.0+38654.0;
const double XSection = 0.01891;

// Used to define root files to read
const string fileName = "../../../data/t_channel/NANOAODJMAR/merged/1.root";

// Define the binning of the different variables to histogram
const map<string, tuple<int, double, double>> ranges = {
  //  "nGoodFatJet"          : {20 , 0 , 20  },
  {"GoodFatJet_pt", {200, 0 , 2000}},
};

// Define the objects and variables to build up the variables to histogram
const map<string, string> objDefinitions = {
  {"GoodFatJet",  "abs(FatJet_eta) < 2.4 && FatJet_pt > 200"},
  //o["nGoodFatJet"]      = "Sum(GoodFatJet)"
  {"GoodFatJet_pt",   "FatJet_pt[GoodFatJet]"},
};


template<class T> double duration(T t0,T t1)
{
  auto elapsed_secs = t1-t0;
  typedef chrono::duration<float> float_seconds;
  auto secs = chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

inline chrono::time_point<std::chrono::steady_clock> now()
{
  return chrono::steady_clock::now();
}

// Book a histogram for a specific variable
RResultPtr<TH1D> bookHistogram(RInterface<ROOT::Detail::RDF::RLoopManager, void> df,
                               string variable, tuple<int, double, double> range_, string weight)
{
  auto histModel = TH1DModel(variable.c_str(), variable.c_str(), get<0>(range_), get<1>(range_), get<2>(range_));
  return df.Histo1D(histModel, variable, weight);
}


// Write a histogram with a given name to the output ROOT file
void writeHistogram(RResultPtr<TH1D> h, string name)
{
  h->SetName(name.c_str());
  h->Write();
}


//void makeHistogramsRDF()
int main()
{
  auto tstart = now();
  ROOT::EnableImplicitMT(8);
  
  // Load dataset
  auto df = RDataFrame("Events", fileName);
  
  string firstKey = objDefinitions.begin()->first;
  string firstValue = objDefinitions.begin()->second;
  
  auto dfAfterCuts = df.Define(firstKey, firstValue);
  bool first = true;
  
  for(auto &[obj, definition] : objDefinitions){
    if(first){
      first = false;
      continue;
    }
    dfAfterCuts = dfAfterCuts.Define(obj, definition);
  }
  
  // Create output file
  TFile *tfile = new TFile("./histograms_cpp.root", "RECREATE");
  tfile->cd();
  
  // Get mean of event weight sums (why?)
  auto dfRun = RDataFrame("Runs", fileName);
  int nGenEvts = dfRun.Mean("genEventSumw_").GetValue();
  
  // Book and fill histograms
  for(auto &[variable, range] : ranges){
    string weight = "genWeight";
    auto hist = bookHistogram(dfAfterCuts, variable, range, weight);
    hist->Scale( XSection*LUMI/nGenEvts );
    writeHistogram(hist, variable);
  }
  
  tfile->Close();
  
  auto elapsed = duration(tstart, now());
  cout<<"Elapsed time: "<<elapsed<<endl;
}


