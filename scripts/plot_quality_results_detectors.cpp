// ROOT: .L plot_quality_results_detectors.cpp+
// plot_quality_results_detectors("quality.root", "quality.pdf");
// Input must be the *_data_quality_check.root histogram file, not GEMCluster tree.
// Arguments: input, output, optional sum-ADC pages, log 1D counts, log 2D counts.
// This only plots saved histograms; it never changes replay cuts or ROOT contents.
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TColor.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace detector_quality {
// EDIT THESE SETTINGS. Limits affect display only; no bins are discarded.
const char* detectorNames[6] = {
    "GEM-#muRWELL tracker (layer 0)", "GEM (layer 1)",
    "DUT (layer 2)", "DUT (layer 3)",
    "GEM-#muRWELL tracker (layer 4)", "GEM-#muRWELL tracker (layer 5)"};
const double positionHalfRange[6] = {60,60,190,190,60,60};
const double peakAdcMax = 2100;
const double clusterSizeMax = 19.5;
const int timeSamples = 9;
const double stripSumMax = 4000;
const double clusterSumMax = 5000;
// Set to "Highest-peak X/Y pair" ONLY for files replayed with c_s = 1.
// Set to "Exactly one X and one Y cluster" only if that replay cut was used.
const char* pairSelection = "X/Y pairs selected during replay";
// Optional annotation, e.g. "ZS: 5 sigma; seed peak >= 50 ADC".
// These are labels only. The macro cannot infer cuts from the histogram file.
const char* cutLabel = "";

struct Panel {
    std::string key, title, x, y;
    double lo, hi, ylo, yhi;
    Panel(std::string k, std::string t, std::string a,
          double l, double h, std::string b="Entries", double yl=0, double yh=0)
      : key(k), title(t), x(a), y(b), lo(l), hi(h), ylo(yl), yhi(yh) {}
};
std::string key(const char* pattern, int d) { return Form(pattern,d); }
void note(double x,double y,const char* s,double size=0.035) {
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.DrawLatex(x,y,s);
}

void page(TCanvas& c,TFile& file,const char* pdf,int d,const char* section,
          const std::vector<Panel>& panels,int pageNo,bool log1,bool log2) {
    c.Clear(); c.cd();
    note(0.04,0.966,detectorNames[d],0.032);
    note(0.04,0.931,section,0.022);
    note(0.04,0.027,Form("%s",gSystem->BaseName(file.GetName())),0.012);
    note(0.93,0.027,Form("%d",pageNo),0.016);
    if (std::string(cutLabel).size()) note(0.04,0.902,cutLabel,0.017);
    // Stack-owned pad removes itself from the canvas on destruction.
    TPad body("dq_body","",0.01,0.05,0.99,0.89);
    body.Draw(); body.Divide(2,3,0.007,0.012);
    std::vector<std::unique_ptr<TH1>> copies;
    for(size_t n=0;n<panels.size();++n) {
        const auto& p=panels[n];
        auto* pad=static_cast<TPad*>(body.cd(n+1));
        pad->SetFillColor(kWhite); pad->SetFrameFillColor(kWhite);
        pad->SetLeftMargin(0.15); pad->SetRightMargin(0.06);
        pad->SetBottomMargin(0.18); pad->SetTopMargin(0.16);
        pad->SetTicks(1,1);
        TH1* original=nullptr; file.GetObject(p.key.c_str(),original);
        note(0.15,0.91,p.title.c_str(),0.053);
        if(!original || original->GetDimension()>2) {
            note(0.15,0.55,"Histogram not available",0.045);
            note(0.08,0.44,p.key.c_str(),0.027);
            std::cerr<<"Missing/unsupported: "<<p.key<<'\n'; continue;
        }
        copies.emplace_back(static_cast<TH1*>(original->Clone(Form("dq_%d_%zu",pageNo,n))));
        auto* h=copies.back().get(); h->SetDirectory(nullptr);
        h->SetTitle(""); h->SetStats(false);
        h->GetXaxis()->SetTitle(p.x.c_str()); h->GetYaxis()->SetTitle(p.y.c_str());
        h->GetXaxis()->SetRangeUser(p.lo,p.hi);
        // Hide the last regular bin and overflow for strip waveform-sum plots.
        if (h->GetDimension() == 1 &&
           (p.key.find("h_raw_xstrip_sum_adc_layer") == 0 ||
            p.key.find("h_raw_ystrip_sum_adc_layer") == 0)) {
        
           const int first = std::max(1, h->GetXaxis()->GetFirst());
           const int last = std::min(
               h->GetNbinsX() - 1,
               h->GetXaxis()->GetLast()
           );
        
           if (first <= last)
               h->GetXaxis()->SetRange(first, last);
        }
	
	for(auto* a:{h->GetXaxis(),h->GetYaxis(),h->GetZaxis()}) {
            a->SetLabelFont(42); a->SetTitleFont(42);
            a->SetLabelSize(0.045); a->SetTitleSize(0.052);
            a->SetNdivisions(505);
        }
        h->GetXaxis()->SetTitleOffset(1.15); h->GetYaxis()->SetTitleOffset(1.25);
        if(h->GetDimension()==2) {
            pad->SetRightMargin(0.19);
            if(p.yhi>p.ylo) h->GetYaxis()->SetRangeUser(p.ylo,p.yhi);
            pad->SetLogz(log2 && h->GetMaximum()>0);
            if(log2 && h->GetMaximum()>0) h->SetMinimum(0.5);
            h->Draw("COLZ");
        } else {
            const int color = (n%2==0) ? TColor::GetColor("#0072B2") : TColor::GetColor("#D55E00");
            h->SetLineColor(color); h->SetLineWidth(2); h->SetFillStyle(0);
            pad->SetLogy(log1 && h->GetMaximum()>0);
            h->SetMinimum(log1 && h->GetMaximum()>0 ? 0.5 : 0.0);
            h->Draw("HIST");
        }
        // Draw annotations after the histogram, which clears the pad on Draw().
        note(0.15,0.91,p.title.c_str(),0.053);
        if(h->GetEntries()==0) note(0.25,0.55,"No entries",0.045);
    }
    c.Modified(); c.Update(); c.Print(pdf);
    // Clear subpads before deleting histogram copies.
    body.Clear();
}
}

void plot_quality_results_detectors(const char* input,const char* output="",
                                    bool sumPages=false,bool logCounts=false,
                                    bool logDensity=true) {
    using namespace detector_quality;
    std::unique_ptr<TFile> f(TFile::Open(input,"READ"));
    if(!f || f->IsZombie()) { std::cerr<<"Cannot open "<<input<<'\n'; return; }
    bool found=false;
    for(int d=0;d<6;++d)
        if(f->Get(key("h_raw_cluster_adc_x_layer%d",d).c_str())) found=true;
    if(!found) {
        std::cerr<<"No expected quality histograms. Use the *_data_quality_check.root file.\n";
        return;
    }
    TString pdf=output;
    if(pdf.IsNull()) { pdf=input; if(pdf.EndsWith(".root")) pdf.Remove(pdf.Length()-5); pdf+="_detectors.pdf"; }
    if(!pdf.EndsWith(".pdf")) { std::cerr<<"Output must end in .pdf\n"; return; }
    TString folder=gSystem->DirName(pdf);
    if(gSystem->AccessPathName(folder) && gSystem->mkdir(folder,true)!=0) {
        std::cerr<<"Cannot create output directory "<<folder<<'\n'; return;
    }
    const bool oldBatch=gROOT->IsBatch(); gROOT->SetBatch(true);
    TStyle saved(*gStyle);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255);
    gStyle->SetCanvasColor(kWhite); gStyle->SetPadColor(kWhite);
    TCanvas canvas("dq_report","Detector quality report",1200,1400);
    canvas.Print(pdf+"[");
    int pages=0;
    for(int d=0;d<6;++d) {
        const double r=positionHalfRange[d];
        page(canvas,*f,pdf,d,"Cluster overview: all surviving clusters; X left, Y right",{
            {key("h_raw_cluster_adc_x_layer%d",d),"X: seed peak ADC","Seed peak ADC",0,peakAdcMax},
            {key("h_raw_cluster_adc_y_layer%d",d),"Y: seed peak ADC","Seed peak ADC",0,peakAdcMax},
            {key("h_raw_cluster_size_x_layer%d",d),"X: cluster size","Strips per cluster",-0.5,clusterSizeMax},
            {key("h_raw_cluster_size_y_layer%d",d),"Y: cluster size","Strips per cluster",-0.5,clusterSizeMax},
            {key("h_raw_cluster_pos_x_layer%d",d),"X: cluster position","X position [mm]",-r,r},
            {key("h_raw_cluster_pos_y_layer%d",d),"Y: cluster position","Y position [mm]",-r,r}
        },++pages,logCounts,logDensity);
        page(canvas,*f,pdf,d,Form("Strips after ZS (top rows); %s (bottom row)",pairSelection),{
            {key("h_raw_xstrip_adc_layer%d",d),"X: strip peak ADC","Strip peak ADC",0,peakAdcMax},
            {key("h_raw_ystrip_adc_layer%d",d),"Y: strip peak ADC","Strip peak ADC",0,peakAdcMax},
            {key("h_raw_xstrip_maxtimebin_layer%d",d),"X: peak time sample","Peak sample index (zero-based)",-0.5,timeSamples-0.5},
            {key("h_raw_ystrip_maxtimebin_layer%d",d),"Y: peak time sample","Peak sample index (zero-based)",-0.5,timeSamples-0.5},
            {key("h_raw_charge_correlation_layer%d",d),"X-Y: peak sample cluster sum ADC","X peak sample cluster sum ADC",0,clusterSumMax,"Y peak sample cluster sum ADC",0,clusterSumMax},
            {key("h_raw_pos_correlation_layer%d",d),"X-Y: cluster position","X position [mm]",-r,r,"Y position [mm]",-r,r}
        },++pages,logCounts,logDensity);
        if(sumPages) page(canvas,*f,pdf,d,"Additional diagnostics: sums over samples; multiplicity as recorded",{
            {key("h_raw_xstrip_sum_adc_layer%d",d),"X: strip waveform sum","Strip sum ADC",0,stripSumMax},
            {key("h_raw_ystrip_sum_adc_layer%d",d),"Y: strip waveform sum","Strip sum ADC",0,stripSumMax},
            {key("h_raw_cluster_sum_adc_x_layer%d",d),"X: cluster peak sample sum","Cluster peak sample sum ADC",0,clusterSumMax},
            {key("h_raw_cluster_sum_adc_y_layer%d",d),"Y: cluster peak sample sum","Cluster peak sample sum ADC",0,clusterSumMax},
            {key("h_raw_x_cluster_multiplicity_layer%d",d),"X: cluster multiplicity","Clusters (zero may not be filled)",-0.5,19.5},
            {key("h_raw_y_cluster_multiplicity_layer%d",d),"Y: cluster multiplicity","Clusters (zero may not be filled)",-0.5,19.5}
        },++pages,logCounts,logDensity);
    }
    canvas.Print(pdf+"]"); saved.Copy(*gStyle); gROOT->SetBatch(oldBatch);
    std::cout<<"Wrote "<<pages<<" pages to "<<pdf<<'\n';
}
