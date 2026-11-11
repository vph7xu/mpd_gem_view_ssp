// gem_dqc_summary.C
//
// Per-detector summary of an mpd_gem_view_ssp data-quality-check ROOT file.
// One page per detector, X and Y planes drawn together (X = blue, Y = red).
//
//   Row 1  ADC     : cluster peak ADC | cluster sum ADC | strip ADC | X-Y charge correlation
//   Row 2  Timing  : strip mean time  | max time sample | X-Y seed time correlation | fired strips
//   Row 3  Position: 1D cluster position | 2D hit map | cluster size | clusters per event
//
// Usage:
//   root -l -b -q 'gem_dqc_summary.C("cluster_0_..._data_quality_check.root")'
//   root -l -b -q 'gem_dqc_summary.C("in.root","out.pdf","GEM0,GEM1,uRWELL0,uRWELL1,GEM2,GEM3")'

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TColor.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace dqc {

const int CX = kBlue + 1;   // X plane
const int CY = kRed + 1;    // Y plane

TFile *gF = nullptr;

TH1 *Get(const char *name)
{
   TH1 *h = (TH1 *)gF->Get(name);
   if (!h || h->GetEntries() <= 0) return nullptr;
   h = (TH1 *)h->Clone(Form("%s_c", name));
   h->SetDirectory(nullptr);
   h->UseCurrentStyle();          // drop axis settings stored in the file
   return h;
}

// Range containing the central (1-2q) of the content of up to two histograms, padded
void QRange(TH1 *a, TH1 *b, double &lo, double &hi, double q = 0.005, double pad = 0.15)
{
   lo = 1e30; hi = -1e30;
   double p[2] = {q, 1 - q}, x[2];
   for (TH1 *h : {a, b}) {
      if (!h || h->Integral() <= 0) continue;
      h->GetQuantiles(2, x, p);
      lo = std::min(lo, x[0]);
      hi = std::max(hi, x[1]);
   }
   if (lo > hi) { lo = 0; hi = 1; return; }
   double d = (hi - lo) * pad;
   lo -= d; hi += d;
}

// Range of non-empty bins of up to two histograms
void NZRange(TH1 *a, TH1 *b, double &lo, double &hi, double pad = 0.03)
{
   lo = 1e30; hi = -1e30;
   for (TH1 *h : {a, b}) {
      if (!h) continue;
      int f = h->FindFirstBinAbove(0), l = h->FindLastBinAbove(0);
      if (f < 1) continue;
      lo = std::min(lo, h->GetXaxis()->GetBinLowEdge(f));
      hi = std::max(hi, h->GetXaxis()->GetBinUpEdge(l));
   }
   if (lo > hi) { lo = 0; hi = 1; return; }
   double d = (hi - lo) * pad;
   lo -= d; hi += d;
}

// Most probable value (skips first bin, light smoothing)
double MPV(TH1 *h)
{
   int best = 2; double bv = -1;
   for (int i = 2; i <= h->GetNbinsX() - 1; i++) {
      double v = h->GetBinContent(i - 1) + h->GetBinContent(i) + h->GetBinContent(i + 1);
      if (v > bv) { bv = v; best = i; }
   }
   return h->GetBinCenter(best);
}

enum Stat { kMean, kPeak };

// Overlay X and Y (and optional dashed extras) in the current pad
void Overlay(TH1 *hx, TH1 *hy, const char *title, const char *xtitle,
             Stat st = kMean, bool norm = false, bool logy = false,
             double xlo = 0, double xhi = 0,
             TH1 *ex = nullptr, TH1 *ey = nullptr)
{
   gPad->SetLogy(logy);
   std::vector<TH1 *> all;
   for (TH1 *h : {hx, hy, ex, ey}) if (h) all.push_back(h);
   if (all.empty()) return;

   if (norm) for (TH1 *h : all) if (h->Integral() > 0) h->Scale(1. / h->Integral());

   double ymax = 0;
   for (TH1 *h : all) {
      if (xhi > xlo) h->GetXaxis()->SetRangeUser(xlo, xhi);
      ymax = std::max(ymax, h->GetMaximum());
   }
   TH1 *first = all[0];
   first->SetTitle(Form("%s;%s;%s", title, xtitle, norm ? "fraction" : "counts"));
   first->SetMaximum(logy ? ymax * 8 : ymax * 1.30);
   if (logy) first->SetMinimum(norm ? 1e-5 : 0.5);
   else first->SetMinimum(0);

   if (hx) { hx->SetLineColor(CX); hx->SetLineWidth(2); }
   if (hy) { hy->SetLineColor(CY); hy->SetLineWidth(2); }
   if (ex) { ex->SetLineColor(CX); ex->SetLineStyle(2); ex->SetLineWidth(1); }
   if (ey) { ey->SetLineColor(CY); ey->SetLineStyle(2); ey->SetLineWidth(1); }

   first->Draw("hist");
   for (size_t i = 1; i < all.size(); i++) all[i]->Draw("hist same");

   TLegend *lg = new TLegend(0.50, 0.72, 0.94, 0.88);
   lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
   auto lab = [&](TH1 *h, const char *p) {
      double n = h->GetEntries();
      if (st == kPeak) return TString::Format("%s: %.0f, MPV %.0f", p, n, MPV(h));
      return TString::Format("%s: %.0f, #mu %.1f, #sigma %.1f", p, n, h->GetMean(), h->GetRMS());
   };
   if (hx) lg->AddEntry(hx, lab(hx, "X"), "l");
   if (hy) lg->AddEntry(hy, lab(hy, "Y"), "l");
   if (ex) lg->AddEntry(ex, "X matched", "l");
   if (ey) lg->AddEntry(ey, "Y matched", "l");
   lg->Draw();
}

// Weighted correlation coefficient and means from a 2D histogram
void Corr(TH2 *h, double &rho, double &mx, double &my)
{
   rho = h->GetCorrelationFactor();
   mx = h->GetMean(1);
   my = h->GetMean(2);
}

void Map2D(TH2 *h, const char *title, const char *xt, const char *yt,
           double xlo, double xhi, double ylo, double yhi,
           bool logz = true, bool diag = false, const char *note = nullptr)
{
   gPad->SetLogz(logz);
   gPad->SetRightMargin(0.15);
   if (!h) {
      TLatex t; t.SetNDC(); t.SetTextAlign(22); t.DrawLatex(0.5, 0.5, "empty");
      return;
   }
   h->SetTitle(Form("%s;%s;%s", title, xt, yt));
   h->GetXaxis()->SetRangeUser(xlo, xhi);
   h->GetYaxis()->SetRangeUser(ylo, yhi);
   if (logz) h->SetMinimum(1);
   h->Draw("colz");
   if (diag) {
      double lo = std::max(xlo, ylo), hi = std::min(xhi, yhi);
      TLine *l = new TLine(lo, lo, hi, hi);
      l->SetLineStyle(2); l->SetLineColor(kBlack); l->Draw();
   }
   if (note) {
      TPaveText *pt = new TPaveText(0.19, 0.80, 0.70, 0.88, "NDC");
      pt->SetFillColor(kWhite); pt->SetFillStyle(1001); pt->SetBorderSize(1);
      pt->SetTextAlign(12); pt->SetTextSize(0.042); pt->SetTextFont(42);
      pt->AddText(note);
      pt->Draw();
   }
}

void DetectorPage(int L, const TString &name, const TString &runlab,
                  const TString &out, bool first, bool last)
{
   auto G = [&](const char *fmt) { return Get(Form(fmt, L)); };

   TH1 *adcx = G("h_raw_cluster_adc_x_layer%d"), *adcy = G("h_raw_cluster_adc_y_layer%d");
   TH1 *madx = G("h_matched_cluster_adc_x_layer%d"), *mady = G("h_matched_cluster_adc_y_layer%d");
   TH1 *sadx = G("h_raw_cluster_sum_adc_x_layer%d"), *sady = G("h_raw_cluster_sum_adc_y_layer%d");
   TH1 *stx = G("h_raw_xstrip_adc_layer%d"), *sty = G("h_raw_ystrip_adc_layer%d");
   TH2 *qc = (TH2 *)G("h_raw_charge_correlation_layer%d");
   TH1 *tx = G("h_raw_strip_mean_time_x_layer%d"), *ty = G("h_raw_strip_mean_time_y_layer%d");
   TH1 *tbx = G("h_raw_xstrip_maxtimebin_layer%d"), *tby = G("h_raw_ystrip_maxtimebin_layer%d");
   TH2 *tc = (TH2 *)G("h_raw_seed_strip_mean_time_corr_layer%d");
   TH1 *fsx = G("h_raw_fired_strip_det%d_plane0"), *fsy = G("h_raw_fired_strip_det%d_plane1");
   TH1 *px = G("h_raw_cluster_pos_x_layer%d"), *py = G("h_raw_cluster_pos_y_layer%d");
   TH2 *pc = (TH2 *)G("h_raw_pos_correlation_layer%d");
   TH1 *csx = G("h_raw_cluster_size_x_layer%d"), *csy = G("h_raw_cluster_size_y_layer%d");
   TH1 *cmx = G("h_raw_x_cluster_multiplicity_layer%d"), *cmy = G("h_raw_y_cluster_multiplicity_layer%d");

   // 1 mm bins for display (raw binning is 1/3 mm)
   if (px) px->Rebin(3);
   if (py) py->Rebin(3);
   if (pc) pc->Rebin2D(3, 3);

   TCanvas *c = new TCanvas(Form("c_layer%d", L), name, 2000, 1150);
   TPad *hdr = new TPad("hdr", "", 0, 0.94, 1, 1);
   TPad *body = new TPad("body", "", 0, 0, 1, 0.94);
   hdr->Draw(); body->Draw();

   hdr->cd();
   TLatex t; t.SetNDC();
   t.SetTextSize(0.55); t.SetTextFont(62); t.DrawLatex(0.01, 0.30, name);
   t.SetTextSize(0.38); t.SetTextFont(42);
   t.DrawLatex(0.16, 0.32, Form("%s   |   X clusters %.0f   |   Y clusters %.0f   |   X-Y matched %.0f",
                                runlab.Data(), adcx ? adcx->GetEntries() : 0.,
                                adcy ? adcy->GetEntries() : 0., qc ? qc->GetEntries() : 0.));
   t.SetTextColor(CX); t.DrawLatex(0.88, 0.32, "X plane");
   t.SetTextColor(CY); t.DrawLatex(0.94, 0.32, "Y plane");

   body->cd();
   body->Divide(4, 3, 0.002, 0.002);
   double lo, hi, lo2, hi2;

   // ---- ADC row
   body->cd(1);
   QRange(adcx, adcy, lo, hi, 0.0, 0.0); QRange(adcx, adcy, lo2, hi, 0.005, 0.0);
   lo = 0;
   Overlay(adcx, adcy, "Cluster peak ADC (dashed: X-Y matched)", "cluster peak ADC",
           kPeak, false, false, lo, hi, madx, mady);
   double qmax = hi;

   body->cd(2);
   QRange(sadx, sady, lo2, hi, 0.005, 0.0);
   Overlay(sadx, sady, "Cluster sum ADC", "cluster sum ADC", kPeak, false, false, 0, hi);

   body->cd(3);
   NZRange(stx, sty, lo, hi);
   Overlay(stx, sty, "Strip ADC (strips in clusters)", "strip ADC", kMean, false, true, 0, hi);

   body->cd(4);
   if (qc) {
      double rho, mx, my; Corr(qc, rho, mx, my);
      Map2D(qc, "Charge correlation (matched)", "X cluster ADC", "Y cluster ADC",
            0, qmax, 0, qmax, true, true,
            Form("#rho = %.2f,  #LTY#GT/#LTX#GT = %.2f", rho, my / mx));
   } else Map2D(nullptr, "Charge correlation", "", "", 0, 1, 0, 1);

   // ---- Timing row
   body->cd(5);
   QRange(tx, ty, lo, hi);
   double tlo = lo, thi = hi;
   Overlay(tx, ty, "Strip mean time", "strip mean time [ns]", kMean, false, false, lo, hi);

   body->cd(6);
   Overlay(tbx, tby, "Time sample of strip max ADC", "time sample", kMean, true, false);

   body->cd(7);
   if (tc) {
      double rho, mx, my; Corr(tc, rho, mx, my);
      Map2D(tc, "Seed-strip time correlation (matched)", "X seed mean time [ns]",
            "Y seed mean time [ns]", tlo, thi, tlo, thi, true, true,
            Form("#rho = %.2f,  #LTt_{Y}-t_{X}#GT = %+.1f ns", rho, my - mx));
   } else Map2D(nullptr, "Seed-strip time correlation", "", "", 0, 1, 0, 1);

   body->cd(8);
   NZRange(fsx, fsy, lo, hi);
   Overlay(fsx, fsy, "Fired strips per event", "fired strips", kMean, false, true, lo, hi);

   // ---- Position row
   body->cd(9);
   NZRange(px, py, lo, hi);
   double plo = lo, phi = hi;
   Overlay(px, py, "Cluster position", "cluster position [mm]", kMean, false, false, lo, hi);

   body->cd(10);
   Map2D(pc, "Hit map, X-Y matched clusters", "X [mm]", "Y [mm]", plo, phi, plo, phi, false);

   body->cd(11);
   NZRange(csx, csy, lo, hi);
   Overlay(csx, csy, "Cluster size", "strips per cluster", kMean, true, false, lo, hi);

   body->cd(12);
   NZRange(cmx, cmy, lo, hi);
   Overlay(cmx, cmy, "Clusters per event", "cluster multiplicity", kMean, true, true, lo, hi);

   c->cd();
   TString opt = first ? (last ? "" : "(") : (last ? ")" : "");
   c->Print(out + opt, Form("Title:%s", name.Data()));
}

} // namespace dqc

void gem_dqc_summary(const char *infile, const char *outfile = "", const char *names = "")
{
   using namespace dqc;

   gStyle->SetOptStat(0);
   gStyle->SetPalette(kRainBow);
   gStyle->SetNumberContours(99);
   gStyle->SetTitleFontSize(0.06);
   gStyle->SetTitleX(0.12); gStyle->SetTitleAlign(13);
   gStyle->SetTitleBorderSize(0); gStyle->SetTitleFillColor(0);
   gStyle->SetPadLeftMargin(0.17); gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadTopMargin(0.11); gStyle->SetPadBottomMargin(0.14);
   gStyle->SetLabelSize(0.045, "xyz");
   gStyle->SetTitleSize(0.05, "xyz");
   gStyle->SetTitleOffset(1.15, "x"); gStyle->SetTitleOffset(1.75, "y");
   gStyle->SetPadGridX(true); gStyle->SetPadGridY(true);
   gStyle->SetGridColor(kGray); gStyle->SetGridStyle(3);
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   gStyle->SetLegendFont(42);
   TH1::AddDirectory(false);

   gF = TFile::Open(infile);
   if (!gF || gF->IsZombie()) { printf("cannot open %s\n", infile); return; }

   std::vector<int> layers;
   for (int L = 0; L < 64; L++)
      if (gF->Get(Form("h_raw_cluster_adc_x_layer%d", L))) layers.push_back(L);
   if (layers.empty()) { printf("no h_raw_cluster_adc_x_layer* histograms found\n"); return; }

   std::vector<TString> nm;
   TObjArray *tok = TString(names).Tokenize(",");
   for (size_t i = 0; i < layers.size(); i++) {
      if (tok && (int)i < tok->GetEntries())
         nm.push_back(((TObjString *)tok->At(i))->GetString().Strip(TString::kBoth));
      else
         nm.push_back(Form("Layer %d", layers[i]));
   }

   TString base = gSystem->BaseName(infile);
   TString runlab = base;
   runlab.ReplaceAll("cluster_0_", "");
   runlab.ReplaceAll(".root_data_quality_check.root", "");
   runlab.ReplaceAll("_data_quality_check.root", "");
   runlab.ReplaceAll(".root", "");

   TString out = outfile;
   if (out == "") {
      out = infile;
      out.ReplaceAll(".root_data_quality_check.root", "");
      out.ReplaceAll(".root", "");
      out += "_summary.pdf";
   }

   for (size_t i = 0; i < layers.size(); i++)
      DetectorPage(layers[i], nm[i], runlab, out, i == 0, i + 1 == layers.size());

   printf("wrote %s (%zu pages)\n", out.Data(), layers.size());
}
