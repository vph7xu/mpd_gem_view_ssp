// ROOT 6 macro: matching detectors on a 3-column, 2-row page, including nested directories.
// root -l -b -q 'PrintAllHistograms.C("input.root","histograms.pdf")'
// Optional arguments: logY (1D), logZ (2D), skipEmpty.
// EDIT PlotSettings() below to change ranges/styles. Ranges affect display only.
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <array>
#include <TStyle.h>
#include <TColor.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TString.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <regex>
#include <limits>
#include <cmath>

namespace GemPdf {
const double AUTO = std::numeric_limits<double>::quiet_NaN();
struct PlotRange {
    const char *pattern; // Full-match regex on histogram NAME (not directory path).
    int dimension;      // 0 = any, 1 = TH1, 2 = TH2.
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int log;            // -1 = use function argument; 0 = linear; 1 = logarithmic.
};

// ===================== USER SETTINGS =====================
// For TH1: x = coordinate, y = bin content/counts.
// For TH2: x/y = coordinates, z = colour scale (bin content).
// AUTO leaves that endpoint automatic. Later matching rows override earlier
// rows, including AUTO; put specific histogram overrides at the END.
// Ranges estimated from the run-479 PDF, not detector cuts or fit windows.
std::vector<PlotRange> PlotSettings() {
    return {
        // pattern                            dim xmin xmax  ymin ymax  zmin zmax log
        {"h_chi2ndf",                          1, 0,   50,   AUTO,AUTO, AUTO,AUTO,-1},
        {"h_.*cluster_adc_.*",                 1, 0, 2200,   AUTO,AUTO, AUTO,AUTO,-1},
        {"h_raw_[xy]strip_adc_layer.*",         1, 0, 2200,   AUTO,AUTO, AUTO,AUTO,-1},
        {"h_raw_charge_correlation_layer.*",    2, 0, 2200,      0,2200, AUTO,AUTO,-1},
        {"h_raw_[xy]_strip_adc_index_layer.*",  2, AUTO,AUTO,    0,2200, AUTO,AUTO,-1},
        {"h_.*cluster_size_.*",                1,-0.5,17.5, AUTO,AUTO, AUTO,AUTO,-1},
        {"h_raw_size_correlation_layer.*",      2,-0.5,17.5,  -0.5,17.5, AUTO,AUTO,-1},
        {"h_raw_[xy]_cluster_multiplicity_.*",  1,-0.5,17.5, AUTO,AUTO, AUTO,AUTO,-1},
        {"h_raw_fired_strip_det.*",             1,-0.5,35.5, AUTO,AUTO, AUTO,AUTO,-1},
        {"h_.*max_timebin_.*",                 1,-0.5, 8.5, AUTO,AUTO, AUTO,AUTO,-1},
        {"h_raw_[xy]strip_maxtimebin_layer.*",  1,-0.5, 8.5, AUTO,AUTO, AUTO,AUTO,-1},
        {"h_raw_strip_mean_time_.*",            1, 0,  250,   AUTO,AUTO, AUTO,AUTO,-1},
        {"h_raw_seed_strip_mean_time_corr_.*",  2, 0,  250,      0, 250, AUTO,AUTO,-1},
        {"h_[xy]resid_gem.*",                  1,-5,    5,   AUTO,AUTO, AUTO,AUTO,-1},
        {"h_[xy]residue_gem.*",                1,-5,    5,   AUTO,AUTO, AUTO,AUTO,-1},
	{"h_raw_cluster_pos_[xy]_layer2",          1,-190,190,AUTO,AUTO,AUTO,AUTO,-1},
	{"h_raw_cluster_pos_[xy]_layer3",          1,-190,190,AUTO,AUTO,AUTO,AUTO,-1},
	{"h_raw_cluster_pos_[xy]_layer0",          1,-80,80,AUTO,AUTO,AUTO,AUTO,-1},
	{"h_raw_cluster_pos_[xy]_layer1",          1,-80,80,AUTO,AUTO,AUTO,AUTO,-1},
	{"h_raw_cluster_pos_[xy]_layer4",          1,-80,80,AUTO,AUTO,AUTO,AUTO,-1},
	{"h_raw_cluster_pos_[xy]_layer5",          1,-80,80,AUTO,AUTO,AUTO,AUTO,-1},
	{"h_raw_pos_correlation_layer2",           2,-190,190,-190,190,AUTO,AUTO,-1},
	{"h_raw_pos_correlation_layer3",           2,-190,190,-190,190,AUTO,AUTO,-1},
	{"h_raw_pos_correlation_layer0",           2,-70,70,-70,70,AUTO,AUTO,-1},
	{"h_raw_pos_correlation_layer1",           2,-70,70,-70,70,AUTO,AUTO,-1},
	{"h_raw_pos_correlation_layer4",           2,-70,70,-70,70,AUTO,AUTO,-1},
	{"h_raw_pos_correlation_layer5",           2,-70,70,-70,70,AUTO,AUTO,-1},
        // Spatial histograms retain their stored ranges unless overridden.
        {"h_2defficiency_xy_gem.*",             2, AUTO,AUTO,AUTO,AUTO,   0,   1, 0},
        // Examples: uncomment and edit to override a single histogram.
        // {"h_raw_cluster_adc_x_layer2", 1, 0,1800, 0,30000, AUTO,AUTO,0},
        // {"h_xresidue_gem2_tracker_exclusive",1,-2,2, 0,6000,AUTO,AUTO,0},
        // {"h_didhit_xy_gem2",2,-100,100,-100,100,0,100,0},
    };
}
const int LINE_COLOR = kBlack;
const int PALETTE = kRainBow; // Change to kRainBow for ROOT's rainbow palette.
const bool SHOW_STATISTICS = false;
const double COUNT_HEADROOM = 1.20; // Automatic linear TH1 ymax = peak * this.
// =========================================================

PlotRange Settings(const TH1 &h) {
    PlotRange result = {"",0,AUTO,AUTO,AUTO,AUTO,AUTO,AUTO,-1};
    for (const auto &r : PlotSettings())
        if ((!r.dimension || r.dimension == h.GetDimension()) &&
            std::regex_match(h.GetName(), std::regex(r.pattern))) result = r;
    return result;
}
void CoordinateRange(TAxis *a, double lo, double hi) {
    if (std::isnan(lo) && std::isnan(hi)) return;
    const double low = std::isnan(lo) ? a->GetXmin() : lo;
    const double high = std::isnan(hi) ? a->GetXmax() : hi;
    if (high > low) a->SetRangeUser(low, high);
    else std::cerr << "Ignoring invalid coordinate range " << low << ", " << high << '\n';
}
void ContentRange(TH1 &h, double lo, double hi, bool logarithmic) {
    if (!std::isnan(lo) && !std::isnan(hi) && hi <= lo) {
        std::cerr << "Ignoring invalid content range.\n"; return;
    }
    // Clear stored content limits, then apply editable endpoints.
    h.SetMinimum(); h.SetMaximum();
    if (!std::isnan(lo)) {
        if (!logarithmic || lo > 0) h.SetMinimum(lo);
        else std::cerr << "Nonpositive content minimum ignored on logarithmic axis.\n";
    } else if (h.GetDimension() == 1 && !logarithmic && h.GetMinimum() >= 0)
        h.SetMinimum(0);
    if (!std::isnan(hi)) h.SetMaximum(hi);
    else if (h.GetDimension() == 1 && !logarithmic && h.GetMaximum() > 0)
        h.SetMaximum(COUNT_HEADROOM * h.GetMaximum());
}
}

namespace GemPdf {
void Find(TDirectory *dir, const std::string &prefix,
          std::vector<std::string> &paths) {
    // Only the latest saved cycle of each object is exported.
    std::map<std::string, TKey *> latest;
    TIter next(dir->GetListOfKeys());
    while (auto *key = dynamic_cast<TKey *>(next())) {
        auto it = latest.find(key->GetName());
        if (it == latest.end() || key->GetCycle() > it->second->GetCycle())
            latest[key->GetName()] = key;
    }
    for (const auto &item : latest) {
        TClass *cls = TClass::GetClass(item.second->GetClassName());
        if (!cls) continue;
        const std::string path = prefix + item.first;
        if (cls->InheritsFrom(TDirectory::Class())) {
            if (auto *sub = dir->GetDirectory(item.first.c_str()))
                Find(sub, path + "/", paths);
        } else if (cls->InheritsFrom(TH1::Class())) paths.push_back(path);
    }
}
bool Empty(const TH1 &h) {
    // Efficiency maps can have bin contents even when GetEntries() is zero.
    for (int i = 0; i < h.GetNcells(); ++i)
        if (h.GetBinContent(i) != 0 || h.GetBinError(i) != 0) return false;
    return true;
}
// Directory remains part of the key: before/after stages never mix.
struct Page {
    std::string label;
    std::array<std::string,6> paths;
    bool detectors = false;
};
std::vector<Page> Group(const std::vector<std::string> &paths) {
    std::vector<Page> pages;
    std::map<std::string,size_t> indices;
    const std::regex detector("(.*(?:gem|layer|det))([0-5])([^0-9].*|$)");
    for (const auto &path : paths) {
        const auto slash = path.find_last_of('/');
        const std::string directory = slash == std::string::npos ? "" : path.substr(0,slash+1);
        const std::string name = slash == std::string::npos ? path : path.substr(slash+1);
        std::smatch m;
        if (std::regex_match(name,m,detector)) {
            const std::string key = directory + m[1].str() + "{0-5}" + m[3].str();
            auto it = indices.find(key);
            if (it == indices.end()) {
                indices[key] = pages.size();
                Page page; page.label=key; page.detectors=true; pages.push_back(page);
            }
            pages[indices[key]].paths[std::stoi(m[2].str())] = path;
        } else {
            Page page; page.label=path; page.paths[0]=path; pages.push_back(page);
        }
    }
    return pages;
}
void Axis(TAxis *a) {
    a->SetLabelFont(42); a->SetTitleFont(42);
    a->SetLabelSize(0.043); a->SetTitleSize(0.050);
    a->SetLabelOffset(0.008); a->SetNdivisions(505);
}
// Draw literal strings: underscores in names should not become subscripts.
void Text(double x, double y, const std::string &s, double size, int align=11) {
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextAlign(align);
    // Fit long labels onto the page without dropping their text.
    t.SetTextSize(std::min(size, 1.35 / std::max(size_t(1), s.size())));
    t.DrawText(x, y, s.c_str());
}
}

void PrintAllHistograms(const char *input, const char *output="histograms.pdf",
                        bool logY=false, bool logZ=false, bool skipEmpty=false) {
    TString pdf(output);
    if (!pdf.EndsWith(".pdf", TString::kIgnoreCase)) {
        std::cerr << "Output must end in .pdf\n"; return;
    }
    std::unique_ptr<TFile> file(TFile::Open(input, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << input << '\n'; return;
    }
    std::vector<std::string> paths;
    GemPdf::Find(file.get(), "", paths);
    if (paths.empty()) { std::cerr << "No TH1/TH2/TH3 histograms found.\n"; return; }
    TString parent = gSystem->DirName(output);
    if (parent != "." && gSystem->AccessPathName(parent) &&
        gSystem->mkdir(parent, true) != 0) {
        std::cerr << "Cannot create output directory " << parent << '\n'; return;
    }

    // Local style; restore the caller's style after export.
    TStyle *previous = gStyle;
    TStyle style("GemPdfStyle", "Histogram PDF style");
    style.SetOptStat(GemPdf::SHOW_STATISTICS ? 1110 : 0); style.SetOptTitle(0); style.SetOptFit(0);
    style.SetCanvasColor(0); style.SetPadColor(0); style.SetFrameFillColor(0);
    style.SetCanvasBorderMode(0); style.SetPadBorderMode(0);
    style.SetFrameBorderMode(0); style.SetFrameLineWidth(1);
    style.SetPaperSize(26, 20); style.SetEndErrorSize(3);
    style.SetPalette(GemPdf::PALETTE); style.SetNumberContours(100); style.cd();
    TCanvas canvas("gem_pdf_canvas", "Histograms", 1500, 1000);
    canvas.Print((pdf + "[").Data());
    int pages = 0, emptyCount = 0, skipped = 0;
    for (const auto &page : GemPdf::Group(paths)) {
        canvas.Clear(); canvas.cd();
        GemPdf::Text(0.03,0.97,page.label,0.022);
        GemPdf::Text(0.97,0.025,"Page " + std::to_string(pages+1),0.016,31);
        // Keep clones alive until every pad has been printed and cleared.
        std::vector<std::unique_ptr<TH1>> drawn;
        bool anyDrawn=false;
        const int panels=page.detectors ? 6 : 1;
        for (int panel=0; panel<panels; ++panel) {
            canvas.cd();
            const int col=panel%3, row=panel/3;
            const double x1=page.detectors ? col/3.0 : 0;
            const double x2=page.detectors ? (col+1)/3.0 : 1;
            const double y1=page.detectors ? 0.06+(1-row)*0.43 : 0.06;
            const double y2=page.detectors ? y1+0.43 : 0.92;
            auto *pad=new TPad(Form("panel_%d",panel),"",x1,y1,x2,y2);
            pad->Draw(); pad->cd();
            const auto &path=page.paths[panel];
            if (path.empty()) {
                GemPdf::Text(0.5,0.5,"Detector " + std::to_string(panel) + ": unavailable",0.045,22);
                continue;
            }
        auto *original = dynamic_cast<TH1 *>(file->Get(path.c_str()));
        if (!original) { ++skipped; continue; }
        // Retain the name so the editable range rules match this clone.
        std::unique_ptr<TH1> h(dynamic_cast<TH1 *>(original->Clone()));
        if (!h) { ++skipped; continue; }
        h->SetDirectory(nullptr);
        const bool empty = GemPdf::Empty(*h);
        if (empty) ++emptyCount;
        if (empty && skipEmpty) { ++skipped; GemPdf::Text(0.5,0.5,"Empty (skipped)",0.045,22); continue; }
        pad->SetLogx(false); pad->SetLogy(false); pad->SetLogz(false);
        pad->SetLeftMargin(0.14); pad->SetBottomMargin(0.19);
        pad->SetTopMargin(0.14); pad->SetRightMargin(h->GetDimension() == 1 ? 0.05 : 0.21);
        pad->SetTicks(1, 1);
        h->SetStats(GemPdf::SHOW_STATISTICS); h->SetLineColor(GemPdf::LINE_COLOR); h->SetLineWidth(2);
        h->SetFillStyle(0); h->SetMarkerColor(GemPdf::LINE_COLOR);
        h->SetMarkerStyle(20); h->SetMarkerSize(0.65);
        GemPdf::Axis(h->GetXaxis()); GemPdf::Axis(h->GetYaxis()); GemPdf::Axis(h->GetZaxis());
        h->GetXaxis()->SetTitleOffset(1.35);
        h->GetYaxis()->SetTitleOffset(1.45);
        h->GetZaxis()->SetTitleOffset(1.35);
        // Preserve axis titles and units supplied by the analysis.
        if (!std::string(h->GetXaxis()->GetTitle()).size()) h->GetXaxis()->SetTitle("Bin coordinate");
        if (h->GetDimension() == 1 && !std::string(h->GetYaxis()->GetTitle()).size())
            h->GetYaxis()->SetTitle(dynamic_cast<TProfile *>(h.get()) ? "Mean value" : "Bin content");
        bool positive = false, negative = false;
        for (int bin=0; bin<h->GetNcells(); ++bin) {
            positive |= h->GetBinContent(bin) > 0;
            negative |= h->GetBinContent(bin) < 0;
        }
        const auto limits = GemPdf::Settings(*h);
        GemPdf::CoordinateRange(h->GetXaxis(), limits.xmin, limits.xmax);
        if (h->GetDimension() > 1)
            GemPdf::CoordinateRange(h->GetYaxis(), limits.ymin, limits.ymax);
        const bool requestedLog = limits.log < 0 ? (h->GetDimension()==1 ? logY : logZ) : limits.log;
        const bool allowLog = positive && !negative;
        const bool useLog = requestedLog && allowLog;
        if (requestedLog && !allowLog)
            std::cerr << "Using linear scale for " << path << " (empty or negative bins).\n";
        if (h->GetDimension() == 1) {
            pad->SetLogy(useLog);
            GemPdf::ContentRange(*h, limits.ymin, limits.ymax, useLog);
            h->Draw(dynamic_cast<TProfile *>(h.get()) ? "E1" : "HIST");
        } else if (h->GetDimension() == 2) {
            pad->SetLogz(useLog);
            GemPdf::ContentRange(*h, limits.zmin, limits.zmax, useLog);
            h->SetContour(100); h->Draw("COLZ");
        } else {
            // 3D histograms remain 3D; no undocumented projection is made.
            h->Draw("BOX2Z");
        }
        const std::string title = std::string(h->GetTitle()).empty() ? path : h->GetTitle();
        GemPdf::Text(0.14, 0.94, (page.detectors ? "Detector " + std::to_string(panel) + " | " : "") + title, 0.040);
        if (empty) GemPdf::Text(0.5, 0.52, "Empty histogram", 0.040, 22);
        anyDrawn=true;
        drawn.push_back(std::move(h));
        }
        if (anyDrawn) {
            canvas.Modified(); canvas.Update(); canvas.Print(pdf.Data());
            std::cout << "Page " << ++pages << ": " << page.label << '\n';
        }
        canvas.Clear(); // Clear pads before destroying histogram clones.
    }
    canvas.Print((pdf + "]").Data());
    previous->cd();
    std::cout << "Wrote " << pages << " pages to " << output
              << "; empty histograms: " << emptyCount << "; skipped: " << skipped << '\n';
}
