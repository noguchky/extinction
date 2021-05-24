
Int_t haddCoin(const std::string& ifilename) {
  std::string ofileprefix;
  {
    TString buff = ifilename.data();
    if (buff.EndsWith(".root")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofileprefix = buff;
  }

  const std::string ofilename = ofileprefix + "_coin.pdf";

  TFile* ifile = new TFile(ifilename.data());
  if (!ifile->IsOpen()) {
    std::cout << "[error] file is not opened, " << ifilename << std::endl;
    return 1;
  }
  
  TH1D* hBh1Timeline   = dynamic_cast<TH1D*>(ifile->Get("hBh1Timeline"));
  TH1D* hBh2Timeline   = dynamic_cast<TH1D*>(ifile->Get("hBh2Timeline"));
  TH1D* hTc1Timeline   = dynamic_cast<TH1D*>(ifile->Get("hTc1Timeline"));
  TH1D* hTc2Timeline   = dynamic_cast<TH1D*>(ifile->Get("hTc2Timeline"));
  TH1D* hExtTimeline   = dynamic_cast<TH1D*>(ifile->Get("hExtTimeline"));
  TH1D* hCoinTdcInSync = dynamic_cast<TH1D*>(ifile->Get("hCoinTdcInSync"));
  TH2D* hCoinMountain  = dynamic_cast<TH2D*>(ifile->Get("hCoinMountain"));

  std::cout << "Draw plots" << std::endl;
  if (!gPad) {
    TCanvas::MakeDefCanvas();
  }
  gPad->SetGrid(true, true);
  gPad->SetLogy(false);
  gPad->Print((ofilename + "[").data());

  std::cout << "hCoinTdcInSync" << std::endl;
  gPad->SetLogy(true);
  {
    hBh1Timeline->SetLineColor(51);
    hBh2Timeline->SetLineColor(61);
    hTc1Timeline->SetLineColor(71);
    hTc2Timeline->SetLineColor(81);
    hExtTimeline->SetLineColor(91);

    hBh1Timeline->Draw();
    hBh2Timeline->Draw("same");
    hTc1Timeline->Draw("same");
    hTc2Timeline->Draw("same");
    hExtTimeline->Draw("same");

    gPad->Print(ofilename.data());
  }
  gPad->SetLogy(false);

  std::cout << "hCoinTdcInSync" << std::endl;
  gPad->SetLogy(true);
  {
    hCoinTdcInSync->Draw("hist");
    hCoinTdcInSync->SetMinimum(0.2);
    gPad->Print(ofilename.data());

    // hCoinTdcInSync->Draw();
    // hCoinTdcInSync->SetMinimum(0.2);
    // gPad->Print(ofilename.data());
  }
  gPad->SetLogy(false);

  std::cout << "hCoinMountain" << std::endl;
  gPad->SetGrid(false, true);
  {
    hCoinMountain->Draw("colz");
    hCoinMountain->SetMinimum(0);
    gPad->Print(ofilename.data());
  }
  gPad->SetGrid(true, true);

  gPad->Print((ofilename + "]").data());

  return 0;
}
