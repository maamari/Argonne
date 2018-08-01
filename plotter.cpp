{
  TCanvas *c = new TCanvas();
  c->SetFillColor(0);
  c->SetGrid();
  
  TGraph *g = new TGraph("output.txt");
  g->SetTitle("18F(d,p)");
  g->GetXaxis()->SetTitle("Angle (Degrees)");
  g->GetYaxis()->SetTitle("KE (MeV)");
  g->SetLineWidth(2);
  g->SetLineColor(2);
  g->Draw("ACP");
  g->GetXaxis()->SetRangeUser(0,180);
  g->Draw("ACP");
}
