// Standard MICE boiler plate for all plots - update the MICE version & dates as appropriate  
TText t1 = TText(0.12,0.785,"MICE Preliminary");
TText t3 = TText(0.12,0.75,"ISIS cycle 2015/04");
TText t2 = TText(0.12,0.715,"LiH, MAUS v3.1.2");
t1.SetNDC(1);
t1.SetTextSize(0.04);
t1.SetTextFont(42);
t2.SetNDC(1);
t2.SetTextSize(0.03);
t2.SetTextFont(42);
t3.SetNDC(1);
t3.SetTextSize(0.04);
t3.SetTextFont(42);
t1.Draw();
t3.Draw();
t2.Draw();

// Add legend to plot - should not be necessary for your work
TLegend *leg = new TLegend(.65,.65,.85,.85);
leg->AddEntry(plot_name, "plot_name","l");
leg->Draw();

// Set name, title & labels for ROOT plot
plot_name->SetName("plot_name");
plot_name->SetTitle("plot_title");
plot_name->GetXaxis()->SetTitle("axis_titlen");
