const int n_pl = 4;
bool logx = false;
TString pls[n_pl] = {"FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN","FTFP_BERT_EMM","QGSP_BERT"};
TString dRs[9] = {"dR < 0.56cm","dR < 1.0cm","dR < 2.0cm","dR < 3.0cm","dR < 5.0cm","dR < 8.0cm","dR < 12.0cm","dR < 18.0cm","dR < 20.0cm"};
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[2] = {"V18","V20"};//{"WJets","TTJets_HT","TTJets","pMSSM_MCMC_473_54451"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"}; 
int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
int line_color[9] = {kRed, kBlue, kRed, kGreen+2, kMagenta, kRed - 3, kBlue + 2 , kCyan + 1 , kGreen + 3 };
//int line_color[11] = {kPink+1, kRed, kBlue,kGray+1 , kGreen+2, kMagenta, kYellow + 2 , kCyan+3,  kBlue + 2 ,kRed+2,kGreen + 3 };
void setLastBinAsOverFlow(TH1F*);
std::map<int,int> rebin_;
std::map<int,float> x_max;
std::map<int,float> x_min;
void rebin_init_layer_pi(){
  rebin_.insert(pair<int, int>(20,2));
  rebin_.insert(pair<int, int>(50,2));
  rebin_.insert(pair<int, int>(100,2));
  rebin_.insert(pair<int, int>(200,4));
    
  x_max.insert(pair<int, int>(20,2500.0));
  x_max.insert(pair<int, int>(50,6000.0));
  x_max.insert(pair<int, int>(100,12000.0));
  x_max.insert(pair<int, int>(200,22000.0));
  
  x_min.insert(pair<int, int>(20,0.0));
  x_min.insert(pair<int, int>(50,0.0));
  x_min.insert(pair<int, int>(100,0.0));
  x_min.insert(pair<int, int>(200,0.0));
  x_min.insert(pair<int, int>(300,0.0));
  
}
void setLastBinAsOverFlow(TH1F* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);

  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);

}

void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", int energy=-1, char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title=""){  
  

  TCanvas *c = new TCanvas(tag_name, tag_name, 700, 600);
  c->SetRightMargin(0.009);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.05);
 c->SetRightMargin(0.009);
  c->SetLeftMargin(0.11);
  c->SetTopMargin(0.05);
c->SetLeftMargin(0.11);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.025);
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.06);
   c->SetBottomMargin(0.12);

  //    gStyle->SetOptStat(1111111);
     gStyle->SetOptStat(0);
  //double pvt_x_min = 0.6;
  double pvt_x_min = 0.75;
  double pvt_x_max = 0.99;
  double pvt_y_min = 0.9;
  //double pvt_dely = 0.18;
  double pvt_dely = 0.15;
  //gStyle->SetOptStat(0);
char* integral_tag = new char[100];

  //gStyle->SetOptFit(0);
  vector<TString> legName;
  //TLegend *legend = new TLegend(0.65,0.95,0.99,0.75);
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.7,0.9,0.96,0.75);  
  legend->SetTextSize(0.035);
  char* lhead = new char[100];
  sprintf(lhead,"#bf{%s:%s} ",leg_head,title);
  legend->SetHeader(lhead);

 auto legend1 = new TLegend(0.6,0.6,0.9,0.73);
  legend1->SetTextSize(0.038);
  legend1->SetLineColor(kWhite);
   legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry1[11];

  TLegendEntry* leg_entry[11];
  float x_label_size = 0.025;
  double ymin = 100000.0;
  double ymax = 0.0;
  cout<<" hist.size() = "<<hist.size()<<endl;


  for(int i = 0; i < (int)hist.size(); i++) {
     hist.at(i)->SetLineWidth(line_width[i]);

    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);

    sprintf(integral_tag,"Integral: %0.2f",hist.at(i)->Integral());
    leg_entry1[i]=legend1->AddEntry(hist.at(i),integral_tag,"");
    leg_entry1[i]->SetTextColor(hist.at(i)->GetLineColor());

    
    if(normalize) {
      hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
          }
    
    hist.at(i)->SetLineWidth(line_width[i]);
    
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" ");
    
    hist.at(i)->GetXaxis()->SetTitleSize(0.04);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.030);
    //    hist.at(i)->GetXaxis()->SetRange(0,1500);
    hist.at(i)->GetYaxis()->SetLabelSize(0.030);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.45);
    // if(DoRebin) {
    //  hist.at(i)->Rebin(4);
    //   //hist.at(i)->Rebin(1);
    // }
    /* hist.at(i)->GetXaxis()->SetRangeUser(x_min[energy],x_max[energy]); */
    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_text[i],"l");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    //  sprintf(integral_tag,"Integral: %0.2f",hist.at(i)->Integral());
    // leg_entry1[i]=legend1->AddEntry(hist.at(i),integral_tag,"");
    // leg_entry1[i]->SetTextColor(hist.at(i)->GetLineColor());

    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    setLastBinAsOverFlow(hist.at(i));
    
    

  }
  if(ymin == 0.0) ymin = 1e-2;
  if(ymin<0.0) ymin = 1e-1;
  //  if(ymax<=10) ymax=10;
  for(int i = 0; i < (int)hist.size(); i++) {
    if(!log_flag) hist.at(i)->GetYaxis()->SetRangeUser(0.0,ymax*1.9);
    else
      {  hist.at(i)->GetYaxis()->SetRangeUser(0.0000000001,ymax*10);
	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
	
      }
  	if(!i) hist.at(i)->Draw("hist");
    else   hist.at(i)->Draw("hist sames");

  }
  
   legend->Draw();
     legend1->Draw();



  if(log_flag) {
      gPad->SetLogy();
    }
  if(logx)
    gPad->SetLogx();

  
  gPad->Update();

 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.035);
  //  textOnTop->DrawLatexNDC(0.12,0.96,"CMS #it{#bf{Preliminary}}");
  
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.035);
  int inlumi=137;
  // sprintf(en_lat,"#bf{%0.1d fb^{-1} (13 TeV)}",inlumi);
  //textOnTop->DrawLatexNDC(0.7,0.96,en_lat);

 

  gPad->Modified();
  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"./Results/%s.png",tag_name);
     c->SaveAs(canvas_name);
    
    sprintf(canvas_name,"./Results/%s.pdf",tag_name);
    c->SaveAs(canvas_name);
    
  }
  
}
const int nfiles=200,nBG=6;                                                                                                                                                              
TFile *f[nfiles];
TFile *f1[nfiles];


void overlayPlots(string pathname, int which_year)
{
  char* hname = new char[200];
  char* hist_name  = new char[200];
  char* hist_name1 = new char[200];
  char* hist_name2 = new char[200];
  char* hist_name3 = new char[200];
  char* hist_name4 = new char[200];
  char* hist_name5 = new char[200];
  char* hist_name6 = new char[200];

  char* full_path = new char[2000];
  char* full_path1 = new char[2000];
  char* full_path2 = new char[2000];
  char* path2 = new char[2000];
  char* title= new char[2000];
  //string filetag;//=new char[20000];                                                                                                                                                                   
  char* full_path3 = new char[2000];
  char* full_path4 = new char[2000];
  char* full_path5 = new char[2000];
  char* full_path6 = new char[2000];
  char* full_path7 = new char[2000];
  char* full_path8 = new char[2000];
  char* full_path9 = new char[2000];
  char* full_path10 = new char[2000];
  char* full_path11= new char[2000];
  char *leg_head = new char[200];
  sprintf(path2,"Results/dPhi_METHadJets");
// TFile *f[nfiles];
//  f[0]=                                                         
    //  /* f[0] = new TFile("Out_TTGJets_v18_Dphi_MetHadJets_corr_v2.r
  char* dataset = new char[100];
  f[0] = new TFile("out.root");
  // f1[0] = new TFile("out.root");
  
  /*if(which_year ==0)
    {
      sprintf (dataset,"Autumn18");
      f[0] = new TFile("Out_Autumn18_WJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
      f[1] = new TFile("Out_Autumn18_TTJets-SiLept_preUL_v18_Electron_wPho_pT_g40_MET200.root");
      f[2] = new TFile("Out_Autumn18_TTJets-DiLept_preUL_v18_Electron_wPho_pT_g40_MET200.root");
      f[3] = new TFile("Out_Autumn18_TTJets-HT_preUL_v18_Electron_wPho_pT_g40_MET200.root");
      f[4] = new TFile("Out_Autumn18_TTGJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");   
      f[5] = new TFile("Out_Autumn18_ZJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
      f[6] = new TFile("Out_Autumn18_GJets_DR_preUL_v18_Electron_wPho_pT_g40_MET200.root");
      f[7] = new TFile("Out_Autumn18_QCD_preUL_v18_Electron_wPho_pT_g40_MET200.root");
      f1[0] = new TFile("Summer20UL18_WJetsToLNu_HT_pt40_Electron_v20.root");
      f1[1] = new TFile("Summer20UL18_TTJets_SingleLeptFromT_pt40_Electron_v20.root");
      f1[2] = new TFile("Summer20UL18_TTJets_DiLept_pt40_Electron_v20.root");
      f1[3] = new TFile("Summer20UL18_TTJets_HT_pt40_Electron_v20.root");
      f1[4] = new TFile("Summer20UL18_TTGJets_pt40_Electron_v20.root");
      f1[5] = new TFile("Summer20UL18_ZJetsToNuNu_HT_pt40_Electron_v20.root");
      f1[6] = new TFile("Summer20UL18_GJets_DR-0p4_HT_pt40_Electron_v20.root");
      f1[7] = new TFile("Summer20UL18_QCD_Jets-HT_pt40_Electron_v20.root");

    }
  else if (which_year==1){
    sprintf (dataset,"Fall17");

    f[0] = new TFile("Out_Fall17_WJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[1] = new TFile("Out_Fall17_TTJets-SiLept_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[2] = new TFile("Out_Fall17_TTJets-DiLept_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[3] = new TFile("Out_Fall17_TTJets-HT_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[4] = new TFile("Out_Fall17_TTGJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[5] = new TFile("Out_Fall17_ZJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[6] = new TFile("Out_Fall17_GJets_DR_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[7] = new TFile("Out_Fall17_QCD_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f1[0] = new TFile("Summer20UL17_WJetsToLNu_HT_pt40_Electron_v20.root");
    f1[1] = new TFile("Summer20UL17_TTJets_SingleLeptFromT_pt40_Electron_v20.root");
    f1[2] = new TFile("Summer20UL17_TTJets_DiLept_pt40_Electron_v20.root");
    f1[3] = new TFile("Summer20UL17_TTJets_HT_pt40_Electron_v20.root");
    f1[4] = new TFile("Summer20UL17_TTGJets_pt40_Electron_v20.root");
    f1[5] = new TFile("Summer20UL17_ZJetsToNuNu_HT_pt40_Electron_v20.root");
    f1[6] = new TFile("Summer20UL17_GJets_DR-0p4_HT_pt40_Electron_v20.root");
    f1[7] = new TFile("Summer20UL17_QCD_Jets-HT_pt40_Electron_v20.root");


  }
  else if(which_year==2){
    sprintf (dataset,"Summer16v3");
    f[0] = new TFile("Out_Summer16v3_WJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[1] = new TFile("Out_Summer16v3_TTJets-SiLept_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[2] = new TFile("Out_Summer16v3_TTJets-DiLept_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[3] = new TFile("Out_Summer16v3_TTJets-HT_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[4] = new TFile("Out_Summer16v3_TTGJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[5] = new TFile("Out_Summer16v3_ZJets_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[6] = new TFile("Out_Summer16v3_GJets_DR_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f[7] = new TFile("Out_Summer16v3_QCD_preUL_v18_Electron_wPho_pT_g40_MET200.root");
    f1[0] = new TFile("Summer20UL16_WJetsToLNu_HT_pt40_Electron_v20.root");
    f1[1] = new TFile("Summer20UL16_TTJets_SingleLeptFromT_pt40_Electron_v20.root");
    f1[2] = new TFile("Summer20UL16_TTJets_DiLept_pt40_Electron_v20.root");
    f1[3] = new TFile("Summer20UL16_TTJets_HT_pt40_Electron_v20.root");
    f1[4] = new TFile("Summer20UL16_TTGJets_pt40_Electron_v20.root");
    f1[5] = new TFile("Summer20UL16_ZJetsToNuNu_HT_pt40_Electron_v20.root");
    f1[6] = new TFile("Summer20UL16_GJets_DR-0p4_HT_pt40_Electron_v20.root");
    f1[7] = new TFile("Summer20UL16_QCD_Jets-HT_pt40_Electron_v20.root");

    }*/

  //  const char* baseline[11] = {"nocut","bkg_comp","veto-lep","veto-iso","Phot-pt_20","NhadJetsCut","dPhi_Met","Met_100","Met_250","st_300_Met250","pt_st_Met_250"};
  vector<string>baseline1;
  //if(which_Lept)
  //baseline1={"h_NhadJets_PreSelection","h_NBJets_PreSelection","h_MET_PreSelection","h_PhoPt_PreSelection","h_Mt_phoMET_PreSelection","h_dPhi_phoMet_PreSelection","h_St_PreSelection"};
  baseline1={"h_MET,h_MET2"};//,"h_elec_pT","h_W_pT","h_Nu_pT","h_mu_pT","h_tau_pT","h_elec_Eta","h_W_Eta","h_Nu_Eta","h_mu_Eta","h_tau_Eta","h_elec_Phi","h_W_Phi","h_Nu_Phi","h_mu_Phi","h_tau_Phi","h_elec_E","h_W_E","h_Nu_E","h_mu_E","h_tau_E","h_GenMET","h_GenHT","h_GenJets","h_recoElec_pT","h_recoMu_pT","h_recoMu_Eta","h_recoElec_Eta","h_recoMu_Phi","h_recoElec_Phi","h_leadJet_pT","h_subleadJet_pT","h_leadJet_Eta","h_subleadJet_Eta"};//
  vector<string> xtitle;
  xtitle = {"N_{Jets}","N_{b-jets}","p^{miss}_{T}","p^{#gamma}_{T}","M^{MET & #gamma}_{T}","d#phi^{MET,#gamma}","Sum of p^{Jets}_{T} & p^{#gamma}_{T}","p^{Elec}_{T}","p^{W}_{T}","p^{Nu}_{T}","p^{Mu}_{T}","p^{Tau}_{T}","Eta coordinate for Elec","Eta coordinate for W","Eta coordinate for Nu","Eta coordinate for Mu","Eta coordinate for Tau","Phi coordinate for Elec","Phi coordinate for W","Phi coordinate for Nu","Phi coordinate for Mu","Phi coordinate for Tau","E_{Elec}","E_{W}","E_{Nu}","E_{Mu}","E_{Tau}","Gen MET","Gen HT","GenJets","reco p^{Elec}_{T}","reco p^{Mu}_{T}","Reco-elec Eta","Reco-mu Eta","Reco-elec Phi","Reco-mu Phi","lead Jet p_{T}","sub-lead Jet p_{T}","lead Jet Eta","sub-lead Jet Eta"};
  cout<<"xtitle.size() "<<"\t"<<xtitle.size()<<endl;
  cout<<"baseline1.size()  "<<baseline1.size()<<endl;  
  vector<int>rebin ={1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2};
  cout<<"rebin.size()  "<<rebin.size()<<endl;
  vector<int> xrange ={16,16,1500,1000,1000,1,1000,1000,1000,1000,1000,1000,1,1,1,1,1,1,1,1,1,1,1000,1000,1000,1000,1000,1500,1500,50,1000,1000,1,1,1,1,1000,1000,1,1};
  cout<<"xrange.size() "<<xrange.size()<<endl;
    vector<string>baseline2 = {"h_madminPhotonDeltaR_after","h_minPho_lep_after"};
  // else
  //   baseline1 = {"Mu_SR","Mu_CR","Mu_failpTcut_SR","Mu_failEtacut_SR","Mu_failAccep_tauHadronic_SR"};

  // const char* baseline[9]= {"Nocut","bkg_comp","veto-lep","veto-iso","dPhi_Met","Met_250","st_300_Met250","pt_st_Met_250","nocut"};
  // const char* baseline1[11]= {"Nocut","nocut_sam","nocut","Met_100","bkg_comp","veto-lep","veto-iso","dPhi_Met","Met_250","st_300_Met250","pt_st_Met_250"};
    const char* filetag[10]={"WJetsToLNu-HT","TTJets_SingleElec","TTJets_DiElec","TTJets-HT","WJets","ZJets","GJets","QCD-Jets","T5bbbbZG_10","T5bbbbZG_50"};//,"T5bbbbZG_200","T5bbbbZG_1500"};
    //    const char* dataset[3] ={"Autumn18","Fall17","Summer16v3"};
   /* vector<TH1F*> hist_list_Njets; */
   /*    vector<TH1F*> hist_list_Bjets; */
   /*    vector<TH1F*> hist_list_MET; */
   /*    vector<TH1F*> hist_list_PhoPt; */
   /*    //vector<TH1F*> hist_list_Mt;                                                                                                                                    */
   /*    vector<TH1F*> hist_list_ST; */
   /*    vector<TH1F*> hist_list_HT; */

  // for(int i_file=0; i_file<2;i_file++)
  //   {
      // vector<TH1F*> hist_list_Njets;
      // vector<TH1F*> hist_list_Bjets;
      // vector<TH1F*> hist_list_MET;
      // vector<TH1F*> hist_list_PhoPt;
      // //vector<TH1F*> hist_list_Mt;
      // vector<TH1F*> hist_list_ST;
      // vector<TH1F*> hist_list_HT;
    //for(int i_cut=0; i_cut<7;i_cut++)
    //{
	  
      for(int i_file=0; i_file<2;i_file++)
	{
	   vector<TH1F*> hist_list_Njets;
	   //sprintf(hist_name,"%s",baseline1[i_cut].c_str());
	   sprintf(hist_name,"%s",baseline1[i_file].c_str());
	  TH1F* hist_Njets_temp = (TH1F*)f[i_file]->Get(hist_name);
	  //TH1F* hist_Bjets_temp = (TH1F*)f1[i_file]->Get(hist_name);
	  //hist_Njets_temp->Rebin(rebin[i_cut]);
	  //hist_Njets_temp->Rebin(rebin[i_cut]);
	  //hist_Bjets_temp->Rebin(rebin[i_cut]);
	  //	  sprintf(title,"%s",
	  //hist_phopt_temp->Rebin(2);
	  // hist_Ht_temp->Rebin(4);
	  // hist_ST_temp->Rebin(4);
	  hist_Njets_temp->GetXaxis()->SetTitle(xtitle[i_cut].c_str());//"mindR(gen-#gamma, q/g)");
	  //sprintf(title,"%s",baseline1[i_cut].c_str());
	  sprintf(title,"%s",baseline1[i_file].c_str());	  
	  //hist_Bjets_temp->GetXaxis()->SetTitle(xtitle[i_cut].c_str());//"mindR(reco-#gamma,gen-l)");
	  if(xrange[i_cut]!=1)
	    { hist_Bjets_temp->GetXaxis()->SetRange(0,xrange[i_cut]); hist_Njets_temp->GetXaxis()->SetRange(0,xrange[i_cut]);}

	  hist_list_Njets.push_back(hist_Njets_temp);
	  hist_list_Njets.push_back(hist_Bjets_temp);
	
      int energy=-1;
      sprintf(full_path,"%s_Overlay_%s_%s",dataset,baseline1[i_cut].c_str(),filetag[i_file]);
      generate_1Dplot(hist_list_Njets,full_path,energy,dataset,true,true,false,true,filetag[i_file]); 
      /* sprintf(full_path6,"Overlay_dPhi_photon_MET_%s",filetag[i_file]); */
      /* generate_1Dplot(hist_list_dPhi,full_path6,energy,leg_head,false,true,false,true); */
      // cout<<hist_list_Njets.at(0)->Integral()<<endl;
      // cout<<hist_list_Njets.at(1)->Integral()<<endl;
      // cout<<"after"<<endl;
      // cout<<hist_list_Bjets.at(0)->Integral()<<endl;
      // cout<<hist_list_Bjets.at(1)->Integral()<<endl;

      // // sprintf(full_path,"Overlay_madmindRphoton_WG_%s_norm",title);//,filetag[i_file]);
      // // generate_1Dplot(hist_list_Njets,full_path,energy,leg_head,true,true,false,true);
      // // sprintf(full_path1,"Overlay_minDrpho_genlep_WG_%s_norm",title);//,filetag[i_file]);
      // // generate_1Dplot(hist_list_Bjets,full_path1,energy,leg_head,true,true,false,true);
      
      //  sprintf(full_path,"Overlay_madmindRphoton_WG_%s_wlog",title);//,filetag[i_file]);                                                       
      //  generate_1Dplot(hist_list_Njets,full_path,energy,leg_head,false,true,false,true);
      //  sprintf(full_path1,"Overlay_minDrpho_genlep_WG_%s_wlog",title);//,filetag[i_file]);                                                        
      //generate_1Dplot(hist_list_Bjets,full_path1,energy,leg_head,false,true,false,true);
	}

      // sprintf(full_path2,"Overlay_MET_%s_wlog",filetag[i_file]);
      // generate_1Dplot(hist_list_MET,full_path2,energy,leg_head,false,true,false,true);
      // sprintf(full_path3,"Overlay_PhoPt_%s_wlog",filetag[i_file]);
      // generate_1Dplot(hist_list_PhoPt,full_path3,energy,leg_head,false,true,false,true);
      // sprintf(full_path4,"Overlay_Ht_%s_wlog",filetag[i_file]);
      // generate_1Dplot(hist_list_HT,full_path4,energy,leg_head,false,true,false,true);
      // sprintf(full_path5,"Overlay_ST_%s_wlog",filetag[i_file]);
      // generate_1Dplot(hist_list_ST,full_path5,energy,leg_head,false,true,false,true);
	

      //       }
      
      //}

}





