#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "TGraphErrors.h"

#define MAX_LINE 2048

void Analysis() {

    const int Ncycles = 10;

    std::vector<double> lambda;
    std::vector<double> lambdaEr;
    std::vector< std::vector<double> > TranspN2(Ncycles);
    std::vector< std::vector<double> > TranspC2F6(Ncycles);

    std::vector<double> meanN2;//(Npoints, 0.0);
    std::vector<double> sigmaN2;//(Npoints, 0.0);
    std::vector<double> meanC2F6;//(Npoints, 0.0);
    std::vector<double> sigmaC2F6;//(Npoints, 0.0);
    char fname[256];
    TGraph * gr[10];
    for (int i = 0; i < Ncycles; i++) {
      gr[i] = new TGraph();
      sprintf(fname,
	      "../nitro_8bar_baseline_05_13/N2_8bar.Sample.Cycle%d.ASC",
	      i + 1);

      std::cout << "Reading: " << fname << std::endl;
      
      FILE *infile = fopen(fname, "r");
      if (!infile) {
	std::cerr << "Cannot open file: " << fname << std::endl;
	continue;
      }

      char line[MAX_LINE];
      bool foundData = false;

      while (fgets(line, sizeof(line), infile)) {
	line[strcspn(line, "\r\n")] = 0;
	if (strcmp(line, "#DATA") == 0) {
	  foundData = true;
	  break;
	}
      }

      if (!foundData) {
	std::cerr << "#DATA not found in " << fname << std::endl;
	fclose(infile);
	continue;
      }

      double x, y;
      int idx = 0;

      while (fscanf(infile, "%lf %lf", &x, &y) == 2) {

	// store wavelength only once (assumes identical grid)
	if (i == 0) {
	  lambda.push_back(x);
	  lambdaEr.push_back(1.0);
	}
	gr[i]->SetPoint(idx,x,y);
	TranspN2[i].push_back(y);
	idx++;
      }
      gr[i]->SetLineColor(30+i);
      if(i==0){
	gr[i]->Draw("AL");
      }
      else
	gr[i]->Draw("L SAMES");
      fclose(infile);
    }

    const int Npoints = lambda.size();
    meanN2.resize(Npoints);
    sigmaN2.resize(Npoints);
    meanC2F6.resize(Npoints);
    sigmaC2F6.resize(Npoints);

    if (Npoints == 0) {
        std::cerr << "No data loaded." << std::endl;
        return;
    }
    // compute mean
    for (int l = 0; l < Npoints; l++) {

        double sum = 0.0;
        int validCycles = 0;

        for (int i = 0; i < Ncycles; i++) {
            if (l < (int)TranspN2[i].size()) {
                sum += TranspN2[i][l];
                validCycles++;
            }
        }

        if (validCycles == 0) continue;

        meanN2[l] = sum / validCycles;

        double var = 0.0;

        for (int i = 0; i < Ncycles; i++) {
            if (l < (int)TranspN2[i].size()) {
                double d = TranspN2[i][l] - meanN2[l];
                var += d * d;
            }
        }

        sigmaN2[l] = std::sqrt(var / validCycles);
    }
    TCanvas *cN2 = new TCanvas();
    // ROOT graph
    TGraphErrors *g = new TGraphErrors(
        Npoints,
        &lambda[0],
        &meanN2[0],
        &lambdaEr[0],
        &sigmaN2[0]
    );
    
    g->SetMarkerStyle(24);
    g->SetMarkerColor(kBlue);
    g->SetTitle("Transmission;Wavelength;Transmission");

    g->Draw("AP");

    ////////////////////////////////////////////
    // Repeat for C2F6
    ////////////////////////////////////////////
    TCanvas *c = new TCanvas();
    TGraph* gr2[10];
    for (int i = 0; i < Ncycles; i++) {
      gr2[i]= new TGraph();
      sprintf(fname,
	      "../nitro_8bar_baseline_05_13/C2F6_3bar.Sample.Cycle%d.ASC",
	      i + 1);

      std::cout << "Reading: " << fname << std::endl;

      FILE *infile = fopen(fname, "r");
      if (!infile) {
	std::cerr << "Cannot open file: " << fname << std::endl;
	continue;
      }

      char line[MAX_LINE];
      bool foundData = false;

      while (fgets(line, sizeof(line), infile)) {
	line[strcspn(line, "\r\n")] = 0;
	if (strcmp(line, "#DATA") == 0) {
	  foundData = true;
	  break;
	}
      }

      if (!foundData) {
	std::cerr << "#DATA not found in " << fname << std::endl;
	fclose(infile);
	continue;
      }

      double x, y;
      int idx = 0;

      while (fscanf(infile, "%lf %lf", &x, &y) == 2) {

	// store wavelength only once (assumes identical grid)
	if (i == 0) {
	  lambda.push_back(x);
	  lambdaEr.push_back(1.0);
	}
	gr2[i]->SetPoint(idx,x,y);
	TranspC2F6[i].push_back(y);
	idx++;
      }
      gr2[i]->SetLineColor(40+i);
      if(i==0){
	gr2[i]->Draw("AL");
      }
      else gr2[i]->Draw("L SAMES");
      fclose(infile);
    }

    //const int Npoints = lambda.size();

    if (Npoints == 0) {
      std::cerr << "No data loaded." << std::endl;
      return;
    }

    // compute mean
    for (int l = 0; l < Npoints; l++) {

      double sum = 0.0;
      int validCycles = 0;

      for (int i = 0; i < Ncycles; i++) {
	if (l < (int)TranspC2F6[i].size()) {
	  sum += TranspC2F6[i][l];
	  validCycles++;
	}
      }

      if (validCycles == 0) continue;

      meanC2F6[l] = sum / validCycles;

      double var = 0.0;

      for (int i = 0; i < Ncycles; i++) {
	if (l < (int)TranspC2F6[i].size()) {
	  double d = TranspC2F6[i][l] - meanC2F6[l];
	  var += d * d;
	}
      }

      sigmaC2F6[l] = std::sqrt(var / validCycles);
    }
  
    TCanvas *cX = new TCanvas();
    // ROOT graph
    TGraphErrors *gC2F6 = new TGraphErrors(
					   Npoints,
					   &lambda[0],
					   &meanC2F6[0],
					   &lambdaEr[0],
					   &sigmaC2F6[0]
					   );

    gC2F6->SetMarkerStyle(24);
    gC2F6->SetMarkerColor(kMagenta);
    gC2F6->SetTitle("Transmission;Wavelength;Transmission");

    gC2F6->Draw("AP");



    /////////////////////////////////////
    ///// Draw Ratio
    ////////////////////////////////////
    vector<double> Ratio;
    vector<double> RatioEr;
    for(int i=0; i<Npoints;i++){
      double r = meanC2F6[i]/meanN2[i];
      double err  = (1/((7.2e+19)*12))*hypot((sigmaC2F6[i]/meanC2F6[i]),(sigmaN2[i]/meanN2[i]));

      double X_sec = -log(r)/((7.2e+19)*12);
      //double X_sec
      Ratio.push_back(X_sec);
      RatioEr.push_back(err);
    }
    TGraphErrors * gRatio =  new TGraphErrors(
	Npoints,
	&lambda[0],
	&Ratio[0],
	&lambdaEr[0],
	&RatioEr[0]
    );
    TCanvas *c3 = new TCanvas();
    gRatio->SetMarkerStyle(20);
    gRatio->SetMarkerColor(20);
    gRatio->SetTitle("Tr(N2/C2F6);Wavelength;Ratio");

    gRatio->Draw("AP");
    TF1 *fit1 = new TF1("fit1",
		       "[0]*TMath::Power(x,[1])",
		       lambda.front(),
		       lambda.back());

    fit1->SetParameters(1e-08, -5.0);

    gRatio->Fit(fit1, "FR","",190,250);

    //gRatio->GetListOfFunctions()->Add(fit);

    //gRatio->Draw("AP");

    //c3->Update();
}
