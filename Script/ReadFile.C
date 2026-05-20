#define MAX_LINE 2048
void ReadFile(){
  vector<double> lambda;
  vector<double> lambdaEr;
  vector<double> TranspN2[10];
  vector<double> centralValueN2;
  char fname[64];
  int sizeAr = 1;
  for(int i=0; i<10; i++){
    //double cv = 0;
    sprintf(fname,"../nitro_8bar_baseline_05_13/C2F6_3bar.Sample.Cycle%d.ASC",(i+1));
    cout<<fname<<endl;
    FILE *infile;
    infile =  fopen(fname, "r");
    int foundData =0;
    char line[MAX_LINE];
    while(fgets(line,sizeof(line),infile)){
      line[strcspn(line, "\r\n")] = 0;

      if (strcmp(line, "#DATA") == 0) {
	foundData = 1;
	break;
      }
    }
    if(!foundData)
      cout<<"No Data"<<endl;
    // Read two-column data
    double x, y;

    while (fscanf(infile, "%lf %lf", &x, &y) == 2) {
      // Process data here
      //printf("%f   %f\n", x, y);
      if(i==0){
	lambda.push_back(x);
	lambdaEr.push_back(1.0);
      }
      TranspN2[i].push_back(y);
    }
    fclose(infile);
  }
  vector<double> mean;
  vector<double> sigma;
  sizeAr = lambda.size();
  for(int l =0; l<sizeAr; l++){
    double cv =0;
    for(int i=0; i<10; i++){
      cv += TranspN2[l][i];
    }
    mean.push_back(cv/10);
    double sig =0;
    for(int i=0; i<10; i++){
      sig += (TranspN2[l][i] - (cv/i));
    }
    sigma.push_back(sqrt(sig/10));
  }
  TGraphErrors * g = new TGraphErrors(sizeAr, &lambda[0], &mean[0], &lambdaEr[0], &sigma[0]);
  g->Draw("AP");
  g->SetMarkerStyle(20);
  g->SetMarkerColor(20);
}
