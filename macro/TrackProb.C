// TrackProb.C - making joint probability / vertex probability histograms

void TrackProb(const char *bbfile, const char *ccfile, const char *qqfile, const char *d0outfile, const char *z0outfile)
{
	TFile *fbb = TFile::Open(bbfile);
	TFile *fcc = TFile::Open(ccfile);
	TFile *fqq = TFile::Open(qqfile);

	TNtuple *ntbb = (TNtuple *)fbb->Get("tracks");
	TNtuple *ntcc = (TNtuple *)fcc->Get("tracks");
	TNtuple *ntqq = (TNtuple *)fqq->Get("tracks");

	cout << "Creating " << d0outfile << endl;
	TFile *fd0 = TFile::Open(d0outfile,"recreate");

	TH1F *hbd0 = new TH1F("hb","hb",38,-1.8,2);
	TH1F *hcd0 = new TH1F("hc","hc",38,-1.8,2);
	TH1F *hqd0 = new TH1F("hq","hq",38,-1.8,2);
	TH1F *ha0d0 = new TH1F("ha0","ha0",38,-1.8,2);

	TH1F *hbpd0 = new TH1F("hbp","hbp",38,-1.8,2);
	TH1F *hcpd0 = new TH1F("hcp","hcp",38,-1.8,2);
	TH1F *hqpd0 = new TH1F("hqp","hqp",38,-1.8,2);
	TH1F *ha0pd0 = new TH1F("ha0p","ha0p",38,-1.8,2);

	TH1F *hbnd0 = new TH1F("hbn","hbn",38,-1.8,2);
	TH1F *hcnd0 = new TH1F("hcn","hcn",38,-1.8,2);
	TH1F *hqnd0 = new TH1F("hqn","hqn",38,-1.8,2);
	TH1F *ha0nd0 = new TH1F("ha0n","ha0n",38,-1.8,2);

	TH1F *hbipd0 = new TH1F("hbip","hbip",200,-5,5);
	TH1F *hcipd0 = new TH1F("hcip","hcip",200,-5,5);
	TH1F *hqipd0 = new TH1F("hqip","hqip",200,-5,5);

	TH1F *hjprobd0 = new TH1F("hjprob","hjprob",50,0,5);
	TH1F *hjprob2d0 = new TH1F("hjprob2","hjprob2",195,5,200);

	// hb/hc/hq/ha0
	ntbb->Project("hb", "log10(abs(sd0))","abs(sd0sig)>5");
	ntcc->Project("hc", "log10(abs(sd0))","abs(sd0sig)>5");
	ntqq->Project("hq", "log10(abs(sd0))","abs(sd0sig)>5");

	ha0d0->Add(hbd0);
	ha0d0->Add(hcd0);
	ha0d0->Add(hqd0);

	// normalize
	hbd0->Divide(ha0d0);
	hcd0->Divide(ha0d0);
	hqd0->Divide(ha0d0);

	cout << "hb/hc/hq written." << endl;

	// hbp/hcp/hqp/ha0p
	ntbb->Project("hbp", "log10(sd0)","abs(sd0sig)>5&&sd0>0");
	ntcc->Project("hcp", "log10(sd0)","abs(sd0sig)>5&&sd0>0");
	ntqq->Project("hqp", "log10(sd0)","abs(sd0sig)>5&&sd0>0");

	ha0pd0->Add(hbpd0);
	ha0pd0->Add(hcpd0);
	ha0pd0->Add(hqpd0);

	// normalize
	hbpd0->Divide(ha0pd0);
	hcpd0->Divide(ha0pd0);
	hqpd0->Divide(ha0pd0);

	cout << "hbp/hcp/hqp written." << endl;

	// hbn/hcn/hqn/ha0n
	ntbb->Project("hbn", "log10(-sd0)","abs(sd0sig)>5&&sd0<0");
	ntcc->Project("hcn", "log10(-sd0)","abs(sd0sig)>5&&sd0<0");
	ntqq->Project("hqn", "log10(-sd0)","abs(sd0sig)>5&&sd0<0");

	ha0nd0->Add(hbnd0);
	ha0nd0->Add(hcnd0);
	ha0nd0->Add(hqnd0);

	// normalize
	hbnd0->Divide(ha0nd0);
	hcnd0->Divide(ha0nd0);
	hqnd0->Divide(ha0nd0);

	cout << "hbn/hcn/hqn written." << endl;

	// hbip/hcip/hqip not normalized
	ntbb->Project("hbip", "sd0sig","abs(sd0sig)<5");
	ntcc->Project("hcip", "sd0sig","abs(sd0sig)<5");
	ntqq->Project("hqip", "sd0sig","abs(sd0sig)<5");

	cout << "hbip/hcip/hqip written." << endl;

	ntqq->Project("hjprob", "abs(sd0sig)","jprobcut>0");
	ntqq->Project("hjprob2", "abs(sd0sig)","jprobcut>0");

	double integ = hjprobd0->Integral(0,50) + hjprob2d0->Integral(0,195);
	hjprobd0->Scale(1./integ);
	hjprob2d0->Scale(1./integ);

	cout << "hjprob/hjprob2 written." << endl;

	fd0->Write();
	fd0->Close();

	cout << "Creating " << z0outfile << endl;
	TFile *fz0 = TFile::Open(z0outfile,"recreate");

	TH1F *hbz0 = new TH1F("hb","hb",38,-1.8,2);
	TH1F *hcz0 = new TH1F("hc","hc",38,-1.8,2);
	TH1F *hqz0 = new TH1F("hq","hq",38,-1.8,2);
	TH1F *ha0z0 = new TH1F("ha0","ha0",38,-1.8,2);

	TH1F *hbpz0 = new TH1F("hbp","hbp",38,-1.8,2);
	TH1F *hcpz0 = new TH1F("hcp","hcp",38,-1.8,2);
	TH1F *hqpz0 = new TH1F("hqp","hqp",38,-1.8,2);
	TH1F *ha0pz0 = new TH1F("ha0p","ha0p",38,-1.8,2);

	TH1F *hbnz0 = new TH1F("hbn","hbn",38,-1.8,2);
	TH1F *hcnz0 = new TH1F("hcn","hcn",38,-1.8,2);
	TH1F *hqnz0 = new TH1F("hqn","hqn",38,-1.8,2);
	TH1F *ha0nz0 = new TH1F("ha0n","ha0n",38,-1.8,2);

	TH1F *hbipz0 = new TH1F("hbip","hbip",200,-5,5);
	TH1F *hcipz0 = new TH1F("hcip","hcip",200,-5,5);
	TH1F *hqipz0 = new TH1F("hqip","hqip",200,-5,5);

	TH1F *hjprobz0 = new TH1F("hjprob","hjprob",50,0,5);
	TH1F *hjprob2z0 = new TH1F("hjprob2","hjprob2",195,5,200);

	// hb/hc/hq/ha0
	ntbb->Project("hb", "log10(abs(sz0))","abs(sz0sig)>5");
	ntcc->Project("hc", "log10(abs(sz0))","abs(sz0sig)>5");
	ntqq->Project("hq", "log10(abs(sz0))","abs(sz0sig)>5");

	ha0z0->Add(hbz0);
	ha0z0->Add(hcz0);
	ha0z0->Add(hqz0);

	// normalize
	hbz0->Divide(ha0z0);
	hcz0->Divide(ha0z0);
	hqz0->Divide(ha0z0);

	cout << "hb/hc/hq written." << endl;

	// hbp/hcp/hqp/ha0p
	ntbb->Project("hbp", "log10(sz0)","abs(sz0sig)>5&&sz0>0");
	ntcc->Project("hcp", "log10(sz0)","abs(sz0sig)>5&&sz0>0");
	ntqq->Project("hqp", "log10(sz0)","abs(sz0sig)>5&&sz0>0");

	ha0pz0->Add(hbpz0);
	ha0pz0->Add(hcpz0);
	ha0pz0->Add(hqpz0);

	// normalize
	hbpz0->Divide(ha0pz0);
	hcpz0->Divide(ha0pz0);
	hqpz0->Divide(ha0pz0);

	cout << "hbp/hcp/hqp written." << endl;

	// hbn/hcn/hqn/ha0n
	ntbb->Project("hbn", "log10(-sz0)","abs(sz0sig)>5&&sz0<0");
	ntcc->Project("hcn", "log10(-sz0)","abs(sz0sig)>5&&sz0<0");
	ntqq->Project("hqn", "log10(-sz0)","abs(sz0sig)>5&&sz0<0");

	ha0nz0->Add(hbnz0);
	ha0nz0->Add(hcnz0);
	ha0nz0->Add(hqnz0);

	// normalize
	hbnz0->Divide(ha0nz0);
	hcnz0->Divide(ha0nz0);
	hqnz0->Divide(ha0nz0);

	cout << "hbn/hcn/hqn written." << endl;

	// hbip/hcip/hqip not normalized
	ntbb->Project("hbip", "sz0sig","abs(sz0sig)<5");
	ntcc->Project("hcip", "sz0sig","abs(sz0sig)<5");
	ntqq->Project("hqip", "sz0sig","abs(sz0sig)<5");

	cout << "hbip/hcip/hqip written." << endl;

	ntqq->Project("hjprob", "abs(sz0sig)","jprobcut>0");
	ntqq->Project("hjprob2", "abs(sz0sig)","jprobcut>0");

	integ = hjprobz0->Integral(0,50) + hjprob2z0->Integral(0,195);
	hjprobz0->Scale(1./integ);
	hjprob2z0->Scale(1./integ);

	cout << "hjprob/hjprob2 written." << endl;

	fz0->Write();
	fz0->Close();

	cout << "All finished." << endl;
}
