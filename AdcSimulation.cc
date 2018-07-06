#include <TFile.h>
#include <cmath>
#include <vector>
#include <TF1.h>
#include <TTree.h>
using namespace std;

double PreampResp(double* x, double* par){
	if(x[0]>0 && x[10]<10)
		return 4.31054*exp(-2.94809*x[0]/par[1])*par[0]-2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*par[0] \
							-2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*cos(2.38722*x[0]/par[1])*par[0] \
										+0.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*par[0] \
													+.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*cos(5.18561*x[0]/par[1])*par[0] \
																+0.762456*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0] \
																			-0.762456*exp(-2.82833*x[0]/par[1])*cos(2.38722*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0] \
																						+0.762456*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0] \
																									-2.6202*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0] \
																												-0.327684*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0] + \
																															+0.327684*exp(-2.40318*x[0]/par[1])*cos(5.18561*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0] \
																																		-0.327684*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0] \
																																					+0.464924*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0];
	else
		return 0;                                                                   
}

class DAC{ // A simple pulser

		double _Vstart;
		double _Vend;
		int _nbit;
		vector<double> _Vdac;

	public:
		DAC(double Vstart=0, double Vend=1.2, int nbit=6){
			_Vstart = Vstart;
			_Vend = Vend;
			_nbit = nbit;
		}

		vector<double> settings(){
			if(_Vdac.size()==0){
				int nsize = int( pow(2,_nbit) );
				for(int i=1; i<nsize; i++){
					double Vp = _Vstart + i * (_Vend - _Vstart)/(nsize-1);
					_Vdac.push_back(Vp);
				}
			}
			return _Vdac;
		}
};

class Preamp{
	double _capacitor;
	double _max_qin;
	double _max_vout;
	double _amp_gain;
	double _peak_time;
	TF1* _wf; // waveform
	double _nsamp; // num of sample in a cycle

public:
	Preamp(double F=4.7, double pt=3.0, double nsamp=100){
		_capacitor = 183;
		_max_qin = 100; // fC
		_max_vout = 1.4; // volts
		_amp_gain = F;
		_peak_time = pt;
		_wf = new TF1("waveform", PreampResp, 0, 10, 2);
		_nsamp = nsamp;
	}

	tuple< vector<double>, vector<double> > GetWaveform(double Vin){
		double Qin = Vin * _capacitor;
		if(Qin>_max_qin){Qin = _max_qin;}
		double Vamp = _amp_gain * Qin; // mV
		Vamp *= 0.001; // volts
		_wf->SetParameter(0, Vamp *10*1.012);//10.12 is for scaling
		_wf->SetParameter(1, _peak_time);
		vector<double> x,y;
		for(int j=0; j<_nsamp; j++){
			double time = 10.0 / _nsamp * j;
			double volt = _wf->Eval(time);
			if(volt > _max_vout) volt = _max_vout;
			cout << "time = " << time << " ; voltage= " << volt << endl;
			x.push_back(time);
			y.push_back(volt);
		}
		return make_tuple(x,y);
	}
};

class FESetting{
	public:
		int _ID;
		double _Vdac;
		double _Gamp;

		FESetting(int id, double v, double g){
			_ID=id;
			_Vdac = v;
			_Gamp = g;
		}
};

double applyNL(double Atrue){
	double dA; // Atrue - Ameas
	if(Atrue<1500){
		dA = Atrue / 75.0;
	}
	else{
		dA = -(Atrue - 4000) / 125.0;
	}
	double Ameas = Atrue - dA;
	if(Ameas>4096) Ameas = 4096;
	return Ameas;
}

double applyNL2(double Atrue){
      double dA;
      dA = -22.0 + Atrue /158.3;

 
	return Atrue - dA; 
}

class SimHeader {
public:
	int Nset;
	double Vdac;
	double Gamp;
	double time;
	double Vtrue;
	double Atrue;
	double Ameas;
};

int main(){
	TFile* ofile = new TFile("protodune-nl-MC.root","recreate");
	TTree* T = new TTree("T","T");
	// TNtuple("T","T","Nset:Vdac:Gamp:time:Vtrue:Atrue:Ameas");
	SimHeader* simheader = new SimHeader();
	T->Branch("Nset", &simheader->Nset, "Nset/I");
	T->Branch("Vdac", &simheader->Vdac, "Vdac/D");
	T->Branch("Gamp", &simheader->Gamp, "Gamp/D");
	T->Branch("time", &simheader->time, "time/D");
	T->Branch("Vtrue", &simheader->Vtrue, "Vtrue/D");
	T->Branch("Atrue", &simheader->Atrue, "Atrue/D");
	T->Branch("Ameas", &simheader->Ameas, "Ameas/D");

	DAC* aDAC = new DAC();
	vector<double> dacSetting = aDAC->settings();
	vector<FESetting> FESets;
	for(int i=0; i<dacSetting.size(); i++){ // fill settings
		FESets.push_back(FESetting(0, dacSetting[i], 4.7));
		FESets.push_back(FESetting(0, dacSetting[i], 7.8));
		FESets.push_back(FESetting(0, dacSetting[i], 14.));
		FESets.push_back(FESetting(0, dacSetting[i], 25.));
	}
	// sort the FE settings in order of waveform amplitude
	sort(FESets.begin(), FESets.end(), [](FESetting& a, FESetting& b) {return (a._Vdac * a._Gamp) < (b._Vdac * b._Gamp);} );
	
	for(int i=0; i<FESets.size(); i++){
		FESetting fe = FESets[i];
		fe._ID = i+1;
		Preamp* aPreamp = new Preamp(fe._Gamp, 3.0, 100); // 3us + nsamp=100
		vector<double> time, wfamp;
		tie(time, wfamp) = aPreamp->GetWaveform(fe._Vdac);
		for(int j=0; j<time.size(); j++){
			//T->Fill(fe._ID, fe._Vdac, fe._Gamp, time[j], wfamp[j], wfamp[j]*2900., applyNL(wfamp[j]*2900.) );
			simheader->Nset = fe._ID;
			simheader->Vdac = fe._Vdac;
			simheader->Gamp = fe._Gamp;
			simheader->time = time[j];
			simheader->Vtrue = wfamp[j];
			simheader->Atrue = wfamp[j]*2900.;
			simheader->Ameas = applyNL(wfamp[j]*2900.);
			T->Fill();
		}
		delete aPreamp;
	}
	delete aDAC;

	T->Write();
	ofile->Close();
}
