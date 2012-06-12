#ifndef __CONFIG_H__
#define __CONFIG_H__

class Config {
	int nreg;
	int maxk;
	int dev;
	int dump;
	int prec;
	int soltype;
	double *kappa;
	double *Ip;
	char meshfn[1024];
	char outfn[1024];
public:	
	explicit Config(const char *cfgfile);
	~Config() {
		delete[] kappa;
		delete[] Ip;
	}
	double getKappa(int region) const {
		return kappa[region];
	}
	double getIp(int region) const {
		return Ip[region];
	}
	const char *getMeshFilename() const {
		return meshfn;
	}
	const char *getOutFilename() const {
		return outfn;
	}
	int getMaxK() const {
		return maxk;
	}
	int getDevice() const {
		return dev;
	}
	bool doDump() const {
		return dump != 0;
	}
	bool usePrec() const {
		return prec != 0;
	}
	int solType() const {
		return soltype;
	}
};

#endif
