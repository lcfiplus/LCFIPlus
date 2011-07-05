// FlavtagCategory.h

#ifndef FlavtagCategory_h
#define FlavtagCategory_h 1

struct FlavtagCategory {
	TString name;
	TString cut;
	std::vector<std::string> vars;
	std::vector<std::string> spec;
	void AddVariable(std::string s, char c) { vars.push_back(s); }
	void AddSpectator(std::string s) { spec.push_back(s); }
};

struct FlavtagType {
	TString name;
	TString cut;
};

#endif
