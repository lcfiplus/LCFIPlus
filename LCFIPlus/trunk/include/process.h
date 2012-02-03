// process.h

#ifndef process_h
#define process_h 1

#include "lcfiplus.h"

namespace lcfiplus{
	struct TrackSelectorConfig;

	class PrimaryVertexFinder : public Algorithm
	{
	public:
		PrimaryVertexFinder(){}
		virtual ~PrimaryVertexFinder(){}

		void init(Parameters *param);
		void process();
		void end();

		ClassDef(PrimaryVertexFinder,1);

	private:
		vector<Vertex *> * _vertex;	//!

		// parameters
		double _chi2th;

		// track cut parameters
		TrackSelectorConfig *_secVtxCfg; //!
	};

	class BuildUpVertex : public Algorithm
	{
	public:
		BuildUpVertex(){}
		virtual ~BuildUpVertex(){}

		void init(Parameters *param);
		void process();
		void end();

		ClassDef(BuildUpVertex,1);

	private:
		std::vector<Vertex *> * _vertices;	//!

		// parameters
		std::string _primvtxcolname;

		// vertex formation limits
		double _chi2thpri;
		double _chi2thsec;
		double _massth;
		double _posth;
		double _chi2orderinglimit;

		// association parameters
		bool _doassoc;
		double _minimumdist;
		double _chi2ratio;

		// track cut parameters
		TrackSelectorConfig *_secVtxCfg; //!
	};

	class JetClustering : public Algorithm
	{
	public:
		JetClustering(){}
		virtual ~JetClustering(){}

		void init(Parameters *param);
		void process();
		void end();

		ClassDef(JetClustering,1);

	private:
		VertexVec * _vertices; //!
		map<double, vector<Jet *> * > _jetsmap; //!
		vector<int> _njets;
		vector<double> _ycut;
		bool _useMuonID;
		double _vsMinDist;
		double _vsMaxDist;
		double _vsK0MassWidth;
		string _vcolname;

	};

}

#endif
