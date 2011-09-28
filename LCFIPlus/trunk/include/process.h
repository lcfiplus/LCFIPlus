// process.h

#ifndef process_h
#define process_h 1

#include "interface.h"
#include "lcfiplus.h"

namespace lcfiplus{
	struct TrackSelectorConfig;

	class PrimaryVertexFinder : public LcfiplusAlgorithm
	{
	public:
		PrimaryVertexFinder(){}
		virtual ~PrimaryVertexFinder(){}

		void init(LcfiplusParameters *param);
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

	class BuildUpVertex : public LcfiplusAlgorithm
	{
	public:
		BuildUpVertex(){}
		virtual ~BuildUpVertex(){}

		void init(LcfiplusParameters *param);
		void process();

		ClassDef(BuildUpVertex,1);

	private:
		std::vector<Vertex *> * _vertices;	//!

		// parameters

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

	class JetClustering : public LcfiplusAlgorithm
	{
	public:
		JetClustering(){}
		virtual ~JetClustering(){}

		void init(LcfiplusParameters *param);
		void process();
		void end();

		ClassDef(JetClustering,1);

	private:
		const std::vector<Vertex *> * _vertices; //!
		std::vector<Jet *> * _jets; //!
		int _njets;
		double _ycut;
		bool _useMuonID;
		double _vsMinDist;
		double _vsMaxDist;
		double _vsK0MassWidth;
	};

}

#endif
