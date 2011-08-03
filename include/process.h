// process.h

#ifndef process_h
#define process_h 1

#include "interface.h"
#include "flavtag.h"
//#include "LcfiInterface.h"

namespace flavtag{
	struct SecondaryVertexConfig;

	class BuildUpVertex : public FlavtagAlgorithm
	{
	public:
		BuildUpVertex(){}
		virtual ~BuildUpVertex(){}

		void init(FlavtagParameters *param);
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
		SecondaryVertexConfig *_secVtxCfg; //!
	};

	class JetClustering : public FlavtagAlgorithm
	{
	public:
		JetClustering(){}
		virtual ~JetClustering(){}

		void init(FlavtagParameters *param);
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
