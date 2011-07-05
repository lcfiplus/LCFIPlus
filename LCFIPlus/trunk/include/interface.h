#ifndef interface_h
#define interface_h 1

#include "flavtag.h"

#include <vector>
#include <sstream>

#include "Rtypes.h"
#include "TObject.h"

using namespace std;

namespace flavtag {

	class FlavtagParameters
	{
	public:
		FlavtagParameters(bool as = true) : _allowstring(as){}
		~FlavtagParameters(){
			// TODO: delete map objects!
		}

		// fetch for non-vector
		template<typename T> void fetch(const char *key, T &ret, T def = T())
		{
			if(_map.find(key) == _map.end()){
				cout << "Parameter " << key << " not found." << endl;
				ret = def;
			}
			else if(_map[key].first == typeid(T).name())ret = *(T *)(_map[key].second);
			else if(_allowstring && _map[key].first == typeid(string).name()){
				istringstream str(*(string *)_map[key].second);
				str >> ret;
			}
			else if(_allowstring && _map[key].first == typeid(vector<string>).name()){
				istringstream str((*(const vector<string> *)_map[key].second)[0]);
				str >> ret;
			}
			else
				throw(Exception("Parameter type invalid."));
			return;
		}

		// fetch for vector
		template<typename T> void fetchArray(const char *key, vector<T> &ret)
		{
			if(_map.find(key) == _map.end()){
				cout << "Parameter " << key << " not found." << endl;
			}
			else if(_map[key].first == typeid(vector<T>).name())ret = *(T *)(_map[key].second);
			else if(_map[key].first == typeid(T).name())ret.push_back(*(T *)(_map[key].second));
			else if(_allowstring && _map[key].first == typeid(string).name()){
				ret.push_back(T());

				istringstream str(*(string *)_map[key].second);
				str >> ret[0];
			}
			else if(_allowstring && _map[key].first == typeid(vector<string>).name()){
				const vector<string> *svec = (const vector<string> *)_map[key].second;
				ret.resize(svec->size());
				for(unsigned int n=0;n<svec->size();n++){
					istringstream str((*svec)[n]);
					str >> ret[n];
				}
			}
			else
				throw(Exception("Parameter type invalid."));
		}

public:
		// for non-vector only (if string parameter)
		template<typename T> T get(const char *key, T def = T()){
			T ret;
			fetch(key, ret, def);
			return ret;
		}

		bool exist(const char *key){return _map.find(key) != _map.end();}

		template<typename T> void add(const char *key, T &data){
			if(_map.find(key) != _map.end())throw(Exception("Double entry."));

			_map[key] = pair<string, void *>(typeid(T).name(), new T(data));
		}

	private:
		map<string, pair<string, void *> > _map;
		bool _allowstring;
	}; 

	class FlavtagAlgorithm
	{
	public:
		FlavtagAlgorithm(){_param = 0;}
		virtual ~FlavtagAlgorithm(){}

		virtual void init(FlavtagParameters *param){
			_param = param;
		}
		virtual void process() = 0;
		virtual void end(){}

	protected:
		FlavtagParameters * GetParam()const{return _param;}
	private:
		FlavtagParameters *_param;

		ClassDef(FlavtagAlgorithm,1);
	};

}

#endif
