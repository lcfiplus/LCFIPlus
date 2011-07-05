#include "EventStore.h"
#include <iostream>

using namespace std;

namespace flavtag {

	// singleton operations
	EventStore * EventStore::_theInstance = NULL;
	EventStore * EventStore::Instance()
	{
		if(_theInstance == NULL) _theInstance = new EventStore;
		return _theInstance;
	}
	
	// ctor/dtor
	EventStore::EventStore(){}
	EventStore::~EventStore()
	{
		// free all void * in the map
		for(_itMap = _objectMap.begin(); _itMap != _objectMap.end(); _itMap ++){
			TClass *cl = TClass::GetClass(_itMap->second.classname.c_str());
			cl->Destructor(_itMap->second.obj);
			_itMap->second.obj = NULL;
		}
		// free the singleton
		_theInstance = NULL;
	}

	// check existence of the named buffer
	bool EventStore::IsExist(const char *name)
	{
		return _objectMap.find(name) != _objectMap.end();
	}

	// obtain class name of the named buffer
	const char * EventStore::GetClassName(const char *name)
	{
		_itMap = _objectMap.find(name);
		if(_itMap == _objectMap.end())return NULL;
		else return _itMap->second.classname.c_str();
	}

	// provide read-only buffer
	void * EventStore::GetObject(const char *name)
	{
		_itMap = _objectMap.find(name);
		if(_itMap == _objectMap.end())return NULL;
		else return _itMap->second.obj;
	}

	// register new object
	void * EventStore::RegisterObject(const char * name, const char *classname, int flags)
	{
		// if already exists: return NULL
		_itMap = _objectMap.find(name);
		if(_itMap != _objectMap.end())return NULL;

		// make the class
		void *newobj = TClass::GetClass(classname)->New();
		if(!newobj)return NULL;

		// register the map
		_objectMap[name] = StoredEntry(string(classname), newobj, flags);

		return newobj;
	}

	// print object map
	void EventStore::Print(){
		cout << "Printing object map in EventStore..." << endl;
		for(_itMap = _objectMap.begin();_itMap != _objectMap.end();_itMap++){
			cout << "Name: " << _itMap->first << ", Class: " << _itMap->second.classname << endl;
		}
		cout << "Printing object map in EventStore finished." << endl;
	}

}
