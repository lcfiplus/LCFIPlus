#include "EventStore.h"
#include "lcfiplus.h"

#include <iostream>

using namespace std;

namespace lcfiplus {

	// ctor/dtor
	EventStore::EventStore(){}
	EventStore::~EventStore()
	{
		map<string, lcfiplus::EventStore::StoredEntry >::iterator itMap;
		// free all void * in the map
		for(itMap = _objectMap.begin(); itMap != _objectMap.end(); itMap ++){
			TClass *cl = TClass::GetClass(itMap->second.classname.c_str());
			cl->Destructor(itMap->second.obj);
			itMap->second.obj = NULL;
		}
	}

	// check existence of the named buffer
	bool EventStore::IsExist(const char *name)const
	{
		return _objectMap.find(name) != _objectMap.end();
	}

	// obtain class name of the named buffer
	const char * EventStore::GetClassName(const char *name)const
	{
		_itMap = _objectMap.find(name);
		if(_itMap == _objectMap.end())return NULL;
		else return _itMap->second.classname.c_str();
	}

	// provide read-only buffer
	void * EventStore::GetObject(const char *name)const
	{
		_itMap = _objectMap.find(name);
		if(_itMap == _objectMap.end())return NULL;
		else return _itMap->second.obj;
	}

	void * const& EventStore::GetObjectRef(const char *name)const
	{
		_itMap = _objectMap.find(name);
		if(_itMap == _objectMap.end())throw(new Exception("EventStore::GetObjectRef: Object not found"));

		return _itMap->second.obj;
	}

	// register new object
	void * EventStore::RegisterObject(const char * name, const char *classname, int flags)
	{
		// if already exists: return NULL
		_itMap = _objectMap.find(name);
		if(_itMap != _objectMap.end())return NULL;

		// make the class
		void *newobj;
		if(classname[strlen(classname)-1] == '*') newobj = new (void *);
		else newobj = TClass::GetClass(classname)->New();

		if(!newobj)return NULL;

		// register the map
		_objectMap[name] = StoredEntry(string(classname), newobj, flags);

		return newobj;
	}

	// print object map
	void EventStore::Print()const{
		cout << "Printing object map in EventStore..." << endl;
		for(_itMap = _objectMap.begin();_itMap != _objectMap.end();_itMap++){
			cout << "Name: " << _itMap->first << ", Class: " << _itMap->second.classname << endl;
		}
		cout << "Printing object map in EventStore finished." << endl;
	}

}
