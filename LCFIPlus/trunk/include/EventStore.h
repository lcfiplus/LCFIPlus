#ifndef EventStore_h
#define EventStore_h 1

#include <map>
#include <string>
#include <iostream>
#include "TClass.h"

#define EVENTSTORE_VERSION 0

using namespace std;

namespace lcfiplus {

 	class EventStore{
			/*
				EventStore: a simple named storage for event data

				To obtain data: use Get(name) (in compiled code) or GetObject(name) (in CINT)
					CAUTION: GetObject() does not have type-check.

				To register data: use Register(name) or RegisterObject(name, class)
					Use a pointer returned (cannot register a pointer from the caller).
					Any classes or vector of classes with ClassDef() can be registered.

				2010/07/07 suehara Initial version
			*/

		public:
			// singleton is tranfered to Event class
//			static EventStore * Instance();

			// register flags
			enum{PERSIST = 1, JET_EXTRACT_VERTEX = 2};

			// Check collection
			bool IsExist(const char *name)const;

			// Object retrieval
			const char * GetClassName(const char *name)const;
			void * GetObject(const char *name)const; // CAUTION: no type check
			template<typename T> bool Get(const char *name, const vector<T*> *& buf)const;	// for pointer-vector classes
			template<typename T> bool Get(const char *name, const vector<T> *& buf)const;	// for vector classes
			template<typename T> bool Get(const char *name, const T*& buf)const;						// for non-vector

			// Object registration
			void * RegisterObject(const char *name, const char *classname, int flags = 0);
			template<typename T> bool Register(const char *name, vector<T*>* &buf, int flags = 0);	// for pointer-vector classes
			template<typename T> bool Register(const char *name, vector<T>* &buf, int flags = 0);	// for vector classes
			template<typename T> bool Register(const char *name, T* &buf, int flags = 0);					// for non-vector

			// print object map
			void Print()const;

			// public dtor
			virtual ~EventStore();

			// for persistency classes /////////////////////////////
			// internal storing structure
			struct StoredEntry{
				string classname;
				void * obj;
				int flag;

				StoredEntry(string cn = string(""), void *ob = 0, int fl = 0) : classname(cn), obj(ob), flag(fl){}
			};

			// map accessor
			const map<string, lcfiplus::EventStore::StoredEntry > & GetObjectMap()const{return _objectMap;}

		protected:
			// reference retrieval: internal use
			void * const& GetObjectRef(const char *name)const; // CAUTION: no type check
			// constructor not public: use Event class
			EventStore();

    private:
			map<string, lcfiplus::EventStore::StoredEntry > _objectMap;
			mutable map<string, lcfiplus::EventStore::StoredEntry >::const_iterator _itMap;

  };


	// template implementations
	template<typename T> bool EventStore::Get(const char *name, const T*&buf)const{
		_itMap = _objectMap.find(name);
		if(_itMap == _objectMap.end() || string(_itMap->second.classname) != TClass::GetClass(typeid(T))->GetName())return false;
		buf = static_cast<const T*>(GetObject(name));
		return true;
	}
	// vector version
	template<typename T> bool EventStore::Get(const char *name, const vector<T>* &buf)const{
		const char *elemclasname = TClass::GetClass(typeid(T))->GetName();
		string vectclasname = "vector<";
		vectclasname += elemclasname;
		vectclasname += ">";

		_itMap = _objectMap.find(name);

		if(_itMap == _objectMap.end() || _itMap->second.classname != vectclasname)return false;
		buf = static_cast<const vector<T>*>(GetObject(name));

		return true;
	}
	// pointer-vector version
	template<typename T> bool EventStore::Get(const char *name, const vector<T*>* &buf)const{
		const char *elemclasname = TClass::GetClass(typeid(T))->GetName();
		string vectclasname = "vector<";
		vectclasname += elemclasname;
		vectclasname += "*>";

		_itMap = _objectMap.find(name);

		if(_itMap == _objectMap.end() || _itMap->second.classname != vectclasname)return false;
		buf = static_cast<const vector<T*>*>(GetObject(name));

		return true;
	}

	template<typename T> bool EventStore::Register(const char *name, T *&buf, int flags){
		string classname = TClass::GetClass(typeid(T))->GetName();
		classname += "*";
		buf = static_cast<T*>(RegisterObject(name, classname.c_str(), flags));
		return buf;
	}
	template<typename T> bool EventStore::Register(const char *name, vector<T> *&buf, int flags){
		const char *elemclasname = TClass::GetClass(typeid(T))->GetName();
		string vectclasname = "vector<";
		vectclasname += elemclasname;
		vectclasname += ">";
		cout << vectclasname << endl;
		buf = static_cast<vector<T>*>(RegisterObject(name, vectclasname.c_str(), flags));
		return buf;
	}
	template<typename T> bool EventStore::Register(const char *name, vector<T*> *&buf, int flags){
		const char *elemclasname = TClass::GetClass(typeid(T))->GetName();
		string vectclasname = "vector<";
		vectclasname += elemclasname;
		vectclasname += "*>";
		cout << vectclasname << endl;
		buf = static_cast<vector<T*>*>(RegisterObject(name, vectclasname.c_str(), flags));
		return buf;
	}

}

#endif
