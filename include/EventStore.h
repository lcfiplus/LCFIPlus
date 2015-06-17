#ifndef EventStore_h
#define EventStore_h 1

#include <map>
#include <list>
#include <string>
#include <iostream>
#include "TClass.h"

#define EVENTSTORE_VERSION 0

using namespace std;

namespace lcfiplus {

// observer class
class EventStoreObserver {
 public:
  EventStoreObserver() {}
  virtual ~EventStoreObserver() {}

  virtual void GetCallback(const char* name, const char* classname) {}
  virtual void RegisterCallback(const char* name, const char* classname, int flags) {}

  ClassDef(EventStoreObserver,1)
};

/**
	A simple named storage for event data.

	To obtain data: use Get(name) (in compiled code) or GetObject(name) (in CINT)

	CAUTION: GetObject() does not have type-check.

	To register data: use Register(name) or RegisterObject(name, class)
	Use a pointer returned (cannot register a pointer from the caller).
	Any classes or vector of classes with ClassDef() can be registered.

	@author T. Suehara, ICEPP, The University of Tokyo
	@version $Id$
 */
class EventStore {

 public:
  // singleton is tranfered to Event class
//			static EventStore * Instance();

  // register flags
  enum {
    DO_NOT_DELETE = 0x0001, // do not delete pointed object at ClearObjects(): to be used in reference collections
    PERSIST = 0x0002, // lcio: to be saved in the output
    JET_WRITE_VERTEX = 0x1000
  };

  // Observer registration
  void RegisterObserver(EventStoreObserver* observer) {
    _observerList.push_back(observer);
  }
  void UnregisterObserver(EventStoreObserver* observer) {
    _observerList.remove(observer);
  }

  // Check collection
  int Count(const char* name)const;
  bool IsExist(const char* name)const {
    return Count(name);
  }
  bool IsExist(const char* name, const char* classname)const;

  // Object retrieval
  const char* GetClassName(const char* name, int idx = 0)const;
  void* GetObject(const char* name, const char* classname = "")const;  // CAUTION: no type check
  template<typename T> bool Get(const char* name, const vector<const T*>*& buf)const;	// for pointer-vector classes, add const
  template<typename T> bool Get(const char* name, const vector<T*>*& buf)const;	// non-const pointer vector prohibited: invoke error
  template<typename T> bool Get(const char* name, const vector<T>*& buf)const;	// for vector classes
  template<typename T> bool Get(const char* name, const T*& buf)const;						// for non-vector

  // Object registration
  void* RegisterObject(const char* name, const char* classname, int flags = 0);
  template<typename T> bool Register(const char* name, vector<T*>*& buf, int flags = 0);	// for pointer-vector classes
  template<typename T> bool Register(const char* name, vector<T>*& buf, int flags = 0);	// for vector classes
  template<typename T> bool Register(const char* name, T*& buf, int flags = 0);					// for non-vector

  // print object map
  void Print()const;

  // clear current objects: to switch event
  void ClearObjects();

  // public dtor
  virtual ~EventStore();

  // for persistency classes /////////////////////////////
  // internal storing structure
  struct StoredEntry {
    string classname;
    void* obj;
    int flag;

    StoredEntry(string cn = string(""), void* ob = 0, int fl = 0) : classname(cn), obj(ob), flag(fl) {}
  };

  // map accessor
  const multimap<string, lcfiplus::EventStore::StoredEntry >& GetObjectMap()const {
    return _objectMap;
  }

 protected:
  // reference retrieval: internal use
  void* const& GetObjectRef(const char* name, const char* classname = "")const;  // CAUTION: no type check
  // constructor not public: use Event class
  EventStore();

 private:
  multimap<string, lcfiplus::EventStore::StoredEntry > _objectMap;
  mutable multimap<string, lcfiplus::EventStore::StoredEntry >::const_iterator _itMap;

  list<EventStoreObserver*> _observerList;
};


// template implementations
template<typename T> bool EventStore::Get(const char* name, const T*& buf)const {
  buf = static_cast<const T*>(GetObject(name, TClass::GetClass(typeid(T))->GetName()));
  return buf;
}
// vector version
template<typename T> bool EventStore::Get(const char* name, const vector<T>*& buf)const {
  const char* elemclasname = TClass::GetClass(typeid(T))->GetName();
  string vectclasname = "vector<";
  vectclasname += elemclasname;
  vectclasname += ">";

  buf = static_cast<const vector<T>*>(GetObject(name, vectclasname.c_str()));
  return buf;
}

// non-const pointer-vector prohibited
template<typename T> bool EventStore::Get(const char* name, const vector<T*>*& buf)const {
  T a = "abc"; // invoke compiler error
  throw ("EventStore::Get: non-const pointer-vector prohibited");
  return false;
}

// pointer-vector version
template<typename T> bool EventStore::Get(const char* name, const vector<const T*>*& buf)const {
  const char* elemclasname = TClass::GetClass(typeid(T))->GetName();
  string vectclasname = "vector<";
  vectclasname += elemclasname;
  vectclasname += "*>";

  buf = static_cast<const vector<const T*>*>(GetObject(name, vectclasname.c_str()));
  return buf;
}

template<typename T> bool EventStore::Register(const char* name, T*& buf, int flags) {
  string classname = TClass::GetClass(typeid(T))->GetName();
  classname += "*";
  buf = static_cast<T*>(RegisterObject(name, classname.c_str(), flags));
  return buf;
}
template<typename T> bool EventStore::Register(const char* name, vector<T>*& buf, int flags) {
  const char* elemclasname = TClass::GetClass(typeid(T))->GetName();
  string vectclasname = "vector<";
  vectclasname += elemclasname;
  vectclasname += ">";
  buf = static_cast<vector<T>*>(RegisterObject(name, vectclasname.c_str(), flags));
  return buf;
}
template<typename T> bool EventStore::Register(const char* name, vector<T*>*& buf, int flags) {
  const char* elemclasname = TClass::GetClass(typeid(T))->GetName();
  string vectclasname = "vector<";
  vectclasname += elemclasname;
  vectclasname += "*>";
//		cout << vectclasname << endl;
  buf = static_cast<vector<T*>*>(RegisterObject(name, vectclasname.c_str(), flags));
  return buf;
}

}

#endif
