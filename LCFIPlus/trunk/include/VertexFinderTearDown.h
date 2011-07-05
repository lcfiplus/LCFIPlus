// VertexFinderTearDown.h

#ifndef VertexFinderTearDown_h
#define VertexFinderTearDown_h 1

#include "flavtag.h"
#include "VertexFitterLCFI.h"
#include "VertexFitterSimple.h"
#include "algoEtc.h"

#include <list>
#include <vector>

using namespace std;

namespace flavtag {

	class SortTracksByChi2{ // decending order
		public:
			bool operator() (const pair<Track *, float> &p1, const pair<Track *, float> &p2){return p1.second > p2.second;}
	};

	// Function for recursive search of vertices using TearDown method
	vector<flavtag::Vertex*> * findTearDownVertices(const Event& evt, const Jet& jet);

	// Primary Vertex finder with TearDown method
	flavtag::Vertex * findPrimaryVertex(const vector<Track *> &tracks, double chi2 = 9.0);
//	flavtag::Vertex * findPrimaryVertex(const vector<Track *> &tracks, const vector<Track *> &beamTracks, double chi2 = 9.0);

	// implementation of TearDown method
	template<template<class T, class Allocator=allocator<T> > class Container = std::vector, template<class Iterator> class VertexFitter = VertexFitterLCFI >
		class VertexFinderTearDown{
		public:
			Vertex * operator () (const Container<Track *> &tracks, const Container<Track *> *fixedTracks = 0,
														double chiSquareThreshold = 9.0, Container<Track *> *residual = 0, Vertex *pointConstraint = 0)
			{
				// copy tracks into a list
				list<Track *>trackList;
				trackList.resize(tracks.size() + (fixedTracks ? fixedTracks->size() : 0));
				list<Track *>::iterator listIt = copy(tracks.begin(), tracks.end(), trackList.begin());
				if(fixedTracks){
					copy(fixedTracks->begin(), fixedTracks->end(), listIt);
				}

				Vertex *resultVertex = 0;
				float worstChi2;
				while(trackList.size() >= 2){
					resultVertex = VertexFitter<list<Track *>::iterator>() (trackList.begin(), trackList.end(), pointConstraint);
					Track *worstTrack = resultVertex->getWorstTrack();
					if(fixedTracks && find(fixedTracks->begin(), fixedTracks->end(), worstTrack) != fixedTracks->end()){
						// sort Chi2Tracks
						vector<pair<Track *, float> > vpair;
						const map<Track *, float> &mpair = resultVertex->getTracksChi2Map();
						vpair.resize(mpair.size());
						partial_sort_copy(mpair.begin(), mpair.end(), vpair.begin(), vpair.end(), SortTracksByChi2());

						unsigned int nworst = 1;
						do{
							worstTrack = vpair[nworst++].first;
						}while(nworst < vpair.size() && find(fixedTracks->begin(), fixedTracks->end(), worstTrack) != fixedTracks->end());

						cout << "The worst track is fixed, " << nworst << "th track will be removed." << endl;
					}

					worstChi2 = resultVertex->getChi2Track(worstTrack);
					if(worstChi2 > chiSquareThreshold){
						trackList.remove(worstTrack);
						if(residual)
							residual->push_back(worstTrack);
						/*
						cout << "Track removed: chi2 = " << worstChi2 << ", vpos = ("
									<< resultVertex->getX() << ", " << resultVertex->getY() << ", " << resultVertex->getZ() << ")" <<  endl;
						 */
						delete resultVertex;
						resultVertex = 0;
					}
					else
						break; // all tracks below chi2 threshold
				}
				
				return resultVertex;
			}
		};
}

#endif
