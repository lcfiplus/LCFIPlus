#!/bin/sh
rm -f dict.cc dict.h

rootcint -f dict.cc -c -I$LCIO/include ../include/flavtag.h ../include/Driver.h ../include/Suehara.h ../include/EventStore.h ../include/LCIOStorer.h ../include/TreeStorer.h ../include/EventNavigator.h ../include/JetFinder.h ../include/interface.h ../include/process.h ../include/MakeNtuple.h ../include/TrainMVA.h LinkDef.h
