#ifndef PREDICT_PEXEL_HH_INCLUDED
#define PREDICT_PEXEL_HH_INCLUDED

#include <string>
#include <utility>
#include <GHMM/ghmm.hh>

#define RLE_PATTERN
#define KLD_PATTERN


// #define VERSION 1
#define VERSION 2

// #define SIGNALP_MODEL
// #define HALDAR_MOTIF

#define ARRAYLEN(x) (sizeof(x) / sizeof(x[0]))
#define ENDOF(x) ((x) + ARRAYLEN(x))

#endif
