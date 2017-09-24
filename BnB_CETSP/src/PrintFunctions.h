#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H

#include "Data.h"

void printDataToFile( Data *, char*, double );
void printDataToFile ( Data *, char*, double, int, double, double, double, double, int, int, double, double, double );
void printDataToFile( Data *, char*, double, int, double, double, double, double, double, double, int, int, double, double, double );
void printDataToMatlab( Data *, int, double, double, vector< int >, vector< vector< double > > );


#endif /* ifndef PRINTFUNCTIONS_H */
