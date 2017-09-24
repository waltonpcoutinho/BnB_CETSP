#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cstdlib>
#include <stdio.h>
#include <cfloat>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>

#include"PrintFunctions.h"
#include"util.h"
#include "Data.h"

void printDataToFile ( Data * dataptr, char* option, double overlap )
{
   ofstream outFile ( "Resultados/resultados_BnB_CETSP_OPTF.txt", ios::app );

   if ( !outFile ) {
      cout << "arquivo nao pode ser criado\n";
      exit (1);
   }

   if ( !outFile ) {
      cout << "arquivo nao pode ser criado\n";
      exit (1);
   }

   outFile << fixed << setiosflags ( ios::showpoint ) << setprecision( 2 );

   outFile  << dataptr->fileName << " "
      << option << " "
      << overlap << " "
      << dataptr->getOverlapRatio() << " ";
   outFile << endl;
   outFile.close();		
}

void printDataToFile ( Data * dataptr, char* option, double overlap, int sizeInst, double bestKnown, double best, double best_lb, double gap_root, int count_SOCP_solved, int itCount, double computationTime, double sbComputationTime, double totalSocpCompTime )
{
   ofstream outFile ( "Resultados/resultados_BnB_CETSP_OPTF.txt", ios::app );

   if ( !outFile ) {
      cout << "arquivo nao pode ser criado\n";
      exit (1);
   }

   outFile << fixed << setiosflags ( ios::showpoint ) << setprecision( 15 );

   outFile  << "OPTF" << " "
      << dataptr->fileName << " "
      << option << " "
      << overlap << " "
      << dataptr->getOverlapRatio() << " "
      << sizeInst << " "
      << bestKnown << " "
      << best << " ";

   if(best_lb >= 999999999999999) outFile << " - ";
   else outFile << best_lb << " ";

   outFile
      << gap_root << " "
      << count_SOCP_solved << " "
      << itCount << " "
      << computationTime << " "
      << sbComputationTime << " "
      << totalSocpCompTime << " ";
   outFile << endl;
   outFile.close();		
}

void printDataToFile( Data * dataptr, char* option, double overlap, int sizeInst, double bestKnown, double ub, double best_lb, double gap_real, double gap_lb_bnb, double gap_root, int count_SOCP_solved, int itCount, double computationTime, double sbComputationTime, double totalSocpCompTime )
{
   ofstream outFile ( "Resultados/resultados_BnB_CETSP_noOPTF.txt", ios::app );

   if ( !outFile ) {
      cout << "arquivo nao pode ser criado\n";
      exit (1);
   }

   outFile << fixed << setiosflags ( ios::showpoint ) << setprecision( 15 );

   outFile  << "noOPTF" << " "
      << dataptr->fileName << " "
      << option << " "
      << overlap << " "
      << dataptr->getOverlapRatio() << " "
      << sizeInst << " "
      << bestKnown << " "
      << ub << " "
      << best_lb << " "
      << gap_real << " "
      << gap_lb_bnb << " "
      << gap_root << " "
      << count_SOCP_solved << " "
      << itCount << " "
      << computationTime << " "
      << sbComputationTime << " "
      << totalSocpCompTime << " ";
   outFile << endl;
   outFile.close();		
}

void printDataToMatlab( Data * dataptr, int sizeInst, double overlap, double best, vector< int > solucao, vector< vector< double > > solucaoXYZ )
{
   char nomeArquivo[ 100 ] = {"Resultados/MATLAB/"};
   strcat (nomeArquivo, dataptr->fileName.c_str());
   strcat (nomeArquivo, ".txt");	
   ofstream outFile2 ( nomeArquivo );

   if ( !outFile2 ) {
      cout << "arquivo nao pode ser criado\n";
      exit (1);
   }

   outFile2 << fixed << setiosflags ( ios::showpoint ) << setprecision( 15 );

   outFile2	<< sizeInst << endl
      << overlap << endl
      << best << endl
      << solucao.size() << endl;
   for ( int i = 0; i < solucao.size(); i++ ){
      outFile2 << solucao[ i ] << endl;
   }
   for ( int j = 0; j < solucaoXYZ[ 0 ].size(); j++ ){
      outFile2 << solucaoXYZ[ 0 ][ j ] << " " << solucaoXYZ[ 1 ][ j ] << " " << solucaoXYZ[ 2 ][ j ] << endl;
   }
   outFile2.close();
}
