/*
 *    BnB_CETSP
 *    main.cpp
 *    Purpose: Solves the CETSP by branch-and-bound
 *
 *    @author Walton Pereira Coutinho
 *    @version 1.0 02/02/2014
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *    
 *    You should have received a copy of the GNU General Public License 
 *    along with this program. If not, see:
 *    <https://github.com/waltonpcoutinho/BnB_CETSP>.
*/

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

#include "util.h"
#include "PrintFunctions.h"
#include "Data.h"
#include "BranchNBound.h"
#include "SolveSocpCplex.h"
#include "structs.h"

using namespace std;

long double somaTeste = 0;

int main(int argc, char** argv) 
{
   // Check input 
   if ( argc != 9 ) {
      cout << "Wrong calling command\n";
      cout << "./exeCETSP [path da Instancia] [OPTIONS] [OVERLAP FACTOR] [TIME LIMIT] [BRANCHING RULE] [BRANCHING STRATEGY] [ROOT SELECTION] [S.B SIZE]" << endl;	
      exit( 1 );
   }

   char *arqInstancia, *option;
   arqInstancia = argv[ 1 ];
   option = argv[ 2 ];
   double overlap = 0;
   overlap = atof( argv[ 3 ] );
   int timeLimit = atoi( argv[ 4 ] );
   int branchingRule = 0;
   int selectingRoot = atoi( argv[ 7 ] );

   if( strcmp( argv[ 5 ], "V1" ) == 0 ) branchingRule = 1;
   if( strcmp( argv[ 5 ], "SB" ) == 0 ) branchingRule = 2;

   int branchingStrategy = 0;
   if( strcmp( argv[ 6 ], "DFS" ) == 0 ) branchingStrategy = 1;
   if( strcmp( argv[ 6 ], "BFS" ) == 0 ) branchingStrategy = 2;
   if( strcmp( argv[ 6 ], "BeFS" ) == 0 ) branchingStrategy = 3;


   Data *dataptr = new Data( arqInstancia, option, overlap, argc, argv );

   int sizeInst = dataptr->getSizeInst();
   cout << "Tamanho da instância: " << sizeInst << endl;
   cout << fixed << setiosflags ( ios::showpoint ) << setprecision( 8 );

   //starting the branch and bound
   double initialTotalTimeBnB = cpuTime();
   list < node* > open;
   list < node* >:: iterator itOpen;

   node * root = new node;

   double bestKnown = dataptr->getUbFromFile();
   cout << "Upper Bound Inicial: " << bestKnown << endl;
   double ub = bestKnown + 1;
   double best = DBL_MAX;
   double best_lb = DBL_MAX;
   double best_ub = DBL_MAX;
   unsigned long int count_SOCP_solved = 0;
   double rootLB = 0;

   int strong_branching_size = atoi( argv[ 8 ] );
   int precision = 4;
   int printHeader = 0;
   int optimalFound = 0;
   int print_counter = floor( sizeInst/2 );

   BranchNBound *bnbPtr = new BranchNBound( dataptr );

   //####################################################################
   vector< mpz_class > levels;
   levels.resize( sizeInst );
   mpz_class quantity = 0;
   mpz_t sizeOfTree;
   bnbPtr->computeSizeTree( sizeInst, sizeOfTree, levels );
   mpf_class temp = 0;
   //####################################################################	

   //### root node selection strategy ###
   if( selectingRoot == 1 ) root->pts = bnbPtr->selectRoot();
   if( selectingRoot == 2 ) root->pts = bnbPtr->selectRoot2();
   if( selectingRoot == 3 ) root->pts = bnbPtr->selectRoot3();

   cout << "Initial root: ";
   for ( int i = 0; i <root->pts.size() ; i++ ){
      cout << root->pts[ i ] << " ";
   }
   cout << endl;
   //####################################

   SolveSocpCplex *solveSocpPtr = new SolveSocpCplex( dataptr, root->pts );
   //solve model
   double totalSocpCompTime = 0;
   double initialSocpCompTime = cpuTime();
   solveSocpPtr->solveSOCP( root->pts );
   double finalSocpCompTime = cpuTime();
   totalSocpCompTime += ( finalSocpCompTime - initialSocpCompTime );
   count_SOCP_solved++;

   root->lb = solveSocpPtr->getF_value();
   rootLB = root->lb;
   root->s_lb = 1;
   solveSocpPtr->printF_value();
   solveSocpPtr->finishSOCP();
   solveSocpPtr->printSolution( root->pts );
   somaTeste += solveSocpPtr->violation;
   
   //get solution
   vector< vector< double > > solutionXYZ;
   vector< double > tempX;
   vector< double > tempY;
   vector< double > tempZ;

   for ( int i = 0; i < root->pts.size(); i++ ){
      tempX.push_back( solveSocpPtr-> getSolutionX( i ) );
      tempY.push_back( solveSocpPtr-> getSolutionY( i ) );
      tempZ.push_back( solveSocpPtr-> getSolutionZ( i ) );
   }

   solutionXYZ.push_back( tempX );
   solutionXYZ.push_back( tempY );
   solutionXYZ.push_back( tempZ );

   //check feasibility
   bool feasibilityTest;
   feasibilityTest = bnbPtr->check_feasibility_Q( root->pts, solutionXYZ );

   cout << endl;
   cout << "Não cobertos na raiz: ";
   for ( int i = 0; i < bnbPtr->notCoveredBalls.size(); i++ ){
      cout << bnbPtr->notCoveredBalls[ i ] << " ";
   }
   cout << endl;
   //check not covered clients
   for ( int j = 0; j < bnbPtr->notCoveredBalls.size(); j++ ){
      root->notCovered.push_back( bnbPtr->notCoveredBalls[ j ] );
   }	

   if ( feasibilityTest == true ){
      cout << "FEASIBLE ROOT" << endl;
      best_lb = rootLB;
      best = rootLB;
      best_ub = rootLB;
      double finalTotalTimeBnB = cpuTime();
      double computationTime = finalTotalTimeBnB - initialTotalTimeBnB;
      double gap_root = ( ( ub - rootLB )/ ub )*100;
      printDataToFile( dataptr, option, overlap, sizeInst, bestKnown, best, best_lb, 0, count_SOCP_solved, 0, computationTime, 0, computationTime );
      return 0;
   }
   else{
      cout << "INFEASIBLE ROOT" << endl;
   }
   cout << endl;

   solutionXYZ.clear();
   tempX.clear();
   tempY.clear();
   tempZ.clear();

   open.push_back( root );
   vector< int > solucao;
   vector< vector< double > > solucaoXYZ;
   int itCount = 0;
   int k = 0;
   double sbComputationTime = 0;

   //set branching rule
   int strongBranchingSize = 0;
   if ( branchingRule == 1 ){
      strongBranchingSize = 1;
   }
   if ( branchingRule == 2 ){
      strongBranchingSize = strong_branching_size;
   }

   int sum = 0;
   vector< sbAuxStruct* > vectorOfChildren;
   list< int >::iterator stBrchit;

   while( !open.empty() && cpuTime() - initialTotalTimeBnB <= timeLimit ){

      node * current;
      //######### Depth First Search ###########	
      if ( branchingStrategy == 1 ){
         current = open.back();
         open.pop_back();
      }
      //######### Depth First Search ###########

      //######### Breadth First Search ###########
      if ( branchingStrategy == 2 ){
         current = open.front();
         open.pop_front();
      }
      //######### Breadth First Search ###########

      //######### Best First Search ###########
      if ( branchingStrategy == 3 ){
         for ( itOpen = open.begin(); itOpen != open.end(); itOpen++ ){
            if ( ( *itOpen )->s_lb == 1 ){
               current = ( *itOpen );
            }
         }
         open.remove( current );
         if ( !open.empty() ){
            bnbPtr->computeLowerBounds( &open, current, &best_lb );
         }
      }
      //######### Best First Search ###########

      if ( current->notCovered.empty() || current->lb > ub ){
         // FATHOMED BY BOUND
         quantity += mpz_class( sizeOfTree )/levels[ current->pts.size() - 3 ] - 1;
         delete current;
      }
      else{
         stBrchit = current->notCovered.begin();

         //control strong branching size
         if ( branchingRule == 2 ){
            strongBranchingSize = strong_branching_size;
            strongBranchingSize -= current->pts.size() - 3;
         }
         if ( strongBranchingSize < 1 ){
            strongBranchingSize = 1;
         }
         //control strong branching size

         //begin strong branching
         double initialTimeSB = cpuTime();
         if ( strongBranchingSize > 1 ){
            cout << "Doing Strong Branching: " << strongBranchingSize << endl;
         }
         for ( int t = 0; t < strongBranchingSize && t < current->notCovered.size(); t++ ){
            k = ( *stBrchit );
            stBrchit++;				
            sbAuxStruct * tempSbStr = new sbAuxStruct;		
            tempSbStr->index = k;
            tempSbStr->sum = 0;

            for ( int i = 0; i < current->pts.size(); i++ ){

               node * child = new node;					
               child->s_lb = 0;
               child->pts = bnbPtr->insert( current->pts, i, k );

               SolveSocpCplex *solveSocpPtr2 = new SolveSocpCplex( dataptr, child->pts );

               initialSocpCompTime = cpuTime();
               solveSocpPtr2->solveSOCP( child->pts );					
               finalSocpCompTime = cpuTime();
               totalSocpCompTime += ( finalSocpCompTime - initialSocpCompTime );

               child->lb = solveSocpPtr2->getF_value();
               solveSocpPtr2->finishSOCP();
               somaTeste += solveSocpPtr2->violation;

               itCount++;

               for ( int i2 = 0; i2 < child->pts.size(); i2++ ){
                  tempX.push_back( solveSocpPtr2-> getSolutionX( i2 ) );
                  tempY.push_back( solveSocpPtr2-> getSolutionY( i2 ) );
                  tempZ.push_back( solveSocpPtr2-> getSolutionZ( i2 ) );
               }

               child->solXYZ.push_back( tempX );
               child->solXYZ.push_back( tempY );
               child->solXYZ.push_back( tempZ );

               //check feasibility
               bool feasibilityTest;
               feasibilityTest = bnbPtr->check_feasibility_Q( child->pts, child->solXYZ );

               if( feasibilityTest ){
                  child->feasible = 1;
               }
               else{
                  child->feasible = 0;
               }
               tempX.clear();
               tempY.clear();
               tempZ.clear();

               child->notCovered.clear();					

               //check not covered clients
               for ( int j = 0; j < bnbPtr->notCoveredBalls.size(); j++ ){
                  child->notCovered.push_back( bnbPtr->notCoveredBalls[ j ] );
               }
               //								
               if ( child->lb < best_ub && child->lb >= ub ){
                  if ( feasibilityTest == true ) best_ub = child->lb;
               }

               //fill the vector of children
               if( child->lb < ub ){
                  tempSbStr->sum += child->lb;
               }
               else{
                  tempSbStr->sum += ub;
               }

               tempSbStr->candidates.push_back( child );
               delete solveSocpPtr2;
            }
            vectorOfChildren.push_back( tempSbStr );
         }

         double finalTimeSB = cpuTime();
         sbComputationTime += finalTimeSB - initialTimeSB;

         //take the best sum
         int bestSum = 0;
         int pos = 0;
         for ( int s0 = 0; s0 < vectorOfChildren.size(); s0++ ){
            if( vectorOfChildren[ s0 ]->sum > bestSum ){
               bestSum = vectorOfChildren[ s0 ]->sum;
               pos = s0;
            }
         }

         //get upper bounds
         for ( int s0 = 0; s0 < vectorOfChildren.size(); s0++ ){
            for ( int s1 = 0; s1 < vectorOfChildren[ s0 ]->candidates.size(); s1++ ){
               if ( vectorOfChildren[ s0 ]->candidates[ s1 ]->lb < ub && vectorOfChildren[ s0 ]->candidates[ s1 ]->feasible == 1 ) ub = vectorOfChildren[ s0 ]->candidates[ s1 ]->lb;
               if ( vectorOfChildren[ s0 ]->candidates[ s1 ]->lb < best && vectorOfChildren[ s0 ]->candidates[ s1 ]->feasible == 1 ){
                  best = vectorOfChildren[ s0 ]->candidates[ s1 ]->lb;
                  solucao = vectorOfChildren[ s0 ]->candidates[ s1 ]->pts;
                  solucaoXYZ = vectorOfChildren[ s0 ]->candidates[ s1 ]->solXYZ;
               }
            }
         }

         for ( int s0 = 0; s0 < vectorOfChildren.size(); s0++ ){
            if( s0 != pos ){
               for ( int s00 = 0; s00 < vectorOfChildren[ s0 ]->candidates.size(); s00++ ){
                  delete vectorOfChildren[ s0 ]->candidates[ s00 ];
               }
               delete vectorOfChildren[ s0 ];
            }
         }

         //select candidates for branching
         for ( int s = 0; s < vectorOfChildren[ pos ]->candidates.size(); s++ ){				
            count_SOCP_solved ++;
            if( vectorOfChildren[ pos ]->candidates[ s ]->lb <= ub ){
               if ( vectorOfChildren[ pos ]->candidates[ s ]->feasible == 1 ){
                  ub = vectorOfChildren[ pos ]->candidates[ s ]->lb;
                  quantity += mpz_class( sizeOfTree )/levels[ vectorOfChildren[ pos ]->candidates[ s ]->pts.size() - 3 ] - 1;
                  if ( vectorOfChildren[ pos ]->candidates[ s ]->lb <  best ){
                     best = vectorOfChildren[ pos ]->candidates[ s ]->lb;
                     solucao = vectorOfChildren[ pos ]->candidates[ s ]->pts;
                     solucaoXYZ = vectorOfChildren[ pos ]->candidates[ s ]->solXYZ;
                  }

                  for ( itOpen = open.begin(); itOpen != open.end(); itOpen++ ){
                     double lowerB = ( *itOpen )->lb;
                     if ( lowerB > ub ){
                        quantity += mpz_class( sizeOfTree )/levels[ ( *itOpen )->pts.size() - 3 ] - 1;
                        open.erase( itOpen );								
                        itOpen = open.begin();
                     }								
                  }
                  delete vectorOfChildren[ pos ]->candidates[ s ];
                  vectorOfChildren[ pos ]->candidates[ s ] = NULL;
               }
               else{
                  open.push_back( vectorOfChildren[ pos ]->candidates[ s ] );
                  //print log
                  if ( print_counter == floor( sizeInst/2 ) ){
                     bnbPtr->printLog(temp, sizeOfTree, vectorOfChildren[ pos ]->candidates[ s ], open, count_SOCP_solved, quantity, best_ub, best_lb, vectorOfChildren[ pos ]->candidates[ s ]->notCovered.size(), &printHeader );
                     print_counter = 0;
                  }
                  print_counter++;
               }
            }
            else{
               quantity += mpz_class( sizeOfTree )/levels[ vectorOfChildren[ pos ]->candidates[ s ]->pts.size() - 3 ] - 1;
               delete vectorOfChildren[ pos ]->candidates[ s ];
               vectorOfChildren[ pos ]->candidates[ s ] = NULL;
            }
         }

         //COMPUTE LOWER BOUNDS
         //---------------------------------------------------------------------
         if ( current->s_lb == 1 && !open.empty() ){
            bnbPtr->computeLowerBounds( &open, current, &best_lb );
         }
         //---------------------------------------------------------------------

         vectorOfChildren.clear();
      }
      temp = mpz_class( sizeOfTree );
      temp = ( quantity/temp )*100;

      if ( cpuTime() - initialTotalTimeBnB > timeLimit ) optimalFound = 1;

   }

   mpf_class prctComp = mpz_class( sizeOfTree );
   prctComp = count_SOCP_solved/prctComp*100;

   double finalTotalTimeBnB = cpuTime();
   double computationTime = finalTotalTimeBnB - initialTotalTimeBnB;

   double gap_root = ( ( ub - rootLB )/ ub )*100;
   double gap_real = ( ( bestKnown - best_lb )/ bestKnown )*100;
   double gap_lb_bnb = ( ( ub - best_lb )/ ub )*100;
   //Finish Branch and Bound

   cout << endl;
   cout << "### Final Log ###" << endl << endl;
   if( optimalFound == 0 ){
      cout << "OPTIMAL SOLUTION FOUND" << endl;
      cout << "Function objective value: " << setiosflags (ios::fixed | ios::showpoint) << setprecision( 15 ) << best << endl;
      cout << "Soma das inviabilidades: " << somaTeste << endl;
      cout << "Sequência: ";
      for ( int i = 0; i < solucao.size(); i++ ){
         cout << solucao[ i ] << " ";
      }
      cout << endl;
      cout << "Solução: \n";
      for ( int j = 0; j < solucaoXYZ[ 0 ].size(); j++ ){
      }
      cout << endl;		

      printDataToFile ( dataptr, option, overlap, sizeInst, bestKnown, best, best_lb, gap_root, count_SOCP_solved, itCount, computationTime, sbComputationTime, totalSocpCompTime );

      printDataToMatlab( dataptr, sizeInst, overlap, best, solucao, solucaoXYZ );
   }
   else{

      cout << "NO OPTIMAL SOLUTION FOUND" << endl;
      cout << setiosflags ( ios::showpoint ) << setprecision( 30 );
      cout << "Lower Bound: " << best_lb << endl;
      cout << "Upper Bound: " << ub << endl;
      cout << "GAP(LB): " << ( (ub - best_lb)/ub )*100 << "% " << endl;

      printDataToFile ( dataptr, option, overlap, sizeInst, bestKnown, ub, best_lb, gap_real, gap_lb_bnb, gap_root, count_SOCP_solved, itCount, computationTime, sbComputationTime, totalSocpCompTime );
   }

   cout << "Número de Nós Resolvidos: " << count_SOCP_solved << endl;
   cout << "Porcentagem da árvore podada: " << temp << "%" << endl;
   cout << "Porcentagem da árvore computada: " << prctComp <<  "%" << endl;
   cout << "Total: " << temp + prctComp << endl;
   cout << "Tempo total: " << computationTime << endl;
   cout << "Tempo S.B: " << sbComputationTime << endl;
   cout << "Tempo SOCP: "<< totalSocpCompTime << endl;

   cout << endl << "#################" << endl;	
   mpz_clear ( sizeOfTree );
   delete root;	
   delete bnbPtr;
   delete solveSocpPtr;
   delete dataptr;

   return 0;
}



















