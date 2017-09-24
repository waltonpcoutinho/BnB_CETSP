#ifndef DATA_H
#define DATA_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<list>
#include<cstdlib>
#include<cstring>
#include<climits>
#include<utility>
#include<algorithm>    // std::sort
#include<cmath>

using namespace std;

class Data{

   public:

      Data( char *, char *, double, int, char** ); //constructor
      ~Data();

      //SET Functions
      void setAllData( ifstream & );
      void setAllData2( ifstream & );
      void deleteZCoordinates();

      //GET Functions
      int getSizeInst();
      double getCoordx( int );
      double getCoordy( int );
      double getCoordz( int );
      double getRadius( int );
      double getDemand( int );
      int getBranchingRuleList( int );
      int getDepotFarthest( int );
      double getUbFromFile();
      double getOverlapRatio();

      //PRINT Functions
      void printInstance();	
      string fileName;

   private:

      int sizeInst;
      double overlapRatio;
      double overlap_for_this_inst;
      double ub;

      vector< double > coordx;
      vector< double > coordy;
      vector< double > coordz;
      vector< double > radius;
      vector< double > demand;

      vector< vector < int > > intersecMatrix;
      vector< int > countIntersecs;
      void setIntersecMatrix();

      void setSizeInst( char * );

      list < int > branchingRuleList;
      list< int >::iterator iterBranchingRuleList;

      void setBranchingRuleList();

      void setDistanceMatrix();

      void setDepotFarthest();

      void setUb( string, string, char *, double );

      vector< vector< double > > distanceMatrix;
      vector< pair< int, double > > farthest;

};

struct pStruct {
   int index;
   int nOfIntersec;
   float radius;
   vector < int > intersections;
   int size;
};

struct coordinates {
   double x;
   double y;
   double z;
};


#endif
