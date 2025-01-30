#ifndef READTREE
#define READTREE

/*
 * Library to read tree from LHCb TupleTools in automatic way which would
 * allow to get variables by names.
 */

#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include "TChain.h"
#include "TTree.h"

typedef enum { INT, FLOAT, DOUBLE, BOOL } Types;

class variable {
 public:
  std::string name;
  std::string bname;
  std::string title;
  float float_value;
  double double_value;
  int int_value;
  bool bool_value;
  float float_array[2000];
  double double_array[2000];
  int arraySize;
  Types type;
  int nGets;
  bool copyToNewTree;
};

class varEq : public std::unary_function< variable*, bool > {
  std::string s;

 public:
  explicit varEq( const std::string& ss ) : s( ss ) {}
  bool operator()( const variable* c ) const {
    return s.compare( c->name ) ? false : true;
  }
};

class TreeReader {
 public:
  TreeReader( std::string treeName );
  TreeReader();
  ~TreeReader();

  void AddFile( std::string filename ) {
    if ( fChain )
      fChain->Add( filename.c_str() );
    std::cout << "Number of trees in chain: " << fChain->GetNtrees()
              << std::endl;
  }
  Long64_t GetEntries() { return fChain->GetEntries(); }
  int GetEntry( Long64_t ientry );

  void Initialize();

  void BranchNewTree( TTree* tree );

  TChain* GetChain() { return fChain; }

  double GetValue( std::string name, int iValue = 0 );
  void setValue( std::string name, double input, int iValue = 0 );

  void setAllOffNewTree();
  void setUseInNewTree( std::string name, bool use = true );
  void setPrintFrequency( long value ) { printFrequency = value; }
  void addFriend( std::string name, std::string file = "" ) {
    fChain->AddFriend( name.c_str(), file.c_str() );
  };

 private:
  void Initialize( TTree* chain, bool topLevel = true );

  TChain* fChain;
  std::vector< variable* > varList;

  int nGets;
  int nSwaps;
  bool continueSorting;
  bool partialSort();
  long printFrequency;
  long nEntries;

  //  ClassDef(TreeReader,2)
};

#endif

