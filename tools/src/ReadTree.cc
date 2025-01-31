#include "ReadTree.hh"

#include <TBranch.h>
#include <TFriendElement.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>

// ClassImp(TreeReader)

TreeReader::TreeReader( std::string treeName ) {
  fChain          = new TChain( treeName.c_str() );
  nGets           = 0;
  nSwaps          = 0;
  printFrequency  = 100000;
  nEntries        = -999;
  continueSorting = true;
}

TreeReader::TreeReader() {
  fChain          = new TChain( "" );
  nGets           = 0;
  nSwaps          = 0;
  printFrequency  = 100000;
  nEntries        = -999;
  continueSorting = true;
}

TreeReader::~TreeReader() {
  // Some statistics for debugging
  std::cout << "Total calls of GetValue function: " << nGets << std::endl;
  std::cout << "Variables accessed with more than 1% frequency and their "
               "positions:\n";
  for ( unsigned int i = 0; i < varList.size(); ++i ) {
    if ( 100. * varList[i]->nGets / nGets >= 1 ) {
      std::cout << "   " << varList[i]->name << " "
                << 100 * varList[i]->nGets / nGets << "  " << i << std::endl;
    }
  }
  std::cout << "Number of calls to sort function: " << nSwaps << std::endl;

  delete fChain;
}

int TreeReader::GetEntry( Long64_t ientry ) {
  if ( ientry % printFrequency == 0 ) {
    std::cout << "Reading entry: " << ientry << "/" << nEntries << std::endl;
  }
  return fChain->GetEntry( ientry );
}

void TreeReader::Initialize() {
  Initialize( fChain );
}

/*
 * In this place we need to get list of all branches and loop over them.
 * For each branch, determine type, fill in variable which gets pushed to
 * varList and tree branch address will be initialized to appropriate
 * variable.
 */
void TreeReader::Initialize( TTree* chain, bool topLevel ) {
  std::cout << "Initialise chain: " << chain->GetName() << std::endl;

  TList* friends = chain->GetListOfFriends();
  if ( friends != nullptr && friends->GetEntries() > 0 ) {
    for ( int i = 0; i < friends->GetEntries(); ++i ) {
      Initialize(
          dynamic_cast< TFriendElement* >( friends->At( i ) )->GetTree(),
          false );
    }
  }

  chain->GetEntry(0);
  TObjArray* branches = chain->GetListOfBranches();
  int nBranches       = branches->GetEntries();

  for ( int i = 0; i < nBranches; ++i ) {
    TBranch* branch = dynamic_cast< TBranch* >( branches->At( i ) );
    TLeaf* leaf     = branch->GetLeaf( branch->GetName() );
    if ( leaf == 0 )  // leaf name is different from branch name
    {
      TObjArray* leafs = branch->GetListOfLeaves();
      leaf             = dynamic_cast< TLeaf* >( leafs->At( 0 ) );
    }
    variable* tmpVar;
    varList.push_back( new variable() );
    tmpVar = ( varList[varList.size() - 1] );

    if ( topLevel ) {
      tmpVar->name = leaf->GetName();
    } else {
      tmpVar->name = std::string( chain->GetName() ) + "_" +
                     std::string( leaf->GetName() );
    }
    tmpVar->bname         = branch->GetName();
    tmpVar->title         = leaf->GetTitle();
    tmpVar->copyToNewTree = true;
    // This does not work before reading event and even then I would have
    // to make assumption.
    // tmpVar->arraySize=leaf->GetLen();
    // Find out whether we have array by inspecting leaf title
    if ( tmpVar->title.find( "[" ) != std::string::npos ) {
      tmpVar->arraySize =
          2000;  // Dirty hack, tree knows actual size, I do not have to
    } else {
      tmpVar->arraySize = 0;
    }
    tmpVar->nGets = 0;

    if ( strcmp( leaf->GetTypeName(), "Float_t" ) == 0 ) {
      tmpVar->type = FLOAT;
    } else if ( strcmp( leaf->GetTypeName(), "Int_t" ) == 0 ) {
      tmpVar->type = INT;
    } else if ( strcmp( leaf->GetTypeName(), "Bool_t" ) == 0 ) {
      tmpVar->type = BOOL;
    } else if ( strcmp( leaf->GetTypeName(), "Double_t" ) == 0 ) {
      tmpVar->type = DOUBLE;
    }

    variable* lastAdded = tmpVar;

    switch ( lastAdded->type ) {
      case FLOAT:
        if ( lastAdded->arraySize > 1 )
          chain->SetBranchAddress( lastAdded->bname.c_str(),
                                   (void*)lastAdded->float_array );
        else
          chain->SetBranchAddress( lastAdded->bname.c_str(),
                                   (void*)&lastAdded->float_value );
        break;

      case DOUBLE:
        if ( lastAdded->arraySize > 1 )
          chain->SetBranchAddress( lastAdded->bname.c_str(),
                                   (void*)lastAdded->double_array );
        else
          chain->SetBranchAddress( lastAdded->bname.c_str(),
                                   (void*)&lastAdded->double_value );
        break;

      case INT:
        chain->SetBranchAddress( lastAdded->bname.c_str(),
                                 (void*)&lastAdded->int_value );
        break;

      case BOOL:
        chain->SetBranchAddress( lastAdded->bname.c_str(),
                                 (void*)&lastAdded->bool_value );
        break;
      default:
        break;
    }
    //    std::cout<<"Branch: "<<lastAdded->name<<"; "<<lastAdded->type<<";
    //    "<<lastAdded->arraySize<<std::endl;
  }

  std::cout << "Set up " << varList.size() << " branches\n";
  if ( topLevel ) {
    nEntries = chain->GetEntries();
  }
}

// This needs update
void TreeReader::BranchNewTree( TTree* tree ) {
  for ( std::vector< variable* >::iterator it = varList.begin();
        it != varList.end(); ++it ) {
    variable* var = *it;
    if ( var->copyToNewTree == false ) {
      continue;
    }
    switch ( var->type ) {
      case FLOAT:
        if ( var->arraySize > 1 )
          tree->Branch( var->name.c_str(), (void*)( var->float_array ),
                        ( var->title + "/F" ).c_str() );
        else
          tree->Branch( var->name.c_str(), (void*)&( var->float_value ),
                        ( var->title + "/F" ).c_str() );
        break;
      case DOUBLE:
        if ( var->arraySize > 1 )
          tree->Branch( var->name.c_str(), (void*)( var->double_array ),
                        ( var->title + "/D" ).c_str() );
        else
          tree->Branch( var->name.c_str(), (void*)&( var->double_value ),
                        ( var->title + "/D" ).c_str() );
        break;
      case INT:
        tree->Branch( var->name.c_str(), (void*)&( var->int_value ),
                      ( var->title + "/I" ).c_str() );
        break;
      case BOOL:
        tree->Branch( var->name.c_str(), (void*)&( var->bool_value ),
                      ( var->title + "/O" ).c_str() );
        break;
      default:
        break;
    }
  }

  return;
}

double TreeReader::GetValue( std::string name, int iValue ) {
  ++nGets;

  //  std::cout<<"DEBUG GetValue: "<<name<<std::endl;

  if ( continueSorting && nGets % 1000 == 0 && varList.size() > 15 ) {
    continueSorting = partialSort();
  }

  variable* myVar =
      ( *( std::find_if( varList.begin(), varList.end(), varEq( name ) ) ) );
  if ( !( ( myVar ) != ( *( varList.end() ) ) ) ) {
    std::cout << "Trying to access non-existing variable with name " << name
              << std::endl;
    std::exit( 1 );
  }
  ++( myVar->nGets );
  //  std::find_if(varList.begin(), varList.end(), varEq(name));

  //  std::cout<<"Accessing branch "<<myVar->name<<" of type "<<myVar->type
  //            <<" with double value "<<myVar->double_value<<std::endl;

  switch ( myVar->type ) {
    case FLOAT:
      if ( myVar->arraySize > 1 )
        return myVar->float_array[iValue];
      else
        return myVar->float_value;
      break;

    case DOUBLE:
      if ( myVar->arraySize > 1 )
        return myVar->double_array[iValue];
      else
        return myVar->double_value;
      break;

    case INT:
      return myVar->int_value;
      break;

    case BOOL:
      return myVar->bool_value;
      break;
    default:
      return -999;
  }

  //  float
  //  value=fChain->GetBranch(name.c_str())->GetLeaf(name.c_str())->GetValue(iValue);
  //  return value;

  return -999;
}

void TreeReader::setAllOffNewTree() {
  for ( std::vector< variable* >::iterator it = varList.begin();
        it != varList.end(); ++it ) {
    variable* var      = *it;
    var->copyToNewTree = false;
  }
}

void TreeReader::setUseInNewTree( std::string name, bool use ) {
  variable* myVar =
      ( *( std::find_if( varList.begin(), varList.end(), varEq( name ) ) ) );
  if ( ( ( myVar ) != ( *( varList.end() ) ) ) ) {
    myVar->copyToNewTree = use;
    return;
  }
}

/*
 * In future I might keep statistics of which variables are accessed often
 * and once in a while swap those to early positions in list to get better
 * than n*ln(n) of find. Ideally if I access only 10 variables, at some
 * point they should occupy first 10 positions while others (order of
 * hundred) ones should be behind. This way find should be more efficient.
 * Have to see whether I can implement and tune idea to actually obtain
 * gain as it will need some statistics collection, checks whether we
 * should still do it and at least some sorting of list.
 *
 * Using chain GetBranch and GetLeaf runs at about twice the speed of my
 * implementation without doing anything special in searching for branch.
 * Should definitely try to be smart.
 */

bool TreeReader::partialSort() {
  bool didSwap = false;
  ++nSwaps;

  // Do not care about first 10 variables
  for ( unsigned int i = 10; i < varList.size(); ++i ) {
    if ( 100.0 * varList[i]->nGets / nGets >=
         1 )  // Care only about cases which are used more often
    {
      for ( unsigned int j = 0; j < i;
            ++j )  // Makes sense to swap only to earlier place and if variable
                   // itself is not used often
      {
        if ( varList[i]->nGets > 2 * varList[j]->nGets &&
             100.0 * varList[j]->nGets / nGets < 1.0 ) {
          variable* tmp = varList[i];
          varList[i]    = varList[j];
          varList[j]    = tmp;
          didSwap       = true;
        }
      }
    }
  }

  return didSwap;
}

void TreeReader::setValue( std::string name, double input, int iValue ) {
  variable* myVar =
      ( *( std::find_if( varList.begin(), varList.end(), varEq( name ) ) ) );
  if ( !( ( myVar ) != ( *( varList.end() ) ) ) ) {
    std::cout << "Trying to access non-existing variable with name " << name
              << std::endl;
    std::exit( 1 );
  }
  switch ( myVar->type ) {
    case FLOAT:
      if ( myVar->arraySize > 1 ) {
        myVar->float_array[iValue] = input;
      } else {
        myVar->float_value = input;
      }
      break;

    case DOUBLE:
      if ( myVar->arraySize > 1 )
        myVar->double_array[iValue] = input;
      else
        myVar->double_value = input;
      break;

    case INT:
      myVar->int_value = input;
      break;

    case BOOL:
      myVar->bool_value = input;
      break;
    default:
      break;
  }

  return;
}

