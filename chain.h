#ifndef CHAIN_H
#define CHAIN_H

//#include "calculate_mean_sd.cpp"
#include "linked_list.h"
#include "aa.h"

class Chain
{
public:
  Chain();
  Chain(char *Input_File);

  char* GetChainType();
  char* GetPrName();
  int   GetPrAaNum();
  int   GetPrAtomNum();
  node<Atom>* GetChainCenter();
  node<Residue>* GetChainHeadAa();
  node<Residue>* GetChainTailAa();
  node<Atom>* GetChainHeadAtom();
  node<Atom>* GetChainTailAtom();


  node<Atom>* GetPrSurfaceHeadCalphaAtom();
  double GetPrAaOmegaAverage(vector<double> input );
  double GetPrAaOmegaSd_Of_Ave(vector<double> input);
  void GetNumberofChainForPr();
  int GetNumberofAaForPr();
  int GetNumberofAtomForPr();


  void SetPrName(char *input);
  void SetPrAaNum(int num);
  void SetPrAtomNum(int num);
  void SetPrCenter();
  void SetPrAtomList(node<Atom> *source_atom);

  void SetChainAaList(Residue& source_aa);

  void SetPrSurfaceAaList();
  void SetPrSurfaceConvexAaList();
  void SetPrSurfaceConcaveAaList();
  void SetPrBurialAaList();
  void SetPrConvexAaList();
  void SetPrConcaveAaList();
  void SetPrBindingSiteAaList();
  
  
  void CalculateOmegaAngleForAaOfThisChain();
  void CalculateDistanceBetweenCaOfAaForPr();
  void CalculateDistanceBetweenAtomForPr();
  void WriteThisChainSummaryToFile();
  void CalculateAndOutputDifferenceInAverageOmegaBetweenBuAndSu();
  void OutputPrSurfaceAalist();
  void OutputPrSurfaceConvexAaList();
  void OutputPrSurfaceConcaveAaList();
  void OutputPrBurialAaList();
  void OutputPrBackBoneAtom();



  //void operator =(Chain*& source);
 

private:

   char *p_name;
   int num_of_aa;
   int num_of_atom;
   node<Atom> *Center_Of_Structure;
   node<Atom> *head_atom_ptr;
   node<Residue> *tail_atom_ptr;
   
   node<Residue> *head_aa_ptr;
   node<Residue> *tail_aa_ptr;

   node<Residue> *surface_head_aa_ptr;
   node<Residue> *surface_Convex_head_aa_ptr;
   node<Residue> *surface_Concave_head_aa_ptr;
   node<Residue> *burial_head_aa_ptr;
   node<Residue> *convex_head_aa_ptr;
   node<Residue> *concave_head_aa_ptr;
   node<Atom> *Binding_Site_Residue_head_ptr;

friend ostream& operator<<(ostream& os,Chain& chain);
};

#endif
