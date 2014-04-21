#ifndef PROTEIN_H
#define PROTEIN_H

//#include "calculate_mean_sd.cpp"
#include "linked_list.h"
#include "aa.h"
#include "chain.h"
#include "model.h"

class Protein:public Model
{
public:
  Protein();
  Protein(char *Input_File);

  char* GetPrName();
  int   GetPrAaNum();
  int   GetPrAtomNum();
  node<Atom>* GetPrCenter();
  node<Residue>* GetPrHeadAa();
  node<Residue>* GetPrTailAa();
  node<Model> *GetPrHeadModel();
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
  void SetPrAaList(node<Residue> *source_aa);
  void SetPrSurfaceAaList();
  void SetPrSurfaceConvexAaList();
  void SetPrSurfaceConcaveAaList();
  void SetPrBurialAaList();
  void SetPrConvexAaList();
  void SetPrConcaveAaList();
  void SetPrBindingSiteAaList();
  void SetPrSs(char *ss_file);
  void SetNcOneAa(char *p_name);
  void SetDifferentAminoAcidTypeInfoForThisProtein();
  void SetPrModelAveOmegaForEachResidue(); 
  
  void CalculateOmegaAngleForAaOfThisProtein();
  void CalculateDistanceBetweenCaOfAaForPr();
  void CalculateDistanceBetweenAtomForPr();
  void WriteThisProteinSummaryToFile();
  void CalculateAndOutputDifferenceInAverageOmegaBetweenBuAndSu();
  void OutputPrSurfaceAalist();
  void OutputPrSurfaceConvexAaList();
  void OutputPrSurfaceConcaveAaList();
  void OutputPrBurialAaList();
  void OutputPrBackBoneAtom();
  void OutputPrAaOmegaBasedOnEachChainCenter(char *input);
  void OutputPrSeSsPo(char *input_file);
 // void OutputPrNameAaPoSsSc();



  //void operator =(Protein*& source);
 

private:

   char *p_name;
   char *p_ss_name;
   char *p_name_chain;
   char *one_letter_aa;
   int num_of_aa;
   int num_of_atom;
   int num_of_chain;
   node<Atom> *head_atom_ptr;
   node<Atom> *tail_atom_ptr;

   node<Residue> *head_aa_ptr;
   node<Residue> *tail_aa_ptr;

   node<Chain> *head_chain_ptr;
   node<Chain> *tail_chain_ptr;

   node<Model> *head_model_ptr;

   node<Aa_Summary>  *aa_sum,*aa_sum_based_on_model;

   node<char> *head_se_ptr;
   node<char> *head_ss_ptr;
   node<int>  *head_pos_ptr;

   node<Atom> *Center_Of_Structure;
   node<Residue> *surface_head_aa_ptr;
   node<Residue> *surface_Convex_head_aa_ptr;
   node<Residue> *surface_Concave_head_aa_ptr;
   node<Residue> *burial_head_aa_ptr;
   node<Residue> *convex_head_aa_ptr;
   node<Residue> *concave_head_aa_ptr;

   node<Atom> *Binding_Site_Residue_head_ptr;

friend ostream& operator<<(ostream& os,Protein& protein);
};

#endif
