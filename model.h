#ifndef MODEL_H
#define MODEL_H

//#include "calculate_mean_sd.cpp"
#include "linked_list.h"
#include "aa.h"
#include "aa_summary.h"

class Model
{
public:
  Model();
  Model(char *Input_File);

  char* GetModelPrname();
  char* GetModelName();
  int   GetModelNum();
  int   GetModelAaNum();
  int   GetModelAtomNum();
  node<Atom>* GetModelCenter();
  node<Residue>* GetModelHeadAa();
  node<Residue>* GetModelTailAa();
  node<Atom>* GetModelHeadAtom();
  node<Atom>* GetModelTailAtom();
 // if(cursor_chain_ptr->next!=NULL)
  node<Aa_Summary>* GetHeadAa_Summary();


  node<Atom>* GetModelSurfaceHeadCalphaAtom();
  double GetModelAaOmegaAverage(vector<double> input );
  double GetModelAaOmegaSd_Of_Ave(vector<double> input);
  int GetNumberofModelForPr();
  void GetNumberofChainForModel();
  int GetNumberofAaForModel();
  int GetNumberofAtomForModel();

  void SetPrName(char *input);
  void SetModelName(char *input);
  void SetModelNum(int num);
  void SetModelAaNum(int num);
  void SetModelAtomNum(int num);
  //void SetModelAaList(Residue &residue)
  void SetModelCenter();
  void SetModelAtomList(node<Atom> *source_atom);

  void SetModelAaList(node<Residue> *source_aa);
  void SetModelSurfaceConcaveAaList();

  void SetModelSurfaceAaList();
  void SetModelSurfaceConvexAaList();
  //void SetModelSurfaceConcaveAaList();
  void SetModelBurialAaList();
  void SetModelConvexAaList();
  void SetModelConcaveAaList();
  void SetModelBindingSiteAaList();
  void SetDifferentAminoAcidTypeInfoForThisModel();
  
  
  void CalculateOmegaAngleForAaOfThisModel();
  void CalculateDistanceBetweenCaOfAaForModel();
  void CalculateDistanceBetweenAtomForModel();
  void WriteThisModelSummaryToFile();
  void CalculateAndOutputDifferenceInAverageOmegaBetweenBuAndSu();
  void OutputModelSurfaceAalist();
  void OutputModelSurfaceConvexAaList();
  void OutputModelSurfaceConcaveAaList();
  void OutputModelBurialAaList();
  void OutputModelBackBoneAtom();



  //void operator =(Chain*& source);
 

private:
   
   char *model_name;
   int   model_no;
   char *p_name;
   int num_of_aa;
   int num_of_atom;
   node<Aa_Summary> *aa_sum;
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

friend ostream& operator<<(ostream& os,Model& model);
};

#endif
