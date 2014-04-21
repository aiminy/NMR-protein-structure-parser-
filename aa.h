#ifndef AA_H
#define AA_H

#include "atom.h"
#include "linked_list.h"


class Residue
{
 public:
    Residue();
    char* GetAaName();    
    int   GetAaNum();
    char* GetAaChain();
    int   GetAaAtomNum();
    int   GetAaSideChainAtomNum();
    void  CalculateAndSetOmegaAngleForAa(Vector3D A,Vector3D B);
    void  CalculateAndSetOmegaMcAngleForAa(Vector3D A,Vector3D B);
    double GetAaMs();
    double GetAaAs();
    double GetAaAverageCurvature();
    double GetOmega();
    double GetOmegaMc();
    double GetAaDistance();
    char*  GetAaPrChain();
    char*  GetAaOneLetterCode();

    node<Residue> *GetAaNeighborList();
 
    node<Atom> *GetAaHeadAtom();
    node<Atom> *GetAaTailAtom();
    node<Atom> *GetAaCalphaAtom();
    node<Atom> *GetAaSideChainGeometricalCenter();
    node<Atom> *GetAaGeometricalCenter();
    node<Atom> *GetAaBetaAtom();
    void SetAaName(char *input);
    void SetAaNum(int num);
    void SetAaChain(char *input);
    void SetAaAtomList(Atom& atom);
    void SetAaDistance(double input);
    void SetAaNeighborAaList(Residue& residue);
    void SetSortedNeighborAaList();
    void SetAaPrChain(char *input);
    void SetAaOneLetterCode(char *input);
    //void operator =(Residue*& source);

    
 private:

   char *aa_name;
   int   aa_num;
   char *aa_chain;
   char *p_name_chain;
   char *one_letter_aa;
   double omega;
   double omegaMc;
   double d;
   node<Atom> *head_atom_ptr;
   node<Atom> *tail_atom_ptr;
   node<Residue> *Nb_head_aa_ptr;
   node<Residue> *Nb_tail_aa_ptr;

 friend ostream& operator<<(ostream& os,Residue& aa);
   
};
#endif
