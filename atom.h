//#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <math.h>


using namespace std;

#ifndef ATOM_H
#define ATOM_H

#include "Vector3D.h"
#include "linked_list.h"


class Atom
{
 public:
    Atom();
   
    char *GetAtomType(); 
    int GetAtom_Index();
    char *GetAtomName();
    double GetAtom_x();
    double GetAtom_y();
    double GetAtom_z();
    int GetResidueNum();
    double GetAtom_radius();
    char *GetResidueName();
    char *GetChainType();
    double GetAtom_occ();
    double GetAtom_B_factor();
    double GetAtom_asa();
    double GetAtomAs();
    double GetAtomMs();
    double GetAtomCurvature();
    double GetAtomDistance();
    node<Atom> *GetAtomNeighborList();

    
    void SetAtomType(char *temp);
    void SetAtomIndex(int input_num);
    void SetAtomName(char *input_name);
    void SetResidueName(char *input_name);
    void SetChainType(char *input_name);
    void SetResidueNum(int input_num);    
    void SetAtom_asa(double input_asa);
    void SetAtomOcc(double Input);
    void SetAtom_x(double a);
    void SetAtom_y(double b);
    void SetAtom_z(double c);
    void SetAtomRadius(double e);
    void SetAtomAs(double e);
    void SetAtomMs(double e);
    void SetAtomCurvature(double e);
    void SetAtomDistance(double d);
    void SetAtomNeighborAtomList(Atom& atom);
    void SetSortedNeighborAtomList();
    void SetAtomBfactor(double Input);
    void SetProteinName(char *temp);
    void SetLineNumber(int input);

    //void operator =(Atom*& atom);

 
    
 private:

   char   *p_name;
   int     LineNum;
   char   *AtomType;
   int    AtomIndex;
   char   *AtomName;
   char   *Re;
   char   *chain_type;
   int    res_num;
   double x;
   double y;
   double z;
   double occ;
   double Bfactor;
   double asa;
   double radius;
   double As;
   double Ms;
   double Curvature;
   double distance;
   node<Atom> *Nb_head_atom_ptr;
   node<Atom> *Nb_tail_atom_ptr;
   

 friend ostream& operator<<(ostream& os,Atom& atom);
   
};
#endif
