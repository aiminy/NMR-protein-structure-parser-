#include "atom.h"

Atom::Atom()
{
  p_name=new char;
  AtomType=new char;
  AtomIndex=-1;
  AtomName=new char;
  x=0;
  y=0;
  z=0;
  radius=0;
  Ms=0;
  As=0;
  Curvature=0;
  occ=0;
  Bfactor=0;
  Nb_head_atom_ptr=NULL;
  Nb_tail_atom_ptr=NULL;
 
  
}

char* Atom::GetAtomType()
{
 return AtomType;
}
int Atom::GetAtom_Index()
{
  return AtomIndex;
}

char* Atom::GetAtomName()
{
return AtomName;
}
 
double Atom::GetAtom_x()
{
 return x;
}

double Atom::GetAtom_y()
{
 return y;
}

double Atom::GetAtom_z()
{
 return z;
}

int Atom::GetResidueNum()
{
 return res_num;
}

double Atom::GetAtom_radius()
{
return radius;
}

char* Atom::GetResidueName()
{
 return Re;
}

char* Atom::GetChainType()
{
return chain_type;
}



double Atom::GetAtom_occ()
{
 return occ;
}

double Atom::GetAtom_B_factor()
{
 return radius;
}

double Atom::GetAtom_asa()
{
return asa; 
}

double Atom::GetAtomAs()
{
return As;
}

double Atom::GetAtomMs()
{
return Ms;
}

double Atom::GetAtomCurvature()
{
return Curvature;
}

double Atom::GetAtomDistance()
{
return distance;
}

node<Atom>* Atom::GetAtomNeighborList()
{
return Nb_head_atom_ptr;
}

void Atom::SetAtomType(char *temp)
{
 AtomType=new char;
 strcpy(AtomType,temp);
}


void Atom::SetAtomIndex(int input_num)
{
AtomIndex=input_num;
}

void Atom::SetAtomName(char *input_atom_name)
{
 AtomName=new char;
 strcpy(AtomName,input_atom_name);
}

void Atom::SetResidueName(char *input_name)
{
 Re=new char;
 strcpy(Re,input_name);
}

void Atom::SetChainType(char *input_name)
{
  chain_type=new char;
  strcpy(chain_type,input_name);
}

void Atom::SetResidueNum(int input_num)
{
 res_num=input_num;
}



void Atom::SetAtom_x(double a)
{
 x=a;
}

void Atom::SetAtom_y(double b)
{
 y=b;
}

void Atom::SetAtom_z(double c)
{
 z=c;
}

void Atom::SetAtomOcc(double Input)
{
 occ=Input;
}

void Atom::SetAtomBfactor(double Input)
{
 Bfactor=Input;
}

void Atom::SetProteinName(char *temp)
{
  p_name=new char;
  strcpy(p_name,temp);

}
    
void Atom::SetLineNumber(int input)
{
  LineNum=input;
}

void Atom::SetAtom_asa(double input_asa)
{
 asa=input_asa;
}

void Atom::SetAtomRadius(double e)
{
radius=e;
}

void Atom::SetAtomAs(double e)
{
As=e;
}

void Atom::SetAtomMs(double e)
{
Ms=e;
}

void Atom::SetAtomCurvature(double e)
{
Curvature=e;
}
    
void Atom::SetAtomDistance(double d)
{
 distance=d;
}

void Atom::SetAtomNeighborAtomList(Atom& atom)
{
 insert_an_node(Nb_head_atom_ptr,Nb_tail_atom_ptr,atom);
}

void Atom::SetSortedNeighborAtomList()
{
 node<Atom> *temp,*cursor_atom_ptr1,*cursor_atom_ptr2;
 cursor_atom_ptr1=Nb_head_atom_ptr;
 cursor_atom_ptr2=Nb_head_atom_ptr;

 while(cursor_atom_ptr1!=NULL)
 {
 cursor_atom_ptr2=cursor_atom_ptr1;

 while(cursor_atom_ptr2!=NULL)
 {
 if(cursor_atom_ptr2->data.GetAtomDistance()<cursor_atom_ptr1->data.GetAtomDistance())
   {
    temp->data=cursor_atom_ptr1->data;
    cursor_atom_ptr1->data=cursor_atom_ptr2->data;
    cursor_atom_ptr2->data=temp->data;
   }
  cursor_atom_ptr2=cursor_atom_ptr2->next;
 }
  cursor_atom_ptr1=cursor_atom_ptr1->next;
 }


}

    
//void Atom::operator =(Atom*& atom)
//{
/*
 AtomIndex=atom->AtomIndex;

 AtomName=new char;
 strcpy(AtomName,atom->AtomName);

 Re=new char;
 strcpy(Re,atom->Re);

 chain_type=new char;
 strcpy(chain_type,atom->chain_type);

 res_num=atom->res_num;
 x=atom->x;
 y=atom->y;
 z=atom->z;
 occ=atom->occ;
 B_factor=atom->B_factor;
 asa=atom->asa;
 radius=atom->radius;
 As=atom->As;
 Ms=atom->Ms;
 Curvature=atom->Curvature;
*/
//}


ostream& operator<<(ostream& os,Atom& atom)
{
   //output format is same as the original pdb file
  os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<atom.AtomType
    <<setiosflags(ios::right)<<setw(5)<<atom.AtomIndex<<setw(2)<<" ";
  os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(4)<<atom.AtomName;
  os<<setw(4)<<atom.Re<<setw(1)<<atom.chain_type;
  os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<atom.res_num<<setiosflags(ios::fixed)
  <<setprecision(3)<<setw(12)<<atom.x
  <<setprecision(3)<<setw(8)<<atom.y
  <<setprecision(3)<<setw(8)<<atom.z
  <<setprecision(2)<<setw(6)<<atom.occ
  <<setprecision(2)<<setw(6)<<atom.Bfactor;
  //os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<" ";
  //os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(4)<<atom.p_name;
 // os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setiosflags(ios::fixed)
 //   <<setw(4)<<atom.LineNum;

 return os;
}
