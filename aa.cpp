#include "aa.h"

Residue::Residue()
{
 aa_name=new char;
 aa_num=0;
 aa_chain=new char;
 omega=999;
 omegaMc=999;
 head_atom_ptr=NULL;
 tail_atom_ptr=NULL;
 d=-1;
 Nb_head_aa_ptr=NULL;
 Nb_tail_aa_ptr=NULL;
 p_name_chain=new char;
 one_letter_aa=new char;
 
}

char* Residue::GetAaName()
{
 return head_atom_ptr->data.GetResidueName();
}

int Residue::GetAaNum()
{
return head_atom_ptr->data.GetResidueNum();
}


char* Residue::GetAaChain()
{
return head_atom_ptr->data.GetChainType();
}

int Residue::GetAaAtomNum()
{
  int length;
  length=linked_list_length(head_atom_ptr);
  return length;
}

int Residue::GetAaSideChainAtomNum()
{
  int length;
  node<Atom> *cursor_atom_ptr;
  cursor_atom_ptr=head_atom_ptr;
  while(strncmp(cursor_atom_ptr->data.GetAtomName(),"N  ",3)==0||
       strncmp(cursor_atom_ptr->data.GetAtomName(),"CA ",3)==0||
       strncmp(cursor_atom_ptr->data.GetAtomName(),"C  ",3)==0||
       strncmp(cursor_atom_ptr->data.GetAtomName(),"O  ",3)==0)
  {cursor_atom_ptr=cursor_atom_ptr->next;}
  length=linked_list_length(cursor_atom_ptr);
  return length;
  
}

void Residue::CalculateAndSetOmegaAngleForAa(Vector3D A,Vector3D B)
{
 
 double A_DOT_B=A.DotProduct(B);
 double A_length=A.CalculateLength();
 double B_length=B.CalculateLength();
 double AB=A_length*B_length;
 double angle=180*acos(A_DOT_B/AB)/3.14;
 
 omega=angle;

}

void Residue::CalculateAndSetOmegaMcAngleForAa(Vector3D A,Vector3D B)
 {
 double A_DOT_B=A.DotProduct(B);
 double A_length=A.CalculateLength();
 double B_length=B.CalculateLength();
 double AB=A_length*B_length;
 double angle=180*acos(A_DOT_B/AB)/3.14;
 
 omegaMc=angle;
}

double Residue::GetAaMs()
{
  node<Atom> *cursor_atom_ptr;
  double Sum_of_Ms;

  Sum_of_Ms=0;
  
  cursor_atom_ptr=head_atom_ptr;

//  cout<<cursor_atom_ptr->data.GetAtomMs()<<endl;

  while(cursor_atom_ptr!=NULL)
   {
    Sum_of_Ms=Sum_of_Ms+cursor_atom_ptr->data.GetAtomMs();
    cursor_atom_ptr=cursor_atom_ptr->next;
   }

  return Sum_of_Ms;

}

double Residue::GetAaAs()
{
  node<Atom> *cursor_atom_ptr;
  double Sum_of_As;

  Sum_of_As=0;
  
  cursor_atom_ptr=head_atom_ptr;
  while(cursor_atom_ptr!=NULL)
   {
    Sum_of_As=Sum_of_As+cursor_atom_ptr->data.GetAtomAs();
    cursor_atom_ptr=cursor_atom_ptr->next;
   }

  return Sum_of_As;

}

double Residue::GetAaAverageCurvature()
{
  node<Atom> *cursor_atom_ptr;
  double AaAverageCurvature;
  double Sum_of_area;
  double Sum_of_Curvature_Area;

  AaAverageCurvature=0;
  Sum_of_area=0;
  Sum_of_Curvature_Area=0;
  
  cursor_atom_ptr=head_atom_ptr;
  while(cursor_atom_ptr!=NULL)
   {
    Sum_of_Curvature_Area=Sum_of_Curvature_Area+cursor_atom_ptr->data.GetAtomMs()*cursor_atom_ptr->data.GetAtomCurvature();
    Sum_of_area=Sum_of_area+cursor_atom_ptr->data.GetAtomMs();
    cursor_atom_ptr=cursor_atom_ptr->next;
   }

   if(Sum_of_area!=0)
   {AaAverageCurvature=Sum_of_Curvature_Area/Sum_of_area;}
   

  return AaAverageCurvature;


}

double Residue::GetOmega()
{
 return omega;
}

double Residue::GetOmegaMc()
{
 return omegaMc;
}

node<Atom>* Residue::GetAaHeadAtom()
{
 return head_atom_ptr;
}

node<Atom>* Residue::GetAaTailAtom()
{
 return tail_atom_ptr;
}


node<Atom>* Residue::GetAaCalphaAtom()
{
 node<Atom> *CalphaAtom_ptr,*cursor_atom_ptr;
 CalphaAtom_ptr=new node<Atom>;

  int AtomId,AaId;
  char *AtomName,*AaName,*Chain;

  AtomName=new char;
  AaName=new char;
  Chain=new char;
  
  double n,av_x,av_y,av_z,av_occ,av_B,av_asa,av_radius;
   n=0;av_x=0;av_y=0;av_z=0;av_asa=0;av_radius=0;
   
  cursor_atom_ptr=head_atom_ptr;
  while(cursor_atom_ptr!=NULL)
   {
     if(strncmp(cursor_atom_ptr->data.GetAtomName(),"CA ",3)==0)
      {
       av_x=av_x+cursor_atom_ptr->data.GetAtom_x();
       av_y=av_y+cursor_atom_ptr->data.GetAtom_y();
       av_z=av_z+cursor_atom_ptr->data.GetAtom_z();
       av_asa=av_asa+cursor_atom_ptr->data.GetAtom_asa();
       av_radius=av_radius+cursor_atom_ptr->data.GetAtom_radius();
       AtomId=cursor_atom_ptr->data.GetResidueNum();
       AaId=cursor_atom_ptr->data.GetResidueNum();
       strcpy(AtomName,cursor_atom_ptr->data.GetAtomName());
       strcpy(AaName,cursor_atom_ptr->data.GetResidueName());
       strcpy(Chain,cursor_atom_ptr->data.GetChainType());
       n=n+1;
      }
      cursor_atom_ptr=cursor_atom_ptr->next;
   }

       if(n!=0)
       {
       av_x=av_x/n;
       av_y=av_y/n;
       av_z=av_z/n;
       av_asa=av_asa/n;
       av_radius=av_radius/n;

    CalphaAtom_ptr->data.SetAtomIndex(AtomId);
    CalphaAtom_ptr->data.SetAtomName(AtomName);
    CalphaAtom_ptr->data.SetResidueName(AaName);
    CalphaAtom_ptr->data.SetChainType(Chain);
    CalphaAtom_ptr->data.SetResidueNum(AaId);    
    CalphaAtom_ptr->data.SetAtom_x(av_x);
    CalphaAtom_ptr->data.SetAtom_y(av_y);
    CalphaAtom_ptr->data.SetAtom_z(av_z);
    CalphaAtom_ptr->data.SetAtom_asa(av_asa);
    CalphaAtom_ptr->data.SetAtomRadius(av_radius);
    }
   else
    {CalphaAtom_ptr=NULL;}

  return CalphaAtom_ptr;

}

node<Atom>* Residue::GetAaBetaAtom()
{
 node<Atom> *CBetaAtom_ptr,*cursor_atom_ptr;
 CBetaAtom_ptr=new node<Atom>;

  int AtomId,AaId;
  char *AtomName,*AaName,*Chain;

  AtomName=new char;
  AaName=new char;
  Chain=new char;
  
  double n,av_x,av_y,av_z,av_occ,av_B,av_asa,av_radius;
   n=0;av_x=0;av_y=0;av_z=0;av_asa=0;av_radius=0;
   
  cursor_atom_ptr=head_atom_ptr;
  while(cursor_atom_ptr!=NULL)
   {
     if(strncmp(cursor_atom_ptr->data.GetAtomName(),"CB ",3)==0)
      {
       av_x=av_x+cursor_atom_ptr->data.GetAtom_x();
       av_y=av_y+cursor_atom_ptr->data.GetAtom_y();
       av_z=av_z+cursor_atom_ptr->data.GetAtom_z();
       av_asa=av_asa+cursor_atom_ptr->data.GetAtom_asa();
       av_radius=av_radius+cursor_atom_ptr->data.GetAtom_radius();
       AtomId=cursor_atom_ptr->data.GetResidueNum();
       AaId=cursor_atom_ptr->data.GetResidueNum();
       strcpy(AtomName,cursor_atom_ptr->data.GetAtomName());
       strcpy(AaName,cursor_atom_ptr->data.GetResidueName());
       strcpy(Chain,cursor_atom_ptr->data.GetChainType());
       n=n+1;
      }
      cursor_atom_ptr=cursor_atom_ptr->next;
   }

       if(n!=0)
       {
       av_x=av_x/n;
       av_y=av_y/n;
       av_z=av_z/n;
       av_asa=av_asa/n;
       av_radius=av_radius/n;

    CBetaAtom_ptr->data.SetAtomIndex(AtomId);
    CBetaAtom_ptr->data.SetAtomName(AtomName);
    CBetaAtom_ptr->data.SetResidueName(AaName);
    CBetaAtom_ptr->data.SetChainType(Chain);
    CBetaAtom_ptr->data.SetResidueNum(AaId);    
    CBetaAtom_ptr->data.SetAtom_x(av_x);
    CBetaAtom_ptr->data.SetAtom_y(av_y);
    CBetaAtom_ptr->data.SetAtom_z(av_z);
    CBetaAtom_ptr->data.SetAtom_asa(av_asa);
    CBetaAtom_ptr->data.SetAtomRadius(av_radius);
    }
   else
    {CBetaAtom_ptr=NULL;}

  return CBetaAtom_ptr;

}
    
    
node<Atom>* Residue::GetAaSideChainGeometricalCenter()
{
 node<Atom> *SideChainCenterAtom_ptr,*cursor_atom_ptr;
 SideChainCenterAtom_ptr=new node<Atom>;

  int AtomId,AaId;
  char *AtomName,*AaName,*Chain;

  AtomName=new char;
  AaName=new char;
  Chain=new char;
  
  double n,av_x,av_y,av_z,av_occ,av_B,av_asa,av_radius;
   n=0;av_x=0;av_y=0;av_z=0;av_asa=0;av_radius=0;
   
  cursor_atom_ptr=head_atom_ptr;

  if(strncmp(cursor_atom_ptr->data.GetResidueName(),"GLY",3)!=0)
  {
  //cout<<cursor_atom_ptr->data.GetResidueName()<<endl; 
  while(cursor_atom_ptr!=NULL)
   {
     if(strncmp(cursor_atom_ptr->data.GetAtomName(),"CA ",3)!=0&&
        strncmp(cursor_atom_ptr->data.GetAtomName(),"N  ",3)!=0&&
        strncmp(cursor_atom_ptr->data.GetAtomName(),"C  ",3)!=0&&
        strncmp(cursor_atom_ptr->data.GetAtomName(),"O  ",3)!=0)
      {
       av_x=av_x+cursor_atom_ptr->data.GetAtom_x();
       av_y=av_y+cursor_atom_ptr->data.GetAtom_y();
       av_z=av_z+cursor_atom_ptr->data.GetAtom_z();
       av_asa=av_asa+cursor_atom_ptr->data.GetAtom_asa();
       av_radius=av_radius+cursor_atom_ptr->data.GetAtom_radius();
       AtomId=cursor_atom_ptr->data.GetResidueNum();
       AaId=cursor_atom_ptr->data.GetResidueNum();
       strcpy(AtomName,"s_c");
       strcpy(AaName,cursor_atom_ptr->data.GetResidueName());
       strcpy(Chain,cursor_atom_ptr->data.GetChainType());
       n=n+1;
      }//end if
      cursor_atom_ptr=cursor_atom_ptr->next;
   }//end while
 
   } 

       if(n!=0)
       {
       av_x=av_x/n;
       av_y=av_y/n;
       av_z=av_z/n;
       av_asa=av_asa/n;
       av_radius=av_radius/n;
       
    SideChainCenterAtom_ptr->data.SetAtomIndex(AtomId);
    SideChainCenterAtom_ptr->data.SetAtomName(AtomName);
    SideChainCenterAtom_ptr->data.SetResidueName(AaName);
    SideChainCenterAtom_ptr->data.SetChainType(Chain);
    SideChainCenterAtom_ptr->data.SetResidueNum(AaId);    
    SideChainCenterAtom_ptr->data.SetAtom_x(av_x);
    SideChainCenterAtom_ptr->data.SetAtom_y(av_y);
    SideChainCenterAtom_ptr->data.SetAtom_z(av_z);
    SideChainCenterAtom_ptr->data.SetAtom_asa(av_asa);
    SideChainCenterAtom_ptr->data.SetAtomRadius(av_radius);
     }
     else
     {SideChainCenterAtom_ptr=NULL;}
 
  return SideChainCenterAtom_ptr;
}

node<Atom>* Residue::GetAaGeometricalCenter()
{
 node<Atom> *AaCenterAtom_ptr,*cursor_atom_ptr;
 AaCenterAtom_ptr=new node<Atom>;

  int AtomId,AaId;
  char *AtomName,*AaName,*Chain;

  AtomName=new char;
  AaName=new char;
  Chain=new char;
  
  double n,av_x,av_y,av_z,av_occ,av_B,av_asa,av_radius;
   n=0;av_x=0;av_y=0;av_z=0;av_asa=0;av_radius=0;
   
  cursor_atom_ptr=head_atom_ptr;
  while(cursor_atom_ptr!=NULL)
   {
       av_x=av_x+cursor_atom_ptr->data.GetAtom_x();
       av_y=av_y+cursor_atom_ptr->data.GetAtom_y();
       av_z=av_z+cursor_atom_ptr->data.GetAtom_z();
       av_asa=av_asa+cursor_atom_ptr->data.GetAtom_asa();
       av_radius=av_radius+cursor_atom_ptr->data.GetAtom_radius();
       AtomId=cursor_atom_ptr->data.GetResidueNum();
       AaId=cursor_atom_ptr->data.GetResidueNum();
       strcpy(AtomName,"a_c");
       strcpy(AaName,cursor_atom_ptr->data.GetResidueName());
       strcpy(Chain,cursor_atom_ptr->data.GetChainType());
       n=n+1;
     // }
      cursor_atom_ptr=cursor_atom_ptr->next;
   }

       av_x=av_x/n;
       av_y=av_y/n;
       av_z=av_z/n;
       av_asa=av_asa/n;
       av_radius=av_radius/n;

    AaCenterAtom_ptr->data.SetAtomIndex(AtomId);
    AaCenterAtom_ptr->data.SetAtomName(AtomName);
    AaCenterAtom_ptr->data.SetResidueName(AaName);
    AaCenterAtom_ptr->data.SetChainType(Chain);
    AaCenterAtom_ptr->data.SetResidueNum(AaId);    
    AaCenterAtom_ptr->data.SetAtom_x(av_x);
    AaCenterAtom_ptr->data.SetAtom_y(av_y);
    AaCenterAtom_ptr->data.SetAtom_z(av_z);
    AaCenterAtom_ptr->data.SetAtom_asa(av_asa);
    AaCenterAtom_ptr->data.SetAtomRadius(av_radius);

  return AaCenterAtom_ptr;

}

double Residue::GetAaDistance()
{
return d;
}
    

char* Residue::GetAaPrChain()
{
   return p_name_chain;

}

char*  Residue::GetAaOneLetterCode()
{
 return one_letter_aa;
}

node<Residue>* Residue::GetAaNeighborList() 
{
  return Nb_head_aa_ptr;
}

void Residue::SetAaName(char *input)
{
 aa_name=new char;
 strcpy(aa_name,input);
}

void Residue::SetAaNum(int num)
{
 aa_num=num;
}

void Residue::SetAaChain(char *input)
{
 aa_chain=new char;
 strcpy(aa_chain,input);
}

void Residue::SetAaAtomList(Atom& atom)
{
   insert_an_node(head_atom_ptr,tail_atom_ptr,atom);

}


void Residue::SetAaDistance(double input)
{
 d=input;
}

void Residue::SetAaNeighborAaList(Residue& residue)
{
   insert_an_node(Nb_head_aa_ptr,Nb_tail_aa_ptr,residue);
}

void Residue::SetSortedNeighborAaList()
{
 node<Residue> *temp,*cursor_aa_ptr1,*cursor_aa_ptr2;
 cursor_aa_ptr1=Nb_head_aa_ptr;
 cursor_aa_ptr2=Nb_head_aa_ptr;

 while(cursor_aa_ptr1!=NULL)
 {
 cursor_aa_ptr2=cursor_aa_ptr1;

 while(cursor_aa_ptr2!=NULL)
 {
 if(cursor_aa_ptr2->data.GetAaDistance()<cursor_aa_ptr1->data.GetAaDistance())
   {
    temp->data=cursor_aa_ptr1->data;
    cursor_aa_ptr1->data=cursor_aa_ptr2->data;
    cursor_aa_ptr2->data=temp->data;
   }
  cursor_aa_ptr2=cursor_aa_ptr2->next;
 }
  cursor_aa_ptr1=cursor_aa_ptr1->next;
 }
}

void Residue::SetAaPrChain(char *input)
{
 p_name_chain=new char;
 strcpy(p_name_chain,input);
  
}

void Residue::SetAaOneLetterCode(char *input)
{
  int i;
 one_letter_aa=new char;
 string aa_list[22]={"CYS","HIS","ILE","MET","SER","VAL","ALA","GLY","LEU","PRO",
                    "THR","PHE","ARG","TYR","TRP","ASP","ASN","ASX","GLU","GLN","GLX","LYS"};

 string one_letter_list[22]={"C","H","I","M","S","V","A","G","L","P","T","F","R","Y","W","D","N","B","E","Q","Z","K"};
    
   for(i=0;i<22;i++)
   {
     //cout<<input<<" "<<aa_list[i].c_str()<<endl;

    if(strncmp(input,aa_list[i].c_str(),4)==0)
    { 
     strcpy(one_letter_aa,one_letter_list[i].c_str());
     break;
    }
   else 
   {
    strcpy(one_letter_aa,"X");
   }
   }

}

//void Residue::operator =(Residue*& source)
//{
// aa_name=new char;
// strcpy(aa_name,source->aa_name);

// aa_num=source->aa_num;

// aa_chain=new char;
// strcpy(aa_chain,source->aa_chain);

// head_atom_ptr=source->head_atom_ptr;

//}


ostream& operator<<(ostream& os,Residue& re)
{

  //if(strcmp(re.GetAaName(),"GLY")!=0&&(re.omega>=0&&re.omega<=180))
  //{
  os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<re.GetAaName()
<<setw(5)<<re.GetAaNum();
  os<<setw(10)<<re.GetAaPrChain();
  os<<setw(10)<<re.GetAaOneLetterCode();
  os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
  //<<setprecision(2)<<setw(7)<<re.GetAaAs()
  //<<setprecision(2)<<setw(7)<<re.GetAaMs()
  //<<setprecision(2)<<setw(7)<<re.GetAaAverageCurvature()
 // <<setprecision(2)<<setw(7)<<re.GetAaCalphaAtom()->data.GetAtom_x()
 // <<setprecision(2)<<setw(7)<<re.GetAaCalphaAtom()->data.GetAtom_y()
 // <<setprecision(2)<<setw(7)<<re.GetAaCalphaAtom()->data.GetAtom_z()
 // <<setprecision(2)<<setw(7)<<re.GetAaSideChainGeometricalCenter()->data.GetAtom_x()
 // <<setprecision(2)<<setw(7)<<re.GetAaSideChainGeometricalCenter()->data.GetAtom_y()
 // <<setprecision(2)<<setw(7)<<re.GetAaSideChainGeometricalCenter()->data.GetAtom_z()
  <<setprecision(2)<<setw(9)<<re.omega;
  //}

  return os; 

}

