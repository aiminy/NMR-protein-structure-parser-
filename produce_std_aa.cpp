#include "cal_tot_asa_of_res.h"

void produce_std_aa(vector <std_aa>&st_aa)
{
  std_aa std_aa_c;

  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"ARG");
  strcpy(std_aa_c.r_atom,"CZ ");
  st_aa.push_back(std_aa_c);
 
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"LYS");
  strcpy(std_aa_c.r_atom,"NZ ");
  st_aa.push_back(std_aa_c);

  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"SER");
  strcpy(std_aa_c.r_atom,"OG ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"ILE");
  strcpy(std_aa_c.r_atom,"CD1");
  st_aa.push_back(std_aa_c);

  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"THR");
  strcpy(std_aa_c.r_atom,"CB ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"ASN");
  strcpy(std_aa_c.r_atom,"CG ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"TYR");
  strcpy(std_aa_c.r_atom,"OH ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"PRO");
  strcpy(std_aa_c.r_atom,"CG ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"LEU");
  strcpy(std_aa_c.r_atom,"CG ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"GLU");
  strcpy(std_aa_c.r_atom,"CD ");
  st_aa.push_back(std_aa_c);

  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"ASP");
  strcpy(std_aa_c.r_atom,"CG ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"ALA");
  strcpy(std_aa_c.r_atom,"CB ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"CYS");
  strcpy(std_aa_c.r_atom,"SG ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"VAL");
  strcpy(std_aa_c.r_atom,"CB ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"MET");
  strcpy(std_aa_c.r_atom,"CE ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"GLN");
  strcpy(std_aa_c.r_atom,"CD ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"HIS");
  strcpy(std_aa_c.r_atom,"CG ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"PHE");
  strcpy(std_aa_c.r_atom,"CZ ");
  st_aa.push_back(std_aa_c);
  
  std_aa_c.res_id=new char[4];
  std_aa_c.r_atom=new char[4];
  strcpy(std_aa_c.res_id,"TRP");
  strcpy(std_aa_c.r_atom,"CH2");
  st_aa.push_back(std_aa_c);
}

