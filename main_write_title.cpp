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

#include "Vector3D.h"


using namespace std;
    
int main()
{

  ofstream aa_file("all_protein_curvature.txt",ios::out);

  aa_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"p_n"
    <<setiosflags(ios::right)<<setw(5)<<"aa_t"<<setw(2)<<" ";
  aa_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "
  <<setw(7)<<"ca_av"
  <<setw(7)<<"ca_sd"
  <<setw(7)<<"ca_n"
  <<setw(7)<<"cx_av"
  <<setw(7)<<"cx_sd"
  <<setw(7)<<"cx_n"
  <<setw(7)<<"su_av"
  <<setw(7)<<"su_sd"
  <<setw(7)<<"su_n"
  <<setw(7)<<"bu_av"
  <<setw(7)<<"bu_sd"
  <<setw(7)<<"bu_n"
  <<setw(7)<<"ov_av"
  <<setw(7)<<"ov_sd"
  <<setw(7)<<"ov_n"<<endl;
  return 0;

}

