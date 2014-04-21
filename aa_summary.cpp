#include "aa_summary.h"

Aa_Summary::Aa_Summary()
{
  Aa_name=NULL;
  average=-1;
  std=-1;
  num=-1;
}

char* Aa_Summary::GetAaName()
{
   return Aa_name;
}

double Aa_Summary::GetAaAverage()
{
   return average;
}

double Aa_Summary::GetAaStd()
{
   return std;
}

int Aa_Summary::GetAaNum()
{
   return num;
}

void Aa_Summary::SetAaName(char *input)
{
 Aa_name=new char;
 strcpy(Aa_name,input);
}

void Aa_Summary::SetAaAverage(double input)
{
 average=input;
}

void Aa_Summary::SetAaStd(double input)
{
 std=input;
}

void Aa_Summary::SetAaNum(int input)
{
  num=input;
}
