#ifndef AA_SUMMARY_H
#define AA_SUMMARY_H

#include "atom.h"
#include "linked_list.h"


class Aa_Summary
{
 public:
    Aa_Summary();

    char* GetAaName();
    double GetAaAverage();
    double GetAaStd();
    int    GetAaNum();
 
    void SetAaName(char *input);
    void SetAaAverage(double input);
    void SetAaStd(double input);
    void SetAaNum(int input);

 private:

   char* Aa_name; 
   double average;
   double std;
   int num;
   
};
#endif
