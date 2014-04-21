#ifndef AMINO_H
#define AMIN0_H

class Aa_type{

public:

Aa_type();
char GetName();
double GetMean();
double GetSd();
double GetNum(); 
void SetName(char *input);
void SetMean(double input);
void SetSd(double input);
void SetNum(int input); 



private:

char   *name;
double mean;
double sd;
int num;


};


