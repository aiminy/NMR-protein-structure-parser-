#include "point.h"
#include <iostream>

using namespace std;

#ifndef TRIANGLE_H
#define TRIANGLE_H

class Triangle
{

 public:
	Triangle();
        void SetTriangleIndex(int num);
        void SetTriangleThreePoint(Point a, Point b, Point c);

 private:

    int index;
    Point first,second,third; 
 
 friend ostream& operator<<(ostream& os, Triangle& triangle);

};
#endif
