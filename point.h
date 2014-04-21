#include <iostream>
#ifndef Point_H
#define Point_H

using namespace std;

class Point {

public:
	Point();
        void SetPointIndex(int num);
        void SetPointCoordinate(double a,double b,double c);
        int GetPointIndex();
        double GetPointX();
        double GetPointY();
        double GetPointZ();
 
        void operator =(Point& source);


private:
        int index;
        double x,y,z;

friend ostream& operator<<(ostream& os,Point& point);

};
#endif
