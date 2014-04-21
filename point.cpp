#include "point.h"

Point::Point()
{
 index=0;
 x=0;
 y=0;
 z=0;
}

void Point::SetPointIndex(int num)
{
 index=num;
}

void Point::SetPointCoordinate(double a, double b, double c)
{
 x=a;
 y=b;
 z=c;
}

int Point::GetPointIndex()
{
 return index;
}

double Point::GetPointX()
{
return x;
}

double Point::GetPointY()
{
return y;
}

double Point::GetPointZ()
{
return z;
}

void Point::operator =(Point& source)
{
 index=source.index; 
 x=source.x;
 y=source.y;
 z=source.z;
}

ostream& operator<<(ostream& os,Point& point)
{
  os<<" PointIndex "<<point.index<<" x "<<point.x<<" y "<<point.y<<" z "<<point.z<<endl;
  return os;
}
 
