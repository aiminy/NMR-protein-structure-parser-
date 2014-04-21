#include "triangle.h"

Triangle::Triangle()
{
 index=0;
 
}

void Triangle::SetTriangleIndex(int num)
{
 index=num;

}

void Triangle::SetTriangleThreePoint(Point a,Point b,Point c)
{
  first=a;
  second=b;
  third=c;
}

ostream& operator<<(ostream& os,Triangle& triangle)
{
  os<<triangle.first;
  os<<triangle.second;
  os<<triangle.third;
 
  //os<<"tets";
 return os;
}
