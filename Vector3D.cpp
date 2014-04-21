/*
 *  Vector3D.cpp
 *  VectorTest
 *
 *  Created by Mark Szymczyk on 8/13/05.
 *
 */

#include "Vector3D.h"

// Constructors

// Default constructor
Vector3D::Vector3D(void)
{
    x = 0.0;
    y = 0.0;
    z = 0.0;    
}

// Constructor with arguments
Vector3D::Vector3D(double xValue, double yValue, double zValue)
{
    x = xValue;
    y = yValue;
    z = zValue;
}

// Copy constructor
Vector3D::Vector3D(const Vector3D& theVector)
{
	x = theVector.x;
	y = theVector.y;
    z = theVector.z;
}

// Destructor
Vector3D::~Vector3D(void)
{
	
}


// Overloaded operators
Vector3D& Vector3D::operator= (const Vector3D& theVector)
{
	if (this != &theVector) {	// Avoid self-assignment
		x = theVector.x;
		y = theVector.y;
        z = theVector.z;
	}
	
	return *this;
}

Vector3D Vector3D::operator+ (const Vector3D& theVector) const
{
    // Returns the sum of two vectors
    
    Vector3D sum;
    
    sum.x = x + theVector.x;
    sum.y = y + theVector.y;
    sum.z = z + theVector.z;
    
    return sum;
}

Vector3D Vector3D::operator- (const Vector3D& theVector) const
{
    // Returns the difference of two vectors
    
    Vector3D diff;
    
    diff.x = x - theVector.x;
    diff.y = y - theVector.y;
    diff.z = z - theVector.z;
    
    return diff;
    
}

Vector3D Vector3D::operator* (double theScalar) const
{
    // Performs multiplication between a vector and a scalar.
    
    Vector3D product;
    
    product.x = x * theScalar;
    product.y = y * theScalar;
    product.z = z * theScalar;
    
    return product;
    
}


// Accessors
double Vector3D::GetXComponent(void)
{
	return x;
}

void Vector3D::SetXComponent(double xValue)
{
	x = xValue;
}

double Vector3D::GetYComponent(void)
{
	return y;
}

void Vector3D::SetYComponent(double yValue)
{
	y = yValue;
}

double Vector3D::GetZComponent(void)
{
	return z;
}

void Vector3D::SetZComponent(double zValue)
{
	z = zValue;
}


// Vector operations
double Vector3D::Vector3D::CalculateLength(void)
{
    // length = square root of (xSquared + ySquared + zSquared)
    
	double length;
	double xSquared = x * x;
	double ySquared = y * y;
    double zSquared = z * z;
    
	length = sqrt(xSquared + ySquared + zSquared);
	return length;
}

double Vector3D::DotProduct(Vector3D v)
{
    // If u = (x1, y1, z1) and v = (x2, y2, z2)
    // dot product = x1x2 + y1y2 + z1z2
    
	double dotProduct = (x * v.x) + (y * v.y) + (z * v.z);
	return dotProduct;
}

Vector3D Vector3D::CrossProduct(Vector3D v)
{
    // If u = (x1, y1, z1) and v = (x2, y2, z2)
    // cross product is a vector (x3, y3, z3)
    // x3 = y1z2 - z1y2
    // y3 = z1x2 - x1z2
    // z3 = x1y2 - y1x2
    
    Vector3D crossProduct;
    
    crossProduct.x = (y * v.z) - (z * v.y);
    crossProduct.y = (z * v.x) - (x * v.z);
    crossProduct.z = (x * v.y) - (y * v.x);
    
    return crossProduct;
}

void Vector3D::Normalize(void)
{
    // Takes the current vector and makes it a vector with length 1.
    // To do this, take each component and divide by its length.
    
    double length = CalculateLength();
    
    // Avoid a possible divide by 0 error.
    if (length == 0) {
        return;
    }
    
    // Multiplying by 1 over the length is faster than
    // dividing by the length.
    double oneOverLength = 1 / length;
    
    x = x * oneOverLength;
    y = y * oneOverLength;
    z = z * oneOverLength;
}

// Equivalents to overloading the +, -, and * operators
Vector3D Vector3D::VectorSum(Vector3D u, Vector3D v)
{
	Vector3D sum;
	sum.x = u.x + v.x;
	sum.y = u.y + v.y;
	sum.z = u.z + v.z;
	
	return sum;
}

Vector3D Vector3D::VectorDifference(Vector3D u, Vector3D v)
{
	Vector3D diff;
	diff.x = u.x - v.x;
	diff.y = u.y - v.y;
	diff.z = u.z - v.z;
	
	return diff;
}

Vector3D Vector3D::ScalarProduct(Vector3D v, double scalar)
{
	// This function multiplies a vector by a scalar quantity
	// Many books refer to the dot product of two vectors as a
	// scalar product, but I can’t think of a better function name
	// for multiplying a vector by a scalar than ScalarProduct.
	Vector3D product;
	product.x = v.x * scalar;
	product.y = v.y * scalar;
	product.z = v.z * scalar;
	
	return product;
}

// Static functions
double Vector3D::DotProduct(Vector3D u, Vector3D v)
{
    // If u = (x1, y1, z1) and v = (x2, y2, z2)
    // dot product = x1x2 + y1y2 + z1z2
    
	double x1 = u.x;
	double y1 = u.y;
    double z1 = u.z;
    
	double x2 = v.x;
	double y2 = v.y;
    double z2 = v.z;
    
	double dotProduct = (x1 * x2) + (y1 * y2) + (z1 * z2);
	return dotProduct;
}

Vector3D Vector3D::CrossProduct(Vector3D u, Vector3D v)
{
    // If u = (x1, y1, z1) and v = (x2, y2, z2)
    // cross product is a vector (x3, y3, z3)
    // x3 = y1z2 - z1y2
    // y3 = z1x2 - x1z2
    // z3 = x1y2 - y1x2
    
	double x1 = u.x;
	double y1 = u.y;
    double z1 = u.z;
    
	double x2 = v.x;
	double y2 = v.y;
    double z2 = v.z;
    
    Vector3D crossProduct;
    
    crossProduct.x = (y1 * z2) - (z1 * y2);
    crossProduct.y = (z1 * x2) - (x1 * z2);
    crossProduct.z = (x1 * y2) - (y1 * x2);
    
    return crossProduct;
}

ostream& operator<<(ostream& os, Vector3D& vec){
  return os << "(" << vec.x << ", "
	    << vec.y << ", "<<vec.z<<")";
}
