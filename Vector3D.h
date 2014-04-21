/*
 *  Vector3D.h
 *  VectorTest
 *
 *  Created by Mark Szymczyk on 8/13/05.
 *
 */

// Include math.h for the square root function
#include <math.h>
#include <iostream> 
using namespace std;

class Vector3D
{
protected:
    double x;
    double y;
    double z;
    
public:
    // Constructors and destructor
    Vector3D(void);
    Vector3D(double xValue, double yValue, double zValue);
    Vector3D(const Vector3D& theVector); // Copy constructor
    virtual ~Vector3D(void);
    
    // Overload operators
    Vector3D& operator= (const Vector3D& theVector);
    Vector3D operator+ (const Vector3D& theVector) const;
    Vector3D operator- (const Vector3D& theVector) const;
    Vector3D operator* (double theScalar) const;
    
    // Accessors
    double GetXComponent(void);
    void SetXComponent(double xValue);
    
    double GetYComponent(void);
    void SetYComponent(double yValue);
    
    double GetZComponent(void);
    void SetZComponent(double zValue);
    
    // Vector operations
    double CalculateLength(void);
    double DotProduct(Vector3D v);
    Vector3D CrossProduct(Vector3D v);
    void Normalize(void);
    
    // These are the equivalents of overloading the
    // +, -, and * operators.  If you wanted to use C 
    // instead of C++, you would use functions like these.
    static Vector3D VectorSum(Vector3D u, Vector3D v);
    static Vector3D VectorDifference(Vector3D u, Vector3D v);
    static Vector3D ScalarProduct(Vector3D v, double scalar);

    // Static functions. These are just different ways
    // of calling some of the vector operation functions.
    static double DotProduct(Vector3D u, Vector3D v);
    static Vector3D CrossProduct(Vector3D u, Vector3D v);

    friend ostream& operator<<(ostream& os, Vector3D& vec);

};

typedef Vector3D* Vector3DPtr;
typedef Vector3D** Vector3DHandle;
