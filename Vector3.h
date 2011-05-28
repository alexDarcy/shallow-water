class Vector3
{
  public:
    double x;
    double y;
    double z;

    Vector3();
    Vector3(double x,double y,double z);
    Vector3(const Vector3 & v);

    Vector3& operator= (const Vector3 & v);

    Vector3& operator+= (const Vector3 & v);
    Vector3&  operator-= (const Vector3 & v);
    Vector3&  operator*= (const double a);
    Vector3&  operator/= (const double a);

    Vector3 operator+ (const Vector3 & v);
    Vector3 operator- (const Vector3 & v);
    Vector3 operator* (const double a);
    Vector3 operator/ (const double a);

    Vector3 crossProduct(const Vector3 & v);
    void normalize();
    double norm();
};
