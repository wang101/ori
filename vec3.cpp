class vec
{
    public:
       double a;
       double b;
       double c;
       vec(double x,double y,double z)
       {
           a=x;
           b=y;
           c=z;
       }
       vec operator+(vec &v1)
       {
           return vec(v1.a+a,v1.b+b,v1.c+c);
       }
       vec operator-(vec &v1)
       {
           return vec(v1.a-a,v1.b-b,v1.c-c);
       }
       vec operator*(vec &v1)
       {
           return vec(v1.a*a,v1.b*b,v1.c*c);
       }
       double dot(vec v1)
       {
           double sum = v1.a*a+v1.b*b+v1.c*c;
           return sum;
       }
};
