#include<bits/stdc++.h>
#include<windows.h>
#include<GL/glut.h>

#define pi (2*acos(0.0))
#define INF numeric_limits<double>::infinity()
#define EPSILON 0.0001


#define AMBIENT 0
#define DIFFUSE 1
#define SPECULAR 2
#define REFLECTION 3

using namespace std;

int windowHeight = 700, windowWidth = 700;
int fovY = 80, aspectRatio = 1;
int nearDist = 1, farDist = 1000;
bool refract = false;

struct Vector
{

    double x,y,z;

    Vector(double x = 0.0,double y = 0.0,double z = 0.0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector operator+(const Vector& other) const
    {
        return Vector(this->x + other.x,this->y + other.y,this->z + other.z);
    }
    Vector operator-(const Vector& other) const
    {
        return Vector(this->x - other.x,this->y - other.y,this->z - other.z);
    }

    Vector operator*(double scalar) const
    {
        return Vector(this->x*scalar,this->y*scalar,this->z*scalar);
    }

    Vector operator/(double scalar) const
    {
        return Vector(this->x/scalar,this->y/scalar,this->z/scalar);
    }

    double dotProduct(const Vector& other)const
    {
        return this->x*other.x + this->y*other.y + this->z*other.z;
    }
    Vector crossProduct(const Vector& other)const
    {
        return Vector(this->y*other.z - this->z*other.y,
                      this->z*other.x - this->x*other.z,
                      this->x*other.y - this->y*other.x
                     );
    }

    double valueOfVector()
    {
        double d=sqrt(x*x + y*y + z*z);
        if(d==0) return EPSILON;
        else return d;
    }
    void normalize()
    {
        *this = *this/this->valueOfVector();
    }

    Vector rotate(double degAngle,const Vector& k)
    {
        Vector v = *this;
        double A = degAngle*pi/180;
        return  v*cos(A) + k.crossProduct(v)*sin(A) +  k * k.dotProduct(v) * (1- cos(A));
    }

    void print()
    {
        cout << "Position(x,y,z): (" << this->x << "," << this->y << "," << this->z << ")" << endl;
    }
};




class Light
{
public:
    Vector lightPosition;
    double color[3];

    Vector lightDirection;
    double cutOffAngle;
    bool   isSpotLight;

    Light(bool isSpotLight)
    {
        this->isSpotLight = isSpotLight;
    }

    void setLightPos(Vector p)
    {
        lightPosition = p;
    }
    void setLightDir(Vector p)
    {
        Vector direction = p - lightPosition;
        direction.normalize();

        if(isSpotLight)
            lightDirection = direction;
    }
    void setCutOffAngle(double angle)
    {
        if(isSpotLight)
            cutOffAngle = angle;
    }
    void setColor(double r, double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void draw()
    {

        glPushMatrix();
        glTranslatef(lightPosition.x, lightPosition.y, lightPosition.z);

        int slices=48,stacks=20;
        int radius = 1;

        Vector points[100][100];
        int i,j;
        double h,r;

        //generate points
        for(i=0; i<=stacks; i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0; j<=slices; j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }

        //draw quads using generated points
        glColor3f(color[0],color[1],color[2]);
        for(i=0; i<stacks; i++)
        {
            for(j=0; j<slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }
                glEnd();
            }
        }
        glPopMatrix();
    }

    void print()
    {
        if(isSpotLight)
        {
            cout<<"\n<Spot Light Source> \n";
            cout<<"Cut Off Angle : "<<cutOffAngle<<endl;
            cout<<"Direcion ";
            lightDirection.print();
        }
        else cout<<"\n<Point Light Source> \n";

        lightPosition.print();
        cout<<"Color(r,g,b): (" << color[0] << "," << color[1] << "," << color[2] << ")" << endl << endl;
    }
};

class Ray
{
public:
    Vector start;
    Vector dir;

    Ray() {}

    Ray(Vector s, Vector d)
    {
        this->start = s;
        this->dir = d;
        this->dir.normalize();
    }
};



class Object
{
public:

    Vector referencePoint;
    double height;
    double width;
    double length;
    double color[3];
    double coeffs[4]; // reflection coefficients
    int shine; // exponent term of specular component

    Object() {}
    virtual void draw() {}
    virtual void print() {}

    virtual double intersect(Ray ray, double color[3], int level)
    {
        return -1;
    }

    void illuminate(Ray ray, Vector intersectionPoint, Vector normal, double color[3], int level);


    virtual Vector normalVector(Vector p)
    {
        return Vector(0, 0, 1);
    }

    void setReferencePoint(Vector p)
    {
        referencePoint = p;
    }

    void setLength(double l)
    {
        length = l;
    }

    void setWidth(double w)
    {
        width = w;
    }

    void setHeight(double h)
    {
        height = h;
    }

    void setColor(double r, double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void setShine(int sh)
    {
        shine = sh;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double reflect)
    {
        coeffs[AMBIENT] = ambient;
        coeffs[DIFFUSE] = diffuse;
        coeffs[SPECULAR] = specular;
        coeffs[REFLECTION] = reflect;
    }

    void printBasics()
    {

        cout<<"Color : (r,g,b) = (";
        for(int i=0; i<3; i++)
        {
            cout<<color[i]<<",";
        }
        cout<<")"<<endl;

        cout<<"Co-Efficients : ";
        cout<<"Ambient("<<coeffs[AMBIENT]<<") ";
        cout<<"Diffuse("<<coeffs[DIFFUSE]<<") ";
        cout<<"Specular("<<coeffs[SPECULAR]<<") ";
        cout<<"Reflection("<<coeffs[REFLECTION]<<")";
        cout<<endl;

        cout<<"Shininess : ";
        cout<<shine<<"\n\n";
    }


};

double myClamp(double low,double high,double currentValue)
{

    if(currentValue < low) return low;
    else if(currentValue > high) return high;
    else return currentValue;
}

extern vector<Light> lightSourceArray;
extern vector<Object*> objectArray;
extern int levelOfRecursion;


void Object::illuminate(Ray ray, Vector intersectionPoint, Vector normal, double color[3], int level)
{

    /* step 3: ambient coloring at intersection point */
    for(int k=0; k<3; k++)
        color[k] = this->color[k]*coeffs[AMBIENT];

    for(int i=0; i<lightSourceArray.size(); i++)
    {

        Vector direction = lightSourceArray[i].lightPosition - intersectionPoint;
        direction.normalize();

        if(lightSourceArray[i].isSpotLight)
        {

            /* Spot Light additional check code starts here*/
            /* adjusting ray direction by reversing it */
            double cos_phi = lightSourceArray[i].lightDirection.dotProduct(direction*(-1.0));

            if(cos_phi <= cos(lightSourceArray[i].cutOffAngle*(pi/180.0)))
                continue;
            //else cout<<"active spotlight"<<endl;

        }

        /* step 4: casting ray in betn light Source and intersection point */
        Ray rayTemp(intersectionPoint + direction*(EPSILON), direction);

        bool inShadow = false;
        double* dummyColor = new double[3];
        double dist = (intersectionPoint - lightSourceArray[i].lightPosition).valueOfVector();

        /* step 5: checking presence of shadow & obscureness by any object*/
        for(int k=0; k<objectArray.size(); k++)
        {
            double tempT = objectArray.at(k)->intersect(rayTemp, dummyColor, 0);

            /* if any object obstructs the ray , it can't travel the max dist i.e. dist betn light & intersect point  */
            if(dist>tempT && tempT > 0)
            {
                inShadow = true;
                break;
            }
        }

        delete[] dummyColor;

        if(!inShadow)
        {

            /* step 6: calculating lambert Value using normal, ray1 (i.e. in betn light Source and intersect point)*/
            double lambertValue = normal.dotProduct(rayTemp.dir); /*  L.N  */

            /* It seems, L is not right as it should come from lightSrc to intersectPoint but
            it is actually from calculation perspective bcz => currently L has direction from
            intersectPoint(slightly above of it) to lightSrc.And from slide we see that we exctly
            need this direction what we have currently to find theta being compatible with N direction
            lambert value is actually equivalent to cos_theta value */

            /*step 7:  finding  reflected ray */
            Vector reflectedRay = normal*(2*lambertValue) - (rayTemp.dir); /*  R = 2(L.N)*N - L  */

            /* step 8: calculating phong Value using normal , reflected ray*/
            double phongValue = reflectedRay.dotProduct(ray.dir*(-1)); /*  R.V  */
            /* (-1) here bcz => V has direction from eye to intersectPoint. But from slide we see that
            we need reverse of it to find phi being compatible with R direction. phong value is actually
            equivalent to cos_phi value */

            if(lambertValue < 0) lambertValue = 0;
            if(phongValue < 0) phongValue = 0;

            for(int k=0; k<3; k++)
            {
                /* step 9: diffuse coloring at intersection point */
                color[k] += lightSourceArray[i].color[k]*coeffs[DIFFUSE]*lambertValue*this->color[k];
                /* step 10: specular coloring at intersection point */
                color[k] += lightSourceArray[i].color[k]*coeffs[SPECULAR]*pow(phongValue, this->shine)*this->color[k];

            }

        }
    }

    for(int k=0; k<3; k++)
        color[k] = (color[k]<0) ? 0 : (color[k] > 1)? 1:color[k];


    if(level >= levelOfRecursion)
        return;

    /* Recursive Reflection implementation starts here*/

    /* (-1) here bcz => old V or new L has direction from eye to intersectPoint. But from slide we see that
    we need reverse of it to find theta being compatible with N direction */
    double rLambertValue = normal.dotProduct(ray.dir*(-1)); /* Now old V => (ray.dir) becomes new L */
    Vector recursiveReflectionDirection = normal*(2*rLambertValue) - (ray.dir*(-1));
    recursiveReflectionDirection.normalize();

    /* step 1: constructing reflected ray from intersection point avoiding self intersection*/
    Ray recursiveReflectionRay(intersectionPoint + recursiveReflectionDirection*(EPSILON), recursiveReflectionDirection);

    double rNearest=-1, rT=INF, rTempT, rTmin;
    double* colorReflected = new double[3];

    /* step 2: find tmin from the nearest intersecting object, using intersect(),as done in the capture()*/
    for(int k=0; k<objectArray.size(); k++)
    {
        rTempT = objectArray.at(k)->intersect(recursiveReflectionRay, colorReflected, 0);
        if(rT>rTempT && rTempT > 0)
        {
            rT = rTempT;
            rNearest = k;
        }
    }

    /* step 3: if found, calling intersect() Recursive Reflection */
    if(rNearest >= 0)
    {
        for(int k=0; k<3; k++)
            colorReflected[k] = 0;

        rTmin = objectArray.at(rNearest)->intersect(recursiveReflectionRay, colorReflected, level+1);
        if(rTmin > 0 && rTmin <INF)
        {

            /*step 4: updating color using the impact of reflection  */
            for(int k=0; k<3; k++)
                color[k] += colorReflected[k]*this->coeffs[REFLECTION];

        }
    }

    delete[] colorReflected;

    /* Recursive Refraction implementation starts here*/

    if(refract)
    {

        double cosi = myClamp(-1, 1,normal.dotProduct(ray.dir)); /* N = normal , I = ray.dir */
        double etai = 1, etat = 1.5; /* glass we assume */
        /*  etai is the index of refraction of the medium the ray is in before entering the second medium (here our objects)
            etat is the refraction index of the material of the physical object the ray has either hit or is about to leave (glass, water, etc.)
        */
        Vector n = normal;
        if (cosi < 0)
        {
            cosi = -cosi;
        }
        else
        {
            swap(etai, etat);
            n = normal*(-1.0);
        }
        double eta = etai / etat;
        double k = 1 - eta * eta * (1 - cosi * cosi);

        /* step 1: constructing refracted ray from intersection point*/
        Vector recursiveRefractionDirection;

        if(k<0)
        {
            recursiveRefractionDirection = Vector(0,0,0);
        }
        else
        {
            recursiveRefractionDirection = ray.dir * eta + n*(eta * cosi - sqrt(k));
            recursiveRefractionDirection.normalize();
        }

        
        Ray recursiveRefractionRay(intersectionPoint + recursiveRefractionDirection*(EPSILON), recursiveRefractionDirection);


        double rrNearest=-1, rrT=INF, rrTempT, rrTmin;
        double* colorRefracted = new double[3];

        /* step 2: find tmin from the nearest intersecting object, using intersect(),as done in the capture()*/
        for(int k=0; k<objectArray.size(); k++)
        {
            rrTempT = objectArray.at(k)->intersect(recursiveRefractionRay, colorRefracted, 0);
            if(rrT>rrTempT && rrTempT > 0)
            {
                rrT = rrTempT;
                rrNearest = k;
            }
        }


        /* step 3: if found, calling intersect() Recursive Refraction */
        if(rrNearest >= 0)
        {
            for(int k=0; k<3; k++)
                colorRefracted[k] = 0;

            rrTmin = objectArray.at(rrNearest)->intersect(recursiveRefractionRay, colorRefracted, level+1);
            if(rrTmin > 0 && rrTmin <INF)
            {

                /*step 4: updating color using the impact of refraction  */
                for(int k=0; k<3; k++)
                    color[k] += colorRefracted[k]*this->coeffs[REFLECTION]/2;

            }
        }

        delete[] colorRefracted;
    }

}






class Sphere: public Object
{
public:

    Sphere(Vector center, double radius)
    {
        referencePoint = center;
        length = radius;
    }

    Vector normalVector(Vector p)
    {
        Vector normal = p - referencePoint;

        normal.normalize();

        return normal;
    }

    void print()
    {
        cout<<"<Sphere>\n";
        cout<<"Center ";
        referencePoint.print();

        cout<<"Radius : ";
        cout<<length<<endl;

        printBasics();
    }

    double intersect(Ray ray, double color[3], int level)
    {

        Vector oc = ray.start - referencePoint; /* ref point is center */
        double a = ray.dir.dotProduct(ray.dir);
        double b = 2 * ray.dir.dotProduct(oc);
        double c = oc.dotProduct(oc) - length*length;  /* length  is radius */

        double minPosT = INF;
        double det = b*b - 4*a*c;

        if(det >=0)
        {
            double d = sqrt(det);
            double t1 = (-b+d)/(2*a);
            double t2 = (-b-d)/(2*a);

            if(t1>0 || t2>0)
            {
                if(t1 <= 0)
                {
                    minPosT = t2;
                }
                else if(t2 <= 0)
                {
                    minPosT = t1;
                }
                else
                {
                    minPosT = min(t1,t2);
                }

                if((minPosT < nearDist) || (minPosT > farDist))
                    return INF;
            }
            else
            {
                return minPosT;
            }
        }
        else
        {
            return minPosT;
        }


        if(level == 0) return minPosT;

        /* step 0: finding intersection point */
        Vector intersectionPoint = ray.start + ray.dir*(minPosT);

        /* step 1: getColorAt intersection point => surface color of SPHERE  => color already set  for SPHERE at beginning  */

        /* step 2: finding normal at intersection point */
        Vector normal =  normalVector(intersectionPoint);
        if(normal.dotProduct(ray.dir) > 0) normal = normal*(-1);

        /* Get Illumination Via (Ambient + Diffuse + Specular) Coloring and also (Reflection + Refraction) */
        illuminate(ray,intersectionPoint,normal,color,level);


        return minPosT;
    }

    void draw()
    {
        glPushMatrix();
        glTranslatef(referencePoint.x, referencePoint.y, referencePoint.z);

        int slices=48,stacks=20;
        int radius = length;

        Vector points[100][100];
        int i,j;
        double h,r;

        //generate points
        for(i=0; i<=stacks; i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0; j<=slices; j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }

        //draw quads using generated points
        glColor3f(color[0],color[1],color[2]);
        for(i=0; i<stacks; i++)
        {
            for(j=0; j<slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }
                glEnd();
            }
        }
        glPopMatrix();
    }

};

class Triangle: public Object
{
public:

    Vector vertices[3];

    Triangle(Vector p[3])
    {
        for(int i=0; i<3; i++)
        {
            vertices[i] = p[i];
        }
    }

    Vector normalVector(Vector p)
    {
        Vector e1 = vertices[1] - vertices[0];
        Vector e2 = vertices[2] - vertices[0];

        Vector normal = e1.crossProduct(e2);

        normal.normalize();

        return normal;
    }
    void print()
    {
        cout<<"<Triangle>\n";
        cout<<"Vertices : \n";
        for(int i=0; i<3; i++)
        {
            vertices[i].print();
        }

        printBasics();
    }

    double determinant(double a[3][3])
    {

        double a1 = a[0][0];
        double a2 = a[1][0];
        double a3 = a[2][0];

        double b1 = a[0][1];
        double b2 = a[1][1];
        double b3 = a[2][1];

        double c1 = a[0][2];
        double c2 = a[1][2];
        double c3 = a[2][2];

        return  a1*(b2*c3 - b3*c2) + b1*(a3*c2 - a2*c3) + c1*(a2*b3 - a3*b2);
    }

    double intersect(Ray ray, double color[3], int level)
    {
        /*
        double matA[3][3]    = {{vertices[0].x - vertices[1].x, vertices[0].x - vertices[2].x, ray.dir.x},
                                {vertices[0].y - vertices[1].y, vertices[0].y - vertices[2].y, ray.dir.y},
                                {vertices[0].z - vertices[1].z, vertices[0].z - vertices[2].z, ray.dir.z}
        };

        double matBeta[3][3] = {{vertices[0].x - ray.start.x, vertices[0].x - vertices[2].x, ray.dir.x},
                                {vertices[0].y - ray.start.y, vertices[0].y - vertices[2].y, ray.dir.y},
                                {vertices[0].z - ray.start.z, vertices[0].z - vertices[2].z, ray.dir.z}
        };


        double matGamma[3][3] ={{vertices[0].x - vertices[1].x, vertices[0].x - ray.start.x, ray.dir.x},
                                {vertices[0].y - vertices[1].y, vertices[0].y - ray.start.y, ray.dir.y},
                                {vertices[0].z - vertices[1].z, vertices[0].z - ray.start.z, ray.dir.z}
        };

        double matT[3][3]     ={{vertices[0].x - vertices[1].x, vertices[0].x - vertices[2].x, vertices[0].x - ray.start.x},
                                {vertices[0].y - vertices[1].y, vertices[0].y - vertices[2].y, vertices[0].y - ray.start.y},
                                {vertices[0].z - vertices[1].z, vertices[0].z - vertices[2].z, vertices[0].z - ray.start.z}
        };

        double minPosT = INF;
        double detA = determinant(matA);
        */

        double detBase, detBeta, detGamma, detT;
        double minPosT = INF;

        Vector a = vertices[0];
        Vector b = vertices[1];
        Vector c = vertices[2];

        detBase = (a.x-b.x)*((a.y-c.y)*ray.dir.z-(a.z-c.z)*ray.dir.y);
        detBase += (a.x-c.x)*((a.z-b.z)*ray.dir.y-(a.y-b.y)*ray.dir.z);
        detBase += ray.dir.x*((a.y-b.y)*(a.z-c.z)-(a.z-b.z)*(a.y-c.y));

        detBeta =  (a.x-ray.start.x)*( (a.y-c.y)*ray.dir.z - (a.z-c.z)*ray.dir.y );
        detBeta += (a.x-c.x)*( (a.z-ray.start.z)*ray.dir.y - (a.y-ray.start.y)*ray.dir.z );
        detBeta += ray.dir.x*( (a.y-ray.start.y)*(a.z-c.z) - (a.z-ray.start.z)*(a.y-c.y) );

        detGamma =  (a.x-b.x)*( (a.y-ray.start.y)*ray.dir.z - (a.z-ray.start.z)*ray.dir.y );
        detGamma += (a.x-ray.start.x)*( (a.z-b.z)*ray.dir.y - (a.y-b.y)*ray.dir.z );
        detGamma += ray.dir.x*( (a.y-b.y)*(a.z-ray.start.z) - (a.z-b.z)*(a.y-ray.start.y) );

        detT =  (a.x-b.x)*( (a.y-c.y)*(a.z-ray.start.z) - (a.z-c.z)*(a.y-ray.start.y) );
        detT += (a.x-c.x)*( (a.z-b.z)*(a.y-ray.start.y) - (a.y-b.y)*(a.z-ray.start.z) );
        detT += (a.x-ray.start.x)*( (a.y-b.y)*(a.z-c.z) - (a.z-b.z)*(a.y-c.y) );

        double detA = detBase;

        if(detA != 0)
        {

            //double beta = determinant(matBeta)/detA;
            //double gamma = determinant(matGamma)/detA;
            //minPosT = determinant(matT)/detA;

            double beta = detBeta/detA;
            double gamma = detGamma/detA;
            minPosT = detT/detA;

            if((beta + gamma < 1) && (beta > 0) && (gamma > 0) && (minPosT >= nearDist) && (minPosT <= farDist))
            {

                if(level == 0) return minPosT;

                /* step 0: finding intersection point */
                Vector intersectionPoint = ray.start + ray.dir*(minPosT);

                /* step 1: getColorAt intersection point => surface color of TRIANGLE  => color already set  for TRIANGLE at beginning  */

                /* step 2: finding normal at intersection point */
                Vector normal = normalVector(Vector(0,0,0));
                if(normal.dotProduct(ray.dir) > 0) normal = normal*(-1);

                /* Get Illumination Via (Ambient + Diffuse + Specular) Coloring and also (Reflection + Refraction) */
                illuminate(ray,intersectionPoint,normal,color,level);

            }
            else
            {
                minPosT = INF;
            }
        }
        return minPosT;
    }

    void draw()
    {
        glPushMatrix();
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(vertices[0].x, vertices[0].y, vertices[0].z);
            glVertex3f(vertices[1].x, vertices[1].y, vertices[1].z);
            glVertex3f(vertices[2].x, vertices[2].y, vertices[2].z);
        }
        glEnd();
        glPopMatrix();
    }
};

class GeneralObject: public Object
{
public:
    double a,b,c,d,e,f,g,h,i,j;

    GeneralObject(double eqnCoeffs[10])
    {

        a = eqnCoeffs[0];
        b = eqnCoeffs[1];
        c = eqnCoeffs[2];
        d = eqnCoeffs[3];
        e = eqnCoeffs[4];
        f = eqnCoeffs[5];
        g = eqnCoeffs[6];
        h = eqnCoeffs[7];
        i = eqnCoeffs[8];
        j = eqnCoeffs[9];

    }


    Vector normalVector(Vector p)
    {

        Vector normal;

        normal.x = 2*a*p.x + d*p.y + e*p.z + g;
        normal.y = 2*b*p.y + d*p.x + f*p.z + h;
        normal.z = 2*c*p.z + e*p.x + f*p.y + i;

        normal.normalize();

        return normal;
    }


    double intersect(Ray ray, double color[3], int level)
    {


        double xd = ray.dir.x;
        double yd = ray.dir.y;
        double zd = ray.dir.z;

        double xo = ray.start.x;
        double yo = ray.start.y;
        double zo = ray.start.z;

        /* co-efficients of t2,t,const */
        double aq = a*xd*xd + b*yd*yd + c*zd*zd + d*xd*yd + e*xd*zd + f*yd*zd;
        double bq = 2*a*xo*xd + 2*b*yo*yd + 2*c*zo*zd + d*(xo*yd + yo*xd) + e*(xo*zd + zo*xd) + f*(yo*zd + yd*zo) + g*xd + h*yd + i*zd;
        double cq = a*xo*xo + b*yo*yo + c*zo*zo + d*xo*yo + e*xo*zo + f*yo*zo + g*xo + h*yo + i*zo + j;

        double det = bq*bq - 4*aq*cq;
        double minPosT = INF;
        double t1,t2;

        if(aq == 0)
        {
            t1 = (-1)*(cq/bq);
            t2 = t1;
        }
        else
        {
            if(det >= 0)
            {
                double d = sqrt(det);
                t1 =(-bq - d)/(2*aq);
                t2 =(-bq + d)/(2*aq);
            }
            else
            {
                return minPosT;
            }
        }

        /* step 0: finding intersection point */
        Vector intersectionPointT1 = ray.start + ray.dir*(t1);
        Vector intersectionPointT2 = ray.start + ray.dir*(t2);

        bool insideCubeT1 = true;
        bool insideCubeT2 = true;

        if(length != 0)
        {
            if(abs(intersectionPointT1.x) < abs(referencePoint.x) || abs(intersectionPointT1.x) > abs(referencePoint.x) + length)
            {
                insideCubeT1 = false;
            }

            if(abs(intersectionPointT2.x) < abs(referencePoint.x) || abs(intersectionPointT2.x) > abs(referencePoint.x) + length)
            {
                insideCubeT2 = false;
            }
        }

        if(width != 0)
        {
            if(abs(intersectionPointT1.y) < abs(referencePoint.y) || abs(intersectionPointT1.y) > abs(referencePoint.y) + width)
            {
                insideCubeT1 = false;
            }

            if(abs(intersectionPointT2.y) < abs(referencePoint.y) || abs(intersectionPointT2.y) > abs(referencePoint.y) + width)
            {
                insideCubeT2 = false;
            }
        }

        if(height != 0)
        {
            if(abs(intersectionPointT1.z) < abs(referencePoint.z) || abs(intersectionPointT1.z) > abs(referencePoint.z) + height)
            {
                insideCubeT1 = false;
            }

            if(abs(intersectionPointT2.z) < abs(referencePoint.z) || abs(intersectionPointT2.z) > abs(referencePoint.z) + height)
            {
                insideCubeT2 = false;
            }
        }

        if(!insideCubeT1 && !insideCubeT2)
        {
            return INF;
        }
        else if(!insideCubeT1)
        {
            minPosT = t2;
        }
        else if(!insideCubeT2)
        {
            minPosT = t1;
        }
        else
        {
            minPosT = min(t1,t2);
        }

        if(minPosT <= 0 || (minPosT < nearDist) || (minPosT > farDist)) return INF;

        if(level == 0) return minPosT;

        Vector intersectionPoint = ray.start + ray.dir*(minPosT);

        /* step 1: getColorAt intersection point => surface color of GENERAL  => color already set  for GENERAL at beginning  */

        /* step 2: finding normal at intersection point */
        Vector normal = normalVector(intersectionPoint);
        if(normal.dotProduct(ray.dir) > 0) normal = normal*(-1);

        /* Get Illumination Via (Ambient + Diffuse + Specular) Coloring and also (Reflection + Refraction) */
        illuminate(ray,intersectionPoint,normal,color,level);

        return minPosT;
    }

    void draw() {}

    void print()
    {
        cout<<"<General>\n";
        cout<<"Equation Co-Efficients : ";

        cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" ";
        cout<<f<<" "<<g<<" "<<h<<" "<<i<<" "<<j<<" ";

        cout<<endl;

        printBasics();
    }
};

class Floor: public Object
{
public:

    Floor(double floorWidth, double tileWidth)
    {

        referencePoint.x = -floorWidth/2;
        referencePoint.y = -floorWidth/2;
        referencePoint.z = 0;

        length = tileWidth;
        width = floorWidth;

        setCoEfficients(0.4, 0.4, 0.4, 0.1);
        setColor(0.5, 0.5, 0.5);
        setShine(3);
    }

    void print()
    {

        cout << "<Floor>\n";
        cout << "Reference Point ";
        referencePoint.print();
        cout << "Length : " << length << endl;
        cout << "Width : " << width << endl << endl;
    }

    double intersect(Ray ray, double* color, int level)
    {

        double d = 0, minPosT = INF;

        /* step 2: finding normal at intersection point */
        Vector normal(0,0,1), intersectionPoint;
        if(normal.dotProduct(ray.dir) > 0) normal = normal*(-1); // upside down

        if(normal.dotProduct(ray.dir) != 0)
        {

            /* step 0: finding intersection point */
            minPosT = -1.0 * (d + normal.dotProduct(ray.start)) / normal.dotProduct(ray.dir);

            intersectionPoint = ray.start + ray.dir*(minPosT);

            if((intersectionPoint.x >= this->referencePoint.x && intersectionPoint.x <= this->referencePoint.x + width) && (intersectionPoint.y >= this->referencePoint.y && intersectionPoint.y <= this->referencePoint.y + width) && (minPosT >= nearDist) && (minPosT <= farDist))
            {

                if(level == 0) return minPosT;

                /* step 1: getColorAt intersection point */
                int pixelX = (int)((referencePoint.x - intersectionPoint.x) / length);
                int pixelY = (int)((referencePoint.y - intersectionPoint.y) / length);

                double colorValue = ((pixelX + pixelY)%2 == 0)? 0 : 1;

                for(int k=0; k<3; k++)
                    this->color[k] = colorValue;

                /* Get Illumination Via (Ambient + Diffuse + Specular) Coloring and also (Reflection + Refraction) */
                illuminate(ray,intersectionPoint,normal,color,level);
            }
            else
            {
                minPosT = INF;
            }
        }
        return minPosT;
    }

    void draw()
    {
        for(double i=referencePoint.x ; i<=width/2 ; i+=length)
        {
            int row = floor((i - referencePoint.x)/length);

            for(double j=referencePoint.y ; j<=width/2 ; j+=length)
            {
                int col = floor((j - referencePoint.y)/length);

                if((row+col)%2 == 0)
                {
                    glColor3f(0.0, 0.0, 0.0);
                }
                else
                {
                    glColor3f(1.0, 1.0, 1.0);
                }


                glBegin(GL_QUADS);
                {
                    glVertex3f(i, j, 0);
                    glVertex3f(i, j+length, 0);
                    glVertex3f(i+length, j+length, 0);
                    glVertex3f(i+length, j, 0);
                }
                glEnd();
            }
        }
    }
};



