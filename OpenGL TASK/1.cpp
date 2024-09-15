#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))


double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double maxLen = 60;
double minLen = 40;

struct point
{

    double x,y,z;

    point(double x = 0.0,double y = 0.0,double z = 0.0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    point operator+(const point& other) const
    {
        return point(this->x + other.x,this->y + other.y,this->z + other.z);
    }
    point operator-(const point& other) const
    {
        return point(this->x - other.x,this->y - other.y,this->z - other.z);
    }

    point operator*(double scalar) const
    {
        return point(this->x*scalar,this->y*scalar,this->z*scalar);
    }

    point operator/(double scalar) const
    {
        return point(this->x/scalar,this->y/scalar,this->z/scalar);
    }

    double dot(const point& other)const
    {
        return this->x*other.x + this->y*other.y + this->z*other.z;
    }
    point cross(const point& other)const
    {
        return point(this->y*other.z - this->z*other.y,
                     this->z*other.x - this->x*other.z,
                     this->x*other.y - this->y*other.x
                    );
    }

    double value()
    {
        return sqrt(x*x + y*y + z*z);
    }

   point rotate(double degAngle,const point& k)
    {
        point v = *this;
        double A = degAngle*pi/180;
        return  v*cos(A) + k.cross(v)*sin(A) +  k * k.dot(v) * (1- cos(A));
    }
};

struct point pos, u, r, l;

void drawAxes()
{
    if(drawaxes==1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f( 100,0,0);
            glVertex3f(-100,0,0);

            glVertex3f(0,-100,0);
            glVertex3f(0, 100,0);

            glVertex3f(0,0, 100);
            glVertex3f(0,0,-100);
        }
        glEnd();
    }
}


void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);
        {
            for(i=-8; i<=8; i++)
            {

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);
            }
        }
        glEnd();
    }
}



void drawSphere(double radius,int slices,int stacks)
{
    glColor3f(0, 1, 0);
    struct point points[100][100];
    int i,j;
    double h,r;
    //generate points
    for(i=0; i<=stacks; i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
        for(j=0; j<=slices; j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*(pi/2.0)); /// one eighth of sphere
            points[i][j].y=r*sin(((double)j/(double)slices)*(pi/2.0));
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0; i<stacks; i++)
    {
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
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
                /*
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                */
            }
            glEnd();

        }
    }
}
void drawSpheres()
{
    double radius = maxLen - minLen ;
    glPushMatrix();            /* at {+x , +y} corner : ABOVE */
    {
        glTranslated(minLen, minLen, minLen);
        drawSphere(radius, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();            /* at {+x , -y} corner : ABOVE */
    {
        glTranslated(minLen, -minLen, minLen);
        glRotated(-90, 0, 0, 1);
        drawSphere(radius, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();            /* at {-x , +y} corner : ABOVE */
    {
        glTranslated(-minLen, minLen, minLen);
        glRotated(90, 0, 0, 1);
        drawSphere(radius, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();            /* at {-x , -y} corner : ABOVE */
    {
        glTranslated(-minLen, -minLen, minLen);
        glRotated(180, 0, 0, 1);
        drawSphere(radius, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();            /* at {+x , +y} corner : BELOW */
    {
        glTranslated(minLen, minLen, -minLen);
        glRotated(90, 0, 1, 0);
        drawSphere(radius, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();            /* at {+x , -y} corner : BELOW */
    {
        glTranslated(minLen, -minLen, -minLen);
        glRotated(180, 1, 0, 0);
        drawSphere(radius, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();            /* at {-x , y} corner : BELOW */
    {
        glTranslated(-minLen, minLen, -minLen);
        glRotated(-180, 0, 1, 0);
        drawSphere(radius, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();            /* at {-x , -y} corner : BELOW */
    {
        glTranslated(-minLen, -minLen, -minLen);
        glRotated(-90, 0, 0, 1);
        glRotated(180, 1, 0, 0);
        drawSphere(maxLen - minLen, 50, 50);
    }
    glPopMatrix();
}

void drawCylinder(double radius, double height, int slices, int stacks)
{
    glColor3f(0.0,0.0,1.0);
    struct point points[100][100];
    int i, j;
    double h, r;
    //generate points
    for (i = 0; i <= stacks; i++)
    {
        h = height * sin(((double)i / (double)stacks)*(pi/2.0));
        //r=radius*cos(((double)i/(double)stacks)*(pi/2));
        //h = height;
        r = radius;
        for (j = 0; j <= slices; j++)
        {
            points[i][j].x = r * cos(((double)j / (double)slices)*(pi/2.0)); /// one fourth of cylinder
            points[i][j].y = r * sin(((double)j / (double)slices)*(pi/2.0));
            points[i][j].z = h;
        }
    }
    //draw quads using generated points
    for (i = 0; i<stacks; i++)
    {
        //glColor3f(0, 1, 0);
        for (j = 0; j<slices; j++)
        {
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
            }
            glEnd();
        }
    }
}

void drawCylinders()
{
    double radius = maxLen - minLen ;
    glPushMatrix();            /* at {+x , +y} corner */
    {
        glTranslated(minLen, minLen, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();           /* at {+x , -y} corner */
    {
        glTranslated(minLen, -minLen, 0);
        glRotated(-90, 0, 0, 1);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();           /* at {-x , +y} corner */
    {
        glTranslated(-minLen, minLen, 0);
        glRotated(90, 0, 0, 1);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();          /* at {-x , -y} corner */
    {
        glTranslated(-minLen, -minLen, 0);
        glRotated(180, 0, 0, 1);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();         /* along {+x , +y} corner To {+x , -y} corner : ABOVE */
    {
        glTranslated(minLen, 0, minLen);
        glRotated(90, 1, 0, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();         /* along {+x , +y} corner To {+x , -y} corner : BELOW */
    {
        glTranslated(minLen, 0, -minLen);
        glRotated(-90, 1, 0, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();         /* along {-x , +y} corner To {-x , -y} corner : ABOVE */
    {
        glTranslated(-minLen, 0, minLen);
        glRotated(90, 0, 0, 1);
        glRotated(-90, 0, 1, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();


    glPushMatrix();         /* along {-x , +y} corner To {-x , -y} corner : BELOW */
    {
        glTranslated(-minLen, 0, -minLen);
        glRotated(90, 0, 0, 1);
        glRotated(90, 0, 1, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();         /* along {+x , +y} corner To {-x , +y} corner : ABOVE */
    {
        glTranslated(0, minLen, minLen); //ok
        glRotated(90, 0, 0, 1);
        glRotated(90, 1, 0, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();         /* along {+x , +y} corner To {+x , -y} corner : BELEOW */
    {
        glTranslated(0, minLen, -minLen);
        glRotated(90, 0, 0, 1);
        glRotated(-90, 1, 0, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();         /* along {+x , -y} corner To {-x , -y} corner : ABOVE */
    {
        glTranslated(0, -minLen, minLen);
        glRotated(-90, 0, 0, 1);
        glRotated(90, 1, 0, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();

    glPushMatrix();         /* along {+x , -y} corner To {-x , -y} corner : BELEOW */
    {
        glTranslated(0, -minLen, -minLen);
        glRotated(-90, 0, 0, 1);
        glRotated(-90, 1, 0, 0);
        drawCylinder(radius, minLen, 50, 50);
    }
    glPopMatrix();
}

void drawSquare(double a)
{
    glColor3f(1.0,0.0,0.0);
    glBegin(GL_QUADS);
    {
        glVertex3f( a, a,0);
        glVertex3f( a,-a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(-a, a,0);
    }
    glEnd();
}


void drawSquares()
{
    glPushMatrix(); /* along +z axis */
    {
        glTranslated(0, 0, maxLen);
        drawSquare(minLen);
    }
    glPopMatrix();

    glPushMatrix(); /* along -z axis */
    {
        glTranslated(0, 0, -maxLen);
        drawSquare(minLen);
    }
    glPopMatrix();

    glPushMatrix(); /* along +y axis */
    {
        glTranslated(0, maxLen, 0);  /// translate third
        glRotated(90, 1, 0, 0);     /// rotate second
        drawSquare(minLen);        /// create first
    }
    glPopMatrix();

    glPushMatrix(); /* along -y axis */
    {
        glTranslated(0, -maxLen, 0);
        glRotated(90, 1, 0, 0);
        drawSquare(minLen);
    }
    glPopMatrix();

    glPushMatrix(); /* along +x axis */
    {
        glTranslated(maxLen, 0, 0);
        glRotated(90, 0, 1, 0);
        drawSquare(minLen);
    }
    glPopMatrix();

    glPushMatrix();  /* along -x axis */
    {
        glTranslated(-maxLen, 0, 0);
        glRotated(90, 0, 1, 0);
        drawSquare(minLen);
    }
    glPopMatrix();
}


void keyboardListener(unsigned char key, int x,int y)
{

    double angle = 2.0;
    switch(key)
    {

    case '1':     /// rotate left
        l = l.rotate(angle,u);
        r = l.cross(u);
        break;
    case '2':     /// rotate right
        l = l.rotate(-angle,u);
        r = l.cross(u);
        break;
    case '3':     /// look up
        l = l.rotate(angle,r);
        u = r.cross(l);
        break;
    case '4':     /// look down
        l = l.rotate(-angle,r);
        u = r.cross(l);
        break;
    case '5':    ///  tilt clockwise
        u = u.rotate(+angle,l);
        r = l.cross(u);
        break;
    case '6':    ///  tilt anticlockwise
        u = u.rotate(-angle,l);
        r = l.cross(u);
        break;

    default:
        break;
    }
}


void specialKeyListener(int key, int x,int y)
{

    double scaler = 3.0 ;
    switch(key)
    {
    case GLUT_KEY_DOWN:		//down arrow key
        pos = pos - l * scaler;
        break;
    case GLUT_KEY_UP:		// up arrow key
        pos = pos + l * scaler;
        break;

    case GLUT_KEY_RIGHT:
        pos = pos + r * scaler;
        break;
    case GLUT_KEY_LEFT:
        pos = pos - r * scaler;
        break;

    case GLUT_KEY_PAGE_UP:
        pos = pos + u * scaler;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos - u * scaler;
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:
        if(minLen > 0) minLen--;
        break;

    case GLUT_KEY_END:
        if(minLen < maxLen) minLen++;
        break;

    default:
        break;
    }
}


void mouseListener(int button, int state, int x, int y) 	//x, y is the x-y of the screen (2D)
{
    switch(button)
    {
    case GLUT_LEFT_BUTTON:
        if(state == GLUT_DOWN) 		// 2 times?? in ONE click? -- solution is checking DOWN or UP
        {
            drawaxes=1-drawaxes;
        }
        break;
    default:
        break;
    }
}



void display()
{

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    //gluLookAt(0,0,200,	0,0,0,	0,1,0);
    gluLookAt(pos.x,pos.y,pos.z,
              pos.x + l.x,pos.y + l.y,pos.z + l.z,
              u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();
    drawSpheres();
    drawCylinders();
    drawSquares();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate()
{
    //angle+=0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    //codes for initialization
    drawgrid=0;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;


    u = point(0,0,1);

    r = point(-1/sqrt(2),1/sqrt(2),0);

    l = point(-1/sqrt(2),-1/sqrt(2),0);

    pos = point(180,180,20);

    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(80,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv)
{
    glutInit(&argc,argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("1705021_Camera");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}
