#include "1705021_Header.h"
#include "bitmap_image.h"

int drawaxes;
double scalar, theta;

int levelOfRecursion, imageDimension, numberOfObjects, numberOfPointLights,numberOfSpotLights;
int floorWidth, tileWidth;

bitmap_image image;
int image_count;

vector<Object*> objectArray;
vector<Light> lightSourceArray;

Vector pos,u,r,l;

void drawAxes()
{
    if(drawaxes==1)
    {
        double len = floorWidth/2.0;
        glColor3f(1.0, 1.0, 0.0);
        glBegin(GL_LINES);
        {
            glVertex3f( len,0,0);
            glVertex3f(-len,0,0);

            glVertex3f(0,-len,0);
            glVertex3f(0, len,0);

            glVertex3f(0,0, len);
            glVertex3f(0,0,-len);
        }
        glEnd();
    }
}


void loadData()
{

    ifstream sceneFile;
    sceneFile.open("scene.txt");

    if(!sceneFile.is_open())
    {
        cout << "failed to open input file" << endl;
        exit(EXIT_FAILURE);
    }
    string lineSegment;

    sceneFile>>levelOfRecursion;
    sceneFile>>imageDimension;
    sceneFile>>numberOfObjects;

    for(int i=0; i<numberOfObjects; i++)
    {

        sceneFile >> lineSegment;

        Object* currentObject;
        double r, g, b;
        double ambient, diffuse, specular, recursiveReflection;
        int shine;


        if(lineSegment == "sphere")
        {

            Vector center;

            sceneFile>>center.x;
            sceneFile>>center.y;
            sceneFile>>center.z;

            double radius;

            sceneFile>>radius;

            currentObject = new Sphere(center, radius);


        }

        if(lineSegment == "triangle")
        {

            Vector vertices[3];

            for(int j=0; j<3; j++)
            {
                sceneFile>>vertices[j].x;
                sceneFile>>vertices[j].y;
                sceneFile>>vertices[j].z;
            }

            currentObject = new Triangle(vertices);

        }

        if(lineSegment == "general")
        {

            double eqnCoeffs[10];

            for(int i=0; i<10; i++)
            {
                sceneFile >> eqnCoeffs[i];
            }
            currentObject = new GeneralObject(eqnCoeffs);

            Vector p;

            sceneFile>>p.x;
            sceneFile>>p.y;
            sceneFile>>p.z;

            double l, w, h;

            sceneFile>>l;
            sceneFile>>w;
            sceneFile>>h;

            currentObject->setReferencePoint(p);

            currentObject->setLength(l);
            currentObject->setWidth(w);
            currentObject->setHeight(h);
        }


        sceneFile>>r;
        sceneFile>>g;
        sceneFile>>b;

        sceneFile>>ambient;
        sceneFile>>diffuse;
        sceneFile>>specular;
        sceneFile>>recursiveReflection;

        sceneFile>>shine;

        currentObject->setColor(r, g, b);
        currentObject->setCoEfficients(ambient, diffuse, specular, recursiveReflection);
        currentObject->setShine(shine);


        objectArray.push_back(currentObject);
    }

    sceneFile >> numberOfPointLights;

    for(int i=0; i< numberOfPointLights; i++)
    {

        Vector p;

        sceneFile>>p.x;
        sceneFile>>p.y;
        sceneFile>>p.z;

        double r, g, b;

        sceneFile>>r;
        sceneFile>>g;
        sceneFile>>b;

        Light light(false);
        light.setLightPos(p);
        light.setColor(r, g, b);

        lightSourceArray.push_back(light);
    }

    sceneFile >> numberOfSpotLights;

    for(int i=0; i< numberOfSpotLights; i++)
    {

        Vector p;

        sceneFile>>p.x;
        sceneFile>>p.y;
        sceneFile>>p.z;

        double r, g, b;

        sceneFile>>r;
        sceneFile>>g;
        sceneFile>>b;

        Vector q;

        sceneFile>>q.x;
        sceneFile>>q.y;
        sceneFile>>q.z;

        double angle;

        sceneFile>>angle;

        Light light(true);
        light.setLightPos(p);
        light.setColor(r, g, b);
        light.setLightDir(q);
        light.setCutOffAngle(angle);

        lightSourceArray.push_back(light);
    }

    sceneFile.close();

    Object* floor = new Floor(floorWidth,tileWidth);

    objectArray.push_back(floor);


    for(int i=0; i<objectArray.size(); i++)
    {
        objectArray.at(i)->print();
    }

    for(int i=0; i<lightSourceArray.size(); i++)
    {
        lightSourceArray.at(i).print();
    }
}

void capture()
{

    cout << "Capturing...\nFrom Eye ";
    pos.print();

    int imageWidth = imageDimension;
    int imageHeight = imageDimension;

    image.setwidth_height(imageWidth, imageHeight, true); /* true for filling pixels with bg color 0x00(black) */

    double viewAngle = fovY;
    double planeDistance = (windowHeight/2.0) / tan((viewAngle * pi) / 360.0);

    Vector topLeft = pos + l * planeDistance - r * (windowWidth/2.0) + u * (windowHeight/2.0);

    double du = windowWidth / (imageWidth * 1.0);
    double dv = windowHeight / (imageHeight * 1.0);

    topLeft = topLeft + r * (0.5 * du) - u * (0.5 * dv);

    int nearest;
    double t, tMin;

    for(int i=0; i<imageWidth; i++) /* left to right = (x,r,du,width) */
    {
        for(int j=0; j<imageHeight; j++) /* top to bottom = (y,u,dv,height) */
        {

            t = INF;
            nearest = -1;

            Vector curPixel = topLeft + r * i * du - u * j * dv;

            Ray ray(pos, curPixel - pos);
            double* color = new double[3];

            for(int k=0; k<objectArray.size(); k++)
            {
                double tempT = objectArray.at(k)->intersect(ray, color, 0);
                if(t>tempT && tempT >0)
                {
                    t = tempT;
                    nearest = k;
                }
            }


            if(nearest >= 0)
            {
                for(int k=0; k<3; k++)
                    color[k] = 0;

                tMin = objectArray.at(nearest)->intersect(ray, color, 1);

                for(int k=0; k<3; k++)
                    color[k] = (color[k]<0) ? 0 : (color[k] > 1)? 1:color[k];

                if(tMin >= nearDist && tMin <=farDist)
                    image.set_pixel(i,j, (int)round(color[0]*255), (int)round(color[1]*255), (int)round(color[2]*255));
            }
            delete[] color;
        }
    }


    stringstream count;
    count << ++image_count;
    string imageName = "Output_" + count.str() + ".bmp";
    image.save_image(imageName);
    cout << imageName<<" : Image Saved..!\n\n";
    image.clear();
}


void clearObjects()
{

    for(int i=0; i<objectArray.size(); i++)
    {
        delete objectArray[i];
    }

    objectArray.clear();
}

void clearLights()
{

    lightSourceArray.clear();
}

void drawObjects()
{
    for(int i=0; i<objectArray.size(); i++)
    {
        objectArray.at(i)->draw();
    }

    for(int i=0; i<lightSourceArray.size(); i++)
    {
        lightSourceArray.at(i).draw();
    }
}

void keyboardListener(unsigned char key, int x,int y)
{

    double angle = theta;
    switch(key)
    {

    case '1':     /// rotate left
        l = l.rotate(angle,u);
        r = l.crossProduct(u);
        break;
    case '2':     /// rotate right
        l = l.rotate(-angle,u);
        r = l.crossProduct(u);
        break;
    case '3':     /// look up
        l = l.rotate(angle,r);
        u = r.crossProduct(l);
        break;
    case '4':     /// look down
        l = l.rotate(-angle,r);
        u = r.crossProduct(l);
        break;
    case '5':    ///  tilt clockwise
        u = u.rotate(+angle,l);
        r = l.crossProduct(u);
        break;
    case '6':    ///  tilt anticlockwise
        u = u.rotate(-angle,l);
        r = l.crossProduct(u);
        break;
    case '0':
        capture();
        break;

    default:
        break;
    }
}
void specialKeyListener(int key, int x,int y)
{

    switch(key)
    {
    case GLUT_KEY_DOWN:		//down arrow key
        pos = pos - l * scalar;
        break;
    case GLUT_KEY_UP:		// up arrow key
        pos = pos + l * scalar;
        break;

    case GLUT_KEY_RIGHT:
        pos = pos + r * scalar;
        break;
    case GLUT_KEY_LEFT:
        pos = pos - r * scalar;
        break;

    case GLUT_KEY_PAGE_UP:
        pos = pos + u * scalar;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos - u * scalar;
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:

        break;

    case GLUT_KEY_END:

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
        }
        break;

    case GLUT_RIGHT_BUTTON:
        if(state == GLUT_DOWN)
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

    gluLookAt(pos.x,pos.y,pos.z,	pos.x+l.x,pos.y+l.y,pos.z+l.z,	u.x,u.y,u.z);

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/

    drawAxes();

    drawObjects();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate()
{
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    //codes for initialization

    scalar = 5;
    theta = 5;

    floorWidth = 1000;
    tileWidth = 20;

    u = Vector(0,0,1);

    r = Vector(-1/sqrt(2),1/sqrt(2),0);

    l = Vector(-1/sqrt(2),-1/sqrt(2),0);

    pos = Vector(120,100,50);

    drawaxes=0;

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
    gluPerspective(fovY, aspectRatio, nearDist, farDist);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}



int main(int argc, char **argv)
{


    glutInit(&argc,argv);
    glutInitWindowSize(windowHeight, windowWidth);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("1705021_Ray_Tracer");


    init();
    loadData();

    if((atexit(clearLights) != 0) || (atexit(clearObjects) != 0))
    {
        cout << "atexit(): atexit() function registration failed" << endl;
        exit(EXIT_FAILURE);
    }


    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL
    return 0;
}
