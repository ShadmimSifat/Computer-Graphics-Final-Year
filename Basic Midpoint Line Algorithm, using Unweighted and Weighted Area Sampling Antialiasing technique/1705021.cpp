#include <iostream>
#include <fstream>
#include <vector>
#include "bitmap_image.h"

using namespace std;

#define VERTICAL_SLOPE 1
#define HORIZONTAL_SLOPE 2
#define MIDPOINT_SLOPE 3     /* m<=1 */
#define NON_MIDPOINT_SLOPE 4 /* m>1 */
#define NEGATIVE_SLOPE_1 5   /* m>=-1 */
#define NEGATIVE_SLOPE_2 6   /* m<-1 */

/* PROTOTYPES */
class Line;
class Rectangle;

/* GLOBAL VARIABLES */
bitmap_image image;
vector<Line> lineList;
ifstream inputFile;

int pixelWidth, pixelHeight;
int totalLines;

/* CLASSES */
class Point
{

public:
    double x, y;

    Point(double x = 0, double y = 0)
    {
        this->x = x;
        this->y = y;
    }

    double dot(const Point &other) const
    {
        return this->x * other.x + this->y * other.y;
    }

    Point operator-(const Point &other) const
    {
        return Point(this->x - other.x, this->y - other.y);
    }

    double value()
    {
        return sqrt(x * x + y * y);
    }

    void print()
    {
        cout << "(" << x << " , " << y << ")" << endl;
    }
};

class Rectangle
{

public:
    Point topLeft, bottomLeft, bottomRight, topRight;

    Rectangle(Point start, Point end, double slope, double t = 1.0)
    {

        double theta = atan(slope);
        // cout << slope << endl;
        topLeft = Point(start.x - t * sin(theta), start.y + t * cos(theta));
        bottomLeft = Point(start.x + t * sin(theta), start.y - t * cos(theta));

        topRight = Point(end.x - t * sin(theta), end.y + t * cos(theta));
        bottomRight = Point(end.x + t * sin(theta), end.y - t * cos(theta));

        // print();
    }

    double getIntensity(int x, int y)
    {

        double intensity;
        int totalOverlaps = 0;

        vector<Point> pixels_16;

        dividePixelBy_16SubPixels(x, y, pixels_16);

        // for(int k = 0;k<pixels_16.size();k++)
        // {
        //     pixels_16[k].print();
        // }

        for (int i = 0; i < pixels_16.size(); i++)
        {
            if (isInside(pixels_16[i]))
            {
                totalOverlaps++;
            }
        }

        intensity = totalOverlaps / 16.0;
        // cout << intensity << endl;

        return intensity;
    }

    double funcValue(Point p1, Point p2, Point m)
    {
        return m.y - p1.y - (p2.y - p1.y) / (p2.x - p1.x) * (m.x - p1.x);
    }

    bool isInside(Point m)
    {

        // m.print();

        double topLineFunc1 = funcValue(topLeft, topRight, m);
        double topLineFunc2 = funcValue(topLeft, topRight, bottomRight);

        double bottomLineFunc1 = funcValue(bottomLeft, bottomRight, m);
        double bottomLineFunc2 = funcValue(bottomLeft, bottomRight, topLeft);

        double leftLineFunc1 = funcValue(topLeft, bottomLeft, m);
        double leftLineFunc2 = funcValue(topLeft, bottomLeft, topRight);

        double rightLineFunc1 = funcValue(topRight, bottomRight, m);
        double rightLineFunc2 = funcValue(topRight, bottomRight, topLeft);

        return (topLineFunc1 * topLineFunc2 >= 0 && bottomLineFunc1 * bottomLineFunc2 >= 0 && leftLineFunc1 * leftLineFunc2 >= 0 && rightLineFunc1 * rightLineFunc2 >= 0);
    }

    void dividePixelBy_16SubPixels(int xx, int yy, vector<Point> &pixels_16)
    {

        double x = xx - 0.5;
        double y = yy - 0.5;

        // cout << x << " " << y << endl;

        double len = 0.25;

        x = x + len / 2;

        for (int i = 0; i < 4; i++)
        {

            double x1 = x + i * len;
            double y1 = y + len / 2;

            // cout << x1 << " " << y1 << endl;
            pixels_16.push_back(Point(x1, y1));

            for (int j = 1; j < 4; j++)
            {

                y1 = y1 + len;

                // cout << x1 << " " << y1 << endl;
                pixels_16.push_back(Point(x1, y1));
            }
            // cout<<endl;
        }
    }

    void print()
    {
        cout << "Top Left: ";
        topLeft.print();
        cout << "Bottom Left: ";
        bottomLeft.print();
        cout << "Bottom Right: ";
        bottomRight.print();
        cout << "Top Right: ";
        topRight.print();
    }
};

class Line
{

public:
    Point start, end;
    int color[3];
    int slopeType;
    int parallelYlineDist;
    double slope;

    Line(Point start, Point end)
    {
        slopeType = slopeTypeOfLine(start, end);

        if (slopeType == MIDPOINT_SLOPE)
        {
            this->start = start;
            this->end = end;
        }
        else if (slopeType == VERTICAL_SLOPE)
        {

            if (start.y < end.y)
            {
                this->start = start;
                this->end = end;
            }
            else
            {
                this->start = end;
                this->end = start;
            }
        }
        else if (slopeType == HORIZONTAL_SLOPE)
        {

            if (start.x < end.x)
            {
                this->start = start;
                this->end = end;
            }
            else
            {
                this->start = end;
                this->end = start;
            }
        }
        else if (slopeType == NON_MIDPOINT_SLOPE)
        {
            /* reflect w.r.to y = x line so that m>1 turns into m<=1 */
            this->start = Point(start.y, start.x);
            this->end = Point(end.y, end.x);
        }
        else if (slopeType == NEGATIVE_SLOPE_1)
        {
            /* reflect w.r.to y = a line so that m>=-1 turns into m<=1  N.B: a = parallelYlineDist here  */

            /* assuming 1st point inserted by tester is always on the lefthand side
             w.r.to other point inserted 2nd ; same Assumption for NEGATIVE_SLOPE_2 */

            this->start = end; /* lower end point */
            parallelYlineDist = end.x;
            this->end = Point(-start.x + 2 * parallelYlineDist, start.y);
        }
        else if (slopeType == NEGATIVE_SLOPE_2) /*  ****here reflect first then swap second****  */
        {

            /* 1st reflect w.r.to y = a line so that  m<-1 turns into m>1  and then*/
            this->start = end;
            parallelYlineDist = end.x;
            this->end = Point(-start.x + 2 * parallelYlineDist, start.y);

            /* 2nd reflect w.r.to y = x line so that m>1 turns into m<=1 */
            this->start = Point(this->start.y, this->start.x);
            this->end = Point(this->end.y, this->end.x);
        }

        setColor(0, 0, 0);
    }

    int slopeTypeOfLine(Point start, Point end)
    {
        double dy = end.y - start.y;
        double dx = end.x - start.x;

        // cout << "m : "<< dy / dx << endl;
        this->slope = dy / dx;

        if (dx == 0)
            return VERTICAL_SLOPE;
        else if (dy == 0)
            return HORIZONTAL_SLOPE;
        else
        {
            double slope = dy / dx;
            if (slope <= 1.0 && slope > 0.0)
                return MIDPOINT_SLOPE; /* m<=1 */
            else if (slope > 1.0)
                return NON_MIDPOINT_SLOPE; /* m>1 */
            else if (slope >= -1.0)
                return NEGATIVE_SLOPE_1; /* m>=-1 */
            else if (slope < -1.0)
                return NEGATIVE_SLOPE_2; /* m<-1 */
        }
    }

    void colorPixel(int x, int y, double intensity = 1.0)
    {

        if (slopeType == NON_MIDPOINT_SLOPE)
            swap(x, y); /* swap back w.r.to y = x line */
        else if (slopeType == NEGATIVE_SLOPE_1)
            x = -x + 2 * parallelYlineDist; /* swap back w.r.to y = a line */
        else if (slopeType == NEGATIVE_SLOPE_2)
        {
            /* this time in reverse order of swapping first ; reflecting second */

            swap(x, y);                     /* swap back w.r.to y = x line */
            x = -x + 2 * parallelYlineDist; /* swap back w.r.to y = a line */
        }

        y = adjustPixelY(y);

        if (x < 0 || x >= pixelWidth || y < 0 || y >= pixelHeight)
            return;
        else
            image.set_pixel(x, y, color[0] * intensity, color[1] * intensity, color[2] * intensity);
    }

    void setColor(int r, int g, int b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    int adjustPixelY(int y)
    {
        return (pixelHeight - 1) - y;
    }

    double getIntensity(double distance)
    {
        // cout<<distance<<endl;

        // return  abs(distance);
        // return round(abs(distance));
        // return 1.0 - round(abs(distance));
        return 1.0 - abs(distance); /* intensity should be up if both centers are closer */
    }

    void drawVerticalLine()
    {

        int yEnd = end.y;
        int yStart = start.y;
        int x = start.x;

        for (int y = yStart; y <= yEnd; y++)
        {
            colorPixel(x, y);
        }
    }

    void drawHorizontalLine()
    {
        int xEnd = end.x;
        int xStart = start.x;
        int y = start.y;

        for (int x = xStart; x <= xEnd; x++)
        {
            colorPixel(x, y);
        }
    }

    void MidPointLineAlgo_R()
    {

        int x0 = start.x;
        int y0 = start.y;

        int x1 = end.x;
        int y1 = end.y;

        int dx = x1 - x0;
        int dy = y1 - y0;

        int d = 2 * dy - dx;
        int incE = 2 * dy;
        int incNE = 2 * (dy - dx);

        int x = x0;
        int y = y0;

        colorPixel(x, y);

        while (x < x1)
        {
            // cout<<x<<" "<<y<<endl;

            if (d < 0)
            {
                d += incE;
            }
            else
            {
                d += incNE;
                y++;
            }

            x++;
            colorPixel(x, y);
        }
    }

    void MidPointLineAlgo_RUA()
    {

        int x0 = start.x;
        int y0 = start.y;

        int x1 = end.x;
        int y1 = end.y;

        int dx = x1 - x0;
        int dy = y1 - y0;

        int d = 2 * dy - dx;
        int incE = 2 * dy;
        int incNE = 2 * (dy - dx);

        int x = x0;
        int y = y0;

        Rectangle rectangle(start, end, slope, 1); /* assumed thickness , t = 1.0 */

        colorPixel(x, y, rectangle.getIntensity(x, y));         /* current pixel */
        colorPixel(x, y + 1, rectangle.getIntensity(x, y + 1)); /* neighbor pixel above */
        colorPixel(x, y - 1, rectangle.getIntensity(x, y - 1)); /* neighbor pixel below */

        while (x < x1)
        {
            if (d < 0)
            {
                d += incE;
            }
            else
            {
                d += incNE;
                y++;
            }

            x++;

            colorPixel(x, y, rectangle.getIntensity(x, y));         /* current pixel */
            colorPixel(x, y + 1, rectangle.getIntensity(x, y + 1)); /* neighbor pixel above */
            colorPixel(x, y - 1, rectangle.getIntensity(x, y - 1)); /* neighbor pixel below */
        }
    }

    void MidPointLineAlgo_RWA()
    {

        int x0 = start.x;
        int y0 = start.y;

        int x1 = end.x;
        int y1 = end.y;

        int dx = x1 - x0;
        int dy = y1 - y0;

        int d = 2 * dy - dx;
        int incE = 2 * dy;
        int incNE = 2 * (dy - dx);

        int two_v_dx = 0;
        double BY_denom = 1.0 / (2.0 * sqrt(dx * dx + dy * dy));
        double two_dx_BY_denom = 2.0 * dx * BY_denom;

        int x = x0;
        int y = y0;

        colorPixel(x, y, getIntensity(0));                   /* current pixel */
        colorPixel(x, y + 1, getIntensity(two_dx_BY_denom)); /* neighbor pixel above */
        colorPixel(x, y - 1, getIntensity(two_dx_BY_denom)); /* neighbor pixel below */

        while (x < x1)
        {
            if (d < 0)
            {

                two_v_dx = d + dx;
                d += incE;
            }
            else
            {
                two_v_dx = d - dx;
                d += incNE;
                y++;
            }

            x++;
            colorPixel(x, y, getIntensity(two_v_dx * BY_denom));
            colorPixel(x, y + 1, getIntensity(two_dx_BY_denom - two_v_dx * BY_denom));
            colorPixel(x, y - 1, getIntensity(two_dx_BY_denom + two_v_dx * BY_denom));
        }
    }

    void print()
    {
        cout << "start = ";
        start.print();
        cout << "end   = ";
        end.print();
        cout << "color = (" << color[0] << " , " << color[1] << " , " << color[2] << ")" << endl;

        cout << "Type = ";
        switch (slopeType)
        {
        case MIDPOINT_SLOPE:
            cout << "MIDPOINT_SLOPE" << endl;
            break;
        case NON_MIDPOINT_SLOPE:
            cout << "NON_MIDPOINT_SLOPE" << endl;
            break;
        case HORIZONTAL_SLOPE:
            cout << "HORIZONTAL_SLOPE" << endl;
            break;
        case VERTICAL_SLOPE:
            cout << "VERTICAL_SLOPE" << endl;
            break;
        case NEGATIVE_SLOPE_1:
            cout << "NEGATIVE_SLOPE_1" << endl;
            break;
        case NEGATIVE_SLOPE_2:
            cout << "NEGATIVE_SLOPE_2" << endl;
            break;
        default:
            cout << "UNKNOWN" << endl;
            break;
        }

        cout << "--------------------------\n";
    }
};

void loadData()
{

    inputFile.open("input.txt");

    if (!inputFile.is_open())
    {
        cout << "failed to open input file" << endl;
        exit(EXIT_FAILURE);
    }

    inputFile >> pixelWidth >> pixelHeight;
    inputFile >> totalLines;

    cout << "pixelWidth = " << pixelWidth << endl;
    cout << "pixelHeight = " << pixelHeight << endl;
    cout << "totalLines = " << totalLines << endl;
    cout << "--------------------------\n";

    int xStart, yStart, xEnd, yEnd;
    int r, g, b;

    for (int i = 0; i < totalLines; i++)
    {
        inputFile >> xStart >> yStart >> xEnd >> yEnd;
        inputFile >> r >> g >> b;

        Line line(Point(xStart, yStart), Point(xEnd, yEnd));
        line.setColor(r, g, b);

        lineList.push_back(line);
    }

    inputFile.close();

    for (int i = 0; i < lineList.size(); i++)
        lineList.at(i).print();
}

void capture(string name)
{

    image.setwidth_height(pixelWidth, pixelHeight, true);

    for (int i = 0; i < lineList.size(); i++)
    {

        if (lineList.at(i).slopeType == VERTICAL_SLOPE)
        {
            lineList.at(i).drawVerticalLine();
        }
        else if (lineList.at(i).slopeType == HORIZONTAL_SLOPE)
        {
            lineList.at(i).drawHorizontalLine();
        }
        else
        {
            if (name == "1_R")
            {
                lineList.at(i).MidPointLineAlgo_R();
            }
            else if (name == "2_RUA")
            {
                lineList.at(i).MidPointLineAlgo_RUA();
            }
            else if (name == "3_RWA")
                lineList.at(i).MidPointLineAlgo_RWA();
        }
    };

    image.save_image(name + ".bmp");
    image.clear();
    cout << "\n";
    cout << name << ".bmp Image Saved..!\n\n";
}

int main()
{
    loadData();
    capture("1_R");
    capture("2_RUA");
    capture("3_RWA");

    return 0;
}
