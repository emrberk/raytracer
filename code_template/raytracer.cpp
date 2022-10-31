#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "utils.h"

typedef unsigned char RGB[3];
parser::Vec3f u;



class Ray {
public:
    parser::Vec3f origin;
    parser::Vec3f direction;

    Ray(parser::Camera camera, int i, int j) {
        // Assuming gaze is normalized
        float left = camera.near_plane.x;
        float right = camera.near_plane.y;
        float bottom = camera.near_plane.z;
        float top = camera.near_plane.w;
        float su = (i + 0.5) * (right - left) / camera.image_width;
        float sv = (j + 0.5) * (top - bottom) / camera.image_height;
        parser::Vec3f m = add(multiply(camera.gaze, camera.near_distance), camera.position);
        parser::Vec3f start = add(m, add(multiply(u,left),multiply(camera.up,top)));
        parser::Vec3f rayPoint = add(start,add(multiply(u,su),multiply(camera.up,-sv)));
        this->origin = camera.position;
        this->direction = add(rayPoint,multiply(camera.position,-1));
    }

};

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };
    /*
    for(auto& camera : scene.cameras) {
        u = cross(camera.gaze, camera.up);
    }
    */

    parser::Camera first = scene.cameras[0];
    u = cross(first.gaze, first.up);
    Ray* r1 = new Ray(first, 0, 0);
    Ray r = *r1;
    std::cout << r.origin.x << " " << r.origin.y << " " << r.origin.z << std::endl;
    std::cout << r.direction.x << " " << r.direction.y << " " << r.direction.z << std::endl;
    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];
    /*
    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }
    

    write_ppm("test.ppm", image, width, height);
   */
}
