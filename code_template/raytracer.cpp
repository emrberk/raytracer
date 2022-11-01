#include <iostream>
#include <climits>
#include "parser.h"
#include "ppm.h"
#include "utils.h"

typedef unsigned char RGB[3];
parser::Vec3f u;
parser::Scene scene;



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

    float intersect(parser::Face face) {
        parser::Vec3f a = scene.vertex_data[face.v0_id - 1];
        parser::Vec3f b = scene.vertex_data[face.v1_id - 1];
        parser::Vec3f c = scene.vertex_data[face.v2_id - 1];

        float A = determinant(a.x - b.x, a.x - c.x, direction.x, \
                    a.y - b.y, a.y - c.y, direction.y, \
                    a.z - b.z, a.z - c.z, direction.z);
    
        float betaD = determinant(a.x - origin.x, a.x - c.x, direction.x, \
                    a.y - origin.y, a.y - c.y, direction.y, \
                    a.z - origin.z, a.z - c.z, direction.z);
        float gammaD = determinant(a.x - b.x, a.x - origin.x, direction.x, \
                    a.y - b.y, a.y - origin.y, direction.y, \
                    a.z - b.z, a.z - origin.z, direction.z);
        float tD = determinant(a.x - b.x, a.x - c.x, a.x -origin.x, \
                    a.y - b.y, a.y - c.y, a.y - origin.y, \
                    a.z - b.z, a.z - c.z, a.z - origin.z);
        float beta = betaD / A;
        float gamma = gammaD / A;
        float t = tD / A;

        if (beta >= 0 && gamma >= 0 && beta + gamma <= 1) {
            return t;
        }

        return -1;
    }

    float intersect(parser::Mesh mesh) {
        float minT = INT_MAX;
        for (auto& face : mesh.faces) {
            float t = intersect(face);
            minT = std::min(t, minT);
        }
        return minT;
    }

    float intersect(parser::Sphere sphere) {
        float A, B, C;
        float delta;
        float result;
        parser::Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];
        //printVec(center);
        float radius = sphere.radius;

        C = pow(origin.x - center.x, 2) + pow(origin.y - center.y, 2) + pow(origin.z - center.z, 2) - pow(radius, 2);
        B = 2 * direction.x * (origin.x - center.x) + 2 * direction.y * (origin.y - center.y) + 2 * direction.z * (origin.z - center.z);
        A = pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2);

        delta = pow(B, 2) - 4 * A * C;

        if (delta < 0){
            return -1;
        } else if (delta == 0) {
            result = -B / (2 * A);
        } else {
            delta = sqrt(delta);
            float t1 = (-B + delta) / (2 * A);
            float t2 = (-B - delta)  / (2 * A);
            result = std::min(t1, t2);
        }
        return result;
    }

    parser::Vec3f computeColor() {
        float minT = INT_MAX;
        parser::Material* material = nullptr;
        parser::Vec3f color  = convert(scene.background_color);
        for (auto& sphere : scene.spheres) {
            float t = intersect(sphere);
            if (t >= 0 && t < minT) {
                minT = t;
                material = &scene.materials[sphere.material_id - 1];
            }
        }
        for (auto& triangle : scene.triangles) {
            float t = intersect(triangle.indices);
            if (t >= 0 && t < minT) {
                minT = t;
                material = &scene.materials[triangle.material_id - 1];
            }
        }
        for (auto& mesh : scene.meshes) {
            for (auto& face : mesh.faces) {
                float t = intersect(face);
                if (t >= 0 && t < minT) {
                    minT = t;
                    material = &scene.materials[mesh.material_id - 1];
                }
            }
        }
        // ambient + diffuse + spectacular calculation

        return material ? material->ambient : color;
    }
};

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file

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

    int i = 0;
    for (auto& camera : scene.cameras) {
        u = cross(camera.gaze, camera.up);
        unsigned char* image = new unsigned char [camera.image_width * camera.image_height * 3];
        for (int y = 0; y < camera.image_height; ++y) {
            for (int x = 0; x < camera.image_width; ++x) {
                /*
                int colIdx = x / columnWidth;
                image[i++] = BAR_COLOR[colIdx][0];
                image[i++] = BAR_COLOR[colIdx][1];
                image[i++] = BAR_COLOR[colIdx][2];
                */

                Ray* r = new Ray(camera, x, y);

                parser::Vec3f color = r->computeColor();
                
                image[(y * camera.image_width + x) * 3] = int(color.x * 255 + 0.5);
                image[(y * camera.image_width + x) * 3 + 1] = int(color.y * 255 + 0.5);;
                image[(y * camera.image_width + x) * 3 + 2] = int(color.z * 255 + 0.5);;
                delete r;
            }
        }
        std::string name = "test";
        name.append(std::to_string(i));
        name.append(".ppm");
        const char* c = name.c_str();
        write_ppm(c, image, camera.image_width, camera.image_height);
        i++;
    }
}
