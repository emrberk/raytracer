#include <iostream>
#include <climits>
#include "parser.h"
#include "ppm.h"
#include "utils.h"

typedef unsigned char RGB[3];
parser::Vec3f u;

typedef enum SurfaceType { TRIANGLE, SPHERE, MESH } SurfaceType;
parser::Scene scene;

typedef struct {
    float t;
    SurfaceType type;
    int id;
    int faceID;
} IntersectData;

parser::Vec3f findSphereNormal(parser::Sphere sphere, parser::Vec3f point) {
    parser::Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];
    parser::Vec3f normal = normalize(subtract(point, center));
    return normal;
}

parser::Vec3f findTriangleNormal(parser::Face face) {
    parser::Vec3f p1 = scene.vertex_data[face.v0_id - 1];
    parser::Vec3f p2 = scene.vertex_data[face.v1_id - 1];
    parser::Vec3f p3 = scene.vertex_data[face.v2_id - 1];
    parser::Vec3f edge1 = subtract(p2, p1);
    parser::Vec3f edge2 = subtract(p3, p1);
    return normalize(cross(edge1, edge2));
}

parser::Vec3f findIrradiance(parser::Vec3f intensity, parser::Vec3f lightDirection, parser::Vec3f intersectionPoint) {
    float lightDistance = length(lightDirection);
    parser::Vec3f irradiance = multiply(intensity, 1.0 / pow(lightDistance,2));
    return irradiance;
}

int clamp(int a) {
    if (a > 255) {
        return 255;
    } else if (a < 0) {
        return 0;
    }
    return a;
}

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
        IntersectData intersectData;
        parser::Vec3f diffuseComponent;

        for (int i = 0; i < scene.spheres.size(); i++) {
            parser::Sphere sphere = scene.spheres[i];
            float t = intersect(sphere);
            if (t >= 0 && t < minT) {
                minT = t;
                material = &scene.materials[sphere.material_id - 1];
                intersectData = { minT, SPHERE, i, -1 };
            }
        }
        for (int i = 0; i < scene.triangles.size(); i++) {
            parser::Triangle triangle = scene.triangles[i];
            float t = intersect(triangle.indices);
            if (t >= 0 && t < minT) {
                minT = t;
                material = &scene.materials[triangle.material_id - 1];
                intersectData = { minT, TRIANGLE, i, -1 };
            }
        }
        for (int i = 0; i < scene.meshes.size(); i++) {
            parser::Mesh mesh = scene.meshes[i];
            for (int j = 0; j < mesh.faces.size(); j++) {
                parser::Face face = mesh.faces[j];
                float t = intersect(face);
                if (t >= 0 && t < minT) {
                    minT = t;
                    material = &scene.materials[mesh.material_id - 1];
                    intersectData = { minT, MESH, i, j };
                }
            }
        }
        //printVec(material ? material->ambient : color);
        // ambient + diffuse + spectacular calculation
        if (!material)
            return color;
        
        for (auto& light : scene.point_lights) {
            parser::Vec3f intersectionPoint = add(origin, multiply(direction, minT));
            parser::Vec3f lightDirection = subtract(light.position, intersectionPoint);
            parser::Vec3f lightNormal = normalize(lightDirection);
            parser::Vec3f surfaceNormal;
            if (intersectData.type == SPHERE) {
                surfaceNormal = findSphereNormal(scene.spheres[intersectData.id], intersectionPoint);
            } else if (intersectData.type == TRIANGLE) {
                surfaceNormal = findTriangleNormal(scene.triangles[intersectData.id].indices);
            } else {
                parser::Face face = scene.meshes[intersectData.id].faces[intersectData.faceID];
                surfaceNormal = findTriangleNormal(face);
            }
            float cosTheta = dot(surfaceNormal,lightNormal);
            parser::Vec3f irradiance = findIrradiance(light.intensity,lightDirection,intersectionPoint);
            parser::Vec3f diffuse = material->diffuse;
            diffuseComponent = scale(multiplyTwo(multiply(irradiance, cosTheta), diffuse));
        }

        parser::Vec3f ambient = scale(material->ambient);
        color = add(ambient, diffuseComponent);
        return color;
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
                
                image[(y * camera.image_width + x) * 3] = clamp(int(color.x));
                image[(y * camera.image_width + x) * 3 + 1] = clamp(int(color.y));
                image[(y * camera.image_width + x) * 3 + 2] = clamp(int(color.z));
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
