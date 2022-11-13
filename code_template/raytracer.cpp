#include <iostream>
#include <float.h>
#include <climits>
#include "parser.h"
#include "ppm.h"
#include "utils.h"

parser::Vec3f u;

typedef enum SurfaceType { TRIANGLE, SPHERE, MESH, NONE } SurfaceType;

parser::Scene scene;

typedef struct {
    float t;
    SurfaceType type;
    int id;
    parser::Material* material;
    int faceID;
} IntersectData;

parser::Vec3f findIrradiance(parser::Vec3f intensity, parser::Vec3f lightDirection) {
    float lightDistance = length(lightDirection);
    parser::Vec3f irradiance = divide(intensity, lightDistance * lightDistance);
    return irradiance;
}

class Ray {
private:
    parser::Vec3f origin;
    parser::Vec3f direction;

    float intersect(parser::Face face) {
        parser::Vec3f a = scene.vertex_data[face.v0_id - 1];
        parser::Vec3f b = scene.vertex_data[face.v1_id - 1];
        parser::Vec3f c = scene.vertex_data[face.v2_id - 1];

        float A = determinant(
            a.x - b.x, a.x - c.x, direction.x, \
            a.y - b.y, a.y - c.y, direction.y, \
            a.z - b.z, a.z - c.z, direction.z
        );
    
        float betaD = determinant(
            a.x - origin.x, a.x - c.x, direction.x, \
            a.y - origin.y, a.y - c.y, direction.y, \
            a.z - origin.z, a.z - c.z, direction.z
        );

        float gammaD = determinant(
            a.x - b.x, a.x - origin.x, direction.x, \
            a.y - b.y, a.y - origin.y, direction.y, \
            a.z - b.z, a.z - origin.z, direction.z
        );
        float tD = determinant(
            a.x - b.x, a.x - c.x, a.x -origin.x, \
            a.y - b.y, a.y - c.y, a.y - origin.y, \
            a.z - b.z, a.z - c.z, a.z - origin.z
        );

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

    parser::Vec3f findSphereNormal(parser::Sphere sphere, parser::Vec3f point) {
        parser::Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];
        parser::Vec3f normal = normalize(subtract(point, center));
        return normal;
    }

    parser::Vec3f findTriangleNormal(parser::Face face) {
        parser::Vec3f p1 = scene.vertex_data[face.v0_id - 1];
        parser::Vec3f p2 = scene.vertex_data[face.v1_id - 1];
        parser::Vec3f p3 = scene.vertex_data[face.v2_id - 1];
        parser::Vec3f edge1 = subtract(p3, p2);
        parser::Vec3f edge2 = subtract(p1, p2);
        parser::Vec3f normal = normalize(cross(edge1, edge2));

        return normal;
    }

    IntersectData intersectWithObjects() {
        float minT = FLT_MAX;
        parser::Material* material = nullptr;
        IntersectData intersectData;
        for (int i = 0; i < scene.spheres.size(); i++) {
            parser::Sphere sphere = scene.spheres[i];
            float t = intersect(sphere);
            if (t >= 0 && t < minT) {
                minT = t;
                material = &scene.materials[sphere.material_id - 1];
                intersectData = { minT, SPHERE, i, material, -1 };
            }
        }
        for (int i = 0; i < scene.triangles.size(); i++) {
            parser::Triangle triangle = scene.triangles[i];
            float t = intersect(triangle.indices);
            if (t >= 0 && t < minT) {
                minT = t;
                material = &scene.materials[triangle.material_id - 1];
                intersectData = { minT, TRIANGLE, i, material, -1 };
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
                    intersectData = { minT, MESH, i, material, j };
                }
            }
        }
        if (!material) {
            intersectData.type = NONE;
        }
        return intersectData;
    }

    parser::Vec3f getLightComponents(parser::PointLight light, parser::Vec3f intersectionPoint, IntersectData intersectData, parser::Vec3f surfaceNormal) {
        parser::Vec3f lightDirection = subtract(light.position, intersectionPoint);
        parser::Vec3f lightNormal = normalize(lightDirection);
        parser::Vec3f irradiance = findIrradiance(light.intensity,lightDirection);
        parser::Vec3f diffuseComponent = { 0, 0, 0 };
        parser::Vec3f specularComponent = { 0, 0, 0 };
        parser::Material* material = intersectData.material;

        // Diffuse
        float cosTheta = std::max(0.0f,dot(surfaceNormal,lightNormal));
        parser::Vec3f diffuse = material->diffuse;
        diffuseComponent = add(diffuseComponent, multiplyTwo(multiply(irradiance, cosTheta), diffuse));

        // Specular
        parser::Vec3f viewNormal = normalize(subtract(this->origin,intersectionPoint));
        parser::Vec3f halfVectorNormal = normalize(add(lightNormal,viewNormal));
        parser::Vec3f specular = material->specular;
        float cosAlpha = std::max(0.0f,dot(surfaceNormal,halfVectorNormal));
        float phongScalar = pow(cosAlpha, material->phong_exponent);
        specularComponent = add(specularComponent, multiplyTwo(multiply(irradiance, phongScalar), specular));

        return add(specularComponent, diffuseComponent);
    }

public:
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
        this->direction = normalize(subtract(rayPoint,camera.position));
    }

    Ray(parser::Vec3f origin, parser::Vec3f direction) {
        this->origin = origin;
        this->direction = direction;
    }

    parser::Vec3f computeColor(int recursionDepth) {
        parser::Vec3f color = convert(scene.background_color);
        IntersectData intersectData = intersectWithObjects();
        if (intersectData.type == NONE) {
            return color;
        }
        parser::Vec3f intersectionPoint = add(origin, multiply(direction, intersectData.t));
        parser::Material* material = intersectData.material;
        parser::Vec3f ambientComponent = multiplyTwo(material->ambient, scene.ambient_light);
        color = ambientComponent;
        parser::Vec3f lightComponents = { 0, 0, 0 };
        parser::Vec3f surfaceNormal;

        if (intersectData.type == SPHERE) {
            surfaceNormal = findSphereNormal(scene.spheres[intersectData.id], intersectionPoint);
        } else if (intersectData.type == TRIANGLE) {
            surfaceNormal = findTriangleNormal(scene.triangles[intersectData.id].indices);
        } else {
            parser::Face face = scene.meshes[intersectData.id].faces[intersectData.faceID];
            surfaceNormal = findTriangleNormal(face);
        }

        for (auto& light : scene.point_lights) {
            parser::Vec3f lightDirection = subtract(light.position, intersectionPoint);
            parser::Vec3f lightNormal = normalize(lightDirection);

            // The light will have no effect on the color in case of a shadow existence.
            parser::Vec3f error = multiply(surfaceNormal, scene.shadow_ray_epsilon);
            parser::Vec3f shadowOrigin = add(intersectionPoint, error);
            Ray shadowRay = Ray(shadowOrigin, lightNormal);
            IntersectData shadowIntersection = shadowRay.intersectWithObjects();

            if (shadowIntersection.type != NONE && shadowIntersection.t <= length(lightDirection)) {
                // Length of lightDirection gives t for the ray from point to light.
                continue;
            }

            lightComponents = add(lightComponents, getLightComponents(light, intersectionPoint, intersectData, surfaceNormal));
        }
        // Mirror reflections
        parser::Vec3f currentColor = add(color, add(ambientComponent, lightComponents));
        if (material->is_mirror && recursionDepth > 0) {
            float cosTheta = std::max(0.0f, dot(surfaceNormal, multiply(direction, -1)));
            parser::Vec3f newDirection = normalize(add(direction, multiply(multiply(surfaceNormal, cosTheta), 2)));
            parser::Vec3f error = multiply(newDirection, scene.shadow_ray_epsilon);
            parser::Vec3f newOrigin = add(intersectionPoint, error);
            parser::Vec3f mirrorCoefficient = material->mirror;
            this->direction = newDirection;
            this->origin = newOrigin;
            parser::Vec3f reflectedColor = this->computeColor(recursionDepth - 1);
            if (reflectedColor != convert(scene.background_color)) {
                parser::Vec3f mirrorComponent = multiplyTwo(reflectedColor, mirrorCoefficient);
                currentColor = add(currentColor, mirrorComponent);
            }
        }

        return clampVec(currentColor);
    }
};

int main(int argc, char* argv[]) {
    scene.loadFromXml(argv[1]);
    int recursionDepth = scene.max_recursion_depth;

    for (auto& camera : scene.cameras) {
        u = cross(camera.gaze, camera.up);
        unsigned char* image = new unsigned char [camera.image_width * camera.image_height * 3];
        for (int y = 0; y < camera.image_height; ++y) {
            for (int x = 0; x < camera.image_width; ++x) {
                Ray* r = new Ray(camera, x, y);

                parser::Vec3f color = r->computeColor(recursionDepth);
                
                image[(y * camera.image_width + x) * 3] = color.x;
                image[(y * camera.image_width + x) * 3 + 1] = color.y;
                image[(y * camera.image_width + x) * 3 + 2] = color.z;
                delete r;
            }
        }

        const char* imageName = camera.image_name.c_str();
        write_ppm(imageName, image, camera.image_width, camera.image_height);
    }
}
