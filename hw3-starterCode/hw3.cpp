/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Your name here>
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glm/glm.hpp>
#include <math.h>


#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;
glm::vec3 topLeft;
glm::vec3 topRight;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
//#define WIDTH 128
// #define HEIGHT 96

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

float epsilon = 0.4f;


struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
bool checkInsideTriangle(glm::vec3 intersectPos, unsigned int triangleIndex);
float triangleIntersection(glm::vec3 rayOrigin, glm::vec3 rayDirection, unsigned int triangle_index, bool& does_intersect);


//holds the "location" in the real world of a pixel
std::vector<std::vector<glm::vec3> > pixelWorldVector;
std::vector<std::vector<glm::vec3> > finalColor;


void generatePixelWorld();

glm::vec3 shootRay(unsigned int pixelX, unsigned int pixelY, unsigned int lightIndex, bool& include_ambient);

//MODIFY THIS FUNCTION
void draw_scene()
{
  generatePixelWorld();
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      //TODO, run this for all lights, add contribution
      bool include_ambient = false;
      glm::vec3 color = shootRay(x,y, 0, include_ambient);
      //printf("color.y = %f\n", color.y);


      if(include_ambient){
        color.x += (ambient_light[0]);
        color.y += ambient_light[1];
        color.z += ambient_light[2];
      } else {
        //printf("no");
      }


      if(color.x > 1.0f){color.x = 1.0f;}
      if(color.y > 1.0f){color.y = 1.0f;}
      if(color.z > 1.0f){color.z = 1.0f;}


      unsigned char red = static_cast<unsigned char>(color.x * 255.0f);
      unsigned char green = static_cast<unsigned char>(color.y * 255.0f);
      unsigned char blue = static_cast<unsigned char>(color.z * 255.0f);




      plot_pixel(x, HEIGHT - y, red, green, blue);

    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}


void clearMyColorBuffer(){
  for(int i = 0; i<HEIGHT; i++){
    for(int j = 0; j<WIDTH; j++){
      finalColor[i][j].x = 0.0f;
      finalColor[i][j].y = 0.0f;
      finalColor[i][j].z = 0.0f;
    }
  }
}


void initVector(){

  finalColor.resize(HEIGHT);
  pixelWorldVector.resize(HEIGHT);
  for(int i = 0; i<HEIGHT; i++){
    pixelWorldVector[i].resize(WIDTH);
    finalColor[i].resize(WIDTH);
  }


}

void generatePixelWorld(){

  float widthF = WIDTH;
  float heightF = HEIGHT;
  float aspect_ratio = widthF/heightF;
  float FOV = 60.0f;
  FOV = FOV * 0.0174533; //convert to rad

  topLeft.x = -aspect_ratio * tan(FOV/2.0f);
  topLeft.y = std::tan(FOV/2.0f);
  topLeft.z = -1.0f;

  topRight.x = aspect_ratio * tan(FOV/2.0f);
  topRight.y = std::tan(FOV/2.0f);
  topRight.z = -1.0f;

  //now we push all the vectors, scanning left to right
  //and top to bottom

  printf("topRight.x = %f\n", topRight.x);
  printf("topRight.y = %f\n", topRight.y);

  float worldWidth = 2 * topRight.x;
  float worldHeight = 2 * topRight.y;

  glm::vec3 travel_pixel = topLeft;
  double rightHop = worldWidth/WIDTH;
  double bottomHop = worldHeight/HEIGHT;

  for(int j = 0; j<HEIGHT; j++){
    for(int i = 0; i<WIDTH; i++){
      /* Add the current scanning pixel */
      pixelWorldVector[j][i] = travel_pixel;
      travel_pixel.x += rightHop;

      //printf("x:%f y:%f z%f \n", travel_pixel.x, travel_pixel.y, travel_pixel.z);

    }
    //printf("travel.x = %f\n", travel_pixel.x);

    travel_pixel.x = topLeft.x; //reset x
    travel_pixel.y -= bottomHop; //go to lower row
  }

}

//given an origin and a ray, check if this ray intersects anything
bool testIntersection(glm::vec3 origin, glm::vec3 ray){
  glm::vec3 viewRay = glm::normalize(ray);

  float epsilon2 = 1.0f;; //greater for intersections

  //check for intersection with all circles 
  int index = -1;
  float t_min = FLT_MAX;
  //printf("max: %f\n", t_min);

  for(int i = 0; i<num_spheres; i++){
    float c_x = static_cast<float>(spheres[i].position[0]);
    float c_y = static_cast<float>(spheres[i].position[1]);
    float c_z = static_cast<float>(spheres[i].position[2]);

    //printf("cx:%f, cy:%f, cz:%f\n", c_x, c_y, c_z);

    float r = static_cast<float>(spheres[i].radius);

    float x0 = origin.x;
    float y0 = origin.y;
    float z0 = origin.z;

    float b = 2*(viewRay.x*(x0 - c_x) + viewRay.y*(y0-c_y) + viewRay.z*(z0 - c_z));
    float c = ((x0 - c_x) * (x0 - c_x)) + ((y0 - c_y) * (y0 - c_y)) + 
    ((z0 - c_z) * (z0 - c_z)) - (r*r);

    float discriminant = (b*b) - (4*c);
    if(discriminant < 0){
      continue; //doesn't intersect
    }

    float term = std::sqrtf(discriminant);

    float t0 = (-b + term)/2;
    float t1 = (-b - term)/2;


    if(t0 > 0.0f && t0 < t_min && t0 > epsilon){
      printf("sphere");
      t_min = t0;
      index = i;
    } 
    if(t1 > 0.0f && t1 < t_min && t1 > epsilon){
      printf("sphere");
      t_min = t1;
      index = i;
    }
  }

  for(int i = 0; i<num_triangles; i++){
    //find intersection with plane of triange
    bool dummy = false;
    float t = triangleIntersection(origin, ray, i, dummy);
    glm::vec3 intersect = origin + (viewRay * t);

    //check if intersect is inside triangle
    bool isInside = checkInsideTriangle(intersect, i);

    //only consider for update if distance is great enough
    glm::vec3 dist = (intersect - origin);
    // if(glm::length(dist) < 0.0000f){
    //   continue;
    // }

    if(isInside){
      if(t > 0.0f && t < t_min && t>epsilon){
        t_min = t;
        index = i;
      }
    }

  }

  if(index == -1){
    return false; //didn't intersect any circles or triangles
  } else {
    return true; //intersected a circle or triangle
  }
}

glm::vec3 getTriangleNormal(unsigned int triangleIndex){
  Triangle triangle = triangles[triangleIndex];
  glm::vec3 p0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  glm::vec3 a = p0;
  glm::vec3 b(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  glm::vec3 c(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
  glm::vec3 A = (b-a);
  glm::vec3 B = (c-a);
  glm::vec3 n = glm::cross(A,B);
  return n;
}

float triangleIntersection(glm::vec3 rayOrigin, glm::vec3 rayDirection, unsigned int triangle_index, bool& does_intersect){

   
    does_intersect = false;
    Triangle triangle = triangles[triangle_index];

    //position on the place
    glm::vec3 p0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
    
    glm::vec3 n = getTriangleNormal(triangle_index);

    //calculate d coeffecient
    //normal and point on a plane

    float d = -glm::dot(n, p0);
    glm::vec3 P1 = rayOrigin;
    glm::vec3 v = rayDirection;

    float n_dot_v = glm::dot(n, v);
    float n_dot_P1 = glm::dot(n, P1);


    float t = -(n_dot_P1 + d)/n_dot_v;
    does_intersect = true;
    return t;
}

float orientation(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3) 
{ 
    float val = (p2.y - p1.y) * (p3.x - p2.x) - 
              (p2.x - p1.x) * (p3.y - p2.y); 
  
    if (val == 0.0f) {printf("\n\n\n\n\n\n\nOHNO%f", 0.0f); return 1.0f;}  // colinear 
  
    return (val > 0.0f)? 1.0f: 1.0f; // clock or counterclock wise 
} 

float GetArea(glm::vec2 A, glm::vec2 B, glm::vec2 C){
  float sign = orientation(A, B, C);
  float area = sign*0.5f*((B.x - A.x)*(C.y-A.y) - (C.x - A.x)*(B.y-A.y));
  return area;
}


bool checkInsideTriangle(glm::vec3 intersectPos, unsigned int triangleIndex){
  Triangle triangle = triangles[triangleIndex];
  glm::vec3 C0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  glm::vec3 C1(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  glm::vec3 C2(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
  //determine how we will collapse the triangle
  glm::vec3 xy(0.0f, 0.0f, 1.0f);
  glm::vec3 yz(1.0f, 0.0f, 0.0f);
  glm::vec3 xz(0.0f, 1.0f, 0.0f);

  glm::vec3 n = getTriangleNormal(triangleIndex); 
  unsigned int planeEnum = 0; //0 = xy, 1 = yz, 2 = xz
  glm::vec3 planeNormal = xy;
  if(fabs(n.y) > fabs(n.x) && fabs(n.y) > fabs(n.z)){
    planeNormal = xz;
    planeEnum = 2;
  } else if(fabs(n.x) > fabs(n.y) && fabs(n.x) > fabs(n.z)){
    planeNormal = yz;
    planeEnum = 1;
  }

  glm::vec2 c0;
  glm::vec2 c1; 
  glm::vec2 c2;
  glm::vec2 intersect;
  //now collapse
  if(planeEnum == 0){ //xy
    c0 = glm::vec2(C0.x, C0.y);
    c1 = glm::vec2(C1.x, C1.y);
    c2 = glm::vec2(C2.x, C2.y);
    intersect = glm::vec2(intersectPos.x, intersectPos.y);
  } else if(planeEnum == 1){ //yz
    c0 = glm::vec2(C0.y, C0.z);
    c1 = glm::vec2(C1.y, C1.z);
    c2 = glm::vec2(C2.y, C2.z);
    intersect = glm::vec2(intersectPos.y, intersectPos.z);
  } else { //xz
    c0 = glm::vec2(C0.x, C0.z);
    c1 = glm::vec2(C1.x, C1.z);
    c2 = glm::vec2(C2.x, C2.z);
    intersect = glm::vec2(intersectPos.x, intersectPos.z);
  }

  float denominator = GetArea(c0,c1,c2);
  //now in anti-clockwise direction, get areas
  float alpha = GetArea(intersect, c1, c2)/denominator;
  float beta = GetArea(c0, intersect, c2)/denominator;
  float gamma = 1.0f - alpha - beta;

  if(alpha < 0.0f || alpha > 1.0f){
    return false;
  }
  if(beta < 0.0f || beta > 1.0f){
    return false;
  }
  if(gamma < 0.0f || gamma > 1.0f){
    return false;
  }

  if(alpha + beta + gamma < 0.99f || alpha + beta + gamma > 1.01f){
    return false;
  }

  return true;






  //choose the one closest to it, i.e. where the dot product is zero, i.e.

}

/* Shoot a ray from the camera location, to the pixel location, and return the 
color for that pixel */
glm::vec3 shootRay(unsigned int pixelX, unsigned int pixelY, unsigned int lightIndex, bool& include_ambient){
  glm::vec3 cameraLocation(0.0f,0.0f,0.0f);
  glm::vec3 pixelLocation = pixelWorldVector[pixelY][pixelX];

  glm::vec3 ray = pixelLocation - cameraLocation;
  glm::vec3 viewRay = glm::normalize(ray);

  //check for intersection with all circles 
  int index = -1;
  float t_min = FLT_MAX;
  int type = 0; //0 = none, 1 = circle, 2 = triangle
  //printf("max: %f\n", t_min);

  for(int i = 0; i<num_spheres; i++){
    float c_x = static_cast<float>(spheres[i].position[0]);
    float c_y = static_cast<float>(spheres[i].position[1]);
    float c_z = static_cast<float>(spheres[i].position[2]);

    //printf("cx:%f, cy:%f, cz:%f\n", c_x, c_y, c_z);
    float r = static_cast<float>(spheres[i].radius);
    float x0 = 0.0f;
    float y0 = 0.0f;
    float z0 = 0.0f;

    float b = 2*(viewRay.x*(x0 - c_x) + viewRay.y*(y0-c_y) + viewRay.z*(z0 - c_z));
    float c = ((x0 - c_x) * (x0 - c_x)) + ((y0 - c_y) * (y0 - c_y)) + 
    ((z0 - c_z) * (z0 - c_z)) - (r*r);

    float discriminant = (b*b) - (4*c);
    if(discriminant < 0){
      continue; //doesn't intersect
    }

    float term = std::sqrtf(discriminant);

    float t0 = (-b + term)/2;
    float t1 = (-b - term)/2;

    if(t0 > 0 && t0 < t_min && t0 > epsilon){
      t_min = t0;
      index = i;
      type = 1;
    } 
    if(t1 > 0 && t1 < t_min && t1 > epsilon){
      t_min = t1;
      index = i;
      type = 1;
    }
  }

  //TODO: Also check for triangles, and update t_min, type, and index
  glm::vec3 triangleIntersectPos;
  for(int i = 0; i<num_triangles; i++){
            
    bool does_intersect = false;
    //d also t?
    float t = triangleIntersection(cameraLocation, viewRay, i, does_intersect);
    if(does_intersect == false){
      continue; //nothing
    }

    if( t > t_min || t < epsilon){
      continue;
    }
    //we have a position of intersection on the plane
    glm::vec3 intersec = viewRay * t;

    //now determine if the intersec was actually in the triangle
    bool inTriangle = checkInsideTriangle(intersec, i);

    if(inTriangle){
      //now check if this can update current t_min
      if(t > 0.0f && t < t_min && t > epsilon){
        t_min = t;
        type = 2;
        index = i;
      }
    }
    //now check if this plane hit is also inside triangle, if it is, test shadow/PHONG
  }

  //after this if still index == -1, intersects nothing
  //else send shadow ray from pos
  if(index == -1){
    include_ambient = false;

    return glm::vec3(1,1,1);
  } else {
    include_ambient = true;
    //calculate the value from lightIndex

    Light currLight = lights[lightIndex];
    glm::vec3 intersectPos = t_min * viewRay;
    glm::vec3 lightPos(currLight.position[0], 
                       currLight.position[1],
                       currLight.position[2]);

    glm::vec3 lightColor(currLight.color[0],
                         currLight.color[1],
                         currLight.color[2]);

    glm::vec3 shadowRay = glm::normalize(lightPos - intersectPos);
    glm::vec3 NORMAL; //init based on circle or triangle
    glm::vec3 kd; //diffuse
    glm::vec3 ks; //specular
    float shiny;

    //now check if this ray intersects anything
    bool intersectTest = testIntersection(intersectPos, shadowRay);

    if(intersectTest == true){
      return glm::vec3(0,0,0);
    } else {

      if(type == 1){
        Sphere curr_sphere = spheres[index];
        float circle_radius = curr_sphere.radius;
        kd.x = curr_sphere.color_diffuse[0];
        kd.y = curr_sphere.color_diffuse[1];
        kd.z = curr_sphere.color_diffuse[2];

        ks.x = curr_sphere.color_specular[0];
        ks.y = curr_sphere.color_specular[1];
        ks.z = curr_sphere.color_specular[2];

        shiny = curr_sphere.shininess;

        //calculate the normal
        NORMAL.x = (intersectPos.x - curr_sphere.position[0]);
        NORMAL.y = (intersectPos.y - curr_sphere.position[1]);
        NORMAL.z = (intersectPos.z - curr_sphere.position[2]);
        NORMAL = glm::normalize(NORMAL);

        //the reflected vector
        glm::vec3 reflected = glm::normalize(glm::reflect(shadowRay, NORMAL));

        //now with this info calc color cont
        float l_dot_n = glm::dot(shadowRay, NORMAL);

        //TODO: FIGURE OUT WHY I NEED TO FLIP THIS
        glm::vec3 v = (-viewRay);
        float r_dot_v = glm::dot(-reflected, v);

        if(l_dot_n < 0.0f){l_dot_n = 0.0f;}
        if(r_dot_v < 0.0f){r_dot_v = 0.0f;}

        float exp_term = std::powf(r_dot_v, shiny);

        glm::vec3 resultColor = lightColor * (kd * l_dot_n) + (ks * exp_term);
        return resultColor;


      } else if(type == 2){
        //interpolate the colors
        return glm::vec3(0,1,0);
      }


    }
  }


}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  initVector();

  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

