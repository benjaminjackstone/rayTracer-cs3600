// Your Name Here!
// CS 3600
// Ray Tracer, crippled
// Modified in 2010 to detect MAC empty lines.
/*
parametric:
	p(t) = e + r*t // constants e come from eye and r comes from ray
	X = e_x + r_x*t 

implicit:
	F(x,y,z) = (x - Xc)^2 + (y - Yc)^2 + (z - Zc)^2 - R^2 = 0 //r is radius

substitution:
	take x,y,z from parametric and substitute inside the implicit
	F(x,y,z) = (e_x + r_x*T - Cx)^2 + (e_y + r_y*T - Cy)^2 + (e_z + r_z*T - Cz)^2 - R^2 = 0
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <time.h>
using namespace std;
#include "glut.h"
#include "graphics.h"
#include "vector.h"

// Turn this to zero to see OpenGL draw the correct picture.
// Turn it to one, and have your code ray-trace it the same way.
#define RAY_TRACE 1


class Sphere
{
public:
	Point3 origin;
	float radius;
	float ambient[4]; // ambient rgb1;
	float diffuse[4]; // diffuse rgb1;
	float specular[4]; // specular rgb1
	float specular_exponent;
};

class Scene
{
public:
	Scene();
	~Scene();
	bool Load(char * filename);
	void CreateRay(int x, int y, Vector3 & ray);
	void CastRay(const Vector3 & ray, float &r,float &g,float &b);
	void CalculateRightAndUp();

	// These values are loaded from the input file:
	int version; // version number of the file format.
	Point3 eye; // eye point
	Point3 at; // at point
	Vector3 up; // up vector
	float near_plane, far_plane; // near and far planes distances
	float fov; // field of view;
	int screen_x, screen_y; // screen size
	int num_point_lights;
	Point4 * pl;
	int num_spheres;
	Sphere * s;

	// These are calculated from the values read from the input file above:
	bool first_time;
	Vector3 forward;
	Vector3 right;
	Point3 ww_tr; // World Windows top right.
	Point3 ww_tl; // World Windows top left.
	Point3 ww_bl; // World Windows bottom left.
	Point3 ww_br; // World Windows bottom right.

	// These values store your ray traced pixels:
	GLubyte *pixels_red; // holds all pixels in the screen
	GLubyte *pixels_green;
	GLubyte *pixels_blue;
};

Scene::Scene()
{
	first_time = true;
	pl=NULL;
	s=NULL;
	pixels_red=NULL;
	pixels_green=NULL;
	pixels_blue=NULL;
}

Scene::~Scene()
{
	if(pl!=NULL)
		delete [] pl;
	if(s!=NULL)
		delete [] s;
	if(pixels_red!=NULL)
		delete [] pixels_red;
	if(pixels_green!=NULL)
		delete [] pixels_green;
	if(pixels_blue!=NULL)
		delete [] pixels_blue;
}
double radians(double d) {
	return d * 3.1416 / 180;
}
// These are the calculations we did in class.
void Scene::CalculateRightAndUp()
{
	first_time = false; // don't need to do this again.


						// Calculate Forward vector
	forward = at - eye;
	// Normalize Forward vector.
	forward.Normalize();
	// Right vector is Forward cross Up
	// Use cross product
	right = CrossProduct(forward, up);
	// Normalize Right vector.
	right.Normalize();
	// Reset Up vector to be Right cross Forward
	up = CrossProduct(right, forward);
	// Set World window corners to be a little in front of eye.
	// look up fov function
	// r is th distance up or down we take divide by screen x
	// u is the steps left or right, divide by screen y
	// tan(radian(fov)/2) = u/f and f = 1;
	double u = tan(radians(fov) / 2);
	double r = (double)screen_x * u / (double)screen_y;
	ww_tr = eye + 1 * forward + u*up + r*right;
	ww_tl = eye + 1 * forward + u*up - r*right;
	ww_br = eye + 1 * forward - u*up + r*right;
	ww_bl = eye + 1 * forward - u*up - r*right;

}

// These are the calculations we did in class.
void Scene::CreateRay(int x, int y, Vector3 &ray)
{
	if(first_time)
	{
		CalculateRightAndUp();
	}

	// Calculate right_ratio and up_ratio.
	double rightRatio = (double)x / (double)(screen_x - 1);
	double upRatio = (double)y / (double)(screen_y - 1);
	// Calculate where x,y on the screen maps to on the WorldWindow (film.)
	Point3 tempTop = AffineSum(ww_tl, ww_tr, rightRatio);
	Point3 tempBottom = AffineSum(ww_bl, ww_br, rightRatio);
	Point3 film = AffineSum(tempBottom, tempTop, upRatio);
	// Generate a ray from the eye to the film position.
	ray = film - eye;
	ray.Normalize();

}

// These are the calculations we did in class.
bool RayIntersectSphere(const Point3 &eye,
						const Vector3 & ray, 
						const Sphere & s,
						float & distance,
						Point3 &intersection_point,
						Vector3 &intersection_normal
						)
{
	Vector3 w = eye - s.origin;

	double A = DotProduct(ray, ray);
	double B = 2 * DotProduct(ray, w);
	double C = DotProduct(w, w) - s.radius*s.radius;

	//check if the sphere is hit
	double lineIntersects = B*B - 4 * A*C;
	if (lineIntersects < 0) {
		//sphere not hit
		return false;
	}
	distance = (-B - sqrt(lineIntersects)) / (2 * A);
	//check if sphere is behind eye
	if (distance < 0) {
		//sphere behind eye
		return false;
	}

	intersection_point = eye + ray * distance;
	intersection_normal = intersection_point - s.origin;
	intersection_normal.Normalize();
	return true;
}

// Some free code for you.
void Scene::CastRay(const Vector3 &ray,float &r,float &g,float &b)
{
	float closest_distance = far_plane; // far plane.
	int closest_hit_sphere = -1;
	Point3 closest_point;
	Vector3 closest_normal;
	for(int i=0; i<num_spheres; i++)
	{
		float distance;
		Point3 intersection_point;
		Vector3 intersection_normal;
		bool hit = RayIntersectSphere(	eye, ray, s[i], 
										distance,
										intersection_point,
										intersection_normal
										);
		if(hit && distance < closest_distance)
		{
			closest_hit_sphere = i;
			closest_distance = distance;
			closest_point = intersection_point;
			closest_normal = intersection_normal;
		}
	}

	// Add correct lighting. This is for Program #11.
		//reflection rule
		//find a vector from intersection of ray from eye to vector from origin to point at which the ray hits the sphere...normalize
		//for(each light source){ find a vector from intersection of ray from eye to vector from origin to point at which the ray hits the sphere...then find the vector to the light source. normalized. then find the difference between the two vectors. use the dot product and if its 1 then its a perfect contribution. all the way to -1 for 0 contribution.}
		//L = vector to the light
		//N = normal vector
		// d = Dot(N, L) if < 0, d = 0;
		/*for each light:
		find L;
		d = Dot(L, closest_normal);
			r += s[closest_hit_sphere].diffuse[0] * d; //clamp it non-negative.
			//specular light
		if d > 0:
		V = vector to the viewer
		*/

		// Add correct lighting. This is for Program #11.
		if (closest_hit_sphere != -1)
		{
			r = g = b = 0.0;

			// AMBIENT
			r += s[closest_hit_sphere].ambient[0];
			g += s[closest_hit_sphere].ambient[1];
			b += s[closest_hit_sphere].ambient[2];

			// DIFFUSE / SPECULAR
			/*
			-------------------------------
			Notes:
			Diffuse and specular are per light
			loop through lights
			//specular
			Specular is applied when the angle between normal and L from intersect is opposite of angle between eye and normal from intersect (aka, the ray)
			To check this, add vector eye to vector L. then normalize the product
			//shadows
			call RayIntersectSphere on vector from intersect to the light. Call on all sphere
			if it hits anything, see if the distance from
			--------------------------------
			intersect = where eye ray strike object
			normal = normal of point i on shape (i - center of shape)
			for each light as lightPoint{
			L = vector from intersect to light source (lightPoint - intersect);
			L.normalize();
			dot = DotProduct(L, normal);
			if(dot <= 0){
			//skip this
			continue;
			}
			for each sphere as sphere{
			hits = Rayintersectsphere(sphere, l);
			if(hits){
			see if sphere is closer than light;
			continue if it is;
			}
			}
			// s is array of spheres
			r += s[closest_hit_sphere].diffuse[0] * dot;
			g += s[closest_hit_sphere].diffuse[1] * dot;
			b += s[closest_hit_sphere].diffuse[2] * dot;
			// specular
			summed = L - ray;
			dot2 = DotProduct(normal of intersect, summed.normalize());
			if(dot2 > 0){
			//specular exponent is from input file
			contribution = pow(dot2, specularexponent);
			r+= s[closest_hit_sphere].specular[0] * contribution;
			g+= s[closest_hit_sphere].specular[1] * contribution;
			b+= s[closest_hit_sphere].specular[2] * contribution;
			}
			}
			if(r > 1) r = 1;
			if(g > 1) g = 1;
			if(b > 1) b = 1;
			*/
			for (int light = 0; light < num_point_lights; light++) {
				Vector3 lightVector = pl[light] - closest_point;
				lightVector.Normalize();

				float tempDistanceForShadow;
				Point3 tempIntersectionForShadow;
				Vector3 tempIntersectionNormal;
				bool blocked = false;
				for (int j = 0; j < num_spheres; j++) {
					bool temp = RayIntersectSphere(closest_point, lightVector, s[j], tempDistanceForShadow, tempIntersectionForShadow, tempIntersectionNormal);
					if (temp) {
						blocked = true;
						break;
					}
				}
				if (!blocked) {
					double dot = DotProduct(lightVector, closest_normal);
					if (dot > 0) {
						double diffuseLight = dot;
						r += diffuseLight*s[closest_hit_sphere].diffuse[0];
						g += diffuseLight*s[closest_hit_sphere].diffuse[1];
						b += diffuseLight*s[closest_hit_sphere].diffuse[2];
					}

					//specular
					//eye + light remember to negate the eye vector!!!
					Vector3 summed = lightVector + ray * -1;
					//normalize them
					summed.Normalize();
					//compare them against the normal
					float dot2 = DotProduct(closest_normal, summed);
					if (dot2 > 0) {
						//specular exponent is from input file
						float contribution = pow(dot2, 42);
						r += s[closest_hit_sphere].specular[0] * contribution;
						g += s[closest_hit_sphere].specular[1] * contribution;
						b += s[closest_hit_sphere].specular[2] * contribution;
					}
				}
			}
			if (r > 1) r = 1;
			if (g > 1) g = 1;
			if (b > 1) b = 1;
		}
		else
		{
			r = g = b = .0;
		}
	}

const int bufsize = 200;
bool ReadActualLine(ifstream & in, char buffer[])
{
	do
	{
		in.getline(buffer, bufsize);
		if(in.eof())
			return false;
	}
	while(strlen(buffer)==0 || buffer[0] == '/' || buffer[0] == '\n' || buffer[0] == '\r'); // ignore empty lines and comment lines.

	return true;
}

bool ReadInteger(ifstream & in, int & x)
{
	char buffer[bufsize];
	if(!ReadActualLine(in,buffer))
		return false;

	if(sscanf_s(buffer, "%i", &x) != 1)
		return false;

	return true;
}

bool ReadFloat(ifstream & in, float & x)
{
	char buffer[bufsize];
	if(!ReadActualLine(in,buffer))
		return false;

	if(sscanf_s(buffer, "%f", &x) != 1)
		return false;

	return true;
}

bool ReadPoint3(ifstream & in, Point3 & p)
{
	char buffer[bufsize];
	if(!ReadActualLine(in,buffer))
		return false;

	if(sscanf_s(buffer, "%f %f %f", &p.p[0], &p.p[1], &p.p[2]) != 3)
		return false;

	return true;
}

// Read 3 values, set last to 1.0;
bool ReadPoint4(ifstream & in, Point4 & p)
{
	char buffer[bufsize];
	if(!ReadActualLine(in,buffer))
		return false;

	if(sscanf_s(buffer, "%f %f %f", &p.p[0], &p.p[1], &p.p[2]) != 3)
		return false;

	p.p[3] = 1.0;

	return true;
}

bool ReadVector(ifstream & in, Vector3 & v)
{
	char buffer[bufsize];
	if(!ReadActualLine(in,buffer))
		return false;

	if(sscanf_s(buffer, "%f %f %f", &v.v[0], &v.v[1], &v.v[2]) != 3)
		return false;

	return true;
}

bool Scene::Load(char * filename)
{
	ifstream in(filename, ios::in);
	if(!in)
	{
		cerr << "Error. Could not open file " << filename << endl;
		return false;
	}

	//
	// Read in all data from non-empty, non-comment lines.
	//
	if(!ReadInteger(in, version))
		return false;

	if(!ReadPoint3(in, eye))
		return false;

	if(!ReadPoint3(in, at))
		return false;

	if(!ReadVector(in, up))
		return false;

	if(!ReadFloat(in, near_plane))
		return false;

	if(!ReadFloat(in, far_plane))
		return false;

	if(!ReadFloat(in, fov))
		return false;

	if(!ReadInteger(in, screen_x))
		return false;

	if(!ReadInteger(in, screen_y))
		return false;

	pixels_red = new GLubyte[screen_x*screen_y];
	if(!pixels_red)
		return false;
	pixels_green = new GLubyte[screen_x*screen_y];
	if(!pixels_green)
		return false;
	pixels_blue = new GLubyte[screen_x*screen_y];
	if(!pixels_blue)
		return false;


	// Read in all point lights.
	if(!ReadInteger(in, num_point_lights))
		return false;

	pl = new Point4[num_point_lights];
	int i;
	for(i=0; i<num_point_lights; i++)
	{
		if(!ReadPoint4(in, pl[i]))
			return false;
	}

	// Read in all Spheres.
	if(!ReadInteger(in, num_spheres))
		return false;

	s = new Sphere[num_spheres];
	for(i=0; i<num_spheres; i++)
	{
		if(!ReadPoint3(in, s[i].origin))
			return false;

		if(!ReadFloat(in, s[i].radius))
			return false;

		char buffer[bufsize];
		if(!ReadActualLine(in,buffer))
			return false;
		if(sscanf_s(buffer, "%f %f %f", &s[i].ambient[0], &s[i].ambient[1], &s[i].ambient[2]) != 3)
			return false;
		s[i].ambient[3] = 1;

		if(!ReadActualLine(in,buffer))
			return false;
		if(sscanf_s(buffer, "%f %f %f", &s[i].diffuse[0], &s[i].diffuse[1], &s[i].diffuse[2]) != 3)
			return false;
		s[i].diffuse[3] = 1;

		if(!ReadActualLine(in,buffer))
			return false;
		if(sscanf_s(buffer, "%f %f %f %f", &s[i].specular[0], &s[i].specular[1], &s[i].specular[2], &s[i].specular_exponent) != 4)
			return false;
		s[i].specular[3] = 1;
	}

	
	return true;
}


Scene the_scene;


void DrawSphere(float x, float y, float z, float radius)
{
	glPushMatrix();
	glTranslatef(x,y,z);
	int slices = 40;
	int stacks = 40;
	glutSolidSphere(radius, slices, stacks);
	glPopMatrix();
}

//
// GLUT callback functions
//

// This callback function gets called by the Glut
// system whenever it decides things need to be redrawn.
void display(void)
{
#if RAY_TRACE
	glClear(GL_COLOR_BUFFER_BIT);
	int offset = 0;
	glBegin(GL_POINTS);
	for(int y=0; y<the_scene.screen_y; y++)
	{
		for(int x=0; x<the_scene.screen_x; x++)
		{
			glColor3ub(the_scene.pixels_red[offset], the_scene.pixels_green[offset], the_scene.pixels_blue[offset]);
			glVertex2d(x+.5,y+.5);
			offset++;
		}
	}
	glEnd();
	glFlush();
#else
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();

	gluLookAt(	the_scene.eye.p[0], the_scene.eye.p[1], the_scene.eye.p[2],
				the_scene.at.p[0], the_scene.at.p[1], the_scene.at.p[2],
				the_scene.up.v[0], the_scene.up.v[1], the_scene.up.v[2]);

	for(int i=0; i<the_scene.num_spheres; i++)
	{
		glMaterialfv(GL_FRONT, GL_AMBIENT, the_scene.s[i].ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, the_scene.s[i].diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, the_scene.s[i].specular);
		glMaterialf(GL_FRONT, GL_SHININESS, the_scene.s[i].specular_exponent);
		Point3 &ori = the_scene.s[i].origin;
		DrawSphere(ori.p[0], ori.p[1], ori.p[2], the_scene.s[i].radius);
	}

	glutSwapBuffers();
#endif
}


// This callback function gets called by the Glut
// system whenever a key is pressed.
void keyboard(unsigned char c, int x, int y)
{
	switch (c) 
	{
		case 27: // escape character means to quit the program
			exit(0);
			break;
		default:
			return; // if we don't care, return without glutPostRedisplay()
	}

	glutPostRedisplay();
}


void SetTopView(int w, int h)
{
	// go into 2D mode
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Depending on the graphics card, either of the next two lines might work.
	gluOrtho2D(0,w,0,h);
	//gluOrtho2D(0,w-1,0,h-1);

	glMatrixMode(GL_MODELVIEW);
}

void SetPerspectiveView(int w, int h)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float aspectRatio = (GLdouble) w/(GLdouble) h;
	gluPerspective( 
	/* field of view in degree */ the_scene.fov,
	/* aspect ratio */ aspectRatio,
	/* Z near */ the_scene.near_plane, /* Z far */ the_scene.far_plane);
	glMatrixMode(GL_MODELVIEW);
}


// This callback function gets called by the Glut
// system whenever the window is resized by the user.
void reshape(int w, int h)
{
	// Set the pixel resolution of the final picture (Screen coordinates).
	glViewport(0, 0, w, h);

#if RAY_TRACE
	SetTopView(w,h);
#else
	SetPerspectiveView(w,h);
#endif

}


// Your initialization code goes here.
void PreInitializeMyStuff()
{
	if(!the_scene.Load("model.txt"))
	{
		exit(1);
	}

#if RAY_TRACE
	int offset = 0;

	Vector3 ray; // represents each ray.
	float r,g,b; // represents the color for each pixel

	for(int y=0; y<the_scene.screen_y; y++)
	{
		for(int x=0; x<the_scene.screen_x; x++)
		{
			the_scene.CreateRay(x, y, ray);
			the_scene.CastRay(ray,r,g,b);
			the_scene.pixels_red[offset] = GLubyte(r*255);
			the_scene.pixels_green[offset] = GLubyte(g*255);
			the_scene.pixels_blue[offset] = GLubyte(b*255);

			offset++;
		}
	}
#endif
}


// Your initialization code goes here.
void InitializeMyStuff()
{
#if RAY_TRACE
#else
	glEnable(GL_LIGHTING);	// enable general lighting

	GLfloat white_light[] = {.8,.8,.8,1};
	GLfloat low_light[] = {.2,.2,.2,1};
	// set light properties
	for(int i=0; i<the_scene.num_point_lights; i++)
	{
		glLightfv(GL_LIGHT0+i, GL_POSITION, the_scene.pl[i].p); // position light i
		glLightfv(GL_LIGHT0+i, GL_DIFFUSE, white_light); // specify light's color
		glLightfv(GL_LIGHT0+i, GL_SPECULAR, white_light);
		glEnable(GL_LIGHT0+i);	// enable the light i.
	}

	glEnable(GL_DEPTH_TEST); // turn on depth buffering
	glShadeModel(GL_SMOOTH);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
#endif
}

int main(int argc, char **argv)
{
	PreInitializeMyStuff();

	glutInit(&argc, argv);
#if RAY_TRACE
	glutInitDisplayMode(GLUT_RGB);
#else
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
#endif
	glutInitWindowSize(the_scene.screen_x, the_scene.screen_y);
	glutInitWindowPosition(0, 0);

	int fullscreen = 0;
	if (fullscreen) 
	{
		glutGameModeString("800x600:32");
		glutEnterGameMode();
	} 
	else 
	{
		glutCreateWindow("Shapes");
	}

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);

	glClearColor(1,1,1,0);	
	InitializeMyStuff();

	glutMainLoop();

	return 0;
}
