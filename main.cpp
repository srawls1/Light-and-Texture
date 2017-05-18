#include<GL/glut.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<list>

using namespace std;

#define ImageW 400
#define ImageH 400

float framebuffer[ImageH][ImageW][3];
float zbuffer[ImageH][ImageW];
float ZMAX = 10000.0;	// NOTE: Assume no point has a Z value greater than 10000.0

char* sourcefile = "triangle.dat";

template<typename T>
inline T clamp(T value, T min, T max)
{
	if (value < min) return min;
	if (value > max) return max;
	return value;
}

template<typename T>
inline T min(T t1, T t2)
{
	return (t1 < t2) ? t1 : t2;
}

template<typename T>
inline T max(T t1, T t2)
{
	return (t1 < t2) ? t2 : t1;
}

struct Vector3
{
	float x, y, z;
	Vector3() {}
	Vector3(float xx, float yy, float zz)
		: x(xx), y(yy), z(zz) {}
};

Vector3 operator+(Vector3 v1, Vector3 v2)
{
	return Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vector3 operator*(float a, Vector3 v)
{
	return Vector3(a * v.x, a * v.y, a * v.z);
}

struct color {
	float r, g, b;		// Color (R,G,B values)
	color() {}
	color(float rr, float gg, float bb)
		: r(rr), g(gg), b(bb) {}
};

color operator*(color c1, color c2)
{
	return color(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b);
}

color operator*(float a, color c)
{
	return color(a * c.r, a * c.g, a * c.b);
}

color operator+(color c1, color c2)
{
	return color(c1.r + c2.r, c1.g + c2.g, c1.b + c2.b);
}

istream& operator>>(istream& in, color& c)
{
	return in >> c.r >> c.g >> c.b;
}

struct vertex {
	float x, y, z;		// x, y, z coordinates
	float nx, ny, nz;		// Normal at the vertex
	float u, v;			// Texture coordinates
};

struct triangle {
	int whichtexture;	// The index number of the corresponding texture to apply
						// Note: Use the color returned by the texture for the
						// ambient, diffuse, and specular color, scaled by the
						// coefficients of ambient, diffuse, and specular reflection
	vertex v[3];		// The three vertices
	float kamb;			// The coefficient of ambient reflection
	float kdiff;		// The coefficient of diffuse reflection
	float kspec;		// The coefficient of specular reflection
	int shininess;		// The exponent to use for Specular Phong Illumination
};

struct light {
	// Note: assume all lights are white
	float x, y, z;		// x, y, z coordinates of light
	color brightness;	// Level of brightness of light (0.0 - 1.0)
};

struct texture {
	// Note access using getTextureRGB provided below
	int xsize, ysize;	// The size of the texture in x and y
	float* elements;	// RGB values
};


int numtriangles;		// The number of triangles in the scene
int numlights;			// The number of lights (not including ambient) in the scene
int numtextures;		// The number of textures used in the scene

color ambientlight;		// The coefficient of ambient light

triangle* trianglelist;	// Array of triangles
light* lightlist;		// Array of lights
texture* texturelist;	// Array of textures

						/* Pass in a pointer to the texture, t, and the texture coordinates, u and v
						Returns (in R,G,B) the color of the texture at those coordinates */
void getTextureRGB(texture* t, float u, float v, float& R, float& G, float& B) {
	int xval, yval;
	if (u<1.0)
		if (u >= 0.0) xval = (int)(u*t->xsize);
		else xval = 0;
	else xval = t->xsize - 1;
	if (v<1.0)
		if (v >= 0.0) yval = (int)(v*t->ysize);
		else yval = 0;
	else yval = t->ysize - 1;

	R = t->elements[3 * (xval*t->ysize + yval)];
	G = t->elements[(3 * (xval*t->ysize + yval)) + 1];
	B = t->elements[(3 * (xval*t->ysize + yval)) + 2];
}


// Draws the scene
void drawit(void)
{
	glDrawPixels(ImageW, ImageH, GL_RGB, GL_FLOAT, framebuffer);
	glFlush();
}

void clearFramebuffer()
{
	for (int i = 0; i<ImageH; i++) {
		for (int j = 0; j<ImageW; j++) {
			framebuffer[j][i][0] = 0.0;
			framebuffer[j][i][1] = 0.0;
			framebuffer[j][i][2] = 0.0;
		}
	}
}

void clearZbuffer()
{
	for (int x = 0; x < ImageW; ++x)
	{
		for (int y = 0; y < ImageH; ++y)
		{
			zbuffer[y][x] = ZMAX;
		}
	}
}

// Sets pixel x,y to the color RGB
void setFramebuffer(int x, int y, color c)
{
	x = clamp(x, 0, ImageW);
	y = clamp(y, 0, ImageH);
	framebuffer[y][x][0] = clamp(c.r, 0.0f, 1.0f);
	framebuffer[y][x][1] = clamp(c.g, 0.0f, 1.0f);
	framebuffer[y][x][2] = clamp(c.b, 0.0f, 1.0f);
}

void setZbuffer(int x, int y, float z)
{
	x = clamp(x, 0, ImageW);
	y = clamp(y, 0, ImageH);
	zbuffer[y][x] = clamp<float>(z, 0, ZMAX);
}

float getZbuffer(int x, int y)
{
	if (x < 0 || x > ImageW || y < 0 || y > ImageH)
	{
		return 0.0f;
	}
	return zbuffer[y][x];
}

// Normalizes the Vector3 passed in
void normalize(Vector3& v) {
	float temp = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	if (temp > 0.0) {
		v.x /= temp;
		v.y /= temp;
		v.z /= temp;
	}
	else {
		v.x = 0.0;
		v.y = 0.0;
		v.z = 0.0;
	}
}

// Returns dot product of two Vector3s
float dot(Vector3 v1, Vector3 v2) {
	return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

// Returns angle between two Vector3s (in radians)
float angle(Vector3 v1, Vector3 v2) {
	normalize(v1);
	normalize(v2);
	return  acos(dot(v1, v2));
}

// Assume v1 and v2 are already normalized
Vector3 slerp(Vector3 v1, Vector3 v2, float t)
{
	float omega = angle(v1, v2);
	float sinomega = sin(omega);
	if (abs(sinomega) < 0.001)
		return v1;
	return sin((1 - t) * omega) / sinomega * v1 + sin(t * omega) / sinomega * v2;
}

struct activeEdge
{
	float x, z;
	float yoff;
	float dyoffdy;
	float dzdy, dxdy;
	float u, v;
	float dudy, dvdy;
	int miny, maxy;
	Vector3 normal;
	Vector3 lowerNorm, upperNorm;
	activeEdge(vertex* v1, vertex* v2) // pass in the lower y value first
	: miny(v1->y), maxy(v2->y),
	  x(v1->x), z(v1->z),
	  u(v1->u), v(v1->v), 
	  lowerNorm(v1->nx, v1->ny, v1->nz), upperNorm(v2->nx, v2->ny, v2->nz),
	  normal(v1->nx, v1->ny, v1->nz), yoff(0)
	{
		float ydif = v1->y - v2->y;
		if (ydif != 0)
		{
			dxdy = (v1->x - v2->x) / ydif;
			dzdy = (v1->z - v2->z) / ydif;
			dudy = (v1->u - v2->u) / ydif;
			dvdy = (v1->v - v2->v) / ydif;
			dyoffdy = 1.0f / (maxy - miny);
			ydif = miny - v1->y;
			if (ydif > 0)
			{
				x = v1->x + dxdy * ydif;
				z = v1->z + dzdy * ydif;
				yoff = dyoffdy * ydif;
				u = v1->u + dudy * ydif;
				v = v1->v + dvdy * ydif;
			}
		}
	}
};

void drawtriangle(triangle tri)
{
	texture* tex = texturelist + tri.whichtexture;

	int minY = min(min(tri.v[0].y, tri.v[1].y), tri.v[2].y);
	int maxY = max(max(tri.v[0].y, tri.v[1].y), tri.v[2].y);
	vector<activeEdge> edges;
	if (tri.v[0].y < tri.v[1].y) edges.push_back(activeEdge(&tri.v[0], &tri.v[1]));
	else edges.push_back(activeEdge(&tri.v[1], &tri.v[0]));
	if (tri.v[1].y < tri.v[2].y) edges.push_back(activeEdge(&tri.v[1], &tri.v[2]));
	else edges.push_back(activeEdge(&tri.v[2], &tri.v[1]));
	if (tri.v[2].y < tri.v[0].y) edges.push_back(activeEdge(&tri.v[2], &tri.v[0]));
	else edges.push_back(activeEdge(&tri.v[0], &tri.v[2]));

	// Create sorted edge table
	vector< list<activeEdge> > sortedEdges(maxY - minY + 1);
	for (activeEdge& l : edges)
	{
		sortedEdges[l.miny - minY].push_back(l);
	}

	// Initialize active edge list
	vector<activeEdge> activeEdges;

	// Start scan lines
	for (int y = minY; y < maxY; ++y)
	{
		// Update active edge list
		for (activeEdge l : sortedEdges[y - minY])
		{
			activeEdges.push_back(l);
		}
		for (auto iter = activeEdges.begin(); iter != activeEdges.end();)
		{
			activeEdge& l = *iter;
			if (y == l.maxy)
			{
				iter = activeEdges.erase(iter);
			}
			else
			{
				l.x += l.dxdy;
				l.u += l.dudy;
				l.v += l.dvdy;
				l.z += l.dzdy;
				l.yoff += l.dyoffdy;
				l.normal = slerp(l.lowerNorm, l.upperNorm, l.yoff);
				normalize(l.normal);
				++iter;
			}
		}

		activeEdge* leftEdge = (activeEdges[0].x < activeEdges[1].x) ? &activeEdges[0] : &activeEdges[1];
		activeEdge* rightEdge = (activeEdges[0].x < activeEdges[1].x) ? &activeEdges[1] : &activeEdges[0];

		float u = leftEdge->u;
		float v = leftEdge->v;
		float z = leftEdge->z;
		Vector3 normal = leftEdge->normal;
		float xoff = 0;
		float xdif = rightEdge->x - leftEdge->x;
		float dudx, dvdx, dzdx, dxoffdx;
		if (xdif > 0)
		{
			dudx = (rightEdge->u - leftEdge->u) / xdif;
			dvdx = (rightEdge->v - leftEdge->v) / xdif;
			dzdx = (rightEdge->z - leftEdge->z) / xdif;
			dxoffdx = 1 / xdif;
		}
		else
		{
			dudx = 0;
			dvdx = 0;
			dzdx = 0;
			dxoffdx = 0;
		}

		// Start Drawing
		for (int x = (int)leftEdge->x; x < rightEdge->x; ++x)
		{
			if (z < getZbuffer(x, y))
			{
				color c;
				getTextureRGB(tex, u, v, c.r, c.g, c.b);
				color cambient = tri.kamb * ambientlight;
				color cdiffuse(0, 0, 0);
				color cspecular(0, 0, 0);
				for (int i = 0; i < numlights; ++i)
				{
					light* l = lightlist + i;
					Vector3 lightdir(l->x - x, l->y - y, l->z - z);
					normalize(lightdir);
					float diffusefactor = clamp<float>(dot(lightdir, normal), 0, 1);
					Vector3 lightparallel = diffusefactor * normal;
					Vector3 lightperpendicular = lightdir + -1 * lightparallel;
					Vector3 reflectdir = lightparallel + -1 * lightperpendicular;
					float specularfactor = -reflectdir.z;
					cdiffuse = cdiffuse + diffusefactor * tri.kdiff * l->brightness;
					cspecular = cspecular + pow(specularfactor, tri.shininess) * tri.kspec * l->brightness;
				}
				color lightcol = cambient + cdiffuse + cspecular;
				setFramebuffer(x, y, c * lightcol);
				setZbuffer(x, y, z);
			}
			u += dudx;
			v += dvdx;
			z += dzdx;
			xoff += dxoffdx;
			normal = slerp(leftEdge->normal, rightEdge->normal, xoff);
			normalize(normal);
		}
	}
}

void display(void)
{
	clearFramebuffer();
	clearZbuffer();
	for (int i = 0; i < numtriangles; ++i)
	{
		triangle& t = trianglelist[i];
		//float nz = t.v[0].nz + t.v[1].nz + t.v[2].nz;
		//if (nz < 0)
		//{
			drawtriangle(t);
		//}
	}
	
	drawit();
}

void init(void)
{
	int i, j, k;

	// Initialize framebuffer to clear
	clearFramebuffer();
	clearZbuffer();

	// Load in data
	ifstream infile(sourcefile);
	if (!infile) {
		cout << "Error! Input file " << sourcefile << " does not exist!" << endl;
		exit(-1);
	}
	infile >> numtriangles >> numlights >> numtextures;

	// First read triangles
	trianglelist = new triangle[numtriangles];
	for (i = 0; i<numtriangles; i++) {
		infile >> trianglelist[i].whichtexture;
		infile >> trianglelist[i].kamb >> trianglelist[i].kdiff >> trianglelist[i].kspec;
		infile >> trianglelist[i].shininess;
		for (j = 0; j<3; j++) {
			infile >> trianglelist[i].v[j].x >> trianglelist[i].v[j].y >> trianglelist[i].v[j].z;
			infile >> trianglelist[i].v[j].nx >> trianglelist[i].v[j].ny >> trianglelist[i].v[j].nz;
			infile >> trianglelist[i].v[j].u >> trianglelist[i].v[j].v;
		}
	}

	// Now read lights
	lightlist = new light[numlights];
	infile >> ambientlight.r >> ambientlight.g >> ambientlight.b;
	for (i = 0; i<numlights; i++) {
		infile >> lightlist[i].x >> lightlist[i].y >> lightlist[i].z;
		infile >> lightlist[i].brightness.r >> lightlist[i].brightness.g >> lightlist[i].brightness.b;
	}

	// Now read textures
	texturelist = new texture[numtextures];
	for (i = 0; i<numtextures; i++) {
		infile >> texturelist[i].xsize >> texturelist[i].ysize;
		texturelist[i].elements = new float[texturelist[i].xsize*texturelist[i].ysize*3];
		for (j = 0; j<texturelist[i].xsize; j++) {
			for (k = 0; k<texturelist[i].ysize; k++) {
				infile >> texturelist[i].elements[3 * (j*texturelist[i].ysize + k)];
				infile >> texturelist[i].elements[3 * (j*texturelist[i].ysize + k) + 1];
				infile >> texturelist[i].elements[3 * (j*texturelist[i].ysize + k) + 2];
			}
		}
	}

	infile.close();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ImageW, ImageH);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Spencer Rawls - Assignment 5");
	init();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}