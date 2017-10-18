#define UNICODE

// WinAPI
#include <Windows.h>
#include <tchar.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// OpenGL
#pragma comment(lib, "opengl32.lib")

// GLUT
#pragma comment(lib, "glut32.lib")
#include "glut.h"

// Other
#include <math.h>

using namespace std;

int xAngle = 0, yAngle = 0, zAngle = 0;

double Width = 0, Height = 0;

// Функции для работы с OpenGL/GLUT
void GLInit();
void GLRenderScene();
void GLKeyDown(unsigned char key, int x, int y);
void Reshape(GLint w, GLint h);
void glEnter2D(void);
void glLeave2D(void);
void glWrite(float x, float y, int *font, char text[256], int kls);

string location = "../МАГА/";

struct nvtr
{
	int uzel[8], numberField;

	nvtr(int number_field, int newNodes[8])
		: numberField(number_field)
	{
		memcpy(uzel, newNodes, 8);
	}

	nvtr()
	{
	}
};

struct point // точка
{
	double x, y, z;

	point(double x, double y, double z)
		: x(x),
		y(y),
		z(z)
	{
	}

	point() : x(0), y(0), z(0)
	{
	}

	friend bool operator==(const point& lhs, const point& rhs)
	{
		return lhs.x == rhs.x
			&& lhs.y == rhs.y
			&& lhs.z == rhs.z;
	}

	friend bool operator!=(const point& lhs, const point& rhs)
	{
		return !(lhs == rhs);
	}
};

vector<nvtr> KE;
vector<point> xyz_points;
int nc = 0;
double xMin, xMax, yMin, yMax, zMin, zMax;
double glxMin, glxMax, glyMin, glyMax, glzMin, glzMax;
double diag;
double zoom = 1, dZoom = 0.1;

void Input()
{
	int i, j, t, k;
	ifstream fileXY(location + "xyz.txt");
	ifstream fileNvtr(location + "nvtr.txt");

	//формируем файл nvtr.txt
	int n = 0;
	fileNvtr >> n;
	KE.resize(n);
	fileNvtr >> xMin >> xMax >> yMin >> yMax >> zMin >> zMax;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 8; j++)
		{
			fileNvtr >> KE[i].uzel[j];
		}
	}

	//формируем файл xyz.txt
	fileXY >> n >> nc;

	xyz_points.resize(n);
	for (i = 0; i < n; i++)
	{
		fileXY >> xyz_points[i].x >> xyz_points[i].y >> xyz_points[i].z;
	}

	diag = sqrt(
		(xMax - xMin)*(xMax - xMin) +
		(yMax - yMin)*(yMax - yMin) +
		(zMax - zMin)*(zMax - zMin))*1.2;
}

int WINAPI wWinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPWSTR lpCmdLine, int nShowCmd)
{
	// Аргументы командной строки (путь к EXE)
	char* argv0 = new char[512];
	GetModuleFileNameA(0, argv0, 512);

	int argc = 1;

	Input();

	// Инициализация GLUT
	glutInit(&argc, &argv0);

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_SINGLE | GLUT_RGB);

	// Координаты и размер окна GLUT
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(800, 600);

	// Создание окна
	glutCreateWindow("");

	// Обработчик рендеринга
	glutDisplayFunc(GLRenderScene);

	glutReshapeFunc(Reshape);

	// Обработчик клавитуры
	glutKeyboardFunc(GLKeyDown);

	// Инициализация OpenGL
	GLInit();

	// Главный цикл
	glutMainLoop();

	return 0;
}

void GLInit()
{
	// Цвет фона - черный
	glClearColor(1, 1, 1, 1);
}

void DrawHalfOfParallelepiped(float x0, float x1, float y0, float y1, float z0, float z1)
{
	glColor3f(0.0, 0.0, 0.0);
	// Передняя грань
	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, z0); // ЛН
	glVertex3f(x0, y1, z0); // ЛВ
	glVertex3f(x1, y1, z0); // ПВ
	glVertex3f(x1, y0, z0); // ПН
	glEnd();

	// Левая грань
	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, z0); // Перед-низ
	glVertex3f(x0, y1, z0); // Перед-верх
	glVertex3f(x0, y1, z1); // Зад-верх
	glVertex3f(x0, y0, z1); // Зад-низ
	glEnd();

	// Нижняя грань
	// Та же верхняя, только, Y = 0
	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, z0); // Перед-лево
	glVertex3f(x1, y0, z0); // Перед-право
	glVertex3f(x1, y0, z1); // Зад-право
	glVertex3f(x0, y0, z1); // Зад-лево
	glEnd();
}

void DrawParallelepiped(float x0, float x1, float y0, float y1, float z0, float z1)
{
	// Передняя грань
	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, z0); // ЛН
	glVertex3f(x0, y1, z0); // ЛВ
	glVertex3f(x1, y1, z0); // ПВ
	glVertex3f(x1, y0, z0); // ПН
	glEnd();

	// Задняя грань
	// Та же передняя, только, Z = 0.25, а, не 0
	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, z1); // ЛН
	glVertex3f(x0, y1, z1); // ЛВ
	glVertex3f(x1, y1, z1); // ПВ
	glVertex3f(x1, y0, z1); // ПН
	glEnd();

	// Левая грань
	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, z0); // Перед-низ
	glVertex3f(x0, y1, z0); // Перед-верх
	glVertex3f(x0, y1, z1); // Зад-верх
	glVertex3f(x0, y0, z1); // Зад-низ
	glEnd();

	// Правая грань
	// Та же левая, только, X = 0.5
	glBegin(GL_POLYGON);
	glVertex3f(x1, y0, z0); // Перед-низ
	glVertex3f(x1, y1, z0); // Перед-верх
	glVertex3f(x1, y1, z1); // Зад-верх
	glVertex3f(x1, y0, z1); // Зад-низ
	glEnd();

	// Верхняя грань
	glBegin(GL_POLYGON);
	glVertex3f(x0, y1, z0); // Перед-лево
	glVertex3f(x1, y1, z0); // Перед-право
	glVertex3f(x1, y1, z1); // Зад-право
	glVertex3f(x0, y1, z1); // Зад-лево
	glEnd();

	// Нижняя грань
	// Та же верхняя, только, Y = 0
	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, z0); // Перед-лево
	glVertex3f(x1, y0, z0); // Перед-право
	glVertex3f(x1, y0, z1); // Зад-право
	glVertex3f(x0, y0, z1); // Зад-лево
	glEnd();
}

void drawText(const char *text, int length, int x, int y)
{
	glMatrixMode(GL_PROJECTION);
	double *matrix = new double[16];
	glGetDoublev(GL_PROJECTION_MATRIX, matrix);
	glLoadIdentity();
	glOrtho(0, 600, 0, 600, -5, 5);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glLoadIdentity();
	glRasterPos2i(x, y);
	for (int i = 0; i<length; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, (int)text[i]);//GLUT_BITMAP_9_BY_15
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(matrix);
	glMatrixMode(GL_MODELVIEW);

}

void GLRenderScene()
{
	// Очищаем буферы
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(-1, -0.5, 0.75, 0, 0, 0, 0, 0, 1);
	//gluPerspective(50, 20, 1, 50);
	glTranslatef(-xMin, -yMin, -zMin);
	glScaled(zoom, zoom, zoom);
	//glOrtho(-3, 3, -3, 3, -3, 3);
	//glTranslatef(-1, 0, 0);
	//glScaled(1, 1, 1);

	// Включаем тест глубины,
	// чтобы грани 3D-фигур - рисовались непрозрачными, и не просвечивали
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// Поворот
	glRotatef(xAngle, 1, 0, 0);
	glRotatef(yAngle, 0, 1, 0);
	glRotatef(zAngle, 0, 0, 1);

	glEnter2D();
	glColor3f(0, 0, 0);
	glWrite(200, 200, (int*)GLUT_BITMAP_8_BY_13, (char*)"qwerty", 6);
	glLeave2D();

	glLineWidth(1.5); 
	glBegin(GL_LINE_STRIP);
	{
		glColor3f(255, 0, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(diag, 0, 0);

		glColor3f(0, 255, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, diag, 0);

		glColor3f(0, 0, 255);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, diag);
	}
	glEnd();

	//glPushMatrix();
	//char*text = new char[5];
	//sprintf(text, "%c", 'X');
	//drawText(text, strlen(text), 3, 0);
	//glPopMatrix();

	// -------- Рисование параллелепипеда ----------
	// Длина x высота x глубина:
	//   0.5    1.0      0.25
	glLineWidth(1);
	glColor3f(0, 0, 0);
	for (int i = 0; i < KE.size(); i++)
	{
		DrawParallelepiped(
			xyz_points[KE[i].uzel[0]].x, xyz_points[KE[i].uzel[1]].x,
			xyz_points[KE[i].uzel[0]].y, xyz_points[KE[i].uzel[2]].y,
			xyz_points[KE[i].uzel[0]].z, xyz_points[KE[i].uzel[4]].z);
	}
	//DrawParallelepiped(0, 0.5, 0, 1, 0, 0.25);

	//DrawHalfOfParallelepiped(-0.5, 0, -1, 0, -0.25, 0);

	glutSwapBuffers();
}

void GLKeyDown(unsigned char key, int x, int y)
{
	if (key == 'w' || key == 'W')
		xAngle += 5;
	if (key == 's' || key == 'S')
		xAngle -= 5;
	if (key == 'q' || key == 'Q')
		yAngle += 5;
	if (key == 'e' || key == 'E')
		yAngle -= 5;
	if (key == 'd' || key == 'D')
		zAngle += 5;
	if (key == 'a' || key == 'A')
		zAngle -= 5;
	if (key == '=' || key == '+')
		zoom += dZoom;
	if (key == '-')
		zoom -= dZoom;

	glutPostRedisplay(); // Перерисовываем окно
}

void Reshape(GLint w, GLint h) // При изменении размеров окна
{
	Width = w;
	Height = h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(0, w, 0, h, -15, 15);
	diag = 20;
	glOrtho(-diag, diag, -diag, diag, -diag, diag);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 20.0);
	//glMatrixMode(GL_MODELVIEW);
}

void glEnter2D(void) {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, Width, 0, Height);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glDisable(GL_DEPTH_TEST);
}

void glLeave2D(void) {
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glEnable(GL_DEPTH_TEST);
}

void glWrite(float x, float y, int *font, char text[256], int kls) {
	int i;
	glRasterPos2f(x, y);
	for (i = 0; i<kls; i++)
		glutBitmapCharacter(font, (int)text[i]);
}