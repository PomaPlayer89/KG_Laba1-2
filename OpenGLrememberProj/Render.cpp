#include "Render.h"

#include <sstream>
#include <iostream>

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include "GUItextRectangle.h"
#include <string>

//-----------------------------------------
#include<vector>
using namespace std;
//-----------------------------------------
bool textureMode = true;
bool lightMode = true;

//-----------------------------------------
//активациа alpha наложения
int alpha = 0;
//смена текстуры
bool texChange = true;
//-----------------------------------------

//класс для настройки камеры
class CustomCamera : public Camera
{
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	
	//значния масеры по умолчанию
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera()
	{
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //создаем объект камеры


//Класс для настройки света
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	
	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //создаем источник света




//старые координаты мыши
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой левой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}
	//----------------------------------------------
	//активация альфа наложения
	if (key == 'W') {
		alpha++;
		if (alpha == 3) {
			alpha = 0;
		}
	}
	//смена текстуры
	if (key == 'Q') {
		texChange = !texChange;
	}
	//-----------------------------------------------
}

void keyUpEvent(OpenGL *ogl, int key)
{
	
}



vector<string> nametex_str = { "Nat.bmp", "picture.bmp"};
GLuint texId[2];

//выполняется перед первым рендером
void initRender(OpenGL *ogl)
{
	//настройка текстур

	//4 байта на хранение пикселя
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//настройка режима наложения текстур
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//включаем текстуры
	glEnable(GL_TEXTURE_2D);
	

	//массив трехбайтных элементов  (R G B)
	RGBTRIPLE *texarray;

	//массив символов, (высота*ширина*4      4, потомучто   выше, мы указали использовать по 4 байта на пиксель текстуры - R G B A)
	char *texCharArray;
	int texW, texH;


	for (int i = 0; i < nametex_str.size(); i++) {
		const char* nametex = nametex_str[i].c_str();
		OpenGL::LoadBMP(nametex, &texW, &texH, &texarray);
		OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

		//генерируем ИД для текстуры
		glGenTextures(1, &texId[i]);
		//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
		glBindTexture(GL_TEXTURE_2D, texId[i]);

		//загружаем текстуру в видеопямять, в оперативке нам больше  она не нужна
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

		//отчистка памяти
		free(texCharArray);
		free(texarray);

		//наводим шмон
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	}


	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH); 


	//   задать параметры освещения
	//  параметр GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  лицевые и изнаночные рисуются одинаково(по умолчанию), 
	//                1 - лицевые и изнаночные обрабатываются разными режимами       
	//                соответственно лицевым и изнаночным свойствам материалов.    
	//  параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение, 
	//                не зависящее от сточников
	// по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}

//-------------------------------------------------------------------------------
class Point {
public:
	double x;
	double y;
	double z;
	Point(double A[]) {
		this->x = A[0];
		this->y = A[1];
		this->z = A[2];
	}
	Point(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	void DrawPoint() {
		glVertex3d(this->x, this->y, this->z);
	}
	void Normal3d(int i = 1) {
		if (i > 0) {
			glNormal3d(x, y, z);
		}
		else {
			//меняем направление нормали
			x *= -1;
			y *= -1;
			z *= -1;
			glNormal3d(x, y, z);
		}
	}
};
class PointXY {
public:
	double x;
	double y;
	PointXY(double x, double y) {
		this->x = x;
		this->y = y;
	}
	void TexCoord2d() {
		glTexCoord2d(x, y);
	}
};
//найти вектор (A - начало вектора, B - конец вектора)
Point SearchVector(Point A, Point B) {
	return Point(B.x - A.x, B.y - A.y, B.z - A.z);
}
//векторное произведение
Point VectorProduct(Point vectorA, Point vectorB) {
	Point result(0, 0, 0);
	result.x = vectorA.y * vectorB.z - vectorB.y * vectorA.z;
	result.y = -1 * vectorA.x * vectorB.z + vectorB.x * vectorA.z;
	result.z = vectorA.x * vectorB.y - vectorB.x * vectorA.y;
	return result;
}
//поиск длины вектора
double SearchVectorLength(Point vector) {
	double length = sqrt(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
	return length;
}
//ищем нормаль к треугольнику и нормализируем ее
 Point SearchNormal(double A1[], double B1[], double C1[], bool flag = true) {

	 Point A(A1), B(B1), C(C1);

	//поиск векторов по заданнымточкам
	Point vector_b = SearchVector(B, C);
	Point vector_a = SearchVector(B, A);

	//векторное произведение
	Point NormalVector = VectorProduct(vector_a, vector_b);

	//поиск длины вектора
	double length = SearchVectorLength(NormalVector);

	if (flag) {
		//нормализация вектора
		NormalVector.x = NormalVector.x / length;
		NormalVector.y = NormalVector.y / length;
		NormalVector.z = NormalVector.z / length;
	}
	return NormalVector;
}
 //находим расстояние между точками
 static double SearchDistancePoints(Point A, Point B) {
	 double X = B.x - A.x;
	 double Y = B.y - A.y;
	 double Z = B.z - A.z;
	 X = pow(X, 2);
	 Y = pow(Y, 2);
	 Z = pow(Z, 2);
	 double r = sqrt(X + Y + Z);
	 return r;
 }
 PointXY UpdatePoint(double A1[], bool set = false, vector<Point> points = { Point(0, 0, 0) }) {
	 Point A(A1);
	 static bool installation = false;
	 static double width = 0;
	 static double height = 0;
	 static double min_x, min_y, max_x, max_y;
	 PointXY newCoord(2, 2);
	 static PointXY new_O(0, 0);
	 if (set == true) {
		 //определяем границы фигуры без учета выпуклости
		 min_x = points[0].x, min_y = points[0].y, max_x = points[0].x, max_y = points[0].y;
		 for (int i = 0; i < points.size(); i++) {
			 if (min_x > points[i].x) {
				 min_x = points[i].x;
			 }
			 if (max_x < points[i].x) {
				 max_x = points[i].x;
			 }
			 if (min_y > points[i].y) {
				 min_y = points[i].y;
			 }
			 if (max_y < points[i].y) {
				 max_y = points[i].y;
			 }
		 }
		 //находим радиус выпуклости
		 double r = SearchDistancePoints(points[1], points[2]) / 2;

		 //расчет длины и ширины прямоугольника куда вписана наша фигура
		 width = abs(min_x) + abs(max_x) + r;
		 height = abs(min_y) + abs(max_y) + r;


		 new_O.x = A.x;
		 new_O.y = A.y;

		 //настройка произведена
		 installation = true;

		 return newCoord;
	 }
	 if (installation) {
		 //расчет новых координат относительно 
		 newCoord.x = (A.x - max_x) / width;
		 newCoord.y = (A.y - max_y) / height;
		 return newCoord;
	 }
	 newCoord.x = 9;
	 newCoord.y = 9;
	 return newCoord;
 }
//-------------------------------------------------------------------------------

void Cic(double z, int i = 1)
{
	double x = 2.5, y = -2, r = 3.64005;
	double A[] = { x, y, z };
	for (double t = 0.27; t < 3.42 - 0.01; t += 0.01)
	{
		double O[] = { x - cos(t) * r, y - sin(t) * r, z };
		double O1[] = { x - cos(t+0.01) * r, y - sin(t+0.01) * r, z };
		SearchNormal(A, O, O1).Normal3d(i);
		UpdatePoint(A).TexCoord2d();
		glVertex3dv(A);
		UpdatePoint(O).TexCoord2d();
		glVertex3dv(O);
		UpdatePoint(O1).TexCoord2d();
		glVertex3dv(O1);
	}
}

void Lateral()
{
	double x = 2.5, y = -2, r = 3.64005;
	//------------------------------------------------
	double number = 1;
	for (double t = 0.27; t < 3.41; t += 0.001) {
		number++;
	}
	//------------------------------------------------
	double i = 0;
	for (double t = 0.27; t < 3.41; t += 0.001)
	{
		double x1 = x - cos(t) * r;
		double y1 = y - sin(t) * r;
		double x2 = x - cos(t + 0.01) * r;
		double y2 = y - sin(t + 0.01) * r;

		double A[] = { x1, y1, 1 };
		double AA[] = { x1, y1, 3 };
		double B[] = { x2, y2, 1 };

		SearchNormal(AA, A, B).Normal3d(-1);
		
		glTexCoord2d(i / number, 0);
 		glVertex3d(x1, y1, 1);

		glTexCoord2d(i / number, 0);
		glVertex3d(x2, y2, 1);

		glTexCoord2d((i + 1.0) / number, 1);
		glVertex3d(x2, y2, 3);

		glTexCoord2d((i + 1.0) / number, 1);
		glVertex3d(x1, y1, 3);

		i++;
	}
}

void Vpyklost()
{
	double X = -4, Y = 2, R = 4.47;
	double number = 1;
	for (double t = 2.675; t <= 4.241; t += 0.001) {
		number++;
	}
	double i = 0;
	for (double t = 2.675; t <= 4.241; t += 0.001)
	{

		double x1 = X - cos(t) * R;
		double x2 = X - cos(t + 0.01) * R;
		double y1 = Y - sin(t) * R;
		double y2 = Y - sin(t + 0.01) * R;

		double A[] = { x1, y1, 1 };
		double AA[] = { x1, y1, 3 };
		double B[] = { x2, y2, 1 };

		SearchNormal(AA, A, B).Normal3d();
		glTexCoord2d(i / number, 0);
		glVertex3d(x1, y1, 1);
		glTexCoord2d(i / number, 0);
		glVertex3d(x2, y2, 1);
		glTexCoord2d((i + 1.0) / number, 1);
		glVertex3d(x2, y2, 3);
		glTexCoord2d((i + 1.0) / number, 1);
		glVertex3d(x1, y1, 3);

		i++;
	}

}

//вверх призмы
void Up() 
{
	//------------------------------------------------------------
	//подключение альфа наложения
	switch (alpha) {
	case 1:
		//включаем режим смешивания
		glEnable(GL_BLEND);
		//задаем опцию для коэффициентов источника и приемника
		glBlendFunc(GL_ONE, GL_ONE);
		break;
	case 2:
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glColor4d(0.7, 0.1, 0.1, 0.6);
		break;
	}
	//------------------------------------------------------------
	double K[] = { 0,0,3 }; 
	double A[] = { -4,0,3 };
	double B[] = { -1,-3,3 };
	double D[] = { 3,1,3 };
	double C[] = { 6,-1,3 };
	double G[] = { -2,6,3 };
	double E[] = { 7,4,3 };
	double F[] = { 4,8,3 };


	glBegin(GL_TRIANGLES);

	SearchNormal(K, A, B).Normal3d(-1);
	UpdatePoint(K).TexCoord2d();
	glVertex3dv(K);
	UpdatePoint(A).TexCoord2d();
	glVertex3dv(A);
	UpdatePoint(B).TexCoord2d();
	glVertex3dv(B);

	SearchNormal(K, B, D).Normal3d(-1);
	UpdatePoint(K).TexCoord2d();
	glVertex3dv(K);
	UpdatePoint(B).TexCoord2d();
	glVertex3dv(B);
	UpdatePoint(D).TexCoord2d();
	glVertex3dv(D);

	SearchNormal(B, C, D).Normal3d(-1);
	UpdatePoint(B).TexCoord2d();
	glVertex3dv(B);
	UpdatePoint(C).TexCoord2d();
	glVertex3dv(C);
	UpdatePoint(D).TexCoord2d();
	glVertex3dv(D);
	/*
	glVertex3dv(K);
	glVertex3dv(D);
	glVertex3dv(G);

	glVertex3dv(D);
	glVertex3dv(G);
	glVertex3dv(E);*/

	SearchNormal(G, F, E).Normal3d();
	UpdatePoint(G).TexCoord2d();
	glVertex3dv(G);
	UpdatePoint(F).TexCoord2d();
	glVertex3dv(F);
	UpdatePoint(E).TexCoord2d();
	glVertex3dv(E);

	SearchNormal(E, K, D).Normal3d(-1);
	UpdatePoint(E).TexCoord2d();
	glVertex3dv(E);
	UpdatePoint(K).TexCoord2d();
	glVertex3dv(K);
	UpdatePoint(D).TexCoord2d();
	glVertex3dv(D);

	double X = -4, Y = 2, R = 4.47;


	for (double t = 2.675; t <= 4.241; t += 0.001)
	{
		double O[] = { X - cos(t) * R , Y - sin(t) * R, 3 };
		double O1[] = { X - cos(t + 0.001) * R , Y - sin(t + 0.001) * R, 3 };
		SearchNormal(E, O, O1).Normal3d();
		UpdatePoint(E).TexCoord2d();
		glVertex3dv(E);
		UpdatePoint(O).TexCoord2d();
		glVertex3dv(O);
		UpdatePoint(O1).TexCoord2d();
		glVertex3dv(O1);
	}

	double O[] = { X - cos(4.242) * R , Y - sin(4.242) * R, 3 };
	SearchNormal(E, G, O).Normal3d(-1);
	UpdatePoint(E).TexCoord2d();
	glVertex3dv(E);
	UpdatePoint(G).TexCoord2d();
	glVertex3dv(G);
	UpdatePoint(O).TexCoord2d();
	glVertex3dv(O);



	Cic(3, -1);


	glEnd();
	//-------------------------------------------------
	//выключаем смешивание
	if (alpha != 0) {
		glDisable(GL_BLEND);
	}
	//-------------------------------------------------
}

//низ призмы
void Bottom() 
{
	double KK[] = { 0,0,1 };
	double AA[] = { -4,0,1 };
	double BB[] = { -1,-3,1 };
	double DD[] = { 3,1,1 };
	double CC[] = { 6,-1,1 };
	double GG[] = { -2,6,1 };
	double EE[] = { 7,4,1 };
	double FF[] = { 4,8,1 };

	glBegin(GL_TRIANGLES);


	SearchNormal(KK, AA, BB).Normal3d();
	UpdatePoint(KK).TexCoord2d();
	glVertex3dv(KK);
	UpdatePoint(AA).TexCoord2d();
	glVertex3dv(AA);
	UpdatePoint(BB).TexCoord2d();
	glVertex3dv(BB);

	SearchNormal(KK, BB, DD).Normal3d();
	UpdatePoint(KK).TexCoord2d();
	glVertex3dv(KK);
	UpdatePoint(BB).TexCoord2d();
	glVertex3dv(BB);
	UpdatePoint(DD).TexCoord2d();
	glVertex3dv(DD);

	SearchNormal(BB, CC, DD).Normal3d();
	UpdatePoint(BB).TexCoord2d();
	glVertex3dv(BB);
	UpdatePoint(CC).TexCoord2d();
	glVertex3dv(CC);
	UpdatePoint(DD).TexCoord2d();
	glVertex3dv(DD);
	/*
	glVertex3dv(KK);
	glVertex3dv(DD);
	glVertex3dv(GG);

	glVertex3dv(DD);
	glVertex3dv(GG);
	glVertex3dv(EE);*/

	SearchNormal(GG, FF, EE).Normal3d(-1);
	UpdatePoint(GG).TexCoord2d();
	glVertex3dv(GG);
	UpdatePoint(FF).TexCoord2d();
	glVertex3dv(FF);
	UpdatePoint(EE).TexCoord2d();
	glVertex3dv(EE);

	SearchNormal(EE, KK, DD).Normal3d();
	UpdatePoint(EE).TexCoord2d();
	glVertex3dv(EE);
	UpdatePoint(KK).TexCoord2d();
	glVertex3dv(KK);
	UpdatePoint(DD).TexCoord2d();
	glVertex3dv(DD);


	double X = -4, Y = 2, R = 4.47;

	for (double t = 2.675; t <= 4.241; t += 0.001)
	{
		double O[] = { X - cos(t) * R , Y - sin(t) * R, 1 };
		double O1[] = { X - cos(t + 0.001) * R , Y - sin(t + 0.001) * R, 1 };
		SearchNormal(EE, O, O1).Normal3d(-1);
		UpdatePoint(EE).TexCoord2d();
		glVertex3dv(EE);
		UpdatePoint(O).TexCoord2d();
		glVertex3dv(O);
		UpdatePoint(O1).TexCoord2d();
		glVertex3dv(O1);
	}

	double O[] = { X - cos(4.242) * R , Y - sin(4.242) * R, 1 };
	SearchNormal(EE, GG, O).Normal3d();
	UpdatePoint(EE).TexCoord2d();
	glVertex3dv(EE);
	UpdatePoint(GG).TexCoord2d();
	glVertex3dv(GG);
	UpdatePoint(O).TexCoord2d();
	glVertex3dv(O);


	Cic(1);

	glEnd();
}

void Side() //боковые поверхности призмы
{
	double K[] = { 0,0,3 };
	double A[] = { -4,0,3 };
	double B[] = { -1,-3,3 };
	double D[] = { 3,1,3 };
	double C[] = { 6,-1,3 };
	double G[] = { -2,6,3 };
	double E[] = { 7,4,3 };
	double F[] = { 4,8,3 };

	double KK[] = { 0,0,1 };
	double AA[] = { -4,0,1 };
	double BB[] = { -1,-3,1 };
	double DD[] = { 3,1,1 };
	double CC[] = { 6,-1,1 };
	double GG[] = { -2,6,1 };
	double EE[] = { 7,4,1 };
	double FF[] = { 4,8,1 };

	glBegin(GL_QUADS);
	vector<PointXY> texCoord = { {1, 1},  {0, 1}, {0, 0}, {1, 0}, {1, 1} };

	SearchNormal(A, B, BB).Normal3d();
	texCoord[0].TexCoord2d();
	glVertex3dv(A);
	texCoord[1].TexCoord2d();
	glVertex3dv(B);
	texCoord[2].TexCoord2d();
	glVertex3dv(BB);
	texCoord[3].TexCoord2d();
	glVertex3dv(AA);

	/*glColor3d(1, 0, 0);
	glVertex3dv(B);
	glVertex3dv(C);
	glVertex3dv(CC);
	glVertex3dv(BB);*/

	SearchNormal(C, D, DD).Normal3d();
	texCoord[0].TexCoord2d();
	glVertex3dv(C);
	texCoord[1].TexCoord2d();
	glVertex3dv(D);
	texCoord[2].TexCoord2d();
	glVertex3dv(DD);
	texCoord[3].TexCoord2d();
	glVertex3dv(CC);

	SearchNormal(D, E, EE).Normal3d();
	texCoord[0].TexCoord2d();
	glVertex3dv(D);
	texCoord[1].TexCoord2d();
	glVertex3dv(E);
	texCoord[2].TexCoord2d();
	glVertex3dv(EE);
	texCoord[3].TexCoord2d();
	glVertex3dv(DD);

	SearchNormal(E, F, FF).Normal3d();
	texCoord[0].TexCoord2d();
	glVertex3dv(E);
	texCoord[1].TexCoord2d();
	glVertex3dv(F);
	texCoord[2].TexCoord2d();
	glVertex3dv(FF);
	texCoord[3].TexCoord2d();
	glVertex3dv(EE);

	SearchNormal(F, G, GG).Normal3d();
	texCoord[0].TexCoord2d();
	glVertex3dv(F);
	texCoord[1].TexCoord2d();
	glVertex3dv(G);
	texCoord[2].TexCoord2d();
	glVertex3dv(GG);
	texCoord[3].TexCoord2d();
	glVertex3dv(FF);
	/*
	glVertex3dv(G);
	glVertex3dv(K);
	glVertex3dv(KK);
	glVertex3dv(GG);*/

	SearchNormal(K, A, AA).Normal3d();
	texCoord[0].TexCoord2d();
	glVertex3dv(K);
	texCoord[1].TexCoord2d();
	glVertex3dv(A);
	texCoord[2].TexCoord2d();
	glVertex3dv(AA);
	texCoord[3].TexCoord2d();
	glVertex3dv(KK);

	glEnd();
}


void Render(OpenGL *ogl)
{



	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);


	//альфаналожение
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//настройка материала
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//фоновая
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//дифузная
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//зеркальная
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//размер блика
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);
	//===================================
	//Прогать тут  
	

	//---------------------------------------
	//смена текстуры
	if (texChange) {
		glBindTexture(GL_TEXTURE_2D, texId[0]);
	}
	else {
		glBindTexture(GL_TEXTURE_2D, texId[1]);
	}
	//---------------------------------------

	glColor3d(0.8, 0.7, 0.4);
	
	Side();

	glBegin(GL_QUADS);
	//выпуклость
	Lateral();
	//впуклость
	Vpyklost();
	glEnd();

	//------------------------------------------------------------
	double K[] = { 0,0,3 };
	double A[] = { -4,0,3 };
	double B[] = { -1,-3,3 };
	double D[] = { 3,1,3 };
	double C[] = { 6,-1,3 };
	double G[] = { -2,6,3 };
	double E[] = { 7,4,3 };
	double F[] = { 4,8,3 };

	//-------------------------------------------------------------
	vector<Point> points = { A, B, C, D, E, F, G, K };
	//производим настройку под нашу фигуру
	UpdatePoint(F, true, points);
	//-------------------------------------------------------------


	glColor3d(1, 0, 0);
	Bottom();

	glColor3d(0.8, 0.7, 0.3);
	//прозрачную рисуем в самом конце
	Up();
	
	
	



   //Сообщение вверху экрана

	
	glMatrixMode(GL_PROJECTION);	//Делаем активной матрицу проекций. 
	                                //(всек матричные операции, будут ее видоизменять.)
	glPushMatrix();   //сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
	glLoadIdentity();	  //Загружаем единичную матрицу
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //врубаем режим ортогональной проекции

	glMatrixMode(GL_MODELVIEW);		//переключаемся на модел-вью матрицу
	glPushMatrix();			  //сохраняем текущую матрицу в стек (положение камеры, фактически)
	glLoadIdentity();		  //сбрасываем ее в дефолт

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //классик моего авторства для удобной работы с рендером текста.
	rec.setSize(300, 200);
	rec.setPosition(10, ogl->getHeight() - 200 - 10);


	std::stringstream ss;
	//--------------------------------------------------------
	ss << "Q - сменить текстуру" << std::endl;
	ss << "W - сменить вид прозрачности" << std::endl;
	//--------------------------------------------------------
	ss << "T - вкл/выкл текстур" << std::endl;
	ss << "L - вкл/выкл освещение" << std::endl;
	ss << "F - Свет из камеры" << std::endl;
	ss << "G - двигать свет по горизонтали" << std::endl;
	ss << "G+ЛКМ двигать свет по вертекали" << std::endl;
	ss << "Коорд. света: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "Коорд. камеры: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "Параметры камеры: R=" << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //восстанавливаем матрицы проекции и модел-вью обратьно из стека.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}