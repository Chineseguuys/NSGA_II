#pragma once

#ifndef _TESTFUNCTION_
#define _TESTFUNCTION_

#ifndef _IOSTREAM_
#include<iostream>
#endif // !_IOSTREAM_

#ifndef _INC_MATH_
#include<math.h>
#endif // !_INC_MATH_

#ifndef _PI_
#define PI  3.141592653589793238462643383279
#endif // !_PI_

#ifdef ZDT1
void test_problem(double* x_real,double* obj)
{
	double f1, f2, g, h;
	int i;
	f1 = x_real[0];
	g = 0.0;
	for (int i = 1; i < 30; i++)
	{
		g += 9 * x_real[i];
	}
	g = g / 29.0;
	g += 1.0;
	h = (1.0 - sqrt(x_real[0] / g));
	f2 = g*h;
	obj[0] = f1;
	obj[1] = f2;
	return;
}

double limits_array[30][2] = { {0.0,1.0},
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
							   { 0.0,1.0 },
	};
#endif // ZDT1

#ifdef ZDT2
void test_problem(double* x_real, double* obj)
{
	double f1, f2, g, h;
	int i;
	f1 = x_real[0];
	g = 0.0;
	for(i=0;i<30;i++)
	{
		g += x_real[i];
	}
	g = 9.0*g / 29.0;
	g += 1.0;
	h = 1.0 - pow((f1 / g), 2);
	f2 = h*g;
	obj[0] = f1;
	obj[1] = f2;
	return;
}

double limits_array[30][2] = { { 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
	};

#endif // ZDT2

#ifdef ZDT3
void test_problem(double* x_real, double* obj)
{
	double f1, f2, g, h;
	int i;
	f1 = x_real[0];
	g = 0.0;
	for (int i = 0; i < 30; i++)
	{
		g += x_real[i];
	}
	g = 9.0*g / 29.0;
	g += 1.0;
	h = 1.0 - sqrt(f1 / g) - (f1 / g)*sin(10.0*PI*f1);
	f2 = h*g;
	obj[0] = f1;
	obj[1] = f2;
	return;
}

double limits_array[30][2] = { { 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
};

#endif // !ZDT3

#ifdef ZDT4
void test_problem(double* x_real, double* obj)
{
	double f1, f2, g, h;
	int i;
	f1 = x_real[0];
	g = 0.0;
	for (i = 1; i < 10; i++)
	{
		g += x_real[i] * x_real[i] - 10.0*cos(4.0*PI*x_real[i]);
	}
	g += 91.0;
	h = 1.0 - sqrt(f1 / g);
	f2 = g*h;
	obj[0] = f1;
	obj[1] = f2;
	return;
}

double limits_array[10][2] = { { 0.0,1.0 },
{ -5.0,5.0 },
{ -5.0,5.0 },
{ -5.0,5.0 },
{ -5.0,5.0 },
{ -5.0,5.0 },
{ -5.0,5.0 },
{ -5.0,5.0 },
{ -5.0,5.0 },
{ -5.0,5.0 }
};

#endif // ZDT4



#ifdef ZDT6
void test_problem(double* x_real, double* obj)
{
	double f1, f2, g, h;
	int i;
	f1 = 1.0 - (exp(-4.0*x_real[0]))*pow((sin(6.0*PI*x_real[0])), 6.0);
	g = 0.0;
	for(i =1;i<10;i++ )
	{
		g += x_real[i];
	}
	g = g / 9.0;
	g = pow(g, 0.25);
	g = 1.0 + 9.0*g;
	h = 1.0 - pow((f1 / g), 2.0);
	f2 = g*h;
	obj[0] = f1;
	obj[1] = f2;
	return;
}
double limits_array[10][2] = { { 0.0,1.0 },
{0.0,1.0},
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 }
};

#endif // ZDT6

#ifdef DTLZ1
void test_problem(double* x_real, double* obj)
{
	double f1, f2, f3, g;
	int i;
	g = 0.0;
	for(i=0;i<10;i++)
	{
		g += (pow((x_real[i] - 0.5), 2) - cos(20.0*PI*(x_real[i] - 0.5)));
	}
	g += 10.0;
	g = g*100.0;
	f1 = 0.5*x_real[0] * x_real[1] * (1 + g);
	f2 = 0.5*x_real[0] * (1 - x_real[1])*(1 + g);
	f3 = 0.5*(1 - x_real[0])*(1 + g);
	obj[0] = f1;
	obj[1] = f2;
	obj[2] = f3;
	return;
}

double limits_array[10][2] = { {0.0,1.0},
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 }
 };
#endif // DTLZ1

#ifdef DTLZ2
void test_problem(double* x_real, double* obj)
{
	double f1, f2, f3, g;
	int i;
	g = 0.0;
	for (i = 0; i < 10; i++)
	{
		g += pow((x_real[i] - 0.5), 2);
	}
	f1 = cos(0.5*PI*x_real[0])*cos(0.5*PI*x_real[1])*(1 + g);
	f2 = cos(0.5*PI*x_real[0])*sin(0.5*PI*x_real[1])*(1 + g);
	f3 = sin(0.5*PI*x_real[0])*(1 + g);
	obj[0] = f1;
	obj[1] = f2;
	obj[2] = f3;
	return;
}

double limits_array[10][2] = { { 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 }
 };


#endif // DTLZ2

#ifdef DTLZ3
void test_problem(double* x_real, double* obj)
{
	double f1, f2, f3, g;
	int i;
	g = 0.0;
	for (i = 0; i < 10; i++)
	{
		g += pow((x_real[i] - 0.5), 2) - cos(20.0*PI*(x_real[i] - 0.5));
	}
	g += 10.0;
	g = g*100;
	f1 = cos(0.5*PI*x_real[0])*cos(0.5*PI*x_real[1])*(1 + g);
	f2 = cos(0.5*PI*x_real[0])*sin(0.5*PI*x_real[1])*(1 + g);
	f3 = sin(0.5*PI*x_real[0])*(1 + g);
	obj[0] = f1;
	obj[1] = f2;
	obj[2] = f3;
	return;
}

double limits_array[10][2] = { { 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 },
{ 0.0,1.0 }
};
#endif // DTLZ3


#endif // !_TESTFUNCTION_
