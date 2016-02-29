#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream> //std
#include <sstream>
#include <fstream> //���� �����
#include <locale.h> //���� � ��� ������, �� ������ ��� ������� ������
#include <time.h>
#include <omp.h>

using namespace std;

static int n;

double generate_double(int a, int b) //�������� ��������� 
{
	double r2 = 0; srand(time(NULL));
	int ifi = (rand() % 4);
	if (ifi == 1)
	{
		r2 = (rand() % (b - a + 1) * a) + ((rand() % (b - a + 1) * a) * 0.1);
	}
	else if (ifi == 2)
	{
		r2 = (rand() % (b - a + 1) * a) + ((rand() % (b - a + 1) * a) * 0.01);
	}
	else if (ifi == 3)
	{
		r2 = (rand() % (b - a + 1) * a) + ((rand() % (b - a + 1) * a) * 0.001);
	}
	else
	{
		r2 = (rand() % (b - a + 1) * a) + ((rand() % (b - a + 1) * a) * 0.0001);
	}
	return r2;
}

double *sort_mass(double *a = new double[n])
{
#pragma omp parallel for
	for (int i = 0; i < n - 1; i++)
	{
		int min_index = i;
		for (int j = i + 1; j < n; j++)
		{
			if (a[j] > a[min_index]) min_index = j;
		}
		if (min_index != i) swap(a[i], a[min_index]);
	}
	return a;
}

double formula(bool tr, double a, double b, double c, double d)
{
	if (tr == true)
	{
		return double((a*(b / c - a)) / (2 * d) + (c / d*a / 2) / (a*b*c*d*0.3));
	}
	else
	{
		return double((a / (b / c - a)) / (2 * d) * 10 + 20 * (c / d*a / 2) / (a*b*c*d*0.3));

	}
}

int main()
{
	srand(time(NULL));
	setlocale(LC_ALL, "Russian"); //���� � ��� ������, �� ������ ��� ������� ������
	cout << "������� ����������� �������� ";
	cin >> n; //����� ����������� ��� ���� ������
	double *a = new double[n], *b = new double[n], *c = new double[n], *d = new double[n];

	double *res1 = new double[n], *res2 = new double[n];
	int proc = 0;
	int cn = omp_get_max_threads();
	cout << "���������� ����� �������� - " << cn << endl;
	cout << "���������� ��������" << endl;
#pragma omp parallel
	{
#pragma omp sections
	{
#pragma omp section
		a[0] = generate_double(1, 1000000);
#pragma omp section
		b[0] = generate_double(1, 3000000);
#pragma omp section
		c[0] = generate_double(1, 2000000);
#pragma omp section
		d[0] = generate_double(1, 4000000);
	}
#pragma omp for schedule(dynamic, 15)
	for (int i = 0;i<n;i++)
	{
		if ((((i * 100) / n) % 10 == 0) && ((i * 100 / n) > proc))
		{
			proc = i * 100 / n;
			cout << proc;
			cout << "% ... ";
		}
		a[i] = generate_double(1, 1000000);
		b[i] = generate_double(1, 3000000);
		c[i] = generate_double(1, 2000000);
		d[i] = generate_double(1, 4000000);
	}
	//}

#pragma omp single nowait
	{
		cout << endl << "������� ���������" << endl;
		/*����� � ���� �������� �������� {
			ofstream out;
			a1 = a; b1 = b; c1 = c; d1 = d;
			out.open("output1_start.txt");
			out << "�������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << a[i]; out << +"   "; }
			out << endl << endl << "�������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << b[i]; out << +"    "; }
			out << endl << endl << "�������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << c[i]; out << +"    "; }
			out << endl << endl << "�������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << d[i]; out << +"    "; }
			out.close(); a1 = NULL; b1 = NULL; c1 = NULL; d1 = NULL;
		}*/
		cout << "���������� ����������" << endl;
	}

#pragma omp sections //����������
	{
#pragma omp section
		a = sort_mass(a);
#pragma omp section 
		b = sort_mass(b);
#pragma omp section 
		c = sort_mass(c);
#pragma omp section 
		d = sort_mass(d);
	}

#pragma omp single nowait
	{
		cout << "���������� ���������" << endl;
		/*����� � ���� ������������� �������� {
			ofstream out;
			out.open("output2_sort.txt"); a2 = a; b2 = b; c2 = c; d2 = d;
			out << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << a[i]; out << +"   "; }
			out << endl << endl << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << b[i]; out << +"    "; }
			out << endl << endl << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << c[i]; out << +"    "; }
			out << endl << endl << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << d[i]; out << +"    "; }
			out.close(); a2 = NULL; b2 = NULL; c2 = NULL; d2 = NULL;
		}*/
		cout << "��������� ��������� �������� ��� ���������� ��������" << endl;
	}

#pragma omp for schedule(dynamic, 1) //���������� ������� �� ���������
	for (int j = 0; j < n; j++)
	{
		if (((int)a[j] % 3) == 0)
		{
			a[j] += a[j] * 0.1;
		}
		if (((int)b[j] % 3) == 0)
		{
			b[j] += b[j] * 0.1;
		}
		if ((int)a[j] % 5 == 0 || (int)a[j] % 5 == 0)
		{
			a[j] += (a[j] * 0.6)*0.31;
		}
		if ((int)b[j] % 5 == 0 || (int)b[j] % 5 == 0)
		{
			b[j] += (b[j] * 0.6)*0.31;
		}
		if ((int)c[j] % 5 == 0 || (int)c[j] % 5 == 0)
		{
			c[j] += (c[j] * 0.6)*0.31;
		}
		if ((int)d[j] % 5 == 0 || (int)d[j] % 5 == 0)
		{
			d[j] += (d[j] * 0.6)*0.31;
		}
	}

#pragma omp single nowait
	{
		cout << "�������� ���������" << endl;
		/*����� � ���� ���������� ������� {
			ofstream out;
			out.open("output3_edit_mass.txt"); a3 = a; b3 = b; c3 = c; d3 = d;
			out << "���������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << a[i]; out << +"   "; }
			out << endl << endl << "���������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << b[i]; out << +"    "; }
			out << endl << endl << "���������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << c[i]; out << +"    "; }
			out << endl << endl << "���������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << d[i]; out << +"    "; }
			out.close(); a3 = NULL; b3 = NULL; c3 = NULL; d3 = NULL;
		}*/
		cout << "��������� ����������" << endl;
	}

#pragma omp sections //����������
	{
#pragma omp section
		a = sort_mass(a);
#pragma omp section 
		b = sort_mass(b);
#pragma omp section 
		c = sort_mass(c);
#pragma omp section 
		d = sort_mass(d);
	}

#pragma omp single nowait
	{
		cout << "��������� ���������� ���������" << endl;
		/*��������� ���������� � ���� {
			ofstream out;
			out.open("output3_sort.txt"); a4 = a; b4 = b; c4 = c; d4 = d;
			out << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << a[i]; out << +"   "; }
			out << endl << endl << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << b[i]; out << +"    "; }
			out << endl << endl << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << c[i]; out << +"    "; }
			out << endl << endl << "������������� ������ '�':" << endl;
			for (int i = 0; i < n; i++) { out << d[i]; out << +"    "; }
			out.close(); a4 = NULL; b4 = NULL; c4 = NULL; d4 = NULL;
		}*/
		cout << "���������� ����� ��������" << endl;
	}

#pragma omp for schedule(dynamic, 5) nowait
	for (int i = 0;i < n;i++)
	{
		res1[i] = formula(true, a[i], b[i], c[i], d[i]);
		while (res1[i] > 499) res1[i] = res1[i] / 2;
	}
#pragma omp for schedule(dynamic, 5)
	for (int i = 0;i < n;i++)
	{
		res2[i] = formula(false, a[i], b[i], c[i], d[i]);
		while (res2[i] > 199) res2[i] = res2[i] / 3;
	}


#pragma omp sections
	{
#pragma omp section
	{
		ofstream out;
		out.open("res1.txt");
		for (int i = 0; i < n; i++)
		{
			out << res1[i]; out << "   ";
		}
		out.close();
	}
#pragma omp section
	{
		ofstream out;
		out.open("res2.txt");
		for (int i = 0; i < n; i++)
		{
			out << res2[i]; out << "    ";
		}
		out.close();
	}
	}
	}
	system("pause");
	return 0;

}

/*
1.	������ �������;
2.	������������� ������� �� �����������;
3.	��������� �������� ������� 3 - � � 1 � 2 �������� �� 10 % ;
4.	�������� ������� 5 ��� 6 ��������� �� 31 % �� 60 % ����� �����.
5.	������ ������������� ������� �� �����������.
6.	�������� 2 �������������� �������;
a.������ ������ res1 = (a*(b / c - a)) / (2 * d) + (c / d*a / 2) / (a*b*c*d*0.3)    if res1[i] > 500  ������ �� 2 ���� �� ����� ������ 500;
b.������ ������ res2 = (a / (b / c - a)) / (2 * d) * 10 + 20 * (c / d*a / 2) / (a*b*c*d*0.3) if res2[i] > 200 ������ �� 3 ���� �� ����� ������ 200;
*/