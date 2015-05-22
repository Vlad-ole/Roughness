using namespace std;
#include <vector>
#include <iostream>
#include <fstream> 
#include <string>

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

struct Point
{
	double x = 0;
	double y = 0;
	double z = 0;
};

double get_alpha(/*struct*/ Point p0, Point p1, Point p2)
{
	Point coeff;

	double x0 = p0.x;
	double y0 = p0.y;
	double z0 = p0.z;

	double x1 = p1.x;
	double y1 = p1.y;
	double z1 = p1.z;

	double x2 = p2.x;
	double y2 = p2.y;
	double z2 = p2.z;

	double A, B, C;

	A = (y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0);
	B = (x1 - x0) * (z2 - z0) - (z1 - z0) * (x2 - x0);
	C = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);

	double alpha = acos( C / sqrt( A*A + B*B + C*C ) );

	
	//cout << "alpha_res = " << alpha << endl;
	return alpha;
}

void test_get_alpha()
{
	Point p0;
	Point p1;
	Point p2;

	double p0.x = 1;
	double p0.y = 1;
	double p0.z = 3;

	double p1.x = 1;
	double p1.y = 0;
	double p1.z = 1;

	double p2.x = 0;
	double p2.y = 1;
	double p2.z = 2;

	cout << get_alpha(p0, p1, p2) << " (radian)" << endl;
	cout << get_alpha(p0, p1, p2) * 180/3.1416 << " (degree)" << endl;
}

void graph2d(string name, double percentage = 100, double xy_scale_nm = 1)
{
	bool write_matrix = true;
	bool write_xyz = false;
	bool draw_surface = true;
	bool calculate = true;
	
	string name_out = name + "__result";

	FILE *f = fopen(name.c_str(), "rb");
	if (f == NULL)
	{
		cout << "Cann't open f file" << endl;
		exit(1);
	}

	
	 
	vector<double> xv;
	vector<double> yv;
	vector<double> zv;

	int x_dimension;
	int y_dimension;
	float z_coord;
	
	vector< vector<double> > vec;

	fread(&x_dimension, 4, 1, f);
	fread(&y_dimension, 4, 1, f);

	cout << "x_dimension = \t" << x_dimension << endl;
	cout << "y_dimension = \t" << y_dimension << endl;

	for (int i = 0; i < x_dimension * y_dimension; i++)
	{
		fread(&z_coord, 4, 1, f);
		 
		zv.push_back(z_coord);
	}

	
	cout << "(int)(y_dimension * percentage/100.0) = \t" << (int)(y_dimension * percentage / 100.0) << endl;

	cout << "There are " << x_dimension * (int)(y_dimension * percentage / 100.0) << " points" << endl;

	char ch;

	//if ((int)(y_dimension * percentage / 100.0) * x_dimension > 1000)
	//{
	//	cout << "WARNING! There are more than 1000 points for your graph. Construction of the graph will take a lot time." << endl;
	//	cout << "Are you sure that you want to continue? (Type 'y' or 'n')" << endl;
	//	
	//	do
	//	{
	//		ch = getchar();

	//	} while ( !((ch == 'y') || (ch == 'n')) );

	//	if (ch == 'n')
	//		exit(1);

	//}


	if (draw_surface)
	{
		TCanvas *c = new TCanvas("c", "Graph2D example", 0, 0, 700, 600);
		TGraph2D *dt = new TGraph2D();
	}
	
	int N = 0;
	//double z_mean = 0;

	

	for (int y = 0; y < (int)(y_dimension * percentage/100.0); y++)
	{
		vector<double> row; // Create an empty row
		for (int x = 0; x < x_dimension; x++)
		{
			if (write_matrix)
				cout << zv[x + y * x_dimension] << "\t";

			if (write_xyz)
				cout << x * xy_scale_nm << "\t" << y * xy_scale_nm << "\t" << zv[x + y * x_dimension] << endl;
			
			row.push_back(zv[x + y * x_dimension]); // Add an element (column) to the row
			
			if (draw_surface)
				dt->SetPoint(N, x * xy_scale_nm, y * xy_scale_nm, zv[x + y * x_dimension]);
			
			//z_mean += zv[x + y * x_dimension];
			
			N++;
		}

		vec.push_back(row); // Add the row to the main vector

		if (write_matrix)
			cout << endl;
	}


	fclose(f);

	if (draw_surface)
	{
		gStyle->SetPalette(1);
		//dt->Draw("P");
		dt->Draw("surf1");
	}

	//dt->Draw("surf1");
	//dt->Draw("COLz");
	//dt->Draw("tri1 p0");


	//for (int j = 0; j < (int)(y_dimension * percentage / 100.0); j++)
	//{
	//	for (int i = 0; i < x_dimension; i++)
	//	{
	//		cout << vec[j][i] << "\t";
	//	}

	//	cout << endl;
	//}

	Point p0;
	Point p1;
	Point p2;

	double alpha;

	ofstream myStream(name_out.c_str());
	
	if (calculate)
	{
		for (int j = 0; j < ((int)(y_dimension * percentage / 100.0) - 1); j++)
		{
			for (int i = 0; i < (x_dimension - 1); i++)
			{
				p0.x = i * xy_scale_nm;
				p0.y = j * xy_scale_nm;
				p0.z = vec[j][i];

				p1.x = (i + 1) * xy_scale_nm;
				p1.y = j * xy_scale_nm;
				p1.z = vec[j][i + 1];

				p2.x = i * xy_scale_nm;
				p2.y = (j + 1) * xy_scale_nm;
				p2.z = vec[j + 1][i];

				//cout << "alpha = \t" << get_alpha(p0, p1, p2) * 180/3.1416 << endl;

				alpha = get_alpha(p0, p1, p2);
				myStream << alpha << endl;

				//--------

				p0.x = (i + 1) * xy_scale_nm;
				p0.y = (j + 1) * xy_scale_nm;
				p0.z = vec[j + 1][i + 1];

				alpha = get_alpha(p0, p1, p2);
				myStream << alpha << endl;


			}

			if (j % 10)
			{
				cout.precision(4);
				cout << (j / ((y_dimension * percentage / 100.0) - 1)) * 100 << " %" << endl;
			}

		}

	}
	
	myStream.close();

	
		
}




void create_test_binary(char name[], int x_dim = 10, int y_dim = 10)
{
	FILE *f = fopen(name, "wb");
	if (f == NULL)
	{
		cout << "Cann't open file" << endl;
		exit(1);
	}


	float z_coord;

	fwrite(&x_dim, 4, 1, f);
	fwrite(&y_dim, 4, 1, f);

	for (int i = 0; i < x_dim * y_dim; i++)
	{
		//z_coord = pow(-1, i);
		z_coord = i*i + i;
	
		fwrite(&z_coord, 4, 1, f);
	}

	fclose(f);
	
}

void hist_sin(char name[], bool IsDivideSin = false)
{	
	TCanvas *c1 = new TCanvas("c1", "A Simple Graph Example", 200, 10, 700, 500);
	c1->SetGrid();


	TH1F *h1 = new TH1F("h1f", "hist", 1000, -5, 5);

	FILE *f = fopen(name, "r");
	if (f == NULL)
	{
		cout << "Cann't open f file" << endl;
		exit(1);
	}


	double x;
	while (!feof(f))
	{
		fscanf(f, "%lf\n", &x);

		if (x > pi/2.0)
			x = x - pi;

		h1->Fill(x * 180 / pi);
	}

	fclose(f);
	// angle has degree dimention now!

	if (IsDivideSin)
	{
		double x_temp;
		double y_temp;
		for (int i = 0; i < h1->GetSize(); i++)
		{
			x_temp = h1->GetBinCenter(i);
			y_temp = h1->GetBinContent(i);

			if (x_temp > 0)
				h1->SetBinContent(i, y_temp / ( sin(x_temp * pi / 180) ) );
			else 
				h1->SetBinContent(i, - y_temp / (sin(x_temp * pi / 180)));

			//cout << x_temp << "\t" << y_temp << endl;
		}
	}

	
	
	h1->SetFillColor(kRed);
	h1->SetFillStyle(3002);
	h1->Draw();	

	c1->Update();
}

//Fitting a TGraph2D
//Author: Olivier Couet

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>

TCanvas* graph2dfit_test()
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();

	TCanvas *c = new TCanvas("c", "Graph2D example", 0, 0, 600, 800);
	
	Double_t rnd, x, y, z;
	Double_t e = 0.3;
	Int_t nd = 400;
	Int_t np = 10000;

	TRandom r;
	Double_t fl = 6;
	TF2  *f2 = new TF2("f2", "1000*(([0]*sin(x)/x)*([1]*sin(y)/y))+200", -fl, fl, -fl, fl);
	f2->SetParameters(1, 1);
	TGraph2D *dt = new TGraph2D();

	// Fill the 2D graph
	Double_t zmax = 0;
	for (Int_t N = 0; N<nd; N++) {
		f2->GetRandom2(x, y);
		// Generate a random number in [-e,e]
		rnd = 2 * r.Rndm()*e - e;
		z = f2->Eval(x, y)*(1 + rnd);
		if (z>zmax) zmax = z;
		dt->SetPoint(N, x, y, z);
	}


	f2->SetParameters(0.5, 1.5);
	dt->Fit(f2);
	TF2 *fit2 = (TF2*)dt->FindObject("f2");

	f2->SetParameters(1, 1);

	for (Int_t N = 0; N<np; N++) {
		f2->GetRandom2(x, y);
		// Generate a random number in [-e,e]
		rnd = 2 * r.Rndm()*e - e;
		z = f2->Eval(x, y)*(1 + rnd);
		h1->Fill(f2->Eval(x, y) - z);
		z = dt->Interpolate(x, y);
		h2->Fill(f2->Eval(x, y) - z);
		z = fit2->Eval(x, y);
		h3->Fill(f2->Eval(x, y) - z);
	}

	gStyle->SetPalette(1);
	
	f2->SetTitle("Original function with Graph2D points on top");
	f2->SetMaximum(zmax);
	gStyle->SetHistTopMargin(0);
	f2->Draw("surf1");
	dt->Draw("same p0");



	return c;
}