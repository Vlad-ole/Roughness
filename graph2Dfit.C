using namespace std;
#include <vector>
#include <iostream>
#include <fstream> 
#include <string>

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;


double surf_func(double p0, double p1, double p2, double p3, double p4, double p5, double x, double y)
{
	return p0*x*x + p1*y*y + 2*p2*x*y + 2*p3*x + 2*p4*y + p5;
}

void graph2dfit(string name, double percentage = 100, double xy_scale_nm = 1)
{
	bool write_matrix = 0;
	bool write_xyz = 0;
	bool draw_surface = true;
	bool calculate = 1;
	bool fit_subtraction = 1;


	//read file
	/////////////////////////////////////////
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


	fread(&x_dimension, 4, 1, f);
	fread(&y_dimension, 4, 1, f);

	cout << "x_dimension = \t" << x_dimension << endl;
	cout << "y_dimension = \t" << y_dimension << endl;

	for (int i = 0; i < x_dimension * y_dimension; i++)
	{
		fread(&z_coord, 4, 1, f);

		zv.push_back(z_coord);
	}
	fclose(f);

	cout << "(int)(y_dimension * percentage/100.0) = \t" << (int)(y_dimension * percentage / 100.0) << endl;
	cout << "There are " << x_dimension * (int)(y_dimension * percentage / 100.0) << " points" << endl;
	/////////////////////////////////////////


	TCanvas *c = new TCanvas("c", "Graph2D example", 0, 0, 1500, 600);
	c->Divide(3, 1);
	TGraph2D *dt = new TGraph2D();
	TGraph2D *dt_supp = new TGraph2D();


	int N = 0;
	for (int y = 0; y < (int)(y_dimension * percentage / 100.0); y++)
	{
		for (int x = 0; x < x_dimension; x++)
		{
			dt->SetPoint(N, x * xy_scale_nm, y * xy_scale_nm, zv[x + y * x_dimension]);
			N++;
		}
	}

	TF2  *f2 = new TF2("f2", " [0]*x*x + [1]*y*y + 2*[2]*x*y + 2*[3]*x + 2*[4]*y + [5] ", 0, x_dimension*xy_scale_nm, 0, y_dimension*xy_scale_nm);
	dt->Fit(f2);
	
	
	c->cd(1);
	gStyle->SetPalette(1);
	dt->Draw("p");

	c->cd(2);
	gStyle->SetPalette(1);
	if (fit_subtraction)
	{
		f2->Draw("surf1");
		dt->Draw("same p");
	}
	else
	{
		dt->Draw("p");
	}
	
		const double par0 = f2->GetParameter(0);
		const double par1 = f2->GetParameter(1);
		const double par2 = f2->GetParameter(2);
		const double par3 = f2->GetParameter(3);
		const double par4 = f2->GetParameter(4);
		const double par5 = f2->GetParameter(5);

	
	double z;
	vector< vector<double> > vec;
	N = 0;
	for (int y = 0; y < (int)(y_dimension * percentage / 100.0); y++)
	{
		vector<double> row; // Create an empty row
		
		for (int x = 0; x < x_dimension; x++)
		{
			if (fit_subtraction)
				z = zv[x + y * x_dimension] - surf_func(par0, par1, par2, par3, par4, par5, x*xy_scale_nm, y*xy_scale_nm);
			else
				z = zv[x + y * x_dimension];
			
			dt_supp->SetPoint(N, x * xy_scale_nm, y * xy_scale_nm, z);
			row.push_back(z); // Add an element (column) to the row
			N++;
		}

		vec.push_back(row); // Add the row to the main vector
	}

	c->cd(3);
	dt_supp->Draw("p");


	if (calculate)
	{

		Point p0;
		Point p1;
		Point p2;

		double alpha;

		ofstream myStream(name_out.c_str());


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

		myStream.close();

	}


}