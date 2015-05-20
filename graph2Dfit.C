using namespace std;
#include <vector>
#include <iostream>
#include <fstream> 
#include <string>

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

void graph2dfit(string name, double percentage = 100, double xy_scale_nm = 1)
{
	bool write_matrix = 0;
	bool write_xyz = 0;
	bool draw_surface = true;
	bool calculate = 0;

	


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


	///////////////////
	//fitting
	TF2  *f2 = new TF2("f2", " [0]*x*x + [1]*y*y + 2*[2]*x*y + 2*[3]*x + 2*[4]*y + [5] ", 0, x_dimension*xy_scale_nm, 0, y_dimension*xy_scale_nm);
	//f2->SetParameters(1, 1);
	/////////////////////


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



	for (int y = 0; y < (int)(y_dimension * percentage / 100.0); y++)
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
		
		//dt->Draw("surf1");

		dt->Fit(f2);
		f2->Draw("surf1");
		dt->Draw("same p");
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