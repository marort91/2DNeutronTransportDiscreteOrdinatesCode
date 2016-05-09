#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//*******************************************************************************************************************************

void input_read(int &N, int &Nx, double &xL, double &xR, int &Ny, double &yB, double &yT,
				int &bc, double &sigt, double &sigs0, double &tol)

//*******************************************************************************************************************************

{

	string line;
	int idx = 0;

	// Ask for input file
	string fid;
	std::cout << "Enter input file with extension: ";
	cin  >> fid;

	ifstream inpfile;
	inpfile.open(fid);

	while ( !inpfile.eof() )
	{
		if ( idx == 0 )
		{
			getline(inpfile, line );
			cout << line << '\n';
			idx++;
		}
		else if ( idx == 2 )
		{
			getline(inpfile,line,'\n');
			N = stoi(line);
			cout << N << '\n';
			idx++;
		}
		else if ( idx == 4 )
		{
			getline(inpfile,line,'\n');
			Nx = stoi(line);
			cout << Nx << '\n';
			idx++;
		}
		else if ( idx == 6 )
		{
			getline(inpfile,line,'\n');
			xL = stod(line);
			cout << xL << '\n';
			idx++;
		}
		else if ( idx == 8 )
		{
			getline(inpfile,line,'\n');
			xR = stod(line);
			cout << xR << '\n';
			idx++;
		}
		else if ( idx == 10 )
		{
			getline(inpfile,line,'\n');
			Ny = stoi(line);
			cout << Ny << '\n';
			idx++;
		}
		else if ( idx == 12 )
		{
			getline(inpfile,line,'\n');
			yB = stod(line);
			cout << yB << '\n';
			idx++;
		}
		else if ( idx == 14 )
		{
			getline(inpfile,line,'\n');
			yT = stod(line);
			cout << yT << '\n';
			idx++;
		}
		else if ( idx == 16 )
		{
			getline(inpfile,line,'\n');
			bc = stoi(line);
			cout << bc << '\n';
			idx++;
		}
		else if ( idx == 18 )
		{
			getline(inpfile,line,'\n');
			sigt = stod(line);
			cout << sigt << '\n';
			idx++;
		}
		else if ( idx == 20 )
		{
			getline(inpfile,line,'\n');
			sigs0 = stod(line);
			cout << sigs0 << '\n';
			idx++;
		}
		else if ( idx == 22 )
		{
			getline(inpfile,line,'\n');
			tol = stod(line);
			cout << tol << '\n';
			idx++;
		}
		else
		{
			getline(inpfile, line);
			cout << line << '\n';
			idx++;
		}
		
	}
	inpfile.close();

}