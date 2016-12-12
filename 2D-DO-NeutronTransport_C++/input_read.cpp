#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//*******************************************************************************************************************************

void input_read(int &N, int &Nx, double &xL, double &xR, int &Ny, double &yB, double &yT,
				int &bc, std::string &sigtfid, std::string &sigsfid, std::string &nusigffid, double &tol, 
				std::string &srcfid, int &Egrp, std::string &calc_mode, std::string &chifid )

//*******************************************************************************************************************************

{

	string line;
	int idx = 0;

	// Ask for input file
	string fid;
	std::cout << "Enter input file with extension: ";
	cin  >> fid;

	cout << "Echoing input file: " << '\n';
	cout << "*******************************************************************************************************************************" << '\n';
	cout << '\n';

	ifstream inpfile;
	inpfile.open(fid);

	if ( !inpfile )
	{
		cout << "File does not exist!!!" << '\n';
		exit(EXIT_FAILURE);
	}

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
			sigtfid = line;
			cout << sigtfid << '\n';
			idx++;
		}
		else if ( idx == 20 )
		{
			getline(inpfile,line,'\n');
			sigsfid = line;
			cout << sigsfid << '\n';
			idx++;
		}
		else if ( idx == 22 )
		{
			getline(inpfile,line,'\n');
			tol = stod(line);
			cout << tol << '\n';
			idx++;
		}
		else if ( idx == 24 )
		{
			getline(inpfile,line,'\n');
			srcfid = line;
			cout << srcfid << '\n';
			idx++;
		}
		else if ( idx == 26 )
		{
			getline(inpfile,line,'\n');
			Egrp = stod(line);
			cout << Egrp << '\n';
			idx++;
		}
		else if ( idx == 28 )
		{
			getline(inpfile,line,'\n');
			nusigffid = line;
			cout << nusigffid << '\n';
			idx++;
		}
		else if ( idx == 30 )
		{
			getline(inpfile,line,'\n');
			calc_mode = line;
			cout << calc_mode << '\n';
			idx++;
		}
		else if (idx == 32 )
		{
			getline(inpfile,line,'\n');
			chifid = line;
			cout << chifid << '\n';
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

	//cout << '\n';

}
