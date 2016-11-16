#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

//*******************************************************************************************************************************

void xs_file_read( std::vector<std::vector<std::vector<double> > > &sigt, std::string sigtfid, 
	               std::vector<std::vector<std::vector<double> > > &sigs, std::string sigsfid, 
	               std::vector<std::vector<std::vector<double> > > &nusigf, std::string nusigffid, int Nx, int Ny )
{

	int Egrp;

	std::string sigtline;
	std::ifstream sigtsrcfile;

	Egrp = 0;

	sigtsrcfile.open(sigtfid);

	if ( sigtsrcfile.is_open() )
	{
		getline(sigtsrcfile,sigtline);

		for ( int x = 0; getline(sigtsrcfile,sigtline); x++ )
		{
			if ( sigtline[0] == 'E')
			{
				Egrp = Egrp + 1;
				getline(sigtsrcfile,sigtline);
				x = 0;
			}

			std::istringstream iss(sigtline);

			for ( int y = 0; y < Ny; y++)
			{
				std::string sub;
				iss >> sub;
				sigt[x][y][Egrp] = ::atof(sub.c_str());
			}
		}
		sigtsrcfile.close();
	}
		
	std::string nusigfline;
	std::ifstream nusigfsrcfile;

	Egrp = 0;

	sigssrcfile.open(sigsfid);

	if ( sigssrcfile.is_open() )
	{
		getline(sigssrcfile,sigsline);

		for ( int x = 0; getline(sigssrcfile,sigsline); x++ )
		{
			if ( sigtline[0] == 'E')
			{
				Egrp = Egrp + 1;
				getline(sigssrcfile,sigsline);
				x = 0;
			}

			std::istringstream iss(sigsline);

			for ( int y = 0; y < Ny; y++)
			{
				std::string sub;
				iss >> sub;
				sigs[x][y][Egrp] = ::atof(sub.c_str());
			}
		}
		sigssrcfile.close();
	}

	std::string nusigfline;
	std::ifstream nusigfsrcfile;

	Egrp = 0;

	nusigfsrcfile.open(sigsfid);

	if ( nusigfsrcfile.is_open() )
	{
		getline(nusigfsrcfile,nusigfline);

		for ( int x = 0; getline(nusigfsrcfile,nusigfline); x++ )
		{
			if ( sigtline[0] == 'E')
			{
				Egrp = Egrp + 1;
				getline(nusigfsrcfile,nusigfline);
				x = 0;
			}

			std::istringstream iss(nusigfline);

			for ( int y = 0; y < Ny; y++)
			{
				std::string sub;
				iss >> sub;
				nusigf[x][y][Egrp] = ::atof(sub.c_str());
			}
		}
		nusigfsrcfile.close();
	}	
}
