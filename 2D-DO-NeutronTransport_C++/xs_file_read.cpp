#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

//*******************************************************************************************************************************

void xs_file_read( std::vector<std::vector<double> > &sigt, std::string sigtfid, 
	               std::vector<std::vector<double> > &sigs, std::string sigsfid, int Nx, int Ny )
{

	std::string sigtline;
	std::ifstream sigtsrcfile;

	sigtsrcfile.open(sigtfid);

	if ( sigtsrcfile.is_open() )
	{
		for ( int x = 0; getline(sigtsrcfile,sigtline); x++ )
		{
			int y = 0;

			std::istringstream iss(sigtline);

			while( iss )
			{
				std::string sub;
				iss >> sub;
				sigt[x][y] = ::atof(sub.c_str());
				y++;
			}
		}
		sigtsrcfile.close();
	}

	std::string sigsline;
	std::ifstream sigssrcfile;

	sigssrcfile.open(sigsfid);

	if ( sigssrcfile.is_open() )
	{
		for ( int x = 0; getline(sigssrcfile,sigsline); x++ )
		{
			int y = 0;

			std::istringstream iss(sigsline);

			while( iss )
			{
				std::string sub;
				iss >> sub;
				sigs[x][y] = ::atof(sub.c_str());
				y++;
			}
		}
		sigssrcfile.close();
	}	
}
