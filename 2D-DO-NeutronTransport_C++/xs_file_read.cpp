#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

//*******************************************************************************************************************************

void xs_file_read( std::vector<std::vector<std::vector<double> > > &sigt, std::string sigtfid, 
	               std::vector<std::vector<double> > &sigs, std::string sigsfid, 
	               std::vector<std::vector<std::vector<double> > > &nusigf, std::string nusigffid, int Nx, int Ny, int E )
{

	int Egrp;

	std::string sigtline;
	std::ifstream sigtsrcfile;

	Egrp = 0;

	sigtsrcfile.open(sigtfid);

	if ( sigtsrcfile.is_open() )
	{

		//for ( int x = 0; x < Nx; x++ )
		for ( int x = 0; getline(sigtsrcfile,sigtline); x++ )
		{
			//getline(sigtsrcfile,sigtline);

			if ( sigtline[0] == 'E')
			{
				Egrp = Egrp + 1 ;
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
		
//	std::string sigsline;
//	std::ifstream sigssrcfile;

	std::ifstream file(sigsfid);

	for ( int i = 0; i < E; i++ )
	{
		for ( int j = 0; j < E; j++ )
		{
			file >> sigs[i][j];
			std::cout << sigs[i][j] << std::endl;
		}
	}


//	sigssrcfile.open(sigsfid);

//	if ( sigssrcfile.is_open() )
//	{

		//for ( int x = 0; x < Nx; x++ )
//		for ( int x = 0; getline(sigssrcfile,sigsline); x++ )
//		{
//			getline(sigssrcfile,sigsline);
//
//			std::istringstream iss(sigsline);

//			for ( int y = 0; y < Egrp; y++)
//			{
//				std::string sub;
//				iss >> sub;
//				std::cout << sub << '\n' << std::endl;
//				sigs[x][y] = ::atof(sub.c_str());
//			}
//		}
//		sigssrcfile.close();
//	}


	std::string nusigfline;
	std::ifstream nusigfsrcfile;

	Egrp = 0;

	nusigfsrcfile.open(sigsfid);

	if ( nusigfsrcfile.is_open() )
	{
		//getline(nusigfsrcfile,nusigfline);

		//for ( int x = 0; x < Nx; x++ )
		for ( int x = 0; getline(nusigfsrcfile,nusigfline); x++ )
		{
			//getline(nusigfsrcfile,nusigfline);

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
