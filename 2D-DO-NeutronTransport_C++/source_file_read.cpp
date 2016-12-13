#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>

//*******************************************************************************************************************************

void source_file_read( std::vector<std::vector<std::vector<double> > > &S, std::string fid, int Nx, int Ny, std::string calc_mode )
{

	int Egrp;

	std::string line;
	std::ifstream srcfile;

	Egrp = 0;

	srcfile.open(fid);

	if ( srcfile.is_open() )
	{
		//getline(srcfile,line);

		//for ( int x = 0; x < Nx; x++ )
		for ( int x = 0; getline(srcfile,line); x++ )
		{
			//getline(srcfile,line);

			if ( line[0] == 'E')
			{
				Egrp = Egrp + 1;
				getline(srcfile,line);
				x = 0;
			}

			std::istringstream iss(line);

			for ( int y = 0; y < Ny; y++)
			{
				std::string sub;
				iss >> sub;
				S[x][y][Egrp] = ::atof(sub.c_str());
			}
		}
		srcfile.close();
	}

	std::cout << "READING SOURCE FILE" << '\n';
	//std::cout << (calc_mode.compare("kEigen") == 1) << '\n';

	//if ( calc_mode.compare("kEigen") == 1 )
	//{
	//	for ( int k = 0; k < Egrp; k++ )
	//	{
	//		for ( int i = 0; i < Nx; i++ )
	//		{
	//			for ( int j = 0; j < Ny; j++ )
	//			{
	//				S[i][j][k] = 0;
	//				std::cout << S[i][j][k] << '\n';
					//std::cout << std::setw(Nx) << S[i][j][k] << " ";
					//if ( ( j + 1 ) % Nx == 0)
					//{
					//	std::cout << std::endl;
					//}
				//}
			//}
		//}
	//}
	
}