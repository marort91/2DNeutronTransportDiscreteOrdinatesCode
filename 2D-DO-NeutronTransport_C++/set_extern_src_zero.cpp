#include <vector>

//*******************************************************************************************************************************

void set_extern_src_zero( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > &S )

//*******************************************************************************************************************************
//
//  SET_EXTERN_SRC_ZERO - External Source Set to Zero for 2D Neutron Transport Solver Criticality Calculations
//
//  Licensing:
//             The MIT License (MIT)
//
//			   Copyright (c) [2016] [Mario I. Ortega]
//
//             Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
//             associated documentation files (the "Software"), to deal in the Software without restriction, including
//             without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
//             copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the 
//             following conditions:
//
//             The above copyright notice and this permission notice shall be included in all copies or substantial 
//             portions of the Software.
//
//             THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
//             LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
//             IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
//             WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
//             OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//  Author: 
//             Mario I. Ortega, University of California, Berkeley
//
//  Modified: 
//             2 May 2016
//
//  Parameters: 
//             Input, int Nx, mesh size in x-direction.
//			   Input, int Ny, mesh size in y-direction.
//			   Input, int Egrp, number of energy groups.
//			   Input, vector<vector<double> > &S, pointer to scalar flux scattering source.
{

	for ( int i = 0; i < Egrp; i++ )
	{
		for ( int j = 0; j < Nx; j++ )
		{
			for ( int k = 0; k < Ny; k++ )
			{
				S[j][k][i] = 0;
			}
		}
	}
}