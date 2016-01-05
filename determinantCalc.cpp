/* 
 * Tushar Maharishi
 * tm5gf
 * 10/02/2015
 * determinantCalc.cpp
 */

#include <iostream>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

int power(int mantissa, int exp)
{

  int ans = 1;
  for(int x = 0; x < exp; x++)
    ans *= mantissa;
  return ans;

}

void printMatrix(vector<double> a, int r, int c)
{

  for(int x = 0; x < r*c; x++)
    {
      if(x % c == 0)
	cout << "\n";
      cout << a[x] << "\t";
    }

  cout << "\n";

}

double determinantCalc(vector<double> a, int n)
{

  if(n == 1) // if 1x1 matrix
    return a[0];

  if(n == 2) // base case, 2x2 matrix
    return a[0] * a[3] - a[1] * a[2]; // ad - bc

  double ans = 0;

  for(int y = 0; y < n; y++) 
    { 
      vector<double> b; // creates a vector b to store the smaller matrices in
      for(int x = 0; x < n*n; x++)
	{
	  if((x)%(n) != (0+y) && x >= n) // ignores the current row and column and 
	    {
	      b.push_back(a[x]); // inserts the remaining rows/cols into the new matrix
	    }
	}
      ans += power(-1,y%2) * a[y] * determinantCalc(b, n-1); 
      // adjusts sign, multiplies cofactor, multiplies new determinant

    }

  return ans;

}

vector<double> rowSwap(vector<double> matrix, int cols, int row1, int row2) // simple swap method, adjusted to work for rows
{
  vector<double> temp; // must save entire row, not just a variable
  for(int x = 0; x < cols; x++)
    {
      temp.push_back(matrix[row1*cols + x]);
      matrix[row1*cols + x] = matrix[row2*cols + x];
      matrix[row2*cols + x] = temp[x];
    }
  return matrix; 
}

int findNonZeroPivot(vector<double> matrix, int r, int c, int col, int failures) // simple find method
{
  for(int row = col-failures; row < r; row++)
    {
      if(matrix[row*c + col + failures] != 0)
	return row;
    }
  return -1; // return -1 if failed (whole column of zeroes has been created)
}

vector<double> rowReductionCalc(vector<double> matrix, int r, int c)
{
  int successes = 0; // successful pivots found per iteration
  int failures = 0; // failed pivots found per iteration
  for(int x = 0; x < c; x++)
    {

      if(successes >= r) // if looking for another pivot when pivot count exceeds rows
	break;

      double tempPivotVal = matrix[x + (x-failures)*c]; // staircase-style pivots
      if(tempPivotVal == 0) // if the current pivot is a 0, perform a row-swap (if possible)
	{
	  int newPiv = findNonZeroPivot(matrix, r, c, x, failures); // find new pivot for row-swap
	  
	  if(newPiv >= 0) // if new pivot has been found, perform row-swap
	    {
	      matrix = rowSwap(matrix, c, newPiv, x-failures);
	    }
	  tempPivotVal = matrix[x + (x-failures)*c]; // re-evaluate the pivot value
	}

      if(tempPivotVal != 0) // must be included to eliminate divide-by-zero errors in the following for-loop
	{
	  successes++;
	  // reduce pivot-locations' values to 1
	  for(int y = (x-failures)*c; y < (x-failures)*c + c; y++) // x-failures means "row of the pivot"
	    //                                                     // x means "col of the pivot"
	    {
	      if(matrix[y] != 0) // gets rid of -0 errors
		matrix[y] /= tempPivotVal;
	    }

	  // subtract pivot row from other rows
	  for(int row = 0; row < r; row++)
	    {
	      if((x-failures) != row) // DO NOT SUBTRACT ROW FROM ITSELF!
		{
		  if(matrix[row*c + x] != 0) // make pivot the only non-zero term in the column
		    {
		      double factor = matrix[row*c + x]; // create scalar for the row combination. IE: R3 = R3 - 4R2
		      for(int col = 0; col < c; col++) // perform arithmetic on each member of the row
			{
			  matrix[row*c + col] -= factor * matrix[(x-failures)*c + col];
			}
		    }
		}
	    }
	}
      else
	failures++;
    }
  
  return matrix;

}

double determinantCalcFaster(vector<double> matrix, int n) // row reduction and multiplying across diagonal, avoids the messy recursive determinant calculator
{
  int determinant = 1;
  matrix = rowReductionCalc(matrix, n, n); // performs row reduction
  for(int i = 0; i < n; i++)
    determinant *= matrix[i*n + i]; // diagonal values
  return determinant;
}

bool isREF(vector<double> matrix, int maxRows, int maxCols)
{

  int leadingTermPosition = -1;
  int numZeroRows = 0;
  bool isZeroRow;

  for(int row = 0; row < maxRows; row++) // CHECK IF PIVOTS ARE IN CORRECT POSITION
    {
      isZeroRow = true;
      for(int col = 0; col < maxCols; col++)
	{
	  if(matrix[row*maxCols + col] != 0) // has a pivot-position been found?
	    {
	      isZeroRow = false;
	      if(col <= leadingTermPosition) // check if leading term follows staircase-format
		{
		  cout << "Value at row " << row+1 << " and column " << col+1 << " is invalid position for pivot." << endl;
		  return false;
		}
	      else
		{
		  leadingTermPosition = col; // update leading term position
		}
	      break;
	    }
	}
      if(isZeroRow) // increment number of zero rows if necessary
	numZeroRows++;
    }

  for(int row = maxRows-numZeroRows; row < maxRows; row++) // CHECK IF ZERO ROWS ARE AT THE BOTTOM
    {
      for(int col = 0; col < maxCols; col++)
	{
	  if(matrix[row*maxCols + col] != 0)
	    {
	      cout << "Zero rows are not properly placed at the bottom of the matrix." << endl;
	      return false;
	    }
	}
    }

  return true;
  
}

bool isRREF(vector<double> matrix, int maxRows, int maxCols)
{

  int leadingTermPosition = -1;
  int numZeroRows = 0;
  bool isZeroRow;

  for(int row = 0; row < maxRows; row++) // CHECK IF PIVOTS ARE IN CORRECT POSITION
    {
      isZeroRow = true;
      for(int col = 0; col < maxCols; col++)
	{
	  if(matrix[row*maxCols + col] != 0) // has a pivot-position been found?
	    {
	      isZeroRow = false;
	      if(col <= leadingTermPosition) // check if leading term follows staircase-format
		{
		  cout << "Value at row " << row+1 << " and column " << col+1 << " is invalid position for pivot." << endl;
		  return false;
		}
	      else if(matrix[row*maxCols + col] != 1) // check if pivot value is a 1
		{
		  cout << "Value at row " << row+1 << " and column " << col+1 << " is not 1. Pivots must have value 1." << endl;
		  return false;
		}
	      else
		{
		  leadingTermPosition = col; // update leading term position
		  for(int r = 0; r < maxRows; r++) // check if pivot is only non-zero value in column
		    {
		      if(r != row)
			if(matrix[r*maxCols + col] != 0)
			  {
			    cout << "Value ar row " << r+1 << " and column " << col+1 << " is not a zero. All non-pivots in a pivot column must be zero." << endl;
			    return false;
			  }
		    }
		}
	      break;
	    }
	}
      if(isZeroRow) // increment number of zero rows if necessary
	numZeroRows++;
    }

  for(int row = maxRows-numZeroRows; row < maxRows; row++) // CHECK IF ZERO ROWS ARE AT THE BOTTOM
    {
      for(int col = 0; col < maxCols; col++)
	{
	  if(matrix[row*maxCols + col] != 0)
	    {
	      cout << "Zero rows are not properly placed at the bottom of the matrix." << endl;
	      return false;
	    }
	}
    }

  return true;
  
}

bool contains(vector<int> list, int element)
{
  for(int i = 0; i < list.size(); i++)
    if(list[i] == element)
      return true;
  return false;
}

void printSolutions(vector<double> matrix, int maxRows, int maxCols) // can only be performed on row-reduced matrices
{
  bool foundPivot;
  vector<int> nonPivots;
  for(int row = 0; row < maxRows; row++) 
    {
      foundPivot = false;	
      for(int col = 0; col < maxCols-1; col++) // maxCols-1 to account for augmented matrix
	{
	  if(matrix[row*maxCols + col] != 0) // has a non-zero value been found?
	    {
	      if(!foundPivot) // is it a pivot?
		{
		  foundPivot = true;
		  cout << "x" << (col+1) << " = ";
		  cout << matrix[row*maxCols + maxCols-1]; // print out augmented value
		}
	      else // if not a pivot
		{
		  if(!contains(nonPivots, col+1)) // add to list of non-pivots
		    nonPivots.push_back(col+1);
		  cout << " + " << (-1*matrix[row*maxCols + col]) << "(x" << (col+1) << ")"; 
		}
	    }
	}
      cout << endl;
    }
  for(int i = 0; i < nonPivots.size(); i++) // print out free variables
    cout << "x" << (nonPivots[i]) << " - free" << endl;
}

bool isConsistent(vector<double> matrix, int maxRows, int maxCols) // must be used on a RREF matrix to give proper answer
{
  bool zeroRow;
  for(int row = 0; row < maxRows; row++)
    {
      zeroRow = true;
      for(int col = 0; col < maxCols-1; col++)
	{
	  if(matrix[row*maxCols + col] != 0)
	    {
	      zeroRow = false;
	      break;
	    }
	}
      if(zeroRow)
	if(matrix[row*maxCols + maxCols-1] != 0)
	  return false;
    }
  return true;
}

int numPivots(vector<double> matrix, int maxRows, int maxCols) // must be used on a RREF matrix to give proper answer
{
  int count = 0;
  for(int row = 0; row < maxRows; row++)
    {
      for(int col = 0; col < maxCols; col++)
	{
	  if(matrix[row*maxCols + col] != 0)
	    {
	      count++;
	      break;
	    }
	}
    }
  return count;
}

bool checkSpan(vector<double> matrix, int maxRows, int maxCols) 
{
  if(maxCols < maxRows)
    return false;
  
  if(numPivots(rowReductionCalc(matrix, maxRows, maxCols), maxRows, maxCols) < maxRows)
    return false;

  return true;
}

bool checkIndependence(vector<double> matrix, int maxRows, int maxCols) // returns true if linearly independent
{
  if(maxCols > maxRows)
    return false;

  if(numPivots(rowReductionCalc(matrix, maxRows, maxCols), maxRows, maxCols) < maxCols)
    return false;

  return true;
}

bool matrixAdditionValidity(int maxRows1, int maxCols1, int maxRows2, int maxCols2)
{
  if(maxRows1 != maxRows2 || maxCols1 != maxCols2)
    return false;
  return true;
}

vector<double> matrixAddition(vector<double> matrix1, int maxRows, int maxCols, vector<double> matrix2, int maxRows2, int maxCols2)
{
  for(int i = 0; i < matrix1.size(); i++) // perform the addition
    {
      matrix1[i] += matrix2[i];
    }
  return matrix1;
}

vector<double> scalarMultiplication(vector<double> matrix, double scalar)
{
  for(int i = 0; i < matrix.size(); i++)
    matrix[i] *= scalar;
  return matrix;
}

bool matrixMultiplicationValidity(int maxCols1, int maxRows2)
{
  if(maxCols1 != maxRows2)
    return false;
  return true;
}

vector<double> matrixMultiplication(vector<double> matrix1, int maxRows, int maxCols, vector<double> matrix2, int maxRows2, int maxCols2)
{
  vector<double> answer;
  double sum;
  for(int i = 0; i < maxRows*maxCols2; i++)
    {
      sum = 0;
      for(int c = 0; c < maxCols; c++)
	{
	  sum += matrix1[(i/maxCols2)*maxCols + c] * matrix2[c*maxCols2 + (i%maxCols2)];
	}
      answer.push_back(sum);
    }
  return answer;
}

vector<double> splitAugmentedMatrix(vector<double> matrix, int colStart, int rowEnd, int colEnd)
{
  vector<double> answer;
  for(int row = 0; row < rowEnd; row++)
    for(int col = colStart; col < colEnd; col++)
      answer.push_back(matrix[row*colEnd + col]);
  return answer;
}

vector< vector<double> > LDU_Factorization(vector<double> matrix, int r, int c)
{

  vector< vector<double> > answer;
  vector<double> L;
  vector<double> D;

  for(int row = 0; row < r; row++) // fill L, and D with 0's, necessary for changing/inserting elements at specific indices in these matrices
    for(int col = 0; col < r; col++)
      {
	if(row == col)
	  {
	    L.push_back(1);
	    D.push_back(1);
	  }
	else
	  {
	    L.push_back(0);
	    D.push_back(0);
	  }
      }

  int successes = 0; // successful pivots found per iteration
  int failures = 0; // failed pivots found per iteration
  for(int x = 0; x < c; x++)
    {

      if(successes >= r) // if looking for another pivot when pivot count exceeds rows
	break;

      double tempPivotVal = matrix[x + (x-failures)*c]; // staircase-style pivots
      if(tempPivotVal == 0) // if the current pivot is a 0, perform a row-swap (if possible)
	{
	  int newPiv = findNonZeroPivot(matrix, r, c, x, failures); // find new pivot for row-swap
	  
	  if(newPiv >= 0) // if new pivot has been found, LDU FACTORIZATION CANNOT BE COMPLETED AS ROW-SWAPS ARE ILLGEAL!
	    {
	      return answer; // return NULL as LDU factorization is not possible for this matrix
	    }
	}

      if(tempPivotVal != 0) // must be included to eliminate divide-by-zero errors in the following for-loop
	{
	  successes++;

	  for(int y = (x-failures); y < r; y++) // x-failures means "row of the pivot", x means "col of the pivot"
	    {
	      if(matrix[y*c + x] != 0) // gets rid of -0 errors)
		L[y*r + x-failures] = matrix[y*c + x] / tempPivotVal; // creates the L matrix
	    }

	  // subtract pivot row from other rows
	  for(int row = x-failures+1; row < r; row++) // start at the row immediately after the pivot row
	    {
	      if(matrix[row*c + x] != 0) // make pivot the only non-zero term in the column
		{
		  double factor = matrix[row*c + x] / tempPivotVal; // create scalar for the row combination. IE: R3 = R3 - 4R2
		  for(int col = 0; col < c; col++) // perform arithmetic on each member of the row
		    {
		      matrix[row*c + col] -= factor * matrix[(x-failures)*c + col]; // creates the U matrix
		    }
		}
	    }
	  int tempPivotPosition = (x-failures)*c + x;
	  D[(x-failures)*r + x-failures] = matrix[tempPivotPosition]; // creates the D matrix

	  for(int y = (x-failures)*c; y < (x-failures)*c + c; y++) // x-failures means "row of the pivot"
	    //                                                     // x means "col of the pivot"
	    {
	      if(matrix[y] != 0) // gets rid of -0 errors
		matrix[y] /= tempPivotVal;
	    }

//matrix[tempPivotPosition] = 1;
	}
      else
	{
	  failures++;
	}


    }
  
  answer.push_back(L);
  answer.push_back(D);
  answer.push_back(matrix); // U

  return answer;

}

vector<double> findBasisVectors(vector<double> matrix, vector<double> basis, int pivotCount, int maxRows, int maxCols)
{
  vector<double> answer;
  int count = 0; // represents number of basis vectors found so far
  for(int i = 0; i < maxRows * pivotCount; i++)
    answer.push_back(0); // fill answer with 0's to allow index-usage later

  for(int row = 0; row < maxRows; row++)
    {
      for(int col = 0; col < maxCols; col++)
	{
	  if(basis[row*maxCols + col] != 0) // if pivot position has been found
	    {
	      for(int r = 0; r < maxRows; r++)
		{
		  answer[r*pivotCount + count] = matrix[r*maxCols + col]; // copies column from original matrix into basis vector matrix
		}
	      count++; // increment basis vectors found
	      break;
	    }
	}
    }
  return answer;
}

bool isStochastic(vector<double> matrix, int maxRows, int maxCols)
{
  double sum;
  for(int c = 0; c < maxCols; c++)
    {
      sum = 0;
      for(int r = 0; r < maxRows; r++)
	{
	  sum += 1000*matrix[r*maxCols + c]; // *1000 included to eliminate 0.1 binary addition errors
	}
      if(sum != 1000)
	return false;
    }
  return true;
}

vector<double> findNullspaceMatrix(vector<double> matrix, int maxRows, int maxCols)
{
  vector<double> answer;
  for(int row = 0; row < maxRows; row++)
    {
      for(int col = 0; col < maxCols; col++)
	{
	  answer.push_back(matrix[row*maxCols + col]); // include all the elements from the original matrix
	  if(col == maxCols-1)
	    answer.push_back(0); // augment the matrix with 0's
	}
    }
  return answer;
}

vector<double> findInverse(vector<double> matrix, int maxRows, int maxCols)
{
  vector<double> matrixInverse;
  int invRows = maxRows;
  int invCols = maxCols*2;
  for(int r = 0; r < maxRows; r++) // creation of the "Super Augmented Matrix" - Monika Abramenko, Fall 2015
    {
      for(int c = 0; c < maxCols; c++)
	{
	  matrixInverse.push_back(matrix[r*maxCols + c]); // place original matrix elements into super augmented matrix
	  if(c == maxCols-1) 
	    {
	      for(int i = 0; i < maxCols; i++) // place the elements of the Identity Matrix into super augmented matrix
		{
		  if(i == r)
		    matrixInverse.push_back(1);
		  else
		    matrixInverse.push_back(0);
		}
	    }
	}
    }

  matrixInverse = rowReductionCalc(matrixInverse, invRows, invCols); // at this point, matrixInverse = the finished super augmented matrix
  matrixInverse = splitAugmentedMatrix(matrixInverse, maxCols, invRows, invCols); // cut off the left-identity matrix
  return matrixInverse;
}

int main()
{

  int n, rows, cols;
  n = -1; // default value for n

  cout << endl;
  cout << "You are creating a matrix." << endl;
  cout << "number of rows = ";
  cin >> rows;
  cout << "number of cols = ";
  cin >> cols;
  if(rows == cols) // only create n if rows = cols
    n = rows;
  
  cout << endl;
  cout << "Please enter each number of the matrix by row (from left to right, top to down)." << endl; // enter elements by row
  cout << "Press the enter-key after each number." << endl;

  cout << endl;
  vector<double> matrix;  
  double temp = 0;
  for(int x = 0; x < rows*cols; x++) // reading in numbers of the matrix
    {
      if(x%cols == 0)
	cout << "Row: " << (x/cols)+1 << endl;
      cin >> temp;
      matrix.push_back(temp);
    }
  
  cout << endl;
  int choice = 0;
  cout << "What would you like to do?" << endl;
  cout << "Press '0' to quit" << endl;
  cout << "Press '1' to print the matrix." << endl;
  cout << "Press '2' for row reduction." << endl;
  cout << "Press '3' to find the determinant." << endl;
  cout << "Press '4' to check if the matrix is in Row Echelon Form (REF)." << endl;
  cout << "Press '5' to check if the matrix is in Reduced Row Echelon Form (RREF)." << endl;
  cout << "Press '6' to solve the system of linear equations (This will treat the matrix as an augmented matrix)." << endl;
  cout << "Press '7' for information on span and linear independence." << endl;
  cout << "Press '8' for information on onto/one-to-one, assuming the matrix is of a linear transformation." << endl;
  cout << "Press '9' for matrix addition. You will be requested for additional information on the matrix to be added." << endl; 
  cout << "Press '10' for scalar multiplication. You will be requested for the scalar." << endl;
  cout << "Press '11' for matrix multiplication. You will be requested for additional information on the matrix to be multiplied." << endl;
  cout << "Press '12' to find the inverse." << endl;
  cout << "Press '13' for the LDU Factorization." << endl;
  cout << "Press '14' for information on Markov Chains. You will be requested the initial state and the Kth state you want to find." << endl;
  cout << "Press '15' for information on basis/dimension/rank/nullity/nullspace." << endl;
  cout << "Press '16' to find the Adjoint matrix. <CURRENTLY BUGGED!>" << endl;

  while(1)
    {
      cout << endl;
      cout << "Choice: ";
      cin >> choice;
      switch(choice)
	{
	case 0:
	  {
	    cout << "The system will now exit." << endl;
	    cout << endl;
	    exit(1);
	  }
	case 1:
	  {
	    cout << "You have chosen to print the matrix." << endl;
	    printMatrix(matrix, rows, cols);
	    break;
	  }
	case 2: 
	  {
	    cout << "You have chosen row reduction." << endl;
	    printMatrix(rowReductionCalc(matrix, rows, cols), rows, cols);
	    break;
	  }
	case 3:
	  {
	    cout << "You have chosen finding the determinant." << endl;
	    if(n == -1) // not an n x n matrix
	      cout << "Invalid dimensions for determinant calculation. Please select another option." << endl;
	    else
	      cout << "Determinant: " << determinantCalc(matrix, n) << endl;
	    break;
	  }
	case 4:
	  {
	    cout << "You have chosen to check REF." << endl;
	    if(isREF(matrix, rows, cols))
	      cout << "The matrix is in REF." << endl;
	    else
	      cout << "The matrix is NOT in REF." << endl;
	    break;
	  }
	case 5:
	  {
	    cout << "You have chosen to check RREF." << endl;
	    if(isRREF(matrix, rows, cols))
	      cout << "The matrix is in RREF." << endl;
	    else
	      cout << "The matrix is NOT in RREF." << endl;
	    break;
	  }
	case 6:
	  {
	    cout << "You have chosen to solve the system of linear equations. This assumes an augmented matrix." << endl;
	    vector<double> solution = rowReductionCalc(matrix, rows, cols);
	    if(isConsistent(solution, rows, cols))
	      printSolutions(solution, rows, cols);
	    else
	      cout << "The matrix is inconsistent. There is no solution." << endl;
	    break;
	  }
	case 7:
	  {
	    bool span;
	    bool indep;
	    span = checkSpan(matrix, rows, cols);
	    if(n == -1)
	      indep = checkIndependence(matrix, rows, cols);
	    else
	      indep = span;
	    
	    cout << "You have chosen information on span and independence." << endl;
	    
	    if(span)
	      cout << "The vectors of the matrix span R^" << rows << "." << endl;
	    else
	      cout << "The vectors of the matrix DO NOT span R^" << rows << "." << endl;
	    if(indep)
	      cout << "The vectors of the matrix are linearly independent." << endl;
	    else
	      cout << "The vectors of the matrix are linearly dependent." << endl;
	    break;
	  }
	case 8:
	  {
	    bool span;
	    bool indep;
	    span = checkSpan(matrix, rows, cols);
	    if(n == -1) // use big-theorem to avoid extra calculation
	      indep = checkIndependence(matrix, rows, cols); // must perform calculation for linear independence if row != col
	    else 
	      indep = span; // can assume independence is true/false if span is true/false due to big-theorem
	    
	    cout << "You have chosen information on onto/one-to-one." << endl;
	    
	    if(span) // span is equivalent to onto
	      cout << "The linear transformation is onto." << endl;
	    else
	      cout << "The linear transformation is NOT onto." << endl;
	    if(indep) // independence is equivalent to one-to-one
	      cout << "The linear transformation is one-to-one." << endl;
	    else
	      cout << "The linear transformation is NOT one-to-one." << endl;
	    break;
	  }
	case 9:
	  {
	    int r, c;
	    cout << "You have chosen matrix addition." << endl;
	    cout << "What are the dimensions of the matrix to be added?" << endl;
	    cout << "Number of rows: ";
	    cin >> r;
	    cout << "Number of cols: ";
	    cin >> c;
	    cout << endl;
	    cout << "Please enter each number of the matrix by row (from left to right, top to down)." << endl;
	    cout << "Press the enter-key after each number." << endl;
	    cout << endl;

	    vector<double> matrix2;  
	    double temp2 = 0;
	    for(int x = 0; x < r*c; x++) // reading in numbers of the matrix
	      {
		if(x%c == 0)
		  cout << "Row: " << (x/c)+1 << endl;
		cin >> temp2;
		matrix2.push_back(temp2);
	      }
	    if(matrixAdditionValidity(rows, cols, r, c))
	      {
		matrix2 = matrixAddition(matrix, rows, cols, matrix2, r, c);
		printMatrix(matrix2, rows, cols);
	      }
	    else
	      cout << "Invalid dimensions for matrix addition." << endl;
	    break;
	  }
	case 10:
	  {
	    double scalar;
	    vector<double> matrix2;
	    cout << "You have chosen scalar multiplication." << endl;
	    cout << "Please provide the scalar to be multiplied." << endl;
	    cout << "Scalar: ";
	    cin >> scalar;
	    matrix2 = scalarMultiplication(matrix, scalar);
	    printMatrix(matrix2, rows, cols);
	    break;
	  }
	case 11:
	  {
	    int r, c;
	    cout << "You have chosen matrix multiplication." << endl;
	    cout << "What are the dimensions of the matrix to be multiplied?" << endl;
	    cout << "Number of rows: ";
	    cin >> r;
	    cout << "Number of cols: ";
	    cin >> c;
	    cout << endl;
	    cout << "Please enter each number of the matrix by row (from left to right, top to down)." << endl;
	    cout << "Press the enter-key after each number." << endl;
	    cout << endl;

	    vector<double> matrix2;  
	    double temp2 = 0;
	    for(int x = 0; x < r*c; x++) // reading in numbers of the matrix
	      {
		if(x%c == 0)
		  cout << "Row: " << (x/c)+1 << endl;
		cin >> temp2;
		matrix2.push_back(temp2);
	      }
	    if(matrixMultiplicationValidity(cols, r))
	      {
		matrix2 = matrixMultiplication(matrix, rows, cols, matrix2, r, c);
		printMatrix(matrix2, rows, c);
	      }
	    else
	      cout << "Invalid dimensions for matrix multiplication." << endl;
	    break;
	  }
	case 12:
	  {
	    cout << "You have chosen to find the inverse." << endl;
	    if(!checkIndependence(matrix, rows, cols) && n != -1) // inverse is allowed if linear independence is allowed due to big-theorem
	      {
		cout << "There is no possible inverse for this matrix." << endl;
	      }
	    else
	      {
		printMatrix(findInverse(matrix, rows, cols), rows, cols);
	      }
	    break;
	  }
	case 13:
	  {
	    cout << "You have chosen LDU Factorization." << endl;
	    vector< vector<double> > temp;
	    temp = LDU_Factorization(matrix, rows, cols);

	    if(temp.size() == 0)
	      {
		cout << "There is no LDU Factorization for this matrix." << endl;
	      }
	    else
	      {
	        cout << "L" << endl;
		printMatrix(temp[0], rows, rows);
		cout << endl;
		cout << "D" << endl;
		printMatrix(temp[1], rows, rows);
		cout << endl;
		cout << "U" << endl;
		printMatrix(temp[2], rows, cols);
		cout << endl;
	      }
	    break;
	  }
	case 14:
	  {
	    cout << "You have chosen information on Markov Chains." << endl;
	    if(!isStochastic(matrix, rows, cols)) // checks if the matrix is legal (stochastic) or not for markov chains
	      {
		cout << "The matrix provided is not stochastic. Information on Markov Chains cannot be performed without a stochastic matrix." << endl;
		break;
	      }
	    cout << "Please enter each number of the initial state." << endl;
	    cout << "Press the enter-key after each number." << endl;

	    vector<double> curr;  
	    double temp2 = 0;
	    for(int x = 0; x < rows; x++) // reading in numbers of the initial matrix
	      {
		cout << "Row " << x+1 << ": ";
		cin >> temp2;
		curr.push_back(temp2);
	      }
	    cout << endl;

	      int kth;
	      cout << "Please enter the Kth state you are looking for: "; // requesting the kth state
	      cin >> kth;
	      cout << endl;

	      vector<double> prev;
	      int counter = 0;
	      bool kthfound = false;
	      while(curr != prev)
		{
		  prev = curr;
		  curr = matrixMultiplication(matrix, rows, cols, curr, cols, 1);
		  counter++;
		  if(kth == counter)
		    {
		      kthfound = true;
		      cout << "The Kth state: " << endl;
		      printMatrix(curr, cols, 1);
		    }
		}
	      if(!kthfound)
		{
		  cout << "The Kth state: " << endl;
		  printMatrix(curr, cols, 1);
		}
	      cout << endl;
	      cout << "The steady state vector: " << endl;
	      printMatrix(curr, cols, 1);
	      
	    break;
	  }
	case 15:
	  {
	    cout << "You have chosen to find the basis/dimension/rank/nullity/nullspace." << endl;
	    vector<double> reduced = rowReductionCalc(matrix, rows, cols);
	    int pivotCount = numPivots(reduced, rows, cols);
	    vector<double> basis = findBasisVectors(matrix, reduced, pivotCount, rows, cols);
	    cout << endl;
	    cout << "Basis: ";
	    printMatrix(basis, rows, pivotCount);
	    cout << "Dimension: " << pivotCount << endl;
	    cout << "Rank: " << pivotCount << endl;
	    int nullity = cols-pivotCount;
	    cout << "Nullity: " << nullity << endl;
	    if(nullity > 0)
	      {
		cout << "Nullspace: " << endl;
		printSolutions(findNullspaceMatrix(reduced, rows, cols), rows, cols+1);
	      }
	    else
	      cout << "There is no nullspace for this subspace." << endl;
	    break;
	  }
	case 16: 
	  {
	    cout << "You have chosen to find the Adjoint matrix." << endl;
	    if(!checkIndependence(matrix, rows, cols) && n != -1) // inverse is allowed if linear independence is allowed due to big-theorem
	      {
		cout << "There is no adjoint for this matrix." << endl;
	      }
	    else
	      {
		printMatrix(scalarMultiplication(findInverse(matrix, rows, cols), determinantCalc(matrix, n)), rows, cols); // adjoint = determinant * inverse
	      }
	    break;
	  }
	default:
	  {
	    cout << "You have chosen an invalid request. Please choose a valid option." << endl;
	  }
	}

      cout << endl << "================================================================================" << endl << endl;
      cout << "Press '0' to quit" << endl;
      cout << "Press '1' to print the matrix." << endl;
      cout << "Press '2' for row reduction." << endl;
      cout << "Press '3' to find the determinant." << endl;
      cout << "Press '4' to check if the matrix is in Row Echelon Form (REF)." << endl;
      cout << "Press '5' to check if the matrix is in Reduced Row Echelon Form (RREF)." << endl;
      cout << "Press '6' to solve the system of linear equations (This will treat the matrix as an augmented matrix)." << endl;
      cout << "Press '7' for information on span and linear independence." << endl;
      cout << "Press '8' for information on onto/one-to-one, assuming the matrix is of a linear transformation." << endl;
      cout << "Press '9' for matrix addition. You will be requested for additional information on the matrix to be added." << endl;
      cout << "Press '10' for scalar multiplication. You will be requested for the scalar." << endl;
      cout << "Press '11' for matrix multiplication. You will be requested for additional information on the matrix to be multiplied." << endl;
      cout << "Press '12' to find the inverse." << endl;
      cout << "Press '13' for the LDU Factorization." << endl;
      cout << "Press '14' for information on Markov Chains. You will be requested the initial state and the Kth state you want to find." << endl;
      cout << "Press '15' for information on basis/dimension/rank/nullity/nullspace." << endl;
      cout << "Press '16' to find the Adjoint matrix. <CURRENTLY BUGGED!>" << endl;

    }
  
  return 0;
}
