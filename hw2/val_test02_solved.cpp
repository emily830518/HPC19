# include <cstdlib>
# include <iostream>

using namespace std;

void junk_data ( );
int main ( );


int main ( )
{
  cout << "\n";
  cout << "TEST02:\n";
  cout << "  C++ version\n";
  cout << "  A sample code for analysis by VALGRIND.\n";

  junk_data ( );

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Normal end of execution.\n";

  return 0;
}

void junk_data ( )
//****************************************************************************80
//
//  Comments:
//  This is the modifed version that fixed memory bugs. 
//
//
{
  int i;
  int *x;

  x = new int[10](); // modifying from new int[10] to new int[10]()
                     // new int[10]() is a way to initialize all elements to 0

//
//  X = { 0, 1, 2, 3, 4, 0, 0, 0, 0, 0 }.
//
  for ( i = 0; i < 5; i++ )
  {
    x[i] = i;
  }
//
//  Copy some values.
//  X = { 0, 1, 0, 3, 4, 0, 0, 0, 0, 0 }.
//
  x[2] = x[7];
  x[5] = x[6];

  for ( i = 0; i < 10; i++ )
  {
    x[i] = 2 * x[i];
  }
//
//  Print X.
//
  for ( i = 0; i < 10; i++ )
  {
    cout << "  " << i << "  " << x[i] << "\n";
  }

  delete [] x;

  return;
}
