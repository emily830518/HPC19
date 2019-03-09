# include <cstdlib>
# include <iostream>

using namespace std;

int main ( );
void f ( int n );


int main ( )
{
  int n = 10;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  C++ version.\n";
  cout << "  A sample code for analysis by VALGRIND.\n";

  f ( n );

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Normal end of execution.\n";

  return 0;
}

void f ( int n )

//****************************************************************************80
//
//  Comments:
//  This is the modifed version that fixed memory bugs. 
//
//
{
  int i;
  int *x;

  x = ( int * ) malloc ( (n+1) * sizeof ( int ) ); // modifying this from n to (n+1)
                                                   // This is the main memory bugs.
                                                   // We should add 1 more room because x array has index from x[0],x[1]....,x[10],
                                                   // there are 11 elements in total. Thus, n should plus 1.

  x[0] = 1;
  cout << "  " << 0 << "  " << x[0] << "\n";

  x[1] = 1;
  cout << "  " << 1 << "  " << x[1] << "\n";

  for ( i = 2; i <= n; i++ )
  {
    x[i] = x[i-1] + x[i-2];
    cout << "  " << i << "  " << x[i] << "\n";
  }

  free(x); // modifying this from delete[] to free. 
           // This is a mismatch problem: malloc is paired with free and new is paired with delete.

  return;
}