#include "test_utils.h"
#include <iostream>
#include <stdexcept>

using namespace std;

void test_1 ()
{
    verify (true);
}

int main (int argc, char **argv)
{
    try
    {
        int *p = new int;
        test_1 ();
        return 0;
    }
    catch (const exception &e)
    {
        cerr << e.what () << endl;
        return -1;
    }
}
