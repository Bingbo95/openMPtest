//// openMPtest.cpp : This file contains the 'main' function. Program execution begins and ends there.
////
//
//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <omp.h>
//
//using namespace std;
//void Test(int n) {
//	for (int i = 0; i < 10000; ++i) {
//		//do nothing, just waste time.
//	}
//	printf("%d, ", n);
//
//}
//
//int main(int argc, char* argv[])
//{
//	cout << omp_get_num_procs() << endl;
//	int i, j;
//	cin >> j;
//	if (j == 0) {
//#pragma omp parallel
//		{
//			cout << "test" << endl;
//		}
//#pragma omp parallel for
//		for (i = 0; i < 10; ++i)
//			Test(i);
//		//system("pause");
//		return 0;
//	}
//	else {
//#pragma omp parallel
//		{
//			cout << "OK" << endl;
//		}
//#pragma omp parallel for
//		for (i = 0; i < 10; ++i)
//			Test(i);
//		//system("pause");
//		return 0;
//	}
//
//
//
//}
//
//// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
//// Debug program: F5 or Debug > Start Debugging menu
//
//// Tips for Getting Started: 
////   1. Use the Solution Explorer window to add/manage files
////   2. Use the Team Explorer window to connect to source control
////   3. Use the Output window to see build output and other messages
////   4. Use the Error List window to view errors
////   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
////   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
