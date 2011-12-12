#include <iostream>
using namespace std;

double GD[7] = {.05, .03, .02, .01, .005, .001, .0001}; //\gamma and \delta
int GDsize = 7;
double EZ[7] = {.05, .03, .02, .01, .005, .001, .0001}; //\epsilon and \zeta
int EZsize = 7;
// double EZ[1] = {.0}; //\epsilon and \zeta
// int EZsize = 1;

int main()
{
	system("rm solution/output.log");
	
	char cmd[1000];
	for (int i=0; i<GDsize; i++) {
		for (int j=0; j<EZsize; j++) {
			cout<<"#*# i: "<<i<<", j: "<<j<<endl;
			memset(cmd, '\0', 1000);
			sprintf(cmd, "./main.out %f %f %f %f", GD[i], GD[i], EZ[j], EZ[j]);
			system(cmd);
		}
	}
	
	return 0;
}
