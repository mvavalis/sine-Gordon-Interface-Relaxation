/*
 * In pseudocode:
 *
 * Init:
 *  -Set each domain's extreme points (pin and pout for the inner and outer rectangle, respectively)
 *  -Initialize the boundary function (uses http://www.dealii.org/developer/doxygen/deal.II/classFunctions_1_1FEFieldFunction.html)
 * Step 2:
 *  -Create a SineGordon_in object, then call it's run function;
 *   supply it with the boundary function.
 *  [Disregard the initial guess inside the domain for u^{(k-1)}_in, use 0 instead.]
 * Step 3:
 *  [Analogous to Step 2.]
 * Step 5:
 *  -Somehow check for convergence.
 */

#include <unistd.h>
#include <assert.h>
#include <signal.h>

#include "more.h"
#include "SineGordon_in.h"
#include "Laplace_out.h"

#define NOF_FIRST_ITERS_TO_PRINT      0
#define PRINT_FIRST_ITERS__THEN_STOP  0

#define STOP_ON_NORM_RAISE 0
double prev_norm = 10e+20;

#define MAX_ITERS 40
double norm_stats[MAX_ITERS];

//used to check for convergence
#define CHK_CONV_SIZE 100
double old_bound_vals[CHK_CONV_SIZE];
double new_bound_vals[CHK_CONV_SIZE];
#define CONV_THRESH ((1e-35) * BETA)

void set_new_bound_vals(InitialGuessFunction &IGFunc, Point<2> &pin1, Point<2> &pin2)
{
	double x, y;
	int i, cnt=0;
	
	for (i=0; i<CHK_CONV_SIZE/4; i++) {
		x = pin1[0];
		x += i * (pin2[0]-pin1[0]) / (CHK_CONV_SIZE/4);
		
		//upper boundary
		y=pin1[1];
		Point<2> tmpp1 (x,y);
		new_bound_vals[cnt] = IGFunc.value(tmpp1); //value
// 		new_bound_vals[cnt] = IGFunc.gradient(tmpp1)[1]; //gradient
		cnt++;
		
		//lower boundary
		y=pin2[1];
		Point<2> tmpp2 (x,y);
		new_bound_vals[cnt] = IGFunc.value(tmpp2); //value
// 		new_bound_vals[cnt] = IGFunc.gradient(tmpp2)[1]; //gradient
		cnt++;
	}
	for (i=0; i<CHK_CONV_SIZE/4; i++) {
		y = pin2[1];
		x += i * (pin1[1]-pin2[1]) / (CHK_CONV_SIZE/4);
		
		//left boundary
		x=pin1[0];
		Point<2> tmpp1 (x,y);
		new_bound_vals[cnt] = IGFunc.value(tmpp1); //value
// 		new_bound_vals[cnt] = IGFunc.gradient(tmpp1)[0]; //gradient
		cnt++;
		
		//right boundary
		x=pin2[0];
		Point<2> tmpp2 (x,y);
		new_bound_vals[cnt] = IGFunc.value(tmpp2); //value
// 		new_bound_vals[cnt] = IGFunc.gradient(tmpp2)[0]; //gradient
		cnt++;
	}
	assert(cnt<=CHK_CONV_SIZE);
}

void cpy_bound_vals()
{
	for (int i=0; i<CHK_CONV_SIZE; i++) {
		old_bound_vals[i] = new_bound_vals[i];
	}
}

bool chk_conv(int iter)
{
	double norm = 0;
	for (int i=0; i<CHK_CONV_SIZE; i++) {
		double tmp = old_bound_vals[i] - new_bound_vals[i];
		norm += tmp*tmp;
	}
	norm = sqrt(norm);
	std::cout<<"### chk_conv() -- norm: "<<norm<<" (CONV_THRESH: "<<CONV_THRESH<<")"<<std::endl;
	
	norm_stats[iter-1] = norm;
	
#if STOP_ON_NORM_RAISE == 1
	//! THIS IS A HACK
	if (norm > prev_norm) {return true;}
	prev_norm = norm;
#endif
	
	if (norm < CONV_THRESH) {return true;}
	else {return false;}
}

void print_output(InitialGuessFunction &IGFunc, SineGordon_in *solver_in, Laplace_out *solver_out, const char *str_append, int k, bool log)
{
	//print_tecplot
	char fname[1000];  memset(fname, '\0', 1000);
// 	sprintf(fname, "%f_%f_%f_%f", IGFunc.GAMMA, IGFunc.DELTA, IGFunc.EPSILON, IGFunc.ZETA);
	sprintf(fname, "%f_%f", IGFunc.GAMMA, IGFunc.DELTA);
	
	char tmp[1000];
	memset(tmp, '\0', 1000);  sprintf(tmp, "%s__solution-Laplace_out--iter_%d%s", fname, k, str_append);    solver_out->output_results(tmp);
	memset(tmp, '\0', 1000);  sprintf(tmp, "%s__solution-SineGordon_in--iter_%d%s", fname, k, str_append);  solver_in->output_results(tmp);
	
	//log
	if (log) {
		FILE *file = fopen("solution/output.log","a+");
		fprintf(file, "%%=============================\n");
		fprintf(file, "%% %s -- iter: %d, norm: %f\n", fname, k, prev_norm);
		int i;
		for (i=1; i<=k; i++) {fprintf(file, "%f ", norm_stats[i-1]);}
		for (; i<=MAX_ITERS; i++) {fprintf(file, "0. ");}
		fprintf(file, "\n");
		fclose(file);
	}
}

bool going_out = false;
void INThandler(int sig)
{
	signal(sig, SIG_IGN);
	going_out = true;
}

int main (int argc, char *argv[])
{
	signal(SIGINT, INThandler);
	
	/* Init */
#if TERSENOV ==1
	/*Tersenov*/
	Point<2> pout1 ( 0., 3.);
	Point<2> pout2 (12., 0.);
	Point<2> pin1 ( 1., 2.);
	Point<2> pin2 (11., 1.);
#else
	/*Mu*/
	Point<2> pout1 (-20., 10.);
	Point<2> pout2 ( 20.,-10.);
	Point<2> pin1 (-10., 2.);
	Point<2> pin2 ( 10.,-2.);
#endif
	
	SineGordon_in *old_solver_in  = NULL;
	Laplace_out   *old_solver_out = NULL;
	std::vector< Functions::FEFieldFunction<2> * > FEFieldFunc_in_vec;
	std::vector< Functions::FEFieldFunction<2> * > FEFieldFunc_out_vec;
	
	InitialGuessFunction IGFunc(FEFieldFunc_in_vec, FEFieldFunc_out_vec, pin1, pin2);
	
	if (argc != 3) {std::cout<<"Wrong input!"<<std::endl; exit(1);}
	IGFunc.GAMMA   = atof(argv[1]); std::cout<<"GAMMA   : "<<IGFunc.GAMMA<<std::endl;
	IGFunc.DELTA   = atof(argv[2]); std::cout<<"DELTA   : "<<IGFunc.DELTA<<std::endl;
// 	IGFunc.EPSILON = atof(argv[3]); std::cout<<"EPSILON : "<<IGFunc.EPSILON<<std::endl;
// 	IGFunc.ZETA    = atof(argv[4]); std::cout<<"ZETA    : "<<IGFunc.ZETA<<std::endl;
	
	set_new_bound_vals(IGFunc, pin1, pin2);
	
// #if NOF_FIRST_ITERS_TO_PRINT != 0
// 	SineGordon_in *solver_in_dummy = new SineGordon_in(IGFunc, pin1, pin2);
// 	Laplace_out *solver_out_dummy = new Laplace_out(IGFunc, pin1, pin2, pout1, pout2);
// 	print_output(IGFunc, solver_in_dummy, solver_out_dummy, "--initial", 0, false);
// #endif
	
	for (unsigned int k=1; ; k++) {
		printf("### main() -- Step 1 (init) (iter: %d)\n", k);
		assert(FEFieldFunc_in_vec.size()  == k-1);
		assert(FEFieldFunc_out_vec.size() == k-1);
		
		/* Step 2: solve inside */
		printf("### main() -- Step 2 (SineGordon in)\n");
		SineGordon_in *solver_in = new SineGordon_in(IGFunc, pin1, pin2);
		FEFieldFunc_in_vec.push_back(solver_in->run(old_solver_in));
		old_solver_in = solver_in;
		
		/* Step 3: solve outside */
		printf("### main() -- Step 3 (Laplace out)\n");
		Laplace_out *solver_out = new Laplace_out(IGFunc, pin1, pin2, pout1, pout2);
		FEFieldFunc_out_vec.push_back(solver_out->run(old_solver_out));
		old_solver_out = solver_out;
		
		IGFunc.current_depth++; //this corresponds to the nof times that the two problems have been solved
		
		/* Step 5: check for convergence */
		printf("### main() -- Step 5 (check conv.)\n");
		cpy_bound_vals();
		set_new_bound_vals(IGFunc, pin1, pin2);
		if (k == MAX_ITERS) {going_out = true;}
		
if (k<=NOF_FIRST_ITERS_TO_PRINT) {print_output(IGFunc, solver_in, solver_out, "", k, false);}
#if PRINT_FIRST_ITERS__THEN_STOP == 1
		else {break;}
#endif
		
#if STOP_ON_NORM_RAISE == 1
		if (chk_conv(k)) {going_out = true;}
#else
		if (chk_conv(k)) {break;}
#endif
		
		if (going_out) {
			print_output(IGFunc, solver_in, solver_out, "--last", k, true);
			break;
		}
	}
}
