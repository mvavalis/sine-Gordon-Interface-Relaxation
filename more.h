#ifndef __MORE_H_
#define __MORE_H_

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>

#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>

#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>

#include <dofs/dof_tools.h>

#include <fe/fe_values.h>
#include <base/quadrature_lib.h>

#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <cmath>

#include <base/logstream.h>

#include <string>
#include <sys/stat.h>

#include <base/convergence_table.h>

#include <numerics/fe_field_function.h>
#include <vector>

using namespace dealii;


/* PARAMS */
#define REFINEMENT_IN  4
#define REFINEMENT_OUT 4

#define BETA .001

#define lambda 1.

/*using constant functions L_{in} and L_{out}*/
// #define L_IN  2.2
// #define L_OUT 1.
/* ****** */

#define min(a,b) (  ((a)<(b))  ?  (a) : (b)  )

#define format__data_out(arg) \
do { \
	outfile << output_dir << filename << "."#arg; \
	std::ofstream output(outfile.str().data()); \
	data_out.write_ ## arg (output); \
} while (0);

#define TERSENOV 0

/* this class has two function that return U_{\varGamma} and u_{in}, respectively */
class InitialGuessFunction : public Function<2>
{
	private:
		std::vector< Functions::FEFieldFunction<2> * > &FEFieldFunc_in_vec;
		std::vector< Functions::FEFieldFunction<2> * > &FEFieldFunc_out_vec;
		
		Point<2> &pin1;
		Point<2> &pin2;
	public:
		unsigned int current_depth;
		
		InitialGuessFunction ( std::vector< Functions::FEFieldFunction<2> * > &FEFieldFunc_in_vec_arg,
							   std::vector< Functions::FEFieldFunction<2> * > &FEFieldFunc_out_vec_arg,
		                       Point<2> &pin1_arg, Point<2> &pin2_arg) :
		   Function<2>(),
		   FEFieldFunc_in_vec  (FEFieldFunc_in_vec_arg),
		   FEFieldFunc_out_vec (FEFieldFunc_out_vec_arg),
		   pin1 (pin1_arg), pin2 (pin2_arg),
		   current_depth (0)
		{}
		
		/*this function is called by the interpolate_boundary_values() function (returns U_{\varGamma})*/
// 		virtual double value(const Point<2> &p, const unsigned int component = 0) const
// 		{
// // 			std::cout<<"InitialGuessFunction.value -- start"<<std::endl;
// 			double val = 4. * atan(exp(p[0]));
// 			
// // 			assert(current_depth == FEFieldFunc_out_vec.size());
// 			for (unsigned int i=0; i<current_depth; i++) {
// 				Tensor<1/*rank*/, 2/*dim*/> tens_in  = FEFieldFunc_in_vec[i]->gradient(p, component);
// 				Tensor<1/*rank*/, 2/*dim*/> tens_out = FEFieldFunc_out_vec[i]->gradient(p, component);
// 				
// 				bool horiz_boundary; //whether we have to do with a horizontal (or a vertical) boundary
// 				assert(p[0]-pin1[0] >= .0);
// 				assert(pin2[0]-p[0] >= .0);
// 				assert(pin1[1]-p[1] >= .0);
// 				assert(p[1]-pin2[1] >= .0);
// 				if ( min(p[0]-pin1[0], pin2[0]-p[0])  <  min(pin1[1]-p[1], p[1]-pin2[1]) ) {
// 					horiz_boundary = false;}
// 				else {
// 					horiz_boundary = true;}
// 				
// 				/*! check out eq. (35), it doesn't matter whether we take positive or negative gradients
// 				in according to where a boundary faces, as the resulting value is squared anyway */
// 				double tmp = /* (1./L_IN )* */ (horiz_boundary ? tens_in[1] :tens_in[0] )/* tens_in  in the direction of v */ +
// 				             /* (1./L_OUT)* */ (horiz_boundary ? tens_out[1]:tens_out[0])/* tens_out in the direction of v */;
// 				
// 				val += BETA * tmp*tmp;
// 			}
// 			
// // 			std::cout<<"InitialGuessFunction.value -- end"<<std::endl;
// 			return val;
// 		}
		
// 		virtual Tensor<1,2> gradient(const Point<2> &p, const unsigned int component = 0) const
// 		{
// 			Tensor<1/*rank*/, 2/*dim*/> grad;
// 			grad[0] = (4. * exp(p[0])) / (exp(2*p[0]) + 1);
// 			grad[1] = 0.;
// 			
// 			for (unsigned int i=0; i<current_depth; i++) {
// 				double val_in  = FEFieldFunc_in_vec[i]->value(p, component);
// 				double val_out = FEFieldFunc_out_vec[i]->value(p, component);
// 				
// 				double tmp = val_in - val_out;
// 				
// 				grad[0] += BETA * tmp;
// 				//grad[0] += BETA * tmp*tmp;
// 			}
// 			
// // 			std::cout<<"$$$  "<<grad[0]<<"   "<<grad[1]<<std::endl;
// 			
// 			return grad;
// 		}
		
		//whether we have to do with a horizontal (or a vertical) boundary
		bool is_horiz_inner(const Point<2> &p) const
		{
			assert(p[0]-pin1[0] >= .0);
			assert(pin2[0]-p[0] >= .0);
			assert(pin1[1]-p[1] >= .0);
			assert(p[1]-pin2[1] >= .0);
			if ( min(p[0]-pin1[0], pin2[0]-p[0])  <  min(pin1[1]-p[1], p[1]-pin2[1]) ) {
				return false;}
			else {
				return true;}
		}
		
// 		//whether we have to do with a horizontal (or a vertical) boundary
// 		bool is_left_inner(const Point<2> &p) const
// 		{
// 			assert(p[0]-pin1[0] >= .0);
// 			assert(pin2[0]-p[0] >= .0);
// 			assert(pin1[1]-p[1] >= .0);
// 			assert(p[1]-pin2[1] >= .0);
// 			if ( p[0]-pin1[0]  <  min(min(pin1[1]-p[1], p[1]-pin2[1]), pin2[0]-p[0]) ) {
// 				return true;}
// 				else {
// 					return false;}
// 		}
// 		
// 		//whether we have to do with a horizontal (or a vertical) boundary
// 		bool is_right_inner(const Point<2> &p) const
// 		{
// 			assert(p[0]-pin1[0] >= .0);
// 			assert(pin2[0]-p[0] >= .0);
// 			assert(pin1[1]-p[1] >= .0);
// 			assert(p[1]-pin2[1] >= .0);
// 			if ( pin2[0]-p[0]  <  min(min(pin1[1]-p[1], p[1]-pin2[1]), p[0]-pin1[0]) ) {
// 				return true;}
// 				else {
// 					return false;}
// 		}
// 		
// 		//whether we have to do with a horizontal (or a vertical) boundary
// 		bool is_upper_inner(const Point<2> &p) const
// 		{
// 			assert(p[0]-pin1[0] >= .0);
// 			assert(pin2[0]-p[0] >= .0);
// 			assert(pin1[1]-p[1] >= .0);
// 			assert(p[1]-pin2[1] >= .0);
// 			if ( pin1[1]-p[1]  <  min(min(pin2[0]-p[0], p[1]-pin2[1]), p[0]-pin1[0]) ) {
// 				return true;}
// 				else {
// 					return false;}
// 		}
		
	public:
		double DELTA;
// 		double ZETA;
		virtual double value(const Point<2> &p, const unsigned int component = 0) const
		{
			double val = 4. * atan(exp(p[0]));
			
			for (unsigned int i=0; i<current_depth; i++) {
				double val_in  = FEFieldFunc_in_vec[i]->value(p, component);
				double val_out = FEFieldFunc_out_vec[i]->value(p, component);
// 				Tensor<1/*rank*/, 2/*dim*/> tens_in  = FEFieldFunc_in_vec[i]->gradient(p, component);
// 				Tensor<1/*rank*/, 2/*dim*/> tens_out = FEFieldFunc_out_vec[i]->gradient(p, component);
				
				val += DELTA * (val_in - val_out);
				
				
// 				if (is_left_inner(p)) { //left (boundary)
// 					val += DELTA * (val_out - val_in);
// 				}
// 				else {
// 					val += DELTA * (val_in - val_out);
// 				}
				
// 				if (!is_horiz_inner(p)) {
// 					val += ZETA * (tens_in[0] - tens_out[0]);
// 				}
// 				else {
// 					val += ZETA * (tens_in[1] - tens_out[1]);
// 				}
			}
			
			return val;
		}
		
		double GAMMA;
// 		double EPSILON;
		//experimental (in -- neumann)
		virtual Tensor<1,2> gradient(const Point<2> &p, const unsigned int component = 0) const
		{
			Tensor<1/*rank*/, 2/*dim*/> grad;
			grad[0] = (4. * exp(p[0])) / (exp(2*p[0]) + 1);
			grad[1] = 0.;
			
			for (unsigned int i=0; i<current_depth; i++) {
// 				double val_in  = FEFieldFunc_in_vec[i]->value(p, component);
// 				double val_out = FEFieldFunc_out_vec[i]->value(p, component);
				Tensor<1/*rank*/, 2/*dim*/> tens_in  = FEFieldFunc_in_vec[i]->gradient(p, component);
				Tensor<1/*rank*/, 2/*dim*/> tens_out = FEFieldFunc_out_vec[i]->gradient(p, component);
				
// 				grad += GAMMA * (tens_out - tens_in);
				grad += GAMMA * (tens_in - tens_out);
// 				if (is_horiz_inner(p)) { //horizontal
// 					std::cout<<"$$$ grad: "/*<< grad[0] <<", "*/<< grad[1] <<std::endl;
// 					std::cout<<"$$$ tens_out: "/*<< tens_out[0] <<", "*/<< tens_out[1] <<std::endl;
// 					std::cout<<"$$$ tens_out: "/*<< tens_out[0] <<", "*/<< tens_in[1] <<std::endl;
// 				}
				
// 				//take under consideration the normal vectors
// 				if (is_left_inner(p)) { //left (boundary)
// 					tens_in = -tens_in;
// 					grad += GAMMA * (tens_out - tens_in);
// // 					grad += GAMMA * (tens_in - tens_out);
// // 					grad += GAMMA * (tens_in + tens_out);
// // 					grad -= GAMMA * (tens_in + tens_out);
// 				}
// 				else if (is_right_inner(p)) { //right
// 					tens_out = -tens_out;
// 					grad += GAMMA * (tens_out - tens_in);
// // 					grad += GAMMA * (tens_in - tens_out);
// // 					grad += GAMMA * (tens_in + tens_out);
// // 					grad -= GAMMA * (tens_in + tens_out);
// 				}
// 				else if (is_upper_inner(p)) { //up
// 					tens_out = -tens_out;
// 					grad += GAMMA * (tens_out - tens_in);
// // 					grad += GAMMA * (tens_in - tens_out);
// // 					grad += GAMMA * (tens_in + tens_out);
// // 					grad -= GAMMA * (tens_in + tens_out);
// 				}
// 				else { //down
// 					tens_in = -tens_in;
// 					grad += GAMMA * (tens_out - tens_in);
// 				}
				
// 				if (!is_horiz_inner(p)) { //vertical
// 					grad[0] += EPSILON * (val_in - val_out);
// // 					grad[0] = 2.;
// 				}
// 				else { //horizontal
// 					grad[1] += EPSILON * (val_in - val_out);
// 					grad[1] = 2.;
// 					grad = -tens_out;
// 				}
// 				grad = -tens_out;
			}
			
			//take under consideration the normal vectors
// 			if (is_left_inner(p)) { //left (boundary)
// 				grad[1] = 0.;
// 				return -grad;
// 			}
// 			else if (is_right_inner(p)) { //right
// 				grad[1] = 0.;
// 				return grad;
// 			}
// 			else if (is_upper_inner(p)) { //up
// 				grad[0] = 0.;
// 				return grad;
// 			}
// 			else { //down
// 				grad[0] = 0.;
// 				return -grad;
// 			}
			return grad;
		}
};

#endif
