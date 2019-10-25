#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/vec.hpp>
#include <boost/qvm/quat.hpp>
#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/mat_traits_array.hpp>
#include <boost/qvm/vec_traits_array.hpp>
//#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/map_vec_mat.hpp>
#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/mat_access.hpp>
#include <fftw3.h>
#include "T2Decay.h"
#include "cmsBField.h"
#include <complex.h>
#include <cerf.h>

class cmMin{
private:
	std::vector<double> gstep;
	std::vector<double> step;
	public:
		std::vector<double> (*functor)(std::vector<double>);
		cmMin(std::vector<double> (*pfun)(std::vector<double>),std::vector<double> stp,std::vector<double> gstp){
			functor=pfun;
			gstep=gstp;
			step=stp;
		};

		std::vector<double> df(std::vector<double> in);
		std::vector<double> d2f(std::vector<double> in);


};

