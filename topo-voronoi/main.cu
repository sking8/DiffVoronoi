//#####################################################################
// Main
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <iostream>
#include "ParseArgs.h"
#include "SoftBodyLinearFemGridDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
    constexpr int d=2;

    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","output","output path");
    parse_args.Add_Integer_Argument("-s",64,"resolution");
    parse_args.Add_Integer_Argument("-test",2,"test");
	parse_args.Add_Integer_Argument("-driver",1,"driver");
	parse_args.Add_Integer_Argument("-lf",200,"last frame");
	parse_args.Add_Integer_Argument("-fr", 25, "frame rate");
	parse_args.Add_Double_Argument("-ts", (real)0.02, "time step size");
	parse_args.Add_Double_Argument("-dp", (real)0.1, "damping");
	parse_args.Add_Option_Argument("-qs", "use quasi-static simulation");
    parse_args.Parse(argc,argv);

    std::string output_dir=parse_args.Get_String_Value("-o");
    const int scale=parse_args.Get_Integer_Value("-s");
	const int driver=parse_args.Get_Integer_Value("-driver");
	const int test=parse_args.Get_Integer_Value("-test");
	const int last_frame=parse_args.Get_Integer_Value("-lf");
	const int frame_rate = parse_args.Get_Integer_Value("-fr");
	const real time_step = parse_args.Get_Double_Value("-ts");
	const bool use_quasi_static = parse_args.Get_Option_Value("-qs");
	const real damping = parse_args.Get_Double_Value("-dp");

	switch(driver){
	case 1:{
		SoftBodyLinearFemGridDriver<d> driver;
		driver.scale=scale;
		driver.output_dir=output_dir;
		driver.test=test;
		driver.last_frame=last_frame;
		driver.Initialize();
		driver.Run();	
	}break;
	}
}

#endif