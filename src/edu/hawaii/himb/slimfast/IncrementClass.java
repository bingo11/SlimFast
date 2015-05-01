package edu.hawaii.himb.slimfast;

import org.la4j.matrix.functor.MatrixFunction;

public class IncrementClass implements MatrixFunction {

	@Override
	public double evaluate(int arg0, int arg1, double arg2) {
		// TODO Auto-generated method stub
		arg2++;
		return arg2;
	}
	
	

}
