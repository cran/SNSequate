depan = function(x)
{
	den=(x>=-1 & x<=1)*(3/4)*(1-x^2)
	return(den)
}

pepan = function(x)
{
	prob=((x>=-1 & x<=1)*((3*x-x^3+2)/4)) + (1*(x>1))
	return(prob)
	
}