# source: https://stackoverflow.com/questions/12088080/how-to-convert-integer-number-into-binary-vector
number2binary <- function(number, noBits) 
{
	binary_vector = rev(as.numeric(intToBits(number)))
	if(missing(noBits)) 
	{
		return(binary_vector)
	} else 
	{
		binary_vector[-(1:(length(binary_vector) - noBits))]
	}
}

# call parameters
source("pars.R")

j=64 # j=64 is for state '1000000' (G1)
number2binary(j, 7) -> initial
	
# Initial condition
Cdh1=initial[1] 
SBF=initial[2]
Cln2=initial[3]
Clb5=initial[4]
Clb2G=initial[5]
Clb2M=initial[6]
Cdc20=initial[7]
Size=0.65
tr=0;
tr_max=600;
	
state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
colnames(state)[10] <- 'Phase'
			
S0=rlnorm(1,log(S0.mean), S0.CV)
	
while(tr<tr_max)
{
    	# Calculate new values by Boolean functions
    	Cdh1new = !( Cln2 || Clb5 || Clb2G || Clb2M ) || ( Cdc20 && !( Cln2 || Clb5 || Clb2M ) )
	SBFnew = ( SBF && !( Clb2G || Clb2M ) ) || 
		( 
			( Cdh1 && !SBF && !Cln2 && !Clb5 && !Clb2G && !Clb2M && !Cdc20 ) &&
		  	( Size > S0 ) && ( runif(1) < (Size-S0)^2 )
		)
	Cln2new = SBF
	Clb5new = ( Clb5 || SBF ) && !( Cdh1 || Cdc20 )
	Clb2Gnew = Clb2M || ( ( Clb5 || Clb2G ) && !Cdh1 )
	Clb2Mnew = ( Clb2G || Clb2M ) && !Cdh1 && ( !Cdc20 || Clb5 ) && !Cln2
	Cdc20new = Clb2M || ( Clb2G && Cdc20 )

    	# Which variables change? del...=0 if no change.
	delCdh1 = Cdh1new - Cdh1
	delSBF  = SBFnew - SBF
	delCln2 = Cln2new - Cln2
	delClb5 = Clb5new - Clb5
	delClb2G = Clb2Gnew - Clb2G
	delClb2M = Clb2Mnew - Clb2M
	delCdc20 = Cdc20new - Cdc20

    	# Propensity
    	x1 = abs(delCdh1)*pCdh1
	x2 = x1 + abs(delSBF)*pSBF
	x3 = x2 + abs(delCln2)*pCln2
	x4 = x3 + abs(delClb5)*pClb5
	x5 = x4 + abs(delClb2G)*pClb2G
	x6 = x5 + abs(delClb2M)*pClb2M
	x7 = x6 + abs(delCdc20)*pCdc20

	# If no variables change, x7=0 so set propensity = pG1
    	if(x7>0) { x8=x7 } else 
    	{ x8=pG1 }

	# Determine which variable actually changes?
    	s2 = runif(1)*x8
	y1 = (s2<x1)
	y2 = (s2>=x1)&&(s2<x2)
	y3 = (s2>=x2)&&(s2<x3)
	y4 = (s2>=x3)&&(s2<x4)
	y5 = (s2>=x4)&&(s2<x5)
	y6 = (s2>=x5)&&(s2<x6)
	y7 = (s2>=x6)&&(s2<x7)
	y8 = (s2>=x7)
	
	# Update real time given the probability per unit time of the possible changes
	delt = -log(runif(1))/x8
	
	# delt during Clb2M or Cdc20 activation are drawn from a lognormal distribution
	# This will overwrite above delt
	if(y6==1 && Clb2M==0)
	{
		delt <- rlnorm(1,log(tM.mean), tM.CV)
	}
	if(y7==1 && Cdc20==0)
	{
		delt <- rlnorm(1,log(tM.mean), tM.CV)
	}
	
	tr = tr + delt   
	
  	# cell division event
    	if (Clb2G==1 && Clb2Gnew==0 && y5==1)
    	{
        	f = rlnorm(1, log(f.mean), f.CV)
     
        	# size at birth of the mother cell
        	Size = Size*exp(mu*delt)*f
        	S0=rlnorm(1,log(S0.mean), S0.CV) ## set a new S0
    	} else
    	{
        	Size = Size*exp(mu*delt)
    	}	

    	# Update variables
    	Cdh1 = Cdh1 + y1*delCdh1
	SBF  = SBF  + y2*delSBF
	Cln2 = Cln2 + y3*delCln2
	Clb5 = Clb5 + y4*delClb5
	Clb2G = Clb2G + y5*delClb2G
	Clb2M = Clb2M + y6*delClb2M
	Cdc20 = Cdc20 + y7*delCdc20
	  
	new.state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
	colnames(new.state)[10] <- 'Phase'
	state <- rbind(state, new.state)	
}
	
