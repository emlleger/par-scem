## Divide pigeons across pigeonholes as fairly as possible.

##  Maintainer: Breanndán Ó Nualláin <bon@science.uva.nl>  

## [lhs rhs] = divvy(n,m,i) divides n objects among m containers and
## returns the interval of n which is assigned to the ith container. 
## All arguments should be positive integers.

## Need to add a condition for m>n.

function [lhs rhs] = divvy(n,m,i)
  if i>m
    error(["Index ", num2str(i), " exceeds number of containers ",...
	   num2str(m)]);
  else
    if m>n
      lhs=i;
      rhs=i;
    else
      per = floor(n/m);  # How many per pigeonhole.
      extras = rem(n,m); # How many pigeonholes get one extra pigeon.
      if (i > extras)
	lhs = 1 + ((i - 1) * per) + extras;
      else
	lhs = 1 + ((i - 1) * (per + 1));
      endif
      if (i > extras)
	inc = 0;
      else
	inc = 1;
      endif
      rhs = min(n, lhs + per + inc - 1);
    endif
  endif
endfunction

# eof.
