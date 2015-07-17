clear 
close
clc
syms v(t) h(t) a betaV betaH MCV MCH
Sol = dsolve(diff(v) == a*MCV*exp(0*(v/MCV-h/MCH)), diff(h) == a*MCH*exp(0*(h/MCH-v/MCV)),v(0)==1.3*MCV,h(0)==1.2*MCH);