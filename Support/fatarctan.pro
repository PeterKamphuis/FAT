FUNCTION FATARCTAN, X, P
  ;C=abs(X[fix(n_elements(x)/3.*1.5)]-X[fix(n_elements(x)/3.)])
  C=X[n_elements(X)-1]*0.1
  RETURN, -1.*ATAN((X-P[0])/(C+abs(P[1])))/!pi*abs(P[2])+P[3] 
END
