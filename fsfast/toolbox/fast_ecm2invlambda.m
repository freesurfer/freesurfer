function invLambda = fast_ecm2invlambda(ecm)

ecm2 = ecm/(diag(ecm) * diag(ecm)'); %'
[u s v] = svd(ecm2);


return;
