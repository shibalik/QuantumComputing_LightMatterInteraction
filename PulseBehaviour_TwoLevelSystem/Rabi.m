for d= -15:15
    for R=0:15
        for t=10:40
           for 0:0.01:t
               H_1=H_0+H*0.01;
               H_0=H_1
        U=exp(-iH/h);
        p= U'*p_0*U;
        P= [0,1]*p*[0;1];
          