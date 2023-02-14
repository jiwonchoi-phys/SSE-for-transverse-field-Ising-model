BEGIN{
    L=4
    M=500
    beta=0.2
    Gmax=3.0440
    Equilibration=1000
    binsize = 100
    repeat = 500
    sidx = 1
    printf("#!/bin/bash\n")
    for (i=0;i<15;++i){
      printf("time build/2d %d %d %.4f %.4f %d %d %d %d\n",L,M,beta,Gmax,Equilibration,binsize,repeat,sidx);
      beta *= sqrt(2);
    }
    #printf("rm sample$idx/l$d.conf\n",L)
}
