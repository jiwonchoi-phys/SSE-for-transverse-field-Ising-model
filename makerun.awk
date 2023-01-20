BEGIN{
    L=4
    M=500
    beta=200
    Gmax=0.5
    Equilibration=1000
    binsize = 500
    repeat = 100
    sidx = 0
    printf("#!/bin/bash\n")
    for (i=0;i<10;++i){
      printf("echo %.4f\n",Gmax)
      printf("time build/2d %d %d %.4f %.4f %d %d %d %d\n",L,M,beta,Gmax,Equilibration,binsize,repeat,sidx);
      Gmax *= sqrt(2);
    }
    #printf("rm sample$idx/l$d.conf\n",L)
}
