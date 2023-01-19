BEGIN{
    L=4
    M=500
    beta=0.01
    Gmax=3
    Equilibration=1000
    binsize = 100
    repeat = 100
    sidx = 0
    printf("#!/bin/bash\n")
    for (i=0;i<12;++i){
      printf("echo %.4f\n",beta)
      printf("./2d %d %d %.4f %.4f %d %d %d %d\n",L,M,beta,Gmax,Equilibration,binsize,repeat,sidx);
      beta *= sqrt(2);
    }
    #printf("rm sample$idx/l$d.conf\n",L)
}
