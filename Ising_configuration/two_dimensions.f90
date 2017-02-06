    !Purpose:    / This program is used to obtain the average value of several quantities in 2D Ising model which adopts wolff update /
    !method:    /Basic method:the metroplice algrithm.In the update,the wolff cluster update has taken the acceptance-rejection rate into
    !            consideration.Therefore if we run Markov chain long enough, the samples, produced by the chain can be
    !            regarded as approximately following the target distribution. Hence we can choose m
    !            points from the chain to calculate the average value./
    !Introduction:/In this program ,i is a Loop Variable;
    !             By calling subroutine Get_confugeration we get the Markov chain and calculate the middle quantities which canbe used in
    !             subroutine date_calculate.
    !             m is the number of points we choose.
    !             nequi is the number of points we skippde at the begining of Markov chain.
    !             n_0 is a interval between the two adjacent from m.
    !             num is the number of points each Markov chain has.
    !             L is the number of Markov chain;
    !             N is the number of grid points of each dimension in 2D ising model.
    !             T is the temperature and Hz is the external magnetic field.
    !
    !          subroutine date_caculation gives the variance of H,M Cv,X .
    !          Hz,T, N,m,num,n_0,nequi is set as commom variables /
    !Author:   / shiqi Liu/
    !Date:     3/1 2016
    program two_dimension
    implicit none
    integer(8)::i
    real(8)::Hz,t1,t2
    integer::N,num,n_0,m,nequi,L,num_T
    common/group1/ Hz,num_T
    common/group2/ N,num,n_0,m,nequi,L

    N=32
    nequi=10000_8
    m=2000_8
    n_0=10_8
    num=m*n_0+Nequi+1
    num_T=15
    !   T=0.1
    Hz=0.0
    L=1

    call RANDOM_SEED()
    call CPU_TIME(t1)
    write(*,*) "start"
    call Get_confugeration()
    call CPU_TIME(t2)
    write(*,*) t2-t1

    stop
    end