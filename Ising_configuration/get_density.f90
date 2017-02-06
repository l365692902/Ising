    subroutine get_density(num_T,m,TER,N)

    implicit none

    integer::i,j,k,start,match_density
    integer::num_T,m,N
    integer::match(num_T,2)
    real(8):: TER(num_T,m,3)
    real(8)::density(m,4)





    !find match points and record them in match array

    do j=1,m
        if(TER(1,j,2)==TER(2,1,2))then
            start=j
            exit
        else
        end if
    end do

    do j=start,m
        if(TER(1,j,3)<=TER(2,j-start+1,3) .and.TER(1,j,3)/=0. )then
            match(1,1)=j
            match(1,2)=j-start+1
            write(*,*) match(1,1),TER(1,match(1,1),2), match(1,2),TER(2,match(1,2),2)
            exit
        else
        end if
    end do


    do i=2,num_T-1
        do j=1,m
            if(TER(i+1,j,2)==TER(i,match(i-1,2),2))then
                start=j
                exit
            else
            end if
        end do
        do j=start,m
            if(TER(i,match(i-1,2)+j-start,3)<=TER(i+1,j,3) .and.TER(i,match(i-1,2)+j-start,3)/=0. )then
                match(i,1)=match(i-1,2)+j-start
                match(i,2)=j
                !write(*,*) TER(i,match(i,1),2),TER(i+1,match(i,2),2)
                exit
            else
            end if
        end do
        write(*,*) match(i,1),TER(i,match(i,1),3), match(i,2),TER(i+1,match(i,2),3)

    end do



    density(1,1)=TER(1,1,1)
    density(1,2)=TER(1,1,2)
    density(1,3)=2.
    density(1,4)=TER(1,1,3)



    k=2
    do i=1,num_T-1
        !write(*,*)match(i,1),match(i-1,2)
        do j=match(i,2),match(i+1,1)
            if(j==match(i,2))then
                match_density=k-1
            else
            end if
            if(TER(i,j,3)==0.)then
                cycle
            end if
            density(k,1)=TER(i,j,1)
            density(k,2)=TER(i,j,2)
            density(k,3)=density(match_density,3)+log(TER(i,j,3)/TER(i,match(i,1),3))-(1./TER(i,match(i,1),1))*(TER(i,match(i,1),2)-TER(i,j,2))
            density(k,4)=TER(i,j,3)
            k=k+1
        end do

    end do







    do i=1,k-1
        write(*,"(2(F7.3),1x,F14.3,F10.6)") density(i,1), density(i,2), density(i,3),density(i,4)
    end do

    !***************************************************************************************************************************
    !******************************************************************************************************************************

    return
    end subroutine