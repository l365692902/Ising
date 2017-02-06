    subroutine Get_confugeration()
    implicit none
    integer::i,j,k,circle,d,before,after,big_L,E_points
    real(8)::Hz,T(num_T),T0
    integer::N,num,n_0,m,nequi,L,num_T
    real(8)::rand,P
    integer::line,colume
    real(8)::E_max,E_min,E_every,M_every,temper
    real(8)::U(L,num_T)
    real(8),allocatable::S(:,:),origin(:,:)
    real(8),allocatable::M_(:),H(:)
    integer,allocatable::change_line(:),change_colume(:)
    integer(8),allocatable::bond_line(:,:),bond_colume(:,:)
    real(8),allocatable::TER(:,:,:)
    real(8),allocatable::P_(:),E(:)

    common/group1/ Hz,num_T
    common/group2/ N,num,n_0,m,nequi,L

    allocate(S(N+1,N+1))
    allocate(origin(N+1,N+1))
    allocate(bond_line(N,N))
    allocate(bond_colume(N,N))
    allocate(change_line(N**3))
    allocate(change_colume(N**3))
    allocate(M_(num))
    allocate(H(num))
    allocate(TER(num_T,m,3))




    T=(/0.1,0.2,0.35,0.5,0.7,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,1000000./)

    !write(*,*)P
    do big_L=1,L
        !generate the 2D confugeration at random
        !*****************************************
        open(unit=10,file="text",status="replace")
        do circle=1,num_T
            T0=T(circle)
            write(*,*) T0

            p=1.-exp(-(1.0/T0/2.0))
            S=0.
            do i=1,N
                do j=1,N
                    call RANDOM_NUMBER(rand)
                    if(rand>0.5)then
                        S(i,j)=0.5
                    else
                        S(i,j)=-0.5
                    end if
                end do
            end do

            ! do i=1,N
            !       write(*,"(8(F5.1))")  S(i,1:N)
            !end do

            S(N+1,1:N)=S(1,1:N)
            S(1:N,N+1)=S(1:N,1)


            !************************wolff update_______wolff update***************************
            do k=1,num
                ! choose a grid point at random,make it the center of the cluster.
                !**********************************************************************************
                call RANDOM_NUMBER(rand)
                line=ceiling(rand*N)
                call RANDOM_NUMBER(rand)
                colume=ceiling(rand*N)
                bond_line(:,:)=0
                bond_colume(:,:)=0
                change_line(:)=0
                change_colume(:)=0
                !write(*,*) "i=",line,"j=",colume

                ! generate the bond( if(S(i,j)=S(i,j+1))then make bond_line(i,j)=1 at the rate of P:  if(S(i,j)==S(i+1,j))then make bond_colume(i,j)=1
                ! at the rate of P),otherwise the bond=0)
                do i=1,N
                    do j=1,N
                        if(S(i,j)==S(i,j+1))then
                            call RANDOM_NUMBER(rand)
                            if(rand<P)then
                                bond_line(i,j)=1
                            end if
                        end if
                    end do
                end do

                do i=1,N
                    do j=1,N
                        if(S(i,j)==S(i+1,j))then
                            call RANDOM_NUMBER(rand)
                            if(rand<P)then
                                bond_colume(i,j)=1
                            end if
                        end if
                    end do
                end do


                !check the first chosen point's around
                change_line(1)=line
                change_colume(1)=colume
                i=1
                j=1
                if(bond_line(change_line(i),change_colume(i))==1)then
                    bond_line(change_line(i),change_colume(i))=2
                    if(change_colume(i)==N)then
                        change_colume(j+1)=1
                    else
                        change_colume(j+1)=change_colume(i)+1
                    end if
                    change_line(j+1)=change_line(i)
                    j=j+1
                    !write(*,*)"右"
                end if



                if(bond_colume(change_line(i),change_colume(i))==1)then
                    bond_colume(change_line(i),change_colume(i))=2
                    if(change_line(i)==N)then
                        change_line(j+1)=1
                    else
                        change_line(j+1)=change_line(i)+1
                    end if
                    change_colume(j+1)=change_colume(i)
                    j=j+1
                    !write(*,*)"下"
                end if



                if(change_line(i)==1)then
                    if(bond_colume(n,change_colume(i))==1)then
                        bond_colume(n,change_colume(i))=2
                        change_line(j+1)=n
                        change_colume(j+1)=change_colume(i)
                        j=j+1
                        !write(*,*)"上"
                    end if
                else
                    if(bond_colume(change_line(i)-1,change_colume(i))==1)then
                        bond_colume(change_line(i)-1,change_colume(i))=2
                        change_line(j+1)=change_line(i)-1
                        change_colume(j+1)=change_colume(i)
                        j=j+1
                        !write(*,*)"上"
                    end if
                end if



                if(change_colume(i)==1)then
                    if(bond_line(change_line(i),n)==1)then
                        bond_line(change_line(i),n)=2
                        change_line(j+1)=change_line(i)
                        change_colume(j+1)=n
                        j=j+1
                        !write(*,*)"左"
                    end if
                else
                    if(bond_line(change_line(i),change_colume(i)-1)==1)then
                        bond_line(change_line(i),change_colume(i)-1)=2
                        change_line(j+1)=change_line(i)
                        change_colume(j+1)=change_colume(i)-1
                        j=j+1
                        !write(*,*)"左"
                    end if
                end if





                !*****************************************************
                !check the cluster outer and record their posotions.
                if(j>=2)then
                    before=2
                    after=j
100                 do i=before,after

                        if(bond_colume(change_line(i),change_colume(i))==1)then
                            bond_colume(change_line(i),change_colume(i))=2

                            if(change_line(i) == N) then
                                change_line(j+1) = 1
                            else
                                change_line(j+1)=change_line(i)+1
                            end if

                            change_colume(j+1)=change_colume(i)
                            j=j+1
                        end if




                        if(bond_line(change_line(i),change_colume(i))==1)then

                            bond_colume(change_line(i),change_colume(i))=2
                            if(change_colume(i)==N)then
                                change_colume(j+1)=1
                            else
                                change_colume(j+1)=change_colume(i)+1
                            end if
                            change_line(j+1)=change_line(i)
                            j=j+1
                        end if



                        if(change_line(i)==1)then
                            if(bond_colume(n,change_colume(i))==1)then
                                bond_colume(n,change_colume(i))=2
                                change_line(j+1)=n
                                change_colume(j+1)=change_colume(i)
                                j=j+1
                            end if
                        else

                            if(bond_colume(change_line(i)-1,change_colume(i))==1)then
                                bond_colume(change_line(i)-1,change_colume(i))=2
                                change_line(j+1)=change_line(i)-1
                                change_colume(j+1)=change_colume(i)
                                j=j+1
                            end if
                        end if



                        if(change_colume(i)==1)then
                            if(bond_line(change_line(i),n)==1)then
                                bond_line(change_line(i),n)=2
                                change_colume(j+1)=n
                                change_line(j+1)=change_line(i)
                                j=j+1
                            end if
                        else
                            if(bond_line(change_line(i),change_colume(i)-1)==1)then

                                bond_line(change_line(i),change_colume(i)-1)=2
                                change_line(j+1)=change_line(i)
                                change_colume(j+1)=change_colume(i)-1
                                j=j+1
                            end if
                        end if

                    end do
                    before=after+1
                    after=j
                    if(before<=after)then
                        goto 100
                    else
                    end if
                else
                end if





                !************cluster_turn_over************************************
                ! write(*,*) "turn it over"
                origin(:,:)=S(:,:)
                S(line,colume)=-origin(line,colume)
                do i=1,N
                    do j=1,N
                        if(bond_line(i,j)==2)then
                            !write(*,*)"bond_line(",i,j,")"
                            !pause
                            S(i,j)=-origin(i,j)
                            S(i,j+1)=-origin(i,j+1)
                        end if
                    end do
                end do

                do j=1,N
                    do i=1,N
                        if(bond_colume(i,j)==2)then
                            !write(*,*)"bond_colume(",i,j,")"
                            !pause
                            S(i,j)=-origin(i,j)
                            S(i+1,j)=-origin(i+1,j)
                        end if
                    end do
                end do


                do i=1,N
                    if( S(N+1,i)/=origin(1,i))then
                        S(1,i)=S(N+1,i)
                    else
                        S(N+1,i)= S(1,i)
                    end if
                    if(S(i,1)/=origin(i,1))then
                        S(i,N+1)=S(i,1)
                    else
                        S(i,1)=S(i,N+1)
                    end if
                end do

!        S store the configurations



                !date calculate
                !************************************************************
                M_every=0.
                do i=1,N
                    do j=1,N
                        M_every=S(i,j)+S(i,j+1)+M_every
                    end do
                end do
                M_(k)=M_every

                E_every=0.
                do i=1,N
                    do j=1,N
                        E_every=S(i,j)*(S(i+1,j)+S(i,j+1))+E_every
                    end do
                end do
                H(k)=-E_every-Hz*M_every
            end do

            !choose m points
            ! this is how to pick m configurations
            ! insert "write to file" around here
            do i=0,m-1
                d=i*n_0+Nequi
                M_(i+1)=M_(d)
                H(i+1)=H(d)

            end do


            !计算能量平均值  U()

            do i=1,m
                U(big_L,circle)=U(big_L,circle)+H(i)
            end do
            U(big_L,circle)=U(big_L,circle)/real(m)



            !*******************************************************************************************

            !reorder
            E_max=MaxVal(H(1:m))
            E_min=MinVal(H(1:m))
            !write(*,*)  E_max,E_min

            E_points=(E_max-E_min)/0.5_8+1
            write(*,*)  E_points
            ! pause

            allocate(P_(E_points))
            allocate(E(E_points))
            P_=0.

            do j=1,E_points
                E(j)=E_min+(j-1)*0.5_8
            end do

            do i=1,m
                do j=1,E_points
                    if(H(i)==E(j))then
                        P_(j)=P_(j)+1.
                        exit
                    else
                    end if
                end do
            end do
            P_=P_/real(m)



            !record in a 3D array which the first colume is temperature,the second colume is energy,the third colume is rate.
            TER(circle,1:E_points,1)=T0
            TER(circle,1:E_points,2)=E(:)
            TER(circle,1:E_points,3)=P_(:)


            do i=1,E_points
                write(10,"(F10.6,1x,2(F14.7))") TER(circle,i,1), TER(circle,i,2), TER(circle,i,3)
            end do
            close(10)

            deallocate(P_)
            deallocate(E)

        end do


        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call get_density(num_T,m,TER,N)



    end do







    deallocate(S)
    deallocate(origin)
    deallocate(M_)
    deallocate(H)
    deallocate(bond_line)
    deallocate(bond_colume)
    deallocate(change_line)
    deallocate(change_colume)
    deallocate(TER)

    return
    end subroutine




