program main
    implicit none
    integer, parameter::L=60
    integer, parameter::N=L*L
    integer::MAT(L,L)                   ! matriz de estado
    integer::NNX(L,L,4)                 ! Nearest neighbors, coordenadas x
    integer::NNY(L,L,4)                 ! Nearest neighbors, coordenadas y
    integer, parameter::term=5000       ! pasos montecarlo termalización
    integer, parameter::tiempo=100000   ! pasos montecarlo proceso
    real::temp              ! temperatura
    integer::t, p 
    real*8::rand1, rand2
    integer::a,b, i, j, k, cosa
    real::m, m0, m1, m2, m4, m8, mm, m2m, m4m, m8m
    real::ra, c
    real::h(-4:4)
    real::varmag, varmag2, varmag4
    real::e0, e, e1, e2, e4, em, e2m, e4m, vare, vare2, errore, ce, errorce, suscep, binder, tau
    real::errorm, errorbinder, errorsuscep
    real::cm
    real::start, sto

    call random_seed()
    ! inicializo el estado y determino nearest neighbors
    do i=1,L
        do j=1,L
            call random_number(rand1)
            if (rand1<0.5) then 
                MAT(i,j)=1
            else 
                MAT(i,j)=-1
            endif

            ! determino los NN con pbc
            if (i>1) then    !caso N
                NNY(i,j,1)=j
                NNX(i,j,1)=i-1
            else 
                NNY(i,j,1)=j
                NNX(i,j,1)=L
            endif
            if (i<L) then    !caso S
                NNY(i,j,2)=j
                NNX(i,j,2)=i+1
            else 
                NNY(i,j,2)=j
                NNX(i,j,2)=1
            endif
            if (j>1) then    !caso W
                NNY(i,j,3)=j-1
                NNX(i,j,3)=i
            else 
                NNY(i,j,3)=L
                NNX(i,j,3)=i
            endif
            if (j<L) then     !caso E
                NNY(i,j,4)=j+1
                NNX(i,j,4)=i
            else 
                NNY(i,j,4)=1
                NNX(i,j,4)=i
            endif
        enddo
    enddo

    open(0, status='unknown', file='last60_2.txt')

    temp=3.5
    do t=1,20
        if (t.eq.1) then
            call cpu_time(start)
        endif
        do j=-4,4,2
            h(j)=min(1.0, exp(-2*j/temp))
        enddo

        ! termalización
        do i=1,term*N
            cosa=0
            ! genero enteros aleatorios (lo he hecho matricial porque quería probar unas cosas)
            call random_number(rand1)
            call random_number(rand2)
            a=1+floor((L)*rand1)        
            b=1+floor((L)*rand2)

            do j=1,4
                cosa=cosa+MAT(NNX(a,b,j),NNY(a,b,j))
            enddo
            cosa=cosa*MAT(a,b)
            call random_number(rand1)
            if (rand1<h(cosa)) then
                MAT(a,b)=-MAT(a,b)
            endif
        enddo

        m1=real(abs((sum(MAT)))/real(N))
        ! inicializo las variables en cada T
        m=0
        m2=0
        m4=0
        m8=0
        c=0
        e=0
        e2=0
        e4=0
        binder=0
        suscep=0
        ce=0

        ! proceso
        do i=1,tiempo
            ! se hace un paso MC
            do k=1,N
                cosa=0
                call random_number(rand1)
                call random_number(rand2)
                a=1+floor((L)*rand1)
                b=1+floor((L)*rand2)

                do j=1,4
                    cosa=cosa+MAT(NNX(a,b,j),NNY(a,b,j))
                enddo
                cosa=cosa*MAT(a,b)
                call random_number(rand1)
                if (rand1<h(cosa)) then
                    MAT(a,b)=-MAT(a,b)
                endif
            enddo

            ! mido y sumo cada 10 pasos MC
            if (mod(i,10)==0) then
                e0=0.0
                ! mido la energía
                do j=1,L
                    do k=1,L
                        e1=0.0
                        ! calculo el campo local
                        do p=1,4
                            e1=e1+MAT(NNX(j,k,p),NNY(j,k,p))
                        enddo
                        e0=e0-e1*MAT(j,k)
                    enddo   
                enddo

                e0=0.5*e0
                e=e+e0
                e2=e2+e0**2.0
                e4=e4+e0**4.0

                m0=real(abs(sum(MAT)))/real(N)
                m=m+m0
                m2=m2+m0**2.0
                m4=m4+m0**4.0
                m8=m8+m0**8.0
                c=c+m0*m1
                m1=m0
                
            endif
            
            ! cada 1000 pasos MC saco una foto
            if (mod(i,1000)==0) then

                em=10.0*e/real(i)             ! <e>
                e2m=10.0*e2/real(i)           ! <e2>
                e4m=10.0*e4/real(i)           ! <e4>

                mm=10*m/real(i)             ! <m>
                m2m=10.0*m2/real(i)           ! <m2>
                m4m=10.0*m4/real(i) 
                m8m=10.0*m8/real(i) 

                varmag=m2m-mm*mm            ! <m2>-<m>2
                vare=e2m-em*em              ! <e2>-<e>2
                varmag2=m4m-m2m*m2m         ! <m4>-<m2>2
                vare2=e4m-e2m*e2m           ! <e4>-<e2>2
                varmag4=m8m-m4m*m4m

                binder=1-real(i)*m4/(30.0*m2*m2)  ! coef binder
                suscep=real(N)*varmag/temp      ! susceptibilidad
                ce=vare/(real(N)*temp**2)       ! calor específico

                cm=((10.0*c/real(i))-mm*mm)/m2m   ! función de autocorrelación

                if (cm.ne.1) then
                    tau=cm/(1-cm)
                endif
                ra=sqrt((2*tau+1)*10/real(i))

                ! calculo los errores
                errorsuscep=ra*(real(N)/temp)*sqrt(varmag2+4.0*m2m*varmag)
                errorce=(ra/(real(N)*temp**2))*sqrt(vare2+4.0*e2m*vare)
                errorbinder=(ra/(sqrt(3.0)*m2m))*sqrt(varmag4+2.0*m4m*varmag2/m2m)
                errorm=varmag*ra
            
                errore=vare*ra/real(N)
                em=em/real(N)

                write(0,*), real(i)/10, temp, mm, errorm, em, errore, binder, errorbinder, suscep, errorsuscep, ce, errorce, tau, cm

            endif
        enddo

        temp=temp-0.1
        print*, 'medida', t, 'terminada'

        ! estimo tiempo de ejecución total en la primera medida
        if (t.eq.1) then
            call cpu_time(sto)
            sto=sto-start
            print*, 'quedan ', sto*19/60, 'minutos'
        endif

    enddo
    close(0)
end program