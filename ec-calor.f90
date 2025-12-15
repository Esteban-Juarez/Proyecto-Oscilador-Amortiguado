program ec_calor
!
! Se resolvera la ecuacion de calor (1D) con condiciones de frontera de Dirichlet homogeneas y condiciones 
! iniciales armonicas por  medio del metodo de diferencias finitas (progresiva y regresiva)
! Los parametros de la constante de calor corresponde a ....
!
! Autor: Juarez, E., E.         Materia: Fisica Computacional           Año: 2025
! 
        implicit none

        integer :: i, j, n, m
        real    :: xmin, xmax, hx, ht, lambda, alfa
        real, allocatable, dimension(:) :: z, w, x
        real    :: T0, DT, l

        open(100,file="datos_entrada.inp")
                read(100,*)
                read(100,*)n, m
                read(100,*)
                read(100,*)xmin, xmax, ht, alfa
                read(100,*)
                read(100,*)T0, DT
        close(100)

        l = xmax - xmin         ! Longitud del intervalo en x
        hx = l/real(n)          ! Longitud de los subintervalos en x

        allocate(x(0:n), z(0:n), w(0:n))        ! Alojamiento0 de arreglos

        do i = 0, n             ! Alojamiento de los arreglos
                x(i) = xmin + real(i)*hx
        end do

        lambda = ht * alfa / hx**2              ! Variable auxiliar para aritmetica

        ! Distribucion de valores iniciales
        ! Condiciones de frontera
        w(0) = T0
        w(n) = T0

        do i = 1, n-1
                w(i) = f(x(i), T0, DT, l)
        end do

        do j = 1, m             ! Ciclo sobre pasos en el tiempo
                ! Condiciones de frontera guardadas en el vector auxiliar
                z(0) = T0
                z(n) = T0

                do i = 1, n-1   ! Solucion al paso siguiente
                        z(i) = (1.0-2.0*lambda)*w(i) + lambda*(w(i-1)+w(i+1))   ! Formula progresiva
                end do

                !open(300,file='resultados.dat')
                do i = 0, n     ! Se vacia la solucion en w
                        w(i) = z(i)
                       ! write(300,*)x(i), w(i)
                        if(mod(j,1000)==0)then  ! Se guarda w cada 1000 pasos
                                write(10000+j,*)x(i), w(i)
                        end if
                end do
                !close(300)
        end do

        contains

                function f(xx, TT0, DDT, LL) result (zz)        ! Define la condicion inicial
                        implicit none
                        real :: xx, zz, TT0, DDT, LL
                        real, parameter :: pi = 3.1415926539
                        ! Condicion inical 
                        zz = TT0 + DDT * sin(pi*xx/LL)
                end function

end program ec_calor
