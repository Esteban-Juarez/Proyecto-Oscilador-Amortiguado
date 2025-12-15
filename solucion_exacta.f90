program solucion_exacta
! Programa para calcular la solucion exacta de la ecuacion de calor
        implicit none

        integer :: i, n
        real    :: alpha2, t, lmax, h
        real, parameter :: pi = 3.141592654
        real, allocatable, dimension(:) :: u

        open(100,file='parametros_solexa.inp')
                read(100,*)
                read(100,*)alpha2, t, n
                read(100,*)
                read(100,*)lmax
        close(100)

        h = lmax / n

        allocate(u(1:n))
        
        open(200, file='exacto.dat')
                do i = 0, n
                
                        u(i) = exp(-pi**2 * alpha2 * t) * sin(pi * h * i)
                        write(200,20)real(i)/10, u(i)
                end do 
        close(200)

        deallocate(u)

        20 format(1F15.2,1F15.8)

end program solucion_exacta

