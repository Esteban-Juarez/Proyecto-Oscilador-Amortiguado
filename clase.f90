! En el siguiente programa se calcular de manera simple la aceleración del oscilador con base al código visto en clase. 

program aceleracion
        
        implicit none

        integer :: i, n
        real    :: m, b, k, a, x, v, h, t, Edis

        m = 0.5
        k = 4.0
        b = 0.5
        h = 0.1

        ! Las condiciones iniciales serviran para usar el Teorema de Existencia y Unicidad para EDO 
        x = 1.0
        v = 0.0

        n = 101
        Edis = 0.0

        open(200, file='clase_resultados.dat')

        do i = 1, n
                a = -(b/m) * v - (k/m) * x

                x = x + v * h + 0.5 * a * h**2
                v = v + a * h
                t = real(i) * h

                write(200, *) t, x, v, a

                if(i == 1 .or. i == n)then
                        Edis = Edis + v**2
                else
                        Edis = Edis + 2.0 * v**2
                end if
        end do

        Edis = b * h * Edis
        write(*,*)Edis
        close(200)
end program aceleracion
