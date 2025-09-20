program talacha
! En este programa se calcularán los valores exactos de la posición, velocidad y aceleración para el oscilador amortiguado. 
! Con el resultado de la velocidad calcularemos la energía disipada.
implicit none 

integer         :: i, n
real            :: b, k, m, h
real            :: A, E, C, D
real            :: decrece, raiz, auxiliar, alfa, cuadrado
real, allocatable, dimension(:)            :: t, pos, vel, acel

b = 3.0
k = 4.0
m = 0.5
h = 0.1
n = 101

decrece = -b / (2.0 * m)
raiz = sqrt((b ** 2) - 4.0 * k * m) / (2.0 * m)
auxiliar = sqrt(-(b ** 2) + 4.0 * k * m) / (2.0 * m)
alfa = 1 / sqrt(1 - ((4.0 * k * m) / (b ** 2)))
cuadrado = sqrt(8.0)

allocate(t(1:n), pos(1:n), vel(1:n), acel(1:n))

open(200, file='exacto.dat')
do i = 1, n
        t(i) = i * h
        ! Sobreamortiguado
        if (b ** 2 > 4.0 * k * m)then
                A = 0.5 * (1.0 + (1.0 / alfa))
                E = 0.5 * (1.0 - (1.0 / alfa))
                pos(i) = exp(decrece * t) * (A * exp(raiz * t) + E * exp(-raiz * t))
                vel(i) = (decrece * pos(i)) + (raiz * exp(decrece * t)) * (A * exp(raiz * t) - E * exp(-raiz * t))
                acel(i) = (decrece * vel(i)) + ((raiz**2) * pos(i))

        ! Amortiguamiento critico
        else if (abs(b - cuadrado) < 1.0e-4)then                ! Dado que el valor crítico se obtiene cuando
                                                                                ! la constante de amortiguamiento b es un número
                                                                                ! racional (sqrt(8)) vamos a acercarnos al valor
                                                                                ! con una diferencia menor que una tolerancia.
                A = 1.0
                E = -decrece
                pos(i) = (A + (B * t)) * exp(decrece * t)
                vel(i) = (decrece * pos(i)) + (E * exp(decrece * t))
                acel(i) = (decrece * vel(i)) + (decrece * E * exp(decrece * t))

        ! Sub-amortiguado
        else
                C = 1.0
                D = b / sqrt(4.0 * k * m - (b ** 2))
                pos(i) = exp(decrece * t) * (C * cos(auxiliar * t) + D * sin(auxiliar * t))
                vel(i) = (decrece * pos(i)) + (auxiliar * exp(decrece * t)) * (C * exp(auxiliar * t) - D * exp(-auxiliar &
                                & * t))
                acel(i) = decrece * vel(i) + (auxiliar**2) * pos(i)
        end if
        write(*,20)t(i), pos(i), vel(i), acel(i)
end do
close(200)

deallocate(t, pos, vel, acel)

20 format(4F15.8)

end program talacha
