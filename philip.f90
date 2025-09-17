program philip

! Propositos: 
!       1. Calcular la velocidad y la aceleracion de un oscilador amortiguado mediante la regla central
!       para aproximar derivadas.
!       2. Utilizar los resultados previos para computar las ecuaciones de MUA del oscilador amortiguado.
!       3. Obtener la aceleracion del oscilador amortiguado mediante la segunda ley de Newton.
!       4. Computar la energía disipada del sistema para distintos intervalos por medio de cuadraturas, 
!       en particular, la cuadratura de Simpson.
!
! Nombre del proyecto                      
!       Operaciones basicas de metodos numericos para resolver la dinamica del oscilador amortiguado.
!       
! Alumno
!       Juarez Espinoza Esteban
!
! Materia                                 Semestre
!       Fisica Computacional              2026-1

        implicit none

! Diccionario de variables
        integer                                 :: i, j, n                                                 
        real                                    :: bfricc, masa, kresor         ! Determinan el tipo de amortiguamiento.                                             
        real                                    :: h, l, factor, t_inicial, t_final, cuadratura                                                 
        real, allocatable, dimension(:)         :: t, x, v, a                   
        real, allocatable, dimension(:)         :: x_unif, v_unif, a_cinec
        real, allocatable, dimension(:)         :: s, f, velocidad                    

! Archivo de entrada con los datos para ocupar los metodos de la regla central y cuadratura.
        open(100, file="input_data.inp")
                read(100, *)
                read(100, *)n, bfricc, masa, kresor, h                          ! El numero de pasos sera el mismo para ambos metodos.
                read(100, *)
                read(100, *)t_inicial, t_final             
        close(100)

! Alojamiento de los arreglos vectoriales para la regla central.
        allocate(t(1:n), x(1:n), v(2:n-1), a(3:n-2))
        allocate(x_unif(3:n-2), v_unif(3:n-2), a_cinec(3:n-2))

! Inicializacion del ciclo para asignar la funcion de posicion
        do i = 1, n
                t(i) = i*h                              
                x(i) = pos(bfricc, masa, kresor, t(i))                          ! Invocacion de la funcion posicion.              
        end do

! Etiqueta para señalar el tipo de amortiguamiento.
        if (bfricc ** 2 > 4.0 * kresor * masa)then
                write(*,*) "Sobre-amortiguamiento."
        else if (abs(bfricc - sqrt(8.0)) < 1.0e-4)then
                write(*,*) "Amortiguamiento critico"
        else
                write(*,*) "Sub-amortiguamiento."
        end if

        open(200, file="Resultados.dat")
! Inicizalicion del ciclo para calcular las derivadas aproximadas por la regla central de 3 puntos.
        do i = 2, n-1
                v(i) = (x(i+1) - x(i-1)) / (2 * h)                              ! Primera derivada.
                if (2 < i .and. i < n-1)then
                        a(i) = (x(i+1) + x(i-1) - 2 * x(i)) / (h ** 2)          ! Segunda derivada.
                end if
                write(200, 20)i, t(i), x(i), v(i), a(i)                         
        end do
        close(200)

        open(300, file="Uniforme.dat")
! Inicializacion del ciclo para resolver las ecuaciones de movimiento uniformemente acelerado.
        do i = 3, n-2
                x_unif(i) = x(i) + v(i) * t(i) + 0.5 * a(i) * (t(i)**2)         ! Posicion.
                v_unif(i) = v(i) + a(i) * t(i)                                  ! Velocidad.
                write(300, 21) i, t(i), x_unif(i), v_unif(i)
        end do
        close(300)

        open(400, file="Cinematica.dat")
! Inicializacion del ciclo para calcular la aceleracion mediente la Segunda Ley de Newton.
        do i = 3, n-2
                a_cinec(i) = -((kresor/masa) * x(i) + (bfricc/masa) * v(i))     ! Aceleracion.
                write(400, 22) i, a_cinec(i)
        end do
        close(400)

        deallocate(t, x, v, a)
        deallocate(x_unif, v_unif, a_cinec)

! Alojamiento de los arreglos vectoriales para la cuadratura de Simpson.
        allocate(s(1:n), f(1:n), velocidad(2:n-1))

        if(mod(n,2)==0)then                                                     ! n debe ser impar para usar la cuadratura.
                write(*,*) "El numero de pasos es par"
                write(*,*) "El programa se detendra"
                stop
        end if

        l = (t_final - t_inicial)/real(n-1)                                     ! Tamaño de los subintervalos.                 

! Inicializacion del ciclo para asignar los pasos para avanzar y la funcion posicion.
        do j = 1, n
                s(j) = t_inicial + real(j-1) * l
                f(j) = pos(bfricc, masa, kresor, s(j)) 
        end do

! Inicializacion del ciclo para asignar el integrando.
        do j = 2, n-1
                velocidad(j) = ((f(j+1) - f(j-1)) / (2.0 * h)) ** 2
        end do

! Se suman por separado los extremos de la cuadratura.
        cuadratura = velocidad(2) + velocidad(n-1)

! Inicializacion del ciclo para sumar los terminos restantes.
        do j = 3, n-2
                if(mod(j,2) == 0)then
                        factor = 4.0                                            ! Integrando par.
                else
                        factor = 2.0                                            ! Integrando impar.
                end if
                cuadratura =  cuadratura + (factor * velocidad(j))
        end do

        write(*,*)"La cuadratura de Simpson nos da La energía disipada del oscilador igual a:"
        write(*, 23) -bfricc * l * cuadratura / 3.0

        deallocate(s, f, velocidad)

        20 format(1I4, 4F15.8)
        21 format(1I4, 3F15.8)
        22 format(1I4, 1F15.8)
        23 format(1F15.8)

        contains

                function pos(b, m, k, t)
                        real :: b, m, k, t
                        real :: A, E, C, D
                        real :: decrece, raiz, auxiliar, alfa, cuadrado
                        real :: pos

                        decrece = -b / (2.0 * m)
                        raiz = sqrt((b ** 2) - 4.0 * k * m) / (2.0 * m)
                        auxiliar = sqrt(-(b ** 2) + 4.0 * k * m) / (2.0 * m)
                        alfa = 1 / sqrt(1 - ((4.0 * k * m) / (b ** 2)))
                        cuadrado = sqrt(8.0) 

                        ! Sobreamortiguado
                        if (b ** 2 > 4.0 * k * m)then
                                A = 0.5 * (1.0 + (1.0 / alfa)) 
                                E = 0.5 * (1.0 - (1.0 / alfa)) 
                                pos = exp(decrece * t) * (A * exp(raiz * t) + E * exp(-raiz * t))
                        
                        ! Amortiguamiento critico
                        else if (abs(b - cuadrado) < 1.0e-4)then                ! Dado que el valor crítico se obtiene cuando
                                                                                ! la constante de amortiguamiento b es un número
                                                                                ! racional (sqrt(8)) vamos a acercarnos al valor 
                                                                                ! con una diferencia menor que una tolerancia.
                                A = 1.0
                                E = -decrece
                                pos = (A + (B * t)) * exp(decrece * t)

                        ! Sub-amortiguado
                        else
                                C = 1.0 
                                D = b / sqrt(4.0 * k * m - (b ** 2))
                                pos = exp(decrece * t) * (C * cos(auxiliar * t) + D * sin(auxiliar * t))
                        end if
                
                end function pos

end program philip
