program pendulo_rk4
! 
! Titulo: Metodo de Runge-Kutta de cuarto orden para resolver la dinámica de sistemas caóticos, pendulo amortiguado y forzado,
! bajo diferentes condiciones (Version Fortran).
! Materia: Fisica Computacional          Autor: Esteban Juarez Espinoza          Semestre: 2026
! Objetivos:
!       1. Implementar el metodo de Runge Kutta de cuarto orden para resolver el sistema.
!       2. Identificar parametros para los que el sistema muestra comportamiento periodico y caotico.
!       3. Realizar graficas, (angulo vs tiempo, espacio fase, diferencia de trayectorias) para analizar visualmente
!       la dinamica del sistema.
!
        implicit none

        integer                                 :: i, j              ! Contadores
        integer                                 :: n                 ! Pasos
        real                                    :: h, t, t0          ! Separacion, tiempo, tiempo inicial
        real                                    :: rho, tau, nu      ! rho = coef fricc, tau = amp externo, nu = frec motor
        real, parameter                         :: pi = 3.1416 
        real                                    :: x1, x2, x3        ! Cambios de variable
        real                                    :: x01, x02, x03     ! Condiciones iniciales (ang y vel ang)
        real                                    :: cambio            ! Diferencia de trayectorias
        real, allocatable, dimension(:,:)       :: x                 ! Soluciones a la edo
        real, allocatable, dimension(:,:)       :: k                 ! Pendientes de Runge - Kutta
        
        ! Archivo de entrada
        open(100, file = "pendulo.inp")
                read(100,*) 
                read(100,*)n, t0, h 
                read(100,*) 
                read(100,*)rho, tau, nu 
                read(100,*) 
                read(100,*)x01, x02 
        close(100)

        allocate(x(1:3,0:n), k(1:3, 1:4))                       ! Arreglos matriciales (renglon, columna)

        ! Asignacion de las condiciones iniciales antes de iniciar el ciclo
        x(1, 0) = x01                                           
        x(2, 0) = x02
        x(3, 0) = t0
        t = t0

        ! Salida de resultados
        open(200, file = "posicion.dat")                        ! Angulos
        open(300, file = "velocidad.dat")                       ! Velocida angular
        open(400, file = "espaciofase.dat")                     ! Velocidad angular vs angulo
        open(500, file = "Trayectorias.dat")                    ! Diferencia de trayectorias

        write(200, *)t, x(1, 0)
        write(300, *)t, x(2, 0)
        write(400, *)x(1, 0), x(2, 0)

        ! Ciclo sobre el numero de pasos
        do i = 1, n
                x1 = x(1, i-1)
                x2 = x(2, i-1)
                x3 = x(3, i-1)

                ! Primeras pendientes
                k(1, 1) = func1(x2)
                k(2, 1) = func2(x1, x2, x3, t, rho, tau, nu)
                k(3, 1) = func3(nu)

                x1 = x(1, i-1) + h * k(1, 1) / 2.0
                x2 = x(2, i-1) + h * k(2, 1) / 2.0
                x3 = x(3, i-1) + h * k(3, 1) / 2.0

                ! Segundas pendientes
                k(1, 2) = func1(x2)
                k(2, 2) = func2(x1, x2, x3, t + h / 2.0, rho, tau, nu)
                k(3, 2) = func3(nu)

                x1 = x(1, i-1) + h * k(1, 2) / 2.0
                x2 = x(2, i-1) + h * k(2, 2) / 2.0
                x3 = x(3, i-1) + h * k(3, 2) / 2.0

                ! Terceras pendientes
                k(1, 3) = func1(x2)
                k(2, 3) = func2(x1, x2, x3, t + h / 2.0, rho, tau, nu)
                k(3, 3) = func3(nu)

                x1 = x(1, i - 1) + h * k(1, 3)
                x2 = x(2, i - 1) + h * k(2, 3)
                x3 = x(3, i - 1) + h * k(3, 3)

                ! Cuartas pendientes
                k(1, 4) = func1(x2)
                k(2, 4) = func2(x1, x2, x3, t + h, rho, tau, nu)
                k(3, 4) = func3(nu)
                  
                ! Regla de Runge - Kutta
                advance: do j = 1, 3                                   
                        x(j, i) = x(j, i - 1) + h * (k(j, 1) + 2.0 * k(j, 2) + 2.0 * k(j, 3) + k(j, 4)) / 6.0
                end do advance

                cambio = abs(x(1, i) - x(2, i))

                t = t + h

                write(200, 20)t, x(1, i)
                write(300, 20)t, x(2, i)
                write(400, 20)x(1, i), x(2, i)
                write(500, 20)t, cambio

        end do 

        close(200)
        close(300)
        close(400)
        close(500)

        20 format(2F15.8)
        contains

                function func1(x)       result(y)
                        implicit none
                        real :: x, y
                        y = x
                end function func1

                function func2(x1, x2, x3, t, rho, tau, nu) result(y)
                        implicit none
                        real :: x1, x2, x3, t, rho, tau, nu, y
                        x3 = 2 * pi * nu * t
                        y = - rho * x2 - sin(x1) + tau * sin(x3)
                end function func2

                function func3(x)       result(y)
                        implicit none
                        real :: x, y
                        y = 2*pi*x
                end function func3

end program pendulo_rk4
