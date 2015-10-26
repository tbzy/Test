!------------------------------------------------------------------------
! lib_bessel.F90  (Last updated: 2014-11-09)
!
! Library for Bessel functions.
!
! Todo:
! * Find a robust spherical Bessel routine valid for both x -> 0 and infty.
!------------------------------------------------------------------------

#include "abbrev.inc"


!////////////////////////////////////////////////////////////////////////
module lib_bessel
    ___
    save

contains


!========================================================================
sub spherical_bessel ( n, x, jn )

!! mask routine.
    ___
    Int, _in_   :: n
    Dble, _in_  :: x
    Dble        :: jn   !! [out]
    !.........................................................

    call spherical_bessel_local ( n, x, jn )

    ! call spherical_bessel_jin ( n, x, jn )   !! not recommended

endsub


!========================================================================
sub spherical_bessel_local ( n, x, jn, jn_eps_ )

!! Local code.
    ___
    Int, _in_   :: n
    Dble, _in_  :: x
    Dble        :: jn  !! [out]
    Dble, _opt_ :: jn_eps_

    Int   m
    Dble  s, c, y, j0, j1, jn_1, jn_2, jn_eps
    !.........................................................

    _optval_( jn_eps, jn_eps_, 1.0d-6 )

    s = sin( x )
    c = cos( x )

!! Use approximate forms for x -> 0 when appropriate.

    if ( n == 0 ) then

        !! j0(x) = sin(x)/x = 1 - x^2/6 + O(x^4)
        jn = 1.0d0 - x**2 / 6.0d0   
    else
        !! jn(x) ~ x^n / (2n+1)!!
        jn = 1.0d0
        do m = 1, n
            _( jn ) * x / dble( 2 * m + 1 )
        enddo
    _fi

    if ( ( n == 0 .and. abs( jn - 1.0d0 ) < jn_eps ) .or. &
         ( n >= 1 .and. abs( jn )         < jn_eps ) ) then

        if ( present( jn_eps_ ) ) put "use approximate form for x = ", x
        return
    _fi

!! Otherwise, use regular forms.

    y = 1.0d0 / x

    j0 = s * y
    j1 = s * y**2 - c * y

    !-------------------------
    selectcase ( n )

    case ( 0 )
        jn = j0

    case ( 1 )
        jn = j1

    case default   !! Use recursion formula.

        jn_2 = j0
        jn_1 = j1

        do m = 2, n
            jn = dble( 2 * m - 1 ) * y * jn_1 - jn_2

            jn_2 = jn_1  !! prepare for next step
            jn_1 = jn
        enddo

    endselect
    !-------------------------

endsub


endmodule


!////////////////////////////////////////////////////////////////////////
! Unit tests.

#ifdef lib_bessel_test


!========================================================================
program main

    call test_spherical_bessel_continuity

    ! call test_spherical_bessel

    ! call plot_spherical_bessel

contains


!========================================================================
sub test_spherical_bessel_continuity ()

    use lib_bessel
    ___
    Int   n, nmax, m

    Dble  dfac, xcut, xleft, xright, shift, inc_rate
    Dble  jn_eps, jn_left, jn_right
    !.........................................................

    put "Input: jn_eps"
    read(*,*) jn_eps

    nmax  = 5
    shift = 1.0d-2

    put "jn_eps = ", jn_eps
    put "shift  = ", shift
    put "nmax   = ", nmax

    !-------------------------
    do n = 0, nmax

        put
        put "n = ", n

        if ( n == 0 ) then

            xcut = sqrt( 6.0d0 * jn_eps )
        else

            !! Double factorial (2n+1)!!
            dfac = 1.0d0
            do m = 1, n
                _( dfac ) * dble( 2 * m + 1 )
            enddo

            xcut = ( jn_eps * dfac )**( 1.0d0 / dble( n ) )
        _fi

        xleft  = xcut * ( 1.0d0 - shift )
        xright = xcut * ( 1.0d0 + shift )

        call spherical_bessel_local &
                ( n, xleft, jn_left, jn_eps_= jn_eps )

        call spherical_bessel_local &
                ( n, xright, jn_right, jn_eps_= jn_eps )

        inc_rate = ( jn_right - jn_left ) / jn_left

        printf(5000) "xcut   = ", xcut
        printf(5000) "xleft  = ", xleft,  "jn_left  = ", jn_left
        printf(5000) "xright = ", xright, "jn_right = ", jn_right
        printf(5100) "inc_rate = ", inc_rate

        if ( n >= 1 ) then
            printf(5100) "inc_rate_ref = ", 2.0d0 * n * shift
        _fi

    enddo
    !-------------------------

5000 format( a10, es20.12, a20, es20.12 )
5100 format( 30x,          a20, es20.12 )

endsub


!========================================================================
sub test_spherical_bessel ()

    use lib_bessel
    ___
    Int   n
    Dble  x, jn_loc, jn_gen
    !.........................................................

    do
        put "Input : n, x"
        read(*,*) n, x

        call spherical_bessel_local( n, x, jn_loc )
        call spherical_bessel      ( n, x, jn_gen )

        put "jn_loc = ", jn_loc
        put "jn_gen = ", jn_gen
    enddo

endsub


!========================================================================
sub plot_spherical_bessel ()

    use lib_bessel
    ___
    Int   ndiv, n, i
    Dble  xmax, x, jn_loc, jn_gen
    !.........................................................

    put "Input : n, xmax, ndiv"
    read(*,*) n, xmax, ndiv

    open( 100, file="jn_loc.dat" )
    open( 200, file="jn_gen.dat" )

    open( 1000, file="diff.dat" )

    do i = 0, ndiv

        x = dble( i ) * xmax / dble( ndiv )

        call spherical_bessel_local ( n, x, jn_loc )
        call spherical_bessel       ( n, x, jn_gen )

        write( 100, "(2es30.20)" ) x, jn_loc
        write( 200, "(2es30.20)" ) x, jn_gen

        write( 1000, "(2es30.20)" ) x, jn_gen - jn_loc

    enddo

endsub


endprogram


#endif /* lib_bessel_test */

