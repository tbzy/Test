!------------------------------------------------------------------------
! fourier.F90  (Last updated: 2015-10-25)
!------------------------------------------------------------------------

#include "abbrev.inc"


!////////////////////////////////////////////////////////////////////////
module fourier_mod

    using( lib_ppfilt ) get_maskfunc, get_kdampfac, alpha, gamma, &
                        lib_ppfilt_init, use_maskfunc_k
    ___
    save

    Char(20) :: func_type = "p0"

    Dble :: rc = 0.64d0   !! cutoff of p(r)
    Dble :: h  = 0.35d0   !! grid spacing in real space

    Int :: bess_order = -1   !! order of spherical Bessel function

    Int :: Ngrid = 1000  !! number of integration points

    Int  :: npow  = 0       !! p(r) ~ r^npow  (r -> 0)
    Dble :: sigma = 1.0d0   !! width param in p(r)

    !! Params for p0(r).
    Dble :: ca( 5 ) = 0.0d0
    Dble :: cx( 5 ) = 0.0d0
    Dble :: cw( 5 ) = 0.1d0

    Bool :: Q_cut_pr = .true.   !! if true, enforce p(r) = 0 for r > rc

    !! Extended integral domain.
    Int :: mult = 3
    
    namelist / fourier_inp / &
            Ngrid, rc, h, sigma, &
            ca, cx, cw, &
            npow, func_type, &
            use_maskfunc_k, bess_order, Q_cut_pr

contains

!------------------------------------------------------------------------
sub test_fourier

    Int   ir, ik

    Dble  r, k, dr, dk, rmask, kmax, D
    Dble  pr, pr_org
    Dble  inte_org, inte_orgrc, inte_step, inte_damp, inte_mask
    Dble  inte_nofilt
    Dble  err_step, err_mask, tmp

    Dble, _alloc_, _dim_(:) :: pr_step, pr_mask, pk_org, pk_mask

    Char(50), _param_ :: fmt = "(10es20.10)"
    !.........................................................

    call lib_ppfilt_init

!! Read params.

    put "file(inp): input"

    open( 10, file="input", status="old" )
    read( 10, nml= fourier_inp )
    close( 10 )

    printf( nml= fourier_inp )

!! Init params.

    rmask = rc * gamma   !! domain of the mask function ("rcut" of Schmid et al)
    kmax  = _pi_ / h

    !! Grid spacing for quadrature.
    dr = rmask / dble( Ngrid )
    dk = kmax  / dble( Ngrid )

    put
    put "rmask = ", rmask
    put "kmax  = ", kmax
    put "kdamp = ", kmax / alpha, "(= kmax/alpha)"  !! "qcut" of Schmid et al

    allocate( pk_org ( 0 : mult * Ngrid ), &
              pk_mask( 0 : mult * Ngrid ), &
              pr_step( 0 : mult * Ngrid ), &
              pr_mask( 0 : mult * Ngrid ) )

!!
!! Forward transform.
!!

    put
    put "file(out): pr_org.dat = original p(r)"
    open( 500, file="pr_org.dat" )

    put
    put "k-space func:"
    put "file(out): pk_org.dat   = Intr[0,rmask] r^2 j(kr) p(r)"
    put "file(out): pk_orgrc.dat = Intr[0,rc]    r^2 j(kr) p(r)"
    put "file(out): pk_mask.dat  = Intr[0,rc]    r^2 j(kr) p(r) / m(r/rmask)"

    open( 10, file="pk_org.dat" )
    open( 20, file="pk_orgrc.dat" )
    open( 30, file="pk_mask.dat" )

    ! Loop over k-grid.
    !-------------------------
    do ik = 1, mult * Ngrid    !! k=0 not calculated

        k = dble( ik ) * dk

        inte_org   = _0_
        inte_orgrc = _0_
        inte_mask  = _0_

        ! Integration over r in [0,rmax].
        !-------------------------
        do ir = 1, Ngrid     !! r = 0 is excluded because of r^2

            r = dble( ir ) * dr

            call calc_radial_func ( r, pr )
              !... pr is supposed to zero for r > rc.

            tmp = dr * r**2 * get_bess( k * r ) * pr

            _( inte_org ) + tmp

            !! integrate over [0,rc].
            if ( r < rc ) then

                _( inte_orgrc ) + tmp
                _( inte_mask  ) + tmp / get_maskfunc( r / rmask )
            _fi

            if ( ik == 1 ) then
                write( 500, fmt ) r, gpcut( pr ), gplog( pr )
            _fi

        enddo
        !-------------------------

        pk_org ( ik ) = inte_org
        pk_mask( ik ) = inte_mask

        write( 10, fmt ) k, gpcut( inte_org ),   gplog( inte_org   )
        write( 20, fmt ) k, gpcut( inte_orgrc ), gplog( inte_orgrc )
        write( 30, fmt ) k, gpcut( inte_mask ),  gplog( inte_mask  )

    enddo
    !-------------------------

    call fclose([ 10, 20, 30, 500 ])

!!
!! Backward transform.
!!

    put
    put "r-space func:"
    put "file(out): pr_nofilt.dat ~ Intk k^2 j(kr) p(k)"
    put "file(out): pr_step.dat   ~ Intk k^2 j(kr) p(k) h(k)"
    put "file(out): pr_damp.dat   ~ Intk k^2 j(kr) p(k) D(k)"
    put "file(out): pr_mask.dat   ~ Intk k^2 j(kr) p_m(k) D(k)" 

    open( 10, file="pr_nofilt.dat" )
    open( 20, file="pr_step.dat" )
    open( 30, file="pr_damp.dat" )
    open( 40, file="pr_mask.dat" )

    put
    put "file(out): pk_damp.dat  = pk_org(k) * D(k)"
    put "file(out): pk_maskD.dat = pk_mask(k) * D(k)"

    open( 500, file="pk_damp.dat" )
    open( 510, file="pk_maskD.dat" )

    err_step = _0_
    err_mask = _0_

    ! Loop over r in [0, rmask].
    !-------------------------
    do ir = 1, Ngrid

        r = dble( ir ) * dr

        call calc_radial_func ( r, pr_org )

        inte_nofilt = _0_
        inte_step   = _0_
        inte_damp   = _0_
        inte_mask   = _0_

        ! Integration over k in [0, mult * kmax].
        !-------------------------
        do ik = 1, mult * Ngrid    !! k = 0 is excluded because of k^2

            k = dble( ik ) * dk

            tmp = (2.0d0 / _pi_) * dk * k**2 * get_bess( k * r )

            D = get_kdampfac( k, kmax/alpha, kmax )

            _( inte_nofilt ) + tmp * pk_org( ik )

            if ( k < kmax ) _( inte_step ) + tmp * pk_org( ik )

            _( inte_damp ) + tmp * pk_org( ik ) * D

            _( inte_mask ) + tmp * pk_mask( ik ) * D
            !...Note: To perform exact cutoff, set use_maskfunc_k = T.

            if ( ir == 1 ) then

                tmp = pk_org( ik ) * D
                write( 500, fmt ) k, gpcut( tmp ), gplog( tmp )

                tmp = pk_mask( ik ) * D
                write( 510, fmt ) k, gpcut( tmp ), gplog( tmp )
            _fi

        enddo
        !-------------------------

        _( inte_mask ) * get_maskfunc( r / rmask )

        write( 10, fmt ) r, gpcut( inte_nofilt ), gplog( inte_nofilt ), &
                            gpcut( inte_nofilt - pr_org )
        write( 20, fmt ) r, gpcut( inte_step  ), gplog( inte_step )
        write( 30, fmt ) r, gpcut( inte_damp  ), gplog( inte_damp )
        write( 40, fmt ) r, gpcut( inte_mask  ), gplog( inte_mask )

        err_step = max( err_step, abs( inte_step - pr_org ) )
        err_mask = max( err_mask, abs( inte_mask - pr_org ) )

        pr_step ( ir ) = inte_step
        pr_mask ( ir ) = inte_mask

    enddo
    !-------------------------

    put
    put "err_step = ", err_step
    put "err_mask = ", err_mask

    call fclose([ 10, 20, 30, 40, 500, 510 ])

!!
!! Again, transform to k-space.
!!

    put
    put "k-space func:"
    put "file(out): pr_step_tok.dat"
    put "file(out): pr_mask_tok.dat"

    open( 10, file="pr_step_tok.dat" )
    open( 20, file="pr_mask_tok.dat" )

    ! Loop over k in [0, mult * kmax].
    !-------------------------
    do ik = 1, mult * Ngrid

        k = dble( ik ) * dk

        inte_step = _0_
        inte_mask = _0_

        ! Integration over r in [0, rmask].
        !-------------------------
        do ir = 1, Ngrid     !! r = 0 is not excluded because of r^2

            r = dble( ir ) * dr

            tmp = dr * r**2 * get_bess( k * r )

            _( inte_step ) + tmp * pr_step( ir )
            _( inte_mask ) + tmp * pr_mask( ir )

        enddo
        !-------------------------

        write( 10, fmt ) k, gpcut( inte_step ), gplog( inte_step )
        write( 20, fmt ) k, gpcut( inte_mask ), gplog( inte_mask )

    enddo
    !-------------------------

    call fclose([ 10, 20 ])

endsub

! Test functions.
!------------------------------------------------------------------------
sub calc_radial_func ( r, fun )

    Dble, _in_ :: r
    Dble       :: fun

    Dble  tmp, rloc
    !.........................................................

    !-------------------------
    selectcase ( func_type )

    case ( "p0" )
        fun = p0( r ) - p0( rc )

    case ( "dvloc_O" )

        rloc = 0.25

        if ( r < 1.0d-6 ) then
            fun = -6.0 * sqrt( 2.0 / _pi_ ) &
                    * ( 1.0 / rloc - 1.0 / sigma )
        else
            fun = -6.0 * (   erf( r / (_sqrt2_ * rloc ) ) / r &
                           - erf( r / (_sqrt2_ * sigma) ) / r )
        _fi

        tmp = r / rloc
        _( fun ) + ( -16.5 + 2.4 * tmp**2 ) * exp( -0.5 * tmp**2 )

    case default
        _stop_( "invalid func_type" )

    endselect
    !-------------------------

    if ( npow > 0 ) _( fun ) * r ** npow

    if ( Q_cut_pr .and. r > rc ) fun = _0_

contains

    function p0( x ) result( val )
        Dble :: x, val

        val = 2.0 * exp( -0.5 * ( x / 0.5 )**2 ) &
                + sum( ca(:) * exp( -( (x-cx(:)) / cw(:) )**2 ) )
               !... ca = -0.5, cx = 0.0, cw = 0.1 are good.
    endfunc

endsub

! Spherical Bessel function. The order is specified by bess_order.
!------------------------------------------------------------------------
function get_bess ( x ) result ( ret )

    use lib_bessel, only: spherical_bessel
    Dble :: x, ret

    if ( bess_order == -1 ) then

        if ( x < 1.0d-6 ) then
            ret = 1.0d0
        else
            ret = sin( x ) / x   !! j0( x )
        _fi
    else
        !! Use a library function.
        call spherical_bessel ( bess_order, x, ret )
    _fi

endfunc

!------------------------------------------------------------------------
sub fclose ( units )

    Int, _in_ :: units(:)
    Int  i

    do i = 1, size( units )
        close( units( i ) )
    enddo
endsub

! For gnuplot.
!------------------------------------------------------------------------
function gplog ( x ) result ( ret )
    Dble :: x, ret

    ret = log10( abs(x) + 1.0d-12 )

endfunc

! For gnuplot.
!------------------------------------------------------------------------
function gpcut ( x ) result ( ret )
    Dble :: x, ret

    if ( abs( x ) < 1.0d-12 ) then
        ret = dsign( 1.0d-12, x )
    else
        ret = x
    _fi

endfunc

end module

