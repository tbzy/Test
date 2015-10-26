!------------------------------------------------------------------------
! lib_ppfilt.F90  (Last updated: 2015-10-25)
!
! Fourier filtering of pseudopotential (PP) based on a mask function.
!
! Ref: Tafipolsky, Schmid, J.Chem.Phys. 124, 174102 (2006)
!      Wang, Phys.Rev.E 64, 201107(R) (2001)
!
! Todo:
! * Receive calc_pfunc() and r as arguments and return the value of
!   pfilt( r ) as output. 
!------------------------------------------------------------------------

#include "abbrev.inc"

#define _checkinit_  if ( .not. Q_modinit ) call lib_ppfilt_init


!////////////////////////////////////////////////////////////////////////
module lib_ppfilt
    ___
    save

    Dble :: alpha = 1.1d0   !! kdamp = kmax / alpha (with kmax = pi/h)
    Dble :: gamma = 2.0d0   !! rmask = gamma * rc(PP)

    !! Table of the mask function.
    Dble, _alloc_ :: mask_tab( : )
    Dble          :: dx_tab
    Int           :: Ngrid_tab

    !! Flag to use Wang's mask function for damping k-components.
    Bool :: use_maskfunc_k = .false.

    Bool :: Q_modinit = .false.

contains

!------------------------------------------------------------------------
sub lib_ppfilt_init

    call init_maskfunc

    Q_modinit = .true.

endsub

! Calc the value of the mask function via linear interpolation.
!------------------------------------------------------------------------
function get_maskfunc ( z, za, zb ) result ( ret )

    Dble, _in_        :: z
    Dble, _in_, _opt_ :: za, zb
    Dble              :: ret
    Int   i
    Dble  x
    !.........................................................

    _checkinit_

    if ( present( za ) .and. present( zb ) ) then
        if ( abs( zb - za ) < 1.0d-12 ) stop "za == zb"

        x = ( z - za ) / ( zb - za )
    else
        x = z
    _fi

    if ( x < 0.0d0 ) then
        ret = 1.0d0 ; return
    _fi
    if ( x > 1.0d0 ) then
        ret = 0.0d0 ; return
    _fi

    i = int( x / dx_tab )

    if ( i == -1 ) then
        i = 0
    elseif ( i == Ngrid_tab ) then
        i = Ngrid_tab - 1
    _fi

    ret = mask_tab( i ) &
            + ( mask_tab( i+1 ) - mask_tab( i ) ) / dx_tab &
            * ( x - dble( i ) * dx_tab )

endfunc

! Calc the damping factor for k-components.
!------------------------------------------------------------------------
function get_kdampfac ( k, ka, kb ) result ( ret )

    Dble, _in_ :: k, ka, kb
    Dble       :: ret

    Dble  tmp
    Dble, _param_ :: eps = 1.0d-5
    !.........................................................

    _checkinit_

    if ( use_maskfunc_k ) then

        ret = get_maskfunc( k, ka, kb )
    else
        if ( k < ka ) then
            ret = 1.0d0
        else
            tmp = -1.0d0 / ( kb - ka )**2 * log( eps )
            ret = exp( - tmp * ( k - ka )**2 )
            !! Note: this is not strictly zero even at k = kmax.
        _fi
    _fi

endfunc

! Table of Wang's mask function with eta = 15.
!------------------------------------------------------------------------
sub init_maskfunc

    Ngrid_tab = 200
    allocate( mask_tab( 0 : Ngrid_tab ) )

    dx_tab = 1.0d0 / dble( Ngrid_tab )

    mask_tab( 0 ) =  0.10000000d+01
    mask_tab( 1 ) =  0.10000000d+01
    mask_tab( 2 ) =  0.99948662d+00
    mask_tab( 3 ) =  0.99863154d+00
    mask_tab( 4 ) =  0.99743557d+00
    mask_tab( 5 ) =  0.99589985d+00
    mask_tab( 6 ) =  0.99402586d+00
    mask_tab( 7 ) =  0.99181538d+00
    mask_tab( 8 ) =  0.98927052d+00
    mask_tab( 9 ) =  0.98639370d+00
    mask_tab( 10 ) =  0.98318766d+00
    mask_tab( 11 ) =  0.97965544d+00
    mask_tab( 12 ) =  0.97580040d+00
    mask_tab( 13 ) =  0.97162618d+00
    mask_tab( 14 ) =  0.96713671d+00
    mask_tab( 15 ) =  0.96233623d+00
    mask_tab( 16 ) =  0.95722924d+00
    mask_tab( 17 ) =  0.95182053d+00
    mask_tab( 18 ) =  0.94611516d+00
    mask_tab( 19 ) =  0.94011842d+00
    mask_tab( 20 ) =  0.93383589d+00
    mask_tab( 21 ) =  0.92727338d+00
    mask_tab( 22 ) =  0.92043693d+00
    mask_tab( 23 ) =  0.91333282d+00
    mask_tab( 24 ) =  0.90596753d+00
    mask_tab( 25 ) =  0.89834777d+00
    mask_tab( 26 ) =  0.89048044d+00
    mask_tab( 27 ) =  0.88237263d+00
    mask_tab( 28 ) =  0.87403161d+00
    mask_tab( 29 ) =  0.86546483d+00
    mask_tab( 30 ) =  0.85667987d+00
    mask_tab( 31 ) =  0.84768450d+00
    mask_tab( 32 ) =  0.83848659d+00
    mask_tab( 33 ) =  0.82909416d+00
    mask_tab( 34 ) =  0.81951535d+00
    mask_tab( 35 ) =  0.80975838d+00
    mask_tab( 36 ) =  0.79983160d+00
    mask_tab( 37 ) =  0.78974340d+00
    mask_tab( 38 ) =  0.77950227d+00
    mask_tab( 39 ) =  0.76911677d+00
    mask_tab( 40 ) =  0.75859548d+00
    mask_tab( 41 ) =  0.74794703d+00
    mask_tab( 42 ) =  0.73718009d+00
    mask_tab( 43 ) =  0.72630334d+00
    mask_tab( 44 ) =  0.71532544d+00
    mask_tab( 45 ) =  0.70425508d+00
    mask_tab( 46 ) =  0.69310092d+00
    mask_tab( 47 ) =  0.68187158d+00
    mask_tab( 48 ) =  0.67057566d+00
    mask_tab( 49 ) =  0.65922170d+00
    mask_tab( 50 ) =  0.64781819d+00
    mask_tab( 51 ) =  0.63637355d+00
    mask_tab( 52 ) =  0.62489612d+00
    mask_tab( 53 ) =  0.61339415d+00
    mask_tab( 54 ) =  0.60187581d+00
    mask_tab( 55 ) =  0.59034914d+00
    mask_tab( 56 ) =  0.57882208d+00
    mask_tab( 57 ) =  0.56730245d+00
    mask_tab( 58 ) =  0.55579794d+00
    mask_tab( 59 ) =  0.54431609d+00
    mask_tab( 60 ) =  0.53286431d+00
    mask_tab( 61 ) =  0.52144984d+00
    mask_tab( 62 ) =  0.51007978d+00
    mask_tab( 63 ) =  0.49876105d+00
    mask_tab( 64 ) =  0.48750040d+00
    mask_tab( 65 ) =  0.47630440d+00
    mask_tab( 66 ) =  0.46517945d+00
    mask_tab( 67 ) =  0.45413176d+00
    mask_tab( 68 ) =  0.44316732d+00
    mask_tab( 69 ) =  0.43229196d+00
    mask_tab( 70 ) =  0.42151128d+00
    mask_tab( 71 ) =  0.41083069d+00
    mask_tab( 72 ) =  0.40025539d+00
    mask_tab( 73 ) =  0.38979038d+00
    mask_tab( 74 ) =  0.37944042d+00
    mask_tab( 75 ) =  0.36921008d+00
    mask_tab( 76 ) =  0.35910371d+00
    mask_tab( 77 ) =  0.34912542d+00
    mask_tab( 78 ) =  0.33927912d+00
    mask_tab( 79 ) =  0.32956851d+00
    mask_tab( 80 ) =  0.31999705d+00
    mask_tab( 81 ) =  0.31056799d+00
    mask_tab( 82 ) =  0.30128436d+00
    mask_tab( 83 ) =  0.29214897d+00
    mask_tab( 84 ) =  0.28316441d+00
    mask_tab( 85 ) =  0.27433307d+00
    mask_tab( 86 ) =  0.26565709d+00
    mask_tab( 87 ) =  0.25713844d+00
    mask_tab( 88 ) =  0.24877886d+00
    mask_tab( 89 ) =  0.24057988d+00
    mask_tab( 90 ) =  0.23254283d+00
    mask_tab( 91 ) =  0.22466884d+00
    mask_tab( 92 ) =  0.21695884d+00
    mask_tab( 93 ) =  0.20941357d+00
    mask_tab( 94 ) =  0.20203357d+00
    mask_tab( 95 ) =  0.19481920d+00
    mask_tab( 96 ) =  0.18777065d+00
    mask_tab( 97 ) =  0.18088790d+00
    mask_tab( 98 ) =  0.17417080d+00
    mask_tab( 99 ) =  0.16761900d+00
    mask_tab( 100 ) =  0.16123200d+00
    mask_tab( 101 ) =  0.15500913d+00
    mask_tab( 102 ) =  0.14894959d+00
    mask_tab( 103 ) =  0.14305240d+00
    mask_tab( 104 ) =  0.13731647d+00
    mask_tab( 105 ) =  0.13174055d+00
    mask_tab( 106 ) =  0.12632327d+00
    mask_tab( 107 ) =  0.12106315d+00
    mask_tab( 108 ) =  0.11595855d+00
    mask_tab( 109 ) =  0.11100775d+00
    mask_tab( 110 ) =  0.10620891d+00
    mask_tab( 111 ) =  0.10156010d+00
    mask_tab( 112 ) =  0.97059268d-01
    mask_tab( 113 ) =  0.92704295d-01
    mask_tab( 114 ) =  0.88492966d-01
    mask_tab( 115 ) =  0.84422989d-01
    mask_tab( 116 ) =  0.80492001d-01
    mask_tab( 117 ) =  0.76697569d-01
    mask_tab( 118 ) =  0.73037197d-01
    mask_tab( 119 ) =  0.69508335d-01
    mask_tab( 120 ) =  0.66108380d-01
    mask_tab( 121 ) =  0.62834685d-01
    mask_tab( 122 ) =  0.59684561d-01
    mask_tab( 123 ) =  0.56655284d-01
    mask_tab( 124 ) =  0.53744102d-01
    mask_tab( 125 ) =  0.50948236d-01
    mask_tab( 126 ) =  0.48264886d-01
    mask_tab( 127 ) =  0.45691239d-01
    mask_tab( 128 ) =  0.43224469d-01
    mask_tab( 129 ) =  0.40861744d-01
    mask_tab( 130 ) =  0.38600231d-01
    mask_tab( 131 ) =  0.36437098d-01
    mask_tab( 132 ) =  0.34369520d-01
    mask_tab( 133 ) =  0.32394681d-01
    mask_tab( 134 ) =  0.30509780d-01
    mask_tab( 135 ) =  0.28712032d-01
    mask_tab( 136 ) =  0.26998673d-01
    mask_tab( 137 ) =  0.25366964d-01
    mask_tab( 138 ) =  0.23814193d-01
    mask_tab( 139 ) =  0.22337676d-01
    mask_tab( 140 ) =  0.20934765d-01
    mask_tab( 141 ) =  0.19602844d-01
    mask_tab( 142 ) =  0.18339338d-01
    mask_tab( 143 ) =  0.17141711d-01
    mask_tab( 144 ) =  0.16007467d-01
    mask_tab( 145 ) =  0.14934157d-01
    mask_tab( 146 ) =  0.13919377d-01
    mask_tab( 147 ) =  0.12960772d-01
    mask_tab( 148 ) =  0.12056034d-01
    mask_tab( 149 ) =  0.11202905d-01
    mask_tab( 150 ) =  0.10399183d-01
    mask_tab( 151 ) =  0.96427132d-02
    mask_tab( 152 ) =  0.89313983d-02
    mask_tab( 153 ) =  0.82631938d-02
    mask_tab( 154 ) =  0.76361106d-02
    mask_tab( 155 ) =  0.70482151d-02
    mask_tab( 156 ) =  0.64976294d-02
    mask_tab( 157 ) =  0.59825322d-02
    mask_tab( 158 ) =  0.55011581d-02
    mask_tab( 159 ) =  0.50517982d-02
    mask_tab( 160 ) =  0.46327998d-02
    mask_tab( 161 ) =  0.42425662d-02
    mask_tab( 162 ) =  0.38795566d-02
    mask_tab( 163 ) =  0.35422853d-02
    mask_tab( 164 ) =  0.32293218d-02
    mask_tab( 165 ) =  0.29392897d-02
    mask_tab( 166 ) =  0.26708663d-02
    mask_tab( 167 ) =  0.24227820d-02
    mask_tab( 168 ) =  0.21938194d-02
    mask_tab( 169 ) =  0.19828122d-02
    mask_tab( 170 ) =  0.17886449d-02
    mask_tab( 171 ) =  0.16102512d-02
    mask_tab( 172 ) =  0.14466132d-02
    mask_tab( 173 ) =  0.12967606d-02
    mask_tab( 174 ) =  0.11597692d-02
    mask_tab( 175 ) =  0.10347601d-02
    mask_tab( 176 ) =  0.92089812d-03
    mask_tab( 177 ) =  0.81739110d-03
    mask_tab( 178 ) =  0.72348823d-03
    mask_tab( 179 ) =  0.63847906d-03
    mask_tab( 180 ) =  0.56169212d-03
    mask_tab( 181 ) =  0.49249371d-03
    mask_tab( 182 ) =  0.43028657d-03
    mask_tab( 183 ) =  0.37450862d-03
    mask_tab( 184 ) =  0.32463165d-03
    mask_tab( 185 ) =  0.28016004d-03
    mask_tab( 186 ) =  0.24062948d-03
    mask_tab( 187 ) =  0.20560566d-03
    mask_tab( 188 ) =  0.17468305d-03
    mask_tab( 189 ) =  0.14748362d-03
    mask_tab( 190 ) =  0.12365560d-03
    mask_tab( 191 ) =  0.10287226d-03
    mask_tab( 192 ) =  0.84830727d-04
    mask_tab( 193 ) =  0.69250769d-04
    mask_tab( 194 ) =  0.55873673d-04
    mask_tab( 195 ) =  0.44461100d-04
    mask_tab( 196 ) =  0.34793983d-04
    mask_tab( 197 ) =  0.26671449d-04
    mask_tab( 198 ) =  0.19909778d-04
    mask_tab( 199 ) =  0.14341381d-04
    mask_tab( 200 ) =  0.98138215d-05

endsub

!------------------------------------------------------------------------
sub plot_maskfunc

    Int   i
    Dble  x
    !.........................................................

    _checkinit_

    put "file(out): mask.dat, mask_fine.dat"

    open( 10, file="mask.dat" )
    open( 20, file="mask_fine.dat" )

    !! Mask func (raw data).
    do i = 0, Ngrid_tab
        write( 10, "(e15.8, 1x, e15.8)" ) dble( i ) * dx_tab, mask_tab( i )
                   !...Format for comparison with the original data file.
    enddo

    !! Mask func on a finer grid.
    do i = 0, Ngrid_tab * 5
        x = dble( i ) * dx_tab / 5.0d0
        write( 20, "(3es20.10)" ) &
                x, &
                get_maskfunc( x ), &
                get_maskfunc( x, 0.0d0, 1.0d0 )
    enddo

    close( 10 )
    close( 20 )

endsub

!------------------------------------------------------------------------
sub plot_kdampfac

    Int   i, Ndiv
    Dble  k, kmax, kdamp
    !.........................................................

    _checkinit_

    put "file(out): kdamp_{ gauss, mask }.dat"

    open( 10, file="kdamp_gauss.dat" )
    open( 20, file="kdamp_mask.dat" )

    kmax = 10.0d0
    Ndiv = 1000

    kdamp = kmax / alpha

    do i = 1, Ndiv * 2

        k = dble( i ) * kmax / dble( Ndiv )

        use_maskfunc_k = .false.
        write( 10, "(2es20.10)" ) k, max( get_kdampfac( k, kdamp, kmax ), 1.0d-15 )

        use_maskfunc_k = .true.
        write( 20, "(2es20.10)" ) k, get_kdampfac( k, kdamp, kmax )
    enddo

    close( 10 )
    close( 20 )

endsub

endmodule


!////////////////////////////////////////////////////////////////////////
! Unit tests.
!////////////////////////////////////////////////////////////////////////

#ifdef lib_ppfilt_test


!========================================================================
program main
    use lib_ppfilt
    ___

    call lib_ppfilt_init

    call plot_maskfunc
    call plot_kdampfac

end program

#endif  /* lib_ppfilt_test */

