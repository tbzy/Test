
!////////////////////////////////////////////////////////////////////////
!! Abbrev words for Fortran (Last updated: 2014-9-2)

#ifndef __ABBREV_INC__
#define __ABBREV_INC__

#define  _bool_         logical
#define  _char_         character
#define  _int_          integer
#define  _real_         real(8)
#define  _dble_         real(8)
#define  _complex_      complex(8)

#define  _Bool_         logical
#define  _Char_         character
#define  _Int_          integer
#define  _Real_         real(8)
#define  _Dble_         real(8)
#define  _Complex_      complex(8)

#define  _sub_          subroutine
#define  _endsub_       end subroutine
#define  _func_         function
#define  _endfunc_      end function
#define  _endmod_       end module
#define  _put_          write(6,*)
#define  _printf_(x)    write(6,x)
#define  _using_(x)     use x, only:
#define  _switch_       selectcase
#define  _endsw_        endselect

#define  ___          implicit none
#define  _in_         intent(in)
#define  _alloc_      allocatable
#define  _param_      parameter
#define  _opt_        optional
#define  _dim_        dimension
#define  _ptr_        pointer
#define  _proc_       procedure
#define  _(x)         x = x
#define  __(x)        type(x)
#define  _fi          endif

#define  _norm_(x)              sqrt(sum((x)**2))
#define  _normalize_(x)         ((x) / _norm_(x))
#define  _proj_(p,x)            ((p) * sum((p)*(x)) / sum((p)**2))
#define  _limit_(x,xmin,xmax)   x = min(max(x,xmin),xmax)
#define  _pbc_(x,box)           x = x - (box) * nint((x) / max(0.1d0,(box)))
! #define  _pbc_(x,box)           x = x - (box) * nint((x) / (box))
#define  _dealloc_(x)           if (allocated(x)) deallocate(x)
#define  _optval_(x,xopt,val)   x = val; if (present(xopt)) x = xopt
#define  _onstep_(i,n,div)      (div > 0 .and. mod(i, max((n)/max(div,1),1)) == 0)

#define  _Re_(x)  dble(x)     /* dreal(x) may also be used */
#define  _Im_(x)  dimag(x)

#define  _line_           repeat("-",70)
#define  _dline_          repeat("=",70)
#define  _section_(x)     put; put _dline_; put x
#define  _routine_(x)     put; put repeat(".",60); put x, ":" ; put
#define  _flush_          flush(6)
#define  _stop_(x)        stop __FILE__ // " : " // x

#define  _0_        0.0d0
#define  _zi_       ( 0.0d0, 1.0d0 )
#define  _pi_       3.1415926535897931d0
#define  _twopi_    6.2831853071795862d0
#define  _fourpi_   12.566370614359172d0
#define  _sqrt2_    1.4142135623730951d0
#define  _sqrt3_    1.7320508075688772d0
#define  _sqrtpi_   1.7724538509055159d0
#define  _sqrt2pi_  2.5066282746310002d0

#define  _atomtype_t        Char(4)   /* Amber=>4, DLPOLY=>8 (8 may be better) */
#define  _atomname_t        Char(4)
#define  _pottype_t         Char(8)   /* potential type */
#define  _filename_t        Char(50)

!! [ Physical units ]
!! #define  _au_to_kcal_         627.509d0     /* latest value */
#define  _au_to_kcal_         627.488d0       /* for comparison with Amber */
#define  _au_to_eV_           27.21138d0
#define  _au_to_Angs_         0.5291772d0
#define  _au_to_amu_          5.4858d-4
#define  _au_to_fsec_         0.0241888d0
#define  _au_to_psec_         0.0241888d-3
#define  _au_to_nsec_         0.0241888d-6
#define  _au_to_wavenum_      219474.6d0
#define  _rad_to_deg_         ( 180.0d0 / _pi_ )
#define  _kB_in_au_           3.166815d-6     /* Boltzmann constant */

#define  _au_to_kcalperAngs_   ( _au_to_kcal_ / _au_to_Angs_ )      /* force */
#define  _au_to_kcalperAngs2_  ( _au_to_kcal_ / _au_to_Angs_**2 )   /* force const */
#define  _au_to_kcalAngs12_    ( _au_to_kcal_ * _au_to_Angs_**12 )  /* LJ A coef */
#define  _au_to_kcalAngs6_     ( _au_to_kcal_ * _au_to_Angs_**6 )   /* LJ B coef */

!! [ Amber internal units ]
#define _amb_to_psec_  ( 1.0d0 / 20.455d0 )
#define _kB_in_amb_    ( 8.31441d-3 / 4.184d0 )

!! [ 1-4 scale factor (default: Amber) ]
#define _scale_es14_  ( 1.0d0 / 1.2d0 )
#define _scale_lj14_  ( 1.0d0 / 2.0d0 )

!! [ Misc ]

#define _f2003_   /* use Fortran2003 features */

!* #define __extern__  /* define this when subroutines are put in external scope */

!! [ Debug ]

#define _debug_ .true.

! #define  __log(x)    put "x", " = ", x ;
! #define  __log(x)      /* no debug output in production run */

#endif /* __abbrev_inc__ */
!////////////////////////////////////////////////////////////////////////

