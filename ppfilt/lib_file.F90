!------------------------------------------------------------------------------
! lib_file.F90  (Last updated: 2014-10-11)
! 
! File handling.
!
! 2014-10-11: xyz reader/writerをアップデートした。Q_EOF_などを追加。
! 2014-10-10: xyzfile_write()にQ_write_box_オプションを追加。
!             これがTだとtitleの欄にbox情報を書く。
! 2014-09-15: xyzfile_t, xyzfile_read(), xyzfile_write() を追加した。
! 2014-09-11: free_unit() => get_free_unit()
! 2011-06-17: status_="new"でも古いファイルを強制的には消さないようにした
!------------------------------------------------------------------------------

#include "abbrev.inc"

#define atomname_t  Char(8)
#define title_t     Char(100)

#define  __writelog  write(logfile,*)


!/////////////////////////////////////////////////////////////////////////////
module lib_file
    ___
    save

    !-------------------------
    type file_t

        Char(200) :: name = "none"   !! file name

        Int :: u = 0      !! file unit

    endtype
    !-------------------------

    Bool :: file_verbose = .true.   !! if true, print file names to stdout

    Int :: logfile = 6     !! a different unit may be used to keep stdout clean

    interface to_str
        module procedure  int_to_str, dble_to_str
    endinterface

    !.........................................................

    !------------------------------
    type xyzfile_t

        Int :: Natom = 0
        
        title_t :: title = ""

        atomname_t, _alloc_ :: atomname(:)

        Dble, _alloc_ :: crd(:,:)   !! xyz coordinates of size (3,Natom)

        Dble :: box( 3 ) = _0_   !! box side lengths

    endtype
    !------------------------------

    !.........................................................

#ifdef _f2003_
    type lib_file_t
        contains
            _proc_, nopass :: file_open
            _proc_, nopass :: get_free_unit
            _proc_, nopass :: read_xyz
            _proc_, nopass :: write_xyz
    endtype

    __(lib_file_t) :: lib_file_
#endif

contains


!========================================================================
sub file_open ( file, name, &
                unit_, status_, form_ )
    ___
    __(file_t), target :: file

    Char(*), _in_         :: name
    Int,     _in_, _opt_  :: unit_
    Char(*), _in_, _opt_  :: status_, form_

    Char(50)  status, form, mode, s
    !...................................................................

    file% name = name

!! File unit.

    if ( present( unit_ ) ) then
        file% u = unit_
    else
        file% u = get_free_unit()
    _fi

!! File status.

    _optval_( status, status_, "unknown" )

    !-------------------------
    select case ( status )

    case ( "old" )
        mode = "in"

    case ( "new" )
        mode = "out"
        s = "rm " // trim( file% name )
        ! put trim( s )
        ! call system ( s )   !! remove an existing file

    case default
        mode = "out"

    endselect
    !-------------------------

    if ( file_verbose ) then
        __writelog
        __writelog "file(", trim( mode ), "): ", trim( file% name )
    _fi

!...Format.

    _optval_( form, form_, "formatted" )

!...Open file.

    open( unit   = file% u,  &
          file   = trim( file% name ),  &
          status = trim( status ),  &
          form   = trim( form )  &
        )

endsub


!========================================================================
function get_free_unit() result( funit )

!! returns a free file unit.
    ___
    Int :: funit

    Bool  opened
    Int, save :: base = 5000
    !...................................................................

    funit = base
    do
        inquire( unit = funit, opened = opened )
        if ( .not. opened ) then
            base = funit + 1   !! to allow consecutive call of this routine
            exit
        _fi
        _( funit ) + 1
    enddo

endfunc


!========================================================================
sub get_last_line( funit, line )
    ___
    Int, _in_ :: funit
    Char(200) :: line
    !.........................................................

    rewind( funit )
    do
        read( funit, "(a200)", end=100 ) line
    enddo

100 continue

endsub


!========================================================================
function int_to_str( n, fmt ) result( str )
    ___
    Int,     _in_        :: n
    Char(*), _in_, _opt_ :: fmt

    Char(30) :: str
    !...................................................................

    if ( present( fmt ) ) then
        write( str, "(" // fmt // ")" ) n
    else
        write( str, "(i30)" ) n
    _fi

    str = adjustL( str )

endfunc


!========================================================================
function dble_to_str( x, fmt ) result( str )
    ___
    Dble,    _in_        :: x
    Char(*), _in_, _opt_ :: fmt

    Char(30) :: str
    !...................................................................

    if ( present( fmt ) ) then
        write( str, "(" // fmt // ")" ) x
    else
        write( str, "(es30.20)" ) x
    _fi

    str = adjustL( str )

endfunc


!========================================================================
function uppercase( str ) result( str2 )
    ___
    Char(*), _in_ :: str
    Char( len( str ) )  str2
    Int  i
    !.........................................................

    str2 = str

    do i = 1, len_trim( str2 )
        if ( "a" <= str2(i:i) .and. str2(i:i) <= "z" ) then
            str2(i:i) = char( ichar( str2(i:i) ) - 32 )  !! primitive way
        _fi
    enddo

endfunc


!========================================================================
function lowercase( str ) result( str2 )
    ___
    Char(*), _in_ :: str
    Char( len( str ) )  str2
    Int  i
    !.........................................................

    str2 = str

    do i = 1, len_trim( str2 )
        if ( "A" <= str2(i:i) .and. str2(i:i) <= "Z" ) then
            str2(i:i) = char( ichar( str2(i:i) ) + 32 )  !! primitive way
        _fi
    enddo

endfunc


!========================================================================
sub read_xyz ( inp, crd, Natom, &
               atomname_, title_, box_, &
               scale_, Q_EOF_ )

!! Low-level xyz reader.
!! Notes:
!! - The input file should be open upon entry.
!! - Natom may be smaller than the total number of atoms in a file.
!! - One should *not* rewind the file below so that
!!   one can read multiple xyz data sequentially.
    ___
    Int, _in_ :: inp              !! input file unit
    Int, _in_ :: Natom            !! number of atoms to read
    Dble      :: crd( 3, Natom )  !! xyz coord [out]

    Char(*), _opt_ :: atomname_( Natom )  !! [out]
    Char(*), _opt_ :: title_              !! [out]
    Dble,    _opt_ :: box_( 3 )           !! [out]
    Dble,    _opt_ :: scale_
    Bool,    _opt_ :: Q_EOF_

    Int         Natom_all, ia
    atomname_t  atomname( Natom )
    title_t     title, stmp
    !.........................................................

    read( inp, *, end=9000 ) Natom_all
    if ( Natom > Natom_all ) _stop_( "read_xyz: Natom too large" )

    read( inp, "(a)" ) title

    do ia = 1, Natom_all
        if ( ia <= Natom ) then
            read( inp, * ) atomname( ia ), crd( :, ia )
        else
            read( inp, * )  !! dummy read
        _fi
    enddo

    if ( present( title_    ) ) title_ = title
    if ( present( atomname_ ) ) atomname_(:) = atomname(:)

    !! Read box side lengths from the title (optional).
    if ( present( box_ ) ) then
        read( title, * ) stmp   !! skip possible blank columns
        if ( stmp( 1:3 ) == "box" ) read( title, * ) stmp, box_( 1:3 )
    _fi

    !! Scale crd and box (optional).
    if ( present( scale_ ) ) then
        _( crd( :, : ) ) * scale_
        if ( present( box_ ) ) _( box_( : ) ) * scale_
    _fi

    if ( present( Q_EOF_ ) ) Q_EOF_ = .false.
    return

9000 continue
    if ( present( Q_EOF_ ) ) Q_EOF_ = .true.

endsub


!========================================================================
sub write_xyz ( out, crd, Natom, &
                atomname_, title_, box_, scale_ )

!! write XYZ data on an already opened file.
    ___
    Int,  _in_  :: out       !! output file unit
    Int,  _in_  :: Natom
    Dble, _in_  :: crd( 3, Natom )

    Char(*), _in_, _opt_ :: atomname_( Natom )
    Char(*), _in_, _opt_ :: title_
    Dble,    _in_, _opt_ :: box_( 3 )
    Dble,    _in_, _opt_ :: scale_

    Int         ia
    Dble        scale
    atomname_t  atomname( Natom )
    title_t     title
    !.........................................................

    _optval_( atomname, atomname_, "X"    )
    _optval_( title,    title_,    ""     )
    _optval_( scale,    scale_,    1.0d0  )

    if ( present( box_ ) ) then
        write( title, "('box  ', 3f20.8)" ) box_( 1 : 3 ) * scale
    _fi

    write( out, "(i0)" ) Natom
    write( out, "(a)" ) trim( adjustL( title ) )

    do ia = 1, Natom
        write( out, "(a8, 2x, 3f20.8)" )  &
                adjustL( atomname( ia ) ), &
                crd( :, ia ) * scale
    enddo

endsub


!========================================================================
sub xyzfile_read ( my, inp, &
                   scale_, Q_EOF_ )

!! xyz file reader (a wrapper for read_xyz()).
!! This routine reads one block of xyz data.
    ___
    __(xyzfile_t) :: my      !! this object (xyz data) [out]
    Int, _in_     :: inp     !! input file unit
    Dble, _opt_   :: scale_
    Bool, _opt_   :: Q_EOF_

    Int   Natom
    Bool  Q_alloc
    !.........................................................

    !! Get Natom for possible memory realloc.
    read( inp, *, end=9000 ) Natom
    backspace( inp )

    my% Natom = Natom

!! Realloc memory if necessary.

    Q_alloc = .false.

    if ( allocated( my% crd ) ) then

        if ( size( my% crd, 2 ) /= Natom ) then
            deallocate( my% crd, my% atomname )
            Q_alloc = .true.
        _fi
    else
        Q_alloc = .true.
    _fi

    if ( Q_alloc ) then
        allocate( my% crd( 3, Natom ), &
                  my% atomname( Natom ) )
    _fi

!! Read one block of xyz data.

    call read_xyz &
            ( inp, my% crd, my% Natom, &
              atomname_ = my% atomname, &
              title_    = my% title, &
              box_      = my% box, &
              scale_    = scale_ )

    if ( present( Q_EOF_ ) ) Q_EOF_ = .false.
    return

9000 continue
    if ( present( Q_EOF_ ) ) Q_EOF_ = .true.

endsub


!========================================================================
sub xyzfile_write ( my, out, &
                    scale_, Q_write_box_ )

!! wrapper for write_xyz().
!! This routine writes one block of xyz data.
    ___
    __(xyzfile_t), _in_ :: my       !! this object (xyz data to be written)
    Int, _in_           :: out      !! output file unit
    Dble, _opt_         :: scale_
    Bool, _opt_         :: Q_write_box_

    Bool  Q_write_box
    !.........................................................

    _optval_( Q_write_box, Q_write_box_, .false. )

    if ( Q_write_box ) then
        call write_xyz &
                ( out, my% crd, my% Natom, &
                  atomname_ = my% atomname, &
                  box_      = my% box, &
                  scale_    = scale_ )
    else
        call write_xyz &
                ( out, my% crd, my% Natom, &
                  atomname_ = my% atomname, &
                  title_    = my% title, &
                  scale_    = scale_ )
    _fi

endsub


endmodule


!/////////////////////////////////////////////////////////////////////////////
! Unit tests.
!/////////////////////////////////////////////////////////////////////////////

#ifdef lib_file_test


!========================================================================
program main

    use lib_file

    ! call test_free_unit
    ! call test_file_open

    ! call test_xyz
    call test_xyzfile

    ! call test_upperlowercase

contains


!========================================================================
sub test_file_open()
    ___
    __(file_t) :: file

    Char(80)  str
    !...................................................................

!...Write.

    call file_open ( file, "test.dat" )
    write( file% u, * ) "this is test.dat"
    close( file% u )

!...Read.

    call file_open ( file, "test.dat", status_= "old" )
    read( file% u, "(a80)" ) str
    close( file% u )

    put "str = ", trim( str )

endsub


!========================================================================
sub test_free_unit()
    ___
    Int      ifile, unit
    Char(1)  answer
    !...................................................................

    do ifile = 1, 5

        unit = get_free_unit()
        put "ifile = ", ifile, " unit = ", unit

        put "open a new file ? [y/n]"
        read(*,*) answer

        if ( answer(1:1) == "y" ) then
            open( unit, file = "test" // trim( to_str( ifile ) ) // ".dat" )
        _fi

    enddo

endsub


!========================================================================
sub write_sample_xyz_file ( filename )
    ___
    Char(*), _in_ :: filename
    !.........................................................

    put "file(out): ", trim(filename)
    open( 10, file=trim(filename) )

    write( 10, "(a)" ) "2"
    write( 10, "(a)" ) "hydrogen molecule"
    write( 10, "(a)" ) "H1    1.0    2.0    3.0"
    write( 10, "(a)" ) "H2    4.0    5.0    6.0"
    write( 10, "(a)" ) "3"
    write( 10, "(a)" ) "box  -777.0  -888.0  -999.0"
    write( 10, "(a)" ) "C1   -1.0   -2.0   -3.0"
    write( 10, "(a)" ) "C2   -4.0   -5.0   -6.0"
    write( 10, "(a)" ) "C3   -7.0   -8.0   -9.0"

    close( 10 )

endsub


!========================================================================
sub test_xyz ()

!! test low-level I/O routines for xyz file.
    ___
    Int, _param_ :: mxatom = 100

    title_t   title1
    Char(4)   atomname1( mxatom )
    Char(8)   atomname2( mxatom )   !! different lengths are used intentionally

    Dble, _dim_( 3, mxatom ) :: crd1, crd2
    Dble   box2( 3 )
    Int    Natom1, Natom2
    !.........................................................

    call write_sample_xyz_file ( "test.xyz" )

!! Read xyz file.

    open( 10, file="test.xyz", status="old" )

    Natom1 = 1   !! could be smaller than the actual number of atoms
    call read_xyz ( 10, crd1, Natom1, &
                    atomname_ = atomname1, &
                    title_    = title1 )

    Natom2 = 3
    call read_xyz ( 10, crd2, Natom2, &
                    atomname_ = atomname2, &
                    box_      = box2 )

    close( 10 )

!! Write xyz file.

    open( 10, file= "test2.xyz" )

    call write_xyz ( 10, crd1, Natom1, &
                     atomname_ = atomname1, &
                     title_    = title1 )

    atomname2( 1 ) = "carbon1"
    call write_xyz ( 10, crd2, Natom2, &
                     atomname_ = atomname2, &
                     box_      = box2, &
                     scale_    = -1.0d0 )

    close( 10 )

!! Check.

    put
    put "=== test.xyz ==="
    call system( "cat test.xyz" )
    put
    put "=== test2.xyz ==="
    call system( "cat test2.xyz" )

endsub


!========================================================================
sub test_xyzfile ()
    ___
    __(xyzfile_t) :: xyz
    Int   inp, out
    Bool  Q_EOF
    !.........................................................

    call write_sample_xyz_file ( "test_inp.xyz" )

    inp = get_free_unit()
    out = get_free_unit()

    open( inp, file="test_inp.xyz", status="old" )
    open( out, file="test_out.xyz" )

    !! Read/write the first block.
    call xyzfile_read  ( xyz, inp )
    call xyzfile_write ( xyz, out )

    !! Read the second block.
    call xyzfile_read  ( xyz, inp )

    !! Write it with different box info.
    xyz% box(:) = [ 1111.0d0, 2222.0d0, 3333.0d0 ]
    call xyzfile_write ( xyz, out, &
                         Q_write_box_=.true. )

    rewind( inp )

    !! Sequential read/write until EOF.
    do
        call xyzfile_read ( xyz, inp, Q_EOF_= Q_EOF )
        if ( Q_EOF ) exit

        call xyzfile_write ( xyz, out )
    enddo

    close( inp )
    close( out )

!! Check.

    put
    put "=== test_inp.xyz ==="
    call system( "cat test_inp.xyz" )
    put
    put "=== test_out.xyz ==="
    call system( "cat test_out.xyz" )

endsub


!========================================================================
sub test_upperlowercase()
    ___
    Char(70) :: str

    str = "abcDEFG hijklMN OPQRstu vwXYZ 1234567890 !#$%&()=-^~\|{}<>?_[]"

    write( *, "(a)" ) str
    write( *, "(a)" ) uppercase( str )
    write( *, "(a)" ) lowercase( str )

endsub


end program

#endif /* unit tests */
!////////////////////////////////////////////////////////////////////////



