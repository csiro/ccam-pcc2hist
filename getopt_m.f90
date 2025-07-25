! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
module getopt_m

   implicit none
   private


!!!  F90 version of getopt, modified from GNU glibc-2.2.2 getopt.c
!!!  Comments from original C code have just a leading !. New comments
!!!  have !!!

!!!  Translated by Martin.Dix@csiro.au.

!!!  Original copyright for getopt.c

!    Copyright (C) 1987, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 2000
!    	Free Software Foundation, Inc.
!
!    The GNU C Library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Library General Public License as
!    published by the Free Software Foundation; either version 2 of the
!    License, or (at your option) any later version.
!
!    The GNU C Library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Library General Public License for more details.
!
!    You should have received a copy of the GNU Library General Public
!    License along with the GNU C Library; see the file COPYING.LIB.  If not,
!    write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!    Boston, MA 02111-1307, USA.

!!!  This version implements only the "POSIXLY_CORRECT" case where the options
!!!  must precede the arguments.

!!!  Note that fortran versions of iargc() all return the number of command
!!!  line arguments. The C version returns the number including the
!!!  program name, so is one larger than the fortran value.
!!!  There are remaining options if optind <= iargc().
!!!  After the end, optind is the index in argv of the first element that's
!!!  not an option, i.e. getarg(optind) returns the first non option.


!  This is also the maximum length of the command line and the optstring
   integer, parameter, public :: MAX_ARGLEN=256

!  For communication from `getopt' to the caller.
!  When `getopt' finds an option that takes an argument,
!  the argument value is returned here.
!   character(len=MAX_ARGLEN), public :: optarg

!  Index in ARGV of the next element to be scanned.
!  This is used for communication to and from the caller
!  and for communication between successive calls to `getopt'.
!
!  On entry to `getopt', zero means this is the first call; initialize.
!
!  When `getopt' returns -1, this is the index of the first of the
!  non-option elements that the caller should itself scan.
!
!!$!  Otherwise, `optind' communicates from one call to the next
!!$!  how much of ARGV has been scanned so far.
!!$!
!!$!  1003.2 says this must be 1 before any call.
!!$   integer, private :: optind = 1

!  Formerly, initialization of getopt depended on optind==0, which
!  causes problems with re-calling getopt as programs generally don't
!  know that.



!  Callers store zero here to inhibit the error message
!  for unrecognized options.

   logical, public :: opterr = .true.

!  Set to an option character which was unrecognized.
!  This must be initialized on some systems to avoid linking in the
!  system's own getopt implementation.

   integer, public :: optopt = ichar("?")

!!! Values for has_arg
   integer, public, parameter :: no_argument = 0,           &
                                 required_argument = 1,     &
                                 optional_argument = 2

!!  Flag specifies how values are returned. If flag=0, then getopt returns
!!  val. Otherwise it returns 0 and flag is set to val. Note that optind
!!  is set to the index of the option too.
   ! This seems unreasonably complicated. It whould be enough just to return
   ! val, which can be zero or a character
   type, public :: loption
      character(len=MAX_ARGLEN) :: name
      integer :: has_arg
      integer :: val
   end type loption

   character(len=*), public, parameter :: &
        getopt_revision = "$Revision: 2.5 $"

   public :: getopt, getcline

! /* Describe how to deal with options that follow non-option ARGV-elements.

!    If the caller did not specify anything,
!    the default is REQUIRE_ORDER if the environment variable
!    POSIXLY_CORRECT is defined, PERMUTE otherwise.

!    REQUIRE_ORDER means don't recognize them as options;
!    stop option processing when the first non-option is seen.
!    This is what Unix does.
!    This mode of operation is selected by either setting the environment
!    variable POSIXLY_CORRECT, or using `+' as the first character
!    of the list of option characters.

!    PERMUTE is the default.  We permute the contents of ARGV as we scan,
!    so that eventually all the non-options are at the end.  This allows options
!    to be given in any order, even with programs that were not written to
!    expect this.

!    RETURN_IN_ORDER is an option available to programs that were written
!    to expect options and other ARGV-elements in any order and that care about
!    the ordering of the two.  We describe each non-option ARGV-element
!    as if it were the argument of an option with character code 1.
!    Using `-' as the first character of the list of option characters
!    selects this mode of operation.

!    The special argument `--' forces an end of option-scanning regardless
!    of the value of `ordering'.  In the case of RETURN_IN_ORDER, only
!    `--' can cause `getopt' to return -1 with `optind' != ARGC.  */

! static enum
! {
!   REQUIRE_ORDER, PERMUTE, RETURN_IN_ORDER
! } ordering;

! /* Value of POSIXLY_CORRECT environment variable.  */
! static char *posixly_correct;
! 
! #ifdef	__GNU_LIBRARY__
! /* We want to avoid inclusion of string.h with non-GNU libraries
!    because there are many ways it can cause trouble.
!    On some systems, it contains special magic macros that don't work
!    in GCC.  */
! # include <string.h>
! # define my_index	strchr
! #else

! # if HAVE_STRING_H
! #  include <string.h>
! # else
! #  include <strings.h>
! # endif

! /* Avoid depending on library functions or files
!    whose names are inconsistent.  */

! #ifndef getenv
! extern char *getenv ();
! #endif

! static char *
! my_index (str, chr)
!      const char *str;
!      int chr;
! {
!   while (*str)
!     {
!       if (*str == chr)
! 	return (char *) str;
!       str++;
!     }
!   return 0;
! }

! /* If using GCC, we can safely declare strlen this way.
!    If not using GCC, it is ok not to declare it.  */
! #ifdef __GNUC__
! /* Note that Motorola Delta 68k R3V7 comes with GCC but not stddef.h.
!    That was relevant to code that was here before.  */
! # if (!defined __STDC__ || !__STDC__) && !defined strlen
! /* gcc with -traditional declares the built-in strlen to return int,
!    and has done so at least since version 2.4.5. -- rms.  */
! extern int strlen (const char *);
! # endif /* not __STDC__ */
! #endif /* __GNUC__ */

! #endif /* not __GNU_LIBRARY__ */
! 
! /* Handle permutation of arguments.  */

! /* Describe the part of ARGV that contains non-options that have
!    been skipped.  `first_nonopt' is the index in ARGV of the first of them;
!    `last_nonopt' is the index after the last of them.  */

! static int first_nonopt;
! static int last_nonopt;

! #ifdef _LIBC
! /* Bash 2.0 gives us an environment variable containing flags
!    indicating ARGV elements that should not be considered arguments.  */

! /* Defined in getopt_init.c  */
! extern char *__getopt_nonoption_flags;

! static int nonoption_flags_max_len;
! static int nonoption_flags_len;

! static int original_argc;
! static char *const *original_argv;

! /* Make sure the environment variable bash 2.0 puts in the environment
!    is valid for the getopt call we must make sure that the ARGV passed
!    to getopt is that one passed to the process.  */
! static void
! __attribute__ ((unused))
! store_args_and_env (int argc, char *const *argv)
! {
!   /* XXX This is no good solution.  We should rather copy the args so
!      that we can compare them later.  But we must not use malloc(3).  */
!   original_argc = argc;
!   original_argv = argv;
! }
! # ifdef text_set_element
! text_set_element (__libc_subinit, store_args_and_env);
! # endif /* text_set_element */

! # define SWAP_FLAGS(ch1, ch2) \
!   if (nonoption_flags_len > 0)						      \
!     {									      \
!       char __tmp = __getopt_nonoption_flags[ch1];			      \
!       __getopt_nonoption_flags[ch1] = __getopt_nonoption_flags[ch2];	      \
!       __getopt_nonoption_flags[ch2] = __tmp;				      \
!     }
! #else	/* !_LIBC */
! # define SWAP_FLAGS(ch1, ch2)
! #endif	/* _LIBC */

! /* Exchange two adjacent subsequences of ARGV.
!    One subsequence is elements [first_nonopt,last_nonopt)
!    which contains all the non-options that have been skipped so far.
!    The other is elements [last_nonopt,optind), which contains all
!    the options processed since those non-options were skipped.

!    `first_nonopt' and `last_nonopt' are relocated so that they describe
!    the new indices of the non-options in ARGV after they are moved.  */

! #if defined __STDC__ && __STDC__
! static void exchange (char **);
! #endif

! static void
! exchange (argv)
!      char **argv;
! {
!   int bottom = first_nonopt;
!   int middle = last_nonopt;
!   int top = optind;
!   char *tem;

!   /* Exchange the shorter segment with the far end of the longer segment.
!      That puts the shorter segment into the right place.
!      It leaves the longer segment in the right place overall,
!      but it consists of two parts that need to be swapped next.  */

! #ifdef _LIBC
!   /* First make sure the handling of the `__getopt_nonoption_flags'
!      string can work normally.  Our top argument must be in the range
!      of the string.  */
!   if (nonoption_flags_len > 0 && top >= nonoption_flags_max_len)
!     {
!       /* We must extend the array.  The user plays games with us and
! 	 presents new arguments.  */
!       char *new_str = malloc (top + 1);
!       if (new_str == NULL)
! 	nonoption_flags_len = nonoption_flags_max_len = 0;
!       else
! 	{
! 	  memset (__mempcpy (new_str, __getopt_nonoption_flags,
! 			     nonoption_flags_max_len),
! 		  '\0', top + 1 - nonoption_flags_max_len);
! 	  nonoption_flags_max_len = top + 1;
! 	  __getopt_nonoption_flags = new_str;
! 	}
!     }
! #endif

!   while (top > middle && middle > bottom)
!     {
!       if (top - middle > middle - bottom)
! 	{
! 	  /* Bottom segment is the short one.  */
! 	  int len = middle - bottom;
! 	  register int i;

! 	  /* Swap it with the top part of the top segment.  */
! 	  for (i = 0; i < len; i++)
! 	    {
! 	      tem = argv[bottom + i];
! 	      argv[bottom + i] = argv[top - (middle - bottom) + i];
! 	      argv[top - (middle - bottom) + i] = tem;
! 	      SWAP_FLAGS (bottom + i, top - (middle - bottom) + i);
! 	    }
! 	  /* Exclude the moved bottom segment from further swapping.  */
! 	  top -= len;
! 	}
!       else
! 	{
! 	  /* Top segment is the short one.  */
! 	  int len = top - middle;
! 	  register int i;

! 	  /* Swap it with the bottom part of the bottom segment.  */
! 	  for (i = 0; i < len; i++)
! 	    {
! 	      tem = argv[bottom + i];
! 	      argv[bottom + i] = argv[middle + i];
! 	      argv[middle + i] = tem;
! 	      SWAP_FLAGS (bottom + i, middle + i);
! 	    }
! 	  /* Exclude the moved top segment from further swapping.  */
! 	  bottom += len;
! 	}
!     }

!   /* Update records for the slots the non-options now occupy.  */

!   first_nonopt += (optind - last_nonopt);
!   last_nonopt = optind;
! }

! Initialize the internal data when the first call is made.
!!$
!!$   function getopt_initialize (argc, argv, optstring)
!!$     int argc;
!!$     char *const *argv;
!!$     const char *optstring;
!!${
!!$  /* Start processing options with ARGV-element 1 (since ARGV-element 0
!!$     is the program name); the sequence of previously skipped
!!$     non-option ARGV-elements is empty.  */
!!$
!!$  first_nonopt = last_nonopt = optind;
!!$
!!$  nextchar = NULL;
!!$
!!$  posixly_correct = getenv ("POSIXLY_CORRECT");
!!$
!!$  /* Determine how to handle the ordering of options and nonoptions.  */
!!$
!!$  if (optstring[0] == '-')
!!$    {
!!$      ordering = RETURN_IN_ORDER;
!!$      ++optstring;
!!$    }
!!$  else if (optstring[0] == '+')
!!$    {
!!$      ordering = REQUIRE_ORDER;
!!$      ++optstring;
!!$    }
!!$  else if (posixly_correct != NULL)
!!$    ordering = REQUIRE_ORDER;
!!$  else
!!$    ordering = PERMUTE;
!!$
!!$#ifdef _LIBC
!!$  if (posixly_correct == NULL
!!$      && argc == original_argc && argv == original_argv)
!!$    {
!!$      if (nonoption_flags_max_len == 0)
!!$	{
!!$	  if (__getopt_nonoption_flags == NULL
!!$	      || __getopt_nonoption_flags[0] == '\0')
!!$	    nonoption_flags_max_len = -1;
!!$	  else
!!$	    {
!!$	      const char *orig_str = __getopt_nonoption_flags;
!!$	      int len = nonoption_flags_max_len = strlen (orig_str);
!!$	      if (nonoption_flags_max_len < argc)
!!$		nonoption_flags_max_len = argc;
!!$	      __getopt_nonoption_flags =
!!$		(char *) malloc (nonoption_flags_max_len);
!!$	      if (__getopt_nonoption_flags == NULL)
!!$		nonoption_flags_max_len = -1;
!!$	      else
!!$		memset (__mempcpy (__getopt_nonoption_flags, orig_str, len),
!!$			'\0', nonoption_flags_max_len - len);
!!$	    }
!!$	}
!!$      nonoption_flags_len = nonoption_flags_max_len;
!!$    }
!!$  else
!!$    nonoption_flags_len = 0;
!!$#endif
!!$
!!$  return optstring;
!!$}

!    Scan elements of ARGV (whose length is ARGC) for option characters
!    given in OPTSTRING.

!    If an element of ARGV starts with '-', and is not exactly "-" or "--",
!    then it is an option element.  The characters of this element
!    (aside from the initial '-') are option characters.  If `getopt'
!    is called repeatedly, it returns successively each of the option characters
!    from each of the option elements.

!    If `getopt' finds another option character, it returns that character,
!    updating `optind' and `nextchar' so that the next call to `getopt' can
!    resume the scan with the following option character or ARGV-element.

!    If there are no more option characters, `getopt' returns -1.
!    Then `optind' is the index in ARGV of the first ARGV-element
!    that is not an option.  (The ARGV-elements have been permuted
!    so that those that are not options now come last.)

!    OPTSTRING is a string containing the legitimate option characters.
!    If an option character is seen that is not listed in OPTSTRING,
!    return '?' after printing an error message.  If you set `opterr' to
!    zero, the error message is suppressed but we still return '?'.

!    If a char in OPTSTRING is followed by a colon, that means it wants an arg,
!    so the following text in the same ARGV-element, or the text of the following
!    ARGV-element, is returned in `optarg'.  Two colons mean an option that
!    wants an optional arg; if there is text in the current ARGV-element,
!    it is returned in `optarg', otherwise `optarg' is set to zero.

!    If OPTSTRING starts with `-' or `+', it requests different methods of
!    handling the non-option ARGV-elements.
!    See the comments about RETURN_IN_ORDER and REQUIRE_ORDER, above.

!    Long-named options begin with `--' instead of `-'.
!    Their names may be abbreviated as long as the abbreviation is unique
!    or is an exact match for some defined option.  If they have an
!    argument, it follows the option name in the same ARGV-element, separated
!    from the option name by a `=', or else the in next ARGV-element.
!    When `getopt' finds a long-named option, it returns 0 if that option's
!    `flag' field is nonzero, the value of the option's `val' field
!    if the `flag' field is zero.

!    The elements of ARGV aren't really const, because we permute them.
!    But we pretend they're const in the prototype to be compatible
!    with other systems.

!    LONGOPTS is a vector of `struct option' terminated by an
!    element containing a name which is zero.

!    LONGIND returns the index in LONGOPT of the long-named option found.
!    It is only valid when a long-named option has been found by the most
!    recent call.

!    If LONG_ONLY is nonzero, '-' as well as '--' can introduce
!    long-named options.

contains

!   subroutine getopt (optstring, nopt, opt, longopts, longind, long_only)
   subroutine getopt (optstring, optind, opt, optarg, longopts, longind, mpi )
      character(len=*), intent(in) :: optstring
      integer, intent(out) :: optind
      integer, intent(out) :: opt
      character(len=*), intent(out) :: optarg
      type(loption), dimension(:), target, intent(in), optional :: longopts
      integer, intent(inout), optional :: longind
!!!   If program is run under mpi, using mpirun, iargc is includes mpirun
!!!   arguments.
      logical, intent(in), optional :: mpi
!      logical, intent(in), optional :: long_only

!  The next char to be scanned in the option-element
!  in which the last option character we returned was found.
!  This allows us to pick up the scan where we left off.
!
!  If this is zero, or a null string, it means resume the scan
!  by advancing to the next ARGV-element.
      integer, save :: nextchar=0

      logical, save :: getopt_initialized = .false.
      logical :: print_errors
      integer, save :: argc, i, nameend
      logical :: hasequals, ambig, exact
      character(len=MAX_ARGLEN), dimension(:), allocatable, save :: argv
      character(len=MAX_ARGLEN) :: optname
      character(len=MAX_ARGLEN), save :: nextstr = ""
      character(len=1) :: c
      type(loption), pointer :: p, pfound
      integer :: indfound
      integer :: temp, temp_p1, temp_p2
      integer :: stderr = 6

      print_errors = opterr
      if (optstring(1:1) == ":") then
         print_errors = .false.
      end if

      if (command_argument_count() < 1) then
         opt = -1
         optind = 1
         return
      end if

      optarg = ""

      if ( .not. getopt_initialized ) then
         optind = 1
         !!! Need to use iargc()+1 to get the same result as with C
         argc = command_argument_count()+1
         if ( present(mpi) ) then
            if ( mpi ) then
               argc = argc - 4 ! Offset for mpirun -np X
            end if
         end if
         allocate ( argv(0:argc-1) )
         do i=0,argc-1
            !call getarg(i,argv(i))
            call get_command_argument(i,argv(i))
         end do
         getopt_initialized = .true.
      end if
!!$  if (optind == 0 || !__getopt_initialized)
!!$    {
!!$      if (optind == 0)
!!$	optind = 1;	/* Don't scan ARGV[0], the program name.  */
!!$      optstring = _getopt_initialize (argc, argv, optstring);
!!$      __getopt_initialized = 1;
!!$    }

!!$!    Test whether ARGV[optind] points to a non-option argument.
!!$!    Either it does not have option syntax, or there is an environment flag
!!$!    from the shell indicating it is not an option.  The later information
!!$!    is only used when the used in the GNU libc.
!!$#ifdef _LIBC
!!$# define NONOPTION_P (argv[optind][0] != '-' || argv[optind][1] == '\0'	      \
!!$		      || (optind < nonoption_flags_len			      \
!!$			  && __getopt_nonoption_flags[optind] == '1'))
!!$#else
!!$# define NONOPTION_P (argv[optind][0] != '-' || argv[optind][1] == '\0')
!!$#endif

!!!   nextchar is an index into nextstring.
!!!   Need to use max(1,nextchar) to avoid bounds error if nextchar=0
!!!   Short circuit of if isn't guaranteed.
      if (nextchar == 0 .or. len_trim(nextstr(max(1,nextchar):)) == 0 ) then

         ! Advance to the next ARGV-element.

         ! Give FIRST_NONOPT & LAST_NONOPT rational values if OPTIND has been
	 ! moved back by the user (who may also have changed the arguments).
!!$         if (last_nonopt > optind) then
!!$            last_nonopt = optind
!!$         end if
!!$         if (first_nonopt > optind) then
!!$            first_nonopt = optind
!!$         end if

!!$      if (ordering == PERMUTE)
!!$	{
!!$	  /* If we have just processed some options following some non-options,
!!$	     exchange them so that the options come first.  */
!!$
!!$	  if (first_nonopt != last_nonopt && last_nonopt != optind)
!!$	    exchange ((char **) argv);
!!$	  else if (last_nonopt != optind)
!!$	    first_nonopt = optind;
!!$
!!$	  /* Skip any additional non-options
!!$	     and extend the range of non-options previously skipped.  */
!!$
!!$	  while (optind < argc && NONOPTION_P)
!!$	    optind++;
!!$	  last_nonopt = optind;
!!$	}

         ! The special ARGV-element `--' means premature end of options.
	 ! Skip it like a null option,
	 ! then exchange with previous non-options as if it were an option,
	 ! then skip everything else like a non-option.

!!!      Split if tests to avoid bounds errors when optind=argc
         if (optind /= argc ) then
            if ( argv(optind) == "--") then
               optind = optind+1

!!$	  if (first_nonopt != last_nonopt && last_nonopt != optind)
!!$	    exchange ((char **) argv);
!!$	  else if (first_nonopt == last_nonopt)
!!$	    first_nonopt = optind;
!!$	  last_nonopt = argc;
!!$
!!$	  optind = argc;
!!$	}
            end if
         end if

      ! If we have done all the ARGV-elements, stop the scan
      ! and back over any non-options that we skipped and permuted.

         if (optind == argc) then
!!$	  /* Set the next-arg-index to point at the non-options
!!$	     that we previously skipped, so the caller will digest them.  */
!!$	  if (first_nonopt != last_nonopt)
!!$	    optind = first_nonopt;
            opt = -1
            return
         end if

      !  If we have come to a non-option and did not permute it,
      !  either stop the scan or describe it to the caller and pass it by.

         if ( argv(optind)(1:1) /= "-" ) then
            opt = -1
            return
         end if

      !  We have found another option-ARGV-element.
      !  Skip the initial punctuation.

         nextchar = 2
         nextstr = argv(optind)
         if (present(longopts) .and. argv(optind)(2:2) == "-") then
            nextchar = 3
         end if

      end if


  ! Decode the current option-ARGV-element.

  !  Check whether the ARGV-element is a long option.

  !  If long_only and the ARGV-element has the form "-f", where f is
  !  a valid short option, don't consider it an abbreviated form of
  !  a long option that starts with f.  Otherwise there would be no
  !  way to give the -f short option.

  !  On the other hand, if there's a long option "fubar" and
  !  the ARGV-element is "-fu", do consider that an abbreviation of
  !  the long option, just like "--fu", and not "-f" with arg "u".

  !  This distinction seems to be the most useful approach.

!!! Note long_only not implemented properly here

      if ( present(longopts) .and. argv(optind)(2:2) == "-") then
!!$  if (longopts != NULL
!!$      && (argv[optind][1] == '-'
!!$	  || (long_only && (argv[optind][2] || !my_index (optstring, argv[optind][1])))))
!!$    {
!!$      char *nameend;
!!$      const struct option *p;
!!$      const struct option *pfound = NULL;
!!$      int exact = 0;
!!$      int ambig = 0;
!!$      int indfound = -1;
!!$      int option_index;

!!!  Find end of option name, either end of string or =
         nameend = index(nextstr,"=")
         if (nameend /= 0) then
            hasequals = .true.
            ! Remove the trailing =
            optname = nextstr(nextchar:nameend-1)
         else
            hasequals = .false.
            nameend = len_trim(argv(optind))
            optname = nextstr(nextchar:nameend)
         end if

         ! Test all long options for either exact match
         ! or abbreviated matches.  */
         nullify(pfound)
         ambig = .false.
         exact = .false.
         indfound = -1
         do i=1,size(longopts)
            p => longopts(i)
            ! To allow partial matching use index rather than ==
            if ( index(p%name, optname) == 1 ) then
               if ( len_trim(p%name) == len_trim(optname) ) then
                  ! Exact match found.
                  pfound => p
                  indfound = i
                  exact = .true.
                  exit
               else if ( .not. associated(pfound) ) then
                  ! First nonexact match found.
                  pfound => p
                  indfound = i
               else if ( pfound%has_arg /= p%has_arg .or. &
                    pfound%val /= p%val ) then
                  ! Second or later nonexact match found
                  ambig = .true.
               end if
            end if
         end do

         if (ambig .and. .not. exact ) then
            if ( print_errors ) then
               write(unit=stderr,fmt="(a,a,a,a)") trim(argv(0)), ": option ", &
                              trim(argv(optind)), " is ambiguous"
               nextchar = len_trim(nextstr)+1 ! Set to end
            end if
            optind = optind + 1
            optopt = 0
            opt = ichar("?")
            return
         end if

         if (associated(pfound)) then
            optind = optind + 1
            if (hasequals) then ! i.e. it has an argument
               if ( pfound%has_arg > 0 ) then
                  ! Go past the =
                  optarg = nextstr(nameend+1:)
               else
                  if (print_errors) then
                     if (argv(optind-1)(2:2) == "-") then
                        ! --option
                        write(unit=stderr,fmt="(a,a,a,a)") trim(argv(0)), &
                 ": option --", trim(pfound%name), " doesn't allow an argument"
                     else
                        ! +option or -option
                        write(unit=stderr,fmt="(a,a,a,a,a)") trim(argv(0)), &
                             ": option --", argv(optind-1)(1:1), &
                              trim(pfound%name), " doesn't allow an argument"
                     end if
                  end if
                  nextchar = len_trim(nextstr)+1 ! Set to end
                  optopt = pfound%val
                  opt = ichar("?")
                  return
               end if
            else if (pfound%has_arg == 1) then
               if (optind < argc) then
                  optarg = argv(optind)
                  optind = optind+1
               else
                  if (print_errors) then
                     write(unit=stderr,fmt="(a,a,a,a)") trim(argv(0)), &
                      ": option ", trim(argv(optind-1)), " requires an argument"
                  end if
                  nextchar = len_trim(nextstr)+1 ! Set to end
                  optopt = pfound%val
                  if ( optstring(1:1) == ":" ) then
                     opt = ichar(":")
                  else
                     opt = ichar("?")
                  end if
                  return
               end if
            end if
            nextchar = len_trim(nextstr)+1 ! Set to end
            if (present(longind)) then
               longind = indfound
            end if
            opt = pfound%val
            return
         end if

      !  Can't find it as a long option.  If this is not getopt_long_only,
      !  or the option starts with '--' or is not a valid short
      !  option, then it's an error.
      !	 Otherwise interpret it as a short option.  */
!!$      if (!long_only || argv[optind][1] == '-'
!!$	  || my_index (optstring, *nextchar) == NULL)
         if ( argv(optind)(2:2) == "-" .or. index(optstring,nextstr(nextchar:nextchar)) == 0 ) then
            if (print_errors) then
               if (argv(optind)(2:2) == "-") then
                  ! --option
                  write(unit=stderr,fmt="(a,a,a)") trim(argv(0)), &
                    ": unrecognized option --", trim(optname)
               else
                  ! +option or -option
                  write(unit=stderr,fmt="(a,a,a,a)") trim(argv(0)), &
                    ": unrecognized option ",argv(optind)(1:1), trim(optname)
               end if
            end if
            nextchar = 0
            optind = optind + 1
            optopt = 0
            opt = ichar("?")
            return
         end if

      end if  ! if (present(longopts)

  ! Look at and handle the next short option-character.
      c = nextstr(nextchar:nextchar)
      nextchar = nextchar+1
      temp = index(optstring, c)

      ! Increment `optind' when we start to process its last character.
      if ( len_trim(nextstr(nextchar:)) == 0 ) then
         optind = optind+1
      end if

      if (temp == 0 .or. c == ":") then
         if (print_errors) then
            write(unit=stderr,fmt="(a,a,a)") trim(argv(0)), ": invalid option -- ", c
         end if
         optopt = ichar(c)
         opt = ichar("?")
         return
      end if

!!$    ! Convenience. Treat POSIX -W foo same as long option --foo
!!$            if (optstring(temp:temp) == "W" .and. &
!!$                optstring(temp+1:temp+1) == ";") then
!!$               nullify(pfound)
!!$               exact = .false.
!!$	int ambig = 0;
!!$	int indfound = 0;
!!$	int option_index;
!!$
!!$	/* This is an option that requires an argument.  */
!!$	if (*nextchar != '\0')
!!$	  {
!!$	    optarg = nextchar;
!!$	    /* If we end this ARGV-element by taking the rest as an arg,
!!$	       we must advance to the next element now.  */
!!$	    optind++;
!!$	  }
!!$	else if (optind == argc)
!!$	  {
!!$	    if (print_errors)
!!$	      {
!!$		/* 1003.2 specifies the format of this message.  */
!!$		fprintf (stderr, _("%s: option requires an argument -- %c\n"),
!!$			 argv[0], c);
!!$	      }
!!$	    optopt = c;
!!$	    if (optstring[0] == ':')
!!$	      c = ':';
!!$	    else
!!$	      c = '?';
!!$	    return c;
!!$	  }
!!$	else
!!$	  /* We already incremented `optind' once;
!!$	     increment it again when taking next ARGV-elt as argument.  */
!!$	  optarg = argv[optind++];
!!$
!!$	/* optarg is now the argument, see if it's in the
!!$	   table of longopts.  */
!!$
!!$	for (nextchar = nameend = optarg; *nameend && *nameend != '='; nameend++)
!!$	  /* Do nothing.  */ ;
!!$
!!$	/* Test all long options for either exact match
!!$	   or abbreviated matches.  */
!!$	for (p = longopts, option_index = 0; p->name; p++, option_index++)
!!$	  if (!strncmp (p->name, nextchar, nameend - nextchar))
!!$	    {
!!$	      if ((unsigned int) (nameend - nextchar) == strlen (p->name))
!!$		{
!!$		  /* Exact match found.  */
!!$		  pfound = p;
!!$		  indfound = option_index;
!!$		  exact = 1;
!!$		  break;
!!$		}
!!$	      else if (pfound == NULL)
!!$		{
!!$		  /* First nonexact match found.  */
!!$		  pfound = p;
!!$		  indfound = option_index;
!!$		}
!!$	      else
!!$		/* Second or later nonexact match found.  */
!!$		ambig = 1;
!!$	    }
!!$	if (ambig && !exact)
!!$	  {
!!$	    if (print_errors)
!!$	      fprintf (stderr, _("%s: option `-W %s' is ambiguous\n"),
!!$		       argv[0], argv[optind]);
!!$	    nextchar += strlen (nextchar);
!!$	    optind++;
!!$	    return '?';
!!$	  }
!!$	if (pfound != NULL)
!!$	  {
!!$	    option_index = indfound;
!!$	    if (*nameend)
!!$	      {
!!$		/* Don't test has_arg with >, because some C compilers don't
!!$		   allow it to be used on enums.  */
!!$		if (pfound->has_arg)
!!$		  optarg = nameend + 1;
!!$		else
!!$		  {
!!$		    if (print_errors)
!!$		      fprintf (stderr, _("\
!!$%s: option `-W %s' doesn't allow an argument\n"),
!!$			       argv[0], pfound->name);
!!$
!!$		    nextchar += strlen (nextchar);
!!$		    return '?';
!!$		  }
!!$	      }
!!$	    else if (pfound->has_arg == 1)
!!$	      {
!!$		if (optind < argc)
!!$		  optarg = argv[optind++];
!!$		else
!!$		  {
!!$		    if (print_errors)
!!$		      fprintf (stderr,
!!$			       _("%s: option `%s' requires an argument\n"),
!!$			       argv[0], argv[optind - 1]);
!!$		    nextchar += strlen (nextchar);
!!$		    return optstring[0] == ':' ? ':' : '?';
!!$		  }
!!$	      }
!!$	    nextchar += strlen (nextchar);
!!$	    if (longind != NULL)
!!$	      *longind = option_index;
!!$	    if (pfound->flag)
!!$	      {
!!$		*(pfound->flag) = pfound->val;
!!$		return 0;
!!$	      }
!!$	    return pfound->val;
!!$	  }
!!$	  nextchar = NULL;
!!$	  return 'W';	/* Let the application handle it.   */
!!$      }

!   Does optstring have something appended to ensure this isn't off the end???
      temp_p1 = min( temp+1, len(optstring) )
      if (optstring(temp_p1:temp_p1) == ":" .and. temp_p1 == temp+1) then
         temp_p2 = min( temp+2, len(optstring) )
         if (optstring(temp_p2:temp_p2) == ":" .and. temp_p2 == temp+2) then
            ! This is an option that accepts an argument optionally.
            if (len_trim(nextstr(nextchar:)) /= 0 ) then
               optarg = trim(nextstr(nextchar:))
               optind = optind + 1
            else
               optarg = ""
               nextchar = 0
            end if
         else
            ! This is an option that requires an argument.
            if (len_trim(nextstr(nextchar:)) /= 0 ) then
               optarg = trim(nextstr(nextchar:))
               ! If we end this ARGV-element by taking the rest as an arg,
               ! we must advance to the next element now.
               optind = optind + 1
            else if (optind == argc) then
               if (print_errors) then
                  write(unit=stderr,fmt="(a,a,a)") trim(argv(0)), ": option requires an argument -- ", c
               end if
               optopt = ichar(c)
               if (optstring(1:1) == ":") then
                  c = ":"
               else
                  c = "?"
               end if
            else
               ! We already incremented `optind' once;
               ! increment it again when taking next ARGV-elt as argument.
               optarg = argv(optind)
               optind = optind + 1
            end if
            nextchar = 0  ! Where does this belong, perhaps on next line?
         end if
      end if
      opt = ichar(c)
      return
   end subroutine getopt

   subroutine getcline ( cline )

!     Get the complete program command line
      character(len=*), intent(out) :: cline
      integer :: iarg
      character(len=MAX_ARGLEN) :: arg

      cline = ''
      do iarg=0,command_argument_count()
         !call getarg(iarg,arg)
         call get_command_argument(iarg,arg) 
!        Use >= here to allow for the extra space
         if ( len_trim(cline) + len_trim(arg) >= len(cline) ) then
            print*, "Error, increase length of command line variable"
            stop
         end if
         cline = cline(1:len_trim(cline)) // " " // trim(arg)
      end do

      !  The loop above adds a leading blank so adjustl
      cline = adjustl(cline)

   end subroutine getcline

end module getopt_m

