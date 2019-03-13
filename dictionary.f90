!Copyright (c) 2008, Arjen Markus
!
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification,
!are permitted provided that the following conditions are met:
!
!Redistributions of source code must retain the above copyright notice,
!this list of conditions and the following disclaimer.
!Redistributions in binary form must reproduce the above copyright notice,
!this list of conditions and the following disclaimer in the documentation
!and/or other materials provided with the distribution.
!Neither the name of the author nor the names of the contributors
!may be used to endorse or promote products derived from this software
!without specific prior written permission.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
!THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! dictionary.f90 --
!     Include file for defining dictionaries:
!     a mapping of strings to some data
!
!     See the example/test program for the way to use this
!
!     Note:
!     Use is made of a hash table. This should speed up most
!     operations. The algorithm for determining the hashkey
!     is taken from Kernighan and Pike: The Practice of Programming
!
!     Note:
!     - Define the length of the strings as
!       parameter "DICT_KEY_LENGTH"
!     - Define a derived type for the data
!       to be stored
!     - Also define a "null" value - DICT_NULL
!       of type DICT_DATA, for use when the
!       key is not found.
!     - Put both in a separate module, that
!       will be used.
!
!     $Id: dictionary.f90,v 1.3 2007/01/26 09:56:43 arjenmarkus Exp $
!
module dictionary

use dict_m, DICT_DATA => MYDATA
use linkedlist
implicit none

private

type HASH_LIST
    type(LINKED_LIST), pointer :: list
end type HASH_LIST

type DICT_STRUCT
    private
    type(HASH_LIST), pointer, dimension(:) :: table
end type DICT_STRUCT

public :: DICT_DATA
public :: DICT_STRUCT
public :: dict_get_key
public :: dict_create
public :: dict_add_key

integer, parameter :: hash_size  = 4993
integer, parameter :: multiplier = 31

contains
!
! Routines and functions specific to dictionaries
!

! dict_create --
!     Create and initialise a dictionary
! Arguments:
!     dict       Pointer to new dictionary
!     key        Key for the first element
!     value      Value for the first element
! Note:
!     This version assumes a shallow copy is enough
!     (that is, there are no pointers within the data
!     to be stored)
!     It also assumes the argument list does not already
!     refer to a list. Use dict_destroy first to
!     destroy up an old list.
!
subroutine dict_create( dict, key, value )
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) ::  key
    type(DICT_DATA), intent(in)  :: value

    type(LIST_DATA)              :: data
    integer                      :: i
    integer                      :: hash

    allocate( dict )
    allocate( dict%table(hash_size) )

    do i = 1,hash_size
        dict%table(i)%list => null()
    enddo

    data%key   = key
    data%value = value

    hash = dict_hashkey( trim(key ) )
    call list_create( dict%table(hash)%list, data )

end subroutine dict_create

! dict_destroy --
!     Destroy an entire dictionary
! Arguments:
!     dict       Pointer to the dictionary to be destroyed
! Note:
!     This version assumes that there are no
!     pointers within the data that need deallocation
!
subroutine dict_destroy( dict )
    type(DICT_STRUCT), pointer  :: dict

    integer                     :: i

    do i = 1,size(dict%table)
        if ( associated( dict%table(i)%list ) ) then
            call list_destroy( dict%table(i)%list )
        endif
    enddo
    deallocate( dict%table )
    deallocate( dict )

end subroutine dict_destroy

! dict_add_key
!     Add a new key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for the new element
!     value      Value for the new element
! Note:
!     If the key already exists, the
!     key's value is simply replaced
!
subroutine dict_add_key( dict, key, value )
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key
    type(DICT_DATA), intent(in)  :: value

    type(LIST_DATA)              :: data
    type(LINKED_LIST), pointer   :: elem
    integer                      :: hash

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
        elem%data%value = value
    else
        data%key   = key
        data%value = value
        hash       = dict_hashkey( trim(key) )
        if ( associated( dict%table(hash)%list ) ) then
            call list_insert( dict%table(hash)%list, data )
        else
            call list_create( dict%table(hash)%list, data )
        endif
    endif

end subroutine dict_add_key

! dict_delete_key
!     Delete a key-value pair from the dictionary
! Arguments:
!     dict       Dictionary in question
!     key        Key to be removed
!
subroutine dict_delete_key( dict, key )
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key

    type(LINKED_LIST), pointer   :: elem
    integer                      :: hash

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
        hash = dict_hashkey( trim(key) )
        call list_delete_element( dict%table(hash)%list, elem )
    endif
end subroutine dict_delete_key

! dict_get_key
!     Get the value belonging to a key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for which the values are sought
!
function dict_get_key( dict, key ) result(value)
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key
    type(DICT_DATA)              :: value

    type(LIST_DATA)              :: data
    type(LINKED_LIST), pointer   :: elem

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
        value = elem%data%value
    else
        value = DICT_NULL
    endif
end function dict_get_key

! dict_has_key
!     Check if the dictionary has a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!
function dict_has_key( dict, key ) result(has)
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key
    logical                      :: has

    type(LINKED_LIST), pointer   :: elem

    elem => dict_get_elem( dict, key )

    has = associated(elem)
end function dict_has_key

! dict_get_elem
!     Find the element with a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!
function dict_get_elem( dict, key ) result(elem)
    type(DICT_STRUCT), pointer   :: dict
    character(len=*), intent(in) :: key

    type(LINKED_LIST), pointer   :: elem
    integer                      :: hash

    hash = dict_hashkey( trim(key) )

    elem => dict%table(hash)%list
    do while ( associated(elem) )
        if ( elem%data%key .eq. key ) then
            exit
        else
            elem => list_next( elem )
        endif
    enddo
end function dict_get_elem

! dict_hashkey
!     Determine the hash value from the string
! Arguments:
!     key        String to be examined
!
integer function dict_hashkey( key )
    character(len=*), intent(in) :: key

    integer                      :: hash
    integer                      :: i
    integer                      :: max_hashkey=(huge(1)-256)/multiplier

    dict_hashkey = 0

    do i = 1,len(key)
        ! Fortran integers are signed, prevent from overflowing to negative
        dict_hashkey = mod(dict_hashkey,max_hashkey)
        dict_hashkey = multiplier * dict_hashkey + ichar(key(i:i))
    enddo

    dict_hashkey = 1 + mod( dict_hashkey-1, hash_size )
end function dict_hashkey

end module dictionary
