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

module dict_m
implicit none

public

!
! The length of the keys
!
integer, parameter :: DICT_KEY_LENGTH = 60

!
! The data that will be stored with each key
type MYDATA
    integer :: index
end type MYDATA

!
! The "null" value for these data
!
type(MYDATA), parameter :: DICT_NULL = mydata( 0 )

type LIST_DATA
    character(len=DICT_KEY_LENGTH) :: key
    type(MYDATA)                   :: value
end type LIST_DATA

end module dict_m
