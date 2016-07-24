module Parameters_mod
#include "config.h"

  implicit none
  !everything public

  integer, parameter :: SP=KIND(1.0)
  integer, parameter :: DP=KIND(1.0D0)
  integer, parameter :: RP=DP

  integer, parameter   :: STRING_LENGTH=1024


end module Parameters_mod
