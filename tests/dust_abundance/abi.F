      program dust_parameter_abi
      use iso_c_binding
      implicit none
#include "grackle.def"
      type(grackle_chemistry_data) chemistry
      integer(c_int) status
      interface
         integer(c_size_t) function grackle_parameter_bytes() bind(C)
            use iso_c_binding
         end function
      end interface
      if (c_sizeof(chemistry) .ne. grackle_parameter_bytes()) stop 1
      status = set_default_chemistry_parameters(chemistry)
      if (status .ne. 1) stop 2
      if (chemistry%dust_chemistry .ne. 0) stop 3
      if (chemistry%h2_on_dust .ne. 0) stop 4
      print *, 'PASS Fortran/C chemistry_data size and defaults'
      end
