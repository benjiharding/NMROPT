module makepar_mod
   implicit none

contains

   subroutine makepar(parfile)
      character(256), intent(in) :: parfile

      open (unit=99, file=parfile, action="write", status="replace")

      write (99, "(A)") '               PARAMETERS FOR NMROPT'
      write (99, "(A)") '               *********************'
      write (99, "(A)") ''
      write (99, "(A)") 'START OF PARAMETERS:'
      write (99, "(A)") 'data.dat                          - file with data'
      write (99, "(A)") '1 4 5 6 7 11                      - columns for dh, x, y, z, var and wt'
      write (99, "(A)") '0                                 - normal score transform var? (0=no, 1=yes)'
      write (99, "(A)") '-1.0e21    1.0e21                 - trimming limits'
      write (99, "(A)") '100                               - number of unconditional realizations'
      write (99, "(A)") '0                                 - simulation type (0=LU (<2500 data), 1=sequential)'
      write (99, "(A)") '5841044                           - random number seed'
      write (99, "(A)") '1 -1                              - debugging level, realization to output (-1 for all)'
      write (99, "(A)") 'nmropt.dbg                        - file for debugging output'
      write (99, "(A)") 'nmropt.out                        - file for network output'
      write (99, "(A)") 'nmrwts.out                        - file for optimized network weights'
      write (99, "(A)") 'nmrobj.out                        - file for objective function value per iteration'
      write (99, "(A)") 'nmr_                              - prefix for target/experimental output files'
      write (99, "(A)") '3                                 - number of network layers (input to output layer)'
      write (99, "(A)") '5 5 1                             - network layer dimensions (input + nugget to output layer)'
      write (99, "(A)") '1 0.1                             - network wt regularization (0=none, 1=L1, 2=L2), lambda'
      write (99, "(A)") 'pool.dat                          - file with covariance structs. of Gaussian pool'
      write (99, "(A)") '0                                 - consider factor precedence? (0=no, 1=yes)'
      write (99, "(A)") '1  1  1  1                        - objective components: varg, ivarg, runs, npoint'
      write (99, "(A)") '1  1  1  1                        - objective weight:     varg, ivarg, runs, npoint'
      write (99, "(A)") '3                                 - number of indicator thresholds'
      write (99, "(A)") '-1.28 0.0 1.28                    - Gaussian indicator thresholds'
      write (99, "(A)") '1  1  1                           - threshold weights'
      write (99, "(A)") '1                                 - runs above or below threshold? (0=below, 1=above, 2=both)'
      write (99, "(A)") '30                                - max number of runs to consider'
      write (99, "(A)") '0                                 - runs target from file? (0=no, 1=yes)'
      write (99, "(A)") 'target_runs.out                   - runs target file'
      write (99, "(A)") '1                                 - npoint above or below threshold? (0=below, 1=above)'
      write (99, "(A)") '30                                - max number of connected steps to consider'
      write (99, "(A)") '0                                 - npoint target from file? (0=no, 1=yes)'
      write (99, "(A)") 'target_npoint.out                 - npoint target file'
      write (99, "(A)") '0.8 0.5 1.0 15 1000               - DE parameters: F, CR lo, CR hi, pop. size, its'
      write (99, "(A)") '0.0 1.0                           - DE bounds: lower, upper'
      write (99, "(A)") 'omega.out                         - file with factor omega bounds'
      write (99, "(A)") '1                                 - num. threads for parallel DE (1=serial, -1 for all)'
      write (99, "(A)") '2                                 - number of experimental variogram directions'
      write (99, "(A)") '0.0 22.5 1000 0.0 22.5 1000 0.0   - dir 01: azm,azmtol,bandhorz,dip,diptol,bandvert,tilt'
      write (99, "(A)") '8  1000.0  500.0                  -        number of lags,lag distance,lag tolerance'
      write (99, "(A)") '90. 22.5 1000 0.0 22.5 1000 0.0   - dir 02: azm,azmtol,bandhorz,dip,diptol,bandvert,tilt'
      write (99, "(A)") '8  1000.0  500.0                  -        number of lags,lag distance,lag tolerance'
      write (99, "(A)") '1                                 - number of target variogram models'
      write (99, "(A)") '3                                 - number of target indicator variogram models'
      write (99, "(A)") '1.0                               - IDW power for variogram optimization weighting'
      write (99, "(A)") '999999                            - max number of exp. variogram pairs (see note)'
      write (99, "(A)") '1                                 - standardize sill? (0=no, 1=yes)'
      write (99, "(A)") '1    0.1                          - nst, nugget effect'
      write (99, "(A)") '1    0.9  0.0   0.0   0.0         - it,cc,ang1,ang2,ang3'
      write (99, "(A)") '        10.0  10.0  10.0          - a_hmax, a_hmin, a_vert'
      write (99, "(A)") '1    0.1                          - inst, nugget effect'
      write (99, "(A)") '1    0.9  0.0   0.0   0.0         - iit,icc,iang1,iang2,iang3'
      write (99, "(A)") '         10.0  10.0  10.0         - ia_hmax, ia_hmin, ia_vert'
      write (99, "(A)") '1    0.1                          - inst, nugget effect'
      write (99, "(A)") '1    0.9  0.0   0.0   0.0         - iit,icc,iang1,iang2,iang3'
      write (99, "(A)") '         10.0  10.0  10.0         - ia_hmax, ia_hmin, ia_vert'
      write (99, "(A)") '1    0.1                          - inst, nugget effect'
      write (99, "(A)") '1    0.9  0.0   0.0   0.0         - iit,icc,iang1,iang2,iang3'
      write (99, "(A)") '         10.0  10.0  10.0         - ia_hmax, ia_hmin, ia_vert'
      close (99)

   end subroutine makepar

   subroutine chknam(str, len)
      !-----------------------------------------------------------------------
      !
      !                   Check for a Valid File Name
      !                   ***************************
      !
      ! This subroutine takes the character string "str" of length "len" and
      ! removes all leading blanks and blanks out all characters after the
      ! first blank found in the string (leading blanks are removed first).
      !
      !
      !
      !-----------------------------------------------------------------------
      character(len=*), intent(inout) :: str
      integer itrim, len
      !
      ! Remove leading blanks:
      str = adjustl(str)
      !
      ! find first two blanks and blank out remaining characters:
      itrim = index(str, '   ')
      if (itrim > 0) str(itrim:) = ' '
      !
      ! Look for "-fi"
      itrim = index(str, '-fi')
      if (itrim > 0) str(itrim:) = ' '
      !
      ! Look for "\fi"
      itrim = index(str, '\fi')
      if (itrim > 0) str(itrim:) = ' '

      return
   end subroutine chknam

end module makepar_mod
